#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TString.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

void CountCE(){
    const TString filename = "output0_t*.root";
    const TString outname = "cut_blocks.root";
    const TString cut = "didChargeExchange>0 && fZ<7000";
    const double gapMm = 300.0;
    const bool Verbose = false; // print the z bins, the density filter, the blocks found and per-block counts

    TChain *ce = new TChain("ChargeExchange_values");
    TChain *info = new TChain("ParticleInfo");
    ce->Add(filename);
    info->Add(filename);
    if (ce->GetEntries() == 0 || info->GetEntries() == 0) return;

    ROOT::EnableImplicitMT(128);

    // CE z positions and (event, track id) of every CE
    std::vector<double> z;
    std::set<std::pair<long,long>> ceKeys;
    {
        // per-thread-slot containers merged afterwards, so no global lock; z is sorted below
        ROOT::RDataFrame dfCE(*ce);
        std::vector<std::vector<double>> slotZ(dfCE.GetNSlots());
        std::vector<std::set<std::pair<long,long>>> slotKeys(dfCE.GetNSlots());
        dfCE.Filter(cut.Data()).ForeachSlot([&](unsigned slot, Int_t ev, Int_t trk, Double_t fz){
            slotZ[slot].push_back(fz);
            slotKeys[slot].insert({ev, trk});
        }, {"fEvent", "TrackID", "fZ"});
        for (auto &v : slotZ) z.insert(z.end(), v.begin(), v.end());
        for (auto &k : slotKeys) ceKeys.insert(k.begin(), k.end());
    }
    if (z.empty()) return;
    std::sort(z.begin(), z.end());

    // CEs in the thin silicon tracker planes between the carbon blocks fill the gaps and chain every block into
    // one. Silicon and carbon bins hold about the same counts, so density cannot tell them apart. Thickness can:
    // a carbon block is a long run of consecutive occupied z bins and a tracker plane is a short one. Keep only
    // runs at least half as long as the longest run.
    {
        const double binMm = 10.0;
        std::map<long,long> bins;
        for (double v : z) bins[std::lround(std::floor(v / binMm))]++;

        std::vector<std::vector<long>> runs; // bin indices of each run of consecutive occupied bins
        for (auto &b : bins) {
            if (runs.empty() || b.first != runs.back().back() + 1) runs.push_back({});
            runs.back().push_back(b.first);
        }
        size_t longest = 0;
        for (auto &r : runs) longest = std::max(longest, r.size());
        const double minRun = 0.5 * longest;
        std::set<long> keptBins;
        for (auto &r : runs) if (r.size() >= minRun) keptBins.insert(r.begin(), r.end());

        if (Verbose) {
            std::cout << "CE (" << cut << "): " << z.size() << " total, z " << z.front() << " to " << z.back() << " mm\n"
                      << "z bins (" << binMm << " mm): " << bins.size() << " occupied in " << runs.size()
                      << " runs, longest " << longest << " bins, keep runs >= " << minRun << " bins\n";
            for (auto &r : runs) {
                std::cout << "  run z " << r.front() * binMm << ".." << (r.back() + 1) * binMm << ": " << r.size()
                          << " bins" << (r.size() >= minRun ? "" : "   DROPPED") << "\n";
            }
        }
        std::vector<double> zKept;
        for (double v : z) if (keptBins.count(std::lround(std::floor(v / binMm)))) zKept.push_back(v);
        z.swap(zKept);
        if (Verbose) std::cout << "after thickness filter: " << z.size() << " CE kept\n";
        if (z.empty()) return;
    }

    // One pass over ParticleInfo, nothing big held in memory:
    //  - tracker z (boundaryType 1) of every particle, to find the tracker planes
    //  - tracker hits and detector entries (boundaryType 2) of CE protons only
    //    (event, parent track id) -> tracker z / detector (x, y)
    std::set<long> trackerZ;
    std::map<std::tuple<long,long,long>, std::pair<double,double>> protonHits; // (event, parent, tracker z) -> (x, y)
    std::map<std::pair<long,long>, std::pair<double,double>> protonDetectorHits;
    {
        // each thread fills its own containers (no lock per row), merged after the pass
        struct Local {
            std::set<long> trackerZ;
            std::map<std::tuple<long,long,long>, std::pair<double,double>> protonHits; // (event, parent, tracker z) -> (x, y)
            std::map<std::pair<long,long>, std::pair<double,double>> protonDetectorHits;
        };
        std::vector<Local> locals(std::max(1u, ROOT::GetThreadPoolSize()));
        ROOT::RDataFrame(*info).ForeachSlot([&](unsigned slot, Int_t ev, Int_t parent, Int_t pdg, Int_t bt,
                                                Double_t fx, Double_t fy, Double_t fz){
            Local &l = locals[slot];
            if (bt == 1) l.trackerZ.insert(std::lround(fz));
            if (pdg != 2212 || !ceKeys.count({ev, parent})) return;
            if (bt == 1) l.protonHits[{ev, parent, std::lround(fz)}] = {fx, fy};
            else if (bt == 2) l.protonDetectorHits[{ev, parent}] = {fx, fy};
        }, {"fEvent", "ParentID", "pdg", "boundaryType", "fX", "fY", "fZ"});

        for (auto &l : locals) {
            trackerZ.insert(l.trackerZ.begin(), l.trackerZ.end());
            protonHits.insert(l.protonHits.begin(), l.protonHits.end());
            for (auto &kv : l.protonDetectorHits) protonDetectorHits[kv.first] = kv.second;
        }
    }

    std::vector<std::pair<size_t,size_t>> blocks;
    size_t start = 0;
    for (size_t i = 1; i < z.size(); i++) {
        if (z[i] - z[i-1] > gapMm) {
            blocks.push_back({start, i - 1});
            start = i;
        }
    }
    blocks.push_back({start, z.size() - 1});
    if (Verbose) {
        std::cout << blocks.size() << " blocks (gap > " << gapMm << " mm splits):\n";
        for (size_t b = 0; b < blocks.size(); b++) {
            std::cout << "  block " << b + 1 << ": z " << z[blocks[b].first] << " to " << z[blocks[b].second] << ", "
                      << blocks[b].second - blocks[b].first + 1 << " CE\n";
        }
        for (size_t i = 1; i < z.size(); i++) {
            if (z[i] - z[i-1] > gapMm) std::cout << "  gap " << z[i-1] << " -> " << z[i] << "\n";
        }
    }

    // group hits into planes: faces of one plane sit mm apart, planes sit hundreds of mm apart.
    // The split size is taken from the data as the biggest jump in the sorted z gaps.
    std::vector<long> zs(trackerZ.begin(), trackerZ.end());
    std::vector<double> gaps;
    for (size_t i = 1; i < zs.size(); i++) gaps.push_back(zs[i] - zs[i-1]);
    std::vector<double> sortedGaps(gaps);
    std::sort(sortedGaps.begin(), sortedGaps.end());
    double split = 0;
    double bestRatio = 5;
    for (size_t i = 1; i < sortedGaps.size(); i++) {
        double ratio = sortedGaps[i] / sortedGaps[i-1];
        if (ratio > bestRatio) { bestRatio = ratio; split = std::sqrt(sortedGaps[i] * sortedGaps[i-1]); }
    }
    
    // [lowest z, highest z] of each tracker plane
    std::vector<std::pair<long,long>> planes;
    for (size_t i = 0; i < zs.size(); i++) {
        if (i == 0 || gaps[i-1] > split) planes.push_back({zs[i], zs[i]});
        else planes.back().second = zs[i];
    }

    if (Verbose) {
        std::cout << trackerZ.size() << " distinct tracker z, plane split " << split << " mm, " << planes.size() << " planes:";
        for (auto &pl : planes) std::cout << " [" << pl.first << "," << pl.second << "]";
        std::cout << "\n";
    }

    //out file
    TFile *out = TFile::Open(outname, "RECREATE");
    for (size_t b = 0; b < blocks.size(); b++) {
        double lo = z[blocks[b].first];
        double hi = z[blocks[b].second];
        TString sel = Form("(%s) && fZ>=%.6f && fZ<=%.6f", cut.Data(), lo - 1.0, hi + 1.0);

        out->cd();
        TTree *t = ce->CopyTree(sel);
        t->SetName(Form("Block%zu", b + 1));
        t->SetTitle(Form("CE in block %zu, z %.0f to %.0f mm", b + 1, lo, hi));
        
        // CE protons (daughters of the CE track) that reach the first tracker past this block
        auto next = std::find_if(planes.begin(), planes.end(),
                                 [&](const std::pair<long,long> &pl){ return pl.first > hi; });
        bool hasNext = next != planes.end();

        // per-CE branches added to the block tree
        Double_t nextTrackerZ = hasNext ? next->first : -1;
        Int_t reachedNextTracker = 0;
        Int_t reachedDetector = 0;   // proton that reached the next tracker also entered the detector (boundaryType 2)
        TBranch *bZ = t->Branch("nextTrackerZ", &nextTrackerZ, "nextTrackerZ/D");
        TBranch *bReached = t->Branch("reachedNextTracker", &reachedNextTracker, "reachedNextTracker/I");
        Double_t trackerX = 0, trackerY = 0;   // proton position at the next tracker plane
        TBranch *bTrkX = t->Branch("trackerX", &trackerX, "trackerX/D");
        TBranch *bTrkY = t->Branch("trackerY", &trackerY, "trackerY/D");
        Double_t detectorX = 0, detectorY = 0;   // proton position entering the detector
        TBranch *bDetX = t->Branch("detectorX", &detectorX, "detectorX/D");
        TBranch *bDetY = t->Branch("detectorY", &detectorY, "detectorY/D");
        TBranch *bDetector = t->Branch("reachedDetector", &reachedDetector, "reachedDetector/I");

        Long64_t nCE = t->GetEntries();
        Long64_t nReachedNext = 0, nReachedDet = 0;
        for (Long64_t i = 0; i < nCE; i++) {
            t->GetEntry(i);
            long ev = std::lround(t->GetLeaf("fEvent")->GetValue());
            long parent = std::lround(t->GetLeaf("TrackID")->GetValue());
            reachedNextTracker = 0;
            trackerX = trackerY = 0;
            if (hasNext) {
                auto hit = protonHits.lower_bound({ev, parent, next->first});
                reachedNextTracker = hit != protonHits.end() && std::get<0>(hit->first) == ev &&
                                     std::get<1>(hit->first) == parent && std::get<2>(hit->first) <= next->second;
                if (reachedNextTracker) { trackerX = hit->second.first; trackerY = hit->second.second; }
            }
            auto det = protonDetectorHits.find({ev, parent});
            reachedDetector = reachedNextTracker && det != protonDetectorHits.end();
            detectorX = reachedDetector ? det->second.first : 0;
            detectorY = reachedDetector ? det->second.second : 0;
            nReachedNext += reachedNextTracker;
            nReachedDet += reachedDetector;
            bZ->Fill();
            bReached->Fill();
            bTrkX->Fill();
            bTrkY->Fill();
            bDetector->Fill();
            bDetX->Fill();
            bDetY->Fill();
        }
        if (Verbose) {
            std::cout << "block " << b + 1 << " (z " << lo << " to " << hi << "): " << nCE << " CE in tree, next tracker z "
                      << nextTrackerZ << ", " << nReachedNext << " reached it, " << nReachedDet << " reached the ZDC\n";
        }
        t->Write();
    }
    out->Close();
}
