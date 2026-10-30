#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include "CutOpen.h"
#include "FitCore.h"

void CountCE(){
    const TString filename = "output0_t*.root";
    const TString outname = "cuts/cut_blocks.root";
    const TString cut = CarbonCeCut();
    const double gapMm = 300.0;
    const bool Verbose = false; // print the z bins, the density filter, the blocks found and per-block counts

    TChain *ce = new TChain("ChargeExchange_values");
    TChain *info = new TChain("ParticleInfo");
    ce->Add(filename);
    info->Add(filename);
    if (ce->GetEntries() == 0 || info->GetEntries() == 0) return;

    ROOT::EnableImplicitMT(NThreads());

    // every carbon-block CE: cuts/cut_CE.root (CarbonCeCut, CutOpen.h), so every row passes. Pass 1 needs no filter, the block copies below only the z window.
    TTree* ceCarbon = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
    if (!ceCarbon){Error("CountCE", "no cuts/cut_CE.root (run CutEvents.cpp)"); return;}

    // CE z positions and (event, track id) of every CE
    std::vector<double> z;
    std::set<std::pair<long,long>> ceKeys;
    {
        // per-thread-slot containers merged afterwards, so no global lock; z is sorted below
        ROOT::RDataFrame dfCE(*ceCarbon);
        std::vector<std::vector<double>> slotZ(dfCE.GetNSlots());
        std::vector<std::set<std::pair<long,long>>> slotKeys(dfCE.GetNSlots());
        dfCE.ForeachSlot([&](unsigned slot, Int_t ev, Int_t trk, Double_t fz){
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
    std::map<std::tuple<long,long,long>, std::array<double,4>> protonHits; // (event, parent, tracker z) -> (x, y, block x, block y)
    std::map<std::pair<long,long>, std::array<double,4>> protonDetectorHits;
    {
        // each thread fills its own containers (no lock per row), merged after the pass
        struct Local {
            std::set<long> trackerZ;
            std::map<std::tuple<long,long,long>, std::array<double,4>> protonHits; // (event, parent, tracker z) -> (x, y, block x, block y)
            std::map<std::pair<long,long>, std::array<double,4>> protonDetectorHits;
        };
        std::vector<Local> locals(std::max(1u, ROOT::GetThreadPoolSize()));
        ROOT::RDataFrame(*info).ForeachSlot([&](unsigned slot, Int_t ev, Int_t parent, Int_t pdg, Int_t bt,
                                                Double_t fx, Double_t fy, Double_t fz, Double_t fdx, Double_t fdy){
            Local &l = locals[slot];
            if (bt == 1) l.trackerZ.insert(std::lround(fz));
            if (pdg != 2212 || !ceKeys.count({ev, parent})) return;
            if (bt == 1) l.protonHits[{ev, parent, std::lround(fz)}] = {fx, fy, fdx, fdy};
            else if (bt == 2) l.protonDetectorHits[{ev, parent}] = {fx, fy, fdx, fdy};
        }, {"fEvent", "ParentID", "pdg", "boundaryType", "fX", "fY", "fZ", "fDx", "fDy"});

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

    //out file: every block tree is one multithreaded RDataFrame pass over cut_CE.root (selection, the per-CE columns below, Snapshot); the row order of a block tree is not the entry order
    gSystem->Unlink(outname);
    std::vector<std::string> ceBranches;
    for (auto *br : *ceCarbon->GetListOfBranches()) ceBranches.push_back(br->GetName());
    for (size_t b = 0; b < blocks.size(); b++) {
        double lo = z[blocks[b].first];
        double hi = z[blocks[b].second];
        TString sel = Form("fZ>=%.6f && fZ<=%.6f", lo - 1.0, hi + 1.0); // the carbon atoms of CarbonCeCut are already applied (cut_CE.root)

        // CE protons (daughters of the CE track) that reach the first tracker past this block
        auto next = std::find_if(planes.begin(), planes.end(),
                                 [&](const std::pair<long,long> &pl){ return pl.first > hi; });
        bool hasNext = next != planes.end();
        const Double_t nextTrackerZ = hasNext ? next->first : -1;
        const long nextLo = hasNext ? next->first : 0, nextHi = hasNext ? next->second : 0;

        // per CE: {reachedNextTracker, trackerX, trackerY, trackerDX, trackerDY, reachedDetector, detectorX, detectorY, detectorDX, detectorDY}
        auto perCE = [&](Int_t evI, Int_t parentI){
            const long ev = evI, parent = parentI;
            std::array<double,10> r = {0,0,0,0,0,0,0,0,0,0};
            if (hasNext) {
                auto hit = protonHits.lower_bound({ev, parent, nextLo});
                if (hit != protonHits.end() && std::get<0>(hit->first) == ev && std::get<1>(hit->first) == parent && std::get<2>(hit->first) <= nextHi) {
                    r[0] = 1; r[1] = hit->second[0]; r[2] = hit->second[1]; r[3] = hit->second[2]; r[4] = hit->second[3];
                }
            }
            auto det = protonDetectorHits.find({ev, parent});
            if (r[0] && det != protonDetectorHits.end()) {
                r[5] = 1; r[6] = det->second[0]; r[7] = det->second[1]; r[8] = det->second[2]; r[9] = det->second[3];
            }
            return r;
        };
        ROOT::RDF::RNode df = ROOT::RDataFrame(*ceCarbon).Filter(sel.Data(), "block z window");
        df = df.Define("ceRow", perCE, {"fEvent", "TrackID"})
               .Define("nextTrackerZ", [nextTrackerZ]{return nextTrackerZ;})
               .Define("reachedNextTracker", [](const std::array<double,10> &r){return (Int_t)r[0];}, {"ceRow"})
               .Define("trackerX", [](const std::array<double,10> &r){return r[1];}, {"ceRow"})
               .Define("trackerY", [](const std::array<double,10> &r){return r[2];}, {"ceRow"})
               .Define("trackerDX", [](const std::array<double,10> &r){return r[3];}, {"ceRow"})
               .Define("trackerDY", [](const std::array<double,10> &r){return r[4];}, {"ceRow"})
               .Define("detectorDX", [](const std::array<double,10> &r){return r[8];}, {"ceRow"})
               .Define("detectorDY", [](const std::array<double,10> &r){return r[9];}, {"ceRow"})
               .Define("detectorX", [](const std::array<double,10> &r){return r[6];}, {"ceRow"})
               .Define("detectorY", [](const std::array<double,10> &r){return r[7];}, {"ceRow"})
               .Define("reachedDetector", [](const std::array<double,10> &r){return (Int_t)r[5];}, {"ceRow"});
        auto nCEr = df.Count();
        auto nReachedNextr = df.Sum<Int_t>("reachedNextTracker");
        auto nReachedDetr = df.Sum<Int_t>("reachedDetector");
        std::vector<std::string> cols = ceBranches;
        for (const char *c : {"nextTrackerZ", "reachedNextTracker", "trackerX", "trackerY", "trackerDX", "trackerDY", "detectorDX", "detectorDY", "detectorX", "detectorY", "reachedDetector"}) cols.push_back(c);
        ROOT::RDF::RSnapshotOptions opts;
        opts.fMode = b == 0 ? "RECREATE" : "UPDATE";
        df.Snapshot(Form("Block%zu", b + 1), outname.Data(), cols, opts);
        Long64_t nCE = *nCEr, nReachedNext = (Long64_t)*nReachedNextr, nReachedDet = (Long64_t)*nReachedDetr;
        {
            TFile f(outname, "UPDATE");
            if (auto *t = (TTree*)f.Get(Form("Block%zu", b + 1))) {
                t->SetTitle(Form("CE in block %zu, z %.0f to %.0f mm", b + 1, lo, hi));
                t->Write("", TObject::kOverwrite);
            }
        }
        if (Verbose) {
            std::cout << "block " << b + 1 << " (z " << lo << " to " << hi << "): " << nCE << " CE in tree, next tracker z "
                      << nextTrackerZ << ", " << nReachedNext << " reached it, " << nReachedDet << " reached the ZDC\n";
        }
    }
}
