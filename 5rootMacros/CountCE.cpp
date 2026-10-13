#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <vector>

void CountCE(){
    const TString filename = "output0_t*.root";
    const TString outname = "cut_blocks.root";
    const TString cut = "didChargeExchange>0 && fZ<7000";
    const double gapMm = 300.0;

    TChain *ce = new TChain("ChargeExchange_values");
    TChain *info = new TChain("ParticleInfo");
    ce->Add(filename);
    info->Add(filename);
    if (ce->GetEntries() == 0 || info->GetEntries() == 0) return;

    // CE z positions and (event, track id) of every CE
    std::vector<double> z;
    std::set<std::pair<long,long>> ceKeys;
    {
        TTreeReader reader(ce);
        TTreeReaderValue<Int_t> ev(reader, "fEvent"), trk(reader, "TrackID");
        TTreeReaderValue<Double_t> did(reader, "didChargeExchange"), fz(reader, "fZ");
        while (reader.Next()) {
            if (*did > 0 && *fz < 7000) {
                z.push_back(*fz);
                ceKeys.insert({*ev, *trk});
            }
        }
    }
    if (z.empty()) return;
    std::sort(z.begin(), z.end());

    // One pass over ParticleInfo, nothing big held in memory:
    //  - tracker z (boundaryType 1) of every particle, to find the tracker planes
    //  - tracker hits and detector entries (boundaryType 2) of CE protons only
    //    (event, parent track id) -> tracker z / detector (x, y)
    std::set<long> trackerZ;
    std::set<std::tuple<long,long,long>> protonHits;
    std::map<std::pair<long,long>, std::pair<double,double>> protonDetectorHits;
    {
        TTreeReader reader(info);
        TTreeReaderValue<Int_t> ev(reader, "fEvent"), parent(reader, "ParentID"),
                                pdg(reader, "pdg"), bt(reader, "boundaryType");
        TTreeReaderValue<Double_t> fx(reader, "fX"), fy(reader, "fY"), fz(reader, "fZ");
        while (reader.Next()) {
            if (*bt == 1) trackerZ.insert(std::lround(*fz));
            if (*pdg != 2212 || !ceKeys.count({*ev, *parent})) continue;
            if (*bt == 1) protonHits.insert({*ev, *parent, std::lround(*fz)});
            else if (*bt == 2) protonDetectorHits[{*ev, *parent}] = {*fx, *fy};
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
        Double_t detectorX = 0, detectorY = 0;   // proton position entering the detector
        TBranch *bDetX = t->Branch("detectorX", &detectorX, "detectorX/D");
        TBranch *bDetY = t->Branch("detectorY", &detectorY, "detectorY/D");
        TBranch *bDetector = t->Branch("reachedDetector", &reachedDetector, "reachedDetector/I");

        Long64_t nCE = t->GetEntries();
        for (Long64_t i = 0; i < nCE; i++) {
            t->GetEntry(i);
            long ev = std::lround(t->GetLeaf("fEvent")->GetValue());
            long parent = std::lround(t->GetLeaf("TrackID")->GetValue());
            reachedNextTracker = 0;
            if (hasNext) {
                auto hit = protonHits.lower_bound({ev, parent, next->first});
                reachedNextTracker = hit != protonHits.end() && std::get<0>(*hit) == ev &&
                                     std::get<1>(*hit) == parent && std::get<2>(*hit) <= next->second;
            }
            auto det = protonDetectorHits.find({ev, parent});
            reachedDetector = reachedNextTracker && det != protonDetectorHits.end();
            detectorX = reachedDetector ? det->second.first : 0;
            detectorY = reachedDetector ? det->second.second : 0;
            bZ->Fill();
            bReached->Fill();
            bDetector->Fill();
            bDetX->Fill();
            bDetY->Fill();
        }
        t->Write();
    }
    out->Close();
}
