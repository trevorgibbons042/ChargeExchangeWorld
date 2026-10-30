#include <TFile.h>
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
    const TString outname = "cuts/cut_blocks.root";
    const TString cut = CarbonCeCut();
    const bool Verbose = false; // print the blocks and tracker planes and the per-block counts


    ROOT::EnableImplicitMT(NThreads());

    // every carbon-block CE: cuts/cut_CE.root (CarbonCeCut, CutOpen.h), so every row passes. Pass 1 needs no filter, the block copies below only the z window.
    TTree* ceCarbon = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
    if (!ceCarbon){Error("CountCE", "no cuts/cut_CE.root (run CutEvents.cpp)"); return;}
    // the CE protons (children of a carbon CE track) are rows of the events that have a carbon CE: ParticleInfo of cut_CE.root, not the sim output
    TTree* info = OpenCut("cuts/cut_CE.root", "ParticleInfo");
    if (!info){Error("CountCE", "no ParticleInfo in cuts/cut_CE.root (run CutEvents.cpp)"); return;}
    if (ceCarbon->GetEntries() == 0 || info->GetEntries() == 0) return;

    // (event, track id) of every carbon-block CE
    std::set<std::pair<long,long>> ceKeys;
    {
        // per-thread-slot containers merged afterwards, so no global lock
        ROOT::RDataFrame dfCE(*ceCarbon);
        std::vector<std::set<std::pair<long,long>>> slotKeys(dfCE.GetNSlots());
        dfCE.ForeachSlot([&](unsigned slot, Int_t ev, Int_t trk){slotKeys[slot].insert({ev, trk});}, {"fEvent", "TrackID"});
        for (auto &k : slotKeys) ceKeys.insert(k.begin(), k.end());
    }
    if (ceKeys.empty()) return;

    // The block and tracker z come from the geometry in settings.csv (2construction.cc): nBlocksTracking tracker planes before each carbon block and after the last one,
    // all elements dz apart around z = 4 m. A carbon block is CarbonThickness thick, a tracker plane detectorThickness (500 um).
    const auto st = ReadSettings();
    const int nBlocks = (int)st.at("nBlocks"), nTrk = (int)st.at("nBlocksTracking");
    const double dzMm = st.at("dz") * 1000.0, carbonMm = st.at("CarbonThickness") * 1000.0, detectorMm = 0.5;
    const int nElements = (nBlocks + 1) * nTrk + nBlocks;
    const double zStartMm = 4000.0 - (nElements - 1) * dzMm / 2.0;
    std::vector<std::pair<double,double>> blocks; // [lowest z, highest z] (mm) of each carbon block
    std::vector<std::pair<long,long>> planes;     // [lowest z, highest z] (mm) of each tracker plane
    for (int i = 0, element = 0; i <= nBlocks; i++){
        for (int j = 0; j < nTrk; j++, element++){const double zc = zStartMm + element * dzMm; planes.push_back({std::lround(zc - detectorMm / 2), std::lround(zc + detectorMm / 2)});}
        if (i == nBlocks) break;
        const double zc = zStartMm + element * dzMm; element++;
        blocks.push_back({zc - carbonMm / 2, zc + carbonMm / 2});
    }
    if (Verbose) {
        std::cout << blocks.size() << " blocks, " << planes.size() << " tracker planes from the geometry:\n";
        for (auto &b : blocks) std::cout << "  block z " << b.first << " to " << b.second << "\n";
        for (auto &pl : planes) std::cout << "  plane z [" << pl.first << "," << pl.second << "]\n";
    }

    // One pass over ParticleInfo, nothing big held in memory:
    //  - tracker hits and detector entries (boundaryType 2) of CE protons only
    //    (event, parent track id) -> tracker z / detector (x, y)
    std::map<std::tuple<long,long,long>, std::array<double,4>> protonHits; // (event, parent, tracker z) -> (x, y, block x, block y)
    std::map<std::pair<long,long>, std::array<double,4>> protonDetectorHits;
    {
        // each thread fills its own containers (no lock per row), merged after the pass
        struct Local {
            std::map<std::tuple<long,long,long>, std::array<double,4>> protonHits; // (event, parent, tracker z) -> (x, y, block x, block y)
            std::map<std::pair<long,long>, std::array<double,4>> protonDetectorHits;
        };
        std::vector<Local> locals(std::max(1u, ROOT::GetThreadPoolSize()));
        ROOT::RDataFrame(*info).ForeachSlot([&](unsigned slot, Int_t ev, Int_t parent, Int_t pdg, Int_t bt,
                                                Double_t fx, Double_t fy, Double_t fz, Double_t fdx, Double_t fdy){
            Local &l = locals[slot];
            if (pdg != 2212 || !ceKeys.count({ev, parent})) return;
            if (bt == 1) l.protonHits[{ev, parent, std::lround(fz)}] = {fx, fy, fdx, fdy};
            else if (bt == 2) l.protonDetectorHits[{ev, parent}] = {fx, fy, fdx, fdy};
        }, {"fEvent", "ParentID", "pdg", "boundaryType", "fX", "fY", "fZ", "fDx", "fDy"});

        for (auto &l : locals) {
            protonHits.insert(l.protonHits.begin(), l.protonHits.end());
            for (auto &kv : l.protonDetectorHits) protonDetectorHits[kv.first] = kv.second;
        }
    }

    //out file: every block tree is one multithreaded RDataFrame pass over cut_CE.root (selection, the per-CE columns below, Snapshot); the row order of a block tree is not the entry order
    gSystem->Unlink(outname);
    std::vector<std::string> ceBranches;
    for (auto *br : *ceCarbon->GetListOfBranches()) ceBranches.push_back(br->GetName());
    for (size_t b = 0; b < blocks.size(); b++) {
        const double lo = blocks[b].first, hi = blocks[b].second;
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
