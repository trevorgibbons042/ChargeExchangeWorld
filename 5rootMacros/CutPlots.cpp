#include <ROOT/RDFHelpers.hxx> // RunGraphs, RResultHandle
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <string>
#include <vector>
#include "FitCore.h"
#include "Plot3DCols.h"
#include <TChain.h>
#include <TString.h>
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
#include <TROOT.h>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include "CutFiles.h"

// Every cut file the plot macros read, in one macro, as a part (CutPlots("<part>")):
//   events   cuts/cut_*.root from the sim output (CutEvents), run_all.sh runs it after ReconAngleAll
//   count    the CE counts per block and tracker plane (CountCE), after events
//   zdc      cuts/cut_Plot3DCE.root, cut_Plot3DAllEv.root and cut_Plot3DAllPart.root, the rows of AnglesRaw by evClass with the Plot3D columns (ZdcColumns, Plot3DCols.h)

// Writes the cut files F1_2_Plot3D.cpp draws from: cut_Plot3DCE.root (AnglesRaw rows with evClass 0), cut_Plot3DAllEv.root (evClass 0 or 1) and cut_Plot3DAllPart.root (all rows),
// each row with the Plot3D columns of its recon fit (isReco, raw_0..26, ok_0..26, Plot3DCols.h).
// Run after ReconAngleAll.cpp, with all the threads: the maps of AnglesRaw are read once (serial), then the files are written in one parallel pass.
static void ZdcColumns(){
    ROOT::EnableImplicitMT(NThreads());
    ROOT::RDF::RSnapshotOptions lazy;
    lazy.fLazy = true;

    const time_t t0 = time(nullptr);
    const ColumnAdder addColumns = Plot3DColumns(); // loads the maps on its first call
    // the three sets of F1_2_Plot3D.cpp are rows of AnglesRaw (one recon fit per track the finder finds) split by evClass: 0 the fit used CE child hits only (CE recon protons), 0 or 1 CE child and other hits
    // (CE + CE mix), any (all recon). The position of a row is the last tracker hit (lastX, lastY, lastXD, lastYD), the KE the ZDC sum of its event (finalKE, finalKED).
    const char* angFile = "cuts/cut_angles_all.root";
    ROOT::RDataFrame dfCE(kAngleTree, angFile), dfAllEv(kAngleTree, angFile), dfAllPart(kAngleTree, angFile);
    std::vector<ROOT::RDF::RResultHandle> handles;
    handles.emplace_back(addColumns(dfCE.Filter([](int c){return c == 0;}, {"evClass"})).Snapshot(kAngleTree, "cuts/cut_Plot3DCE.root", "", lazy));
    const time_t t1 = time(nullptr); // the first addColumns call has loaded the maps
    printf(">>> CutPlots zdc maps took %.0f s\n", difftime(t1, t0)); fflush(stdout);
    handles.emplace_back(addColumns(dfAllEv.Filter([](int c){return c == 0 || c == 1;}, {"evClass"})).Snapshot(kAngleTree, "cuts/cut_Plot3DAllEv.root", "", lazy));
    handles.emplace_back(addColumns(dfAllPart).Snapshot(kAngleTree, "cuts/cut_Plot3DAllPart.root", "", lazy));
    ROOT::RDF::RunGraphs(handles);
    printf(">>> CutPlots zdc write took %.0f s\n", difftime(time(nullptr), t1)); fflush(stdout);

}

// ===== CutEvents (part "events", run after ReconAngleAll) =====
void CutEvents(){
    ROOT::EnableImplicitMT(NThreads());
    
    TChain *events = new TChain("ParticleInfo");
    TChain *ChargeExchange = new TChain("ChargeExchange_values");
    
    TString filename = "output0_t*.root";
    events->Add(filename);
    ChargeExchange->Add(filename);
    
    // (fEvent, TrackID) of every carbon-block CE, read on all threads as one key so the pairs stay together.
    // TChain::Draw + GetV1 kept only part of them (SetEstimate(-1) on a wildcard TChain): a cut file came out with 1 row.
    ROOT::RDataFrame dfCESel(*ChargeExchange);
    auto keys = dfCESel.Filter(CarbonCeCut())
        .Define("ceKey", [](int ev, int id){return ((long long)ev << 32) | (unsigned int)id;}, {"fEvent", "TrackID"})
        .Take<long long>("ceKey");

    // Lookup tables instead of one huge "a || b || ..." string per cut: those strings were JIT-compiled by cling
    // (serial, slow) and tested term by term for every entry.
    std::unordered_map<int, std::vector<int>> ceTracksInEvent;
    std::unordered_set<long long> ceEventParent;
    for (long long k : *keys){
        ceEventParent.insert(k);
        ceTracksInEvent[int(k >> 32)].push_back(int(k & 0xffffffff));
    }
    printf(">>> CutEvents: %zu carbon CE in %zu events\n", keys->size(), ceTracksInEvent.size());

    // event has a charge exchange
    auto inEvent = [&](int ev){return ceTracksInEvent.count(ev) > 0;};
    // track is a child of a charge exchange track in this event
    auto isParent = [&](int ev, int parent){return ceEventParent.count(((long long)ev << 32) | (unsigned int)parent) > 0;};

    ROOT::RDataFrame df(*events);
    ROOT::RDataFrame dfCEAll(*ChargeExchange);

    // The cuts chain off each other (see CutFiles.h): the big sim output is read once, in stage 1, for the two base files; the other cut is a filter of the small cut_CE.
    //   stage 1 (sim output):  cut_CE (events with a carbon CE), cut_CEall (every CE)
    //   stage 2 (cut_CE):      cut_ParentTracker (the protons that are children of a carbon CE track, at the first tracker)
    //   stage 3 (cut_CEall):   the carbon CE rows, added as the second tree (UPDATE) of cut_CE
    // Each stage is lazy: all the cuts of a stage are written in one pass over their input. The result handles own the open output files, so they are cleared (files closed)
    // before the next stage reads those files; left open, cut_CE.root came out without ChargeExchange_values.
    ROOT::RDF::RSnapshotOptions lazy;
    lazy.fLazy = true;
    const std::string pInfo = events->GetName(), ceTree = ChargeExchange->GetName();

    {
        std::vector<ROOT::RDF::RResultHandle> handles;
        handles.emplace_back(df.Filter(inEvent, {"fEvent"}).Snapshot(pInfo, "cuts/cut_CE.root", "", lazy));
        handles.emplace_back(dfCEAll.Snapshot(ceTree, "cuts/cut_CEall.root", "", lazy)); // every CE of the sim (the whole ChargeExchange_values tree, one row per CE)
        ROOT::RDF::RunGraphs(handles);
    }
    {
        std::vector<ROOT::RDF::RResultHandle> handles;
        ROOT::RDataFrame dfCEfile(pInfo, "cuts/cut_CE.root"); // stays alive until RunGraphs
        auto children = dfCEfile.Filter(isParent, {"fEvent", "ParentID"});
        handles.emplace_back(children.Filter("boundaryType == 0 && pdg == 2212").Snapshot(pInfo, "cuts/cut_ParentTracker.root", "", lazy));
        ROOT::RDF::RunGraphs(handles);
    }
    ROOT::RDF::RSnapshotOptions update;
    update.fMode = "UPDATE";
    ROOT::RDataFrame dfCEallRows(ceTree, "cuts/cut_CEall.root");
    dfCEallRows.Filter(CarbonCeCut()).Snapshot(ceTree, "cuts/cut_CE.root", "", update); // the carbon-block CE rows, beside ParticleInfo
}

// ===== CountCE (part "count", run after events) =====
void CountCE(){
    const TString outname = "cuts/cut_blocks.root";

    ROOT::EnableImplicitMT(NThreads());

    // every carbon-block CE: cuts/cut_CE.root (CarbonCeCut, CutFiles.h), so every row passes. Pass 1 needs no filter, the block copies below only the z window.
    TTree* ceCarbon = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
    if (!ceCarbon){Error("CountCE", "no cuts/cut_CE.root (run CutPlots.cpp (events))"); return;}
    // the CE protons (children of a carbon CE track) are rows of the events that have a carbon CE: ParticleInfo of cut_CE.root, not the sim output
    TTree* info = OpenCut("cuts/cut_CE.root", "ParticleInfo");
    if (!info){Error("CountCE", "no ParticleInfo in cuts/cut_CE.root (run CutPlots.cpp (events))"); return;}
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

    // The block and tracker z come from the geometry in settings.csv (2construction.cc): carbon blocks evenly spread from z = 0 to the ZDC face, a tracker plane either side of each block (none in front of the first) and one in front of the ZDC, each a gap dz from the block or ZDC edge
    // (an old run: nBlocksTracking planes per block, dz apart around z = 4 m). A carbon block is CarbonThickness thick, a tracker plane detectorThickness (500 um).
    const auto st = ReadSettings();
    const double carbonMm = st.at("CarbonThickness") * 1000.0, detectorMm = 0.5;
    const StackLayout stack = StackZ(st);
    std::vector<std::pair<double,double>> blocks; // [lowest z, highest z] (mm) of each carbon block
    std::vector<std::pair<long,long>> planes;     // [lowest z, highest z] (mm) of each tracker plane
    for (const double zc : stack.planeZ) planes.push_back({std::lround(zc - detectorMm / 2), std::lround(zc + detectorMm / 2)});
    for (const double zc : stack.blockZ) blocks.push_back({zc - carbonMm / 2, zc + carbonMm / 2});

    // One pass over ParticleInfo, nothing big held in memory:
    //  - tracker hits (entry rows, boundaryType 0) and detector entries (boundaryType 2) of CE protons only
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
            if (bt == 0) l.protonHits[{ev, parent, std::lround(fz)}] = {fx, fy, fdx, fdy};
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
               .Define("reachedNextTracker", [](const std::array<double,10> &r){return (Int_t)r[0];}, {"ceRow"})
               .Define("reachedDetector", [](const std::array<double,10> &r){return (Int_t)r[5];}, {"ceRow"});
        std::vector<std::string> cols = ceBranches;
        cols.push_back("reachedNextTracker");
        // the double columns of ceRow, one Define each
        const auto addCol = [&](const char *name, int idx){
            df = df.Define(name, [idx](const std::array<double,10> &r){return r[idx];}, {"ceRow"});
            cols.push_back(name);
        };
        addCol("trackerX", 1); addCol("trackerY", 2); addCol("trackerDX", 3); addCol("trackerDY", 4);
        addCol("detectorDX", 8); addCol("detectorDY", 9); addCol("detectorX", 6); addCol("detectorY", 7);
        cols.push_back("reachedDetector");
        ROOT::RDF::RSnapshotOptions opts;
        opts.fMode = b == 0 ? "RECREATE" : "UPDATE";
        df.Snapshot(Form("Block%zu", b + 1), outname.Data(), cols, opts);
        {
            TFile f(outname, "UPDATE");
            if (auto *t = (TTree*)f.Get(Form("Block%zu", b + 1))) {
                t->SetTitle(Form("CE in block %zu, z %.0f to %.0f mm", b + 1, lo, hi));
                t->Write("", TObject::kOverwrite);
            }
        }
    }
}

void CutPlots(const char* partArg = "zdc"){ // a const char*: calling a macro with a std::string argument fails on the Mac cling (Failed to materialize symbols)
    const std::string part = partArg;
    if (part == "events"){CutEvents(); return;} // the cut files of the sim output, after ReconAngleAll in run_all.sh
    if (part == "count"){CountCE(); return;} // the CE counts per block, after CutEvents
    if (part != "zdc"){Error("CutPlots", "part must be events, count or zdc"); return;}
    const time_t t0 = time(nullptr);
    ZdcColumns();
    printf(">>> CutPlots zdc took %.0f s\n", difftime(time(nullptr), t0)); fflush(stdout);
}
