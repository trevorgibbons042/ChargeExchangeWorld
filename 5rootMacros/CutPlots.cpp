#include <ROOT/RDFHelpers.hxx> // RunGraphs, RResultHandle
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>
#include <TEntryList.h>
#include <TFile.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "CutFiles.h"
#include "FitCore.h"
#include "Plot3DCols.h"
#include "QualityCut.h"

// Every cut file the plot macros read besides the ones CutEvents.cpp and ReconAngleAll.cpp write, in one macro (the order matters: the samples read what the zdc part writes):
//   trees    cuts/ke_gt80.root, the KE cut entry lists of cut_APandT.root and cut_angles_all.root (CutTrees, CutFiles.h)
//   zdc      cuts/cut_ZdcCE.root and cut_ZdcAllEv.root, the ZDC protons with the Plot3D columns (ZdcColumns, Plot3DCols.h)
//   samples  one cut file per cut of 6_cuts and cut_hitZDC_AllEv.root (CutSamples)
// CutPlots() runs all three; CutPlots("samples") only the last one (plot_cut.sh). Only the samples that have a ZDC KE cut (megacut2, megacut3, the KE scan, hitZDC) read the zdc files, so Plot.pbs runs trees, zdc and the other samples
// side by side, then the ZDC KE samples as several processes (CUT_ONLY).

// Applies every shared cut (KeCuts, CutFiles.h) once to the already cut trees of cut_APandT.root (CalcAPandT.cpp) and to the proton tree AnglesRaw of cut_angles_all.root
// (ReconAngleAll.cpp), true and detector energy, and writes the entries that pass as TEntryLists to cuts/<ke dir>.root (no residual cut): one list per tree,
// named like the tree. Run after CalcAPandT.cpp.
// AnglesRaw gets one list per view (Angles, AnglesCE: the rows AngleSelection picks) and energy, named view + suffix, with the KE cut on that energy's KE branch.
static void CutTrees(){
    TFile* fAP = TFile::Open("cuts/cut_APandT.root", "READ");
    TFile* fAng = TFile::Open("cuts/cut_angles_all.root", "READ");
    if (!fAP || !fAng){Error("CutTrees", "cuts/cut_APandT.root or cut_angles_all.root missing"); return;}
    // One list to build: the tree it indexes, its name, the rows it covers (sel) and the energy suffix of its KE cut.
    struct Job {TTree* t; std::string name, base, sel, sfx;};
    std::vector<Job> jobs;
    TTree* raw = (TTree*)fAng->Get(kAngleTree);
    for (const std::string sfx : {"", "_DE"}){
        for (const auto& n : ApTrees()) jobs.push_back({(TTree*)fAP->Get((n + sfx).c_str()), n + sfx, n, "", sfx});
        for (const auto& n : AngleTrees()) jobs.push_back({raw, n + sfx, n, AngleSelection(n), sfx});
    }
    for (const auto& j : jobs) if (!j.t){Error("CutTrees", "a tree is missing (rerun CalcAPandT.cpp / ReconAngleAll.cpp)"); return;}

    gSystem->mkdir("cuts", true);
    ROOT::EnableImplicitMT(NThreads());
    for (const auto& ke : KeCuts()){
        // One event loop per tree on all the threads books every list of that tree at once: the entries that pass are collected per thread (with their entry number, rdfentry_),
        // sorted, and entered into the TEntryLists afterwards.
        struct Want {std::string name, title; int view; bool detector;}; // view: 0 none, 1 Angles (every row), 2 AnglesCE (inCE)
        struct Hit {TEntryList* list; std::vector<std::vector<ULong64_t>> perSlot;};
        TFile out(("cuts/" + ke.dir + ".root").c_str(), "RECREATE");
        const unsigned nSlots = std::max(1u, NThreads()) + 1;
        auto finish = [&](TTree* t, std::vector<Hit>& hits){
            out.cd();
            for (auto& h : hits){
                std::vector<ULong64_t> all;
                for (auto& v : h.perSlot) all.insert(all.end(), v.begin(), v.end());
                std::sort(all.begin(), all.end());
                for (ULong64_t e : all) h.list->Enter(e);
                h.list->SetDirectory(&out);
            }
        };
        // cut_APandT.root: one tree per name and energy, every list on finalKE (the _DE trees hold the detector KE there)
        for (const std::string sfx : {"", "_DE"}) for (const auto& n : ApTrees()){
            TTree* t = (TTree*)fAP->Get((n + sfx).c_str());
            std::vector<Hit> hits(1);
            out.cd();
            hits[0] = {new TEntryList((n + sfx).c_str(), ke.cut.c_str(), t), std::vector<std::vector<ULong64_t>>(nSlots)};
            ROOT::RDataFrame df0(n + sfx, "cuts/cut_APandT.root");
            auto df = WithEntryNumber(df0, df0.GetNSlots()); // the true entry numbers (rdfentry_ is wrong under MT on a tree with several clusters)
            const double minKE = ke.minKE;
            df.ForeachSlot([&](unsigned slot, ULong64_t e, double keVal){
                if (!(keVal > minKE)) return;
                hits[0].perSlot[slot].push_back(e);
            }, {"entryNo", "finalKE"});
            finish(t, hits);
        }
        // AnglesRaw: the views Angles and AnglesCE of both energies in one loop
        {
            std::vector<Want> wants;
            for (const std::string sfx : {"", "_DE"}) for (const auto& v : AngleTrees())
                wants.push_back({v + sfx, "", v == "AnglesCE" ? 2 : 1, !sfx.empty()});
            std::vector<Hit> hits(wants.size());
            out.cd();
            for (size_t i = 0; i < wants.size(); i++){
                const std::string cut = AngleSelection(wants[i].view == 2 ? "AnglesCE" : "Angles") + " && " + KeCutFor(ke.cut, wants[i].detector ? "_DE" : "");
                hits[i] = {new TEntryList(wants[i].name.c_str(), cut.c_str(), raw), std::vector<std::vector<ULong64_t>>(nSlots)};
            }
            ROOT::RDataFrame df0(kAngleTree, "cuts/cut_angles_all.root");
            auto df = WithEntryNumber(df0, df0.GetNSlots());
            const double minKE = ke.minKE;
            df.ForeachSlot([&](unsigned slot, ULong64_t e, int inCE, double keTrue, double keDet){
                for (size_t i = 0; i < wants.size(); i++){
                    const Want& w = wants[i];
                    if (w.view == 2 && !inCE) continue;
                    if (!((w.detector ? keDet : keTrue) > minKE)) continue;
                    hits[i].perSlot[slot].push_back(e);
                }
            }, {"entryNo", "inCE", "finalKE", "finalKED"});
            finish(raw, hits);
        }
        out.Write();
        printf("[CutTrees] cuts/%s.root\n", ke.dir.c_str());
    }
}

// Writes the two ZDC proton cut files Plot3D.cpp draws from: cut_ZdcCE.root (the ZDC protons of the CE children, cut_ParentZdc.root) and cut_ZdcAllEv.root (the ZDC protons
// of every event, cut_AllEvProtons.root), each row with the Plot3D columns of its reconstructed track (isReco, raw_0..16, ok_0..16, Plot3DCols.h).
// Run after CutEvents.cpp and ReconAngleAll.cpp, with all the threads: the maps of AnglesRaw are read once (serial), then both files are written in one parallel pass.
static void ZdcColumns(){
    ROOT::EnableImplicitMT(NThreads());
    ROOT::RDF::RSnapshotOptions lazy;
    lazy.fLazy = true;

    const time_t t0 = time(nullptr);
    const ColumnAdder addColumns = Plot3DColumns(); // loads the maps on its first call
    ROOT::RDataFrame dfCE("ParticleInfo", "cuts/cut_ParentZdc.root"), dfAllEv("ParticleInfo", "cuts/cut_AllEvProtons.root"), dfAllPart("ParticleInfo", "cuts/cut_AllEvParticles.root");
    std::vector<ROOT::RDF::RResultHandle> handles;
    handles.emplace_back(addColumns(dfCE).Snapshot("ParticleInfo", "cuts/cut_ZdcCE.root", "", lazy));
    const time_t t1 = time(nullptr); // the first addColumns call has loaded the maps
    printf(">>> CutPlots zdc maps took %.0f s\n", difftime(t1, t0)); fflush(stdout);
    handles.emplace_back(addColumns(dfAllEv).Snapshot("ParticleInfo", "cuts/cut_ZdcAllEv.root", "", lazy));
    handles.emplace_back(addColumns(dfAllPart).Snapshot("ParticleInfo", "cuts/cut_ZdcAllPart.root", "", lazy)); // the same for every recorded particle, any pdg (the "all recon particles" plots)
    ROOT::RDF::RunGraphs(handles);
    printf(">>> CutPlots zdc write took %.0f s\n", difftime(time(nullptr), t1)); fflush(stdout);
}

// Writes one cut file per cut of 6_cuts (Plot3D.cpp), each holding only the rows that pass, so Plot3D reads the file instead of cutting all of AnglesRaw with a Draw string.
// AnglesRaw cuts (every fitted proton, no ZDC requirement, cut_angles_all.root) keep the branches Plot3D reads:
//   cut_ResCuts       resDRms < 4 mm (3_ResCuts)
//   cut_nParticles1/2 at most 1 / 2 reconstructed particles in the event
//   cut_nZdc1         at most 1 charged particle at the ZDC in the event (8_nZdcParticles, EventSummary of the sim)
//   cut_nHits         2, 4, 6 or 8 tracker hits
//   cut_megacut_1     6_megaCut/1_1particlecut: at most 1 reconstructed particle in the event, 2, 4, 6 or 8 tracker hits, recon theta and recon block theta 0 to 1 degree
//   cut_megacut2_1    7_megaCut2/1_1particlecut: the same, and a reconstructed ZDC proton with ZDC KE > 80% of the beam KE
//   cut_megacut3_1    9_megaCut3/1_1particlecut: the cut_megacut2 cuts, 2, 4, 6 or 8 charged tracker hits in the event (nTrackerHitsEv, all particles added up) and exactly 1 charged particle at the ZDC (nZdcCharged)
// ZDC cut (cut_ZdcAllEv.root, the ZDC protons with the Plot3D columns of ZdcColumns.cpp):
//   cut_hitZDC_AllEv           the proton has a reconstructed track (isReco)
// Run after ZdcColumns.cpp, with all the threads: every input is read once and all the cut files of it are written in the same pass.
// CUT_ONLY (environment, comma separated, e.g. megacut3 or megacut3,ResCuts): CutSamples writes only the cut files whose name holds one of these (cut_megacut3_1.root, cut_megacut3_1_ke85.root, ...),
// and CutPlots("all") then runs only the samples part (the trees and zdc files stay as they are). Empty: every file.
static bool CutWanted(const std::string& file){
    const char* e = std::getenv("CUT_ONLY");
    if (!e || !*e) return true;
    std::stringstream ss(e);
    for (std::string w; std::getline(ss, w, ',');) if (!w.empty() && file.find(w) != std::string::npos) return true;
    return false;
}
static void CutSamples(){
    ROOT::EnableImplicitMT(NThreads());
    ROOT::RDF::RSnapshotOptions lazy;
    lazy.fLazy = true;
    std::vector<ROOT::RDF::RResultHandle> handles;

    const std::vector<std::string> rawCols = {"fEvent", "TrackID", "PDG", "isFromCE", "inCE", "finalKE", "resRms", "resDRms", "impact", "impactD", "r2", "r2D", "nProtonsEv", "nZdcCharged", "nTrackerHitsEv", "nTrackerTracksEv",
                                              "resHalfRatioD", "outTheta", "outDTheta", "outPhi", "outDPhi", "nHits", "firstXD", "firstYD", "z_bx_reconVertexOut", "ParentID"};
    ROOT::RDataFrame raw("AnglesRaw", "cuts/cut_angles_all.root");
    // baseline of every cut of 6_cuts: recon theta and recon block theta from 0 to 1 degree (NaN fails)
    auto base = raw.Filter([](double theta, double thetaBlock){const double deg = theta * TMath::RadToDeg(), degBlock = thetaBlock * TMath::RadToDeg(); return deg >= 0 && deg <= 1 && degBlock >= 0 && degBlock <= 1;}, {"outTheta", "outDTheta"});
    auto snap = [&](ROOT::RDF::RNode node, const char* file){if (CutWanted(file)) handles.emplace_back(node.Snapshot("AnglesRaw", file, rawCols, lazy));};
    snap(base.Filter([](double resD){return resD < kRes4Mm;}, {"resDRms"}), "cuts/cut_ResCuts.root");
    snap(base.Filter([](int n){return n <= 1;}, {"nProtonsEv"}), "cuts/cut_nParticles1.root");
    snap(base.Filter([](int n){return n <= 2;}, {"nProtonsEv"}), "cuts/cut_nParticles2.root");
    snap(base.Filter([](int n){return n >= 0 && n <= 1;}, {"nZdcCharged"}), "cuts/cut_nZdc1.root"); // -1: no EventSummary, fails
    snap(base.Filter([](int n){return n == 2 || n == 4 || n == 6 || n == 8;}, {"nHits"}), "cuts/cut_nHits.root");

    // 7_megaCut2: the same cuts as 6_megaCut, and the proton must also be a reconstructed ZDC proton with ZDC KE > 80% of the beam KE
    using KeyKeys = std::set<std::pair<int, int>>;
    const auto keKeys = [](double frac){ // (fEvent, TrackID) of the reconstructed ZDC protons with ZDC KE > frac of the beam KE
        KeyKeys keys;
        ROOT::RDataFrame zdcKe("ParticleInfo", "cuts/cut_ZdcAllEv.root");
        const double cutVal = kBeamKE * frac;
        auto pass = zdcKe.Filter([cutVal](int r, double ke){return r == 1 && ke > cutVal;}, {"isReco", "fDEnergy"});
        auto ev = pass.Take<int>("fEvent"), tid = pass.Take<int>("TrackID");
        for (size_t i = 0; i < ev->size(); i++) keys.emplace_hint(keys.end(), (*ev)[i], (*tid)[i]);
        return keys;
    };
    // only the keys a wanted file needs are read (CUT_ONLY): the jobs of Plot.pbs run the sample files in several processes, and the files without a ZDC KE cut do not wait for the zdc part
    const bool need80 = CutWanted("cuts/cut_megacut2_1.root") || CutWanted("cuts/cut_megacut3_1.root");
    const KeyKeys ke80Keys = need80 ? keKeys(0.80) : KeyKeys();
    // KE scan of 7_megaCut2 and 9_megaCut3: the same cuts as cut_megacut2_1 and cut_megacut3_1 (ZDC KE > 80%) with ZDC KE > 85, 90 and 95% of the beam KE
    // (cut_megacut2_1_ke85.root to cut_megacut3_1_ke95.root; the phi fits are Plot3D.cpp part 4, the overlay with the 80% folder is PlotKEOverlay.cpp). 6_megaCut + KE > x is the megacut2 sample, so it has no files of its own.
    std::vector<KeyKeys> scanKeys;
    scanKeys.reserve(6); // the filters keep a reference into it
    for (const int family : {2, 3}) for (const int maxP : {1}) for (const int pct : {85, 90, 95}){
        const std::string scanFile = Form("cuts/cut_megacut%d_%d_ke%d.root", family, maxP, pct);
        scanKeys.push_back(CutWanted(scanFile) ? keKeys(pct / 100.0) : KeyKeys());
        const KeyKeys& keys = scanKeys.back();
        if (family == 2)
            snap(base.Filter([&keys, maxP](int np, int n, int e, int t){return np <= maxP && (n == 2 || n == 4 || n == 6 || n == 8) && keys.count({e, t}) > 0;},
                             {"nProtonsEv", "nHits", "fEvent", "TrackID"}), scanFile.c_str());
        else
            snap(base.Filter([&keys, maxP](int np, int n, int e, int t, int nTh, int nZ){return np <= maxP && (n == 2 || n == 4 || n == 6 || n == 8) && keys.count({e, t}) > 0 && (nTh == 2 || nTh == 4 || nTh == 6 || nTh == 8) && nZ == 1;},
                             {"nProtonsEv", "nHits", "fEvent", "TrackID", "nTrackerHitsEv", "nZdcCharged"}), scanFile.c_str());
    }
    for (const int maxP : {1}){ // 1_1particlecut of 6_megaCut, 7_megaCut2 and 9_megaCut3: at most 1 reconstructed particle in the event
        snap(base.Filter([maxP](int np, int n){return np <= maxP && (n == 2 || n == 4 || n == 6 || n == 8);},
                         {"nProtonsEv", "nHits"}), Form("cuts/cut_megacut_%d.root", maxP));
        snap(base.Filter([&ke80Keys, maxP](int np, int n, int e, int t){return np <= maxP && (n == 2 || n == 4 || n == 6 || n == 8) && ke80Keys.count({e, t}) > 0;},
                         {"nProtonsEv", "nHits", "fEvent", "TrackID"}), Form("cuts/cut_megacut2_%d.root", maxP));
        snap(base.Filter([&ke80Keys, maxP](int np, int n, int e, int t, int nTh, int nZ){return np <= maxP && (n == 2 || n == 4 || n == 6 || n == 8) && ke80Keys.count({e, t}) > 0
                                                                                                          && (nTh == 2 || nTh == 4 || nTh == 6 || nTh == 8) && nZ == 1;},
                         {"nProtonsEv", "nHits", "fEvent", "TrackID", "nTrackerHitsEv", "nZdcCharged"}), Form("cuts/cut_megacut3_%d.root", maxP)); // 9_megaCut3: megaCut2, 2, 4, 6 or 8 charged tracker hits in the event and exactly 1 charged particle at the ZDC
    }

    // opened only when their file is wanted: cut_ZdcAllEv.root and cut_ZdcAllPart.root are written by the zdc part, which Plot.pbs runs beside the samples without a ZDC KE cut; both stay alive until RunGraphs
    std::unique_ptr<ROOT::RDataFrame> zdcAllEv, zdcAllPart;
    if (CutWanted("cuts/cut_hitZDC_AllEv.root")){
        zdcAllEv.reset(new ROOT::RDataFrame("ParticleInfo", "cuts/cut_ZdcAllEv.root"));
        handles.emplace_back(zdcAllEv->Filter([](int r){return r == 1;}, {"isReco"}).Snapshot("ParticleInfo", "cuts/cut_hitZDC_AllEv.root", "", lazy));
    }
    if (CutWanted("cuts/cut_hitZDC_AllPart.root")){ // every recorded particle at the ZDC that has a reconstructed track, any pdg
        zdcAllPart.reset(new ROOT::RDataFrame("ParticleInfo", "cuts/cut_ZdcAllPart.root"));
        handles.emplace_back(zdcAllPart->Filter([](int r){return r == 1;}, {"isReco"}).Snapshot("ParticleInfo", "cuts/cut_hitZDC_AllPart.root", "", lazy));
    }
    ROOT::RDF::RunGraphs(handles);
}

void CutPlots(const char* partArg = "all"){ // a const char*: calling a macro with a std::string argument fails on the Mac cling (Failed to materialize symbols)
    const std::string part = partArg;
    if (part != "all" && part != "trees" && part != "zdc" && part != "samples"){Error("CutPlots", "part must be all, trees, zdc or samples"); return;}
    // one timing line per part, flushed: the job log then shows where it stopped, and the buffered prints of a part are out before the next one starts
    auto timed = [](const char* name, void (*f)()){
        const time_t t0 = time(nullptr);
        f();
        printf(">>> CutPlots %s took %.0f s\n", name, difftime(time(nullptr), t0)); fflush(stdout);
    };
    const char* only = std::getenv("CUT_ONLY");
    const bool samplesOnly = only && *only; // CUT_ONLY: one cut file is redone, the trees and zdc parts are not
    if ((part == "all" || part == "trees") && !samplesOnly) timed("trees", CutTrees);
    if ((part == "all" || part == "zdc") && !samplesOnly) timed("zdc", ZdcColumns);
    if (part == "all" || part == "samples") timed("samples", CutSamples);
}
