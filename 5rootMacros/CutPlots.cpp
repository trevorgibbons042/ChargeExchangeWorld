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
// CutPlots() runs all three; CutPlots("samples") only the last one (plot_cut.sh).

// Applies every shared cut (KeCuts, CutFiles.h) once to the already cut trees of cut_APandT.root (CalcAPandT.cpp) and to the proton tree AnglesRaw of cut_angles_all.root
// (ReconAngleAll.cpp), true and detector energy, and writes the entries that pass as TEntryLists to cuts/<ke dir>.root (the residual threshold from RES_CUT_MM, QualityCut.h): one list per tree,
// named like the tree, plus <tree>_block with the block residual cut (kBlockCut) on top for the trees that have it. Run after CalcAPandT.cpp.
// AnglesRaw gets one list per view (Angles, AnglesCE: the rows AngleSelection picks) and energy, named view + suffix, with the KE cut on that energy's KE branch.
static void CutTrees(){
    TFile* fAP = TFile::Open("cuts/cut_APandT.root", "READ");
    TFile* fAng = TFile::Open("cuts/cut_angles_all.root", "READ");
    if (!fAP || !fAng){Error("CutTrees", "cuts/cut_APandT.root or cut_angles_all.root missing"); return;}
    // One list to build: the tree it indexes, its name (without _block), the rows it covers (sel) and the energy suffix of its KE cut.
    struct Job {TTree* t; std::string name, base, sel, sfx;};
    std::vector<Job> jobs;
    TTree* raw = (TTree*)fAng->Get(kAngleTree.c_str());
    for (const std::string sfx : {"", "_DE"}){
        for (const auto& n : kApTrees) jobs.push_back({(TTree*)fAP->Get((n + sfx).c_str()), n + sfx, n, "", sfx});
        for (const auto& n : kAngleTrees) jobs.push_back({raw, n + sfx, n, AngleSelection(n), sfx});
    }
    for (const auto& j : jobs) if (!j.t){Error("CutTrees", "a tree is missing (rerun CalcAPandT.cpp / ReconAngleAll.cpp)"); return;}

    gSystem->mkdir("cuts", true);
    ROOT::EnableImplicitMT(NThreads());
    for (const auto& ke : KeCuts()){
        // One event loop per tree on all the threads books every list of that tree at once: the entries that pass are collected per thread (with their entry number, rdfentry_),
        // sorted, and entered into the TEntryLists afterwards.
        struct Want {std::string name, title; bool block; int view; bool detector;}; // view: 0 none, 1 Angles (every row), 2 AnglesCE (inCE)
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
        // cut_APandT.root: one tree per name and energy, every list on finalKE (the _DE trees hold the detector KE there) and resDRms
        for (const std::string sfx : {"", "_DE"}) for (const auto& n : kApTrees){
            TTree* t = (TTree*)fAP->Get((n + sfx).c_str());
            const bool hasBlock = std::find(kBlockTrees.begin(), kBlockTrees.end(), n) != kBlockTrees.end();
            std::vector<Hit> hits(hasBlock ? 2 : 1);
            out.cd();
            hits[0] = {new TEntryList((n + sfx).c_str(), ke.cut.c_str(), t), std::vector<std::vector<ULong64_t>>(nSlots)};
            if (hasBlock) hits[1] = {new TEntryList((n + sfx + "_block").c_str(), (ke.cut + " && " + kBlockCut).c_str(), t), std::vector<std::vector<ULong64_t>>(nSlots)};
            ROOT::RDataFrame df0(n + sfx, "cuts/cut_APandT.root");
            auto df = WithEntryNumber(df0, df0.GetNSlots()); // the true entry numbers (rdfentry_ is wrong under MT on a tree with several clusters)
            const double minKE = ke.minKE;
            df.ForeachSlot([&](unsigned slot, ULong64_t e, double keVal, double res){
                if (!(keVal > minKE)) return;
                hits[0].perSlot[slot].push_back(e);
                if (hasBlock && PassesBlockResidual(res)) hits[1].perSlot[slot].push_back(e);
            }, {"entryNo", "finalKE", "resDRms"});
            finish(t, hits);
        }
        // AnglesRaw: the views Angles and AnglesCE of both energies in one loop
        {
            std::vector<Want> wants;
            for (const std::string sfx : {"", "_DE"}) for (const auto& v : kAngleTrees) for (const bool block : {false, true})
                wants.push_back({v + sfx + (block ? "_block" : ""), "", block, v == "AnglesCE" ? 2 : 1, !sfx.empty()});
            std::vector<Hit> hits(wants.size());
            out.cd();
            for (size_t i = 0; i < wants.size(); i++){
                const std::string cut = AngleSelection(wants[i].view == 2 ? "AnglesCE" : "Angles") + " && " + KeCutFor(ke.cut, wants[i].detector ? "_DE" : "") + (wants[i].block ? " && " + kBlockCut : "");
                hits[i] = {new TEntryList(wants[i].name.c_str(), cut.c_str(), raw), std::vector<std::vector<ULong64_t>>(nSlots)};
            }
            ROOT::RDataFrame df0(kAngleTree, "cuts/cut_angles_all.root");
            auto df = WithEntryNumber(df0, df0.GetNSlots());
            const double minKE = ke.minKE;
            df.ForeachSlot([&](unsigned slot, ULong64_t e, int inCE, double keTrue, double keDet, double res){
                const bool blockOk = PassesBlockResidual(res);
                for (size_t i = 0; i < wants.size(); i++){
                    const Want& w = wants[i];
                    if (w.view == 2 && !inCE) continue;
                    if (!((w.detector ? keDet : keTrue) > minKE)) continue;
                    if (w.block && !blockOk) continue;
                    hits[i].perSlot[slot].push_back(e);
                }
            }, {"entryNo", "inCE", "finalKE", "finalKED", "resDRms"});
            finish(raw, hits);
        }
        out.Write();
        printf("[CutTrees] cuts/%s.root\n", ke.dir.c_str());
    }
}

// Writes the two ZDC proton cut files Plot3D.cpp draws from: cut_ZdcCE.root (the ZDC protons of the CE children, cut_ParentZdc.root) and cut_ZdcAllEv.root (the ZDC protons
// of every event, cut_AllEvProtons.root), each row with the Plot3D columns of its reconstructed track (isReco, isRecoD, isRes2, raw_0..16, ok_0..16, Plot3DCols.h).
// Run after CutEvents.cpp and ReconAngleAll.cpp, with all the threads: the maps of AnglesRaw are read once (serial), then both files are written in one parallel pass.
static void ZdcColumns(){
    ROOT::EnableImplicitMT(NThreads());
    ROOT::RDF::RSnapshotOptions lazy;
    lazy.fLazy = true;

    const ColumnAdder addColumns = Plot3DColumns(); // loads the maps on its first call
    ROOT::RDataFrame dfCE("ParticleInfo", "cuts/cut_ParentZdc.root"), dfAllEv("ParticleInfo", "cuts/cut_AllEvProtons.root"), dfAllPart("ParticleInfo", "cuts/cut_AllEvParticles.root");
    std::vector<ROOT::RDF::RResultHandle> handles;
    handles.emplace_back(addColumns(dfCE).Snapshot("ParticleInfo", "cuts/cut_ZdcCE.root", "", lazy));
    handles.emplace_back(addColumns(dfAllEv).Snapshot("ParticleInfo", "cuts/cut_ZdcAllEv.root", "", lazy));
    handles.emplace_back(addColumns(dfAllPart).Snapshot("ParticleInfo", "cuts/cut_ZdcAllPart.root", "", lazy)); // the same for every recorded particle, any pdg (the "all recon particles" plots)
    ROOT::RDF::RunGraphs(handles);
}

// Writes one cut file per cut of 6_cuts (Plot3D.cpp), each holding only the rows that pass, so Plot3D reads the file instead of cutting all of AnglesRaw with a Draw string.
// AnglesRaw cuts (every fitted proton, no ZDC requirement, cut_angles_all.root) keep the branches Plot3D reads:
//   cut_ResCuts       resDRms < 4 mm (3_ResCuts)
//   cut_nParticles1/2 at most 1 / 2 reconstructed particles in the event
//   cut_nZdc1         at most 1 charged particle at the ZDC in the event (8_nZdcParticles, EventSummary of the sim)
//   cut_nHits         2, 4, 6 or 8 tracker hits
//   cut_megacut_1..4  6_megaCut/1_particleCut to 4_particleCut: at most 1 (1, 3) / 2 (2, 4) reconstructed particles in the event, 2, 4, 6 or 8 tracker hits, recon theta and recon block theta 0 to 1 degree, and for 3 and 4 also resDRms < 4 mm
//   cut_megacut2_1..4 7_megaCut2/1_particleCut to 4_particleCut: the same, and a reconstructed ZDC proton with ZDC KE > 80% of the beam KE
//   cut_megacut3_1..4 9_megaCut3/1_1particlecut to 4_2particleREScut: the cut_megacut2 cuts, and 2, 4, 6 or 8 charged tracker hits in the event (nTrackerHitsEv, all particles added up)
// ZDC cut (cut_ZdcAllEv.root, the ZDC protons with the Plot3D columns of ZdcColumns.cpp):
//   cut_hitZDC_AllEv           the proton has a reconstructed track (isReco)
// Run after ZdcColumns.cpp, with all the threads: every input is read once and all the cut files of it are written in the same pass.
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
    auto snap = [&](ROOT::RDF::RNode node, const char* file){handles.emplace_back(node.Snapshot("AnglesRaw", file, rawCols, lazy));};
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
    const KeyKeys ke80Keys = keKeys(0.80);
    // KE scan of 9_megaCut3/1_1particlecut: the same cuts with ZDC KE > 85, 90 and 95% of the beam KE (cut_megacut3_1_ke85.root to _ke95.root, phi fits in Plot3D.cpp part 4)
    std::vector<KeyKeys> scanKeys;
    scanKeys.reserve(3); // the filters keep a reference into it
    for (const int pct : {85, 90, 95}){
        scanKeys.push_back(keKeys(pct / 100.0));
        const KeyKeys& keys = scanKeys.back();
        snap(base.Filter([&keys](int np, int n, int e, int t, int nTh){return np <= 1 && (n == 2 || n == 4 || n == 6 || n == 8) && keys.count({e, t}) > 0 && (nTh == 2 || nTh == 4 || nTh == 6 || nTh == 8);},
                         {"nProtonsEv", "nHits", "fEvent", "TrackID", "nTrackerHitsEv"}), Form("cuts/cut_megacut3_1_ke%d.root", pct));
    }
    for (const int i : {1, 2, 3, 4}){ // 1_particleCut to 4_particleCut of 6_megaCut and 7_megaCut2: at most 1 (1, 3) / 2 (2, 4) reconstructed particles in the event, 3 and 4 also resDRms < 4 mm
        const int maxP = i % 2 ? 1 : 2;
        const bool res4 = i > 2;
        snap(base.Filter([maxP, res4](int np, int n, double resD){return np <= maxP && (n == 2 || n == 4 || n == 6 || n == 8) && (!res4 || resD < kRes4Mm);},
                         {"nProtonsEv", "nHits", "resDRms"}), Form("cuts/cut_megacut_%d.root", i));
        snap(base.Filter([&ke80Keys, maxP, res4](int np, int n, double resD, int e, int t){return np <= maxP && (n == 2 || n == 4 || n == 6 || n == 8) && (!res4 || resD < kRes4Mm) && ke80Keys.count({e, t}) > 0;},
                         {"nProtonsEv", "nHits", "resDRms", "fEvent", "TrackID"}), Form("cuts/cut_megacut2_%d.root", i));
        snap(base.Filter([&ke80Keys, maxP, res4](int np, int n, double resD, int e, int t, int nTh){return np <= maxP && (n == 2 || n == 4 || n == 6 || n == 8) && (!res4 || resD < kRes4Mm) && ke80Keys.count({e, t}) > 0
                                                                                                          && (nTh == 2 || nTh == 4 || nTh == 6 || nTh == 8);},
                         {"nProtonsEv", "nHits", "resDRms", "fEvent", "TrackID", "nTrackerHitsEv"}), Form("cuts/cut_megacut3_%d.root", i)); // 9_megaCut3: megaCut2 and 2, 4, 6 or 8 charged tracker hits in the event
    }

    ROOT::RDataFrame zdcAllEv("ParticleInfo", "cuts/cut_ZdcAllEv.root"); // stays alive until RunGraphs
    handles.emplace_back(zdcAllEv.Filter([](int r){return r == 1;}, {"isReco"}).Snapshot("ParticleInfo", "cuts/cut_hitZDC_AllEv.root", "", lazy));
    ROOT::RDataFrame zdcAllPart("ParticleInfo", "cuts/cut_ZdcAllPart.root"); // every recorded particle at the ZDC that has a reconstructed track, any pdg
    handles.emplace_back(zdcAllPart.Filter([](int r){return r == 1;}, {"isReco"}).Snapshot("ParticleInfo", "cuts/cut_hitZDC_AllPart.root", "", lazy));
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
    if (part == "all" || part == "trees") timed("trees", CutTrees);
    if (part == "all" || part == "zdc") timed("zdc", ZdcColumns);
    if (part == "all" || part == "samples") timed("samples", CutSamples);
}
