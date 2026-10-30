#include <TEntryList.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

#include "CutFiles.h"

// Applies every shared cut (KeCuts, CutFiles.h) once to the already cut trees of cut_APandT.root (CalcAPandT.cpp) and to the proton tree AnglesRaw of cut_angles_all.root
// (ReconAngleAll.cpp), true and detector energy, and writes the entries that pass as TEntryLists to cuts/<ke dir>.root (the residual threshold from RES_CUT_MM, QualityCut.h): one list per tree,
// named like the tree, plus <tree>_block with the block residual cut (kBlockCut) on top for the trees that have it. Run after CalcAPandT.cpp.
// AnglesRaw gets one list per view (Angles, AnglesCE: the rows AngleSelection picks) and energy, named view + suffix, with the KE cut on that energy's KE branch.
void CutTrees(){
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
    if (!raw->GetBranch("passQuality")){Error("CutTrees", "AnglesRaw has no passQuality branch (rerun ReconAngleAll.cpp)"); return;}

    gSystem->mkdir("cuts", true);
    ROOT::EnableImplicitMT(NThreads());
    for (const auto& ke : KeCuts()){
        // One event loop per tree on all the threads books every list of that tree at once: the entries that pass are collected per thread (with their entry number, rdfentry_),
        // sorted, and entered into the TEntryLists afterwards.
        struct Want {std::string name, title; bool block; int view; bool detector;}; // view: 0 none, 1 Angles (passQuality), 2 AnglesCE (passQuality and inCE)
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
            ROOT::RDataFrame df(n + sfx, "cuts/cut_APandT.root");
            const double minKE = ke.minKE;
            df.ForeachSlot([&](unsigned slot, ULong64_t e, double keVal, double res){
                if (!(keVal > minKE)) return;
                hits[0].perSlot[slot].push_back(e);
                if (hasBlock && PassesBlockResidual(res)) hits[1].perSlot[slot].push_back(e);
            }, {"rdfentry_", "finalKE", "resDRms"});
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
            ROOT::RDataFrame df(kAngleTree, "cuts/cut_angles_all.root");
            const double minKE = ke.minKE;
            df.ForeachSlot([&](unsigned slot, ULong64_t e, int pq, int inCE, double keTrue, double keDet, double res){
                if (!pq) return;
                const bool blockOk = PassesBlockResidual(res);
                for (size_t i = 0; i < wants.size(); i++){
                    const Want& w = wants[i];
                    if (w.view == 2 && !inCE) continue;
                    if (!((w.detector ? keDet : keTrue) > minKE)) continue;
                    if (w.block && !blockOk) continue;
                    hits[i].perSlot[slot].push_back(e);
                }
            }, {"rdfentry_", "passQuality", "inCE", "finalKE", "finalKED", "resDRms"});
            finish(raw, hits);
        }
        out.Write();
        printf("[CutTrees] cuts/%s.root\n", ke.dir.c_str());
    }
}
