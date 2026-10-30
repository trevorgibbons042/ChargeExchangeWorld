#include <TEntryList.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cstdlib>
#include <string>

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
    for (const auto& ke : KeCuts()){
        TFile out(("cuts/" + ke.dir + ".root").c_str(), "RECREATE");
        for (const auto& j : jobs){
            for (const bool block : {false, true}){
                if (block && std::find(kBlockTrees.begin(), kBlockTrees.end(), j.base) == kBlockTrees.end()) continue;
                std::string cut = KeCutFor(ke.cut, j.t == raw ? j.sfx : ""); // the _DE trees of cut_APandT.root hold the detector KE in finalKE, AnglesRaw in finalKED
                if (block) cut = cut.empty() ? kBlockCut : cut + " && " + kBlockCut;
                if (!j.sel.empty()) cut = cut.empty() ? j.sel : j.sel + " && " + cut;
                const std::string name = j.name + (block ? "_block" : "");
                out.cd();
                if (cut.empty()){ // no cut: every entry passes, so no scan of the tree is needed
                    TEntryList* l = new TEntryList(name.c_str(), "no cut", j.t);
                    const Long64_t n = j.t->GetEntries();
                    for (Long64_t i = 0; i < n; i++) l->Enter(i);
                    l->SetDirectory(&out);
                    continue;
                }
                j.t->Draw((">>" + name).c_str(), cut.c_str(), "entrylist");
                TEntryList* l = (TEntryList*)out.Get(name.c_str());
                if (l) l->SetTitle(cut.c_str());
            }
        }
        out.Write();
        printf("[CutTrees] cuts/%s.root\n", ke.dir.c_str());
    }
}
