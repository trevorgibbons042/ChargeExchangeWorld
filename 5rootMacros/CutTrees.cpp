#include <TEntryList.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <string>

#include "CutFiles.h"

// Applies every shared cut (KeCuts, CutFiles.h) once to the already cut trees of cut_APandT.root (CalcAPandT.cpp) and cut_angles_all.root
// (ReconAngleAll.cpp), true and detector energy, and writes the entries that pass as TEntryLists to cuts/<ke dir>.root (the residual threshold from RES_CUT_MM, QualityCut.h): one list per tree,
// named like the tree, plus <tree>_block with the block residual cut (kBlockCut) on top for the trees that have it. Run after CalcAPandT.cpp.
void CutTrees(){
    TFile* fAP = TFile::Open("cuts/cut_APandT.root", "READ");
    TFile* fAng = TFile::Open("cuts/cut_angles_all.root", "READ");
    if (!fAP || !fAng){Error("CutTrees", "cuts/cut_APandT.root or cut_angles_all.root missing"); return;}
    std::vector<TTree*> trees;
    for (const std::string sfx : {"", "_DE"}){
        for (const auto& n : kApTrees) trees.push_back((TTree*)fAP->Get((n + sfx).c_str()));
        for (const auto& n : kAngleTrees) trees.push_back((TTree*)fAng->Get((n + sfx).c_str()));
    }
    for (TTree* t : trees) if (!t){Error("CutTrees", "a tree is missing (rerun CalcAPandT.cpp / ReconAngleAll.cpp)"); return;}

    gSystem->mkdir("cuts", true);
    for (const auto& ke : KeCuts()){
        TFile out(("cuts/" + ke.dir + ".root").c_str(), "RECREATE");
        for (TTree* t : trees){
            for (const bool block : {false, true}){
                std::string base = t->GetName();
                if (base.size() > 3 && base.compare(base.size() - 3, 3, "_DE") == 0) base.resize(base.size() - 3);
                if (block && std::find(kBlockTrees.begin(), kBlockTrees.end(), base) == kBlockTrees.end()) continue;
                std::string cut = ke.cut;
                if (block) cut = cut.empty() ? kBlockCut : cut + " && " + kBlockCut;
                const std::string name = std::string(t->GetName()) + (block ? "_block" : "");
                out.cd();
                t->Draw((">>" + name).c_str(), cut.c_str(), "entrylist");
                TEntryList* l = (TEntryList*)out.Get(name.c_str());
                if (l) l->SetTitle(cut.empty() ? "no cut" : cut.c_str());
            }
        }
        out.Write();
        printf("[CutTrees] cuts/%s.root\n", ke.dir.c_str());
    }
}
