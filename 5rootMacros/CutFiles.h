#pragma once
#include <TEntryList.h>
#include <TError.h>
#include <TFile.h>
#include <TTree.h>

#include <map>
#include <string>
#include <vector>

#include "FitCore.h"
#include "QualityCut.h"

// The cuts shared by the plot macros (PlotAngle.cpp, PlotPolarization.cpp, PlotFomPolarization.cpp, PlotDiffPhiCuts.cpp): every KE cut. CutTrees.cpp applies them once to the trees of cut_APandT.root and
// cut_angles_all.root and writes the entries that pass as TEntryLists to cuts/<dir>.root; the plot macros set those lists on the trees
// (ApplyCut) instead of cutting with a TTree::Draw string.
struct KeCut {std::string dir, cut, label;}; // cut: the selection CutTrees.cpp applies, as a TTree::Draw string (finalKE is the KE of the tree's energy)

// No KE cut, finalKE above 50 to 99% of the beam KE (kBeamKE, FitCore.h), then each again with the 2 mm block residual cut on top (QualityCut.h).
// keLabel names the energy in the plot titles ("ZDC KE" or "true KE").
static std::vector<KeCut> KeCuts(const std::string& keLabel = "KE"){
    std::vector<KeCut> kes = {{"ke_all", "", "no " + keLabel + " cut"}};
    for (int pct : {50, 60, 70, 80, 90, 95, 99}){
        const std::string minKE = std::to_string((int)(kBeamKE * pct / 100.0));
        kes.push_back({"ke_gt" + std::to_string(pct), "finalKE > " + minKE, keLabel + " > " + minKE + " MeV (" + std::to_string(pct) + "%)"});
    }
    return kes;
}

// Trees CutTrees.cpp writes lists for, per energy suffix ("" or "_DE"). The block lists (name + "_block", block residual cut too) only exist
// for the trees with block position quantities (kBlockTrees).
static const std::vector<std::string> kApTrees = {"CE", "CalcValsRe", "CalcValsReD", "CalcValsReDIFF", "CalcValsReDDIFF"}; // cut_APandT.root
static const std::vector<std::string> kAngleTrees = {"Angles", "AnglesCE"};                                               // cut_angles_all.root
static const std::vector<std::string> kBlockTrees = {"CalcValsReD", "CalcValsReDDIFF", "Angles", "AnglesCE"};

// The entries of t that pass KE cut dir (block: and the block residual cut), from cuts/<dir>.root.
static TEntryList* CutList(const std::string& dir, TTree* t, bool block = false){
    static std::map<std::string, TFile*> files;
    TFile*& f = files[dir];
    if (!f) f = TFile::Open(("cuts/" + dir + ".root").c_str(), "READ");
    const std::string name = std::string(t->GetName()) + (block ? "_block" : "");
    TEntryList* l = f ? (TEntryList*)f->Get(name.c_str()) : nullptr;
    if (!l){Error("CutList", "%s missing in cuts/%s.root (run CutTrees.cpp)", name.c_str(), dir.c_str()); return nullptr;}
    // the list belongs to the tree as this macro opened it, so TTree::SetEntryList uses it as is
    l->SetTreeName(t->GetName());
    if (t->GetCurrentFile()) l->SetFileName(t->GetCurrentFile()->GetName());
    return l;
}

// From here on Draw, GetEntries(selection) and Process on t (friends follow) only see the entries of the cut; t->SetEntryList(nullptr) clears it.
// GetEntries() without a selection still counts every entry.
static void ApplyCut(TTree* t, const std::string& dir, bool block = false){t->SetEntryList(CutList(dir, t, block));}
