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

// The one KE cut: finalKE above 80% of the beam KE (kBeamKE, FitCore.h). keLabel names the energy in the plot titles ("ZDC KE" or "true KE").
static std::vector<KeCut> KeCuts(const std::string& keLabel = "KE"){
    const std::string minKE = std::to_string((int)(kBeamKE * 80 / 100.0));
    return {{"ke_gt80", "finalKE > " + minKE, keLabel + " > " + minKE + " MeV (80%)"}};
}

// Trees CutTrees.cpp writes lists for, per energy suffix ("" or "_DE"). The block lists (name + "_block", block residual cut too) only exist
// for the trees with block position quantities (kBlockTrees).
// cut_angles_all.root holds one tree, AnglesRaw (every fitted proton, ReconAngleAll.cpp). The "trees" Angles and AnglesCE are views of it: the rows with
// AngleSelection(view) set. Their lists are named view + suffix, so the lists of both energies live side by side in one cuts/<ke dir>.root.
static const std::string kAngleTree = "AnglesRaw";
static const std::vector<std::string> kApTrees = {"CE", "CalcValsRe", "CalcValsReD", "CalcValsReDIFF", "CalcValsReDDIFF"}; // cut_APandT.root
static const std::vector<std::string> kAngleTrees = {"Angles", "AnglesCE"};                                               // views of kAngleTree
static const std::vector<std::string> kBlockTrees = {"CalcValsReD", "CalcValsReDDIFF", "Angles", "AnglesCE"};

// The rows of AnglesRaw a view stands for: the protons that pass the quality cut (passQuality), for AnglesCE those of a CE event too (inCE).
static std::string AngleSelection(const std::string& view){return view == "AnglesCE" ? "passQuality == 1 && inCE == 1" : "passQuality == 1";}
// KE branch of an energy suffix: finalKE (true) or finalKED (detector).
static std::string KeVar(const std::string& sfx){return sfx.empty() ? "finalKE" : "finalKED";}
// A KeCuts() cut for the energy of sfx: the cut strings say finalKE, the detector KE is finalKED.
static std::string KeCutFor(std::string cut, const std::string& sfx){
    if (!sfx.empty()) for (size_t p = cut.find("finalKE"); p != std::string::npos; p = cut.find("finalKE", p + 8)) cut.insert(p + 7, "D");
    return cut;
}

// The entries of t that pass KE cut dir (block: and the block residual cut), from cuts/<dir>.root. The list is the one named like t, or name (a view + suffix, AnglesRaw).
static TEntryList* CutList(const std::string& dir, TTree* t, bool block = false, const std::string& listName = ""){
    static std::map<std::string, TFile*> files;
    TFile*& f = files[dir];
    if (!f) f = TFile::Open(("cuts/" + dir + ".root").c_str(), "READ");
    const std::string name = (listName.empty() ? std::string(t->GetName()) : listName) + (block ? "_block" : "");
    TEntryList* l = f ? (TEntryList*)f->Get(name.c_str()) : nullptr;
    if (!l){Error("CutList", "%s missing in cuts/%s.root (run CutTrees.cpp)", name.c_str(), dir.c_str()); return nullptr;}
    // the list belongs to the tree as this macro opened it, so TTree::SetEntryList uses it as is
    l->SetTreeName(t->GetName());
    if (t->GetCurrentFile()) l->SetFileName(t->GetCurrentFile()->GetName());
    return l;
}

// From here on Draw, GetEntries(selection) and Process on t (friends follow) only see the entries of the cut; t->SetEntryList(nullptr) clears it.
// GetEntries() without a selection still counts every entry.
static void ApplyCut(TTree* t, const std::string& dir, bool block = false, const std::string& listName = ""){t->SetEntryList(CutList(dir, t, block, listName));}

// The selection string CutTrees.cpp built the list of (dir, block) from: DrawPlan's listCut for a tree that has that list set. dir is the KeCuts() dir ("ke_gt80").
static std::string ListCut(const std::string& dir, bool block){
    std::string cut;
    for (const auto& ke : KeCuts()) if (ke.dir == dir) cut = ke.cut;
    if (block) cut = cut.empty() ? std::string(kBlockCut) : cut + " && " + kBlockCut;
    return cut;
}
