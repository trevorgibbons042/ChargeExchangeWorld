#pragma once
#include <TEntryList.h>
#include <TError.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <map>
#include <string>
#include <vector>

#include "FitCore.h"
#include "QualityCut.h"

// The cuts shared by the plot macros (PlotAngle.cpp, PlotPolarization.cpp, PlotFomPolarization.cpp, PlotDiffPhiCuts.cpp): every KE cut. CutPlots.cpp (trees) applies them once to the trees of cut_APandT.root and
// cut_angles_all.root and writes the entries that pass as TEntryLists to cuts/<dir>.root; the plot macros set those lists on the trees
// (ApplyCut) instead of cutting with a TTree::Draw string.
struct KeCut {std::string dir, cut, label; double minKE;}; // cut: the selection as a TTree::Draw string (finalKE is the KE of the tree's energy), minKE: the same cut as a number (CutPlots.cpp (trees))

// The one KE cut: finalKE above 80% of the beam KE (kBeamKE, FitCore.h). keLabel names the energy in the plot titles ("ZDC KE" or "true KE").
static std::vector<KeCut> KeCuts(const std::string& keLabel = "KE"){
    const int minKE = (int)(kBeamKE * 80 / 100.0);
    const std::string min = std::to_string(minKE);
    return {{"ke_gt80", "finalKE > " + min, keLabel + " > " + min + " MeV (80%)", (double)minKE}};
}

// Trees CutPlots.cpp (trees) writes lists for, per energy suffix ("" or "_DE"). There is no residual cut: each list is the KE cut only.
// cut_angles_all.root holds one tree, AnglesRaw (every fitted proton, ReconAngleAll.cpp). The "trees" Angles and AnglesCE are views of it: the rows with
// AngleSelection(view) set. Their lists are named view + suffix, so the lists of both energies live side by side in one cuts/<ke dir>.root.
// No namespace-scope globals with a runtime initializer (std::string, std::vector): ACLiC cannot load such a macro on LCRC (Failed to materialize symbols).
static constexpr const char* kAngleTree = "AnglesRaw";
inline const std::vector<std::string>& ApTrees(){static const std::vector<std::string> v = {"CE", "CalcValsRe", "CalcValsReD", "CalcValsReDIFF", "CalcValsReDDIFF"}; return v;} // cut_APandT.root
inline const std::vector<std::string>& AngleTrees(){static const std::vector<std::string> v = {"Angles", "AnglesCE"}; return v;}                                               // views of kAngleTree

// The rows of AnglesRaw a view stands for: every fitted proton for Angles, only those of a CE event (inCE) for AnglesCE.
static std::string AngleSelection(const std::string& view){return view == "AnglesCE" ? "inCE == 1" : "1";}
// KE branch of an energy suffix: finalKE (true) or finalKED (detector).
static std::string KeVar(const std::string& sfx){return sfx.empty() ? "finalKE" : "finalKED";}
// A KeCuts() cut for the energy of sfx: the cut strings say finalKE, the detector KE is finalKED.
static std::string KeCutFor(std::string cut, const std::string& sfx){
    if (!sfx.empty()) for (size_t p = cut.find("finalKE"); p != std::string::npos; p = cut.find("finalKE", p + 8)) cut.insert(p + 7, "D");
    return cut;
}

// The entries of t that pass KE cut dir, from cuts/<dir>.root. The list is the one named like t, or listName (a view + suffix, AnglesRaw).
static TEntryList* CutList(const std::string& dir, TTree* t, const std::string& listName = ""){
    static std::map<std::string, TFile*> files;
    TFile*& f = files[dir];
    if (!f) f = TFile::Open(("cuts/" + dir + ".root").c_str(), "READ");
    const std::string name = listName.empty() ? std::string(t->GetName()) : listName;
    TEntryList* l = f ? (TEntryList*)f->Get(name.c_str()) : nullptr;
    if (!l){Error("CutList", "%s missing in cuts/%s.root (run CutPlots.cpp)", name.c_str(), dir.c_str()); return nullptr;}
    // the list belongs to the tree as this macro opened it, so TTree::SetEntryList uses it as is
    l->SetTreeName(t->GetName());
    // absolute path: TTree::SetEntryList matches the list to the tree's file by name, and a relative name makes RDataFrame (Columns, TreeColumns.h) see an empty list
    if (t->GetCurrentFile()){
        const std::string fn = t->GetCurrentFile()->GetName();
        l->SetFileName(gSystem->IsAbsoluteFileName(fn.c_str()) ? fn.c_str() : gSystem->ConcatFileName(gSystem->WorkingDirectory(), fn.c_str()));
    }
    return l;
}

// From here on Draw, GetEntries(selection) and Process on t (friends follow) only see the entries of the cut; t->SetEntryList(nullptr) clears it.
// GetEntries() without a selection still counts every entry.
static void ApplyCut(TTree* t, const std::string& dir, const std::string& listName = ""){t->SetEntryList(CutList(dir, t, listName));}

