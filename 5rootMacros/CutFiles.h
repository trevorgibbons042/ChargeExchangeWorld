#pragma once
#include <TError.h>
#include <TFile.h>
#include <TTree.h>

#include <map>
#include <string>

#include "FitCore.h"

// cut_angles_all.root holds one tree, AnglesRaw (every fitted proton, ReconAngleAll.cpp). The "trees" Angles and AnglesCE are views of it: the rows with
// AngleSelection(view) set.
// No namespace-scope globals with a runtime initializer (std::string, std::vector): ACLiC cannot load such a macro on LCRC (Failed to materialize symbols).
static constexpr const char* kAngleTree = "AnglesRaw";

// The rows of AnglesRaw a view stands for: every fitted proton for Angles, only those of a CE event (inCE) for AnglesCE.
static std::string AngleSelection(const std::string& view){return view == "AnglesCE" ? "inCE == 1" : "1";}
// KE branch of an energy suffix: zdcKEAll (true) or zdcKEAllD (detector), the ZDC KE sum of all particles of the event.
static std::string KeVar(const std::string& sfx){return sfx.empty() ? "zdcKEAll" : "zdcKEAllD";}

// The fixed cut files CutPlots.cpp (events), ReconAngleAll.cpp and CutPlots.cpp write into cuts/: each holds only the rows that pass the cut in its name, so the plot macros
// read them as they are, with no cut string. CutPlots.cpp (events) chains them: cut_CE and cut_CEall come from the sim output, cut_ParentTracker from cut_CE.
//   cut_CE.root            ChargeExchange_values: the carbon-block CE (CarbonCeCut, FitCore.h)         ParticleInfo: every row of the events that have one
//   cut_CEall.root         ChargeExchange_values: every CE of the sim
//   cut_ParentTracker.root ParticleInfo: the children of a carbon CE track that are protons at the first tracker (boundaryType == 0)
//   cut_Plot3DCE.root, cut_Plot3DAllEv.root, cut_Plot3DAllPart.root   AnglesRaw rows with evClass 0, 0 or 1, any, plus the Plot3D columns (Plot3DCols.h, CutPlots.cpp (zdc))

// The tree of a cut file, opened once and kept open (READ); nullptr with an Error if the file or tree is missing.
inline TTree* OpenCut(const std::string& file, const std::string& tree){
    static std::map<std::string, TFile*> files;
    TFile*& f = files[file];
    if (!f) f = TFile::Open(file.c_str(), "READ");
    TTree* t = f && !f->IsZombie() ? (TTree*)f->Get(tree.c_str()) : nullptr;
    if (!t) Error("OpenCut", "tree %s missing in %s (run CutPlots.cpp (events) / ReconAngleAll.cpp / CutPlots.cpp)", tree.c_str(), file.c_str());
    return t;
}
