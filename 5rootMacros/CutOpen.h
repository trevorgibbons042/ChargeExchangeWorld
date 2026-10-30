#pragma once
#include <TError.h>
#include <TFile.h>
#include <TTree.h>

#include <map>
#include <string>

// The fixed cut files CutEvents.cpp, ReconAngleAll.cpp and CutPlots.cpp write into cuts/: each holds only the rows that pass the cut in its name, so the plot macros
// read them as they are, with no cut string. CutEvents.cpp chains them: cut_CE, cut_AllEvProtons and cut_CEall come from the sim output, cut_ParentZdc and cut_ParentTracker from cut_CE.
//   cut_CE.root            ChargeExchange_values: the carbon-block CE (CarbonCeCut, FitCore.h)         ParticleInfo: every row of the events that have one
//   cut_CEall.root         ChargeExchange_values: every CE of the sim
//   cut_ParentZdc.root     ParticleInfo: the children of a carbon CE track that are protons at the ZDC (boundaryType == 2, pdg == 2212)
//   cut_ParentTracker.root ParticleInfo: the children of a carbon CE track that are protons at the first tracker (boundaryType == 0)
//   cut_AllEvProtons.root  ParticleInfo: the protons at the ZDC of every event
//   cut_ZdcCE.root, cut_ZdcAllEv.root   cut_ParentZdc.root / cut_AllEvProtons.root plus the Plot3D columns of the reconstructed track (Plot3DCols.h)
//   cut_ResCuts, cut_nProtons1, cut_nProtons2, cut_nHits, cut_megacut .root  AnglesRaw: the fitted protons that pass one cut of 6_cuts (CutPlots.cpp (samples))
//   cut_hitZDC_AllEv.root   ParticleInfo: the ZDC protons of every event with a reconstructed track (CutPlots.cpp (samples))

// The tree of a cut file, opened once and kept open (READ); nullptr with an Error if the file or tree is missing.
inline TTree* OpenCut(const std::string& file, const std::string& tree){
    static std::map<std::string, TFile*> files;
    TFile*& f = files[file];
    if (!f) f = TFile::Open(file.c_str(), "READ");
    TTree* t = f && !f->IsZombie() ? (TTree*)f->Get(tree.c_str()) : nullptr;
    if (!t) Error("OpenCut", "tree %s missing in %s (run CutEvents.cpp / ReconAngleAll.cpp / CutPlots.cpp)", tree.c_str(), file.c_str());
    return t;
}
