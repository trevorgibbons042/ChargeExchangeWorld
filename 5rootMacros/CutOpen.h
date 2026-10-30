#pragma once
#include <TError.h>
#include <TFile.h>
#include <TTree.h>

#include <map>
#include <string>

// The fixed cut files CutEvents.cpp, ReconAngleAll.cpp and ZdcColumns.cpp write into cuts/: each holds only the rows that pass the cut in its name, so the plot macros
// read them as they are, with no cut string.
//   cut_CE.root            ChargeExchange_values: the carbon-block CE (CarbonCeCut, FitCore.h)         ParticleInfo: every row of the events that have one
//   cut_CEtracker.root     ChargeExchange_values: the CE outside the carbon blocks (the silicon trackers)
//   cut_CEall.root         ChargeExchange_values: every CE of the sim
//   cut_ParentID.root      ParticleInfo: the children of a carbon CE track
//   cut_ParentZdc.root     ParticleInfo: the children of a carbon CE track that are protons at the ZDC (boundaryType == 2, pdg == 2212)
//   cut_ParentTracker.root ParticleInfo: the children of a carbon CE track that are protons at the first tracker (boundaryType == 0)
//   cut_AllEvProtons.root  ParticleInfo: the protons at the ZDC of every event
//   cut_ZdcCE.root, cut_ZdcAllEv.root   cut_ParentZdc.root / cut_AllEvProtons.root plus the Plot3D columns of the reconstructed track (Plot3DCols.h)
//   cut_ResCuts, cut_nProtons1, cut_nProtons2, cut_zOrder, cut_resHalf, cut_nHits, cut_firstHit, cut_thetacut .root   AnglesRaw: the fitted protons that pass one cut of 6_cuts (CutSamples.cpp)
//   cut_hitZDC_CE/_AllEv.root, cut_KECUT80_CE/_AllEv.root   ParticleInfo: the ZDC protons with a reconstructed track, and those above 80% of the beam KE (CutSamples.cpp)

// The tree of a cut file, opened once and kept open (READ); nullptr with an Error if the file or tree is missing.
inline TTree* OpenCut(const std::string& file, const std::string& tree){
    static std::map<std::string, TFile*> files;
    TFile*& f = files[file];
    if (!f) f = TFile::Open(file.c_str(), "READ");
    TTree* t = f && !f->IsZombie() ? (TTree*)f->Get(tree.c_str()) : nullptr;
    if (!t) Error("OpenCut", "tree %s missing in %s (run CutEvents.cpp / ReconAngleAll.cpp / ZdcColumns.cpp)", tree.c_str(), file.c_str());
    return t;
}
