#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx> // RunGraphs, RResultHandle
#include <ROOT/RSnapshotOptions.hxx>

#include <string>
#include <vector>

#include "FitCore.h"
#include "Plot3DCols.h"

// Writes the two ZDC proton cut files Plot3D.cpp draws from: cut_ZdcCE.root (the ZDC protons of the CE children, cut_ParentZdc.root) and cut_ZdcAllEv.root (the ZDC protons
// of every event, cut_AllEvProtons.root), each row with the Plot3D columns of its reconstructed track (isReco, isRecoD, isRes2, thetaOk, raw_0..16, ok_0..16, Plot3DCols.h).
// Run after CutEvents.cpp and ReconAngleAll.cpp, with all the threads: the maps of AnglesRaw are read once (serial), then both files are written in one parallel pass.
void ZdcColumns(){
    ROOT::EnableImplicitMT(NThreads());
    ROOT::RDF::RSnapshotOptions lazy;
    lazy.fLazy = true;

    const ColumnAdder addColumns = Plot3DColumns(); // loads the maps on its first call
    ROOT::RDataFrame dfCE("ParticleInfo", "cuts/cut_ParentZdc.root"), dfAllEv("ParticleInfo", "cuts/cut_AllEvProtons.root");
    std::vector<ROOT::RDF::RResultHandle> handles;
    handles.emplace_back(addColumns(dfCE).Snapshot("ParticleInfo", "cuts/cut_ZdcCE.root", "", lazy));
    handles.emplace_back(addColumns(dfAllEv).Snapshot("ParticleInfo", "cuts/cut_ZdcAllEv.root", "", lazy));
    ROOT::RDF::RunGraphs(handles);
}
