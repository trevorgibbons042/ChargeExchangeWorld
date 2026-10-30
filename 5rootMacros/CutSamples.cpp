#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx> // RunGraphs, RResultHandle
#include <ROOT/RSnapshotOptions.hxx>
#include <TMath.h>

#include <string>
#include <vector>

#include "BlockLimits.h"
#include "FitCore.h"
#include "QualityCut.h"

// Writes one cut file per cut of 6_cuts (Plot3D.cpp), each holding only the rows that pass, so Plot3D reads the file instead of cutting all of AnglesRaw with a Draw string.
// AnglesRaw cuts (every fitted proton, no ZDC requirement, cut_angles_all.root) keep the branches Plot3D reads:
//   cut_ResCuts       passQuality and resDRms < 2 mm (3_ResCuts/2mmCutResD)
//   cut_nProtons1/2   at most 1 / 2 protons with a track in the event
//   cut_zOrder        planes hit in increasing z
//   cut_resHalf       residual of the hits after the first 2 over the first 2 hits <= 2
//   cut_nHits         3 to 16 tracker hits
//   cut_firstHit      first tracker hit within +-5 mm in x and y
//   cut_thetacut      recon block theta below the ZDC corner limit of its block (BlockLimits.h)
// ZDC cuts (cut_ZdcCE.root and cut_ZdcAllEv.root, the ZDC protons with the Plot3D columns of ZdcColumns.cpp):
//   cut_hitZDC_CE / _AllEv     the proton has a reconstructed track (isReco)
//   cut_KECUT80_CE / _AllEv    also ZDC KE (fDEnergy) above 80% of the beam KE
// Run after ZdcColumns.cpp, with all the threads: every input is read once and all the cut files of it are written in the same pass.
void CutSamples(){
    ROOT::EnableImplicitMT(NThreads());
    FillBlockLimits(ReadSettings());
    ROOT::RDF::RSnapshotOptions lazy;
    lazy.fLazy = true;
    std::vector<ROOT::RDF::RResultHandle> handles;

    const std::vector<std::string> rawCols = {"fEvent", "TrackID", "isFromCE", "passQuality", "inCE", "finalKE", "resRms", "resDRms", "impact", "impactD", "r2", "r2D", "nProtonsEv", "zOrderOK",
                                              "resHalfRatioD", "outTheta", "outDTheta", "outPhi", "outDPhi", "nHits", "firstXD", "firstYD", "z_bx_reconVertexOut"};
    ROOT::RDataFrame raw("AnglesRaw", "cuts/cut_angles_all.root");
    auto snap = [&](ROOT::RDF::RNode node, const char* file){handles.emplace_back(node.Snapshot("AnglesRaw", file, rawCols, lazy));};
    snap(raw.Filter([](int pq, double resD){return pq == 1 && resD < kRes2Mm;}, {"passQuality", "resDRms"}), "cuts/cut_ResCuts.root");
    snap(raw.Filter([](int n){return n <= 1;}, {"nProtonsEv"}), "cuts/cut_nProtons1.root");
    snap(raw.Filter([](int n){return n <= 2;}, {"nProtonsEv"}), "cuts/cut_nProtons2.root");
    snap(raw.Filter([](int ok){return ok > 0;}, {"zOrderOK"}), "cuts/cut_zOrder.root");
    snap(raw.Filter([](double r){return r <= 2;}, {"resHalfRatioD"}), "cuts/cut_resHalf.root"); // NaN fails the test, as in the Draw cut
    snap(raw.Filter([](int n){return n >= 3 && n <= 16;}, {"nHits"}), "cuts/cut_nHits.root");
    snap(raw.Filter([](double x, double y){return std::fabs(x) <= 5 && std::fabs(y) <= 5;}, {"firstXD", "firstYD"}), "cuts/cut_firstHit.root");
    snap(raw.Filter([](double thetaD, double z){return BlockThetaOk(thetaD * TMath::RadToDeg(), z);}, {"outDTheta", "z_bx_reconVertexOut"}), "cuts/cut_thetacut.root");

    const double minKE = kBeamKE * 0.80;
    ROOT::RDataFrame zdcCE("ParticleInfo", "cuts/cut_ZdcCE.root"), zdcAllEv("ParticleInfo", "cuts/cut_ZdcAllEv.root"); // both stay alive until RunGraphs
    auto zdcCuts = [&](ROOT::RDF::RNode zdc, const std::string& sample){
        auto reco = zdc.Filter([](int r){return r == 1;}, {"isReco"});
        handles.emplace_back(reco.Snapshot("ParticleInfo", "cuts/cut_hitZDC_" + sample + ".root", "", lazy));
        handles.emplace_back(reco.Filter([minKE](double e){return e > minKE;}, {"fDEnergy"}).Snapshot("ParticleInfo", "cuts/cut_KECUT80_" + sample + ".root", "", lazy));
    };
    zdcCuts(zdcCE, "CE");
    zdcCuts(zdcAllEv, "AllEv");
    ROOT::RDF::RunGraphs(handles);
}
