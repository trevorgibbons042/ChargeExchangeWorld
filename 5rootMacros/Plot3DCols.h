#pragma once
// The Plot3D maps and the cached "p3d_v1" columns (CutCache.h tag) that replace the per-entry isReco/rawVal/rawOk/thetaOk calls of the TTree::Draw cuts.
// Shared by Plot3D.cpp and BuildCaches.cpp.
#include <TFile.h>
#include <TError.h>
#include <TMath.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <utility>

#include "BlockLimits.h"
#include "CutCache.h"
#include "FitCore.h"
#include "QualityCut.h"

// (fEvent, TrackID) of every proton with a reconstructed track (cut_angles_all.root, Angles). The position plots only show these protons,
// the ones that can actually be drawn on the ZDC.
// gRecoTracksD: the ones that also pass the block residual cut (PassesBlockResidual, QualityCut.h), for the block position plots only.
// gRes2Tracks: block residual under 2 mm, the _res2mm twin of every KE cut (kRes2Mm, QualityCut.h).
inline std::set<std::pair<int, int>> gRecoTracks, gRecoTracksD, gRes2Tracks;
inline bool isReco(int fEvent, int trackID){return gRecoTracks.count({fEvent, trackID}) > 0;}
inline bool isRes2(int fEvent, int trackID){return gRes2Tracks.count({fEvent, trackID}) > 0;}
inline bool isRecoD(int fEvent, int trackID){return gRecoTracksD.count({fEvent, trackID}) > 0;}

// Cut candidates of every fitted proton (cut_angles_all.root, AnglesRaw, before any cut) by (fEvent, TrackID), for the 1_additionalinfo plots of the KE cut folders.
// Index: see kInfoVars in Plot3D (14, 15: recon phi and recon block phi in degrees, 0 to 360). rawOk is the Draw cut that keeps the protons that have a finite value.
inline std::map<std::pair<int, int>, std::array<double, 17>> gRawVals; // 16: recon vertex z (mm), from the block position line fit
inline bool rawOk(int fEvent, int trackID, int k){
    auto it = gRawVals.find({fEvent, trackID});
    return it != gRawVals.end() && std::isfinite(it->second[k]);
}
inline double rawVal(int fEvent, int trackID, int k){return gRawVals.at({fEvent, trackID})[k];}

// True if the recon block theta (k 11) is within the limit of the block nearest the recon vertex z (k 16), BlockLimits.h.
inline bool thetaOk(int fEvent, int trackID){
    const auto it = gRawVals.find({fEvent, trackID});
    return it != gRawVals.end() && BlockThetaOk(it->second[11], it->second[16]);
}

// Fills the maps above from cuts/cut_angles_all.root and the block limits from settings.csv. False (with an Error) if the file or tree is missing.
inline bool LoadPlot3DMaps(){
    TFile* fReco = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree* reco = fReco ? (TTree*)fReco->Get("Angles") : nullptr;
    if (!reco){Error("Plot3D", "missing Angles tree in cut_angles_all.root (run ReconAngleAll.cpp first)"); return false;}
    int ev, tid;
    double resDRms = 0;
    reco->SetBranchAddress("fEvent", &ev); reco->SetBranchAddress("TrackID", &tid);
    reco->SetBranchAddress("resDRms", &resDRms);
    for (Long64_t i = 0; i < reco->GetEntries(); i++){
        reco->GetEntry(i);
        gRecoTracks.insert({ev, tid});
        if (PassesBlockResidual(resDRms)) gRecoTracksD.insert({ev, tid});
        if (resDRms < kRes2Mm) gRes2Tracks.insert({ev, tid});
    }
    TTree* raw = (TTree*)fReco->Get("AnglesRaw");
    if (raw && raw->GetBranch("nProtonsEv")){
        double v[17]; double thetaRad, thetaBlockRad, phiRad, phiBlockRad; int nP, zOk, nH;
        raw->SetBranchAddress("fEvent", &ev); raw->SetBranchAddress("TrackID", &tid);
        const char* dn[6] = {"resRms", "resDRms", "impact", "impactD", "r2", "r2D"};
        for (int i = 0; i < 6; i++) raw->SetBranchAddress(dn[i], &v[i]);
        raw->SetBranchAddress("nProtonsEv", &nP); raw->SetBranchAddress("zOrderOK", &zOk); raw->SetBranchAddress("nHits", &nH);
        raw->SetBranchAddress("outPhi", &phiRad); raw->SetBranchAddress("outDPhi", &phiBlockRad);
        raw->SetBranchAddress("firstXD", &v[12]); raw->SetBranchAddress("firstYD", &v[13]);
        raw->SetBranchAddress("z_bx_reconVertexOut", &v[16]); raw->SetBranchAddress("resHalfRatioD", &v[8]); raw->SetBranchAddress("outTheta", &thetaRad); raw->SetBranchAddress("outDTheta", &thetaBlockRad);
        for (Long64_t i = 0; i < raw->GetEntries(); i++){
            raw->GetEntry(i);
            v[14] = (phiRad < 0 ? phiRad + 2 * TMath::Pi() : phiRad) * TMath::RadToDeg(); v[15] = (phiBlockRad < 0 ? phiBlockRad + 2 * TMath::Pi() : phiBlockRad) * TMath::RadToDeg();
            v[6] = nP; v[7] = zOk; v[9] = thetaRad * TMath::RadToDeg(); v[10] = nH; v[11] = thetaBlockRad * TMath::RadToDeg();
            std::array<double, 17> a; std::copy(v, v + 17, a.begin());
            gRawVals[{ev, tid}] = a;
        }
    }
    fReco->Close();
    FillBlockLimits(ReadSettings());
    return true;
}

inline const char* kPlot3DTag = "p3d_v1"; // bump when the column logic below changes
// A change to the cuts or limits the columns are built from makes the caches stale: the macros are rsynced into the run folder, which updates their mtime.
inline const std::vector<std::string> kPlot3DTagDeps = {"cuts/cut_angles_all.root", "QualityCut.h", "BlockLimits.h", "Plot3DCols.h"};

// One map lookup per row and thread: the 40 columns below are defined in a row and all ask for the same (fEvent, TrackID).
inline const std::array<double, 17>* P3Row(int ev, int tid){
    thread_local int lastEv = 0, lastTid = 0; thread_local const std::array<double, 17>* lastRow = nullptr; thread_local bool valid = false;
    if (valid && lastEv == ev && lastTid == tid) return lastRow;
    auto it = gRawVals.find({ev, tid});
    lastEv = ev; lastTid = tid; lastRow = it == gRawVals.end() ? nullptr : &it->second; valid = true;
    return lastRow;
}

// Adds isReco, isRecoD, isRes2, thetaOk (Int_t 0/1), raw_0..raw_16 (Double_t, NaN if the key is missing) and ok_0..ok_16 (Int_t, the rawOk logic).
// The maps are read only while the event loop runs, so the lambdas are safe on several threads.
inline ColumnAdder Plot3DColumns(){
    return [](ROOT::RDF::RNode n){
        static const bool mapsLoaded = LoadPlot3DMaps(); // only when a cache is built: a cache hit never reads the 213M rows
        if (!mapsLoaded) throw std::runtime_error("Plot3DColumns: the maps could not be loaded (cuts/cut_angles_all.root)");
        n = n.Define("isReco", [](int ev, int tid){return (Int_t)isReco(ev, tid);}, {"fEvent", "TrackID"})
             .Define("isRecoD", [](int ev, int tid){return (Int_t)isRecoD(ev, tid);}, {"fEvent", "TrackID"})
             .Define("isRes2", [](int ev, int tid){return (Int_t)isRes2(ev, tid);}, {"fEvent", "TrackID"})
             .Define("thetaOk", [](int ev, int tid){
                 const auto* r = P3Row(ev, tid);
                 return (Int_t)(r && BlockThetaOk((*r)[11], (*r)[16]));
             }, {"fEvent", "TrackID"});
        for (int k = 0; k < 17; k++){
            n = n.Define("raw_" + std::to_string(k), [k](int ev, int tid){
                const auto* r = P3Row(ev, tid);
                return r ? (*r)[k] : std::numeric_limits<double>::quiet_NaN();
            }, {"fEvent", "TrackID"});
            n = n.Define("ok_" + std::to_string(k), [k](int ev, int tid){
                const auto* r = P3Row(ev, tid);
                return (Int_t)(r && std::isfinite((*r)[k]));
            }, {"fEvent", "TrackID"});
        }
        return n;
    };
}

// Rewrites the isReco/isRecoD/isRes2/thetaOk/rawVal/rawOk calls of a TTree::Draw expression or cut into the cached columns.
inline std::string P3(const std::string& s){
    static const std::regex fn(R"(\b(isRecoD|isReco|isRes2|thetaOk)\s*\(\s*fEvent\s*,\s*TrackID\s*\))");
    static const std::regex val(R"(\brawVal\s*\(\s*fEvent\s*,\s*TrackID\s*,\s*(\d+)\s*\))");
    static const std::regex ok(R"(\brawOk\s*\(\s*fEvent\s*,\s*TrackID\s*,\s*(\d+)\s*\))");
    std::string o = std::regex_replace(s, fn, "($1==1)");
    o = std::regex_replace(o, val, "raw_$1");
    o = std::regex_replace(o, ok, "(ok_$1==1)");
    return o;
}

// The cached trees of Plot3D (spec A.7 C7 to C9), see BuildCaches.cpp. The C7/C8 sets are the ZDC protons of the CE children and of every event.
inline const std::vector<std::string> kZdcProtonAtoms = {"boundaryType == 2 && pdg == 2212"};
inline CacheTree Plot3DCE(){ // C7: CE protons at the ZDC (cut_ParentID) + columns
    return GetCutTree("output0_t*.root", "ParticleInfo", {"@isParentCE", kZdcProtonAtoms[0]}, kPlot3DTag, Plot3DColumns(), kPlot3DTagDeps);
}
inline CacheTree Plot3DAllEv(){ // C8: ZDC protons of every event (cut_AllEvProtons) + columns
    return GetCutTree("output0_t*.root", "ParticleInfo", {kZdcProtonAtoms[0]}, kPlot3DTag, Plot3DColumns(), kPlot3DTagDeps);
}
inline CacheTree Plot3DRaw(){ // C9: AnglesRaw + columns
    return GetCutTree("cuts/cut_angles_all.root", "AnglesRaw", {}, kPlot3DTag, Plot3DColumns(), kPlot3DTagDeps);
}
