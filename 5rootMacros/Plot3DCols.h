#pragma once
// The Plot3D maps, the columns that replace the per-entry isReco/rawVal/rawOk/thetaOk calls of the TTree::Draw cuts (P3), and where each tree of Plot3D gets them:
// the ZDC proton cut files cut_ZdcCE.root / cut_ZdcAllEv.root hold them as columns (CutPlots.cpp (zdc)), AnglesRaw gets them as TTree aliases of its own branches (Plot3DRawAliases).
// Shared by Plot3D.cpp and CutPlots.cpp.
#include <TFile.h>
#include <TError.h>
#include <TMath.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadExecutor.hxx>

#include <algorithm>
#include <functional>
#include <iterator>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <utility>

#include "FitCore.h"
#include "QualityCut.h"

using ColumnAdder = std::function<ROOT::RDF::RNode(ROOT::RDF::RNode)>; // adds the Plot3D columns to a dataframe node (Define with compiled lambdas)

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

// Fills the maps above from cuts/cut_angles_all.root and the block limits from settings.csv. False (with an Error) if the file or tree is missing.
// one fitted proton of AnglesRaw while LoadPlot3DMaps builds the maps: its key, which sets it belongs to (bits: reco, reco block, res2) and its 17 values
struct P3Item {std::pair<int, int> key; unsigned char flags; std::array<double, 17> a;};
inline bool LoadPlot3DMaps(){
    TFile* fReco = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree* raw = fReco ? (TTree*)fReco->Get("AnglesRaw") : nullptr;
    if (!raw){Error("Plot3D", "missing AnglesRaw tree in cut_angles_all.root (run ReconAngleAll.cpp first)"); return false;}
    // one pass over AnglesRaw, each entry read once and on all the threads (RDataFrame): every row is a reconstructed track (the Angles tree), gRawVals has every row.
    // Every slot collects its rows in a vector, the vectors are sorted and merged on all the threads, and the maps are built in key order (insertion at the end, no tree search), one map per thread.
    const bool withRaw = raw->GetBranch("nProtonsEv") != nullptr;
    const unsigned char kReco = 1, kRecoD = 2, kRes2 = 4;
    ROOT::RDataFrame df(*raw);
    std::vector<std::vector<P3Item>> parts(df.GetNSlots());
    auto flagsOf = [&](double resDRms){
        return (unsigned char)(kReco | (PassesBlockResidual(resDRms) ? kRecoD : 0) | (resDRms < kRes2Mm ? kRes2 : 0));
    };
    if (withRaw){
        const double rad = TMath::RadToDeg(), twoPi = 2 * TMath::Pi();
        df.ForeachSlot([&](unsigned slot, int ev, int tid, double resRms, double resDRms, double impact, double impactD, double r2, double r2D, int nP, int nH,
                           double phiRad, double phiBlockRad, double firstXD, double firstYD, double zRecon, double resHalfRatioD, double thetaRad, double thetaBlockRad){
            P3Item it; it.key = {ev, tid}; it.flags = flagsOf(resDRms);
            double* v = it.a.data();
            v[0] = resRms; v[1] = resDRms; v[2] = impact; v[3] = impactD; v[4] = r2; v[5] = r2D;
            v[6] = nP; v[7] = 1; // slot 7 is unused (it held zOrderOK), kept so the other indexes stay put
            v[8] = resHalfRatioD; v[9] = thetaRad * rad; v[10] = nH; v[11] = thetaBlockRad * rad; v[12] = firstXD; v[13] = firstYD;
            v[14] = (phiRad < 0 ? phiRad + twoPi : phiRad) * rad; v[15] = (phiBlockRad < 0 ? phiBlockRad + twoPi : phiBlockRad) * rad; v[16] = zRecon;
            parts[slot].push_back(it);
        }, {"fEvent", "TrackID", "resRms", "resDRms", "impact", "impactD", "r2", "r2D", "nProtonsEv", "nHits",
            "outPhi", "outDPhi", "firstXD", "firstYD", "z_bx_reconVertexOut", "resHalfRatioD", "outTheta", "outDTheta"});
    } else {
        df.ForeachSlot([&](unsigned slot, int ev, int tid, double resDRms){
            P3Item it; it.key = {ev, tid}; it.flags = flagsOf(resDRms); it.a.fill(0); it.a[1] = resDRms;
            parts[slot].push_back(it);
        }, {"fEvent", "TrackID", "resDRms"});
    }
    fReco->Close();

    ROOT::TThreadExecutor pool;
    const auto byKey = [](const P3Item& x, const P3Item& y){return x.key < y.key;};
    const unsigned n = (unsigned)parts.size();
    pool.Foreach([&](unsigned s){std::sort(parts[s].begin(), parts[s].end(), byKey);}, ROOT::TSeqU(n));
    for (unsigned step = 1; step < n; step *= 2){
        pool.Foreach([&](unsigned k){
            const unsigned i = k * 2 * step;
            if (i + step >= n) return;
            std::vector<P3Item> merged;
            merged.reserve(parts[i].size() + parts[i + step].size());
            std::merge(parts[i].begin(), parts[i].end(), parts[i + step].begin(), parts[i + step].end(), std::back_inserter(merged), byKey);
            parts[i].swap(merged);
            std::vector<P3Item>().swap(parts[i + step]);
        }, ROOT::TSeqU((n + 2 * step - 1) / (2 * step)));
    }
    const std::vector<P3Item>& all = parts[0];
    pool.Foreach([&](unsigned which){
        if (which == 0){ if (withRaw) for (const P3Item& it : all) gRawVals.emplace_hint(gRawVals.end(), it.key, it.a); return; }
        std::set<std::pair<int, int>>& target = which == 1 ? gRecoTracks : which == 2 ? gRecoTracksD : gRes2Tracks;
        const unsigned char bit = which == 1 ? kReco : which == 2 ? kRecoD : kRes2;
        for (const P3Item& it : all) if (it.flags & bit) target.emplace_hint(target.end(), it.key);
    }, ROOT::TSeqU(4u));
    std::vector<P3Item>().swap(parts[0]);
    return true;
}

// One map lookup per row and thread: the 40 columns below are defined in a row and all ask for the same (fEvent, TrackID).
inline const std::array<double, 17>* P3Row(int ev, int tid){
    thread_local int lastEv = 0, lastTid = 0; thread_local const std::array<double, 17>* lastRow = nullptr; thread_local bool valid = false;
    if (valid && lastEv == ev && lastTid == tid) return lastRow;
    auto it = gRawVals.find({ev, tid});
    lastEv = ev; lastTid = tid; lastRow = it == gRawVals.end() ? nullptr : &it->second; valid = true;
    return lastRow;
}

// Adds isReco, isRecoD, isRes2 (Int_t 0/1), raw_0..raw_16 (Double_t, NaN if the key is missing) and ok_0..ok_16 (Int_t, the rawOk logic).
// The maps are read only while the event loop runs, so the lambdas are safe on several threads.
inline ColumnAdder Plot3DColumns(){
    return [](ROOT::RDF::RNode n){
        static const bool mapsLoaded = LoadPlot3DMaps(); // the 213M rows of AnglesRaw are read once, when CutPlots.cpp (zdc) builds the columns
        if (!mapsLoaded) throw std::runtime_error("Plot3DColumns: the maps could not be loaded (cuts/cut_angles_all.root)");
        n = n.Define("isReco", [](int ev, int tid){return (Int_t)isReco(ev, tid);}, {"fEvent", "TrackID"})
             .Define("isRecoD", [](int ev, int tid){return (Int_t)isRecoD(ev, tid);}, {"fEvent", "TrackID"})
             .Define("isRes2", [](int ev, int tid){return (Int_t)isRes2(ev, tid);}, {"fEvent", "TrackID"});
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

// Rewrites the isReco/isRecoD/isRes2/rawVal/rawOk calls of a TTree::Draw expression or cut into the columns (or, on AnglesRaw, the aliases).
inline std::string P3(const std::string& s){
    static const std::regex fn(R"(\b(isRecoD|isReco|isRes2)\s*\(\s*fEvent\s*,\s*TrackID\s*\))");
    static const std::regex val(R"(\brawVal\s*\(\s*fEvent\s*,\s*TrackID\s*,\s*(\d+)\s*\))");
    static const std::regex ok(R"(\brawOk\s*\(\s*fEvent\s*,\s*TrackID\s*,\s*(\d+)\s*\))");
    std::string o = std::regex_replace(s, fn, "($1==1)");
    o = std::regex_replace(o, val, "raw_$1");
    o = std::regex_replace(o, ok, "(ok_$1==1)");
    return o;
}

// The same names on AnglesRaw, as TTree aliases of its own branches (no extra columns in the file): a row of AnglesRaw is one fitted proton, so every row is a reco proton (isReco = 1),
// raw_k its k-th value of gRawVals (LoadPlot3DMaps) and ok_k whether that value is finite. P3() rewrites the Draw strings to these names for every tree.
inline void Plot3DRawAliases(TTree* raw){
    const char* val[17] = {"resRms", "resDRms", "impact", "impactD", "r2", "r2D", "nProtonsEv", "1", "resHalfRatioD", "outTheta*TMath::RadToDeg()", "nHits",
                           "outDTheta*TMath::RadToDeg()", "firstXD", "firstYD", "outPhi*TMath::RadToDeg()", "outDPhi*TMath::RadToDeg()", "z_bx_reconVertexOut"}; // phi in 0 to 360 (ReconAngleAll.cpp)
    for (int k = 0; k < 17; k++){
        raw->SetAlias(Form("raw_%d", k), val[k]);
        raw->SetAlias(Form("ok_%d", k), Form("TMath::Finite(raw_%d)", k));
    }
    raw->SetAlias("isReco", "1");
    raw->SetAlias("isRecoD", ("(" + std::string(kBlockCut) + ")").c_str());
    raw->SetAlias("isRes2", Form("(resDRms < %g)", kRes2Mm));
}
