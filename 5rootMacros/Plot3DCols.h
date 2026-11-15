#pragma once
// The Plot3D maps, the columns that replace the per-entry isReco/rawVal/rawOk calls of the TTree::Draw cuts (P3), and where each tree of Plot3D gets them:
// the cut files cut_Plot3DCE.root / cut_Plot3DAllEv.root (rows of AnglesRaw) hold them as columns (CutPlots.cpp (zdc)).
// Shared by F1_2_Plot3D.cpp and CutPlots.cpp.
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
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

#include <unordered_map>

#include "FitCore.h"

using ColumnAdder = std::function<ROOT::RDF::RNode(ROOT::RDF::RNode)>; // adds the Plot3D columns to a dataframe node (Define with compiled lambdas)

// Every recon event of AnglesRaw (cut_angles_all.root, before any cut) by fEvent, as three flat arrays in key order (LoadPlot3DMaps), found by a binary search of gP3Keys.
// A std::map and three std::sets of the same rows took about 100 GB at 213M rows, the arrays 32 GB.
// gP3Flags: kP3Reco the event has a reconstructed fit (cut_angles_all.root, Angles): the position plots only show the ZDC protons of these events, the ones that can actually be drawn on the ZDC.
// There is no residual flag: the block position plots use the same isReco as the true position ones.
// gP3Vals: the cut candidates for the 1_additionalinfo plots of the KE cut folders. Index: 14, 15: RECONPHI_ESTNEU and RECONPHIBL_ESTNEU in degrees, 0 to 360; 16: recon vertex z (mm), from the block position line fit.
// Slots 17 and 18 are nTrackerHitsEv (tracker entries of every charged particle in the event) and nTrackerTracksEv (charged tracks hitting 2+ planes), NaN without them.
// Slot 7 is nZdcCharged, the charged particles at the ZDC in the event (EventSummary of the sim, via ReconAngleAll.cpp), NaN without it.
// Empty when AnglesRaw has no nProtonsEv column.
// Slot 19 is nMaxChargedPerPlaneEv, the most charged tracks that entered one tracker plane in the event (TrackerHitsAll, ReconAngleAll.cpp), NaN without it.
// Slots 20, 21, 22: the est neutron range of the block recon proton (RECONTHETABL_ESTNEU): thetaRelNDMin, thetaRelNDMax (deg) and phiRelNDHalf (rad), NaN without them (the error by source of the phi fits).
// Slots 23, 24, 25: the block grid range alone of that proton (thetaRelNDBMin, thetaRelNDBMax in deg, phiRelNDBHalf in rad): the slope error of the block position fit, NaN without them.
// Slot 26 is the average charged tracker plane hits per hit plane in the event, nTrackerHitsEv / nPlanesHitEv (5_avgTrackerHitsPerPlane), NaN without a hit plane or without the columns.
inline constexpr int kP3Slots = 27; // values per recon event
inline constexpr long long kP3Stride = 1024; // a row of AnglesRaw is a track of an event: its key is fEvent * kP3Stride + trackNo
inline std::vector<long long> gP3Keys;
inline std::vector<unsigned char> gP3Flags;
inline std::vector<std::array<double, kP3Slots>> gP3Vals;
inline constexpr unsigned char kP3Reco = 1;
inline long P3Index(int fEvent){ // position of the first recon track of the event in gP3Keys, -1 if the event has none (the event-level lookups: isReco, rawVal of an event)
    const auto it = std::lower_bound(gP3Keys.begin(), gP3Keys.end(), (long long)fEvent * kP3Stride);
    return (it != gP3Keys.end() && *it / kP3Stride == fEvent) ? (long)(it - gP3Keys.begin()) : -1;
}
inline long P3IndexTrack(int fEvent, int trackNo){ // position of the track of the event in gP3Keys, -1 if it is not there (the lookups of an AnglesRaw row)
    const long long key = (long long)fEvent * kP3Stride + trackNo;
    const auto it = std::lower_bound(gP3Keys.begin(), gP3Keys.end(), key);
    return (it != gP3Keys.end() && *it == key) ? (long)(it - gP3Keys.begin()) : -1;
}
inline bool isReco(int fEvent){const long i = P3Index(fEvent); return i >= 0 && (gP3Flags[i] & kP3Reco) != 0;}

// Fills the arrays above from cuts/cut_angles_all.root and the block limits from settings.csv. False (with an Error) if the file or tree is missing.
// one event of AnglesRaw while LoadPlot3DMaps builds the arrays: its key, which sets it belongs to (bit: reco) and its kP3Slots values
struct P3Item {long long key; unsigned char flags; std::array<double, kP3Slots> a;};
inline bool LoadPlot3DMaps(){
    TFile* fReco = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree* raw = fReco ? (TTree*)fReco->Get("AnglesRaw") : nullptr;
    if (!raw){Error("Plot3D", "missing AnglesRaw tree in cut_angles_all.root (run ReconAngleAll.cpp first)"); return false;}
    // one pass over AnglesRaw, each entry read once and on all the threads (RDataFrame): every row is a reconstructed event (the Angles tree), gP3Vals has every row.
    // Every slot collects its rows in a vector, the vectors are sorted and merged on all the threads, and the arrays are filled from the merged vector in key order.
    const bool withRaw = raw->GetBranch("nProtonsEv") != nullptr;
    const unsigned char kReco = kP3Reco;
    ROOT::RDataFrame df0(*raw);
    ROOT::RDF::RNode df = df0;
    for (const char* col : {"thetaRelNDMin", "thetaRelNDMax", "phiRelNDHalf", "thetaRelNDBMin", "thetaRelNDBMax", "phiRelNDBHalf"}){
        if (!raw->GetBranch(col)) df = df.Define(col, []{return std::numeric_limits<double>::quiet_NaN();}); // an AnglesRaw file made before the range columns
    }
    for (const char* col : {"nZdcCharged", "nTrackerHitsEv", "nTrackerTracksEv", "nMaxChargedPerPlaneEv", "nPlanesHitEv"}){
        if (!raw->GetBranch(col)) df = df.Define(col, []{return (int)-1;}); // an AnglesRaw file made before this column
    }
    if (!raw->GetBranch("trackNo")) df = df.Define("trackNo", []{return (int)0;}); // an AnglesRaw file made before the row per track
    std::vector<std::vector<P3Item>> parts(df.GetNSlots());
    if (withRaw){
        const double rad = TMath::RadToDeg(), twoPi = 2 * TMath::Pi();
        df.ForeachSlot([&](unsigned slot, int ev, double resRms, double resDRms, double impact, double impactD, double r2, double r2D, int nP, int nZdc, int nH,
                           double phiRad, double phiBlockRad, double firstXD, double firstYD, double zRecon, double resHalfRatioD, double thetaRad, double thetaBlockRad, int nTH, int nTT, int nMP, double thMinDRad, double thMaxDRad, double phHalfD, double thMinBRad, double thMaxBRad, double phHalfB, int nPH, int trackNo){
            P3Item it; it.key = (long long)ev * kP3Stride + trackNo; it.flags = kReco;
            double* v = it.a.data();
            v[0] = resRms; v[1] = resDRms; v[2] = impact; v[3] = impactD; v[4] = r2; v[5] = r2D;
            v[6] = nP; v[7] = nZdc >= 0 ? nZdc : NAN;
            v[8] = resHalfRatioD; v[9] = thetaRad * rad; v[10] = nH; v[11] = thetaBlockRad * rad; v[12] = firstXD; v[13] = firstYD;
            v[14] = (phiRad < 0 ? phiRad + twoPi : phiRad) * rad; v[15] = (phiBlockRad < 0 ? phiBlockRad + twoPi : phiBlockRad) * rad; v[16] = zRecon; v[17] = nTH >= 0 ? nTH : NAN; v[18] = nTT >= 0 ? nTT : NAN; v[19] = nMP >= 0 ? nMP : NAN;
            v[20] = thMinDRad * rad; v[21] = thMaxDRad * rad; v[22] = phHalfD;
            v[23] = thMinBRad * rad; v[24] = thMaxBRad * rad; v[25] = phHalfB;
            v[26] = (nPH > 0 && nTH >= 0) ? double(nTH) / nPH : NAN;
            parts[slot].push_back(it);
        }, {"fEvent", "resRms", "resDRms", "impact", "impactD", "r2", "r2D", "nProtonsEv", "nZdcCharged", "nHits",
            "phiRelN", "phiRelND", "firstXD", "firstYD", "z_bx_reconVertexOut", "resHalfRatioD", "thetaRelN", "thetaRelND", "nTrackerHitsEv", "nTrackerTracksEv", "nMaxChargedPerPlaneEv", "thetaRelNDMin", "thetaRelNDMax", "phiRelNDHalf", "thetaRelNDBMin", "thetaRelNDBMax", "phiRelNDBHalf", "nPlanesHitEv", "trackNo"});
    } else {
        df.ForeachSlot([&](unsigned slot, int ev, double resDRms){
            P3Item it; it.key = (long long)ev * kP3Stride; it.flags = kReco; it.a.fill(0); it.a[1] = resDRms;
            parts[slot].push_back(it);
        }, {"fEvent", "resDRms"});
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
    gP3Keys.clear(); gP3Flags.clear(); gP3Vals.clear();
    gP3Keys.reserve(all.size()); gP3Flags.reserve(all.size());
    if (withRaw) gP3Vals.reserve(all.size());
    for (const P3Item& it : all){
        if (!gP3Keys.empty() && gP3Keys.back() == it.key){gP3Flags.back() |= it.flags; continue;} // a repeated key is one entry holding the flags of all its rows, as the sets did
        gP3Keys.push_back(it.key); gP3Flags.push_back(it.flags);
        if (withRaw) gP3Vals.push_back(it.a);
    }
    std::vector<P3Item>().swap(parts[0]);
    return true;
}

// One lookup per row and thread: the 56 columns below are defined in a row and all ask for the same fEvent (and trackNo; -1 for a row that is not a track, the first track of the event).
inline const std::array<double, kP3Slots>* P3Row(int ev, int trackNo){
    thread_local int lastEv = 0, lastTrack = 0; thread_local const std::array<double, kP3Slots>* lastRow = nullptr; thread_local bool valid = false;
    if (valid && lastEv == ev && lastTrack == trackNo) return lastRow;
    const long i = trackNo < 0 ? P3Index(ev) : P3IndexTrack(ev, trackNo);
    lastEv = ev; lastTrack = trackNo; lastRow = (i < 0 || gP3Vals.empty()) ? nullptr : &gP3Vals[i]; valid = true;
    return lastRow;
}

// Adds isReco (Int_t 0/1), raw_0..raw_26 (Double_t, NaN if the key is missing) and ok_0..ok_26 (Int_t, 1 if the value is finite).
// The arrays are read only while the event loop runs, so the lambdas are safe on several threads.
inline ColumnAdder Plot3DColumns(){
    return [](ROOT::RDF::RNode n){
        static const bool mapsLoaded = LoadPlot3DMaps(); // the 213M rows of AnglesRaw are read once, when CutPlots.cpp (zdc) builds the columns
        if (!mapsLoaded) throw std::runtime_error("Plot3DColumns: the maps could not be loaded (cuts/cut_angles_all.root)");
        // a row of AnglesRaw (it has trackNo) takes the values of its own track; any other tree (the ZDC hits of an event) the first track of its event
        const auto cols = n.GetColumnNames();
        const bool ownTrack = std::find(cols.begin(), cols.end(), "trackNo") != cols.end();
        n = n.Define("p3Track", [ownTrack](int tn){return ownTrack ? tn : -1;}, {ownTrack ? "trackNo" : "fEvent"});
        n = n.Define("isReco", [](int ev){return (Int_t)isReco(ev);}, {"fEvent"});
        for (int k = 0; k < kP3Slots; k++){
            n = n.Define("raw_" + std::to_string(k), [k](int ev, int tn){
                const auto* r = P3Row(ev, tn);
                return r ? (*r)[k] : std::numeric_limits<double>::quiet_NaN();
            }, {"fEvent", "p3Track"});
            n = n.Define("ok_" + std::to_string(k), [k](int ev, int tn){
                const auto* r = P3Row(ev, tn);
                return (Int_t)(r && std::isfinite((*r)[k]));
            }, {"fEvent", "p3Track"});
        }
        return n;
    };
}

// Rewrites the isReco/rawVal/rawOk calls of a TTree::Draw expression or cut into the columns.
inline std::string P3(const std::string& s){
    thread_local std::unordered_map<std::string, std::string> cache; // the same expressions and cuts come back in many calls
    if (const auto hit = cache.find(s); hit != cache.end()) return hit->second;
    static const std::regex fn(R"(\b(isReco)\s*\(\s*fEvent\s*\))");
    static const std::regex val(R"(\brawVal\s*\(\s*fEvent\s*,\s*(\d+)\s*\))");
    static const std::regex ok(R"(\brawOk\s*\(\s*fEvent\s*,\s*(\d+)\s*\))");
    std::string o = std::regex_replace(s, fn, "($1==1)");
    o = std::regex_replace(o, val, "raw_$1");
    o = std::regex_replace(o, ok, "(ok_$1==1)");
    cache.emplace(s, o);
    return o;
}
