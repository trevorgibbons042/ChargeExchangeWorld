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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <unordered_map>

#include "FitCore.h"
#include "QualityCut.h"

using ColumnAdder = std::function<ROOT::RDF::RNode(ROOT::RDF::RNode)>; // adds the Plot3D columns to a dataframe node (Define with compiled lambdas)

// Every fitted proton of AnglesRaw (cut_angles_all.root, before any cut) by (fEvent, TrackID), as three flat arrays in key order (LoadPlot3DMaps), found by a binary search of gP3Keys.
// A std::map and three std::sets of the same rows took about 100 GB at 213M rows, the arrays 32 GB.
// gP3Flags: kP3Reco the proton has a reconstructed track (cut_angles_all.root, Angles): the position plots only show these protons, the ones that can actually be drawn on the ZDC.
// There is no residual flag: the block position plots use the same isReco as the true position ones.
// gP3Vals: the cut candidates for the 1_additionalinfo plots of the KE cut folders. Index: see kInfoVars in Plot3D (14, 15: RECONPHI_ESTNEU and RECONPHIBL_ESTNEU in degrees, 0 to 360); 16: recon vertex z (mm), from the block position line fit.
// Slots 17 and 18 are nTrackerHitsEv (tracker entries of every charged particle in the event) and nTrackerTracksEv (charged tracks hitting 2+ planes), NaN without them.
// Slot 7 is nZdcCharged, the charged particles at the ZDC in the event (EventSummary of the sim, via ReconAngleAll.cpp), NaN without it.
// Empty when AnglesRaw has no nProtonsEv column. rawOk is the Draw cut that keeps the protons that have a finite value.
// Slot 19 is nMaxChargedPerPlaneEv, the most charged tracks that entered one tracker plane in the event (ParticleInfo, ReconAngleAll.cpp), NaN without it.
inline constexpr int kP3Slots = 20; // values per fitted proton
inline std::vector<std::pair<int, int>> gP3Keys;
inline std::vector<unsigned char> gP3Flags;
inline std::vector<std::array<double, kP3Slots>> gP3Vals;
inline constexpr unsigned char kP3Reco = 1;
inline long P3Index(int fEvent, int trackID){ // position of the key in gP3Keys, -1 if it is not there
    const std::pair<int, int> key{fEvent, trackID};
    const auto it = std::lower_bound(gP3Keys.begin(), gP3Keys.end(), key);
    return (it != gP3Keys.end() && *it == key) ? (long)(it - gP3Keys.begin()) : -1;
}
inline bool P3HasFlag(int fEvent, int trackID, unsigned char bit){const long i = P3Index(fEvent, trackID); return i >= 0 && (gP3Flags[i] & bit) != 0;}
inline bool isReco(int fEvent, int trackID){return P3HasFlag(fEvent, trackID, kP3Reco);}

inline bool rawOk(int fEvent, int trackID, int k){
    const long i = P3Index(fEvent, trackID);
    return i >= 0 && !gP3Vals.empty() && std::isfinite(gP3Vals[i][k]);
}
inline double rawVal(int fEvent, int trackID, int k){
    const long i = P3Index(fEvent, trackID);
    if (i < 0 || gP3Vals.empty()) throw std::out_of_range("rawVal: key not in AnglesRaw");
    return gP3Vals[i][k];
}

// Fills the arrays above from cuts/cut_angles_all.root and the block limits from settings.csv. False (with an Error) if the file or tree is missing.
// one fitted proton of AnglesRaw while LoadPlot3DMaps builds the arrays: its key, which sets it belongs to (bit: reco) and its kP3Slots values
struct P3Item {std::pair<int, int> key; unsigned char flags; std::array<double, kP3Slots> a;};
inline bool LoadPlot3DMaps(){
    TFile* fReco = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree* raw = fReco ? (TTree*)fReco->Get("AnglesRaw") : nullptr;
    if (!raw){Error("Plot3D", "missing AnglesRaw tree in cut_angles_all.root (run ReconAngleAll.cpp first)"); return false;}
    // Only the rows of AnglesRaw that a ZDC input of CutPlots.cpp (zdc) asks for are kept: the (fEvent, TrackID) keys of those inputs are read first (small files), sorted, and the AnglesRaw
    // rows are filtered on them before any other branch is read. A key that is not kept is never looked up, so the columns do not change.
    std::vector<std::pair<int, int>> want;
    for (const char* keyFile : {"cuts/cut_ParentZdc.root", "cuts/cut_AllEvProtons.root", "cuts/cut_AllEvParticles.root"}){
        ROOT::RDataFrame dfKeys("ParticleInfo", keyFile);
        std::vector<std::vector<std::pair<int, int>>> slotKeys(dfKeys.GetNSlots());
        dfKeys.ForeachSlot([&slotKeys](unsigned slot, int ev, int tid){slotKeys[slot].emplace_back(ev, tid);}, {"fEvent", "TrackID"});
        for (const auto& v : slotKeys) want.insert(want.end(), v.begin(), v.end());
    }
    std::sort(want.begin(), want.end());
    want.erase(std::unique(want.begin(), want.end()), want.end());
    // one pass over AnglesRaw, each entry read once and on all the threads (RDataFrame): every row is a reconstructed track (the Angles tree), gP3Vals has every row.
    // Every slot collects its rows in a vector, the vectors are sorted and merged on all the threads, and the arrays are filled from the merged vector in key order.
    const bool withRaw = raw->GetBranch("nProtonsEv") != nullptr;
    const unsigned char kReco = kP3Reco;
    ROOT::RDataFrame df0(*raw);
    ROOT::RDF::RNode df = df0;
    for (const char* col : {"nZdcCharged", "nTrackerHitsEv", "nTrackerTracksEv", "nMaxChargedPerPlaneEv"}){
        if (!raw->GetBranch(col)) df = df.Define(col, []{return (int)-1;}); // an AnglesRaw file made before this column
    }
    df = df.Filter([&want](int ev, int tid){return std::binary_search(want.begin(), want.end(), std::make_pair(ev, tid));}, {"fEvent", "TrackID"}); // before the ForeachSlot: its other branches are read for the kept rows only
    std::vector<std::vector<P3Item>> parts(df.GetNSlots());
    if (withRaw){
        const double rad = TMath::RadToDeg(), twoPi = 2 * TMath::Pi();
        df.ForeachSlot([&](unsigned slot, int ev, int tid, double resRms, double resDRms, double impact, double impactD, double r2, double r2D, int nP, int nZdc, int nH,
                           double phiRad, double phiBlockRad, double firstXD, double firstYD, double zRecon, double resHalfRatioD, double thetaRad, double thetaBlockRad, int nTH, int nTT, int nMP){
            P3Item it; it.key = {ev, tid}; it.flags = kReco;
            double* v = it.a.data();
            v[0] = resRms; v[1] = resDRms; v[2] = impact; v[3] = impactD; v[4] = r2; v[5] = r2D;
            v[6] = nP; v[7] = nZdc >= 0 ? nZdc : NAN;
            v[8] = resHalfRatioD; v[9] = thetaRad * rad; v[10] = nH; v[11] = thetaBlockRad * rad; v[12] = firstXD; v[13] = firstYD;
            v[14] = (phiRad < 0 ? phiRad + twoPi : phiRad) * rad; v[15] = (phiBlockRad < 0 ? phiBlockRad + twoPi : phiBlockRad) * rad; v[16] = zRecon; v[17] = nTH >= 0 ? nTH : NAN; v[18] = nTT >= 0 ? nTT : NAN; v[19] = nMP >= 0 ? nMP : NAN;
            parts[slot].push_back(it);
        }, {"fEvent", "TrackID", "resRms", "resDRms", "impact", "impactD", "r2", "r2D", "nProtonsEv", "nZdcCharged", "nHits",
            "phiRelN", "phiRelND", "firstXD", "firstYD", "z_bx_reconVertexOut", "resHalfRatioD", "thetaRelN", "thetaRelND", "nTrackerHitsEv", "nTrackerTracksEv", "nMaxChargedPerPlaneEv"});
    } else {
        df.ForeachSlot([&](unsigned slot, int ev, int tid, double resDRms){
            P3Item it; it.key = {ev, tid}; it.flags = kReco; it.a.fill(0); it.a[1] = resDRms;
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

// One lookup per row and thread: the 40 columns below are defined in a row and all ask for the same (fEvent, TrackID).
inline const std::array<double, kP3Slots>* P3Row(int ev, int tid){
    thread_local int lastEv = 0, lastTid = 0; thread_local const std::array<double, kP3Slots>* lastRow = nullptr; thread_local bool valid = false;
    if (valid && lastEv == ev && lastTid == tid) return lastRow;
    const long i = P3Index(ev, tid);
    lastEv = ev; lastTid = tid; lastRow = (i < 0 || gP3Vals.empty()) ? nullptr : &gP3Vals[i]; valid = true;
    return lastRow;
}

// Adds isReco (Int_t 0/1), raw_0..raw_16 (Double_t, NaN if the key is missing) and ok_0..ok_16 (Int_t, the rawOk logic).
// The arrays are read only while the event loop runs, so the lambdas are safe on several threads.
inline ColumnAdder Plot3DColumns(){
    return [](ROOT::RDF::RNode n){
        static const bool mapsLoaded = LoadPlot3DMaps(); // the 213M rows of AnglesRaw are read once, when CutPlots.cpp (zdc) builds the columns
        if (!mapsLoaded) throw std::runtime_error("Plot3DColumns: the maps could not be loaded (cuts/cut_angles_all.root)");
        n = n.Define("isReco", [](int ev, int tid){return (Int_t)isReco(ev, tid);}, {"fEvent", "TrackID"});
        for (int k = 0; k < kP3Slots; k++){
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

// Rewrites the isReco/rawVal/rawOk calls of a TTree::Draw expression or cut into the columns (or, on AnglesRaw, the aliases).
inline std::string P3(const std::string& s){
    thread_local std::unordered_map<std::string, std::string> cache; // the same expressions and cuts come back in many calls
    if (const auto hit = cache.find(s); hit != cache.end()) return hit->second;
    static const std::regex fn(R"(\b(isReco)\s*\(\s*fEvent\s*,\s*TrackID\s*\))");
    static const std::regex val(R"(\brawVal\s*\(\s*fEvent\s*,\s*TrackID\s*,\s*(\d+)\s*\))");
    static const std::regex ok(R"(\brawOk\s*\(\s*fEvent\s*,\s*TrackID\s*,\s*(\d+)\s*\))");
    std::string o = std::regex_replace(s, fn, "($1==1)");
    o = std::regex_replace(o, val, "raw_$1");
    o = std::regex_replace(o, ok, "(ok_$1==1)");
    cache.emplace(s, o);
    return o;
}

// The same names on AnglesRaw, as TTree aliases of its own branches (no extra columns in the file): a row of AnglesRaw is one fitted proton, so every row is a reco proton (isReco = 1),
// raw_k its k-th value of gP3Vals (LoadPlot3DMaps) and ok_k whether that value is finite. P3() rewrites the Draw strings to these names for every tree.
inline void Plot3DRawAliases(TTree* raw){
    const char* val[kP3Slots] = {"resRms", "resDRms", "impact", "impactD", "r2", "r2D", "nProtonsEv", "nZdcCharged", "resHalfRatioD", "thetaRelN*TMath::RadToDeg()", "nHits",
                           "thetaRelND*TMath::RadToDeg()", "firstXD", "firstYD", "phiRelN*TMath::RadToDeg()", "phiRelND*TMath::RadToDeg()", "z_bx_reconVertexOut", "nTrackerHitsEv", "nTrackerTracksEv", "nMaxChargedPerPlaneEv"}; // phi in 0 to 360 (ReconAngleAll.cpp)
    for (int k = 0; k < kP3Slots; k++){
        raw->SetAlias(Form("raw_%d", k), val[k]);
        raw->SetAlias(Form("ok_%d", k), Form("TMath::Finite(raw_%d)", k));
    }
    raw->SetAlias("isReco", "1");
}

// The start-theta scan of Plot3D.cpp part 4 in one pass over a sample instead of one pass per start value. th, phi (deg) are the recon (block) theta and phi of the protons of a sample, apScale the
// APScale of settings.csv. For the start values b * 0.01 deg (b 0 to 99) it returns the phi histogram of the protons with start <= theta <= 1 deg (kScanBins bins: 0 and kScanBins - 1 are the under
// and overflow of a TH1D(50, 0, 360), the others its bins) and the sum of their model A (AEffOfTheta, the one at their own theta), and aMeanAll, the mean A over every proton of the sample.
// Row ranges run on all threads (MapRows, FitCore.h) and are added in range order; the bin counts are exact, the A sums can differ from a serial sum in the last digits.
inline constexpr int kScanStarts = 100, kScanBins = 52;
// fit[b]: the unbinned P (UnbinnedPFromCos, FitCore.h) of the protons of start b, A fixed to their mean model A (aSum[b] / count), phi normalised over loDeg to hiDeg; ok is false below 10 protons.
struct StartScan {std::vector<std::array<double, kScanBins>> phi = std::vector<std::array<double, kScanBins>>(kScanStarts); std::vector<double> aSum = std::vector<double>(kScanStarts, 0.0), aSum2 = std::vector<double>(kScanStarts, 0.0); AMean aAll; std::vector<PFit> fit = std::vector<PFit>(kScanStarts);};
inline StartScan ScanStartTheta(const std::vector<double>& th, const std::vector<double>& ph, double apScale, double loDeg = 0, double hiDeg = 360, double thHi = 1.0){ // thHi: the upper theta of the protons that count (deg), 1 unless the folder has no upper limit
    const auto parts = MapRows((Long64_t)th.size(), [&](Long64_t lo, Long64_t hi){
        StartScan p;
        for (auto& a : p.phi) a.fill(0.0);
        for (Long64_t i = lo; i < hi; i++){
            const double a = AEffOfTheta(th[i], apScale);
            p.aAll.Add(a);
            if (th[i] < 0 || th[i] > thHi) continue;
            int b = std::min(kScanStarts - 1, (int)std::floor(th[i] / 0.01)); // the last start value with b * 0.01 <= theta, as the comparison of the plain loop
            while (b < kScanStarts - 1 && (b + 1) * 0.01 <= th[i]) b++;
            while (b > 0 && b * 0.01 > th[i]) b--;
            const double x = ph[i];
            const int bin = x < 0 ? 0 : !(x < 360.0) ? kScanBins - 1 : 1 + (int)(50.0 * x / 360.0); // TH1::FindBin of a fixed-width axis
            p.phi[b][bin]++; p.aSum[b] += a; p.aSum2[b] += a * a;
        }
        return p;
    });
    StartScan s;
    for (auto& a : s.phi) a.fill(0.0);
    for (const StartScan& p : parts){
        for (int b = 0; b < kScanStarts; b++){
            for (int j = 0; j < kScanBins; j++) s.phi[b][j] += p.phi[b][j];
            s.aSum[b] += p.aSum[b]; s.aSum2[b] += p.aSum2[b];
        }
        s.aAll += p.aAll;
    }
    for (int b = kScanStarts - 2; b >= 0; b--){ // start b takes every proton of the buckets b and above
        for (int j = 0; j < kScanBins; j++) s.phi[b][j] += s.phi[b + 1][j];
        s.aSum[b] += s.aSum[b + 1]; s.aSum2[b] += s.aSum2[b + 1];
    }
    { // the unbinned P of every start: the protons of the buckets b and above (th 0 to 1 deg), sorted by bucket so that start b is one stretch of the cos(phi) array
        const auto bucketOf = [](double t){
            int b = std::min(kScanStarts - 1, (int)std::floor(t / 0.01)); // the same bucket as the loop above
            while (b < kScanStarts - 1 && (b + 1) * 0.01 <= t) b++;
            while (b > 0 && b * 0.01 > t) b--;
            return b;
        };
        std::vector<Long64_t> off(kScanStarts + 1, 0);
        std::vector<int> bk(th.size(), -1);
        for (size_t i = 0; i < th.size(); i++) if (th[i] >= 0 && th[i] <= 1.0){bk[i] = bucketOf(th[i]); off[bk[i] + 1]++;}
        for (int b = 0; b < kScanStarts; b++) off[b + 1] += off[b];
        std::vector<double> cosSorted(off[kScanStarts]);
        std::vector<Long64_t> pos(off.begin(), off.end() - 1);
        for (size_t i = 0; i < th.size(); i++) if (bk[i] >= 0) cosSorted[pos[bk[i]]++] = std::cos(ph[i] * TMath::Pi() / 180);
        std::vector<int> starts(kScanStarts);
        std::iota(starts.begin(), starts.end(), 0);
        ROOT::TThreadExecutor ex(NThreads());
        s.fit = ex.Map([&](int b){
            const Long64_t n = off[kScanStarts] - off[b];
            return n >= 10 ? UnbinnedPFromCos(cosSorted.data() + off[b], n, nullptr, AMean{(double)n, s.aSum[b], s.aSum2[b]}, loDeg, hiDeg, false) : PFit();
        }, starts);
    }
    return s;
}
