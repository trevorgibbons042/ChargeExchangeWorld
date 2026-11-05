#include <TFile.h>
#include <ROOT/TBufferMerger.hxx>
#include <TChain.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "QualityCut.h"
#include "CutOpen.h"
#include "FitCore.h"
#include <chrono>

// One proton point as the read pass sees it. The read keeps these flat per slot; the tracks are built from them later, in parallel (a std::map per slot made the read as slow as the fits).
struct Point {int eventID, trackID, pdg, parentID; double x, y, z, xD, yD, zD; ULong64_t entry;};
// The ZDC row of a proton (boundaryType 2), the only place a KE is taken from: key is the event in the high 32 bits and the trackID in the low 32 bits, KE the true KE and KED the smeared one.
struct ZdcHit {ULong64_t key; double KE, KED;};

struct FitResult {double x0, y0, bx, by;
    double x0D, y0D, bxD, byD;};
// The points of one track as a struct of arrays: point i is (x[i], y[i], z[i], ...), one array per quantity.
// entry: row number in ParticleInfo (rows of a track are written in the order it hit them).
struct Track {
    std::vector<double> x, y, z, xD, yD, zD;
    std::vector<int> ParentID;
    std::vector<ULong64_t> entry;
    std::size_t size() const {return z.size();}
    void clear(){ // empties the track and keeps the capacity, so one Track serves every track of a shard without new allocations
        x.clear(); y.clear(); z.clear(); xD.clear(); yD.clear(); zD.clear();
        ParentID.clear(); entry.clear();
    }
    void push(double x_, double y_, double z_, double xD_, double yD_, double zD_, int parent, ULong64_t entry_){
        x.push_back(x_); y.push_back(y_); z.push_back(z_); xD.push_back(xD_); yD.push_back(yD_); zD.push_back(zD_);
        ParentID.push_back(parent); entry.push_back(entry_);
    }
    void append(const Track& o){
        for (std::size_t i = 0; i < o.size(); i++) push(o.x[i], o.y[i], o.z[i], o.xD[i], o.yD[i], o.zD[i], o.ParentID[i], o.entry[i]);
    }
};

// Least-squares line x = x0 + bx z, y = y0 + by z through the points of t: all of them, or with block the detector block points (zD != 0) at their block positions.
// Leaves the four outputs alone (NAN from the caller) with fewer than 2 points or no spread in z.
static void LineFit(const Track& t, bool block, double& x0, double& y0, double& bx, double& by){
    const std::vector<double>& X = block ? t.xD : t.x; const std::vector<double>& Y = block ? t.yD : t.y; const std::vector<double>& Z = block ? t.zD : t.z;
    double meanX = 0.0, meanY = 0.0, meanZ = 0.0;
    std::size_t N = 0;
    for (std::size_t i = 0; i < t.size(); i++) {
        if (block && t.zD[i] == 0.0) {continue;}
        meanX += X[i]; meanY += Y[i]; meanZ += Z[i];
        N++;
    }
    if (N < 2) {return;}
    meanX /= N; meanY /= N; meanZ /= N;

    double Szz = 0.0, Szx = 0.0, Szy = 0.0;
    for (std::size_t i = 0; i < t.size(); i++) {
        if (block && t.zD[i] == 0.0) {continue;}
        const double dz = Z[i] - meanZ;
        Szz += dz * dz;
        Szx += dz * (X[i] - meanX);
        Szy += dz * (Y[i] - meanY);
    }
    if (std::abs(Szz) >= 1e-20) {
        bx = Szx / Szz; by = Szy / Szz;
        x0 = meanX - bx * meanZ; y0 = meanY - by * meanZ;
    }
}

FitResult fitTrack(const Track& t){
    FitResult fit{NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN};
    LineFit(t, false, fit.x0, fit.y0, fit.bx, fit.by);
    LineFit(t, true, fit.x0D, fit.y0D, fit.bxD, fit.byD);
    return fit;
}

// RMS transverse distance (mm) of the points to the fitted line x = x0 + bx z, y = y0 + by z, and the R^2 of the same line fit, x and y together:
// 1 - sum of squared residuals / sum of squared deviations from the mean hit position. block: only the detector block points (zD != 0), at their block positions.
// rms is NAN with no points; r2 is NAN with fewer than 3 points or no spread. One function for both: the mean position and the residuals are one pass each.
static void FitQuality(const Track& t, double x0, double y0, double bx, double by, bool block, double& rms, double& r2){
    double sx = 0, sy = 0; int n = 0;
    for (std::size_t i = 0; i < t.size(); i++) {
        if (block && t.zD[i] == 0.0) {continue;}
        sx += block ? t.xD[i] : t.x[i]; sy += block ? t.yD[i] : t.y[i]; n++;
    }
    if (n == 0) {rms = NAN; r2 = NAN; return;}
    const double mx = sx / n, my = sy / n;
    double ssRes = 0, ssTot = 0;
    for (std::size_t i = 0; i < t.size(); i++) {
        if (block && t.zD[i] == 0.0) {continue;}
        const double x = block ? t.xD[i] : t.x[i], y = block ? t.yD[i] : t.y[i], z = block ? t.zD[i] : t.z[i];
        const double rx = x - (x0 + bx * z), ry = y - (y0 + by * z);
        ssRes += rx * rx + ry * ry;
        ssTot += (x - mx) * (x - mx) + (y - my) * (y - my);
    }
    rms = std::sqrt(ssRes / n);
    r2 = (n >= 3 && ssTot > 0) ? 1.0 - ssRes / ssTot : NAN;
}

// RMS transverse distance (mm) of the block hits hits[lo] to hits[hi - 1] (indices into t, all with zD != 0) to the line x = x0 + bx z, y = y0 + by z, at their block positions. NAN with no hits.
static double ResidualRmsHits(const Track& t, const std::vector<std::size_t>& hits, std::size_t lo, std::size_t hi, double x0, double y0, double bx, double by){
    double sum = 0.0;
    for (std::size_t k = lo; k < hi; k++) {
        const std::size_t i = hits[k];
        const double rx = t.xD[i] - (x0 + bx * t.zD[i]), ry = t.yD[i] - (y0 + by * t.zD[i]);
        sum += rx * rx + ry * ry;
    }
    return hi > lo ? std::sqrt(sum / (hi - lo)) : NAN;
}

// Closest approach (mm) of the line x = x0 + bx z, y = y0 + by z to the beam axis (x = y = 0).
static double ImpactToAxis(double x0, double y0, double bx, double by){
    const double slope = std::hypot(bx, by);
    return slope > 0 ? std::fabs(x0 * by - y0 * bx) / slope : std::hypot(x0, y0);
}

// One pass fits every proton and writes the tree AnglesRaw with the true (finalKE) and the smeared (finalKED) KE at the ZDC, NAN for a track that never hit the ZDC.
void ReconAngleAll(){
    using Clock = std::chrono::steady_clock;
    const auto tStart = Clock::now();
    auto secs = [](Clock::time_point a, Clock::time_point b){return std::chrono::duration<double>(b - a).count();};
    ROOT::EnableImplicitMT(NThreads());
    ROOT::TTreeProcessorMT::SetTasksPerWorkerHint(4); // a raw file is split into 4 tasks per thread, so one file keeps all the threads busy

    // CE neutrons (same selection as CutEvents.cpp): event -> CE track IDs.
    // the carbon-block CE come from cuts/cut_CE.root (CarbonCeCut, CutOpen.h): every row passes, no filter
    TTree* ceCarbon = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
    if (!ceCarbon){Error("ReconAngleAll", "no cuts/cut_CE.root (run CutEvents.cpp)"); return;}
    ROOT::RDataFrame dfCE(*ceCarbon);
    std::vector<std::map<int, std::set<int>>> slotCeTracks(dfCE.GetNSlots());
    dfCE.ForeachSlot([&](unsigned slot, int eventID, int trackID){
        slotCeTracks[slot][eventID].insert(trackID);
    }, {"fEvent", "TrackID"});
    std::map<int, std::set<int>> ceTracks;
    for (auto &m : slotCeTracks) {
        for (auto &[ev, ids] : m) ceTracks[ev].insert(ids.begin(), ids.end());
    }

    // Per-event counts of the sim (EventSummary: every event, no reconstruction needed): event -> charged particles at the ZDC, tracker entries of all charged particles, charged tracks that hit 2+ tracker planes.
    // -1 for an event without a row, a sim without the tree or without that column.
    struct EvCounts {int zdcCharged = -1, trackerHits = -1, trackerTracks = -1;};
    std::vector<EvCounts> evCountsOf;
    {
        TChain es("EventSummary");
        es.Add("output0_t*.root");
        if (es.LoadTree(0) >= 0){
            ROOT::RDataFrame dfEsAll(es);
            ROOT::RDF::RNode dfEs = dfEsAll;
            for (const char* col : {"nZdcCharged", "nTrackerHitsEv", "nTrackerTracksEv"}){
                if (!es.GetBranch(col)) dfEs = dfEs.Define(col, []{return (int)-1;}); // a sim run before that column
            }
            std::vector<std::vector<std::pair<int, EvCounts>>> slotRows(dfEsAll.GetNSlots());
            dfEs.ForeachSlot([&](unsigned slot, int eventID, int nZ, int nH, int nT){slotRows[slot].push_back({eventID, EvCounts{nZ, nH, nT}});},
                             {"fEvent", "nZdcCharged", "nTrackerHitsEv", "nTrackerTracksEv"});
            int maxEvent = -1;
            for (const auto& rows : slotRows) for (const auto& r : rows) maxEvent = std::max(maxEvent, r.first);
            evCountsOf.assign(maxEvent + 1, EvCounts());
            for (const auto& rows : slotRows) for (const auto& r : rows) evCountsOf[r.first] = r.second;
        }
    }

    // One proton tree: every recon proton before the quality cut (QualityStudy.cpp). The old trees Angles, AnglesCE and their _DE copies are views of it
    // (AngleSelection, CutFiles.h): inCE picks the rows, finalKE / finalKED the true / smeared ZDC KE.
    // The file is written through a TBufferMerger: the shards fill their own trees on all the threads (see the write below).
    ROOT::TBufferMerger merger("cuts/cut_angles_all.root", "RECREATE");

    // A group of raw files at a time (ANA_FILES_PER_PASS, default 8): an event, and so every track, lives in a single output0_t<thread>.root, so the tracks of a group are complete
    // and only that group's proton points are held in memory (all files at once does not fit at 100M events). One file alone cannot keep all the threads busy.
    TChain rawFiles("ParticleInfo");
    rawFiles.Add("output0_t*.root");
    const unsigned filesPerPass = std::getenv("ANA_FILES_PER_PASS") ? std::max(1, std::atoi(std::getenv("ANA_FILES_PER_PASS"))) : 8;
    std::vector<std::string> rawNames;
    for (const TObject* rawFile : *rawFiles.GetListOfFiles()) rawNames.push_back(rawFile->GetTitle());
    const unsigned nShards = NThreads();
    for (size_t first = 0; first < rawNames.size(); first += filesPerPass) {
        TChain group("ParticleInfo");
        for (size_t i = first; i < std::min(rawNames.size(), first + filesPerPass); i++) group.Add(rawNames[i].c_str());
        ROOT::RDataFrame df0(group);
        auto df = WithEntryNumber(df0, df0.GetNSlots()); // entryNo: the row of ParticleInfo in its file (rdfentry_ is wrong under MT on a tree with several clusters); orders the points of a track

        const auto tRead0 = Clock::now();
        // The read puts every point straight into the bucket of its shard (event % nShards) of its slot: all the tracks of an event land in one shard, so the per-event
        // proton count is complete per shard, and a shard reads only its own buckets (a shard that looked through every point of every slot would scan all of them, nShards times over).
        std::vector<std::vector<std::vector<Point>>> buckets(df.GetNSlots(), std::vector<std::vector<Point>>(nShards));
        // boundaryType 1 (tracker exit, older sim runs only) is left out. A ZDC row (boundaryType 2) is not a fit point: it goes to zdcBuckets, the only source of a KE (the trackers give none),
        // so a tracker hit is its entry row (boundaryType 0) only; particles that never reach the ZDC are reconstructed too, with KE NAN.
        // The sim writes only charged particles (FillParticleInfo in 6detector.hh); the pdg test below is a guard so a neutral row never becomes a recon particle.
        // A nucleus (pdg >= 1000000000) is charged when its Z is above 0, an unknown pdg is kept. The PDG table is loaded once here so the threads only read it.
        TDatabasePDG::Instance()->GetParticle(2212);
        const auto isCharged = [](int pdg){
            const int a = pdg < 0 ? -pdg : pdg;
            if (a >= 1000000000) return (a / 10000) % 1000 > 0;
            const TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(pdg);
            return p == nullptr || p->Charge() != 0.;
        };
        // The selection is a typed Filter lambda, not a Filter string: a string Filter is JIT-compiled, which fails on a ROOT whose cling cannot materialize libc++ symbols (the Mac).
        // A Filter node, unlike a test inside the ForeachSlot lambda, reads the other branches only for rows that pass, not for every particle row.
        std::vector<std::vector<std::vector<ZdcHit>>> zdcBuckets(df.GetNSlots(), std::vector<std::vector<ZdcHit>>(nShards));
        df.Filter([isCharged](int boundaryType, int pdg){return boundaryType != 1 && isCharged(pdg);}, {"boundaryType", "pdg"})
          .ForeachSlot([&](unsigned slot,
            int boundaryType, int eventID, int trackID, int pdg,
            double x, double y, double z, double xD, double yD, double zD, double KE, double KED, int ParentID, ULong64_t entry){
            const unsigned shard = (unsigned)eventID % nShards;
            if (boundaryType == 2) {zdcBuckets[slot][shard].push_back({(ULong64_t(std::uint32_t(eventID)) << 32) | std::uint32_t(trackID), KE, KED}); return;}
            buckets[slot][shard].push_back({eventID, trackID, pdg, ParentID, x, y, z, xD, yD, zD, entry});
        }, {"boundaryType", "fEvent", "TrackID", "pdg", "fX", "fY", "fZ", "fDx", "fDy", "fDz", "fKEnergy", "fDEnergy", "ParentID", "entryNo"});
        // The fits run on every thread: each shard builds its tracks, fits them and keeps the results as rows. Only the TTree fills below are serial.
        struct Row {
            int trackID, pdg, eventID, parentID, isFromCE, inCE, nProtonsEv, nZdcCharged, nTrackerHitsEv, nTrackerTracksEv, nHits, nMaxChargedPerPlaneEv;
            double theta, phi, dTheta, dPhi, x0, y0, bx, by, x0D, y0D, bxD, byD, finalKE, finalKED;
            double xVertex, yVertex, zVertex, zBxRecon, zByRecon, zBxReconDIFF, zByReconDIFF;
            double resRms, resDRms, impact, impactD, r2, r2D, firstXD, firstYD, resHalfRatioD;
        };
        std::vector<std::vector<Row>> shardRows(nShards);
        const auto tFit0 = Clock::now(); // the event loop above is the read
        ROOT::TThreadExecutor pool;
        pool.Foreach([&](unsigned shard){
            // Points are sorted and de-duplicated per track below, so the merge order does not matter.
            // Every point of an event is in the buckets of one shard, so the shard builds its tracks alone: its points go in one flat vector, sorted by track
            // and then by z (then x, y, ...), so a track is a run of rows with its points already in order. That is cheaper than a std::map of growing vectors.
            // pdg is not part of the sort or of the track key: a trackID is unique within an event, so it fixes the particle.
            std::vector<Point> pts;
            std::size_t nPts = 0;
            for (const auto &slotBuckets : buckets) nPts += slotBuckets[shard].size();
            pts.reserve(nPts);
            for (auto &slotBuckets : buckets) {
                pts.insert(pts.end(), slotBuckets[shard].begin(), slotBuckets[shard].end());
                slotBuckets[shard].clear(); slotBuckets[shard].shrink_to_fit();
            }
            // The ZDC rows of the shard, sorted by track key: a track's KE is looked up here, a track without a row never hit the ZDC. A particle is killed at the ZDC, so a track has one row.
            std::vector<ZdcHit> zdc;
            for (auto &slotZdc : zdcBuckets) {
                zdc.insert(zdc.end(), slotZdc[shard].begin(), slotZdc[shard].end());
                slotZdc[shard].clear(); slotZdc[shard].shrink_to_fit();
            }
            std::sort(zdc.begin(), zdc.end(), [](const ZdcHit &a, const ZdcHit &b){return a.key < b.key;});
            // The sort moves 24-byte keys (event and track packed in one word, then z), not the Points. The other position fields only break ties,
            // so points of equal z keep the order of the full comparison.
            struct Key {ULong64_t track; double z; std::uint32_t idx;};
            std::vector<Key> keys(pts.size());
            for (std::size_t i = 0; i < pts.size(); i++) {
                keys[i] = {(ULong64_t(std::uint32_t(pts[i].eventID)) << 32) | std::uint32_t(pts[i].trackID), pts[i].z, std::uint32_t(i)};
            }
            std::sort(keys.begin(), keys.end(), [&](const Key &a, const Key &b){
                if (a.track != b.track) return a.track < b.track;
                if (a.z != b.z) return a.z < b.z;
                const Point &p = pts[a.idx], &q = pts[b.idx];
                return std::tie(p.x, p.y, p.zD, p.xD, p.yD, p.entry) < std::tie(q.x, q.y, q.zD, q.xD, q.yD, q.entry);
            });
            auto at = [&](std::size_t i) -> const Point& {return pts[keys[i].idx];}; // the i-th point in sorted order
            // the runs of one track (begin, end), and the number of tracks of the event of each: every run counts, also one that is dropped below
            std::vector<std::pair<std::size_t, std::size_t>> runs;
            for (std::size_t i = 0, j; i < keys.size(); i = j) {
                for (j = i + 1; j < keys.size() && keys[j].track == keys[i].track; j++) {}
                runs.emplace_back(i, j);
            }
            // per event (a group of runs): the tracks of any pdg, and the most distinct charged tracks that entered one tracker plane (a plane is its centre z, 1 um steps;
            // the tracks ParticleInfo holds, every one of them hit 2+ planes or is a CE proton)
            std::vector<int> nProtonsInEvent(runs.size()), nMaxPerPlaneInEvent(runs.size());
            std::vector<std::pair<long, ULong64_t>> planeTracks;
            for (std::size_t i = 0, j; i < runs.size(); i = j) {
                for (j = i + 1; j < runs.size() && (keys[runs[j].first].track >> 32) == (keys[runs[i].first].track >> 32); j++) {}
                planeTracks.clear();
                for (std::size_t k = i; k < j; k++) for (std::size_t p = runs[k].first; p < runs[k].second; p++) {
                    if (at(p).zD != 0.0) planeTracks.emplace_back(std::lround(at(p).zD * 1000.0), keys[p].track);
                }
                std::sort(planeTracks.begin(), planeTracks.end());
                planeTracks.erase(std::unique(planeTracks.begin(), planeTracks.end()), planeTracks.end());
                int best = 0;
                for (std::size_t a = 0, b; a < planeTracks.size(); a = b) {
                    for (b = a + 1; b < planeTracks.size() && planeTracks[b].first == planeTracks[a].first; b++) {}
                    best = std::max(best, int(b - a));
                }
                for (std::size_t k = i; k < j; k++) {nProtonsInEvent[k] = int(j - i); nMaxPerPlaneInEvent[k] = best;}
            }

            std::vector<Row> &rows = shardRows[shard];
            Track t; // the points of the track in z order, repeats of the same position dropped; one Track for the whole shard, emptied for every track
            std::vector<std::size_t> blockHits; // tracker hits by z, reused the same way
            for (std::size_t run = 0; run < runs.size(); run++) {
                const auto key = std::make_tuple(at(runs[run].first).eventID, at(runs[run].first).trackID, at(runs[run].first).pdg);
                t.clear();
                for (std::size_t i = runs[run].first; i < runs[run].second; i++) {
                    const Point &p = at(i);
                    if (i > runs[run].first) {
                        const Point &q = at(i - 1);
                        if (p.x == q.x && p.y == q.y && p.z == q.z && p.xD == q.xD && p.yD == q.yD && p.zD == q.zD) continue;
                    }
                    t.push(p.x, p.y, p.z, p.xD, p.yD, p.zD, p.parentID, p.entry);
                }
                Row r;
                if (t.size() < 2) {continue;}

                const std::size_t endPoint   = std::max_element(t.z.begin(), t.z.end()) - t.z.begin();
                const std::size_t startPoint = std::min_element(t.z.begin(), t.z.end()) - t.z.begin();

                FitResult fit = fitTrack(t);
                if (!std::isfinite(fit.bx)  || !std::isfinite(fit.by))  {continue;}
                if (!std::isfinite(fit.bxD) || !std::isfinite(fit.byD)) {continue;}

                r.theta = std::atan2(std::hypot(fit.bx, fit.by), 1.0);
                r.phi   = std::atan2(fit.by, fit.bx);
                if (r.phi < 0) {r.phi += 2.0 * M_PI;}

                r.dTheta = std::atan2(std::hypot(fit.bxD, fit.byD), 1.0);
                r.dPhi   = std::atan2(fit.byD, fit.bxD);
                if (r.dPhi < 0) {r.dPhi += 2.0 * M_PI;}

                r.eventID = std::get<0>(key);
                r.trackID = std::get<1>(key);
                r.pdg     = std::get<2>(key);

                r.xVertex = t.x[startPoint];
                r.yVertex = t.y[startPoint];
                r.zVertex = t.z[startPoint];

                r.x0 = fit.x0;   r.y0 = fit.y0;   r.bx = fit.bx;   r.by = fit.by;
                r.x0D = fit.x0D; r.y0D = fit.y0D; r.bxD = fit.bxD; r.byD = fit.byD;

                // KE only from the ZDC row of this track (true KE and the smeared one); NAN when the track never hit the ZDC
                const ULong64_t trackKey = (ULong64_t(std::uint32_t(r.eventID)) << 32) | std::uint32_t(r.trackID);
                const auto zdcRow = std::lower_bound(zdc.begin(), zdc.end(), trackKey, [](const ZdcHit &h, ULong64_t k){return h.key < k;});
                const bool hitZdc = zdcRow != zdc.end() && zdcRow->key == trackKey;
                r.finalKE  = hitZdc ? zdcRow->KE : NAN;
                r.finalKED = hitZdc ? zdcRow->KED : NAN;
                r.parentID = t.ParentID[endPoint];

                r.zBxRecon = (fit.bxD != 0) ? -fit.x0D/fit.bxD : 0;
                r.zByRecon = (fit.byD != 0) ? -fit.y0D/fit.byD : 0;
                r.zBxReconDIFF = r.zVertex - r.zBxRecon;
                r.zByReconDIFF = r.zVertex - r.zByRecon;

                auto ce = ceTracks.find(r.eventID);
                r.isFromCE = (ce != ceTracks.end() && ce->second.count(r.parentID)) ? 1 : 0;
                r.inCE = (ce != ceTracks.end()) ? 1 : 0;

                FitQuality(t, fit.x0, fit.y0, fit.bx, fit.by, false, r.resRms, r.r2);
                FitQuality(t, fit.x0D, fit.y0D, fit.bxD, fit.byD, true, r.resDRms, r.r2D);
                r.impact  = ImpactToAxis(fit.x0, fit.y0, fit.bx, fit.by);
                r.impactD = ImpactToAxis(fit.x0D, fit.y0D, fit.bxD, fit.byD);

                r.nProtonsEv = nProtonsInEvent[run];
                r.nMaxChargedPerPlaneEv = nMaxPerPlaneInEvent[run];
                const EvCounts ec = (r.eventID >= 0 && (size_t)r.eventID < evCountsOf.size()) ? evCountsOf[r.eventID] : EvCounts();
                r.nZdcCharged = ec.zdcCharged; r.nTrackerHitsEv = ec.trackerHits; r.nTrackerTracksEv = ec.trackerTracks;
                {
                    blockHits.clear();
                    for (std::size_t i = 0; i < t.size(); i++) {if (t.zD[i] != 0.0) {blockHits.push_back(i);}}
                    std::sort(blockHits.begin(), blockHits.end(), [&](std::size_t a, std::size_t b){return t.zD[a] < t.zD[b];});
                    r.nHits = (int)blockHits.size();
                    r.firstXD = r.nHits ? t.xD[blockHits.front()] : NAN; r.firstYD = r.nHits ? t.yD[blockHits.front()] : NAN;
                    r.resHalfRatioD = NAN;
                    if (r.nHits > 2) {
                        const double r1 = ResidualRmsHits(t, blockHits, 0, 2, fit.x0D, fit.y0D, fit.bxD, fit.byD), r2h = ResidualRmsHits(t, blockHits, 2, blockHits.size(), fit.x0D, fit.y0D, fit.bxD, fit.byD);
                        if (r1 > 0) {r.resHalfRatioD = r2h / r1;}
                    }
                }
                rows.push_back(r);
            }
        }, ROOT::TSeqU(nShards));

        // TTree::Fill is not thread-safe: every shard fills a TTree of its own in a file the TBufferMerger gives it, in parallel, and the merger appends them into AnglesRaw.
        // A file's rows come from the shards in the order they finish, nothing reads them by position.
        const auto tWrite0 = Clock::now();
        ROOT::TThreadExecutor writePool;
        writePool.Foreach([&](unsigned shard){
            auto outFile = merger.GetFile();
            TTree *tree = new TTree("AnglesRaw", "Reconstructed Angles, every recon proton before the quality cut");
            // one branch per Row field; finalKE / finalKED are the true and smeared KE at the ZDC (NAN without a hit), the quality variables are in QualityCut.h, r2 / r2D are plotted by QualityStudy.cpp, not cut on,
            // firstXD / firstYD the block position (mm) of the first tracker hit (lowest z), resHalfRatioD the block residual of the hits after the first 2 over the first 2 (NAN with 2 or fewer),
            // nProtonsEv the tracks with a fit in the event (any pdg), nTrackerHitsEv / nTrackerTracksEv from the sim EventSummary (all charged particles), nMaxChargedPerPlaneEv the most charged tracks
            // (ParticleInfo) that entered one tracker plane in the event, inCE: the event has a carbon CE
            Row cur;
            tree->Branch("TrackID", &cur.trackID);
            tree->Branch("PDG", &cur.pdg);
            tree->Branch("fEvent", &cur.eventID);
            tree->Branch("outTheta", &cur.theta);
            tree->Branch("outPhi", &cur.phi);
            tree->Branch("outDTheta", &cur.dTheta);
            tree->Branch("outDPhi", &cur.dPhi);
            tree->Branch("x0Out", &cur.x0);
            tree->Branch("y0Out", &cur.y0);
            tree->Branch("bxOut", &cur.bx);
            tree->Branch("byOut", &cur.by);
            tree->Branch("x0DOut", &cur.x0D);
            tree->Branch("y0DOut", &cur.y0D);
            tree->Branch("bxDOut", &cur.bxD);
            tree->Branch("byDOut", &cur.byD);
            tree->Branch("finalKE", &cur.finalKE);
            tree->Branch("finalKED", &cur.finalKED);
            tree->Branch("inCE", &cur.inCE);
            tree->Branch("ParentID", &cur.parentID);
            tree->Branch("xVertexOut", &cur.xVertex);
            tree->Branch("yVertexOut", &cur.yVertex);
            tree->Branch("zVertexOut", &cur.zVertex);
            tree->Branch("z_bx_reconVertexOut", &cur.zBxRecon);
            tree->Branch("z_by_reconVertexOut", &cur.zByRecon);
            tree->Branch("z_bx_reconVertexDIFFOut", &cur.zBxReconDIFF);
            tree->Branch("z_by_reconVertexDIFFOut", &cur.zByReconDIFF);
            tree->Branch("isFromCE", &cur.isFromCE);
            tree->Branch("resRms", &cur.resRms);
            tree->Branch("resDRms", &cur.resDRms);
            tree->Branch("impact", &cur.impact);
            tree->Branch("impactD", &cur.impactD);
            tree->Branch("r2", &cur.r2);
            tree->Branch("r2D", &cur.r2D);
            tree->Branch("nProtonsEv", &cur.nProtonsEv);
            tree->Branch("nZdcCharged", &cur.nZdcCharged);
            tree->Branch("nTrackerHitsEv", &cur.nTrackerHitsEv);
            tree->Branch("nTrackerTracksEv", &cur.nTrackerTracksEv);
            tree->Branch("nHits", &cur.nHits);
            tree->Branch("nMaxChargedPerPlaneEv", &cur.nMaxChargedPerPlaneEv);
            tree->Branch("firstXD", &cur.firstXD);
            tree->Branch("firstYD", &cur.firstYD);
            tree->Branch("resHalfRatioD", &cur.resHalfRatioD);
            for (const Row &r : shardRows[shard]) {
                cur = r;
                tree->Fill();
            }
            outFile->Write();
            shardRows[shard].clear(); shardRows[shard].shrink_to_fit();
        }, ROOT::TSeqU(nShards));
        printf(">>> ReconAngleAll group %zu: read %.1f s, sort and fit %.1f s, write %.1f s\n", first, secs(tRead0, tFit0), secs(tFit0, tWrite0), secs(tWrite0, Clock::now())); fflush(stdout);
    }
    printf(">>> ReconAngleAll loop total %.1f s (the merger closes after this)\n", secs(tStart, Clock::now())); fflush(stdout);
}
