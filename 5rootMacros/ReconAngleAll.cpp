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
#include <fstream>
#include <thread>

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

// Weight of a point in every line fit, relative to the 2 points with the lowest z (weight 1): the points near the vertex have scattered the least, the later ones more.
// Reads the raw files names[lo] to names[hi - 1] once through and throws the bytes away: the next group of the loop below then finds them in the page cache. On LCRC the
// first read of a group of 8 files took 11.6 s against 3.3 s when they were cached; this runs beside the fit and the write of the group before.
static void PrefetchFiles(const std::vector<std::string>& names, std::size_t lo, std::size_t hi){
    std::vector<char> buf(8u << 20);
    for (std::size_t i = lo; i < hi; i++){
        std::ifstream f(names[i], std::ios::binary);
        while (f) f.read(buf.data(), (std::streamsize)buf.size());
    }
}
static constexpr double kWeightFar = 0.25;

// Weighted least-squares line x = x0 + bx z, y = y0 + by z through the points of t: all of them, or with block the detector block points (zD != 0) at their block positions.
// The 2 points with the lowest z of the set have weight 1, every other point kWeightFar.
// Leaves the four outputs alone (NAN from the caller) with fewer than 2 points or no spread in z.
static void LineFit(const Track& t, bool block, double& x0, double& y0, double& bx, double& by){
    const std::vector<double>& X = block ? t.xD : t.x; const std::vector<double>& Y = block ? t.yD : t.y; const std::vector<double>& Z = block ? t.zD : t.z;
    const std::size_t none = t.size();
    std::size_t first = none, second = none; // the 2 points with the lowest z
    std::size_t N = 0;
    for (std::size_t i = 0; i < t.size(); i++) {
        if (block && t.zD[i] == 0.0) {continue;}
        N++;
        if (first == none || Z[i] < Z[first]) {second = first; first = i;}
        else if (second == none || Z[i] < Z[second]) {second = i;}
    }
    if (N < 2) {return;}
    const auto weight = [&](std::size_t i){return (i == first || i == second) ? 1.0 : kWeightFar;};
    double sw = 0.0, meanX = 0.0, meanY = 0.0, meanZ = 0.0;
    for (std::size_t i = 0; i < t.size(); i++) {
        if (block && t.zD[i] == 0.0) {continue;}
        const double w = weight(i);
        sw += w; meanX += w * X[i]; meanY += w * Y[i]; meanZ += w * Z[i];
    }
    meanX /= sw; meanY /= sw; meanZ /= sw;

    double Szz = 0.0, Szx = 0.0, Szy = 0.0;
    for (std::size_t i = 0; i < t.size(); i++) {
        if (block && t.zD[i] == 0.0) {continue;}
        const double w = weight(i), dz = Z[i] - meanZ;
        Szz += w * dz * dz;
        Szx += w * dz * (X[i] - meanX);
        Szy += w * dz * (Y[i] - meanY);
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

    // Carbon block centres (z, mm, ascending) and half thickness: the vertex tube of a track is its fitted line across the block just upstream of its first tracker hit.
    const auto settings = ReadSettings();
    const std::vector<double> blockZ = StackZ(settings).blockZ;
    const double halfBlockMm = settings.at("CarbonThickness") * 1000.0 / 2.0;
    // Position error of one tracker hit (mm): the tracker is 400 mm over nCols_blocks cells (2construction.cc), a hit sits at a cell centre, so the error is pitch / sqrt(12).
    const double hitSigmaMm = 400.0 / settings.at("nCols_blocks") / std::sqrt(12.0);

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
        std::thread prefetch; // the files of the next group, read into the page cache while this group is processed
        if (first + filesPerPass < rawNames.size()) prefetch = std::thread(PrefetchFiles, std::cref(rawNames), first + filesPerPass, std::min(rawNames.size(), first + 2 * (size_t)filesPerPass));
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
            int tubeBlock, inTubeZ;
            double tubeZ0, tubeZ1, tubeDist, tubeDistD, tubeRadius, tubeX0, tubeY0, tubeX1, tubeY1, thetaNMin, thetaNMax, phiNMid, phiNHalf;
            double thetaRelN, thetaRelNMin, thetaRelNMax, phiRelN, phiRelNHalf, thetaRelND, thetaRelNDMin, thetaRelNDMax, phiRelND, phiRelNDHalf; // proton angles in the frame of the estimated neutron
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
                // Vertex tube: the fitted lines (true and block-position) extended back across the carbon block just upstream of the first tracker hit.
                // tubeBlock -1 (all tube values NAN, inTubeZ 0) with no tracker hit or a first hit upstream of every block.
                r.tubeBlock = -1; r.inTubeZ = 0;
                r.tubeZ0 = r.tubeZ1 = r.tubeDist = r.tubeDistD = r.tubeRadius = NAN;
                r.thetaNMin = r.thetaNMax = r.phiNMid = r.phiNHalf = r.tubeX0 = r.tubeY0 = r.tubeX1 = r.tubeY1 = NAN;
                r.thetaRelN = r.thetaRelNMin = r.thetaRelNMax = r.phiRelN = r.phiRelNHalf = NAN;
                r.thetaRelND = r.thetaRelNDMin = r.thetaRelNDMax = r.phiRelND = r.phiRelNDHalf = NAN;
                if (r.nHits > 0) {
                    const double zFirst = t.zD[blockHits.front()];
                    const int idx = int(std::upper_bound(blockZ.begin(), blockZ.end(), zFirst) - blockZ.begin()) - 1;
                    if (idx >= 0) {
                        r.tubeBlock = idx;
                        r.tubeZ0 = blockZ[idx] - halfBlockMm; r.tubeZ1 = blockZ[idx] + halfBlockMm;
                        r.inTubeZ = (r.zVertex >= r.tubeZ0 && r.zVertex <= r.tubeZ1) ? 1 : 0;
                        r.tubeDist  = std::hypot(r.xVertex - (fit.x0 + fit.bx * r.zVertex),   r.yVertex - (fit.y0 + fit.by * r.zVertex));
                        // The tube line is the weighted block-position fit itself (x0D, bxD, y0D, byD, LineFit): the 2 hits right after the CE block (the first 2 in z, the first hit is already
                        // behind the block) have weight 1, every later hit kWeightFar, since the later hits have scattered more. A hit of weight w has error hitSigmaMm / sqrt(w), so
                        // sigma^2 = s^2 (1/W + (z - zMean)^2 / Sxx) with W = sum w, Sxx = sum w (z - zMean)^2, taken at the block centre. NAN with fewer than 2 hits or no spread in zD.
                        double sw = 0.0, zMean = 0.0;
                        for (std::size_t k = 0; k < blockHits.size(); k++) {
                            const double w = k < 2 ? 1.0 : kWeightFar;
                            sw += w; zMean += w * t.zD[blockHits[k]];
                        }
                        zMean /= sw;
                        double sxx = 0.0;
                        for (std::size_t k = 0; k < blockHits.size(); k++) {
                            const double w = k < 2 ? 1.0 : kWeightFar, dzh = t.zD[blockHits[k]] - zMean;
                            sxx += w * dzh * dzh;
                        }
                        if (blockHits.size() >= 2 && sxx > 0.0) {
                            const auto lineX = [&](double zz){return fit.x0D + fit.bxD * zz;};
                            const auto lineY = [&](double zz){return fit.y0D + fit.byD * zz;};
                            const double dz = blockZ[idx] - zMean;
                            constexpr double kTubeSigmas = 2.0; // holds 92% of the CE vertices of the 801x801 run (1.5: 67%, 2.25: 98%, 2.45: 99%); PlotTubeDebug.cpp repeats this number
                            r.tubeRadius = kTubeSigmas * hitSigmaMm * std::sqrt(1.0 / sw + dz * dz / sxx);
                            r.tubeDistD = std::hypot(r.xVertex - lineX(r.zVertex), r.yVertex - lineY(r.zVertex));
                            r.tubeX0 = lineX(r.tubeZ0); r.tubeY0 = lineY(r.tubeZ0);
                            r.tubeX1 = lineX(r.tubeZ1); r.tubeY1 = lineY(r.tubeZ1);
                            // Initial neutron direction range: the neutron flies straight from the beam origin (0, 0, 0) to its CE vertex, so every point of the tube (the weighted line over the block
                            // z range, radius tubeRadius) gives a direction. Sampled along the line: at a sample (x, y, z) with rho = hypot(x, y), theta runs from atan((rho - R) / z) (0 when rho < R) to
                            // atan((rho + R) / z), phi covers phiCentre +- asin(R / rho), all of 2 pi when rho < R. thetaNMin / thetaNMax in rad, phiNMid / phiNHalf: centre in [0, 2 pi) and half width of the phi range.
                            const double R = r.tubeRadius, phiRef = std::atan2(lineY(blockZ[idx]), lineX(blockZ[idx]));
                            constexpr int kTubeSamples = 41;
                            // The neutron does not always start at the origin: with randPos the gun draws its start from a gaussian (5generator.cc: sigma x 0.0055 mm, y 0.003 mm, z 15 mm), not stored
                            // per event. The start counts as an error like the tube: kTubeSigmas sigma of the beam spot widens the tube radius in x, y (rBeam) and moves the start by +-zBeam along z.
                            constexpr double kBeamSigmaX = 0.0055, kBeamSigmaY = 0.003, kBeamSigmaZ = 15.0;
                            const double rBeam = kTubeSigmas * std::hypot(kBeamSigmaX, kBeamSigmaY), zBeam = kTubeSigmas * kBeamSigmaZ, Rn = R + rBeam;
                            double thLo = INFINITY, thHi = 0.0, phLo = INFINITY, phHi = -INFINITY;
                            for (int s = 0; s < kTubeSamples; s++) {
                                const double zs = r.tubeZ0 + (r.tubeZ1 - r.tubeZ0) * s / (kTubeSamples - 1);
                                const double xs = lineX(zs), ys = lineY(zs), rho = std::hypot(xs, ys);
                                thLo = std::min(thLo, std::atan2(std::max(0.0, rho - Rn), zs + zBeam)); // the start behind the sample is the longest flight, the smallest angle
                                thHi = std::max(thHi, std::atan2(rho + Rn, zs - zBeam));
                                const double half = rho > Rn ? std::asin(Rn / rho) : M_PI;
                                const double d = std::remainder(std::atan2(ys, xs) - phiRef, 2.0 * M_PI); // phi of the sample relative to the centre of the tube, in (-pi, pi]
                                phLo = std::min(phLo, d - half); phHi = std::max(phHi, d + half);
                            }
                            r.thetaNMin = thLo; r.thetaNMax = thHi;
                            r.phiNHalf = std::min(M_PI, (phHi - phLo) / 2.0);
                            r.phiNMid = std::fmod(phiRef + (phLo + phHi) / 2.0 + 4.0 * M_PI, 2.0 * M_PI);

                            // The proton direction of each fit (true position: bx, by; block position: bxD, byD) in the frame where the estimated neutron points along z. The estimate is the direction
                            // from the beam origin to the tube centre at the middle of the block; the Min / Max (theta) and Half (phi, the largest offset from the estimate) come from every neutron
                            // the start and the tube allow: the start z at 0 and +-zBeam (the x, y spread of the start is inside Rn), the end at the block front, middle and back, each at the tube centre and
                            // on a circle of radius Rn around it (8 points).
                            const double pTrue[3] = {fit.bx, fit.by, 1.0}, pBlock[3] = {fit.bxD, fit.byD, 1.0};
                            const double zMid = blockZ[idx], nCentre[3] = {lineX(zMid), lineY(zMid), zMid};
                            double thC[2], phC[2], thMin[2], thMax[2], phHalf[2];
                            // The block direction is the slope of the block-position line fit, and a block hit sits on a grid (error hitSigmaMm = pitch / sqrt(12) per hit): the slope error of the weighted fit is
                            // kTubeSigmas * hitSigmaMm / sqrt(Sxx) in x and in y. The block range (thetaRelNDMin/Max, phiRelNDHalf) also covers the direction at the four corners of that slope error.
                            const double slopeSig = kTubeSigmas * hitSigmaMm / std::sqrt(sxx);
                            const double vx[5] = {0.0, 1.0, 1.0, -1.0, -1.0}, vy[5] = {0.0, 1.0, -1.0, 1.0, -1.0};
                            RelToDir(pTrue, nCentre, thC[0], phC[0]); RelToDir(pBlock, nCentre, thC[1], phC[1]);
                            for (int f = 0; f < 2; f++) {thMin[f] = thC[f]; thMax[f] = thC[f]; phHalf[f] = 0.0;}
                            for (const double zStart : {-zBeam, 0.0, zBeam}) for (const double zEnd : {r.tubeZ0, zMid, r.tubeZ1}) for (int k = 0; k < 9; k++) {
                                const double a = (k - 1) * M_PI / 4.0, ex = lineX(zEnd) + (k ? Rn * std::cos(a) : 0.0), ey = lineY(zEnd) + (k ? Rn * std::sin(a) : 0.0);
                                const double nS[3] = {ex, ey, zEnd - zStart};
                                for (int f = 0; f < 2; f++) for (int v = 0; v < (f ? 5 : 1); v++) { // the block direction also moves by +-slopeSig in both slopes (the corners v 1 to 4)
                                    const double pV[3] = {pBlock[0] + (f ? vx[v] * slopeSig : 0.0), pBlock[1] + (f ? vy[v] * slopeSig : 0.0), 1.0};
                                    double th, ph; RelToDir(f ? pV : pTrue, nS, th, ph);
                                    thMin[f] = std::min(thMin[f], th); thMax[f] = std::max(thMax[f], th);
                                    phHalf[f] = std::max(phHalf[f], std::fabs(std::remainder(ph - phC[f], 2.0 * M_PI)));
                                }
                            }
                            r.thetaRelN = thC[0]; r.thetaRelNMin = thMin[0]; r.thetaRelNMax = thMax[0]; r.phiRelN = phC[0]; r.phiRelNHalf = phHalf[0];
                            r.thetaRelND = thC[1]; r.thetaRelNDMin = thMin[1]; r.thetaRelNDMax = thMax[1]; r.phiRelND = phC[1]; r.phiRelNDHalf = phHalf[1];
                        }
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
            // vertex tube: tubeBlock the index of the carbon block just upstream of the first tracker hit (-1 none), tubeZ0 / tubeZ1 its front and back face (mm),
            // inTubeZ the true vertex z lies in the block, tubeDist / tubeDistD the transverse distance (mm) of the true vertex from the true / D line at the vertex z;
            // every line fit (x0Out.. and x0DOut..) is weighted: the 2 lowest-z points weight 1, the rest kWeightFar (0.25). The tube centre at a face is x0DOut + bxDOut * tubeZ (tubeX0.. hold it),
            // tubeDist / tubeDistD the distance of the true vertex from the true / block fit, tubeRadius the 2 sigma fit error (mm) of the block fit at the block centre (tubeDistD <= tubeRadius: vertex inside the tube)
            tree->Branch("tubeBlock", &cur.tubeBlock);
            tree->Branch("inTubeZ", &cur.inTubeZ);
            tree->Branch("tubeZ0", &cur.tubeZ0);
            tree->Branch("tubeZ1", &cur.tubeZ1);
            // initial neutron direction range of the tube, neutron from the beam origin (0, 0, 0): thetaNMin / thetaNMax (rad), phiNMid (centre, [0, 2 pi) rad) and phiNHalf (half width, pi = every phi, the tube covers the axis)
            tree->Branch("thetaNMin", &cur.thetaNMin);
            tree->Branch("thetaNMax", &cur.thetaNMax);
            tree->Branch("phiNMid", &cur.phiNMid);
            tree->Branch("phiNHalf", &cur.phiNHalf);
            // tubeX0 / tubeY0 and tubeX1 / tubeY1: the centre of the tube (the weighted line) at the front and back face of the block (mm)
            tree->Branch("tubeX0", &cur.tubeX0);
            tree->Branch("tubeY0", &cur.tubeY0);
            tree->Branch("tubeX1", &cur.tubeX1);
            tree->Branch("tubeY1", &cur.tubeY1);
            // proton angles (rad) in the frame of the estimated neutron (RelToDir): thetaRelN / phiRelN for the true-position fit (outTheta, outPhi), thetaRelND / phiRelND for the block fit (outDTheta, outDPhi),
            // the neutron estimate is the direction from the beam origin to the tube centre at the block middle; Min / Max the theta range and Half the phi half width over every neutron the start spread and the tube allow
            tree->Branch("thetaRelN", &cur.thetaRelN);
            tree->Branch("thetaRelNMin", &cur.thetaRelNMin);
            tree->Branch("thetaRelNMax", &cur.thetaRelNMax);
            tree->Branch("phiRelN", &cur.phiRelN);
            tree->Branch("phiRelNHalf", &cur.phiRelNHalf);
            tree->Branch("thetaRelND", &cur.thetaRelND);
            tree->Branch("thetaRelNDMin", &cur.thetaRelNDMin);
            tree->Branch("thetaRelNDMax", &cur.thetaRelNDMax);
            tree->Branch("phiRelND", &cur.phiRelND);
            tree->Branch("phiRelNDHalf", &cur.phiRelNDHalf);
            tree->Branch("tubeRadius", &cur.tubeRadius);
            tree->Branch("tubeDist", &cur.tubeDist);
            tree->Branch("tubeDistD", &cur.tubeDistD);
            for (const Row &r : shardRows[shard]) {
                cur = r;
                tree->Fill();
            }
            outFile->Write();
            shardRows[shard].clear(); shardRows[shard].shrink_to_fit();
        }, ROOT::TSeqU(nShards));
        printf(">>> ReconAngleAll group %zu: read %.1f s, sort and fit %.1f s, write %.1f s\n", first, secs(tRead0, tFit0), secs(tFit0, tWrite0), secs(tWrite0, Clock::now())); fflush(stdout);
        if (prefetch.joinable()) prefetch.join();
    }
    printf(">>> ReconAngleAll loop total %.1f s (the merger closes after this)\n", secs(tStart, Clock::now())); fflush(stdout);
}
