#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "QualityCut.h"
#include "CutCache.h"
#include "FitCore.h"

struct FitResult {double x0, y0, bx, by;
    double x0D, y0D, bxD, byD;};
// The points of one track as a struct of arrays: point i is (x[i], y[i], z[i], ...), one array per quantity. KE: true KE, KED: detector KE,
// entry: row number in ParticleInfo (rows of a track are written in the order it hit them).
struct Track {
    std::vector<double> x, y, z, xD, yD, zD, KE, KED;
    std::vector<int> ParentID;
    std::vector<ULong64_t> entry;
    std::size_t size() const {return z.size();}
    void push(double x_, double y_, double z_, double xD_, double yD_, double zD_, double KE_, int parent, double KED_, ULong64_t entry_){
        x.push_back(x_); y.push_back(y_); z.push_back(z_); xD.push_back(xD_); yD.push_back(yD_); zD.push_back(zD_);
        KE.push_back(KE_); ParentID.push_back(parent); KED.push_back(KED_); entry.push_back(entry_);
    }
    void append(const Track& o){
        for (std::size_t i = 0; i < o.size(); i++) push(o.x[i], o.y[i], o.z[i], o.xD[i], o.yD[i], o.zD[i], o.KE[i], o.ParentID[i], o.KED[i], o.entry[i]);
    }
    Track select(const std::vector<std::size_t>& idx) const { // the points idx, in that order
        Track t;
        for (std::size_t i : idx) t.push(x[i], y[i], z[i], xD[i], yD[i], zD[i], KE[i], ParentID[i], KED[i], entry[i]);
        return t;
    }
};

FitResult fitTrack(const Track& t){
    FitResult fit{NAN, NAN, NAN, NAN,
                NAN, NAN, NAN, NAN};

    const std::size_t N = t.size();
    if (N < 2) {return fit;}

    double meanX = 0.0, meanY = 0.0, meanZ = 0.0;
    for (std::size_t i = 0; i < N; i++) {
        meanX += t.x[i];
        meanY += t.y[i];
        meanZ += t.z[i];
    }
    meanX /= N;
    meanY /= N;
    meanZ /= N;

    double Szz = 0.0, Szx = 0.0, Szy = 0.0;
    for (std::size_t i = 0; i < N; i++) {
        double dz = t.z[i] - meanZ;
        Szz += dz * dz;
        Szx += dz * (t.x[i] - meanX);
        Szy += dz * (t.y[i] - meanY);
    }

    if (std::abs(Szz) >= 1e-20) {
        fit.bx = Szx / Szz;
        fit.by = Szy / Szz;
        fit.x0 = meanX - fit.bx * meanZ;
        fit.y0 = meanY - fit.by * meanZ;
    }

    double meanXD = 0.0, meanYD = 0.0, meanZD = 0.0;
    std::size_t ND = 0;
    for (std::size_t i = 0; i < N; i++) {
        if (t.zD[i] == 0.0) {continue;}
        meanXD += t.xD[i];
        meanYD += t.yD[i];
        meanZD += t.zD[i];
        ND++;
    }

    if (ND < 2) {return fit;}

    meanXD /= ND;
    meanYD /= ND;
    meanZD /= ND;

    double SzzD = 0.0, SzxD = 0.0, SzyD = 0.0;
    for (std::size_t i = 0; i < N; i++) {
        if (t.zD[i] == 0.0) {continue;}
        double dzD = t.zD[i] - meanZD;
        SzzD += dzD * dzD;
        SzxD += dzD * (t.xD[i] - meanXD);
        SzyD += dzD * (t.yD[i] - meanYD);
    }

    if (std::abs(SzzD) >= 1e-20) {
        fit.bxD = SzxD / SzzD;
        fit.byD = SzyD / SzzD;
        fit.x0D = meanXD - fit.bxD * meanZD;
        fit.y0D = meanYD - fit.byD * meanZD;
    }

    return fit;
}

// RMS transverse distance (mm) of the points to the fitted line x = x0 + bx z, y = y0 + by z. block: only the detector block points (zD != 0), at their block positions.
static double ResidualRms(const Track& t, double x0, double y0, double bx, double by, bool block){
    double sum = 0.0; int n = 0;
    for (std::size_t i = 0; i < t.size(); i++) {
        const double x = block ? t.xD[i] : t.x[i], y = block ? t.yD[i] : t.y[i], z = block ? t.zD[i] : t.z[i];
        if (block && z == 0.0) {continue;}
        sum += std::pow(x - (x0 + bx * z), 2) + std::pow(y - (y0 + by * z), 2);
        n++;
    }
    return n > 0 ? std::sqrt(sum / n) : NAN;
}

// R^2 of the same line fit, x and y together: 1 - sum of squared residuals / sum of squared deviations from the mean hit position.
// Same points as ResidualRms. NAN with fewer than 3 points or no spread.
static double RSquared(const Track& t, double x0, double y0, double bx, double by, bool block){
    double sx = 0, sy = 0; int n = 0;
    for (std::size_t i = 0; i < t.size(); i++) {
        if (block && t.zD[i] == 0.0) {continue;}
        sx += block ? t.xD[i] : t.x[i]; sy += block ? t.yD[i] : t.y[i]; n++;
    }
    if (n < 3) {return NAN;}
    const double mx = sx / n, my = sy / n;
    double ssRes = 0, ssTot = 0;
    for (std::size_t i = 0; i < t.size(); i++) {
        const double x = block ? t.xD[i] : t.x[i], y = block ? t.yD[i] : t.y[i], z = block ? t.zD[i] : t.z[i];
        if (block && z == 0.0) {continue;}
        ssRes += std::pow(x - (x0 + bx * z), 2) + std::pow(y - (y0 + by * z), 2);
        ssTot += std::pow(x - mx, 2) + std::pow(y - my, 2);
    }
    return ssTot > 0 ? 1.0 - ssRes / ssTot : NAN;
}

// Closest approach (mm) of the line x = x0 + bx z, y = y0 + by z to the beam axis (x = y = 0).
static double ImpactToAxis(double x0, double y0, double bx, double by){
    const double slope = std::hypot(bx, by);
    return slope > 0 ? std::fabs(x0 * by - y0 * bx) / slope : std::hypot(x0, y0);
}

// One pass fits every proton and writes the tree AnglesRaw with the true (finalKE) and the detector (finalKED) KE of the track end.
void ReconAngleAll(){
    ROOT::EnableImplicitMT(NThreads());

    // CE neutrons (same selection as CutEvents.cpp): event -> CE track IDs.
    // the carbon-block CE come from the cut cache of CarbonCeCut (the seed cut_CE.root, CutCache.h): every row passes, no filter
    const CacheTree ceCarbon = GetCutTree("output0_t*.root", "ChargeExchange_values", {CarbonCeCut()});
    if (!ceCarbon.tree){Error("ReconAngleAll", "no cut cache for the carbon CE (run CutEvents.cpp)"); return;}
    ROOT::RDataFrame dfCE(*ceCarbon.tree);
    std::vector<std::map<int, std::set<int>>> slotCeTracks(dfCE.GetNSlots());
    dfCE.ForeachSlot([&](unsigned slot, int eventID, int trackID){
        slotCeTracks[slot][eventID].insert(trackID);
    }, {"fEvent", "TrackID"});
    std::map<int, std::set<int>> ceTracks;
    for (auto &m : slotCeTracks) {
        for (auto &[ev, ids] : m) ceTracks[ev].insert(ids.begin(), ids.end());
    }

    TFile *outFile = new TFile("cuts/cut_angles_all.root", "RECREATE");
    // The one proton tree: every recon proton before the quality cut (QualityStudy.cpp). The old trees Angles, AnglesCE and their _DE copies are views of it
    // (AngleSelection, CutFiles.h): passQuality, inCE pick the rows, finalKE / finalKED the true / detector KE.
    TTree *outTree = new TTree("AnglesRaw", "Reconstructed Angles, every recon proton before the quality cut");

    int outTrackID, outPDG, outEventID, fParentID;
    double outTheta, outPhi, outDTheta, outDPhi;
    double x0Out, y0Out, bxOut, byOut;
    double x0DOut, y0DOut, bxDOut, byDOut;
    double xVertexOut, yVertexOut, zVertexOut;
    double z_bx_reconVertexOut, z_by_reconVertexOut;
    double z_bx_reconVertexDIFFOut, z_by_reconVertexDIFFOut;
    double finalKE, finalKED; // true and detector KE of the track end
    double resRms, resDRms, impact, impactD; // quality variables, see QualityCut.h
    double r2, r2D; // R^2 of the line fit (true and block position), plotted by QualityStudy.cpp, not cut on
    double firstXD = NAN, firstYD = NAN; // block position (mm) of the first tracker hit (lowest z)
    int nProtonsEv, zOrderOK, nHits; // protons with a track in the event, planes hit in increasing z (1/0), tracker hits (block points)
    double resHalfRatioD; // block residual of the tracker hits after the first 2 over the first 2 hits, NAN with 2 or fewer hits
    int passResD = 1; // block residual cut, cut on only by the block plots (QualityCut.h)
    int isFromCE;
    int passQuality, inCE; // z order and quality cut passed (1/0), event has a carbon CE (1/0)

    {
        TTree *tree = outTree;
        tree->Branch("TrackID", &outTrackID);
        tree->Branch("PDG", &outPDG);
        tree->Branch("fEvent", &outEventID);
        tree->Branch("outTheta", &outTheta);
        tree->Branch("outPhi", &outPhi);
        tree->Branch("outDTheta", &outDTheta);
        tree->Branch("outDPhi", &outDPhi);
        tree->Branch("x0Out", &x0Out);
        tree->Branch("y0Out", &y0Out);
        tree->Branch("bxOut", &bxOut);
        tree->Branch("byOut", &byOut);
        tree->Branch("x0DOut", &x0DOut);
        tree->Branch("y0DOut", &y0DOut);
        tree->Branch("bxDOut", &bxDOut);
        tree->Branch("byDOut", &byDOut);
        tree->Branch("finalKE", &finalKE);
        tree->Branch("finalKED", &finalKED);
        tree->Branch("passQuality", &passQuality);
        tree->Branch("inCE", &inCE);
        tree->Branch("ParentID", &fParentID);
        tree->Branch("xVertexOut", &xVertexOut);
        tree->Branch("yVertexOut", &yVertexOut);
        tree->Branch("zVertexOut", &zVertexOut);
        tree->Branch("z_bx_reconVertexOut", &z_bx_reconVertexOut);
        tree->Branch("z_by_reconVertexOut", &z_by_reconVertexOut);
        tree->Branch("z_bx_reconVertexDIFFOut", &z_bx_reconVertexDIFFOut);
        tree->Branch("z_by_reconVertexDIFFOut", &z_by_reconVertexDIFFOut);
        tree->Branch("isFromCE", &isFromCE);
        tree->Branch("resRms", &resRms);
        tree->Branch("resDRms", &resDRms);
        tree->Branch("impact", &impact);
        tree->Branch("impactD", &impactD);
        tree->Branch("r2", &r2);
        tree->Branch("r2D", &r2D);
        tree->Branch("passResD", &passResD);
        tree->Branch("nProtonsEv", &nProtonsEv);
        tree->Branch("zOrderOK", &zOrderOK);
        tree->Branch("nHits", &nHits);
        tree->Branch("firstXD", &firstXD); tree->Branch("firstYD", &firstYD);
        tree->Branch("resHalfRatioD", &resHalfRatioD);
    }

    // One raw file at a time: an event, and so every track, lives in a single output0_t<thread>.root, so the tracks of a file are complete
    // and only that file's proton points are held in memory (all files at once does not fit at 100M events).
    TChain rawFiles("ParticleInfo");
    rawFiles.Add("output0_t*.root");
    for (const TObject* rawFile : *rawFiles.GetListOfFiles()) {
        ROOT::RDataFrame df("ParticleInfo", rawFile->GetTitle());

        std::vector<std::map<std::tuple<int, int, int>, Track>> slotTracks(df.GetNSlots());
        // boundaryType 2 is the ZDC row, left out of the fit points; protons that never reach the ZDC are reconstructed too
        // The selection is a typed Filter lambda, not a Filter string: a string Filter is JIT-compiled, which fails on a ROOT whose cling cannot materialize libc++ symbols (the Mac).
        // A Filter node, unlike a test inside the ForeachSlot lambda, reads the other branches only for rows that pass (protons), not for every particle row.
        df.Filter([](int pdg, int boundaryType){return pdg == 2212 && boundaryType != 2;}, {"pdg", "boundaryType"})
          .ForeachSlot([&](unsigned slot,
            int eventID, int trackID, int pdg,
            double x, double y, double z, double xD, double yD, double zD, double KE, double KED, int ParentID, ULong64_t entry){
            slotTracks[slot][{eventID, trackID, pdg}].push(x, y, z, xD, yD, zD, KE, ParentID, KED, entry);
        }, {"fEvent", "TrackID", "pdg", "fX", "fY", "fZ", "fDx", "fDy", "fDz", "fKEnergy", "fDEnergy", "ParentID", "rdfentry_"});
        // The merge and the fits run on every thread: tracks are sharded by event (all tracks of an event land in one shard, so the per-event
        // proton count is complete per shard), each shard builds its tracks, fits them and keeps the results as rows. Only the TTree fills below are serial.
        struct Row {
            int trackID, pdg, eventID, parentID, isFromCE, inCE, passResD, nProtonsEv, zOrderOK, nHits;
            bool passQuality;
            double theta, phi, dTheta, dPhi, x0, y0, bx, by, x0D, y0D, bxD, byD, finalKE, finalKED;
            double xVertex, yVertex, zVertex, zBxRecon, zByRecon, zBxReconDIFF, zByReconDIFF;
            double resRms, resDRms, impact, impactD, r2, r2D, firstXD, firstYD, resHalfRatioD;
        };
        const unsigned nShards = NThreads();
        std::vector<std::vector<Row>> shardRows(nShards);
        ROOT::TThreadExecutor pool;
        pool.Foreach([&](unsigned shard){
            // Points are sorted and de-duplicated per track below, so the merge order does not matter.
            // Each key is taken by exactly one shard, so its Track is moved out of the slot map, not copied.
            std::map<std::tuple<int, int, int>, Track> tracks;
            for (auto &slotMap : slotTracks) {
                for (auto &[key, t] : slotMap) {
                    if ((unsigned)std::get<0>(key) % nShards != shard) continue;
                    auto it = tracks.find(key);
                    if (it == tracks.end()) tracks.emplace(key, std::move(t)); else it->second.append(t);
                }
            }

            std::map<int, int> nProtonsInEvent;
            for (const auto &kv : tracks) nProtonsInEvent[std::get<0>(kv.first)]++;

            std::vector<Row> &rows = shardRows[shard];
            for (auto &kv : tracks) {
                const auto &key = kv.first;
                Track &t = kv.second; // not a structured binding: the sort lambdas below capture it, which C++17 does not allow for a binding
                Row r;
                // Tracker hit order: rows in ParticleInfo come in the order the track hit them, fDz (zD) is the tracker plane z.
                // Planes must be hit with strictly increasing z (repeat rows of the same plane are fine); 4 to 2 to 3 or 4 to 3 is dropped.
                {
                    std::vector<std::size_t> hits;
                    for (std::size_t i = 0; i < t.size(); i++) {if (t.zD[i] != 0.0) {hits.push_back(i);}}
                    std::sort(hits.begin(), hits.end(), [&](std::size_t a, std::size_t b){return t.entry[a] < t.entry[b];});
                    bool outOfOrder = false;
                    for (std::size_t i = 1; i < hits.size() && !outOfOrder; i++) {
                        if (t.zD[hits[i]] < t.zD[hits[i - 1]]) {outOfOrder = true;}
                    }
                    r.zOrderOK = outOfOrder ? 0 : 1; // out of order tracks stay in AnglesRaw for the plots, the other trees drop them below
                }
                { // points sorted by z (then x, y), repeats of the same position dropped
                    std::vector<std::size_t> order(t.size());
                    for (std::size_t i = 0; i < order.size(); i++) order[i] = i;
                    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b){
                        if (t.z[a] != t.z[b]) return t.z[a] < t.z[b];
                        if (t.x[a] != t.x[b]) return t.x[a] < t.x[b];
                        return t.y[a] < t.y[b];
                    });
                    order.erase(std::unique(order.begin(), order.end(), [&](std::size_t a, std::size_t b){
                        return t.x[a]  == t.x[b]  && t.y[a]  == t.y[b]  && t.z[a]  == t.z[b] &&
                               t.xD[a] == t.xD[b] && t.yD[a] == t.yD[b] && t.zD[a] == t.zD[b];
                    }), order.end());
                    t = t.select(order);
                }

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

                r.finalKE  = t.KE[endPoint];
                r.finalKED = t.KED[endPoint];
                r.parentID = t.ParentID[endPoint];

                r.zBxRecon = (fit.bxD != 0) ? -fit.x0D/fit.bxD : 0;
                r.zByRecon = (fit.byD != 0) ? -fit.y0D/fit.byD : 0;
                r.zBxReconDIFF = r.zVertex - r.zBxRecon;
                r.zByReconDIFF = r.zVertex - r.zByRecon;

                auto ce = ceTracks.find(r.eventID);
                r.isFromCE = (ce != ceTracks.end() && ce->second.count(r.parentID)) ? 1 : 0;
                r.inCE = (ce != ceTracks.end()) ? 1 : 0;

                r.resRms  = ResidualRms(t, fit.x0, fit.y0, fit.bx, fit.by, false);
                r.resDRms = ResidualRms(t, fit.x0D, fit.y0D, fit.bxD, fit.byD, true);
                r.impact  = ImpactToAxis(fit.x0, fit.y0, fit.bx, fit.by);
                r.impactD = ImpactToAxis(fit.x0D, fit.y0D, fit.bxD, fit.byD);
                r.r2  = RSquared(t, fit.x0, fit.y0, fit.bx, fit.by, false);
                r.r2D = RSquared(t, fit.x0D, fit.y0D, fit.bxD, fit.byD, true);

                r.passResD = PassesBlockResidual(r.resDRms) ? 1 : 0; // block plots only, see QualityCut.h
                r.passQuality = PassesQuality(r.resRms, r.impact, r.impactD);

                r.nProtonsEv = nProtonsInEvent[r.eventID];
                {
                    std::vector<std::size_t> blockHits; // tracker hits by z
                    for (std::size_t i = 0; i < t.size(); i++) {if (t.zD[i] != 0.0) {blockHits.push_back(i);}}
                    std::sort(blockHits.begin(), blockHits.end(), [&](std::size_t a, std::size_t b){return t.zD[a] < t.zD[b];});
                    r.nHits = (int)blockHits.size();
                    r.firstXD = r.nHits ? t.xD[blockHits.front()] : NAN; r.firstYD = r.nHits ? t.yD[blockHits.front()] : NAN;
                    r.resHalfRatioD = NAN;
                    if (r.nHits > 2) {
                        const Track first = t.select({blockHits.begin(), blockHits.begin() + 2}), second = t.select({blockHits.begin() + 2, blockHits.end()});
                        const double r1 = ResidualRms(first, fit.x0D, fit.y0D, fit.bxD, fit.byD, true), r2h = ResidualRms(second, fit.x0D, fit.y0D, fit.bxD, fit.byD, true);
                        if (r1 > 0) {r.resHalfRatioD = r2h / r1;}
                    }
                }
                rows.push_back(r);
            }
        }, ROOT::TSeqU(nShards));
        slotTracks.clear();

        // TTree::Fill is not thread-safe: the branch variables are set and the trees filled one row at a time.
        for (auto &rows : shardRows) {
            for (const Row &r : rows) {
                outTheta = r.theta; outPhi = r.phi; outDTheta = r.dTheta; outDPhi = r.dPhi;
                outEventID = r.eventID; outTrackID = r.trackID; outPDG = r.pdg;
                xVertexOut = r.xVertex; yVertexOut = r.yVertex; zVertexOut = r.zVertex;
                x0Out = r.x0;   y0Out = r.y0;   bxOut = r.bx;   byOut = r.by;
                x0DOut = r.x0D; y0DOut = r.y0D; bxDOut = r.bxD; byDOut = r.byD;
                finalKE = r.finalKE; finalKED = r.finalKED; fParentID = r.parentID;
                z_bx_reconVertexOut = r.zBxRecon; z_by_reconVertexOut = r.zByRecon;
                z_bx_reconVertexDIFFOut = r.zBxReconDIFF; z_by_reconVertexDIFFOut = r.zByReconDIFF;
                isFromCE = r.isFromCE;
                resRms = r.resRms; resDRms = r.resDRms; impact = r.impact; impactD = r.impactD; r2 = r.r2; r2D = r.r2D;
                passResD = r.passResD; nProtonsEv = r.nProtonsEv; zOrderOK = r.zOrderOK; nHits = r.nHits;
                firstXD = r.firstXD; firstYD = r.firstYD; resHalfRatioD = r.resHalfRatioD;

                passQuality = r.zOrderOK && r.passQuality ? 1 : 0; // out of order tracks and tracks that fail the quality cut stay in the tree for the plots, the views drop them
                inCE = r.inCE;
                outTree->Fill();
            }
            rows.clear(); rows.shrink_to_fit();
        }
    }

    outFile->cd();
    outTree->Write();
    outFile->Close();
}
