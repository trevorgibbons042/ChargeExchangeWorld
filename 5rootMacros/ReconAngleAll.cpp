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

// One pass fits every proton and writes the true-KE trees (Angles, AnglesCE) and the detector-KE trees (Angles_DE, AnglesCE_DE); they differ only in finalKE.
void ReconAngleAll(){
    ROOT::EnableImplicitMT(128);

    // CE neutrons (same selection as CutEvents.cpp): event -> CE track IDs.
    ROOT::RDataFrame dfCE("ChargeExchange_values", "output0_t*.root");
    std::vector<std::map<int, std::set<int>>> slotCeTracks(dfCE.GetNSlots());
    dfCE.Filter(CarbonCeCut())
        .ForeachSlot([&](unsigned slot, int eventID, int trackID){
        slotCeTracks[slot][eventID].insert(trackID);
    }, {"fEvent", "TrackID"});
    std::map<int, std::set<int>> ceTracks;
    for (auto &m : slotCeTracks) {
        for (auto &[ev, ids] : m) ceTracks[ev].insert(ids.begin(), ids.end());
    }

    TFile *outFile = new TFile("cuts/cut_angles_all.root", "RECREATE");
    TTree *outTree = new TTree("Angles", "Reconstructed Angles, all recon protons");
    TTree *outTreeCE = new TTree("AnglesCE", "Reconstructed Angles, recon protons in CE events");
    TTree *outTreeDE = new TTree("Angles_DE", "Reconstructed Angles, all recon protons, detector KE");
    TTree *outTreeCEDE = new TTree("AnglesCE_DE", "Reconstructed Angles, recon protons in detector KE");
    TTree *outTreeRaw = new TTree("AnglesRaw", "Reconstructed Angles, every recon proton before the quality cut (QualityStudy.cpp)");

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

    for (TTree *tree : {outTree, outTreeCE, outTreeDE, outTreeCEDE, outTreeRaw}) {
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
        tree->Branch("finalKE", (tree == outTreeDE || tree == outTreeCEDE) ? &finalKED : &finalKE);
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
        df.Filter("pdg == 2212 && boundaryType != 2").ForeachSlot([&](unsigned slot,
            int eventID, int trackID, int pdg,
            double x, double y, double z, double xD, double yD, double zD, double KE, double KED, int ParentID, ULong64_t entry){
            slotTracks[slot][{eventID, trackID, pdg}].push(x, y, z, xD, yD, zD, KE, ParentID, KED, entry);
        }, {"fEvent", "TrackID", "pdg", "fX", "fY", "fZ", "fDx", "fDy", "fDz", "fKEnergy", "fDEnergy", "ParentID", "rdfentry_"});
        // Points are sorted and de-duplicated per track below, so the merge order does not matter.
        std::map<std::tuple<int, int, int>, Track> tracks;
        for (auto &slotMap : slotTracks) {
            for (auto &[key, t] : slotMap) tracks[key].append(t);
        }

        std::map<int, int> nProtonsInEvent;
        for (const auto &kv : tracks) nProtonsInEvent[std::get<0>(kv.first)]++;

        for (auto &kv : tracks) {
            const auto &key = kv.first;
            Track &t = kv.second; // not a structured binding: the sort lambdas below capture it, which C++17 does not allow for a binding
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
                zOrderOK = outOfOrder ? 0 : 1; // out of order tracks stay in AnglesRaw for the plots, the other trees drop them below
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

            outTheta = std::atan2(std::hypot(fit.bx, fit.by), 1.0);
            outPhi   = std::atan2(fit.by, fit.bx);
            if (outPhi < 0) {outPhi += 2.0 * M_PI;}

            outDTheta = std::atan2(std::hypot(fit.bxD, fit.byD), 1.0);
            outDPhi   = std::atan2(fit.byD, fit.bxD);
            if (outDPhi < 0) {outDPhi += 2.0 * M_PI;}

            outEventID = std::get<0>(key);
            outTrackID = std::get<1>(key);
            outPDG     = std::get<2>(key);

            xVertexOut = t.x[startPoint];
            yVertexOut = t.y[startPoint];
            zVertexOut = t.z[startPoint];

            x0Out = fit.x0;   y0Out = fit.y0;   bxOut = fit.bx;   byOut = fit.by;
            x0DOut = fit.x0D; y0DOut = fit.y0D; bxDOut = fit.bxD; byDOut = fit.byD;

            finalKE   = t.KE[endPoint];
            finalKED  = t.KED[endPoint];
            fParentID = t.ParentID[endPoint];

            z_bx_reconVertexOut = (fit.bxD != 0) ? -fit.x0D/fit.bxD : 0;
            z_by_reconVertexOut = (fit.byD != 0) ? -fit.y0D/fit.byD : 0;
            z_bx_reconVertexDIFFOut = zVertexOut - z_bx_reconVertexOut;
            z_by_reconVertexDIFFOut = zVertexOut - z_by_reconVertexOut;

            auto ce = ceTracks.find(outEventID);
            isFromCE = (ce != ceTracks.end() && ce->second.count(fParentID)) ? 1 : 0;

            resRms  = ResidualRms(t, fit.x0, fit.y0, fit.bx, fit.by, false);
            resDRms = ResidualRms(t, fit.x0D, fit.y0D, fit.bxD, fit.byD, true);
            impact  = ImpactToAxis(fit.x0, fit.y0, fit.bx, fit.by);
            impactD = ImpactToAxis(fit.x0D, fit.y0D, fit.bxD, fit.byD);
            r2  = RSquared(t, fit.x0, fit.y0, fit.bx, fit.by, false);
            r2D = RSquared(t, fit.x0D, fit.y0D, fit.bxD, fit.byD, true);

            passResD = PassesBlockResidual(resDRms) ? 1 : 0; // block plots only, see QualityCut.h

            outEventID = std::get<0>(key);
            nProtonsEv = nProtonsInEvent[outEventID];
            {
                std::vector<std::size_t> blockHits; // tracker hits by z
                for (std::size_t i = 0; i < t.size(); i++) {if (t.zD[i] != 0.0) {blockHits.push_back(i);}}
                std::sort(blockHits.begin(), blockHits.end(), [&](std::size_t a, std::size_t b){return t.zD[a] < t.zD[b];});
                nHits = (int)blockHits.size();
                firstXD = nHits ? t.xD[blockHits.front()] : NAN; firstYD = nHits ? t.yD[blockHits.front()] : NAN;
                resHalfRatioD = NAN;
                if (nHits > 2) {
                    const Track first = t.select({blockHits.begin(), blockHits.begin() + 2}), second = t.select({blockHits.begin() + 2, blockHits.end()});
                    const double r1 = ResidualRms(first, fit.x0D, fit.y0D, fit.bxD, fit.byD, true), r2h = ResidualRms(second, fit.x0D, fit.y0D, fit.bxD, fit.byD, true);
                    if (r1 > 0) {resHalfRatioD = r2h / r1;}
                }
            }

            outTreeRaw->Fill();
            if (!zOrderOK) {continue;}
            if (!PassesQuality(resRms, impact, impactD)) {continue;}
            outTree->Fill();
            outTreeDE->Fill();
            if (ce != ceTracks.end()) {outTreeCE->Fill(); outTreeCEDE->Fill();}
        }
    }

    outFile->cd();
    outTree->Write();
    outTreeCE->Write();
    outTreeDE->Write();
    outTreeCEDE->Write();
    outTreeRaw->Write();
    outFile->Close();
}
