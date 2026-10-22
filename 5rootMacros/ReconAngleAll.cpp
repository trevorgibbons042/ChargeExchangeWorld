#include <TFile.h>
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

struct Pos {double x, y, z; double xD, yD, zD; double KE; int ParentID; double KED;}; // KE: true KE, KED: detector KE
struct FitResult {double x0, y0, bx, by;
    double x0D, y0D, bxD, byD;};
struct Track {std::vector<Pos> points;};

FitResult fitTrack(const std::vector<Pos>& points){
    FitResult fit{NAN, NAN, NAN, NAN,
                NAN, NAN, NAN, NAN};

    const std::size_t N = points.size();
    if (N < 2) {return fit;}

    double meanX = 0.0, meanY = 0.0, meanZ = 0.0;
    for (const auto& p : points) {
        meanX += p.x;
        meanY += p.y;
        meanZ += p.z;
    }
    meanX /= N;
    meanY /= N;
    meanZ /= N;

    double Szz = 0.0, Szx = 0.0, Szy = 0.0;
    for (const auto& p : points) {
        double dz = p.z - meanZ;
        Szz += dz * dz;
        Szx += dz * (p.x - meanX);
        Szy += dz * (p.y - meanY);
    }

    if (std::abs(Szz) >= 1e-20) {
        fit.bx = Szx / Szz;
        fit.by = Szy / Szz;
        fit.x0 = meanX - fit.bx * meanZ;
        fit.y0 = meanY - fit.by * meanZ;
    }

    double meanXD = 0.0, meanYD = 0.0, meanZD = 0.0;
    std::size_t ND = 0;
    for (const auto& p : points) {
        if (p.zD == 0.0) {continue;}
        meanXD += p.xD;
        meanYD += p.yD;
        meanZD += p.zD;
        ND++;
    }

    if (ND < 2) {return fit;}

    meanXD /= ND;
    meanYD /= ND;
    meanZD /= ND;

    double SzzD = 0.0, SzxD = 0.0, SzyD = 0.0;
    for (const auto& p : points) {
        if (p.zD == 0.0) {continue;}
        double dzD = p.zD - meanZD;
        SzzD += dzD * dzD;
        SzxD += dzD * (p.xD - meanXD);
        SzyD += dzD * (p.yD - meanYD);
    }

    if (std::abs(SzzD) >= 1e-20) {
        fit.bxD = SzxD / SzzD;
        fit.byD = SzyD / SzzD;
        fit.x0D = meanXD - fit.bxD * meanZD;
        fit.y0D = meanYD - fit.byD * meanZD;
    }

    return fit;
}

// One pass fits every proton and writes the true-KE trees (Angles, AnglesCE) and the detector-KE trees (Angles_DE, AnglesCE_DE); they differ only in finalKE.
void ReconAngleAll(){
    ROOT::EnableImplicitMT(128);

    ROOT::RDataFrame df("ParticleInfo", "output0_t*.root");

    // One container per thread slot, merged after each loop, so the entry loops need no global lock.
    std::vector<std::set<std::pair<int, int>>> slotReachedEnd(df.GetNSlots());
    df.Filter("boundaryType == 2 && pdg == 2212")
      .ForeachSlot([&](unsigned slot, int eventID, int trackID){
        slotReachedEnd[slot].insert({eventID, trackID});
      }, {"fEvent", "TrackID"});
    std::set<std::pair<int, int>> reachedEnd;
    for (auto &s : slotReachedEnd) reachedEnd.insert(s.begin(), s.end());

    // CE neutrons (same selection as CutEvents.cpp): event -> CE track IDs.
    ROOT::RDataFrame dfCE("ChargeExchange_values", "output0_t*.root");
    std::vector<std::map<int, std::set<int>>> slotCeTracks(dfCE.GetNSlots());
    dfCE.Filter("didChargeExchange>0 && fZ<7000")
        .ForeachSlot([&](unsigned slot, int eventID, int trackID){
        slotCeTracks[slot][eventID].insert(trackID);
    }, {"fEvent", "TrackID"});
    std::map<int, std::set<int>> ceTracks;
    for (auto &m : slotCeTracks) {
        for (auto &[ev, ids] : m) ceTracks[ev].insert(ids.begin(), ids.end());
    }

    std::vector<std::map<std::tuple<int, int, int>, Track>> slotTracks(df.GetNSlots());
    df.Filter("pdg == 2212").ForeachSlot([&](unsigned slot,
        int eventID, int trackID, int pdg,
        double x, double y, double z, double xD, double yD, double zD, double KE, double KED, int ParentID){
        if (!reachedEnd.count({eventID, trackID})) {return;}
        slotTracks[slot][{eventID, trackID, pdg}].points.push_back({x, y, z, xD, yD, zD, KE, ParentID, KED});
    }, {"fEvent", "TrackID", "pdg", "fX", "fY", "fZ", "fDx", "fDy", "fDz", "fKEnergy", "fDEnergy", "ParentID"});
    // Points are sorted and de-duplicated per track below, so the merge order does not matter.
    std::map<std::tuple<int, int, int>, Track> tracks;
    for (auto &slotMap : slotTracks) {
        for (auto &[key, t] : slotMap) {
            auto &dst = tracks[key].points;
            dst.insert(dst.end(), t.points.begin(), t.points.end());
        }
    }

    TFile *outFile = new TFile("cut_angles_all.root", "RECREATE");
    TTree *outTree = new TTree("Angles", "Reconstructed Angles, all protons");
    TTree *outTreeCE = new TTree("AnglesCE", "Reconstructed Angles, all protons in CE events");
    TTree *outTreeDE = new TTree("Angles_DE", "Reconstructed Angles, all protons, detector KE");
    TTree *outTreeCEDE = new TTree("AnglesCE_DE", "Reconstructed Angles, all protons in CE events, detector KE");

    int outTrackID, outPDG, outEventID, fParentID;
    double outTheta, outPhi, outDTheta, outDPhi;
    double x0Out, y0Out, bxOut, byOut;
    double x0DOut, y0DOut, bxDOut, byDOut;
    double xVertexOut, yVertexOut, zVertexOut;
    double z_bx_reconVertexOut, z_by_reconVertexOut;
    double z_bx_reconVertexDIFFOut, z_by_reconVertexDIFFOut;
    double finalKE, finalKED; // true and detector KE of the track end
    int isFromCE;

    for (TTree *tree : {outTree, outTreeCE, outTreeDE, outTreeCEDE}) {
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
    }

    for (auto &[key, t] : tracks) {
        std::sort(t.points.begin(), t.points.end(), [](const Pos &a, const Pos &b){
            if (a.z != b.z) return a.z < b.z;
            if (a.x != b.x) return a.x < b.x;
            return a.y < b.y;
        });

        t.points.erase(std::unique(t.points.begin(), t.points.end(), [](const Pos &a, const Pos &b){
            return a.x  == b.x  && a.y  == b.y  && a.z  == b.z &&
                   a.xD == b.xD && a.yD == b.yD && a.zD == b.zD;
        }), t.points.end());

        if (t.points.size() < 2) {continue;}

        auto byZ = [](const Pos &a, const Pos &b){return a.z < b.z;};
        auto endPoint   = std::max_element(t.points.begin(), t.points.end(), byZ);
        auto startPoint = std::min_element(t.points.begin(), t.points.end(), byZ);

        FitResult fit = fitTrack(t.points);
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

        xVertexOut = startPoint->x;
        yVertexOut = startPoint->y;
        zVertexOut = startPoint->z;

        x0Out = fit.x0;   y0Out = fit.y0;   bxOut = fit.bx;   byOut = fit.by;
        x0DOut = fit.x0D; y0DOut = fit.y0D; bxDOut = fit.bxD; byDOut = fit.byD;

        finalKE   = endPoint->KE;
        finalKED  = endPoint->KED;
        fParentID = endPoint->ParentID;

        z_bx_reconVertexOut = (fit.bxD != 0) ? -fit.x0D/fit.bxD : 0;
        z_by_reconVertexOut = (fit.byD != 0) ? -fit.y0D/fit.byD : 0;
        z_bx_reconVertexDIFFOut = zVertexOut - z_bx_reconVertexOut;
        z_by_reconVertexDIFFOut = zVertexOut - z_by_reconVertexOut;

        auto ce = ceTracks.find(outEventID);
        isFromCE = (ce != ceTracks.end() && ce->second.count(fParentID)) ? 1 : 0;

        outTree->Fill();
        outTreeDE->Fill();
        if (ce != ceTracks.end()) {outTreeCE->Fill(); outTreeCEDE->Fill();}
    }

    outFile->cd();
    outTree->Write();
    outTreeCE->Write();
    outTreeDE->Write();
    outTreeCEDE->Write();
    outFile->Close();
}
