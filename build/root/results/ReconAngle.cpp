#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <memory>
#include <filesystem>
#include <sstream>
#include <TStyle.h>
#include <limits>
#include <map>
#include <mutex>
#include <cmath>
#include <tuple>
#include <vector>

struct Pos {double x, y, z;
    double xD, yD, zD;};
struct FitResult {double x0, y0, bx, by;
    double x0D, y0D, bxD, byD;};
struct Track {std::vector<Pos> points;};

FitResult fitTrack(const std::vector<Pos>& points){
    FitResult fit;
    const std::size_t N = points.size();

    if (N < 2) {
        fit.x0 = NAN;
        fit.y0 = NAN;
        fit.bx = NAN;
        fit.by = NAN;
        fit.x0D = NAN;
        fit.y0D = NAN;
        fit.bxD = NAN;
        fit.byD = NAN;
        return fit;
    }

    double meanX = 0.0;
    double meanY = 0.0;
    double meanZ = 0.0;
    double meanXD = 0.0;
    double meanYD = 0.0;
    double meanZD = 0.0;

    for (const auto& p : points) {
        meanX += p.x;
        meanY += p.y;
        meanZ += p.z;
        meanXD += p.xD;
        meanYD += p.yD;
        meanZD += p.zD;
    }

    meanX /= N;
    meanY /= N;
    meanZ /= N;
    meanXD /= N;
    meanYD /= N;
    meanZD /= N;

    double Szz = 0.0;
    double Szx = 0.0;
    double Szy = 0.0;
    double SzzD = 0.0;
    double SzxD = 0.0;
    double SzyD = 0.0;


    for (const auto& p : points) {
        const double dz = p.z - meanZ;
        const double dx = p.x - meanX;
        const double dy = p.y - meanY;
        Szz += dz * dz;
        Szx += dz * dx;
        Szy += dz * dy;
        
        double dzD = p.zD - meanZD;
        double dxD = p.xD - meanXD;
        double dyD = p.yD - meanYD;

        SzzD += dzD * dzD;
        SzxD += dzD * dxD;
        SzyD += dzD * dyD;
    }

    if (std::abs(Szz) < 1e-20) {
        fit.x0 = NAN;
        fit.y0 = NAN;
        fit.bx = NAN;
        fit.by = NAN;
        return fit;
    } else{
        fit.bx = Szx / Szz;
        fit.by = Szy / Szz;

        fit.x0 = meanX - fit.bx * meanZ;
        fit.y0 = meanY - fit.by * meanZ;
    }
    
    if (std::abs(SzzD) < 1e-20) {
        fit.x0D = NAN;
        fit.y0D = NAN;
        fit.bxD = NAN;
        fit.byD = NAN;
    }
    else {
        fit.bxD = SzxD / SzzD;
        fit.byD = SzyD / SzzD;

        fit.x0D = meanXD - fit.bxD * meanZD;
        fit.y0D = meanYD - fit.byD * meanZD;
    }

    return fit;
}

struct CalcValPrint {
    int eventID;
    int trackID;
    int pdg;
    double theta;
    double phi;
};


void ReconAngle(){
    ROOT::EnableImplicitMT();
    
    bool Verbose = false;
    
    ROOT::RDataFrame df("ParticleInfo", "cut_EOL.root");
    TFile *outFile = new TFile("cut_angles.root", "RECREATE");
    TTree *outTree = new TTree("Angles", "Reconstructed Angles");

    int outTrackID;
    int outPDG;
    int outEventID;
    
    double outTheta;
    double outPhi;
    double outDTheta;
    double outDPhi;
    
    double x0Out;
    double y0Out;
    double bxOut;
    double byOut;
    
    double x0DOut;
    double y0DOut;
    double bxDOut;
    double byDOut;
    
    outTree->Branch("TrackID", &outTrackID);
    outTree->Branch("PDG", &outPDG);
    outTree->Branch("fEvent", &outEventID);
    outTree->Branch("outTheta", &outTheta);
    outTree->Branch("outPhi", &outPhi);
    outTree->Branch("outDTheta", &outDTheta);
    outTree->Branch("outDPhi", &outDPhi);
    outTree->Branch("x0Out", &x0Out);
    outTree->Branch("y0Out", &y0Out);
    outTree->Branch("bxOut", &bxOut);
    outTree->Branch("byOut", &byOut);
    outTree->Branch("x0DOut", &x0DOut);
    outTree->Branch("y0DOut", &y0DOut);
    outTree->Branch("bxDOut", &bxDOut);
    outTree->Branch("byDOut", &byDOut);
    
    
    std::map<std::tuple<int, int, int>, Track> tracks;
    std::mutex m;
    df.Foreach([&](
        int eventID, int trackID, int pdg, int boundary,
        double x, double y, double z, double xD, double yD, double zD,
        ULong64_t entry){
        std::lock_guard<std::mutex> lock(m);
        
        auto &t = tracks[{eventID, trackID, pdg}];
        t.points.push_back({x, y, z, xD, yD, zD});
        },{"fEvent","TrackID", "pdg", "boundaryType", "fX", "fY", "fZ", "fDx", "fDy", "fDz", "rdfentry_"});

    for (auto &[key, t] : tracks) {
        std::sort(
            t.points.begin(),
            t.points.end(),
            [](const Pos &a, const Pos &b)
            {
                if (a.z != b.z) return a.z < b.z;
                if (a.x != b.x) return a.x < b.x;
                return a.y < b.y;
            }
        );

        t.points.erase(std::unique(t.points.begin(),t.points.end(), [](const Pos &a, const Pos &b){
            return a.x  == b.x  &&
                   a.y  == b.y  &&
                   a.z  == b.z  &&
                   a.xD == b.xD &&
                   a.yD == b.yD &&
                   a.zD == b.zD;
        }), t.points.end());
        
        if (t.points.size() < 2){continue;}
        
        auto endPoint = std::max_element(t.points.begin(),t.points.end(),[](const Pos &a, const Pos &b){
                return a.z < b.z;
            });

        FitResult fit = fitTrack(t.points);
        double zEnd = endPoint->z;
        double xEnd = fit.x0 + fit.bx * zEnd;
        double yEnd = fit.y0 + fit.by * zEnd;
        
        
        if (!std::isfinite(fit.bx) || !std::isfinite(fit.by)){continue;}
        if (!std::isfinite(fit.bxD) || !std::isfinite(fit.byD)){continue;}
        
        double theta = std::atan2(std::hypot(xEnd, yEnd),zEnd);
        double phi = std::atan2(yEnd, xEnd);
        if (phi < 0) {phi += 2.0 * M_PI;}
        
        outEventID = std::get<0>(key);
        outTrackID = std::get<1>(key);
        outPDG = std::get<2>(key);
        outTheta = theta;
        outPhi = phi;
        x0Out = fit.x0;
        y0Out = fit.y0;
        bxOut = fit.bx;
        byOut = fit.by;
        
        
        auto endPointD = std::max_element(t.points.begin(), t.points.end(),[](const Pos &a, const Pos &b){
                return a.zD < b.zD;
            });
        
        double zEndD = endPointD->zD;
        double xEndD = fit.x0D + fit.bxD * zEndD;
        double yEndD = fit.y0D + fit.byD * zEndD;
        
        double thetaD = std::atan2(std::hypot(xEndD, yEndD), zEndD);
        double phiD = std::atan2(yEndD, xEndD);
        if (phiD < 0) {phiD += 2.0 * M_PI;}
        
        outDTheta = thetaD;
        outDPhi = phiD;
        x0DOut = fit.x0D;
        y0DOut = fit.y0D;
        bxDOut = fit.bxD;
        byDOut = fit.byD;
        
        if (Verbose == true){
            for (size_t i = 0; i < t.points.size(); i++) {
                std::cout
                    << "Point " << i
                    << ": x = " << t.points[i].x
                    << ", y = " << t.points[i].y
                    << ", z = " << t.points[i].z
                    << " | xD = " << t.points[i].xD
                    << ", yD = " << t.points[i].yD
                    << ", zD = " << t.points[i].zD
                    << std::endl;
                std::cout
                    << "Event: "   << std::get<0>(key)
                    << " TrackID: " << std::get<1>(key)
                    << " PDG: "     << std::get<2>(key)
                    << "\n";
                
                std::cout
                    << "theta: "   << theta
                    << " phi: " << phi
                    << " thetaD: " << thetaD
                    << " phiD: " << phiD
                    << "\n";
            }
        }
        
        outTree->Fill();
    }
    outFile->cd();
    outTree->Write();
    outFile->Close();

    ROOT::RDataFrame dfCV("CalcVals", "cut_EOL.root");
    std::vector<CalcValPrint> calcValues;
    std::mutex mCV;

    dfCV.Foreach([&](int eventID,int trackID,int pdg, double theta, double phi){
    std::lock_guard<std::mutex> lock(mCV);
    calcValues.push_back({eventID,trackID,pdg,theta,phi});
    },{"fEvent","TrackID","pdg","thetaCalcXY0detector0Pos","phiCalcXY0detector0Pos"});

    std::sort(
        calcValues.begin(),
        calcValues.end(),
        [](const CalcValPrint &a, const CalcValPrint &b)
        {
            if (a.eventID != b.eventID)
                return a.eventID < b.eventID;

            if (a.trackID != b.trackID)
                return a.trackID < b.trackID;

            return a.pdg < b.pdg;
        }
    );

    if (Verbose == true) {
        for (const auto &v : calcValues) {
            std::cout
                << "Event: " << v.eventID
                << " TrackID: " << v.trackID
                << " PDG: " << v.pdg
                << std::endl;
            std::cout
                << "thetaCalc: " << v.theta
                << ", phiCalc: " << v.phi
                << std::endl;
        }
    }

}
