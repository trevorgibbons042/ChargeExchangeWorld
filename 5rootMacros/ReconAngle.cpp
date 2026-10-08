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

#include <TCanvas.h>
#include <TGraph2D.h>
#include <TPolyLine3D.h>

struct Pos {double x, y, z; double xD, yD, zD; double KE; int ParentID;};
struct FitResult {double x0, y0, bx, by;
    double x0D, y0D, bxD, byD;};
struct Track {std::vector<Pos> points;};

FitResult fitTrack(const std::vector<Pos>& points){
    FitResult fit{NAN, NAN, NAN, NAN,
                NAN, NAN, NAN, NAN};

    const std::size_t N = points.size();

    if (N < 2) {return fit;}

    double meanX = 0.0;
    double meanY = 0.0;
    double meanZ = 0.0;

    for (const auto& p : points) {
        meanX += p.x;
        meanY += p.y;
        meanZ += p.z;
    }

    meanX /= N;
    meanY /= N;
    meanZ /= N;

    double Szz = 0.0;
    double Szx = 0.0;
    double Szy = 0.0;

    for (const auto& p : points) {

        double dz = p.z - meanZ;
        double dx = p.x - meanX;
        double dy = p.y - meanY;

        Szz += dz * dz;
        Szx += dz * dx;
        Szy += dz * dy;
    }

    if (std::abs(Szz) >= 1e-20) {

        fit.bx = Szx / Szz;
        fit.by = Szy / Szz;

        fit.x0 = meanX - fit.bx * meanZ;
        fit.y0 = meanY - fit.by * meanZ;
    }

    double meanXD = 0.0;
    double meanYD = 0.0;
    double meanZD = 0.0;

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

    double SzzD = 0.0;
    double SzxD = 0.0;
    double SzyD = 0.0;

    for (const auto& p : points) {
        if (p.zD == 0.0) {continue;}

        double dzD = p.zD - meanZD;
        double dxD = p.xD - meanXD;
        double dyD = p.yD - meanYD;

        SzzD += dzD * dzD;
        SzxD += dzD * dxD;
        SzyD += dzD * dyD;
    }

    if (std::abs(SzzD) >= 1e-20) {

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


void PlotFit3D(
    const std::vector<Pos>& points,
    const FitResult& fit,
    int eventID,
    int trackID
);

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
    
    double finalKE;
    int fParentID;
    
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
    outTree->Branch("finalKE", &finalKE);
    outTree->Branch("ParentID", &fParentID);
    
    
    
    std::map<std::tuple<int, int, int>, Track> tracks;
    std::mutex m;
    df.Foreach([&](
        int eventID, int trackID, int pdg, int boundary,
        double x, double y, double z, double xD, double yD, double zD, double KE, int ParentID, ULong64_t entry){
        std::lock_guard<std::mutex> lock(m);
        
        auto &t = tracks[{eventID, trackID, pdg}];
        t.points.push_back({x, y, z, xD, yD, zD, KE, ParentID});
        },{"fEvent","TrackID", "pdg", "boundaryType", "fX", "fY", "fZ", "fDx", "fDy", "fDz", "fKEnergy", "ParentID", "rdfentry_"});

    bool plotted = false;

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
        
        auto startPoint = std::min_element(t.points.begin(),t.points.end(),[](const Pos &a, const Pos &b){
                return a.z < b.z;
            });

        FitResult fit = fitTrack(t.points);
        
        /*
        double zEnd = endPoint->z;
        double zStart = startPoint->z;
        double dz = zEnd - zStart;
        
        double xStart = fit.x0 + fit.bx * zStart;
        double xEnd = fit.x0 + fit.bx * zEnd;
        double dx = xEnd - xStart;
        
        double yStart = fit.y0 + fit.by * zStart;
        double yEnd = fit.y0 + fit.by * zEnd;
        double dy = yEnd - yStart;
         
        double theta = std::atan2(std::hypot(dx, dy),dz);
        double phi = std::atan2(dy, dx);
        */
        
        if (!std::isfinite(fit.bx) || !std::isfinite(fit.by)){continue;}
        if (!std::isfinite(fit.bxD) || !std::isfinite(fit.byD)){continue;}
        
        double theta = std::atan2(std::hypot(fit.bx, fit.by), 1.0);
        double phi = std::atan2(fit.by, fit.bx);
        if (phi < 0) {phi += 2.0 * M_PI;}
        
        int pdg = std::get<2>(key);
        
        int eventID = std::get<0>(key);
        int trackID = std::get<1>(key);

        if (!plotted && pdg == 2212) {
            std::cout
                << "PLOTTING event " << eventID
                << " track " << trackID
                << " with " << t.points.size()
                << " points"
                << std::endl;

            PlotFit3D(t.points,fit,eventID,trackID);
            plotted = true;
        }
        
        if (Verbose == true){
            if (theta > 0.01 && pdg == 2212) {
                double zEnd = endPoint->z;
                double zStart = startPoint->z;
                double yEnd = endPoint->y;
                double yStart = startPoint->y;
                double xEnd = endPoint->x;
                double xStart = startPoint->x;
                int eventID = std::get<0>(key);
                int trackID = std::get<1>(key);
                
                std::cout
                    << "Event = " << std::get<0>(key)
                    << " Track = " << std::get<1>(key)
                    << " pdg = " << std::get<2>(key)
                    << " xStart = " << xStart << " yStart = " << yStart << " zStart = " << zStart
                    << " xEnd = " << xEnd << " yEnd = " << yEnd << " zEnd = " << zEnd
                    << " theta = " << theta
                    << std::endl;
            }
        }
        
        
        
        
        outEventID = std::get<0>(key);
        outTrackID = std::get<1>(key);
        outPDG = std::get<2>(key);
        outTheta = theta;
        outPhi = phi;
        x0Out = fit.x0;
        y0Out = fit.y0;
        bxOut = fit.bx;
        byOut = fit.by;
        
        finalKE = endPoint->KE;
        fParentID = endPoint->ParentID;
        
        
        auto endPointD = std::max_element(t.points.begin(), t.points.end(),[](const Pos &a, const Pos &b){
                return a.zD < b.zD;
            });
        
        double zEndD = endPointD->zD;
        double xEndD = fit.x0D + fit.bxD * zEndD;
        double yEndD = fit.y0D + fit.byD * zEndD;
                
        double thetaD = std::atan2(std::hypot(fit.bxD, fit.byD), 1.0);
        double phiD = std::atan2(fit.byD, fit.bxD);
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
}

void PlotFit3D(
    const std::vector<Pos>& points,
    const FitResult& fit,
    int eventID,
    int trackID){
    std::vector<double> x, y, z;
    std::vector<double> xD, yD, zD;

    double zMin  =  1e99;
    double zMax  = -1e99;
    double zDMin =  1e99;
    double zDMax = -1e99;

    for (const auto& p : points) {

        x.push_back(p.x);
        y.push_back(p.y);
        z.push_back(p.z);

        zMin = std::min(zMin, p.z);
        zMax = std::max(zMax, p.z);

        if (p.zD != 0.0) {

            xD.push_back(p.xD);
            yD.push_back(p.yD);
            zD.push_back(p.zD);

            zDMin = std::min(zDMin, p.zD);
            zDMax = std::max(zDMax, p.zD);
        }
    }

    if (x.size() < 2 || xD.size() < 2)
        return;


    auto* c = new TCanvas(
        Form("fit_%d_%d", eventID, trackID),
        Form("Event %d Track %d", eventID, trackID),
        1400,
        700
    );

    c->Divide(2,1);


    // =====================================================
    // X(Z)
    // =====================================================

    c->cd(1);

    auto* gX = new TGraph(z.size());

    for (size_t i = 0; i < z.size(); ++i) {
        gX->SetPoint(i, z[i], x[i]);
    }

    gX->SetTitle(
        Form(
            "Event %d Track %d : x(z);z;x",
            eventID,
            trackID
        )
    );

    gX->SetMarkerStyle(20);
    gX->SetMarkerSize(0.7);

    gX->GetXaxis()->SetLimits(0,26000);
    gX->SetMinimum(-100);
    gX->SetMaximum(100);

    gX->Draw("AP");


    // detector points
    auto* gXD = new TGraph(zD.size());

    for (size_t i = 0; i < zD.size(); ++i) {
        gXD->SetPoint(i, zD[i], xD[i]);
    }

    gXD->SetMarkerStyle(24);
    gXD->SetMarkerSize(0.7);
    gXD->Draw("P SAME");


    // real-position fit
    auto* lineX = new TGraph(2);

    lineX->SetPoint(
        0,
        zMin,
        fit.x0 + fit.bx*zMin
    );

    lineX->SetPoint(
        1,
        zMax,
        fit.x0 + fit.bx*zMax
    );

    lineX->SetLineWidth(3);

    // dashed real fit
    lineX->SetLineStyle(2);

    lineX->Draw("L SAME");


    // detector fit
    auto* lineXD = new TGraph(2);

    lineXD->SetPoint(
        0,
        zDMin,
        fit.x0D + fit.bxD*zDMin
    );

    lineXD->SetPoint(
        1,
        zDMax,
        fit.x0D + fit.bxD*zDMax
    );

    lineXD->SetLineWidth(3);
    lineXD->Draw("L SAME");


    // =====================================================
    // Y(Z)
    // =====================================================

    c->cd(2);

    auto* gY = new TGraph(z.size());

    for (size_t i = 0; i < z.size(); ++i) {
        gY->SetPoint(i, z[i], y[i]);
    }

    gY->SetTitle(
        Form(
            "Event %d Track %d : y(z);z;y",
            eventID,
            trackID
        )
    );

    gY->SetMarkerStyle(20);
    gY->SetMarkerSize(0.7);

    gY->GetXaxis()->SetLimits(0,26000);
    gY->SetMinimum(-100);
    gY->SetMaximum(100);

    gY->Draw("AP");


    // detector points
    auto* gYD = new TGraph(zD.size());

    for (size_t i = 0; i < zD.size(); ++i) {
        gYD->SetPoint(i, zD[i], yD[i]);
    }

    gYD->SetMarkerStyle(24);
    gYD->SetMarkerSize(0.7);
    gYD->Draw("P SAME");


    // real-position fit
    auto* lineY = new TGraph(2);

    lineY->SetPoint(
        0,
        zMin,
        fit.y0 + fit.by*zMin
    );

    lineY->SetPoint(
        1,
        zMax,
        fit.y0 + fit.by*zMax
    );

    lineY->SetLineWidth(3);

    // dashed real fit
    lineY->SetLineStyle(2);

    lineY->Draw("L SAME");


    // detector fit
    auto* lineYD = new TGraph(2);

    lineYD->SetPoint(
        0,
        zDMin,
        fit.y0D + fit.byD*zDMin
    );

    lineYD->SetPoint(
        1,
        zDMax,
        fit.y0D + fit.byD*zDMax
    );

    lineYD->SetLineWidth(3);
    lineYD->Draw("L SAME");


    c->Modified();
    c->Update();
}
