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


namespace
{
  const double AP_A = 0.7079;
  const double AP_k = 0.3860;
  const double z_end = 23000;
}

using Key = std::pair<int, int>;

struct Angles {double LastPhi, LastTheta, LastMomentum, LastT, LastAP; int PDG;};
struct Track {std::vector<Angles> points;};

void CalcAPandT(){
    ROOT::EnableImplicitMT();
    
    bool Verbose = false;
    
    ROOT::RDataFrame dfCEV("ChargeExchange_values", "cut_CE.root");
    ROOT::RDataFrame dfAngles("Angles", "cut_angles.root");
    
    TFile *outFile = new TFile("cut_APandT.root", "RECREATE");
    TTree *outTree = new TTree("CE", "values from ChargeExchange directly");
    TTree *outTree2 = new TTree("CalcValsRe", "Calculated AP and t from recon");
    TTree *outTree3 = new TTree("CalcValsReD", "Calculated AP and t from recon Dectector Blocks");
    TTree *outTree4 = new TTree("CalcValsReDIFF", "Difference Calculated AP and t from recon");
    TTree *outTree5 = new TTree("CalcValsReDDIFF", "Difference Calculated AP and t from recon Dectector Blocks");
    
    int outTrackID;
    int outPDG;
    int outEventID;
    int outType;
    
    double outThetaCE;
    double outPhiCE;
    double outMomCE;
    double outTCE;
    double outAPCE;
    double outFOMCE;
    
    double outThetaRe;
    double outPhiRe;
    double outMomRe;
    double outTRe;
    double outAPRe;
    double outFOMRe;
    
    double outThetaReD;
    double outPhiReD;
    double outMomReD;
    double outTReD;
    double outAPReD;
    double outFOMReD;
    
    double outThetaReDIFF;
    double outPhiReDIFF;
    double outTReDIFF;
    double outAPReDIFF;
    double outFOMReDIFF;
    
    double outThetaReDDIFF;
    double outPhiReDDIFF;
    double outTReDDIFF;
    double outAPReDDIFF;
    double outFOMReDDIFF;
    
    outTree->Branch("TrackID", &outTrackID);
    outTree->Branch("PDG", &outPDG);
    outTree->Branch("fEvent", &outEventID);
    outTree->Branch("fType", &outType);
    outTree->Branch("thetaCE", &outThetaCE);
    outTree->Branch("phiCE", &outPhiCE);
    outTree->Branch("tCE", &outTCE);
    outTree->Branch("apCE", &outAPCE);
    outTree->Branch("fomCE", &outFOMCE);
    outTree->Branch("momCE",  &outMomCE);
    
    outTree2->Branch("TrackID", &outTrackID);
    outTree2->Branch("PDG", &outPDG);
    outTree2->Branch("fEvent", &outEventID);
    outTree2->Branch("fType", &outType);
    outTree2->Branch("thetaRe", &outThetaRe);
    outTree2->Branch("phiRe", &outPhiRe);
    outTree2->Branch("tRe", &outTRe);
    outTree2->Branch("apRe", &outAPRe);
    outTree2->Branch("fomRe", &outFOMRe);
    outTree2->Branch("momRe",  &outMomRe);
    
    outTree3->Branch("TrackID", &outTrackID);
    outTree3->Branch("PDG", &outPDG);
    outTree3->Branch("fEvent", &outEventID);
    outTree3->Branch("fType", &outType);
    outTree3->Branch("thetaReD", &outThetaReD);
    outTree3->Branch("phiReD", &outPhiReD);
    outTree3->Branch("tReD", &outTReD);
    outTree3->Branch("apReD", &outAPReD);
    outTree3->Branch("fomReD", &outFOMReD);
    outTree3->Branch("momReD", &outMomReD);
    
    outTree4->Branch("TrackID", &outTrackID);
    outTree4->Branch("PDG", &outPDG);
    outTree4->Branch("fEvent", &outEventID);
    outTree4->Branch("fType", &outType);
    outTree4->Branch("thetaReD", &outThetaReDIFF);
    outTree4->Branch("phiReD", &outPhiReDIFF);
    outTree4->Branch("tReD", &outTReDIFF);
    outTree4->Branch("apReD", &outAPReDIFF);
    outTree4->Branch("fomReD", &outFOMReDIFF);
    
    outTree5->Branch("TrackID", &outTrackID);
    outTree5->Branch("PDG", &outPDG);
    outTree5->Branch("fEvent", &outEventID);
    outTree5->Branch("fType", &outType);
    outTree5->Branch("thetaReD", &outThetaReDDIFF);
    outTree5->Branch("phiReD", &outPhiReDDIFF);
    outTree5->Branch("tReD", &outTReDDIFF);
    outTree5->Branch("apReD", &outAPReDDIFF);
    outTree5->Branch("fomReD", &outFOMReDDIFF);
    
    
    std::map<std::tuple<int, int, int>, Track> tracksCE;
    std::mutex m;
    
    dfCEV.Foreach([&](
    int eventID, int trackID, double LastPhi, double LastTheta, double LastMomentum, double LastT, double LastAP, double z, ULong64_t entry){
        std::lock_guard<std::mutex> lock(m);
                          
        auto &t = tracksCE[{eventID, trackID, 0}];
    
        t.points.push_back({LastPhi, LastTheta, LastMomentum, LastT, LastAP, 2112});
    
     
    },{"fEvent", "TrackID", "LastPhi", "LastTheta", "LastMomentum", "LastT", "LastAP", "fZ", "rdfentry_"});

    std::map<std::tuple<int, int, int>, Track> tracksAngles;
    std::mutex mAngles;
    
    dfAngles.Foreach([&](
    int eventID, int ParentID, double outTheta, double outPhi, double outDTheta, double outDPhi, double finalKE, int pdg, ULong64_t entry){
        
        std::lock_guard<std::mutex> lock(mAngles);
        if (pdg != 2212) {return;}
         
        auto &tRecon = tracksAngles[{eventID, ParentID, 1}];
        
        double mom_recon = std::sqrt(finalKE*finalKE+2*finalKE*938);
        double t_recon = std::acos(1-outTheta);
        double AP_recon = AP_A*(1-std::exp(-(AP_k*mom_recon)))*std::sqrt(t_recon);
                     
        tRecon.points.push_back({outPhi, outTheta, mom_recon, t_recon, AP_recon, pdg});
        
        auto &tReconBlock = tracksAngles[{eventID, ParentID, 2}];
        double t_reconBlock = std::acos(1-outDTheta);
        double AP_reconBlock = (AP_A)*(1-std::exp(-(AP_k*mom_recon)))*std::sqrt(t_reconBlock);
        
        tReconBlock.points.push_back({outDPhi, outDTheta, mom_recon, t_reconBlock, AP_reconBlock});
        },{"fEvent", "ParentID","outTheta", "outPhi", "outDTheta", "outDPhi", "finalKE", "PDG","rdfentry_"});
    
    if (Verbose == true){
        for (auto &[key, tRecon] : tracksAngles) {
            for (const auto& p : tRecon.points) {
                if (std::get<2>(key) == 2){continue;}
                std::cout
                    << "Event = " << std::get<0>(key)
                    << " TrackID = " << std::get<1>(key)
                    << " Type = " << std::get<2>(key)
                    << " theta = " << p.LastTheta
                    << " pdg = " << p.PDG
                    << std::endl;
            }
        }
    }

    

    
    for (auto &[key, tCE] : tracksCE) {
        if (tCE.points.empty()) {continue;}

        int eventID = std::get<0>(key);
        int CETrackID = std::get<1>(key);
        
        auto itRecon = tracksAngles.find({eventID, CETrackID, 1});
        auto itReconD =tracksAngles.find({eventID, CETrackID, 2});
        
        if (itRecon == tracksAngles.end()) {continue;}
        if (itReconD == tracksAngles.end()) {continue;}
        auto &tRecon  = itRecon->second;
        auto &tReconD = itReconD->second;

        if (tRecon.points.empty() || tReconD.points.empty()) {continue;}

        const auto &ce = tCE.points.back();
        const auto &re = tRecon.points.back();
        const auto &reD = tReconD.points.back();

        outEventID = eventID;
        outTrackID = CETrackID;

        outPhiCE   = ce.LastPhi;
        outThetaCE = ce.LastTheta;
        outMomCE   = ce.LastMomentum;
        outTCE     = ce.LastT;
        outAPCE    = ce.LastAP;

        outPhiRe   = re.LastPhi;
        outThetaRe = re.LastTheta;
        outMomRe   = re.LastMomentum;
        outTRe     = re.LastT;
        outAPRe    = re.LastAP;

        outPhiReD   = reD.LastPhi;
        outThetaReD = reD.LastTheta;
        outMomReD   = reD.LastMomentum;
        outTReD     = reD.LastT;
        outAPReD    = reD.LastAP;
        
        outPDG = re.PDG;
        
        outThetaReDIFF = outThetaRe - outThetaCE;
        outPhiReDIFF = outPhiRe - outPhiCE;
        outTReDIFF = outTRe - outThetaCE;
        outAPReDIFF = outAPRe - outAPCE;
        outFOMReDIFF = 0;
        
        outThetaReDDIFF = outThetaReD - outThetaCE;
        outPhiReDDIFF = outPhiReD - outPhiCE;
        outTReDDIFF = outTReD - outThetaCE;
        outAPReDDIFF = outAPReD - outAPCE;
        outFOMReDDIFF = 0;
        

        outTree->Fill();
        outTree2->Fill();
        outTree3->Fill();
        outTree4->Fill();
        outTree5->Fill();
    }

    outFile->cd();
    outTree->Write();
    outTree2->Write();
    outTree3->Write();
    outTree4->Write();
    outTree5->Write();
    outFile->Close();
}
