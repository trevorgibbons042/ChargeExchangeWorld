#include "Settings.h"
#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TROOT.h>
#include <TMath.h>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <memory>
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
  const double AP_A = 0.7079;   // same fit as G4ChargeExchange.cc (there per GeV, here momentum is converted to GeV)
  const double AP_k = 0.3860;
  const double z_end = 23000;

  // Lab momentum of the 25 GeV KE beam neutron, MeV. LastMomentum is the CMS momentum of the model, so it cannot be
  // used as the lab momentum in t and in the AP exponent (the model uses the lab pin.mag() there).
  const double beamKE = 25000.0;
  const double neutronMass = 939.565;
  const double pBeam = std::sqrt(beamKE*beamKE + 2.0*beamKE*neutronMass);

  // Transverse-only t in GeV^2: (p_n sin(theta))^2, theta = proton lab polar angle (rad), p_n = lab beam momentum
  // (MeV, pBeam). The full invariant also has an energy term, which the measured final KE (after energy loss) overwhelms.
  double TFromAngle(double pn, double theta){
    double pT = pn*std::sin(theta)/1000.0;
    return pT*pT;
  }

  // Same formula as the model, including its minus sign. pbeam in MeV, t in GeV^2.
  // AP times apScale (the sim's /CE/APScale), clamped to -1.0 at the low end like the model.
  double APFromT(double pbeam, double t, double apScale){
    return std::max(-1.0, apScale*-AP_A*(1-std::exp(-AP_k*pbeam/1000.0))*std::sqrt(std::max(t, 0.0)));
  }
}

using Key = std::pair<int, int>;
struct Angles {double LastPhi, LastTheta, LastMomentum, LastT, LastAP; int PDG; double FinalKE; ULong64_t entry = 0;};

struct Track {std::vector<Angles> points;};

void CalcAPandT(std::string energy = "fKEnergy"){ // energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the cut_angles file
    const double apScale = kAPScale; // Settings.h
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    ROOT::EnableImplicitMT(128);
    
    bool Verbose = false;
    
    ROOT::RDataFrame dfCEV("ChargeExchange_values", "cut_CE.root");
    ROOT::RDataFrame dfAngles("Angles" + suffix, "cut_angles.root");
    
    TFile *outFile = new TFile("cut_APandT.root", energy == "fDEnergy" ? "UPDATE" : "RECREATE"); // detector-energy pass appends its _DE trees to the true-KE file
    TTree *outTree = new TTree(("CE" + suffix).c_str(), "values from ChargeExchange directly");
    TTree *outTree2 = new TTree(("CalcValsRe" + suffix).c_str(), "Calculated AP and t from recon");
    TTree *outTree3 = new TTree(("CalcValsReD" + suffix).c_str(), "Calculated AP and t from recon Dectector Blocks");
    TTree *outTree4 = new TTree(("CalcValsReDIFF" + suffix).c_str(), "Difference Calculated AP and t from recon");
    TTree *outTree5 = new TTree(("CalcValsReDDIFF" + suffix).c_str(), "Difference Calculated AP and t from recon Dectector Blocks");
    
    int outTrackID;
    int outPDG;
    int outEventID;
    int outType;
    double outFinalKE;
    
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
    outTree->Branch("finalKE", &outFinalKE);
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
    outTree2->Branch("finalKE", &outFinalKE);
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
    outTree3->Branch("finalKE", &outFinalKE);
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
    outTree4->Branch("finalKE", &outFinalKE);
    outTree4->Branch("thetaReD", &outThetaReDIFF);
    outTree4->Branch("phiReD", &outPhiReDIFF);
    outTree4->Branch("tReD", &outTReDIFF);
    outTree4->Branch("apReD", &outAPReDIFF);
    outTree4->Branch("fomReD", &outFOMReDIFF);
    
    outTree5->Branch("TrackID", &outTrackID);
    outTree5->Branch("PDG", &outPDG);
    outTree5->Branch("fEvent", &outEventID);
    outTree5->Branch("fType", &outType);
    outTree5->Branch("finalKE", &outFinalKE);
    outTree5->Branch("thetaReD", &outThetaReDDIFF);
    outTree5->Branch("phiReD", &outPhiReDDIFF);
    outTree5->Branch("tReD", &outTReDDIFF);
    outTree5->Branch("apReD", &outAPReDDIFF);
    outTree5->Branch("fomReD", &outFOMReDDIFF);
    
    
    // One map per thread slot, merged after the loop, so the entry loops need no global lock. The last point of a
    // track is used below, so each merged track is sorted back into entry order.
    auto mergeSlots = [](std::vector<std::map<std::tuple<int, int, int>, Track>>& slotMaps){
        std::map<std::tuple<int, int, int>, Track> merged;
        for (auto &slotMap : slotMaps) {
            for (auto &[key, t] : slotMap) {
                auto &dst = merged[key].points;
                dst.insert(dst.end(), t.points.begin(), t.points.end());
            }
        }
        for (auto &[key, t] : merged) {
            std::stable_sort(t.points.begin(), t.points.end(),
                [](const Angles &a, const Angles &b){return a.entry < b.entry;});
        }
        return merged;
    };

    std::vector<std::map<std::tuple<int, int, int>, Track>> slotTracksCE(dfCEV.GetNSlots());
    dfCEV.ForeachSlot([&](unsigned slot,
    int eventID, int trackID, double LastPhi, double LastTheta, double LastMomentum, double LastT, double LastAP, double z, ULong64_t entry){
        auto &t = slotTracksCE[slot][{eventID, trackID, 0}];
        t.points.push_back({LastPhi, LastTheta, LastMomentum, LastT, LastAP, 2112, 0.0, entry});
    },{"fEvent", "TrackID", "LastPhi", "LastTheta", "LastMomentum", "LastT", "LastAP", "fZ", "rdfentry_"});
    std::map<std::tuple<int, int, int>, Track> tracksCE = mergeSlots(slotTracksCE);

    std::vector<std::map<std::tuple<int, int, int>, Track>> slotTracksAngles(dfAngles.GetNSlots());
    dfAngles.ForeachSlot([&](unsigned slot,
    int eventID, int ParentID, double outTheta, double outPhi, double outDTheta, double outDPhi, double finalKE, int pdg, ULong64_t entry){
        if (pdg != 2212) {return;}
         
        auto &tRecon = slotTracksAngles[slot][{eventID, ParentID, 1}];
        
                     
        tRecon.points.push_back({outPhi, outTheta, 0.0, 0.0, 0.0, pdg, finalKE, entry});
        
        auto &tReconBlock = slotTracksAngles[slot][{eventID, ParentID, 2}];
        
        tReconBlock.points.push_back({outDPhi, outDTheta, 0.0, 0.0, 0.0, pdg, finalKE, entry});
        },{"fEvent", "ParentID","outTheta", "outPhi", "outDTheta", "outDPhi", "finalKE", "PDG","rdfentry_"});
    std::map<std::tuple<int, int, int>, Track> tracksAngles = mergeSlots(slotTracksAngles);
    
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

        outPhiCE   = ce.LastPhi*180/TMath::Pi();
        outThetaCE = ce.LastTheta*180/TMath::Pi();
        outMomCE   = ce.LastMomentum;
        outTCE     = ce.LastT;
        outAPCE    = ce.LastAP;

        outPhiRe   = re.LastPhi*180/TMath::Pi();
        outThetaRe = re.LastTheta*180/TMath::Pi();
        // recon t and AP use the lab beam momentum
        outMomRe   = pBeam;
        outTRe     = TFromAngle(pBeam, re.LastTheta);
        outAPRe    = APFromT(pBeam, outTRe, apScale);

        outPhiReD   = reD.LastPhi*180/TMath::Pi();
        outThetaReD = reD.LastTheta*180/TMath::Pi();
        outMomReD   = pBeam;
        outTReD     = TFromAngle(pBeam, reD.LastTheta);
        outAPReD    = APFromT(pBeam, outTReD, apScale);
        
        outPDG = re.PDG;
        outFinalKE = re.FinalKE;
        
        outThetaReDIFF = outThetaRe - outThetaCE;
        outPhiReDIFF = std::remainder(outPhiRe - outPhiCE, 360.0);
        outTReDIFF = outTRe - outTCE;
        outAPReDIFF = outAPRe - outAPCE;
        outFOMReDIFF = 0;
        
        outThetaReDDIFF = outThetaReD - outThetaCE;
        outPhiReDDIFF = std::remainder(outPhiReD - outPhiCE, 360.0);
        outTReDDIFF = outTReD - outTCE;
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
