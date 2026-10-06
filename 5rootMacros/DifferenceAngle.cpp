#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

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
#include <iostream>

struct AngleInfo { double theta; double phi;};

struct DifferenceResult {int eventID;int trackID;int pdg;
    double thetaEOL; double phiEOL; double thetaRecon; double phiRecon;
    double deltaTheta; double deltaPhi;
    double thetaReconD; double phiReconD; double deltaThetaD;
    double deltaPhiD;};

bool Verbose = false;

void DifferenceAngle(){
    ROOT::EnableImplicitMT();
    
    ROOT::RDataFrame dfEOL("CalcVals", "cut_EOL.root");
    ROOT::RDataFrame dfAngles("Angles", "cut_angles.root");
    
    using AngleMap = std::map<std::tuple<int, int, int>, AngleInfo>;
    const unsigned int nSlotsEOL = dfEOL.GetNSlots();
    std::vector<AngleMap> slotMaps(nSlotsEOL);
    
    dfEOL.ForeachSlot([&](
    unsigned int slot,
    int eventID, int trackID, int pdg,
    double theta, double phi) {
        
        if (Verbose == true){std::cout << "slot: " << slot << ", eventID: " << eventID << ", trackID: " << trackID << ", pdg: " << pdg << ", theta: " << theta << ", phi: " << phi << std::endl;}
        
        auto key = std::make_tuple(eventID, trackID, pdg);
        slotMaps[slot][key] = {theta, phi};
    },{"fEvent","TrackID","pdg","thetaCalcXY0detector0Pos","phiCalcXY0detector0Pos"});
    
    AngleMap EOLAngles;
    for (auto &slotMap : slotMaps) {
        for (auto &[key, angle] : slotMap) {
            EOLAngles[key] = angle;
            }
        }
    if (Verbose == true){std::cout << "Tracks in CalcVals: "<< EOLAngles.size() << std::endl;}
    
    const unsigned int nSlotsAngles = dfAngles.GetNSlots();
    std::vector<std::vector<DifferenceResult>>slotResults(nSlotsAngles);
    
    dfAngles.ForeachSlot([&](unsigned int slot,
    int eventID, int trackID, int pdg,
    double outTheta, double outPhi,
    double outDTheta,double outDPhi){
        if (Verbose == true){
            std::cout << "slot: " << slot << ", eventID: " << eventID << ", trackID: " << trackID << ", pdg: "<< pdg << ", outTheta: " << outTheta << ", outPhi: " << outPhi << "outDTheta: " << outDTheta << ", outDPhi: " << outDPhi << std::endl;
        }
        
        auto key = std::make_tuple(eventID,trackID, pdg);
        auto it = EOLAngles.find(key);
         
        // Track not found in CalcVals
        if (it == EOLAngles.end()) {return;}
        double thetaEOL = it->second.theta;
        double phiEOL = it->second.phi;
        
        double deltaTheta = outTheta - thetaEOL;
        double deltaPhi = outPhi - phiEOL;
        
        double deltaThetaD = outDTheta - thetaEOL;
        double deltaPhiD = outDPhi - phiEOL;
         
        slotResults[slot].push_back({
             eventID, trackID, pdg,
             thetaEOL, phiEOL,
             outTheta, outPhi,
             deltaTheta, deltaPhi,
             outDTheta, outDPhi,
             deltaThetaD, deltaPhiD});
    }, {"fEvent", "TrackID", "PDG", "outTheta", "outPhi", "outDTheta", "outDPhi"});
    
    
    if (Verbose == true){
        for(const auto &slot: slotResults){
            for(const auto &r: slot){
                std::cout <<
                "Event: "      << r.eventID <<
                " TrackID: "   << r.trackID <<
                " PDG: "       << r.pdg <<
                " thetaEOL: "  << r.thetaEOL <<
                " phiEOL: "    << r.phiEOL <<
                " outTheta: "  << r.thetaRecon <<
                " outPhi: "    << r.phiRecon <<
                " deltaTheta: "    << r.deltaTheta <<
                " deltaPhi: "      << r.deltaPhi <<
                
                " outDTheta: "  << r.thetaReconD <<
                " outDPhi: "    << r.phiReconD <<
                " deltaDTheta: "    << r.deltaThetaD <<
                " deltaDPhi: "      << r.deltaPhiD <<
                std::endl;
            }
        }
    std::size_t totalResults = 0;
    for (const auto &results : slotResults) {totalResults += results.size();}
    std::cout << "Tracks in Difference: " << totalResults<< std::endl;
    std::cout << "CalcVales entries = " << *dfEOL.Count() << "\n";
    std::cout << "CalcVales unique tracks = " << EOLAngles.size() << "\n";
    std::cout << "Angles entries = " << *dfAngles.Count() << "\n";
    }
    
    TFile *outFile = new TFile("angle_differences.root", "RECREATE");
    TTree *outTree =new TTree("AngleDifference","Angle Differences");
    
    int outEventID;
    int outTrackID;
    int outPDG;

    double thetaCalc;
    double phiCalc;

    double thetaRecon;
    double phiRecon;
    double deltaTheta;
    double deltaPhi;
    
    double thetaReconD;
    double phiReconD;
    double deltaThetaD;
    double deltaPhiD;
    
    outTree->Branch("fEvent",&outEventID);
    outTree->Branch("TrackID",&outTrackID);
    outTree->Branch("PDG",&outPDG);
    outTree->Branch("thetaCalcXY0detector0Pos",&thetaCalc);
    outTree->Branch("phiCalcXY0detector0Pos",&phiCalc);
    outTree->Branch("outTheta",&thetaRecon);
    outTree->Branch("outPhi",&phiRecon);
    
    outTree->Branch("deltaTheta",&deltaTheta);
    outTree->Branch("deltaPhi",&deltaPhi);
    
    outTree->Branch("outDTheta", &thetaReconD);
    outTree->Branch("outDPhi", &phiReconD);
    outTree->Branch("deltaThetaD", &deltaThetaD);
    outTree->Branch("deltaPhiD", &deltaPhiD);
    
    for (const auto &results : slotResults) {
        for (const auto &r : results) {

            outEventID = r.eventID;
            outTrackID = r.trackID;
            outPDG     = r.pdg;

            thetaCalc = r.thetaEOL;
            phiCalc   = r.phiEOL;

            thetaRecon = r.thetaRecon;
            phiRecon   = r.phiRecon;
            deltaTheta = r.deltaTheta;
            deltaPhi   = r.deltaPhi;
            
            thetaReconD = r.thetaReconD;
            phiReconD   = r.phiReconD;
            deltaThetaD = r.deltaThetaD;
            deltaPhiD   = r.deltaPhiD;

            outTree->Fill();
        }
    }
    outFile->cd();
    outTree->Write();
    outFile->Close();
    
}
