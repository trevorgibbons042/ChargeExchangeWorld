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


struct Pos {double x, y, z;};

struct Track {Pos p0, p2;
    bool has0 = false;
    bool has2 = false;
    ULong64_t first0 = std::numeric_limits<ULong64_t>::max();
    ULong64_t first2 = std::numeric_limits<ULong64_t>::max();
};

void ReconAngle(){
    ROOT::EnableImplicitMT();
    
    ROOT::RDataFrame df("ParticleInfo", "cut_EOL.root");
    TFile *outFile = new TFile("cut_angles.root", "RECREATE");
    TTree *outTree = new TTree("Angles", "Reconstructed Angles");

    int outTrackID;
    int outPDG;
    double outAngle;
    int outEventID;
    
    outTree->Branch("TrackID", &outTrackID);
    outTree->Branch("PDG", &outPDG);
    outTree->Branch("Angle", &outAngle);
    outTree->Branch("fEvent", &outEventID);
    
    std::map<std::tuple<int, int, int>, Track> tracks;
    std::mutex m;
    df.Foreach([&](
        int eventID, int trackID, int pdg, int boundary,
        double x, double y, double z, ULong64_t entry){
        std::lock_guard<std::mutex> lock(m);
        
        auto &t = tracks[{eventID, trackID, pdg}];
        if ((boundary == 1 || boundary == 0) && entry < t.first2) {
            t.p0 = {x, y, z};
            t.has0 = true;
            t.first2 = entry;
        } if (boundary == 2) {
            t.p2 = {x, y, z};
            t.has2 = true;
            }
        },{"fEvent","TrackID", "pdg", "boundaryType", "fX", "fY", "fZ", "rdfentry_"});

    for (auto &[key, t] : tracks) {
        if (!t.has0 || !t.has2){continue;}
        
        double dx = t.p2.x - t.p0.x;
        double dy = t.p2.y - t.p0.y;
        double dz = t.p2.z - t.p0.z;
        
        double length = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (length == 0){continue;}
        
        double theta = std::acos(dz / length);
        
        outEventID = std::get<0>(key);
        outTrackID = std::get<1>(key);
        outPDG = std::get<2>(key);
        outAngle = theta;
        outTree->Fill();
    }
    outFile->cd();
    outTree->Write();
    outFile->Close();
}
