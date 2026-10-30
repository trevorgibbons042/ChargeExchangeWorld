#include "FitCore.h"
#include "CutFiles.h"
#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TROOT.h>
#include <TMath.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/TBufferMerger.hxx>
#include <ROOT/TThreadExecutor.hxx>

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
  const double z_end = 23000;

  // Lab momentum of the beam neutron, MeV (kBeamMomentum, FitCore.h). LastMomentum is the CMS momentum of the model, so it cannot be
  // used as the lab momentum in t and in the AP exponent (the model uses the lab pin.mag() there).
  const double pBeam = kBeamMomentum;
}

using Key = std::pair<int, int>;
struct Angles {double LastPhi, LastTheta, LastMomentum, LastT, LastAP; int PDG; double FinalKE; ULong64_t entry = 0; int passResD = 1; double resDRms = 0;};

struct Track {std::vector<Angles> points;};

void CalcAPandT(std::string energy = "fKEnergy"){ // energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the Angles tree
    const double apScale = ReadSettings().at("APScale"); // settings.csv from the sim
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    ROOT::EnableImplicitMT(NThreads());
    
    bool Verbose = false;
    
    ROOT::RDataFrame dfCEV("ChargeExchange_values", "cuts/cut_CE.root");
    ROOT::RDataFrame dfAnglesAll(kAngleTree, "cuts/cut_angles_all.root");
    auto dfAngles = dfAnglesAll.Filter("passQuality==1 && isFromCE==1"); // the fitted CE protons that pass the quality cut, ReconAngleAll.cpp fits every proton
    
    using Key3 = std::tuple<int, int, int>;
    using TrackMap = std::map<Key3, Track>;
    // One map per thread slot, filled without a lock. Merged here on all the threads: the tracks of every slot go into the bucket of their shard (event % nShards, one slot per thread),
    // then each shard merges its own buckets. The last point of a track is used below, so each merged track is sorted back into entry order. The result is one map per shard.
    const unsigned nShards = std::max(1u, NThreads());
    ROOT::TThreadExecutor pool;
    auto mergeSlots = [&](std::vector<TrackMap>& slotMaps){
        std::vector<std::vector<std::vector<std::pair<Key3, Track>>>> buckets(slotMaps.size(), std::vector<std::vector<std::pair<Key3, Track>>>(nShards));
        pool.Foreach([&](unsigned slot){
            for (auto &[key, t] : slotMaps[slot]) buckets[slot][(unsigned)std::get<0>(key) % nShards].emplace_back(key, std::move(t));
            slotMaps[slot].clear();
        }, ROOT::TSeqU((unsigned)slotMaps.size()));
        std::vector<TrackMap> merged(nShards);
        pool.Foreach([&](unsigned shard){
            for (auto &slotBuckets : buckets) {
                for (auto &kt : slotBuckets[shard]) {
                    auto &dst = merged[shard][kt.first].points;
                    dst.insert(dst.end(), kt.second.points.begin(), kt.second.points.end());
                }
                slotBuckets[shard].clear(); slotBuckets[shard].shrink_to_fit();
            }
            for (auto &[key, t] : merged[shard]) {
                std::stable_sort(t.points.begin(), t.points.end(),
                    [](const Angles &a, const Angles &b){return a.entry < b.entry;});
            }
        }, ROOT::TSeqU(nShards));
        return merged;
    };


    std::vector<TrackMap> slotTracksCE(dfCEV.GetNSlots());
    dfCEV.ForeachSlot([&](unsigned slot,
    int eventID, int trackID, double LastPhi, double LastTheta, double LastMomentum, double LastT, double LastAP, double z, ULong64_t entry){
        auto &t = slotTracksCE[slot][{eventID, trackID, 0}];
        // the CE phi is the lab phi of the proton (LastPhiLab, same secondary as LastTheta), wrapped to 0 to 2 pi like the recon phi
        t.points.push_back({LastPhi < 0 ? LastPhi + 2.0 * M_PI : LastPhi, LastTheta, LastMomentum, LastT, LastAP, 2112, 0.0, entry});
    },{"fEvent", "TrackID", "LastPhiLab", "LastTheta", "LastMomentum", "LastT", "LastAP", "fZ", "rdfentry_"});
    std::vector<TrackMap> tracksCE = mergeSlots(slotTracksCE);


    std::vector<TrackMap> slotTracksAngles(dfAngles.GetNSlots());
    dfAngles.ForeachSlot([&](unsigned slot,
    int eventID, int ParentID, double outTheta, double outPhi, double outDTheta, double outDPhi, double finalKE, int pdg, ULong64_t entry, int passResD, double resDRms){
        if (pdg != 2212) {return;}

        auto &tRecon = slotTracksAngles[slot][{eventID, ParentID, 1}];


        tRecon.points.push_back({outPhi, outTheta, 0.0, 0.0, 0.0, pdg, finalKE, entry});

        auto &tReconBlock = slotTracksAngles[slot][{eventID, ParentID, 2}];

        tReconBlock.points.push_back({outDPhi, outDTheta, 0.0, 0.0, 0.0, pdg, finalKE, entry, passResD, resDRms});
        },{"fEvent", "ParentID","outTheta", "outPhi", "outDTheta", "outDPhi", KeVar(suffix), "PDG","rdfentry_", "passResD", "resDRms"});
    std::vector<TrackMap> tracksAngles = mergeSlots(slotTracksAngles);
    

    if (Verbose == true){
        for (auto &shardMap : tracksAngles) for (auto &[key, tRecon] : shardMap) {
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


    // The five trees are filled in the same loop, so their entries line up (the plots use them as friends). TTree::Fill is not thread-safe: the CE tracks are cut into chunks, every chunk is
    // filled by one thread into five trees of its own in a file the TBufferMerger gives it, and a chunk's five trees reach the output file together, so the entries stay aligned.
    // The rows come in the order the chunks finish; nothing reads them by position.
    ROOT::TBufferMerger merger("cuts/cut_APandT.root", energy == "fDEnergy" ? "UPDATE" : "RECREATE"); // detector-energy pass appends its _DE trees to the true-KE file
    std::vector<std::pair<const Key3*, const Track*>> ceTracks;
    for (auto &shardMap : tracksCE) for (auto &kt : shardMap) ceTracks.emplace_back(&kt.first, &kt.second);
    const std::size_t nChunks = std::min<std::size_t>(ceTracks.size(), nShards * 4);
    pool.Foreach([&](unsigned chunk){
        auto outFile = merger.GetFile();
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
        int outPassResD; // block residual cut (QualityCut.h), only on the block trees: the block plots cut on passResD == 1
        double outResDRms; // block fit residual of the track, on every tree: CutTrees.cpp cuts it into the block lists (kBlockCut, QualityCut.h)

        outTree->Branch("TrackID", &outTrackID);
        outTree->Branch("PDG", &outPDG);
        outTree->Branch("fEvent", &outEventID);
        outTree->Branch("fType", &outType);
        outTree->Branch("finalKE", &outFinalKE);
        outTree->Branch("resDRms", &outResDRms);
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
        outTree2->Branch("resDRms", &outResDRms);
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
        outTree3->Branch("resDRms", &outResDRms);
        outTree3->Branch("thetaReD", &outThetaReD);
        outTree3->Branch("phiReD", &outPhiReD);
        outTree3->Branch("tReD", &outTReD);
        outTree3->Branch("apReD", &outAPReD);
        outTree3->Branch("fomReD", &outFOMReD);
        outTree3->Branch("momReD", &outMomReD);
        outTree3->Branch("passResD", &outPassResD);
    
        outTree4->Branch("TrackID", &outTrackID);
        outTree4->Branch("PDG", &outPDG);
        outTree4->Branch("fEvent", &outEventID);
        outTree4->Branch("fType", &outType);
        outTree4->Branch("finalKE", &outFinalKE);
        outTree4->Branch("resDRms", &outResDRms);
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
        outTree5->Branch("resDRms", &outResDRms);
        outTree5->Branch("thetaReD", &outThetaReDDIFF);
        outTree5->Branch("phiReD", &outPhiReDDIFF);
        outTree5->Branch("tReD", &outTReDDIFF);
        outTree5->Branch("apReD", &outAPReDDIFF);
        outTree5->Branch("fomReD", &outFOMReDDIFF);
        outTree5->Branch("passResD", &outPassResD);
        for (std::size_t idx = chunk * ceTracks.size() / nChunks; idx < (chunk + 1) * ceTracks.size() / nChunks; idx++) {
            const Key3 &key = *ceTracks[idx].first;
            const Track &tCE = *ceTracks[idx].second;
                if (tCE.points.empty()) {continue;}

                int eventID = std::get<0>(key);
                int CETrackID = std::get<1>(key);
        
                auto itRecon = tracksAngles[(unsigned)eventID % nShards].find({eventID, CETrackID, 1});
                auto itReconD =tracksAngles[(unsigned)eventID % nShards].find({eventID, CETrackID, 2});
        
                if (itRecon == tracksAngles[(unsigned)eventID % nShards].end()) {continue;}
                if (itReconD == tracksAngles[(unsigned)eventID % nShards].end()) {continue;}
                auto &tRecon  = itRecon->second;
                auto &tReconD = itReconD->second;

                if (tRecon.points.empty() || tReconD.points.empty()) {continue;}

                const auto &ce = tCE.points.back();
                const auto &re = tRecon.points.back();
                const auto &reD = tReconD.points.back();
                outPassResD = reD.passResD;
                outResDRms = reD.resDRms;

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
        outFile->Write();
    }, ROOT::TSeqU((unsigned)nChunks));
}
