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
#include <array>
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

  // Lab momentum of the beam neutron, MeV (kBeamMomentum, FitCore.h). LastMomentum is the CMS momentum of the model, so it cannot be
  // used as the lab momentum in t and in the AP exponent (the model uses the lab pin.mag() there).
  const double pBeam = kBeamMomentum;
}

using Key = std::pair<int, int>;
struct Angles {double LastPhi, LastTheta, LastMomentum, LastT, LastAP; int PDG; double FinalKE; ULong64_t entry = 0; double resDRms = 0; double FinalKED = 0;}; // FinalKED: the detector KE, for the _DE trees

struct Track {std::vector<Angles> points;};

static void CalcAPandTBoth(bool both){ // both: the true-KE trees and the detector-KE (_DE) trees from ONE read of the inputs (only finalKE differs between them); else the true-KE trees only
    const double apScale = ReadSettings().at("APScale"); // settings.csv from the sim
    const std::vector<std::string> suffixes = both ? std::vector<std::string>{"", "_DE"} : std::vector<std::string>{""};
    ROOT::EnableImplicitMT(NThreads());
    
    bool Verbose = false;
    
    // entryNo: the true entry number of each row (rdfentry_ is wrong under MT on a tree with several clusters), taken before the Filter below
    ROOT::RDataFrame dfCEV0("ChargeExchange_values", "cuts/cut_CE.root");
    auto dfCEV = WithEntryNumber(dfCEV0, dfCEV0.GetNSlots());
    ROOT::RDataFrame dfAnglesAll0(kAngleTree, "cuts/cut_angles_all.root");
    auto dfAnglesAll = WithEntryNumber(dfAnglesAll0, dfAnglesAll0.GetNSlots());
    auto dfAngles = dfAnglesAll.Filter("isFromCE==1"); // the fitted CE protons, ReconAngleAll.cpp fits every proton
    
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
    },{"fEvent", "TrackID", "LastPhiLab", "LastTheta", "LastMomentum", "LastT", "LastAP", "fZ", "entryNo"});
    std::vector<TrackMap> tracksCE = mergeSlots(slotTracksCE);


    std::vector<TrackMap> slotTracksAngles(dfAngles.GetNSlots());
    dfAngles.ForeachSlot([&](unsigned slot,
    int eventID, int ParentID, double outTheta, double outPhi, double outDTheta, double outDPhi, double finalKE, int pdg, ULong64_t entry, double resDRms, double finalKED){
        if (pdg != 2212) {return;}

        auto &tRecon = slotTracksAngles[slot][{eventID, ParentID, 1}];


        tRecon.points.push_back({outPhi, outTheta, 0.0, 0.0, 0.0, pdg, finalKE, entry, 0.0, finalKED});

        auto &tReconBlock = slotTracksAngles[slot][{eventID, ParentID, 2}];

        tReconBlock.points.push_back({outDPhi, outDTheta, 0.0, 0.0, 0.0, pdg, finalKE, entry, resDRms, finalKED});
        },{"fEvent", "ParentID","outTheta", "outPhi", "outDTheta", "outDPhi", KeVar(""), "PDG","entryNo", "resDRms", KeVar("_DE")});
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
    ROOT::TBufferMerger merger("cuts/cut_APandT.root", "RECREATE");
    std::vector<std::pair<const Key3*, const Track*>> ceTracks;
    for (auto &shardMap : tracksCE) for (auto &kt : shardMap) ceTracks.emplace_back(&kt.first, &kt.second);
    const std::size_t nChunks = std::min<std::size_t>(ceTracks.size(), nShards * 4);
    pool.Foreach([&](unsigned chunk){
        auto outFile = merger.GetFile();
        int outTrackID, outPDG, outEventID, outType;
        double outFinalKE[2] = {0, 0}; // [0] true KE, [1] detector KE (the _DE trees); every other value is the same in both
        double outResDRms; // block fit residual of the track, on every tree (no cut on it)
        double outThetaCE, outPhiCE, outMomCE, outTCE, outAPCE;
        double outThetaRe, outPhiRe, outMomRe, outTRe, outAPRe;
        double outThetaReD, outPhiReD, outMomReD, outTReD, outAPReD;
        double outThetaReDIFF, outPhiReDIFF, outTReDIFF, outAPReDIFF;
        double outThetaReDDIFF, outPhiReDDIFF, outTReDDIFF, outAPReDDIFF;

        // the five trees of one energy (suffix "" true KE, "_DE" detector KE), finalKE bound to keVar
        auto makeSet = [&](const std::string& suffix, double* keVar) -> std::array<TTree*, 5> {
        // the six branches every tree starts with
        auto makeTree = [&](const std::string& name, const char* title){
            TTree* t = new TTree((name + suffix).c_str(), title);
            t->Branch("TrackID", &outTrackID);
            t->Branch("PDG", &outPDG);
            t->Branch("fEvent", &outEventID);
            t->Branch("fType", &outType);
            t->Branch("finalKE", keVar);
            t->Branch("resDRms", &outResDRms);
            return t;
        };
        TTree *outTree = makeTree("CE", "values from ChargeExchange directly");
        TTree *outTree2 = makeTree("CalcValsRe", "Calculated AP and t from recon");
        TTree *outTree3 = makeTree("CalcValsReD", "Calculated AP and t from recon Dectector Blocks");
        TTree *outTree4 = makeTree("CalcValsReDIFF", "Difference Calculated AP and t from recon");
        TTree *outTree5 = makeTree("CalcValsReDDIFF", "Difference Calculated AP and t from recon Dectector Blocks");

        outTree->Branch("thetaCE", &outThetaCE);
        outTree->Branch("phiCE", &outPhiCE);
        outTree->Branch("tCE", &outTCE);
        outTree->Branch("apCE", &outAPCE);
        outTree->Branch("momCE",  &outMomCE);

        outTree2->Branch("thetaRe", &outThetaRe);
        outTree2->Branch("phiRe", &outPhiRe);
        outTree2->Branch("tRe", &outTRe);
        outTree2->Branch("apRe", &outAPRe);
        outTree2->Branch("momRe",  &outMomRe);

        outTree3->Branch("thetaReD", &outThetaReD);
        outTree3->Branch("phiReD", &outPhiReD);
        outTree3->Branch("tReD", &outTReD);
        outTree3->Branch("apReD", &outAPReD);
        outTree3->Branch("momReD", &outMomReD);

        outTree4->Branch("thetaReD", &outThetaReDIFF);
        outTree4->Branch("phiReD", &outPhiReDIFF);
        outTree4->Branch("tReD", &outTReDIFF);
        outTree4->Branch("apReD", &outAPReDIFF);

        outTree5->Branch("thetaReD", &outThetaReDDIFF);
        outTree5->Branch("phiReD", &outPhiReDDIFF);
        outTree5->Branch("tReD", &outTReDDIFF);
        outTree5->Branch("apReD", &outAPReDDIFF);
        return {outTree, outTree2, outTree3, outTree4, outTree5};
        };
        std::vector<std::array<TTree*, 5>> trees; // one set per energy, filled in the same row loop so the entries of a set line up
        for (std::size_t s = 0; s < suffixes.size(); s++) trees.push_back(makeSet(suffixes[s], &outFinalKE[s]));
        for (std::size_t idx = chunk * ceTracks.size() / nChunks; idx < (chunk + 1) * ceTracks.size() / nChunks; idx++) {
            const Key3 &key = *ceTracks[idx].first;
            const Track &tCE = *ceTracks[idx].second;
                if (tCE.points.empty()) {continue;}

                int eventID = std::get<0>(key);
                int CETrackID = std::get<1>(key);
        
                const auto& shard = tracksAngles[(unsigned)eventID % nShards];
                auto itRecon = shard.find({eventID, CETrackID, 1});
                auto itReconD = shard.find({eventID, CETrackID, 2});
                if (itRecon == shard.end() || itReconD == shard.end()) {continue;}
                auto &tRecon  = itRecon->second;
                auto &tReconD = itReconD->second;

                if (tRecon.points.empty() || tReconD.points.empty()) {continue;}

                const auto &ce = tCE.points.back();
                const auto &re = tRecon.points.back();
                const auto &reD = tReconD.points.back();
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
                outFinalKE[0] = re.FinalKE; outFinalKE[1] = re.FinalKED;
        
                outThetaReDIFF = outThetaRe - outThetaCE;
                outPhiReDIFF = std::remainder(outPhiRe - outPhiCE, 360.0);
                outTReDIFF = outTRe - outTCE;
                outAPReDIFF = outAPRe - outAPCE;
        
                outThetaReDDIFF = outThetaReD - outThetaCE;
                outPhiReDDIFF = std::remainder(outPhiReD - outPhiCE, 360.0);
                outTReDDIFF = outTReD - outTCE;
                outAPReDDIFF = outAPReD - outAPCE;
        

                for (auto& set : trees) for (TTree* tr : set) tr->Fill();
        }
        outFile->Write();
    }, ROOT::TSeqU((unsigned)nChunks));
}

// Both energies in one pass: the true-KE trees and the detector-KE _DE trees are written to cut_APandT.root together (the inputs are read and the tracks built once).
// The beam spread (settings.csv BeamMomentumSD) as relative errors on t and A_P, the same for every event (FitCore.h BeamTFrac, BeamAFrac); the P errors of the plots include
// the A_P one. Printed and written to beam_errors.csv in the run folder.
static void WriteBeamErrors(){
    printf("beam p = %.1f MeV/c, SD = %.1f MeV/c: sigma_t/t = %.3f %%, sigma_A/A = %.3f %% (also the relative beam error on P)\n", kBeamMomentum, kBeamMomentumSD, 100 * kBeamTFrac, 100 * kBeamAFrac);
    if (FILE* f = fopen("beam_errors.csv", "w")){
        fprintf(f, "BeamMomentum,BeamMomentumSD,sigma_t_over_t,sigma_A_over_A\n%.6f,%.6f,%.8g,%.8g\n", kBeamMomentum, kBeamMomentumSD, kBeamTFrac, kBeamAFrac);
        fclose(f);
    }
}

// pass "all" (default): both energies. pass "true": the true-KE pass only (cut_APandT.root without the _DE trees), enough for the 6_cuts plots (Submit6Cuts.sh), which read CalcValsReDDIFF only.
void CalcAPandT(const char* pass = "all"){ // a const char*: a std::string macro argument fails on the Mac cling
    WriteBeamErrors();
    CalcAPandTBoth(std::string(pass) != "true");
}
