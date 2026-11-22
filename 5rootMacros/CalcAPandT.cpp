#include "FitCore.h"
#include "CutFiles.h"
#include <TTree.h>
#include <TROOT.h>
#include <TMath.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/TBufferMerger.hxx>
#include <ROOT/TThreadExecutor.hxx>

#include <algorithm>
#include <array>
#include <map>
#include <cmath>
#include <tuple>
#include <vector>


namespace
{

  // Lab momentum of the beam neutron, MeV (kBeamMomentum, FitCore.h), used in t and in the AP exponent (the model uses the lab pin.mag() there).
  const double pBeam = kBeamMomentum;
}

// One row of a track: the recon rows carry the est neutron range (ThMin, ThMax, PhHalf: rad, ReconAngleAll.cpp) and the lab angles (LabPhi, LabTheta: about the lab z axis, no neutron estimate), ZdcKEAllD the detector KE.
// ThB*, PhB*: the block grid range in the est neutron frame; ThL*, PhL*: the same about the lab z axis (block only).
struct Angles {double LastPhi, LastTheta; ULong64_t entry = 0; double ZdcKEAllD = 0; double LabPhi = 0, LabTheta = 0; double ThMin = 0, ThMax = 0, PhHalf = 0; double ThBMin = NAN, ThBMax = NAN, PhBHalf = NAN, ThLMin = NAN, ThLMax = NAN, PhLHalf = NAN;};

struct Track {std::vector<Angles> points;};

// Writes cut_APandT.root: the detector-KE trees (suffix _DE) CE_DE, CalcValsRe_DE, CalcValsReD_DE, CalcValsReDIFF_DE and CalcValsReDDIFF_DE (the inputs are read and the tracks built once).
static void CalcAPandTBoth(){
    const double apScale = ReadSettings().at("APScale"); // settings.csv from the sim
    ROOT::EnableImplicitMT(NThreads());
    
    // entryNo: the true entry number of each row (rdfentry_ is wrong under MT on a tree with several clusters), taken before the Filter below
    ROOT::RDataFrame dfCEV0("ChargeExchange_values", "cuts/cut_CE.root");
    auto dfCEV = WithEntryNumber(dfCEV0, dfCEV0.GetNSlots());
    TTree* anglesRaw = OpenCut("cuts/cut_angles_all.root", kAngleTree); // with the Events tree as a friend: zdcKEAllD and the other event columns by name
    if (!anglesRaw){Error("CalcAPandT", "no AnglesRaw in cuts/cut_angles_all.root (run ReconAngleAll.cpp)"); return;}
    ROOT::RDataFrame dfAnglesAll0(*anglesRaw);
    auto dfAnglesAll = WithEntryNumber(dfAnglesAll0, dfAnglesAll0.GetNSlots());
    auto dfAngles = dfAnglesAll.Filter("isFromCE==1"); // the recon tracks whose fit used a CE child hit (evClass 0 or 1), ReconAngleAll.cpp fits every track the finder finds
    
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
    int eventID, int trackID, double LastPhiLab, double LastThetaLab, ULong64_t entry, double dx, double dy, double dz){
        auto &t = slotTracksCE[slot][{eventID, trackID, 0}];
        // the CE angles are the lab angles of the proton (LastPhiLab, same secondary as LastTheta) relative to the TRUE neutron direction (CENDirX/Y/Z): CETHETA_NEUFRAME, CEPHI_NEUFRAME (rad, phi 0 to 2 pi)
        const double LastTheta = CEThetaNeu(LastThetaLab, LastPhiLab, dx, dy, dz), LastPhi = CEPhiNeu(LastThetaLab, LastPhiLab, dx, dy, dz);
        t.points.push_back({LastPhi, LastTheta, entry});
    },{"fEvent", "TrackID", "LastPhiLab", "LastTheta", "entryNo", "CENDirX", "CENDirY", "CENDirZ"});
    std::vector<TrackMap> tracksCE = mergeSlots(slotTracksCE);


    // An event can have several recon tracks that used a CE child hit (ReconAngleAll.cpp writes a row per track the finder finds). The CE of the event is paired with the one that used the most CE child hits, the lowest trackNo on a tie.
    std::map<int, std::pair<int, int>> pickOf; // event -> (CE child hits used, trackNo) of that track
    {
        std::vector<std::map<int, std::pair<int, int>>> slotPick(dfAngles.GetNSlots());
        dfAngles.ForeachSlot([&](unsigned slot, int eventID, int trackNo, int nCEUsed){
            auto it = slotPick[slot].find(eventID);
            if (it == slotPick[slot].end() || nCEUsed > it->second.first || (nCEUsed == it->second.first && trackNo < it->second.second)) slotPick[slot][eventID] = {nCEUsed, trackNo};
        }, {"fEvent", "trackNo", "nCEUsed"});
        for (const auto &m : slotPick) for (const auto &[ev, p] : m) {
            auto it = pickOf.find(ev);
            if (it == pickOf.end() || p.first > it->second.first || (p.first == it->second.first && p.second < it->second.second)) pickOf[ev] = p;
        }
    }

    std::vector<TrackMap> slotTracksAngles(dfAngles.GetNSlots());
    dfAngles.ForeachSlot([&](unsigned slot,
    int eventID, double thetaRelN, double phiRelN, double thetaRelND, double phiRelND, ULong64_t entry, double zdcKEAllD, double labPhi, double labTheta, double labPhiD, double labThetaD, double thMin, double thMax, double phHalf, double thMinD, double thMaxD, double phHalfD, double thBMinD, double thBMaxD, double phBHalfD, double thLMinD, double thLMaxD, double phLHalfD, int trackNo){
        if (trackNo != pickOf.at(eventID).second) return; // the recon track of the event the CE is paired with (pickOf)
        // the key has no CE track (-1): every CE track of the event gets the recon track picked above
        auto &tRecon = slotTracksAngles[slot][{eventID, -1, 1}];


        tRecon.points.push_back({phiRelN, thetaRelN, entry, zdcKEAllD, labPhi, labTheta, thMin, thMax, phHalf});

        auto &tReconBlock = slotTracksAngles[slot][{eventID, -1, 2}];

        tReconBlock.points.push_back({phiRelND, thetaRelND, entry, zdcKEAllD, labPhiD, labThetaD, thMinD, thMaxD, phHalfD, thBMinD, thBMaxD, phBHalfD, thLMinD, thLMaxD, phLHalfD});
        },{"fEvent", "thetaRelN", "phiRelN", "thetaRelND", "phiRelND", "entryNo", KeVar("_DE"), "outPhi", "outTheta", "outDPhi", "outDTheta", "thetaRelNMin", "thetaRelNMax", "phiRelNHalf", "thetaRelNDMin", "thetaRelNDMax", "phiRelNDHalf", "thetaRelNDBMin", "thetaRelNDBMax", "phiRelNDBHalf", "thetaLabDBMin", "thetaLabDBMax", "phiLabDBHalf", "trackNo"});
    std::vector<TrackMap> tracksAngles = mergeSlots(slotTracksAngles);
    

    // The five trees are filled in the same loop, so their entries line up (the plots use them as friends). TTree::Fill is not thread-safe: the CE tracks are cut into chunks, every chunk is
    // filled by one thread into five trees of its own in a file the TBufferMerger gives it, and a chunk's five trees reach the output file together, so the entries stay aligned.
    // The rows come in the order the chunks finish; nothing reads them by position.
    ROOT::TBufferMerger merger("cuts/cut_APandT.root", "RECREATE", 106); // zlib level 6 like the sim output (7run.cc)
    std::vector<std::pair<const Key3*, const Track*>> ceTracks;
    for (auto &shardMap : tracksCE) for (auto &kt : shardMap) ceTracks.emplace_back(&kt.first, &kt.second);
    const std::size_t nChunks = std::min<std::size_t>(ceTracks.size(), nShards * 4);
    pool.Foreach([&](unsigned chunk){
        auto outFile = merger.GetFile();
        int outEventID;
        double outZdcKEAllD; // the detector KE of the event (zdcKEAllD of AnglesRaw)
        double outThetaCE, outPhiCE, outTCE, outAPCE;
        double outThetaRe, outPhiRe, outTRe, outAPRe;
        double outThetaReD, outPhiReD, outTReD, outAPReD;
        double outThetaReLab, outPhiReLab, outThetaReDLab, outPhiReDLab; // recon angles about the lab z axis (degrees), no neutron estimate
        double outThetaReMin, outThetaReMax, outAPReMin, outAPReMax, outPhiReHalf; // est neutron range of the recon proton: theta (deg) and the AP at its ends, phi half width (rad)
        double outThetaReDMin, outThetaReDMax, outAPReDMin, outAPReDMax, outPhiReDHalf;
        double outThetaReDBMin, outThetaReDBMax, outAPReDBMin, outAPReDBMax, outPhiReDBHalf; // the block grid range alone (est neutron frame)
        double outThetaReDLabMin, outThetaReDLabMax, outPhiReDLabHalf; // the block grid range about the lab z axis
        double outThetaReDIFF, outPhiReDIFF, outTReDIFF, outAPReDIFF;
        double outThetaReDDIFF, outPhiReDDIFF, outTReDDIFF, outAPReDDIFF;

        // the five trees, zdcKEAllD bound to the detector KE
        auto makeSet = [&]() -> std::array<TTree*, 5> {
        const std::string suffix = "_DE";
        // the two branches every tree starts with
        auto makeTree = [&](const std::string& name, const char* title){
            TTree* t = new TTree((name + suffix).c_str(), title);
            t->Branch("fEvent", &outEventID);
            t->Branch("zdcKEAllD", &outZdcKEAllD);
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

        outTree2->Branch("thetaRe", &outThetaRe);
        outTree2->Branch("phiRe", &outPhiRe);
        outTree2->Branch("tRe", &outTRe);
        outTree2->Branch("apRe", &outAPRe);
        outTree2->Branch("thetaReLab", &outThetaReLab);
        outTree2->Branch("phiReLab", &outPhiReLab);
        outTree2->Branch("thetaReMin", &outThetaReMin);
        outTree2->Branch("thetaReMax", &outThetaReMax);
        outTree2->Branch("apReMin", &outAPReMin);
        outTree2->Branch("apReMax", &outAPReMax);
        outTree2->Branch("phiReHalf", &outPhiReHalf);

        outTree3->Branch("thetaReD", &outThetaReD);
        outTree3->Branch("phiReD", &outPhiReD);
        outTree3->Branch("tReD", &outTReD);
        outTree3->Branch("apReD", &outAPReD);
        outTree3->Branch("thetaReDLab", &outThetaReDLab);
        outTree3->Branch("phiReDLab", &outPhiReDLab);
        outTree3->Branch("thetaReDMin", &outThetaReDMin);
        outTree3->Branch("thetaReDMax", &outThetaReDMax);
        outTree3->Branch("apReDMin", &outAPReDMin);
        outTree3->Branch("apReDMax", &outAPReDMax);
        outTree3->Branch("phiReDHalf", &outPhiReDHalf);
        outTree3->Branch("thetaReDBMin", &outThetaReDBMin);
        outTree3->Branch("thetaReDBMax", &outThetaReDBMax);
        outTree3->Branch("apReDBMin", &outAPReDBMin);
        outTree3->Branch("apReDBMax", &outAPReDBMax);
        outTree3->Branch("phiReDBHalf", &outPhiReDBHalf);
        outTree3->Branch("thetaReDLabMin", &outThetaReDLabMin);
        outTree3->Branch("thetaReDLabMax", &outThetaReDLabMax);
        outTree3->Branch("phiReDLabHalf", &outPhiReDLabHalf);

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
        const std::array<TTree*, 5> trees = makeSet(); // filled in the same row loop so the entries line up
        for (std::size_t idx = chunk * ceTracks.size() / nChunks; idx < (chunk + 1) * ceTracks.size() / nChunks; idx++) {
            const Key3 &key = *ceTracks[idx].first;
            const Track &tCE = *ceTracks[idx].second;
                if (tCE.points.empty()) {continue;}

                int eventID = std::get<0>(key);
        
                const auto& shard = tracksAngles[(unsigned)eventID % nShards];
                auto itRecon = shard.find({eventID, -1, 1});
                auto itReconD = shard.find({eventID, -1, 2});
                if (itRecon == shard.end() || itReconD == shard.end()) {continue;}
                auto &tRecon  = itRecon->second;
                auto &tReconD = itReconD->second;

                if (tRecon.points.empty() || tReconD.points.empty()) {continue;}

                const auto &ce = tCE.points.back();
                const auto &re = tRecon.points.back();
                const auto &reD = tReconD.points.back();
                outEventID = eventID;

                outPhiCE   = ce.LastPhi*180/TMath::Pi();
                outThetaCE = ce.LastTheta*180/TMath::Pi();
                // CE t and AP from the CETHETA_NEUFRAME and the lab beam momentum, the same formula as the recon ones below (the sim LastT, LastAP are about the lab z axis)
                outTCE     = TFromAngle(pBeam, ce.LastTheta);
                outAPCE    = APFromT(pBeam, outTCE, apScale);

                outPhiRe   = re.LastPhi*180/TMath::Pi();
                outThetaRe = re.LastTheta*180/TMath::Pi();
                // recon t and AP use the lab beam momentum
                outTRe     = TFromAngle(pBeam, re.LastTheta);
                outAPRe    = APFromT(pBeam, outTRe, apScale);
                outPhiReLab = re.LabPhi*180/TMath::Pi(); outThetaReLab = re.LabTheta*180/TMath::Pi();
                outPhiReDLab = reD.LabPhi*180/TMath::Pi(); outThetaReDLab = reD.LabTheta*180/TMath::Pi();

                outPhiReD   = reD.LastPhi*180/TMath::Pi();
                outThetaReD = reD.LastTheta*180/TMath::Pi();
                outTReD     = TFromAngle(pBeam, reD.LastTheta);
                outAPReD    = APFromT(pBeam, outTReD, apScale);
        
                // the est neutron range: theta at both ends, and the AP there (the same t and AP formulas as the central value); phi half width as it is
                outThetaReMin = re.ThMin*180/TMath::Pi(); outThetaReMax = re.ThMax*180/TMath::Pi(); outPhiReHalf = re.PhHalf;
                outAPReMin = APFromT(pBeam, TFromAngle(pBeam, re.ThMin), apScale); outAPReMax = APFromT(pBeam, TFromAngle(pBeam, re.ThMax), apScale);
                outThetaReDMin = reD.ThMin*180/TMath::Pi(); outThetaReDMax = reD.ThMax*180/TMath::Pi(); outPhiReDHalf = reD.PhHalf;
                outAPReDMin = APFromT(pBeam, TFromAngle(pBeam, reD.ThMin), apScale); outAPReDMax = APFromT(pBeam, TFromAngle(pBeam, reD.ThMax), apScale);
                outThetaReDBMin = reD.ThBMin*180/TMath::Pi(); outThetaReDBMax = reD.ThBMax*180/TMath::Pi(); outPhiReDBHalf = reD.PhBHalf;
                outAPReDBMin = APFromT(pBeam, TFromAngle(pBeam, reD.ThBMin), apScale); outAPReDBMax = APFromT(pBeam, TFromAngle(pBeam, reD.ThBMax), apScale);
                outThetaReDLabMin = reD.ThLMin*180/TMath::Pi(); outThetaReDLabMax = reD.ThLMax*180/TMath::Pi(); outPhiReDLabHalf = reD.PhLHalf;

                outZdcKEAllD = re.ZdcKEAllD;
        
                outThetaReDIFF = outThetaRe - outThetaCE;
                outPhiReDIFF = std::remainder(outPhiRe - outPhiCE, 360.0);
                outTReDIFF = outTRe - outTCE;
                outAPReDIFF = outAPRe - outAPCE;
        
                outThetaReDDIFF = outThetaReD - outThetaCE;
                outPhiReDDIFF = std::remainder(outPhiReD - outPhiCE, 360.0);
                outTReDDIFF = outTReD - outTCE;
                outAPReDDIFF = outAPReD - outAPCE;
        

                for (TTree* tr : trees) tr->Fill();
        }
        outFile->Write();
    }, ROOT::TSeqU((unsigned)nChunks));
}

// The beam spread (settings.csv BeamMomentumSD) as relative errors on t and A_P, the same for every event (FitCore.h BeamTFrac, BeamAFrac); the P errors of the plots include
// the A_P one. Printed and written to beam_errors.csv in the run folder.
static void WriteBeamErrors(){
    printf("beam p = %.1f MeV/c, SD = %.1f MeV/c: sigma_t/t = %.3f %%, beam (A_max - A_min) / 2A = %.3f %% (also the relative beam error on P); model coefficients (A, k fit) (A_max - A_min) / 2A = %.3f %%\n", kBeamMomentum, kBeamMomentumSD, 100 * kBeamTFrac, 100 * kBeamAFrac, 100 * kModelAFrac);
    if (FILE* f = fopen("beam_errors.csv", "w")){
        fprintf(f, "BeamMomentum,BeamMomentumSD,sigma_t_over_t,sigma_A_over_A,sigma_A_over_A_model\n%.6f,%.6f,%.8g,%.8g,%.8g\n", kBeamMomentum, kBeamMomentumSD, kBeamTFrac, kBeamAFrac, kModelAFrac);
        fclose(f);
    }
}

// Writes beam_errors.csv, then the detector-KE trees of cuts/cut_APandT.root.
void CalcAPandT(){
    WriteBeamErrors();
    CalcAPandTBoth();
}
