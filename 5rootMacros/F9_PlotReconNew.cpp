#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLine.h>

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "FitCore.h"
#include "CutFiles.h"
#include <map>
#include <fstream>
#include <TH2D.h>
#include <TRandom3.h>

// 6_ReconNew: the event recon of ReconAngleAll.cpp (cuts/cut_angles_all.root): one fit per track the seed and extend finder finds in an event (hits taken off and the finder run again), no cut of any kind.
// Every plot overlays the three kinds of event (evClass): CE only (every fitted hit a CE child proton), CE mix (CE child hits and other particles), non-CE (no CE child hit).
// Output: six folders under plots/6_ReconNew_1/ (1_noCuts, 2_KECUT80, 3_nPlanesUsedEven, 4_Theta0_1, 5_avgTrackerHitsPerPlane1, 6_megaCut), the pdfs <n>_<name>.pdf titled with the cuts of the folder.
// Only with USE_BDT=1 in the environment, and with cuts/tmva_scores.root (TmvaApply.cpp, one row per AnglesRaw row in the same order) there are 36 BDT score cuts, -0.50 to 0.20 in steps of 0.02, each with the cut (>=) as its only cut:
// every cut writes cuts.txt and scan_stats.csv to plots/10_BDTcuts/scan_csv/BDT<cut>/ (F10_BDTCuts.cpp plots the csv files); only the cut whose total recon fit P (pAll) is closest to the true P (PPValue) is run again (mode 42 + its index) into plots/10_BDTcuts/1_BDT<cut>/ with the CE, CE + CE mix and total recon phi fits. No variable plots.
namespace {
struct NewRow {
    double tubeX0, tubeY0, tubeZ0, tubeX1, tubeY1, tubeZ1, thMin, thMax, phHalf, th2Min, th2Max, ph2Half, outTheta, outDTheta, outPhi, outDPhi, thetaRelND, phiRelND, impactD, r2D, resDRms, zdcKEAllD, bdt;
    int eventID, trackNo, nRecoEv, nTrackerTracks, nMaxPlane, nTrackerHits, nPlanesHit, nHits, nPlanesUsed, evClass, nCEUsed, nOtherUsed, nProtonUsed, nPionUsed, nKaonUsed, nNucleusUsed, nCEPart, nOtherPart, nCEReco, nOtherReco;
};
struct NewVar {
    const char* name; const char* title; int bins; double lo, hi; bool logY;
    std::function<double(const NewRow&)> value;
};
constexpr double kMegaKePercent = 90.0; // the KE cut of 6_megaCut: tighter than the 80% of 2_KECUT80, it costs about 4% of the CE and takes out much of the background
constexpr double kKeCutPercent = 80.0; // the KE cut: the smeared all-particle ZDC total (zdcKEAllD, charged and neutral) above this percent of the beam KE (80% of 25 GeV = 20000 MeV)
// the three curves are cumulative: CE only, CE only + CE mix, and every event (total recon particles)
const char* ClassLabel(int c){return c == 0 ? "CE recon" : c == 1 ? "CE + CE mix recon" : "total recon particles";}
int ClassColor(int c){return c == 0 ? CbCE() : c == 1 ? CbProton() : CbAll();}
double DegDiffPhi(double a, double b){return std::remainder(a - b, 2.0 * M_PI) * 180.0 / M_PI;}
}

void F9_PlotReconNew(){
    TFile* f = TFile::Open("cuts/cut_angles_all.root");
    TTree* t = f ? (TTree*)f->Get("AnglesRaw") : nullptr;
    if (!t){Error("F9_PlotReconNew", "no cuts/cut_angles_all.root (run ReconAngleAll.cpp)"); return;}
    if (!t->GetBranch("nNucleusUsed")){Error("F9_PlotReconNew", "cut_angles_all.root has no nNucleusUsed (rerun ReconAngleAll.cpp)"); return;}
    if (!t->GetBranch("zdcKEAllD")){Error("F9_PlotReconNew", "cut_angles_all.root has no zdcKEAllD, the smeared total ZDC KE of the KE cut (rerun ReconAngleAll.cpp)"); return;}

    const double rad = 180.0 / M_PI;
    const double keAxisMax = 1.2 * ReadSettings().at("BeamKE"); // the KE axes run to 120% of the beam KE (settings.csv), so the plots follow the sim energy
    const std::vector<NewVar> vars = {
        {"nPlanesUsed", "planes in the fit (exactly 1 charged hit each);planes", 9, 1.5, 10.5, true, [](const NewRow& x){return (double)x.nPlanesUsed;}},
        {"thetaBlock", "block fit theta in the est neutron frame;theta (degrees)", 55, 0.0, 2.2, true, [=](const NewRow& x){return x.thetaRelND * rad;}},
        {"thetaLabBlock", "block fit theta, lab frame;theta (degrees)", 60, 0.0, 6.0, true, [=](const NewRow& x){return x.outDTheta * rad;}},
        {"phiBlock", "block fit phi in the est neutron frame;phi (degrees)", 36, 0.0, 360.0, false, [=](const NewRow& x){return x.phiRelND * rad;}},
        {"phiLabBlock", "block fit phi, lab frame;phi (degrees)", 36, 0.0, 360.0, false, [=](const NewRow& x){return x.outDPhi * rad;}},
        {"blockDiff_theta", "block fit theta - true position fit theta;#Delta theta (degrees)", 80, -0.4, 0.4, true, [=](const NewRow& x){return (x.outDTheta - x.outTheta) * rad;}},
        {"blockDiff_phi", "block fit phi - true position fit phi;#Delta phi (degrees)", 80, -40.0, 40.0, true, [](const NewRow& x){return DegDiffPhi(x.outDPhi, x.outPhi);}},
        {"impactD", "impact of the block fit on the beam axis;impact (mm)", 60, 0.0, 60.0, true, [](const NewRow& x){return x.impactD;}},
        {"r2D", "R^{2} of the block fit;R^{2}", 60, 0.0, 1.0, true, [](const NewRow& x){return x.r2D;}},
        {"resDRms", "RMS distance of the block hits to the fit;resDRms (mm)", 60, 0.0, 12.0, true, [](const NewRow& x){return x.resDRms;}},
        {"KEEnd", "smeared total ZDC KE of the event (all particles, charged and neutral);KE (MeV)", 60, 0.0, keAxisMax, false, [](const NewRow& x){return x.zdcKEAllD;}},
        {"nMaxChargedPerPlaneEv", "most charged tracks in any one tracker plane of the event (all planes, not only the fit planes);tracks", 10, 0.5, 10.5, true, [](const NewRow& x){return (double)x.nMaxPlane;}},
        {"avgTrackerHitsPerPlane", "average charged tracker plane hits per hit plane in the event (nTrackerHitsEv / planes hit, all planes, not only the fit planes);nTrackerHitsEv / planes hit", 100, 0.5, 5.5, true, [](const NewRow& x){return x.nPlanesHit > 0 ? (double)x.nTrackerHits / x.nPlanesHit : NAN;}},
        {"nReconstructedTracksInEvent", "reconstructed tracks in the event (every fit the finder makes, whether or not they are the same particle), one entry per event (its first fit, so the legend N counts every fit);reconstructed tracks in the event", 10, 0.5, 10.5, true, [](const NewRow& x){return x.trackNo == 0 ? (double)x.nRecoEv : NAN;}},
        {"nUnreconstructedCEParticles", "CE particles (children of the carbon CE) that entered a tracker and have no hit in any fit of the event, one entry per event (its first fit, so the legend N counts every fit); events with no fit are not in AnglesRaw;CE particles not reconstructed", 8, -0.5, 7.5, true, [](const NewRow& x){return x.trackNo == 0 ? (double)std::max(0, x.nCEPart - x.nCEReco) : NAN;}},
        {"nUnreconstructedOtherParticles", "non-CE particles (every other particle) that entered a tracker and have no hit in any fit of the event, one entry per event (its first fit, so the legend N counts every fit); events with no fit are not in AnglesRaw;non-CE particles not reconstructed", 11, -0.5, 10.5, true, [](const NewRow& x){return x.trackNo == 0 ? (double)std::max(0, x.nOtherPart - x.nOtherReco) : NAN;}}};

    gStyle->SetOptStat(0);
    TCanvas c("cReconNew", "", 900, 650);
    c.SetRightMargin(0.04);
    const Long64_t n = t->GetEntries();
    // one multithreaded read (ColumnsUnordered, no ordered merge: the rows only need their columns aligned) of the 19 branches used here (not every branch of the tree), then every plot loop below runs over this vector in memory: 48 passes of GetEntry over all branches took hours at 100M events
    std::vector<NewRow> rows(n);
    ROOT::EnableImplicitMT(NThreads()); // implicit MT: ColumnsUnordered reads the tree on all the threads (NThreads: ANA_NTHREADS, FitCore.h, only returns the count)
    // the BDT score of every row (cuts/tmva_scores.root, TmvaApply.cpp, one row per AnglesRaw row in the same order). Not a friend tree: that made the read of AnglesRaw single-threaded.
    // The scores are read on all the threads with their entry numbers (bdtByEntry), AnglesRaw is read with its own entry number, and the two meet on it. Without that file, or with another row count, there are no BDT folders.
    bool haveBdt = false;
    std::vector<double> bdtByEntry;
    TFile* scoresFile = (gSystem->AccessPathName("cuts/tmva_scores.root") || !gSystem->Getenv("USE_BDT")) ? nullptr : TFile::Open("cuts/tmva_scores.root");
    if (scoresFile && !scoresFile->IsZombie()){
        TTree* st = (TTree*)scoresFile->Get("TMVAScores");
        haveBdt = st && st->GetEntries() == t->GetEntries();
        if (st && !haveBdt) Warning("F9_PlotReconNew", "cuts/tmva_scores.root has %lld rows, AnglesRaw %lld: no BDT folders (rerun TmvaApply.cpp)", (long long)st->GetEntries(), (long long)t->GetEntries());
        if (haveBdt){
            const auto sc = ColumnsUnordered(st, {"entryNo", "BDT"});
            bdtByEntry.assign(st->GetEntries(), NAN);
            for (size_t i = 0; i < sc[0].size(); i++) bdtByEntry[(size_t)sc[0][i]] = sc[1][i];
        }
    }
    int tmvaEvery = 1; // the events with fEvent % tmvaEvery == 0 trained the BDT (TmvaSplit.cpp, cuts/tmva_prescale.txt): the BDT cut folders leave them out, the fits and efficiencies use the other events only
    if (haveBdt){
        std::ifstream pf("cuts/tmva_prescale.txt");
        if (!(pf >> tmvaEvery) || tmvaEvery < 1){tmvaEvery = 1; Warning("F9_PlotReconNew", "no cuts/tmva_prescale.txt (rerun TmvaSplit.cpp): the BDT cut folders keep the events the BDT trained on");}
    }
    std::vector<std::string> colNames = {"outTheta", "outDTheta", "outPhi", "outDPhi", "thetaRelND", "phiRelND", "impactD", "r2D", "resDRms", "zdcKEAllD", "nHits", "nPlanesUsed", "evClass", "nCEUsed", "nOtherUsed", "nProtonUsed", "nPionUsed", "nKaonUsed", "nNucleusUsed", "nMaxChargedPerPlaneEv", "nTrackerHitsEv", "nPlanesHitEv", "thetaRelNDMin", "thetaRelNDMax", "phiRelNDHalf", "thetaRelNDBMin", "thetaRelNDBMax", "phiRelNDBHalf", "fEvent", "tubeX0", "tubeY0", "tubeZ0", "tubeX1", "tubeY1", "tubeZ1", "trackNo", "nTrackerTracksEv", "nCEPartEv", "nOtherPartEv", "nCERecoEv", "nOtherRecoEv"};
    if (haveBdt) colNames.push_back("entryNo"); // the entry number of the row in AnglesRaw, the key of bdtByEntry
    // the columns are read in groups of 9 and put into rows by entry number: all 41 at once (and the per-thread copies the merge makes) needed about 240 GB at 272M rows and the job died with exit 9 and no plots
    const size_t nCol = colNames.size();
    // the tree is also read in chunks of entries (the filter on entryNo): a 230 GB tree (700M rows) held its 10 columns of every row next to the rows and the job was killed with no message
    constexpr Long64_t kChunk = 40000000;
    for (Long64_t c0 = 0; c0 < n; c0 += kChunk)
    for (size_t g0 = 0; g0 < nCol; g0 += 9){
        const size_t g1 = std::min(nCol, g0 + 9);
        std::vector<std::string> names = {"entryNo"}; // the key that puts a row of this group into the same rows[] entry as in the other groups
        names.insert(names.end(), colNames.begin() + g0, colNames.begin() + g1);
        const auto cols = ColumnsUnordered(t, names, Form("entryNo >= %lld && entryNo < %lld", (long long)c0, (long long)std::min(n, c0 + kChunk)));
        MapRows((Long64_t)cols[0].size(), [&](Long64_t a, Long64_t z){
        for (Long64_t i = a; i < z; i++){
            NewRow& x = rows[(size_t)cols[0][i]];
            if (g0 == 0 && !haveBdt) x.bdt = NAN;
            for (size_t c = g0; c < g1; c++){
                const double v = cols[c - g0 + 1][i];
                switch (c){
                case 0: x.outTheta = v; break;
                case 1: x.outDTheta = v; break;
                case 2: x.outPhi = v; break;
                case 3: x.outDPhi = v; break;
                case 4: x.thetaRelND = v; break;
                case 5: x.phiRelND = v; break;
                case 6: x.impactD = v; break;
                case 7: x.r2D = v; break;
                case 8: x.resDRms = v; break;
                case 9: x.zdcKEAllD = v; break;
                case 10: x.nHits = (int)v; break;
                case 11: x.nPlanesUsed = (int)v; break;
                case 12: x.evClass = (int)v; break;
                case 13: x.nCEUsed = (int)v; break;
                case 14: x.nOtherUsed = (int)v; break;
                case 15: x.nProtonUsed = (int)v; break;
                case 16: x.nPionUsed = (int)v; break;
                case 17: x.nKaonUsed = (int)v; break;
                case 18: x.nNucleusUsed = (int)v; break;
                case 19: x.nMaxPlane = (int)v; break;
                case 20: x.nTrackerHits = (int)v; break;
                case 21: x.nPlanesHit = (int)v; break;
                case 22: x.thMin = v * rad; break;
                case 23: x.thMax = v * rad; break;
                case 24: x.phHalf = v; break;
                case 25: x.th2Min = v * rad; break;
                case 26: x.th2Max = v * rad; break;
                case 27: x.ph2Half = v; break;
                case 28: x.eventID = (int)v; break;
                case 29: x.tubeX0 = v; break;
                case 30: x.tubeY0 = v; break;
                case 31: x.tubeZ0 = v; break;
                case 32: x.tubeX1 = v; break;
                case 33: x.tubeY1 = v; break;
                case 34: x.tubeZ1 = v; break;
                case 35: x.trackNo = (int)v; break;
                case 36: x.nTrackerTracks = (int)v; break;
                case 37: x.nCEPart = (int)v; break;
                case 38: x.nOtherPart = (int)v; break;
                case 39: x.nCEReco = (int)v; break;
                case 40: x.nOtherReco = (int)v; break;
                case 41: x.bdt = bdtByEntry[(size_t)v]; break; // only there with haveBdt (entryNo is the last column)
                }
            }
        }
        return 0;
        });
    }
    // fits per event: the number of rows with the same fEvent (one row per fit), counted in an array indexed by the event number; the nReconstructedTracksInEvent plot (reconstructed particles per event, not only protons) reads it on the first fit of every event
    {
        int maxEv = 0;
        for (const auto& r : rows) if (r.eventID > maxEv) maxEv = r.eventID;
        std::vector<unsigned char> perEvent((size_t)maxEv + 1, 0);
        for (const auto& r : rows) if (r.eventID >= 0 && perEvent[r.eventID] < 255) perEvent[r.eventID]++;
        for (auto& r : rows) r.nRecoEv = r.eventID >= 0 ? perEvent[r.eventID] : 0;
    } // 35 columns of every row (36 GB at 128M rows) are not needed once rows is filled
    const auto settings = ReadSettings(); // settings.csv of the sim, read once
    if (settings.count("PPValue")) gPTrue = settings.at("PPValue"); // the true P, for the pull on the phi fit pages
    const double apScale = settings.at("APScale");
    const double keCut = kKeCutPercent / 100.0 * settings.at("BeamKE"); // whole events whose smeared total ZDC KE (zdcKEAllD, all particles) is not above kKeCutPercent of the beam KE are not drawn or counted
    const std::string kePct = Form("%.0f", kKeCutPercent);
    int made = 0;
    const char* scanEnv = gSystem->Getenv("SCAN_ONLY");
    const bool scanOnly = scanEnv && *scanEnv && std::string(scanEnv) != "0"; // SCAN_ONLY set: write only cuts.txt and scan_stats.csv of every folder, no pdf (SubmitPlotNew.sh)
    // plots/6_ReconNew_1/1_noCuts (every event; its KE plot alone shows the events above the KE cut) and plots/6_ReconNew_1/2_KECUT80 (every plot with the 80% KE cut).
    const int maxPlanes = (int)std::lround(ReadSettings().at("nBlocks") * (2.0 + ReadSettings().at("nAdditionalTrackers"))); // every tracker of the stack: 2 per block plus nAdditionalTrackers in each gap behind a block (12 for 1 block with 10 extra); 2 * nBlocks cut nearly every CE fit once there are additional trackers
    // Every folder after 2_ has one cut of its own, 6_megaCut has them all:
    // 3_nPlanesUsedEven: a whole-event cut on the planes in the fit, an even number and at most 2 * nBlocks (EvenHits, FitCore.h)
    // 4_Theta0_1: the block fit theta in the est neutron frame (the thetaBlock plot) from 0 to 1 degree
    // 5_avgTrackerHitsPerPlane1: on average one charged tracker hit per hit plane (nTrackerHitsEv == planes hit, the avgTrackerHitsPerPlane plot at 1)
    // 6_megaCut: the KE cut (90%), theta 0 to 1 degree and the average of 1 hit per plane and exactly 1 recon track in the event and at least 3 planes in the fit, all together
    // modes 6 to 41: the BDT score cut -0.50 to 0.20 in steps of 0.02 (plots/10_BDTcuts), only when the scores are there
    std::vector<int> modes = {0, 1, 2, 3, 4, 5, 100}; // mode 100: 6_nReconstructedTracksInEvent1, the events with exactly 1 reconstructed track
    if (haveBdt) for (int m = 6; m < 42; m++) modes.push_back(m);
    for (size_t mi = 0; mi < modes.size(); mi++){
    const int mode = modes[mi]; // modes 42 to 77: the BDT cut with the P closest to the true P, run again with its fit pdfs (pushed below after the last BDT cut)
    const bool mega = mode == 5;
    const bool oneTrack = mode == 100; // only events with exactly one reconstructed track (nRecoEv, the rows of the event)
    const bool bdtMode = mode >= 6 && !oneTrack;
    const int bdtHund = bdtMode ? -50 + 2 * (mode >= 42 ? mode - 42 : mode - 6) : 0; // the cut in hundredths, so a cut of 0 is exactly 0
    const double bdtCut = bdtHund / 100.0;
    const int bdtFullNo = (mode >= 42 && !oneTrack) ? 1 : 0; // the one cut with a plot folder, the one whose total recon P is closest to the true P
    const bool bdtFull = bdtFullNo > 0;
    const std::string bdtStat = Form("plots/10_BDTcuts/scan_csv/BDT%.2f", bdtCut); // cuts.txt and scan_stats.csv of every BDT cut
    const bool ke80 = mode == 1 || mega, evenPlanes = mode == 2, // not in the megaCut: it cost about 40% of the CE for almost no purity
         theta01 = mode == 3 || mega, avgOne = mode == 4 || mega, firstTrack = mega;
    const double keCutMode = mega ? kMegaKePercent / 100.0 * settings.at("BeamKE") : keCut; // this folder's KE cut (the megaCut has its own percent)
    const std::string kePctMode = mega ? Form("%.0f", kMegaKePercent) : kePct;
    const std::string dir = oneTrack ? "plots/6_ReconNew_1/6_nReconstructedTracksInEvent1" : mode == 0 ? "plots/6_ReconNew_1/1_noCuts" : mode == 1 ? "plots/6_ReconNew_1/2_KECUT" + kePct : mode == 2 ? "plots/6_ReconNew_1/3_nPlanesUsedEven" : mode == 3 ? "plots/6_ReconNew_1/4_Theta0_1" : mode == 4 ? "plots/6_ReconNew_1/5_avgTrackerHitsPerPlane1" : mode == 5 ? "plots/6_ReconNew_1/6_megaCut" : bdtFull ? std::string(Form("plots/10_BDTcuts/%d_BDT%.2f", bdtFullNo, bdtCut)) : bdtStat;
    const std::string statDir = bdtMode ? bdtStat : dir; // where cuts.txt and scan_stats.csv go
    const std::string infoDir = dir + "/1_additionalinfo"; // the main plot (KEEnd) in the folder, the other plots in 1_additionalinfo, the polarization fits of the megaCut in 2_phi_cosfit
    gSystem->mkdir((bdtMode ? dir : infoDir).c_str(), true); gSystem->mkdir(statDir.c_str(), true); // the BDT folders have no variable plots, so no 1_additionalinfo
    // the whole-event cuts of this folder
    const auto passes = [&](const NewRow& x){
        if (bdtMode && !(std::isfinite(x.bdt) && x.bdt >= bdtCut)) return false;
        if (bdtMode && tmvaEvery > 1 && x.eventID % tmvaEvery == 0) return false; // an event the BDT trained on
        if (ke80 && !(std::isfinite(x.zdcKEAllD) && x.zdcKEAllD > keCutMode)) return false;
        if (evenPlanes && !EvenHits(x.nPlanesUsed, maxPlanes)) return false;
        if (theta01 && !(x.thetaRelND * rad >= 0.0 && x.thetaRelND * rad <= 1.0)) return false; // NAN (no tube) fails too
        if (avgOne && !(x.nPlanesHit > 0 && x.nTrackerHits == x.nPlanesHit)) return false;
        if (firstTrack && x.nPlanesUsed < 3) return false; // no 2-plane tracks (with the even planes cut, 4 or more): two points always fit a line, nothing checks them
        if (firstTrack && x.nRecoEv != 1) return false; // only events with exactly 1 recon track (nRecoEv: the fits of the event)
        if (oneTrack && x.nRecoEv != 1) return false;
        return true;
    };
    std::string cutText; // every cut of the folder, named in the title; "no cuts" when there is none
    if (ke80) cutText += ", smeared total ZDC KE (all particles) above " + kePctMode + "% of the beam KE";
    if (evenPlanes) cutText += ", planes in the fit even and at most " + std::to_string(maxPlanes);
    if (theta01) cutText += ", theta 0 to 1 degree";
    if (avgOne) cutText += ", average tracker hits per hit plane = 1";
    if (firstTrack) cutText += ", exactly 1 reconstructed track in the event, at least 3 planes in the fit";
    if (oneTrack) cutText += ", exactly 1 reconstructed track in the event";
    if (bdtMode) cutText += Form(", BDT score >= %.2f (TmvaApply.cpp), no other cut", bdtCut);
    if (bdtMode && tmvaEvery > 1) cutText += Form(", events the BDT trained on (fEvent %% %d == 0) left out", tmvaEvery);
    if (cutText.empty()) cutText = ", no cuts";
    {std::ofstream cutsTxt(statDir + "/cuts.txt"); cutsTxt << "Folder " << statDir << ": every recon track of cuts/cut_angles_all.root (AnglesRaw, one row per track the finder finds), cuts" << cutText << "\n";}
    int number = 0;
    for (const auto& v : vars){
        if (scanOnly || bdtMode) break;
        number++;
        const bool keVar = std::string(v.name) == "KEEnd";
        TH1D* h[3];
        for (int k = 0; k < 3; k++){
            h[k] = new TH1D(Form("hRN_%s_%d", v.name, k), "", v.bins, keVar && ke80 ? keCutMode : v.lo, v.hi);
            h[k]->SetLineColor(ClassColor(k)); h[k]->SetLineWidth(2);
        }
        double nClass[3] = {0, 0, 0};
        // the rows run on all the threads (MapRows, FitCore.h): each range counts its own bins (no TH1 in the threads), the ranges are added up here; same bins as TH1::Fill (x below lo or from hi on is left out)
        const double hLo = h[0]->GetXaxis()->GetXmin(), hHi = h[0]->GetXaxis()->GetXmax();
        const int nb = v.bins;
        struct VarPart {std::vector<double> bins[3]; double nClass[3];};
        const auto varParts = MapRows(n, [&](Long64_t a, Long64_t z){
            VarPart p;
            for (int k = 0; k < 3; k++){p.bins[k].assign(nb, 0.0); p.nClass[k] = 0;}
            for (Long64_t i = a; i < z; i++){
                const NewRow& r = rows[i];
                if (r.evClass < 0 || r.evClass > 2) continue;
                const double x = v.value(r);
                if (!passes(r)) continue;
                for (int k = r.evClass; k < 3; k++){ // cumulative: CE only events are in all three, CE mix events in the last two, non-CE events in the last only
                    p.nClass[k]++;
                    if (std::isfinite(x) && x >= hLo && x < hHi) p.bins[k][std::min(nb - 1, (int)((x - hLo) / (hHi - hLo) * nb))]++;
                }
            }
            return p;
        });
        for (const auto& p : varParts) for (int k = 0; k < 3; k++){
            nClass[k] += p.nClass[k];
            for (int j = 0; j < nb; j++){h[k]->AddBinContent(j + 1, p.bins[k][j]); h[k]->SetEntries(h[k]->GetEntries() + p.bins[k][j]);}
        }
        double top = 0;
        for (int k = 0; k < 3; k++) top = std::max(top, h[k]->GetMaximum());
        const bool logY = v.logY || keVar; // the KE plot is log in every folder
        c.Clear(); c.SetLogy(logY);
        const std::string path = keVar ? Form("%s/1_KEEnd.pdf", dir.c_str()) : Form("%s/%d_%s.pdf", infoDir.c_str(), number, v.name);
        // PlotTitle numbers the title from the file name; "no cuts" in the title: the sample is every event with 2 or more planes of exactly 1 hit
        const std::string raw = v.title; const size_t semi = raw.find(';');
        h[0]->SetTitle(PlotTitle(path, raw.substr(0, semi) + cutText + raw.substr(semi)).c_str());
        h[0]->SetMaximum(logY ? top * ((std::string(v.name) == "nMaxChargedPerPlaneEv" || std::string(v.name) == "avgTrackerHitsPerPlane") ? 60.0 : 3.0) : top * 1.25);
        h[0]->SetMinimum(logY ? 0.5 : 0.0);
        h[0]->Draw("hist");
        h[2]->Draw("hist same"); h[1]->Draw("hist same"); h[0]->Draw("hist same"); // CE (orange) is drawn last, on top, then CE + CE mix (blue), then the total (green) at the back; h[0] is drawn first too, it owns the frame
        TLegend leg(keVar ? 0.42 : 0.62, 0.74, 0.96, 0.90);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        for (int k = 0; k < 3; k++) leg.AddEntry(h[k], Form("%s, N = %.0f", ClassLabel(k), nClass[k]), "l");
        leg.Draw();
        c.SaveAs(path.c_str());
        for (int k = 0; k < 3; k++) delete h[k];
        made++;
    }
    // the last plot: the fitted hits by particle type, one group of bars per kind of event. CE only events hold CE protons only; the CE mix events show which particles are mixed in
    // (a proton that is not a CE child is counted apart from the CE proton); the non-CE events are all "other" hits.
    if (!scanOnly && !bdtMode){
        number++;
        const int nCat = 6;
        const char* cats[nCat] = {"CE proton", "proton (not CE)", "pi+ / pi-", "K+ / K-", "nuclei (d, t, He)", "other"};
        const char* kinds[3] = {"CE only events", "CE mix events", "non-CE events"};
        TH1D* b[3];
        for (int k = 0; k < 3; k++){
            b[k] = new TH1D(Form("hRN_types_%d", k), "", nCat, 0.0, nCat);
            b[k]->SetFillColor(ClassColor(k)); b[k]->SetLineColor(ClassColor(k)); b[k]->SetBarWidth(0.28); b[k]->SetBarOffset(0.06 + 0.3 * k);
            for (int j = 0; j < nCat; j++) b[k]->GetXaxis()->SetBinLabel(j + 1, cats[j]);
        }
        double nEv[3] = {0, 0, 0};
        struct TypePart {double nEv[3] = {0, 0, 0}, cnt[3][6] = {{0}};};
        const auto typeParts = MapRows(n, [&](Long64_t a, Long64_t z){
            TypePart p;
            for (Long64_t i = a; i < z; i++){
                const NewRow& r = rows[i];
                if (r.evClass < 0 || r.evClass > 2 || !passes(r)) continue;
                p.nEv[r.evClass]++;
                double* q = p.cnt[r.evClass];
                q[0] += r.nCEUsed; q[1] += r.nProtonUsed; q[2] += r.nPionUsed; q[3] += r.nKaonUsed; q[4] += r.nNucleusUsed;
                q[5] += r.nOtherUsed - r.nProtonUsed - r.nPionUsed - r.nKaonUsed - r.nNucleusUsed;
            }
            return p;
        });
        for (const auto& p : typeParts) for (int k = 0; k < 3; k++){nEv[k] += p.nEv[k]; for (int j = 0; j < nCat; j++) b[k]->AddBinContent(j + 1, p.cnt[k][j]);}
        double top = 0;
        for (int k = 0; k < 3; k++) top = std::max(top, b[k]->GetMaximum());
        c.Clear(); c.SetLogy(true);
        const std::string path = Form("%s/%d_particleTypes.pdf", infoDir.c_str(), number);
        b[0]->SetTitle(PlotTitle(path, "fitted hits by particle type" + cutText + ";particle of the fitted hit;fitted hits").c_str());
        b[0]->SetMaximum(top * 60.0); b[0]->SetMinimum(0.5); // room above the bars for the legend
        b[0]->Draw("bar"); b[1]->Draw("bar same"); b[2]->Draw("bar same");
        TLegend leg(0.62, 0.74, 0.96, 0.90);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        for (int k = 0; k < 3; k++) leg.AddEntry(b[k], Form("%s, N = %.0f", kinds[k], nEv[k]), "f");
        leg.Draw();
        c.SaveAs(path.c_str());
        for (int k = 0; k < 3; k++) delete b[k];
        made++;
    }
    // The polarization fits of the three samples (PhiCosFit in FitCore.h, the unbinned fit): the block phi in the est neutron frame, A_eff from the recon block theta of each event,
    // the est neutron range and block grid range errors, the phi spike at 0 left out. Same cuts as the folder.
    // Every folder gets the CE and CE + CE mix fits, the megaCut also the total recon fit, its 50% subsamples and the matrix pages.
    {
        if (!bdtMode || bdtFull) gSystem->mkdir((dir + "/2_phi_cosfit").c_str(), true);
        // the whole-event cuts of this folder (passes) and the fit sample, applied to the rows in memory; the phi spike at 0 is left out and theta must exist (NAN fails)
        const auto inFit = [&](const NewRow& x, int sample){
            const double ph = x.phiRelND * rad;
            return (sample == 0 ? x.evClass == 0 : x.evClass <= sample) && passes(x) && ph > kBlockPhiLo && ph < kBlockPhiHi && x.thetaRelND == x.thetaRelND;
        };
        // A_eff at the recon block theta of a row and at the ends of its two est neutron theta ranges (the centre value where a range is missing); r1, r2 say a range exists
        struct AEffIn {double a, a1lo, a1hi, a2lo, a2hi; bool r1, r2;};
        const auto aeffOf = [&](const NewRow& x){
            const double a = AEffOfTheta(x.thetaRelND * rad, apScale);
            const bool r1 = std::isfinite(x.thMin) && std::isfinite(x.thMax), r2 = std::isfinite(x.th2Min) && std::isfinite(x.th2Max);
            return AEffIn{a, r1 ? AEffOfTheta(x.thMin, apScale) : a, r1 ? AEffOfTheta(x.thMax, apScale) : a, r2 ? AEffOfTheta(x.th2Min, apScale) : a, r2 ? AEffOfTheta(x.th2Max, apScale) : a, r1, r2};
        };
        c.SetLogy(false); // the canvas still has the log scale of the plot before
        gStyle->SetOptStat(1111); gStyle->SetOptFit(11); // fit parameters with errors (the fits are unbinned: no chi2 or probability) // PhiCosFit moves and fills the stats box with the fit, so it needs one
                const char* sampleName[3] = {"CE recon protons", "CE + CE mix recon protons", "total recon particles"};
        const char* sampleFile[3] = {"CE", "CEandMix", "all"};
        // The CE sample is fitted at the actual CE angle (ChargeExchange_values: the true phi and theta of the CE proton in the true neutron frame, CEPHI_NEUFRAME) and its own LastAP, not at the recon angle in the est neutron frame:
        // the events are the CE recon events that pass this folder's cuts (evClass 0 and cutExpr, matched by fEvent), only the angle of each is the true one. The rows go to a small file so Columns reads them on all the threads.
        TTree* ceT = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
        gROOT->cd();
        std::vector<std::vector<double>> ce;
        std::map<long long, size_t> byEvent;
        double smearBias = 0; // the signed A_eff bias of the recon theta: the CE recon events of this folder, mean AP at the recon block theta minus mean AP at the true CE theta; it is corrected in A_eff and P of the two recon fits and quoted apart from the uncertainty (the CE fit uses the true angle)
        if (ceT){
            ce = Columns(ceT, {"fEvent", "CEPhiNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "CEThetaNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "LastAP", "CENDirX", "CENDirY", "CENDirZ"});
            for (size_t i = 0; i < ce[0].size(); i++) byEvent.emplace((long long)ce[0][i], i);
            // on all the threads (MapRows): each row range sums its own rows, the sums are added in range order
            struct SmearPart {double n = 0, sRec = 0, sTrue = 0;};
            const auto smearParts = MapRows(n, [&](Long64_t a, Long64_t z){
                SmearPart p;
                for (Long64_t i = a; i < z; i++){
                    const NewRow& x = rows[i];
                    if (!inFit(x, 0)) continue;
                    const auto it = byEvent.find((long long)x.eventID);
                    if (it == byEvent.end()) continue;
                    p.n++; p.sRec += AEffOfTheta(x.thetaRelND * rad, apScale); p.sTrue += AEffOfTheta(ce[2][it->second], apScale);
                }
                return p;
            });
            double nA = 0, sRec = 0, sTrue = 0;
            for (const auto& p : smearParts){nA += p.n; sRec += p.sRec; sTrue += p.sTrue;}
            if (nA > 0) smearBias = (sRec - sTrue) / nA;
            // 1_additionalinfo/2_matrix2 and 3_matrix3 of the megaCut (the CE recon events of this folder, every cut on):
            // matrix2 M_recB_est_to_lab: the recon block angle in the est neutron frame (x, the "truth" of the page) against the recon block angle in the lab frame (y); matrix3 N_true_to_est: the true neutron direction (CENDir) against the est one (the tube centre at the block middle).
            // theta, phi and the flattened (theta slice, phi bin); the probability pages are normalised per column of x, the counts pages are the counts. Bins (0.1 deg theta bins, 36 phi bins, 10 x 50 flat bins).
            if (mega && !scanOnly){
                gStyle->SetPalette(kViridis);
                const double tMaxM = 2.3, tMaxN = 2.0;
                const int kSlice = 10, kFlat = 50;
                const auto wrap360 = [](double a){a = std::fmod(a, 360.0); return a < 0 ? a + 360.0 : a;};
                const auto flatBin = [&](double th, double ph, double tm){return std::min(kSlice - 1, (int)(th / (tm / kSlice))) * kFlat + std::min(kFlat - 1, (int)(ph / (360.0 / kFlat)));};
                TH2D* m2[3] = {new TH2D("hM2_theta", "", 23, 0, tMaxM, 23, 0, tMaxM), new TH2D("hM2_phi", "", 36, 0, 360, 36, 0, 360), new TH2D("hM2_flat", "", kSlice * kFlat, 0, kSlice * kFlat, kSlice * kFlat, 0, kSlice * kFlat)};
                TH2D* m3[3] = {new TH2D("hN_theta", "", 50, 0, tMaxN, 50, 0, tMaxN), new TH2D("hN_phi", "", 36, 0, 360, 36, 0, 360), new TH2D("hN_flat", "", kSlice * kFlat, 0, kSlice * kFlat, kSlice * kFlat, 0, kSlice * kFlat)};
                for (const NewRow& x : rows){
                    if (x.evClass != 0 || !passes(x) || !std::isfinite(x.thetaRelND) || !std::isfinite(x.outDTheta)) continue;
                    const double thE = x.thetaRelND * rad, phE = wrap360(x.phiRelND * rad), thL = x.outDTheta * rad, phL = wrap360(x.outDPhi * rad);
                    if (thE < tMaxM && thL < tMaxM){
                        m2[0]->Fill(thE, thL); m2[1]->Fill(phE, phL); m2[2]->Fill(flatBin(thE, phE, tMaxM) + 0.5, flatBin(thL, phL, tMaxM) + 0.5);
                    }
                    const auto it = byEvent.find((long long)x.eventID);
                    if (it == byEvent.end() || !std::isfinite(x.tubeZ0) || !std::isfinite(x.tubeZ1)) continue;
                    const double tx = ce[4][it->second], ty = ce[5][it->second], tz = ce[6][it->second];
                    const double ex = 0.5 * (x.tubeX0 + x.tubeX1), ey = 0.5 * (x.tubeY0 + x.tubeY1), ez = 0.5 * (x.tubeZ0 + x.tubeZ1);
                    const double thT = std::acos(std::max(-1.0, std::min(1.0, tz / std::sqrt(tx * tx + ty * ty + tz * tz)))) * rad, thN = std::acos(std::max(-1.0, std::min(1.0, ez / std::sqrt(ex * ex + ey * ey + ez * ez)))) * rad;
                    const double phT = wrap360(std::atan2(ty, tx) * rad), phN = wrap360(std::atan2(ey, ex) * rad);
                    if (thT < tMaxN && thN < tMaxN){
                        m3[0]->Fill(thT, thN); m3[1]->Fill(phT, phN); m3[2]->Fill(flatBin(thT, phT, tMaxN) + 0.5, flatBin(thN, phN, tMaxN) + 0.5);
                    }
                }
                const auto save2D = [&](TH2D* h, const std::string& path, const std::string& title, const std::string& axes, const char* zTitle, bool norm){
                    TH2D* pr = (TH2D*)h->Clone(Form("%s_%s", h->GetName(), norm ? "p" : "n"));
                    if (norm){ // P(est | truth): each column of x (the truth) sums to 1
                        for (int bx = 1; bx <= pr->GetNbinsX(); bx++){
                            double col = 0; for (int by = 1; by <= pr->GetNbinsY(); by++) col += h->GetBinContent(bx, by);
                            for (int by = 1; by <= pr->GetNbinsY(); by++) pr->SetBinContent(bx, by, col > 0 ? h->GetBinContent(bx, by) / col : 0.0);
                        }
                    }
                    double lo = 1e300; for (int bx = 1; bx <= pr->GetNbinsX(); bx++) for (int by = 1; by <= pr->GetNbinsY(); by++){const double v = pr->GetBinContent(bx, by); if (v > 0) lo = std::min(lo, v);}
                    c.Clear(); c.SetLogy(false); c.SetLogz(true); c.SetRightMargin(0.18);
                    pr->SetStats(false);
                    pr->SetMinimum(lo < 1e300 ? 0.5 * lo : 1e-3); if (norm) pr->SetMaximum(1.0);
                    pr->SetTitle(PlotTitle(path, title + cutText + axes).c_str());
                    pr->GetZaxis()->SetTitle(zTitle);
                    pr->Draw("COLZ");
                    c.SaveAs(path.c_str());
                    c.SetLogz(false); c.SetRightMargin(0.04);
                    delete pr;
                };
                const char* pg[3] = {"theta", "phi", "thetaphi"};
                const std::string ax2[3] = {";recon block theta, neutron estimate frame (deg);recon block theta, lab frame (deg)", ";recon block phi, neutron estimate frame (deg);recon block phi, lab frame (deg)", ";recon block (theta, phi) bin number, neutron estimate frame;recon block (theta, phi) bin number, lab frame"};
                const std::string ax3[3] = {";true neutron theta, lab (deg);estimated neutron theta, lab (deg)", ";true neutron phi, lab (deg);estimated neutron phi, lab (deg)", ";true neutron (theta, phi) bin number;estimated neutron (theta, phi) bin number"};
                gSystem->mkdir((infoDir + "/2_matrix2").c_str(), true); gSystem->mkdir((infoDir + "/3_matrix3").c_str(), true);
                for (int a = 0; a < 3; a++){
                    save2D(m2[a], Form("%s/2_matrix2/%d_M_recB_est_to_lab_%s.pdf", infoDir.c_str(), a + 1, pg[a]), "M_recB_est_to_lab, CE recon protons", ax2[a], "P(recon | true)", true);
                    save2D(m2[a], Form("%s/2_matrix2/%d_M_recB_est_to_lab_%s_counts.pdf", infoDir.c_str(), a + 4, pg[a]), "M_recB_est_to_lab counts, CE recon protons", ax2[a], "Protons per bin", false);
                    save2D(m3[a], Form("%s/3_matrix3/%d_N_true_to_est_%s.pdf", infoDir.c_str(), a + 1, pg[a]), "N_true_to_est, CE recon protons", ax3[a], "P(recon | true)", true);
                    save2D(m3[a], Form("%s/3_matrix3/%d_N_true_to_est_%s_counts.pdf", infoDir.c_str(), a + 4, pg[a]), "N_true_to_est counts, CE recon protons", ax3[a], "Neutrons per bin", false);
                    made += 4;
                }
                for (int a = 0; a < 3; a++){delete m2[a]; delete m3[a];}
            }
        }
        double nFit[3] = {0, 0, 0}, pFitV[3] = {0, 0, 0}, eFitV[3] = {0, 0, 0}, pRawV[3] = {0, 0, 0}; // pFitV, eFitV: P with the biases corrected and its uncertainty; pRawV: P of the fit before the correction
        { // every folder
            // scan_stats.csv of the folder, the numbers of the block scan plots (PlotBlockScan.cpp), one row: signal = events with a CE hit (evClass 0 and 1, CE + CE mix), background = non-CE events (evClass 2).
            // sAll, bAll count every recon event with no cut, nS, nB the ones that pass the cuts of this folder (passes); meanA is the mean A_eff at the recon block theta of the signal events in the cuts;
            // nCETrue is the number of CE events in the carbon blocks (cut_CE.root); fit 1 and fit 2 are the P fits of the CE + CE mix and the total recon samples, exactly as the plots of 2_phi_cosfit.
            // the counts and the three fit samples below are built on all the threads (MapRows): each row range returns its own partial result, the partials are merged in range order
            struct CountPart {double sAll = 0, bAll = 0, nS = 0, nB = 0, sumA = 0, nA = 0, s0All = 0, n0 = 0, sumA0 = 0, nA0 = 0;}; // s0All, n0, sumA0, nA0: the CE only tracks (evClass 0), the signal of F10_BDTCuts.cpp
            const auto countParts = MapRows(n, [&](Long64_t a, Long64_t z){
                CountPart p;
                for (Long64_t i = a; i < z; i++){
                    const NewRow& x = rows[i];
                    if (x.evClass < 0 || x.evClass > 2) continue;
                    if (bdtMode && tmvaEvery > 1 && x.eventID % tmvaEvery == 0) continue; // the efficiencies of the BDT folders are of the events the BDT did not train on
                    if (x.evClass <= 1) p.sAll++; else p.bAll++;
                    if (x.evClass == 0) p.s0All++;
                    if (!passes(x)) continue;
                    if (x.evClass <= 1){
                        p.nS++;
                        if (std::isfinite(x.thetaRelND)){p.sumA += AEffOfTheta(x.thetaRelND * rad, apScale); p.nA++;}
                    } else p.nB++;
                    if (x.evClass == 0){
                        p.n0++;
                        if (std::isfinite(x.thetaRelND)){p.sumA0 += AEffOfTheta(x.thetaRelND * rad, apScale); p.nA0++;}
                    }
                }
                return p;
            });
            double sAll = 0, bAll = 0, nS = 0, nB = 0, sumA = 0, nA = 0, s0All = 0, n0 = 0, sumA0 = 0, nA0 = 0;
            for (const auto& p : countParts){sAll += p.sAll; bAll += p.bAll; nS += p.nS; nB += p.nB; sumA += p.sumA; nA += p.nA; s0All += p.s0All; n0 += p.n0; sumA0 += p.sumA0; nA0 += p.nA0;}
            for (int s = 0; s < 3; s++){ // 0 the CE sample, at the true CE angle and its LastAP like 1_phiBlock_cosfit_CE.pdf; 1 CE + CE mix and 2 total recon, at the recon block angle
                struct FitPart {std::vector<double> ph; AMean a; PhiRange q1, q2;};
                auto fitParts = MapRows(n, [&](Long64_t lo, Long64_t hi){
                    FitPart p;
                    for (Long64_t i = lo; i < hi; i++){
                        const NewRow& x = rows[i];
                        if (!inFit(x, s)) continue;
                        if (s == 0){
                            const auto it = byEvent.find((long long)x.eventID);
                            if (it == byEvent.end()) continue;
                            const double a0 = ce[3][it->second];
                            p.ph.push_back(ce[1][it->second]); p.a.Add(a0, a0, a0, a0, a0);
                            continue;
                        }
                        const AEffIn e = aeffOf(x);
                        p.ph.push_back(x.phiRelND * rad);
                        p.a.Add(e.a, e.a1lo, e.a1hi, e.a2lo, e.a2hi);
                        if (e.r1) p.q1.Add(x.phHalf);
                        if (e.r2) p.q2.Add(x.ph2Half);
                    }
                    return p;
                });
                std::vector<double> phF; AMean aF; PhiRange q1, q2;
                size_t nPh = 0;
                for (const auto& p : fitParts) nPh += p.ph.size();
                phF.reserve(nPh);
                for (auto& p : fitParts){phF.insert(phF.end(), p.ph.begin(), p.ph.end()); std::vector<double>().swap(p.ph); aF += p.a; q1 += p.q1; q2 += p.q2;}
                nFit[s] = phF.size();
                if (phF.size() < 10) continue;
                aF.smearBias = s == 0 ? 0.0 : smearBias; // the CE fit uses the true angle: no recon bias term
                const PFit uf = UnbinnedP(phF, aF, kBlockPhiLo, kBlockPhiHi, nullptr, true, &q1, &q2);
                pFitV[s] = uf.p; eFitV[s] = uf.err; pRawV[s] = uf.pRaw;
            }
            std::ofstream scan(statDir + "/scan_stats.csv");
            scan.precision(10); // counts of millions and 250M neutrons stay whole numbers
            scan << "nBlocks,TotalNeutrons,PPValue,nCETrue,sAll,bAll,nS,nB,meanA,nFitCEMix,pCEMix,errCEMix,nFitAll,pAll,errAll,nFitCE,pCE,errCE,sCEOnly,nCEOnly,meanACEOnly,pRawAll,pRawCE\n"; // pRawAll, pRawCE are the P of the total recon and CE fits before the bias correction (pAll, pCE are corrected); // sCEOnly, nCEOnly, meanACEOnly (CE only tracks, evClass 0) are after the CE fit columns; the CE columns are last so older readers (PlotBlockScan.cpp reads the first 15, columns 16 to 18 only if present) still work
            scan << settings.at("nBlocks") << "," << settings.at("TotalNeutrons") << "," << settings.at("PPValue") << "," << byEvent.size() << "," << sAll << "," << bAll << "," << nS << "," << nB << "," << (nA > 0 ? sumA / nA : 0.0) << ","
                 << nFit[1] << "," << pFitV[1] << "," << eFitV[1] << "," << nFit[2] << "," << pFitV[2] << "," << eFitV[2] << "," << nFit[0] << "," << pFitV[0] << "," << eFitV[0] << "," << s0All << "," << n0 << "," << (nA0 > 0 ? sumA0 / nA0 : 0.0) << "," << pRawV[2] << "," << pRawV[0] << "\n";
        }
        if (bdtMode && mode == 41 && !scanOnly){ // every BDT cut is scanned: the cut with the total recon P (pAll) closest to the true P gets its fit pdfs
            int best = -1; double bestDiff = 1e300;
            for (int j = 0; j < 36; j++){
                std::ifstream in(Form("plots/10_BDTcuts/scan_csv/BDT%.2f/scan_stats.csv", (-50 + 2 * j) / 100.0));
                std::string head, line, tok; std::vector<double> v;
                if (!std::getline(in, head) || !std::getline(in, line)) continue;
                std::stringstream ss(line);
                while (std::getline(ss, tok, ',')) v.push_back(std::atof(tok.c_str()));
                if (v.size() < 15 || v[12] < 10) continue; // no total recon fit
                if (std::fabs(v[13] - v[2]) < bestDiff){bestDiff = std::fabs(v[13] - v[2]); best = j;}
            }
            if (best >= 0) modes.push_back(42 + best);
        }
        if (scanOnly && mega && !bdtMode){ // SCAN_ONLY has no pdf, but the 100 random 50% subsamples of the total recon fit (the page of 4_phiBlock_cosfit_all_50pct_subsamples.pdf) still go to subsamples.csv, same seeds and fit
            struct SubC {double ph, a, a1lo, a1hi, a2lo, a2hi, h1, h2; int ce;};
            std::vector<SubC> all;
            for (const NewRow& x : rows){
                if (!inFit(x, 2)) continue;
                const AEffIn e = aeffOf(x);
                all.push_back({x.phiRelND * rad, e.a, e.a1lo, e.a1hi, e.a2lo, e.a2hi, e.r1 ? x.phHalf : NAN, e.r2 ? x.ph2Half : NAN, x.evClass == 0 ? 1 : 0});
            }
            gSystem->mkdir((dir + "/2_phi_cosfit").c_str(), true);
            std::ofstream subCsv(dir + "/2_phi_cosfit/subsamples.csv");
            subCsv.precision(10);
            subCsv << "subsample,nAll,nCE,purity,P,err\n";
            if (all.size() >= 10) for (int s = 0; s < 100; s++){
                std::vector<size_t> idx(all.size());
                for (size_t i = 0; i < idx.size(); i++) idx[i] = i;
                TRandom3 rng(12345 + 200 + s);
                for (size_t i = idx.size() - 1; i > 0; i--) std::swap(idx[i], idx[rng.Integer(i + 1)]);
                const size_t m = size_t(std::lround(0.5 * idx.size()));
                std::vector<double> phS; AMean aM; PhiRange r1, r2; size_t ceS = 0;
                for (size_t q = 0; q < m; q++){
                    const SubC& e = all[idx[q]];
                    phS.push_back(e.ph); ceS += e.ce;
                    aM.Add(e.a, e.a1lo, e.a1hi, e.a2lo, e.a2hi);
                    if (std::isfinite(e.h1)) r1.Add(e.h1);
                    if (std::isfinite(e.h2)) r2.Add(e.h2);
                }
                aM.smearBias = smearBias;
                const PFit uf = UnbinnedP(phS, aM, kBlockPhiLo, kBlockPhiHi, nullptr, true, &r1, &r2);
                subCsv << s + 1 << "," << m << "," << ceS << "," << (m ? (double)ceS / m : 0.0) << "," << uf.p << "," << uf.err << "\n";
            }
        }
        if (scanOnly || (bdtMode && !bdtFull)) continue; // the plots below are not wanted
        for (int k = 0; k < (mega || bdtFull ? 3 : 2); k++){
            number++;
            const std::string phiHist = Form("hPhiNew_%d_%d", mode, k), phiPdf = Form("%s/2_phi_cosfit/%d_phiBlock_cosfit_%s.pdf", dir.c_str(), k + 1, sampleFile[k]);
            if (k == 0 && ceT){
                const std::string trueFile = Form("cuts/cut_ceTrueNew_%d.root", mode);
                {
                    TFile out(trueFile.c_str(), "RECREATE");
                    double ph = 0, th = 0, ap = 0;
                    TTree tt("CETrue", "CE recon events at the true CE angle");
                    tt.Branch("CEPhiNeuTrue", &ph); tt.Branch("CEThetaNeuTrue", &th); tt.Branch("LastAP", &ap);
                    for (const NewRow& x : rows){ // the CE recon events of this folder
                        if (!inFit(x, 0)) continue;
                        const auto it = byEvent.find((long long)x.eventID);
                        if (it == byEvent.end()) continue; // no carbon-block CE row for this event
                        ph = ce[1][it->second]; th = ce[2][it->second]; ap = ce[3][it->second];
                        tt.Fill();
                    }
                    tt.Write();
                }
                TFile* tf = TFile::Open(trueFile.c_str());
                TTree* tt = tf ? (TTree*)tf->Get("CETrue") : nullptr;
                if (tt){
                    PhiCosFit(&c, tt, "CEPhiNeuTrue", "CEThetaNeuTrue", "", phiHist,
                              std::string(sampleName[k]) + ", CEPHI_NEUFRAME Distribution, the actual CE angle, A_{eff} from the CE #theta and LastAP" + cutText, phiPdf,
                              {}, "LastAP", {}, {}, ClassColor(k), 0, 33);
                    made++;
                    tf->Close();
                    continue;
                }
            }
            // the rows of this fit go to a small file (the est neutron frame phi and theta in degrees, the two range errors), so PhiCosFit does not read the 23M rows of the tree again
            const std::string fitFile = Form("cuts/cut_phiFitNew_%d_%d.root", mode, k);
            {
                TFile out(fitFile.c_str(), "RECREATE");
                double ph = 0, th = 0, r1[3] = {0, 0, 0}, r2[3] = {0, 0, 0};
                TTree tp("PhiFit", "fitted events of one sample");
                tp.Branch("phiRelNDdeg", &ph); tp.Branch("thetaDeg", &th);
                tp.Branch("rThMin", &r1[0]); tp.Branch("rThMax", &r1[1]); tp.Branch("rPhHalf", &r1[2]);
                tp.Branch("r2ThMin", &r2[0]); tp.Branch("r2ThMax", &r2[1]); tp.Branch("r2PhHalf", &r2[2]);
                // the rows of the sample are found on all the threads (MapRows, in row order); the TTree is filled from them here, a TTree writer is not thread-safe
                const auto sampleRows = MapRows(n, [&](Long64_t a, Long64_t z){
                    std::vector<Long64_t> v;
                    for (Long64_t i = a; i < z; i++) if (inFit(rows[i], k)) v.push_back(i);
                    return v;
                });
                for (const auto& v : sampleRows) for (const Long64_t i : v){
                    const NewRow& x = rows[i];
                    ph = x.phiRelND * rad; th = x.thetaRelND * rad;
                    r1[0] = x.thMin; r1[1] = x.thMax; r1[2] = x.phHalf; r2[0] = x.th2Min; r2[1] = x.th2Max; r2[2] = x.ph2Half;
                    tp.Fill();
                }
                tp.Write();
            }
            TFile* ff = TFile::Open(fitFile.c_str());
            TTree* tp = ff ? (TTree*)ff->Get("PhiFit") : nullptr;
            if (tp){
                PhiCosFit(&c, tp, "phiRelNDdeg", "thetaDeg", "", phiHist,
                          std::string(sampleName[k]) + ", RECONPHIBL_ESTNEU Distribution, A_{eff} at the recon block #theta" + cutText, phiPdf,
                          {}, "", {"rThMin", "rThMax", "rPhHalf"}, {"r2ThMin", "r2ThMax", "r2PhHalf"}, ClassColor(k), smearBias, 33);
                ff->Close();
            }
            gSystem->Unlink(fitFile.c_str());
            made++;
            if (k == 2 && !bdtMode){ // the all recon fit again for 10 random 50% subsamples (no replacement, fixed seeds), one page each, then a page of the 10 P values, to see the effect of the statistics
                struct SubIn {double ph, a, a1lo, a1hi, a2lo, a2hi, h1, h2; int ce;}; // ce: 1 for a CE only recon track (evClass 0), for the purity of each subsample
                std::vector<SubIn> all;
                const double pTrue = settings.at("PPValue");
                for (const NewRow& x : rows){
                    if (!inFit(x, 2)) continue;
                    const AEffIn e = aeffOf(x);
                    all.push_back({x.phiRelND * rad, e.a, e.a1lo, e.a1hi, e.a2lo, e.a2hi, e.r1 ? x.phHalf : NAN, e.r2 ? x.ph2Half : NAN, x.evClass == 0 ? 1 : 0});
                }
                const int nSub = 100, nShow = 5, nb = 33; const double frac = 0.5; // 100 subsamples are fitted, the phi plots of the first 5 are drawn, the P page has all 100
                const std::string path = dir + "/2_phi_cosfit/4_phiBlock_cosfit_all_50pct_subsamples.pdf";
                std::vector<double> pSub, eSub, purSub; double pFull = 0, eFull = 0; size_t ceAll = 0; for (const SubIn& e : all) ceAll += e.ce;
                std::ofstream subCsv(dir + "/2_phi_cosfit/subsamples.csv"); // P, error and purity of each of the 100 subsamples (row full: the whole sample), for the plots over runs
                subCsv.precision(10);
                subCsv << "subsample,nAll,nCE,purity,P,err\n";
                const double purFull = all.empty() ? 0.0 : (double)ceAll / all.size();
                if (all.size() >= 10){
                    pFull = pFitV[2]; eFull = eFitV[2]; // the fit of the full sample, already done for scan_stats.csv (same rows, same A_eff and ranges)
                    for (int s = 0; s < nSub; s++){
                        std::vector<size_t> idx(all.size());
                        for (size_t i = 0; i < idx.size(); i++) idx[i] = i;
                        TRandom3 rng(12345 + 200 + s);
                        for (size_t i = idx.size() - 1; i > 0; i--) std::swap(idx[i], idx[rng.Integer(i + 1)]);
                        const size_t m = size_t(std::lround(frac * idx.size()));
                        std::vector<double> phS; AMean aM; PhiRange r1, r2; size_t ceS = 0;
                        TH1D hS(Form("hPhiSubNew_%d", s), "", nb, 0, 360);
                        for (size_t q = 0; q < m; q++){
                            const SubIn& e = all[idx[q]];
                            phS.push_back(e.ph); hS.Fill(e.ph); ceS += e.ce;
                            aM.Add(e.a, e.a1lo, e.a1hi, e.a2lo, e.a2hi);
                            if (std::isfinite(e.h1)) r1.Add(e.h1);
                            if (std::isfinite(e.h2)) r2.Add(e.h2);
                        }
                        aM.smearBias = smearBias;
                        const PFit uf = UnbinnedP(phS, aM, kBlockPhiLo, kBlockPhiHi, nullptr, true, &r1, &r2);
                        pSub.push_back(uf.p); eSub.push_back(uf.err); purSub.push_back(m ? (double)ceS / m : 0.0);
                        subCsv << s + 1 << "," << m << "," << ceS << "," << purSub.back() << "," << uf.p << "," << uf.err << "\n";
                        if (s >= nShow) continue;
                        c.Clear(); c.cd();
                        c.SetLeftMargin(0.16); c.SetRightMargin(0.25); c.SetTopMargin(0.17); c.SetBottomMargin(0.13);
                        hS.SetStats(0); hS.SetLineColor(ClassColor(2)); hS.SetLineWidth(2); hS.SetMinimum(0); hS.SetMaximum(1.1 * hS.GetMaximum());
                        hS.SetTitle(Form("random 50%% subsample %d of %d (first %d shown), %s, RECONPHIBL_ESTNEU Distribution, %d phi bins of %.1f deg%s;RECONPHIBL_ESTNEU (degrees);Count", s + 1, nSub, nShow, sampleName[2], nb, 360.0 / nb, cutText.c_str()));
                        hS.Draw("hist E"); // error bars sqrt(N) on every phi bin
                        TF1 fd("fSubDrawNew", "[0] * (1 + [1] * [2] * [3] * cos(x * TMath::Pi() / 180))", 0, 360); // corrected P and A_eff, D the phi dilution: the curve stays on the data
                        fd.SetParameters(hS.GetEntries() / nb, uf.p, aM.Corrected(), uf.dil); fd.SetLineColor(kBlack); fd.SetLineStyle(2); fd.SetLineWidth(3); fd.Draw("same");
                        TLine meanBin(0, hS.GetEntries() / nb, 360, hS.GetEntries() / nb); meanBin.SetLineStyle(3); meanBin.SetLineWidth(2); meanBin.Draw(); // dotted line at the mean bin content
                        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026);
                        tx.DrawLatex(0.18, 0.93, "Unbinned fit  N(#phi) #propto 1 + P #upoint A_{eff} cos#phi, dashed = fit, A_{eff} fixed at the recon block #theta");
                        tx.DrawLatex(0.18, 0.895, Form("%.0f entries (%.0f%% of %zu), %s = %.4f #pm %.4f,  P = %.3g #pm %.2g", hS.GetEntries(), 100.0 * hS.GetEntries() / all.size(), all.size(), AEffEq(), aM.Corrected(), aM.Err(), uf.p, uf.err));
                        tx.DrawLatex(0.18, 0.86, BlockPhiNote());
                        DrawErrorBySource(aM, uf, true);
                        c.Modified(); c.Update();
                        c.Print((path + (s == 0 ? "(" : "")).c_str());
                    }
                    c.Clear(); c.cd();
                    c.SetLeftMargin(0.12); c.SetRightMargin(0.04); c.SetTopMargin(0.17); c.SetBottomMargin(0.12);
                    double lo = pSub[0], hi = pSub[0], mean = 0, rms = 0;
                    for (size_t i = 0; i < pSub.size(); i++){lo = std::min(lo, pSub[i]); hi = std::max(hi, pSub[i]); mean += pSub[i] / pSub.size();}
                    for (const double p : pSub) rms += (p - mean) * (p - mean);
                    rms = pSub.size() > 1 ? std::sqrt(rms / (pSub.size() - 1)) : 0;
                    lo = std::min(lo, std::min(pFull, pTrue)); hi = std::max(hi, std::max(pFull, pTrue));
                    const double pad = 0.1 * (hi - lo > 0 ? hi - lo : 0.1);
                    TH1D hP("hSubPNew", Form("P of %d random 50%% subsamples, %s%s;Polarization fit P;Subsamples", nSub, sampleName[2], cutText.c_str()), 50, lo - pad, hi + pad); // 50 bins
                    for (const double p : pSub) hP.Fill(p);
                    hP.SetStats(0); hP.SetLineColor(ClassColor(2)); hP.SetLineWidth(2); hP.SetFillColorAlpha(ClassColor(2), 0.35); hP.SetMinimum(0); hP.SetMaximum(1.3 * hP.GetMaximum());
                    hP.Draw("hist");
                    TLine lf(pFull, 0, pFull, hP.GetMaximum()); lf.SetLineColor(ClassColor(2)); lf.SetLineWidth(3); lf.Draw();
                    TLine lt(pTrue, 0, pTrue, hP.GetMaximum()); lt.SetLineStyle(3); lt.SetLineWidth(3); lt.Draw();
                    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026);
                    tx.DrawLatex(0.14, 0.93, "solid line = P of the full sample, dotted black = true P (PPValue)");
                    tx.DrawLatex(0.14, 0.89, Form("full P = %.3g #pm %.2g, subsample mean %.3g, RMS %.2g, expected RMS %.2g (full error #times #sqrt{(1-f)/f}, f = 0.5)", pFull, eFull, mean, rms, eFull * std::sqrt((1 - frac) / frac)));
                    c.Modified(); c.Update();
                    c.Print(path.c_str());
                    // last page: the purity (CE only recon tracks / all recon tracks) of each of the subsamples
                    c.Clear(); c.cd();
                    double plo = purSub[0], phi = purSub[0], pmean = 0, prms = 0;
                    for (const double v : purSub){plo = std::min(plo, v); phi = std::max(phi, v); pmean += v / purSub.size();}
                    for (const double v : purSub) prms += (v - pmean) * (v - pmean);
                    prms = purSub.size() > 1 ? std::sqrt(prms / (purSub.size() - 1)) : 0;
                    plo = std::min(plo, purFull); phi = std::max(phi, purFull);
                    const double ppad = 0.1 * (phi - plo > 0 ? phi - plo : 0.01);
                    TH1D hQ("hSubPurNew", Form("purity of %d random 50%% subsamples, %s%s;purity (CE recon / all recon tracks);Subsamples", nSub, sampleName[2], cutText.c_str()), 50, plo - ppad, phi + ppad);
                    for (const double v : purSub) hQ.Fill(v);
                    hQ.SetStats(0); hQ.SetLineColor(ClassColor(2)); hQ.SetLineWidth(2); hQ.SetFillColorAlpha(ClassColor(2), 0.35); hQ.SetMinimum(0); hQ.SetMaximum(1.3 * hQ.GetMaximum());
                    hQ.Draw("hist");
                    TLine lq(purFull, 0, purFull, hQ.GetMaximum()); lq.SetLineColor(ClassColor(2)); lq.SetLineWidth(3); lq.Draw();
                    TLatex tq; tq.SetNDC(); tq.SetTextSize(0.026);
                    tq.DrawLatex(0.14, 0.93, "solid line = purity of the full sample");
                    tq.DrawLatex(0.14, 0.89, Form("full purity = %.3f, subsample mean %.3f, RMS %.3f", purFull, pmean, prms));
                    c.Modified(); c.Update();
                    c.Print((path + ")").c_str());
                    c.Clear(); c.SetLeftMargin(0.13); c.SetRightMargin(0.30); c.SetTopMargin(0.1); c.SetBottomMargin(0.13);
                    made += nShow + 2;
                }
            }
        }
        gStyle->SetOptStat(0);
    }
    if (bdtMode) continue; // the BDT folders stop after the fits
    // the last plot of the folder: the fraction of the recon events that are CE (CE only, and CE + CE mix) in each bin of the smeared total ZDC KE (the KEEnd plot), same cuts as the folder.
    // Per bin, fraction = events of that kind / all recon events, with the binomial error. The bins hold the events counted on all the threads (MapRows).
    {
        number++;
        const double kLo = ke80 ? keCutMode : 0.0, kHi = keAxisMax;
        const int nb = 60;
        struct FracPart {std::vector<double> ce, mix, all;};
        const auto fracParts = MapRows(n, [&](Long64_t a, Long64_t z){
            FracPart p; p.ce.assign(nb, 0.0); p.mix.assign(nb, 0.0); p.all.assign(nb, 0.0);
            for (Long64_t i = a; i < z; i++){
                const NewRow& r = rows[i];
                if (r.evClass < 0 || r.evClass > 2 || !passes(r)) continue;
                const double x = r.zdcKEAllD;
                if (!std::isfinite(x) || x < kLo || x >= kHi) continue;
                const int bin = std::min(nb - 1, (int)((x - kLo) / (kHi - kLo) * nb));
                p.all[bin]++;
                if (r.evClass == 0) p.ce[bin]++;
                if (r.evClass <= 1) p.mix[bin]++;
            }
            return p;
        });
        TH1D* fr[2];
        for (int k = 0; k < 2; k++){fr[k] = new TH1D(Form("hRN_keFrac_%d_%d", mode, k), "", nb, kLo, kHi); fr[k]->SetLineColor(ClassColor(k)); fr[k]->SetLineWidth(2);}
        for (int j = 0; j < nb; j++){
            double ce = 0, mix = 0, all = 0;
            for (const auto& p : fracParts){ce += p.ce[j]; mix += p.mix[j]; all += p.all[j];}
            if (all <= 0) continue;
            const double f0 = ce / all, f1 = mix / all;
            fr[0]->SetBinContent(j + 1, f0); fr[0]->SetBinError(j + 1, std::sqrt(f0 * (1.0 - f0) / all));
            fr[1]->SetBinContent(j + 1, f1); fr[1]->SetBinError(j + 1, std::sqrt(f1 * (1.0 - f1) / all));
        }
        c.Clear(); c.SetLogy(false);
        const std::string path = Form("%s/%d_KEEndFraction.pdf", infoDir.c_str(), number);
        fr[1]->SetTitle(PlotTitle(path, "fraction of the recon events that are CE, per bin of the smeared total ZDC KE" + cutText + ";KE (MeV);fraction of the recon events in the bin").c_str());
        fr[1]->SetMinimum(0.0); fr[1]->SetMaximum(1.05);
        fr[1]->Draw("hist"); fr[0]->Draw("hist same"); // CE only on top, CE + CE mix behind
        TLegend leg(0.12, 0.78, 0.60, 0.90);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(fr[0], "CE recon / total recon particles", "l");
        leg.AddEntry(fr[1], "CE + CE mix recon / total recon particles", "l");
        leg.Draw();
        c.SaveAs(path.c_str());
        for (int k = 0; k < 2; k++) delete fr[k];
        made++;
    }
    }
    printf(">>> F9_PlotReconNew: %d plots in plots/6_ReconNew_1 (1_noCuts, 2_KECUT%s, 3_ to 5_ one cut each, 6_megaCut)%s\n", made, kePct.c_str(), haveBdt ? ", BDT cut folders in plots/10_BDTcuts" : ", no cuts/tmva_scores.root so no BDT folders");
}
