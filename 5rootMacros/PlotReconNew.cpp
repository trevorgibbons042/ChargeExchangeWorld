#include <TFile.h>
#include "Palette.h"
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
#include "PhiCosFit.h"
#include "CutOpen.h"
#include <map>
#include <fstream>
#include <TH2D.h>
#include <TRandom3.h>
#include <TGraphErrors.h>

// 9_ReconNew: the event recon of ReconAngleAll.cpp(true) (cuts/cut_angles_new.root): one fit per event through the tracker planes with exactly 1 charged hit, no cut of any kind.
// Every plot overlays the three kinds of event (evClass): CE only (every fitted hit a CE child proton), CE mix (CE child hits and other particles), non-CE (no CE child hit).
// Output: plots/9_ReconNew/1_noCuts/<n>_<name>.pdf, titled "no cuts".
namespace {
struct NewRow {
    double tubeX0, tubeY0, tubeZ0, tubeX1, tubeY1, tubeZ1, thMin, thMax, phHalf, th2Min, th2Max, ph2Half, outTheta, outDTheta, outPhi, outDPhi, thetaRelND, phiRelND, impactD, r2D, resDRms, finalKED;
    int eventID, nMaxPlane, nTrackerHits, nPlanesHit, nHits, nPlanesUsed, evClass, nCEUsed, nOtherUsed, nProtonUsed, nPionUsed, nKaonUsed, nNucleusUsed;
};
struct NewVar {
    const char* name; const char* title; int bins; double lo, hi; bool logY;
    std::function<double(const NewRow&)> value;
};
constexpr double kKeCutPercent = 80.0; // the KE cut: summed (smeared) KE on the last tracker above this percent of the beam KE (80% of 25 GeV = 20000 MeV)
// the three curves are cumulative: CE only, CE only + CE mix, and every event (total recon particles)
const char* ClassLabel(int c){return c == 0 ? "CE recon" : c == 1 ? "CE + CE mix recon" : "total recon particles";}
int ClassColor(int c){return c == 0 ? CbCE() : c == 1 ? CbProton() : CbAll();}
double DegDiffPhi(double a, double b){return std::remainder(a - b, 2.0 * M_PI) * 180.0 / M_PI;}
}

void PlotReconNew(){
    TFile* f = TFile::Open("cuts/cut_angles_new.root");
    TTree* t = f ? (TTree*)f->Get("AnglesRaw") : nullptr;
    if (!t){Error("PlotReconNew", "no cuts/cut_angles_new.root (run ReconAngleAll.cpp(true))"); return;}
    NewRow r;
    t->SetBranchAddress("outTheta", &r.outTheta); t->SetBranchAddress("outDTheta", &r.outDTheta);
    t->SetBranchAddress("outPhi", &r.outPhi); t->SetBranchAddress("outDPhi", &r.outDPhi);
    t->SetBranchAddress("thetaRelND", &r.thetaRelND); t->SetBranchAddress("phiRelND", &r.phiRelND); t->SetBranchAddress("impactD", &r.impactD);
    t->SetBranchAddress("r2D", &r.r2D); t->SetBranchAddress("resDRms", &r.resDRms); t->SetBranchAddress("finalKED", &r.finalKED);
    t->SetBranchAddress("nHits", &r.nHits); t->SetBranchAddress("nPlanesUsed", &r.nPlanesUsed); t->SetBranchAddress("evClass", &r.evClass);
    if (!t->GetBranch("nNucleusUsed")){Error("PlotReconNew", "cut_angles_new.root has no nProtonUsed (rerun ReconAngleAll.cpp(true))"); return;}
    t->SetBranchAddress("nCEUsed", &r.nCEUsed); t->SetBranchAddress("nOtherUsed", &r.nOtherUsed);
    t->SetBranchAddress("nProtonUsed", &r.nProtonUsed); t->SetBranchAddress("nPionUsed", &r.nPionUsed);
    t->SetBranchAddress("nKaonUsed", &r.nKaonUsed); t->SetBranchAddress("nNucleusUsed", &r.nNucleusUsed);

    const double rad = 180.0 / M_PI;
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
        {"KEEnd", "summed KE of the charged hits on the last tracker, ZDC energy smearing;KE (MeV)", 60, 0.0, 30000.0, false, [](const NewRow& x){return x.finalKED;}},
        {"nMaxChargedPerPlaneEv", "most charged tracks in any one tracker plane of the event (all planes, not only the fit planes);tracks", 10, 0.5, 10.5, true, [](const NewRow& x){return (double)x.nMaxPlane;}},
        {"avgTrackerHitsPerPlane", "average charged tracker plane hits per hit plane in the event (nTrackerHitsEv / planes hit, all planes, not only the fit planes);nTrackerHitsEv / planes hit", 100, 0.5, 5.5, true, [](const NewRow& x){return x.nPlanesHit > 0 ? (double)x.nTrackerHits / x.nPlanesHit : NAN;}}};

    gStyle->SetOptStat(0);
    TCanvas c("cReconNew", "", 900, 650);
    c.SetRightMargin(0.04);
    const Long64_t n = t->GetEntries();
    // one multithreaded read (ColumnsUnordered, no ordered merge: the rows only need their columns aligned) of the 19 branches used here (not every branch of the tree), then every plot loop below runs over this vector in memory: 48 passes of GetEntry over all branches took hours at 100M events
    std::vector<NewRow> rows(n);
    NThreads(); // ROOT::EnableImplicitMT: Columns reads the tree on all the threads
    const auto cols = ColumnsUnordered(t, {"outTheta", "outDTheta", "outPhi", "outDPhi", "thetaRelND", "phiRelND", "impactD", "r2D", "resDRms", "finalKED", "nHits", "nPlanesUsed", "evClass", "nCEUsed", "nOtherUsed", "nProtonUsed", "nPionUsed", "nKaonUsed", "nNucleusUsed", "nMaxChargedPerPlaneEv", "nTrackerHitsEv", "nPlanesHitEv", "thetaRelNDMin", "thetaRelNDMax", "phiRelNDHalf", "thetaRelNDBMin", "thetaRelNDBMax", "phiRelNDBHalf", "fEvent", "tubeX0", "tubeY0", "tubeZ0", "tubeX1", "tubeY1", "tubeZ1"});
    for (Long64_t i = 0; i < n; i++){
        NewRow& x = rows[i]; size_t k = 0;
        x.outTheta = cols[k++][i]; x.outDTheta = cols[k++][i]; x.outPhi = cols[k++][i]; x.outDPhi = cols[k++][i]; x.thetaRelND = cols[k++][i]; x.phiRelND = cols[k++][i]; x.impactD = cols[k++][i]; x.r2D = cols[k++][i]; x.resDRms = cols[k++][i]; x.finalKED = cols[k++][i];
        x.nHits = (int)cols[k++][i]; x.nPlanesUsed = (int)cols[k++][i]; x.evClass = (int)cols[k++][i]; x.nCEUsed = (int)cols[k++][i]; x.nOtherUsed = (int)cols[k++][i]; x.nProtonUsed = (int)cols[k++][i]; x.nPionUsed = (int)cols[k++][i]; x.nKaonUsed = (int)cols[k++][i]; x.nNucleusUsed = (int)cols[k++][i]; x.nMaxPlane = (int)cols[k++][i]; x.nTrackerHits = (int)cols[k++][i]; x.nPlanesHit = (int)cols[k++][i];
        x.thMin = cols[k++][i] * rad; x.thMax = cols[k++][i] * rad; x.phHalf = cols[k++][i]; x.th2Min = cols[k++][i] * rad; x.th2Max = cols[k++][i] * rad; x.ph2Half = cols[k++][i]; x.eventID = (int)cols[k++][i];
        x.tubeX0 = cols[k++][i]; x.tubeY0 = cols[k++][i]; x.tubeZ0 = cols[k++][i]; x.tubeX1 = cols[k++][i]; x.tubeY1 = cols[k++][i]; x.tubeZ1 = cols[k++][i];
    }
    const double keCut = kKeCutPercent / 100.0 * ReadSettings().at("BeamKE"); // whole events whose summed (smeared) KE on the last tracker is not above kKeCutPercent of the beam KE are not drawn or counted
    const std::string kePct = Form("%.0f", kKeCutPercent);
    int made = 0;
    // plots/9_ReconNew_1/1_noCuts (every event; its KE plot alone shows the events above the KE cut) and plots/9_ReconNew_1/2_KECUT80 (every plot with the 80% KE cut).
    const int maxPlanes = MaxTrackerHits(); // 2 * nBlocks, the most planes a track can cross
    // Every folder after 2_ has one cut of its own, 6_megaCut has them all:
    // 3_nPlanesUsedEven: a whole-event cut on the planes in the fit, an even number and at most 2 * nBlocks (EvenHits, FitCore.h)
    // 4_Theta0_1: the block fit theta in the est neutron frame (the thetaBlock plot) from 0 to 1 degree
    // 5_avgTrackerHitsPerPlane1: on average one charged tracker hit per hit plane (nTrackerHitsEv == planes hit, the avgTrackerHitsPerPlane plot at 1)
    // 6_megaCut: the KE cut, the even planes, theta 0 to 1 degree and the average of 1 hit per plane, all together
    for (const int mode : {0, 1, 2, 3, 4, 5}){
    const bool mega = mode == 5;
    const bool ke80 = mode == 1 || mega, evenPlanes = mode == 2 || mega, theta01 = mode == 3 || mega, avgOne = mode == 4 || mega;
    const std::string dir = mode == 0 ? "plots/9_ReconNew_1/1_noCuts" : mode == 1 ? "plots/9_ReconNew_1/2_KECUT" + kePct : mode == 2 ? "plots/9_ReconNew_1/3_nPlanesUsedEven" : mode == 3 ? "plots/9_ReconNew_1/4_Theta0_1" : mode == 4 ? "plots/9_ReconNew_1/5_avgTrackerHitsPerPlane1" : "plots/9_ReconNew_1/6_megaCut";
    const std::string infoDir = dir + "/1_additionalinfo"; // like plots/6_cuts/<cut>: the main plot (KEEnd) in the folder, the other plots in 1_additionalinfo, the polarization fits of the megaCut in 2_phi_cosfit
    gSystem->mkdir(infoDir.c_str(), true);
    // the whole-event cuts of this folder
    const auto passes = [&](const NewRow& x){
        if (ke80 && !(std::isfinite(x.finalKED) && x.finalKED > keCut)) return false;
        if (evenPlanes && !EvenHits(x.nPlanesUsed, maxPlanes)) return false;
        if (theta01 && !(x.thetaRelND * rad >= 0.0 && x.thetaRelND * rad <= 1.0)) return false; // NAN (no tube) fails too
        if (avgOne && !(x.nPlanesHit > 0 && x.nTrackerHits == x.nPlanesHit)) return false;
        return true;
    };
    std::string cutText; // every cut of the folder, named in the title; "no cuts" when there is none
    if (ke80) cutText += ", KE cut above " + kePct + "% of the beam KE";
    if (evenPlanes) cutText += ", planes in the fit even and at most " + std::to_string(maxPlanes);
    if (theta01) cutText += ", theta 0 to 1 degree";
    if (avgOne) cutText += ", average tracker hits per hit plane = 1";
    if (cutText.empty()) cutText = ", no cuts";
    {std::ofstream cutsTxt(dir + "/cuts.txt"); cutsTxt << "Folder " << dir << ": every recon event of cuts/cut_angles_new.root (AnglesRaw, planes with exactly 1 charged hit), cuts" << cutText << "\n";}
    int number = 0;
    for (const auto& v : vars){
        number++;
        const bool keVar = std::string(v.name) == "KEEnd";
        const bool cutOn = ke80;
        TH1D* h[3];
        for (int k = 0; k < 3; k++){
            h[k] = new TH1D(Form("hRN_%s_%d", v.name, k), "", v.bins, keVar && ke80 ? keCut : v.lo, v.hi);
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
    {
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
    // The polarization fits of the three samples (PhiCosFit.h, the unbinned fit of 3_ReconCalcVals and 6_cuts): the block phi in the est neutron frame, A_eff from the recon block theta of each event,
    // the est neutron range and block grid range errors, the phi spike at 0 left out. Same cuts as the folder.
    // Every folder gets the CE and CE + CE mix fits, the megaCut also the total recon fit, its 50% subsamples and the matrix pages.
    {
        gSystem->mkdir((dir + "/2_phi_cosfit").c_str(), true);
        // the whole-event cuts of this folder (passes) and the fit sample, applied to the rows in memory; the phi spike at 0 is left out and theta must exist (NAN fails)
        const auto inFit = [&](const NewRow& x, int sample){
            const double ph = x.phiRelND * rad;
            return (sample == 0 ? x.evClass == 0 : x.evClass <= sample) && passes(x) && ph > kBlockPhiLo && ph < kBlockPhiHi && x.thetaRelND == x.thetaRelND;
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
        double smearErr = 0; // the A_eff bias of the recon theta: the CE recon events of this folder, mean AP at the recon block theta minus mean AP at the true CE theta; it is an error term of the two recon fits (the CE fit uses the true angle)
        if (ceT){
            ce = Columns(ceT, {"fEvent", "CEPhiNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "CEThetaNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "LastAP", "CENDirX", "CENDirY", "CENDirZ"});
            for (size_t i = 0; i < ce[0].size(); i++) byEvent.emplace((long long)ce[0][i], i);
            const double apScale = ReadSettings().at("APScale");
            double nA = 0, sRec = 0, sTrue = 0;
            for (const NewRow& x : rows) if (inFit(x, 0)){
                const auto it = byEvent.find((long long)x.eventID);
                if (it == byEvent.end()) continue;
                nA++; sRec += AEffOfTheta(x.thetaRelND * rad, apScale); sTrue += AEffOfTheta(ce[2][it->second], apScale);
            }
            if (nA > 0) smearErr = std::fabs(sRec - sTrue) / nA;
            // 1_additionalinfo/2_matrix2 and 3_matrix3 of the megaCut, as plots/6_cuts/<cut>/1_additionalinfo of SmearMatrix.cpp (the CE recon events of this folder, every cut on):
            // matrix2 M_recB_est_to_lab: the recon block angle in the est neutron frame (x, the "truth" of the page) against the recon block angle in the lab frame (y); matrix3 N_true_to_est: the true neutron direction (CENDir) against the est one (the tube centre at the block middle).
            // theta, phi and the flattened (theta slice, phi bin); the probability pages are normalised per column of x, the counts pages are the counts. Same bins as SmearMatrix.cpp (0.1 deg theta bins, 36 phi bins, 10 x 50 flat bins).
            if (mega){
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
        for (int k = 0; k < (mega ? 3 : 2); k++){
            number++;
            if (k == 0 && ceT){
                std::vector<double> ev; // the CE recon events of this folder, in memory
                for (const NewRow& x : rows) if (inFit(x, 0)) ev.push_back(x.eventID);
                const std::string trueFile = Form("cuts/cut_ceTrueNew_%d.root", mode);
                {
                    TFile out(trueFile.c_str(), "RECREATE");
                    double ph = 0, th = 0, ap = 0;
                    TTree tt("CETrue", "CE recon events at the true CE angle");
                    tt.Branch("CEPhiNeuTrue", &ph); tt.Branch("CEThetaNeuTrue", &th); tt.Branch("LastAP", &ap);
                    for (const double e : ev){
                        const auto it = byEvent.find((long long)e);
                        if (it == byEvent.end()) continue; // no carbon-block CE row for this event
                        ph = ce[1][it->second]; th = ce[2][it->second]; ap = ce[3][it->second];
                        tt.Fill();
                    }
                    tt.Write();
                }
                TFile* tf = TFile::Open(trueFile.c_str());
                TTree* tt = tf ? (TTree*)tf->Get("CETrue") : nullptr;
                if (tt){
                    PhiCosFit(&c, tt, "CEPhiNeuTrue", "CEThetaNeuTrue", "", Form("hPhiNew_%d_%d", mode, k),
                              std::string(sampleName[k]) + ", CEPHI_NEUFRAME Distribution, the actual CE angle, A_{eff} from the CE #theta and LastAP" + cutText, Form("%s/2_phi_cosfit/%d_phiBlock_cosfit_%s.pdf", dir.c_str(), k + 1, sampleFile[k]),
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
                for (const NewRow& x : rows){
                    if (!inFit(x, k)) continue;
                    ph = x.phiRelND * rad; th = x.thetaRelND * rad;
                    r1[0] = x.thMin; r1[1] = x.thMax; r1[2] = x.phHalf; r2[0] = x.th2Min; r2[1] = x.th2Max; r2[2] = x.ph2Half;
                    tp.Fill();
                }
                tp.Write();
            }
            TFile* ff = TFile::Open(fitFile.c_str());
            TTree* tp = ff ? (TTree*)ff->Get("PhiFit") : nullptr;
            if (tp){
                PhiCosFit(&c, tp, "phiRelNDdeg", "thetaDeg", "", Form("hPhiNew_%d_%d", mode, k),
                          std::string(sampleName[k]) + ", RECONPHIBL_ESTNEU Distribution, A_{eff} at the recon block #theta" + cutText, Form("%s/2_phi_cosfit/%d_phiBlock_cosfit_%s.pdf", dir.c_str(), k + 1, sampleFile[k]),
                          {}, "", {"rThMin", "rThMax", "rPhHalf"}, {"r2ThMin", "r2ThMax", "r2PhHalf"}, ClassColor(k), smearErr, 33);
                ff->Close();
            }
            gSystem->Unlink(fitFile.c_str());
            made++;
            if (k == 2){ // the all recon fit again for 10 random 50% subsamples (no replacement, fixed seeds), one page each, then a page of the 10 P values, to see the effect of the statistics
                struct SubIn {double ph, a, a1lo, a1hi, a2lo, a2hi, h1, h2;};
                std::vector<SubIn> all;
                const double apScale = ReadSettings().at("APScale"), pTrue = ReadSettings().at("PPValue");
                for (const NewRow& x : rows){
                    if (!inFit(x, 2)) continue;
                    const double th = x.thetaRelND * rad;
                    const bool r1 = std::isfinite(x.thMin) && std::isfinite(x.thMax), r2 = std::isfinite(x.th2Min) && std::isfinite(x.th2Max);
                    const double a = AEffOfTheta(th, apScale);
                    all.push_back({x.phiRelND * rad, a, r1 ? AEffOfTheta(x.thMin, apScale) : a, r1 ? AEffOfTheta(x.thMax, apScale) : a, r2 ? AEffOfTheta(x.th2Min, apScale) : a, r2 ? AEffOfTheta(x.th2Max, apScale) : a, r1 ? x.phHalf : NAN, r2 ? x.ph2Half : NAN});
                }
                const int nSub = 100, nShow = 5, nb = 33; const double frac = 0.5; // 100 subsamples are fitted, the phi plots of the first 5 are drawn, the P page has all 100
                const std::string path = dir + "/2_phi_cosfit/4_phiBlock_cosfit_all_50pct_subsamples.pdf";
                std::vector<double> pSub, eSub; double pFull = 0, eFull = 0;
                if (all.size() >= 10){
                    std::vector<double> phA; AMean aA; PhiRange pr1, pr2;
                    for (const SubIn& s : all){phA.push_back(s.ph); aA.Add(s.a, s.a1lo, s.a1hi, s.a2lo, s.a2hi); if (std::isfinite(s.h1)) pr1.Add(s.h1); if (std::isfinite(s.h2)) pr2.Add(s.h2);}
                    aA.smearErr = smearErr;
                    const PFit uf0 = UnbinnedP(phA, aA, kBlockPhiLo, kBlockPhiHi, nullptr, true, &pr1, &pr2);
                    pFull = uf0.p; eFull = uf0.err;
                    for (int s = 0; s < nSub; s++){
                        std::vector<size_t> idx(all.size());
                        for (size_t i = 0; i < idx.size(); i++) idx[i] = i;
                        TRandom3 rng(12345 + 200 + s);
                        for (size_t i = idx.size() - 1; i > 0; i--) std::swap(idx[i], idx[rng.Integer(i + 1)]);
                        const size_t m = size_t(std::lround(frac * idx.size()));
                        std::vector<double> phS; AMean aM; PhiRange r1, r2;
                        TH1D hS(Form("hPhiSubNew_%d", s), "", nb, 0, 360);
                        for (size_t q = 0; q < m; q++){
                            const SubIn& e = all[idx[q]];
                            phS.push_back(e.ph); hS.Fill(e.ph);
                            aM.Add(e.a, e.a1lo, e.a1hi, e.a2lo, e.a2hi);
                            if (std::isfinite(e.h1)) r1.Add(e.h1);
                            if (std::isfinite(e.h2)) r2.Add(e.h2);
                        }
                        aM.smearErr = smearErr;
                        const PFit uf = UnbinnedP(phS, aM, kBlockPhiLo, kBlockPhiHi, nullptr, true, &r1, &r2);
                        pSub.push_back(uf.p); eSub.push_back(uf.err);
                        if (s >= nShow) continue;
                        c.Clear(); c.cd();
                        c.SetLeftMargin(0.16); c.SetRightMargin(0.25); c.SetTopMargin(0.17); c.SetBottomMargin(0.13);
                        hS.SetStats(0); hS.SetLineColor(ClassColor(2)); hS.SetLineWidth(2); hS.SetMinimum(0); hS.SetMaximum(1.1 * hS.GetMaximum());
                        hS.SetTitle(Form("random 50%% subsample %d of %d (first %d shown), %s, RECONPHIBL_ESTNEU Distribution, %d phi bins of %.1f deg%s;RECONPHIBL_ESTNEU (degrees);Count", s + 1, nSub, nShow, sampleName[2], nb, 360.0 / nb, cutText.c_str()));
                        hS.Draw("hist");
                        TF1 fd("fSubDrawNew", "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                        fd.SetParameters(hS.GetEntries() / nb, uf.p, aM.Mean()); fd.SetLineColor(kBlack); fd.SetLineStyle(2); fd.SetLineWidth(3); fd.Draw("same");
                        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026);
                        tx.DrawLatex(0.18, 0.93, "Unbinned fit  N(#phi) #propto 1 + P #upoint A_{eff} cos#phi, dashed = fit, A_{eff} fixed at the recon block #theta");
                        tx.DrawLatex(0.18, 0.895, Form("%.0f entries (%.0f%% of %zu), %s = %.4f #pm %.4f,  P = %.3g #pm %.2g", hS.GetEntries(), 100.0 * hS.GetEntries() / all.size(), all.size(), AEffEq(), aM.Mean(), aM.Err(), uf.p, uf.err));
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
                    TH1D hP("hSubPNew", Form("P of %d random 50%% subsamples, %s%s;Polarization fit P;Subsamples", nSub, sampleName[2], cutText.c_str()), 20, lo - pad, hi + pad);
                    for (const double p : pSub) hP.Fill(p);
                    hP.SetStats(0); hP.SetLineColor(ClassColor(2)); hP.SetLineWidth(2); hP.SetFillColorAlpha(ClassColor(2), 0.35); hP.SetMinimum(0); hP.SetMaximum(1.3 * hP.GetMaximum());
                    hP.Draw("hist");
                    TLine lf(pFull, 0, pFull, hP.GetMaximum()); lf.SetLineColor(ClassColor(2)); lf.SetLineWidth(3); lf.Draw();
                    TLine lt(pTrue, 0, pTrue, hP.GetMaximum()); lt.SetLineStyle(3); lt.SetLineWidth(3); lt.Draw();
                    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026);
                    tx.DrawLatex(0.14, 0.93, "solid line = P of the full sample, dotted black = true P (PPValue)");
                    tx.DrawLatex(0.14, 0.89, Form("full P = %.3g #pm %.2g, subsample mean %.3g, RMS %.2g, expected RMS %.2g (full error #times #sqrt{(1-f)/f}, f = 0.5)", pFull, eFull, mean, rms, eFull * std::sqrt((1 - frac) / frac)));
                    c.Modified(); c.Update();
                    c.Print((path + ")").c_str());
                    c.Clear(); c.SetLeftMargin(0.13); c.SetRightMargin(0.30); c.SetTopMargin(0.1); c.SetBottomMargin(0.13);
                    made += nShow + 1;
                }
            }
        }
        gStyle->SetOptStat(0);
    }
    // the last plot of the folder: the fraction of the recon events that are CE (CE only, and CE + CE mix) in each bin of the summed KE on the last tracker (the KEEnd plot), same cuts as the folder.
    // Per bin, fraction = events of that kind / all recon events, with the binomial error. The bins hold the events counted on all the threads (MapRows).
    {
        number++;
        const double kLo = ke80 ? keCut : 0.0, kHi = 30000.0;
        const int nb = 60;
        struct FracPart {std::vector<double> ce, mix, all;};
        const auto fracParts = MapRows(n, [&](Long64_t a, Long64_t z){
            FracPart p; p.ce.assign(nb, 0.0); p.mix.assign(nb, 0.0); p.all.assign(nb, 0.0);
            for (Long64_t i = a; i < z; i++){
                const NewRow& r = rows[i];
                if (r.evClass < 0 || r.evClass > 2 || !passes(r)) continue;
                const double x = r.finalKED;
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
        fr[1]->SetTitle(PlotTitle(path, "fraction of the recon events that are CE, per bin of the summed KE on the last tracker" + cutText + ";KE (MeV);fraction of the recon events in the bin").c_str());
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
    printf(">>> PlotReconNew: %d plots in plots/9_ReconNew_1 (1_noCuts, 2_KECUT%s, 3_ to 5_ one cut each, 6_megaCut)\n", made, kePct.c_str());
}
