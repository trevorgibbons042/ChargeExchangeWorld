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

// 9_ReconNew: the event recon of ReconAngleAll.cpp(true) (cuts/cut_angles_new.root): one fit per event through the tracker planes with exactly 1 charged hit, no cut of any kind.
// Every plot overlays the three kinds of event (evClass): CE only (every fitted hit a CE child proton), CE mix (CE child hits and other particles), non-CE (no CE child hit).
// Output: plots/9_ReconNew/1_noCuts/<n>_<name>.pdf, titled "no cuts".
namespace {
struct NewRow {
    double thMin, thMax, phHalf, th2Min, th2Max, ph2Half, outTheta, outDTheta, outPhi, outDPhi, thetaRelND, phiRelND, impactD, r2D, resDRms, finalKED;
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
    const auto cols = ColumnsUnordered(t, {"outTheta", "outDTheta", "outPhi", "outDPhi", "thetaRelND", "phiRelND", "impactD", "r2D", "resDRms", "finalKED", "nHits", "nPlanesUsed", "evClass", "nCEUsed", "nOtherUsed", "nProtonUsed", "nPionUsed", "nKaonUsed", "nNucleusUsed", "nMaxChargedPerPlaneEv", "nTrackerHitsEv", "nPlanesHitEv", "thetaRelNDMin", "thetaRelNDMax", "phiRelNDHalf", "thetaRelNDBMin", "thetaRelNDBMax", "phiRelNDBHalf", "fEvent"});
    for (Long64_t i = 0; i < n; i++){
        NewRow& x = rows[i]; size_t k = 0;
        x.outTheta = cols[k++][i]; x.outDTheta = cols[k++][i]; x.outPhi = cols[k++][i]; x.outDPhi = cols[k++][i]; x.thetaRelND = cols[k++][i]; x.phiRelND = cols[k++][i]; x.impactD = cols[k++][i]; x.r2D = cols[k++][i]; x.resDRms = cols[k++][i]; x.finalKED = cols[k++][i];
        x.nHits = (int)cols[k++][i]; x.nPlanesUsed = (int)cols[k++][i]; x.evClass = (int)cols[k++][i]; x.nCEUsed = (int)cols[k++][i]; x.nOtherUsed = (int)cols[k++][i]; x.nProtonUsed = (int)cols[k++][i]; x.nPionUsed = (int)cols[k++][i]; x.nKaonUsed = (int)cols[k++][i]; x.nNucleusUsed = (int)cols[k++][i]; x.nMaxPlane = (int)cols[k++][i]; x.nTrackerHits = (int)cols[k++][i]; x.nPlanesHit = (int)cols[k++][i];
        x.thMin = cols[k++][i] * rad; x.thMax = cols[k++][i] * rad; x.phHalf = cols[k++][i]; x.th2Min = cols[k++][i] * rad; x.th2Max = cols[k++][i] * rad; x.ph2Half = cols[k++][i]; x.eventID = (int)cols[k++][i];
    }
    const double keCut = kKeCutPercent / 100.0 * ReadSettings().at("BeamKE"); // whole events whose summed (smeared) KE on the last tracker is not above kKeCutPercent of the beam KE are not drawn or counted
    const std::string kePct = Form("%.0f", kKeCutPercent);
    int made = 0;
    // plots/9_ReconNew_1/1_noCuts (every event; its KE plot alone shows the events above the KE cut) and plots/9_ReconNew_1/2_KECUT80 (every plot with the 80% KE cut).
    const int maxPlanes = MaxTrackerHits(); // 2 * nBlocks, the most planes a track can cross
    // 3_KECUTnPlanesUsedEven: the KE cut and a whole-event cut on the planes in the fit, an even number and at most 2 * nBlocks (EvenHits, FitCore.h)
    // 4_KECUTnPlanesUsedEvenTheta0_1: those two cuts and the block fit theta in the est neutron frame (the thetaBlock plot) from 0 to 1 degree
    for (const int mode : {0, 1, 2, 3, 4}){
    const bool ke80 = mode >= 1, evenPlanes = mode >= 2, theta01 = mode >= 3, avgOne = mode == 4; // 5_: the cuts of 4_ and on average one charged tracker hit per hit plane (nTrackerHitsEv == planes hit, the avgTrackerHitsPerPlane plot at 1)
    const std::string dir = mode == 0 ? "plots/9_ReconNew_1/1_noCuts" : mode == 1 ? "plots/9_ReconNew_1/2_KECUT" + kePct : mode == 2 ? "plots/9_ReconNew_1/3_KECUTnPlanesUsedEven" : mode == 3 ? "plots/9_ReconNew_1/4_KECUTnPlanesUsedEvenTheta0_1" : "plots/9_ReconNew_1/5_KECUTnPlanesUsedEvenTheta0_1avgTrackerHitsPerPlane1";
    gSystem->mkdir(dir.c_str(), true);
    // the whole-event cuts of this folder
    const auto passes = [&](const NewRow& x){
        if (ke80 && !(std::isfinite(x.finalKED) && x.finalKED > keCut)) return false;
        if (evenPlanes && !EvenHits(x.nPlanesUsed, maxPlanes)) return false;
        if (theta01 && !(x.thetaRelND * rad >= 0.0 && x.thetaRelND * rad <= 1.0)) return false; // NAN (no tube) fails too
        if (avgOne && !(x.nPlanesHit > 0 && x.nTrackerHits == x.nPlanesHit)) return false;
        return true;
    };
    const std::string cutText = avgOne ? ", KE cut above " + kePct + "% of the beam KE, planes in the fit even and at most " + std::to_string(maxPlanes) + ", theta 0 to 1 degree, average tracker hits per hit plane = 1"
        : theta01 ? ", KE cut above " + kePct + "% of the beam KE, planes in the fit even and at most " + std::to_string(maxPlanes) + ", theta 0 to 1 degree"
        : evenPlanes ? ", KE cut above " + kePct + "% of the beam KE, planes in the fit even and at most " + std::to_string(maxPlanes)
        : ke80 ? ", KE cut above " + kePct + "% of the beam KE" : ", no cuts";
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
        const bool logY = v.logY || (keVar && mode == 0); // the KE plot without cuts is log, in the KE cut folders it is linear
        c.Clear(); c.SetLogy(logY);
        const std::string path = Form("%s/%d_%s.pdf", dir.c_str(), number, v.name);
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
        const std::string path = Form("%s/%d_particleTypes.pdf", dir.c_str(), number);
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
    {
        // the whole-event cuts of this folder (passes) and the fit sample, applied to the rows in memory; the phi spike at 0 is left out and theta must exist (NAN fails)
        const auto inFit = [&](const NewRow& x, int sample){
            const double ph = x.phiRelND * rad;
            return (sample == 0 ? x.evClass == 0 : x.evClass <= sample) && passes(x) && ph > kBlockPhiLo && ph < kBlockPhiHi && x.thetaRelND == x.thetaRelND;
        };
        c.SetLogy(false); // the canvas still has the log scale of the plot before
        gStyle->SetOptStat(1111); gStyle->SetOptFit(111); // PhiCosFit moves and fills the stats box with the fit, so it needs one
                const char* sampleName[3] = {"CE recon protons", "CE + CE mix recon protons", "total recon particles"};
        const char* sampleFile[3] = {"CE", "CEandMix", "all"};
        // The CE sample is fitted at the actual CE angle (ChargeExchange_values: the true phi and theta of the CE proton in the true neutron frame, CEPHI_NEUFRAME) and its own LastAP, not at the recon angle in the est neutron frame:
        // the events are the CE recon events that pass this folder's cuts (evClass 0 and cutExpr, matched by fEvent), only the angle of each is the true one. The rows go to a small file so Columns reads them on all the threads.
        TTree* ceT = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
        gROOT->cd();
        std::vector<std::vector<double>> ce;
        std::map<long long, size_t> byEvent;
        double phiD = 1; // the dilution of P by the smeared recon phi: mean cos(recon phi - true phi) of those CE events, an error term of P of the two recon fits
        double smearErr = 0; // the A_eff bias of the recon theta: the CE recon events of this folder, mean AP at the recon block theta minus mean AP at the true CE theta; it is an error term of the two recon fits (the CE fit uses the true angle)
        if (ceT){
            ce = Columns(ceT, {"fEvent", "CEPhiNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "CEThetaNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "LastAP"});
            for (size_t i = 0; i < ce[0].size(); i++) byEvent.emplace((long long)ce[0][i], i);
            const double apScale = ReadSettings().at("APScale");
            double nA = 0, sRec = 0, sTrue = 0, sCos = 0;
            for (const NewRow& x : rows) if (inFit(x, 0)){
                const auto it = byEvent.find((long long)x.eventID);
                if (it == byEvent.end()) continue;
                nA++; sRec += AEffOfTheta(x.thetaRelND * rad, apScale); sTrue += AEffOfTheta(ce[2][it->second], apScale); sCos += std::cos((x.phiRelND * rad - ce[1][it->second]) * M_PI / 180.0);
            }
            if (nA > 0){smearErr = std::fabs(sRec - sTrue) / nA; phiD = sCos / nA;}
        }
        for (int k = 0; k < 3; k++){
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
                              std::string(sampleName[k]) + ", CEPHI_NEUFRAME Distribution, the actual CE angle, A_{eff} from the CE #theta and LastAP" + cutText, Form("%s/%d_phiBlock_cosfit_%s.pdf", dir.c_str(), number, sampleFile[k]),
                              {}, "LastAP", {}, {}, ClassColor(k));
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
                          std::string(sampleName[k]) + ", RECONPHIBL_ESTNEU Distribution, A_{eff} at the recon block #theta" + cutText, Form("%s/%d_phiBlock_cosfit_%s.pdf", dir.c_str(), number, sampleFile[k]),
                          {}, "", {"rThMin", "rThMax", "rPhHalf"}, {"r2ThMin", "r2ThMax", "r2PhHalf"}, ClassColor(k), smearErr, phiD);
                ff->Close();
            }
            gSystem->Unlink(fitFile.c_str());
            made++;
        }
        gStyle->SetOptStat(0);
    }
    }
    printf(">>> PlotReconNew: %d plots in plots/9_ReconNew_1 (1_noCuts, 2_KECUT%s, 3_, 4_ and 5_)\n", made, kePct.c_str());
}
