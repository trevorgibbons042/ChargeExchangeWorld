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
    double outTheta, outDTheta, outPhi, outDPhi, thetaRelND, phiRelND, impactD, r2D, resDRms, finalKED;
    int nMaxPlane, nHits, nPlanesUsed, evClass, nCEUsed, nOtherUsed, nProtonUsed, nPionUsed, nKaonUsed, nNucleusUsed;
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
        {"nMaxChargedPerPlaneEv", "most charged tracks in any one tracker plane of the event (all planes, not only the fit planes);tracks", 10, 0.5, 10.5, true, [](const NewRow& x){return (double)x.nMaxPlane;}}};

    gStyle->SetOptStat(0);
    TCanvas c("cReconNew", "", 900, 650);
    c.SetRightMargin(0.04);
    const Long64_t n = t->GetEntries();
    // one multithreaded read (Columns) of the 19 branches used here (not every branch of the tree), then every plot loop below runs over this vector in memory: 48 passes of GetEntry over all branches took hours at 100M events
    std::vector<NewRow> rows(n);
    NThreads(); // ROOT::EnableImplicitMT: Columns reads the tree on all the threads
    const auto cols = Columns(t, {"outTheta", "outDTheta", "outPhi", "outDPhi", "thetaRelND", "phiRelND", "impactD", "r2D", "resDRms", "finalKED", "nHits", "nPlanesUsed", "evClass", "nCEUsed", "nOtherUsed", "nProtonUsed", "nPionUsed", "nKaonUsed", "nNucleusUsed", "nMaxChargedPerPlaneEv"});
    for (Long64_t i = 0; i < n; i++){
        NewRow& x = rows[i]; size_t k = 0;
        x.outTheta = cols[k++][i]; x.outDTheta = cols[k++][i]; x.outPhi = cols[k++][i]; x.outDPhi = cols[k++][i]; x.thetaRelND = cols[k++][i]; x.phiRelND = cols[k++][i]; x.impactD = cols[k++][i]; x.r2D = cols[k++][i]; x.resDRms = cols[k++][i]; x.finalKED = cols[k++][i];
        x.nHits = (int)cols[k++][i]; x.nPlanesUsed = (int)cols[k++][i]; x.evClass = (int)cols[k++][i]; x.nCEUsed = (int)cols[k++][i]; x.nOtherUsed = (int)cols[k++][i]; x.nProtonUsed = (int)cols[k++][i]; x.nPionUsed = (int)cols[k++][i]; x.nKaonUsed = (int)cols[k++][i]; x.nNucleusUsed = (int)cols[k++][i]; x.nMaxPlane = (int)cols[k++][i];
    }
    const double keCut = kKeCutPercent / 100.0 * ReadSettings().at("BeamKE"); // whole events whose summed (smeared) KE on the last tracker is not above kKeCutPercent of the beam KE are not drawn or counted
    const std::string kePct = Form("%.0f", kKeCutPercent);
    int made = 0;
    // plots/9_ReconNew_1/1_noCuts (every event; its KE plot alone shows the events above the KE cut) and plots/9_ReconNew_1/2_KECUT80 (every plot with the 80% KE cut).
    const int maxPlanes = MaxTrackerHits(); // 2 * nBlocks, the most planes a track can cross
    // 3_KECUTnPlanesUsedEven: the KE cut and a whole-event cut on the planes in the fit, an even number and at most 2 * nBlocks (EvenHits, FitCore.h)
    // 4_KECUTnPlanesUsedEvenTheta0_2: those two cuts and the block fit theta in the est neutron frame (the thetaBlock plot) from 0 to 2 degrees
    for (const int mode : {0, 1, 2, 3}){
    const bool ke80 = mode >= 1, evenPlanes = mode >= 2, theta02 = mode == 3;
    const std::string dir = mode == 0 ? "plots/9_ReconNew_1/1_noCuts" : mode == 1 ? "plots/9_ReconNew_1/2_KECUT" + kePct : mode == 2 ? "plots/9_ReconNew_1/3_KECUTnPlanesUsedEven" : "plots/9_ReconNew_1/4_KECUTnPlanesUsedEvenTheta0_2";
    gSystem->mkdir(dir.c_str(), true);
    // the whole-event cuts of this folder
    const auto passes = [&](const NewRow& x){
        if (ke80 && !(std::isfinite(x.finalKED) && x.finalKED > keCut)) return false;
        if (evenPlanes && !EvenHits(x.nPlanesUsed, maxPlanes)) return false;
        if (theta02 && !(x.thetaRelND * rad >= 0.0 && x.thetaRelND * rad <= 2.0)) return false; // NAN (no tube) fails too
        return true;
    };
    const std::string cutText = theta02 ? ", KE cut above " + kePct + "% of the beam KE, planes in the fit even and at most " + std::to_string(maxPlanes) + ", theta 0 to 2 degrees"
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
        for (Long64_t i = 0; i < n; i++){
            const NewRow& r = rows[i];
            if (r.evClass < 0 || r.evClass > 2) continue;
            const double x = v.value(r);
            if (!passes(r)) continue;
            for (int k = r.evClass; k < 3; k++){ // cumulative: CE only events are in all three, CE mix events in the last two, non-CE events in the last only
                nClass[k]++;
                if (std::isfinite(x)) h[k]->Fill(x);
            }
        }
        double top = 0;
        for (int k = 0; k < 3; k++) top = std::max(top, h[k]->GetMaximum());
        const bool logY = v.logY || (keVar && mode == 0); // the KE plot without cuts is log, in the KE cut folders it is linear
        c.Clear(); c.SetLogy(logY);
        const std::string path = Form("%s/%d_%s.pdf", dir.c_str(), number, v.name);
        // PlotTitle numbers the title from the file name; "no cuts" in the title: the sample is every event with 2 or more planes of exactly 1 hit
        const std::string raw = v.title; const size_t semi = raw.find(';');
        h[0]->SetTitle(PlotTitle(path, raw.substr(0, semi) + (theta02 ? Form(", KE cut above %s%% of the beam KE, planes in the fit even and at most %d, theta 0 to 2 degrees", kePct.c_str(), maxPlanes) : evenPlanes ? Form(", KE cut above %s%% of the beam KE, planes in the fit even and at most %d", kePct.c_str(), maxPlanes) : ke80 ? (", KE cut above " + kePct + "% of the beam KE").c_str() : ", no cuts") + raw.substr(semi)).c_str());
        h[0]->SetMaximum(logY ? top * (std::string(v.name) == "nMaxChargedPerPlaneEv" ? 60.0 : 3.0) : top * 1.25);
        h[0]->SetMinimum(logY ? 0.5 : 0.0);
        h[0]->Draw("hist");
        h[1]->Draw("hist same"); h[2]->Draw("hist same");
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
        for (Long64_t i = 0; i < n; i++){
            const NewRow& r = rows[i];
            if (r.evClass < 0 || r.evClass > 2 || !passes(r)) continue;
            nEv[r.evClass]++;
            b[r.evClass]->AddBinContent(1, r.nCEUsed); b[r.evClass]->AddBinContent(2, r.nProtonUsed);
            b[r.evClass]->AddBinContent(3, r.nPionUsed); b[r.evClass]->AddBinContent(4, r.nKaonUsed); b[r.evClass]->AddBinContent(5, r.nNucleusUsed);
            b[r.evClass]->AddBinContent(6, r.nOtherUsed - r.nProtonUsed - r.nPionUsed - r.nKaonUsed - r.nNucleusUsed);
        }
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
        std::string cutExpr = Form("phiRelND * TMath::RadToDeg() > %g && phiRelND * TMath::RadToDeg() < %g && thetaRelND == thetaRelND", kBlockPhiLo, kBlockPhiHi);
        if (ke80) cutExpr += Form(" && finalKED > %.3f", keCut);
        if (evenPlanes) cutExpr += " && " + EvenHitsDraw("nPlanesUsed");
        if (theta02) cutExpr += " && thetaRelND * TMath::RadToDeg() >= 0 && thetaRelND * TMath::RadToDeg() <= 2";
        c.SetLogy(false); // the canvas still has the log scale of the plot before
        gStyle->SetOptStat(1111); gStyle->SetOptFit(111); // PhiCosFit moves and fills the stats box with the fit, so it needs one
        const char* sampleCut[3] = {"evClass == 0", "evClass <= 1", "evClass <= 2"};
        const char* sampleName[3] = {"CE recon protons", "CE + CE mix recon protons", "total recon particles"};
        const char* sampleFile[3] = {"CE", "CEandMix", "all"};
        const std::vector<std::string> range = {"thetaRelNDMin * TMath::RadToDeg()", "thetaRelNDMax * TMath::RadToDeg()", "phiRelNDHalf"};
        const std::vector<std::string> range2 = {"thetaRelNDBMin * TMath::RadToDeg()", "thetaRelNDBMax * TMath::RadToDeg()", "phiRelNDBHalf"};
        // The CE sample is fitted at the actual CE angle (ChargeExchange_values: the true phi and theta of the CE proton in the true neutron frame, CEPHI_NEUFRAME) and its own LastAP, not at the recon angle in the est neutron frame:
        // the events are the CE recon events that pass this folder's cuts (evClass 0 and cutExpr, matched by fEvent), only the angle of each is the true one. The rows go to a small file so Columns reads them on all the threads.
        TTree* ceT = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
        gROOT->cd();
        for (int k = 0; k < 3; k++){
            number++;
            if (k == 0 && ceT){
                const auto ce = Columns(ceT, {"fEvent", "CEPhiNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "CEThetaNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "LastAP"});
                std::map<long long, size_t> byEvent;
                for (size_t i = 0; i < ce[0].size(); i++) byEvent.emplace((long long)ce[0][i], i);
                const auto ev = Columns(t, {"fEvent"}, std::string(sampleCut[0]) + " && " + cutExpr)[0];
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
            PhiCosFit(&c, t, "phiRelND * TMath::RadToDeg()", "thetaRelND * TMath::RadToDeg()", std::string(sampleCut[k]) + " && " + cutExpr, Form("hPhiNew_%d_%d", mode, k),
                      std::string(sampleName[k]) + ", RECONPHIBL_ESTNEU Distribution, A_{eff} at the recon block #theta" + cutText, Form("%s/%d_phiBlock_cosfit_%s.pdf", dir.c_str(), number, sampleFile[k]),
                      {}, "", range, range2, ClassColor(k));
            made++;
        }
        gStyle->SetOptStat(0);
    }
    }
    printf(">>> PlotReconNew: %d plots in plots/9_ReconNew_1 (1_noCuts, 2_KECUT%s, 3_ and 4_)\n", made, kePct.c_str());
}
