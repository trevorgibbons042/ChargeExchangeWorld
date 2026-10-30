#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"

#include "BlockLimits.h"
#include "CutOpen.h"
#include "FitCore.h"
#include "QualityCut.h"
#include "TreeColumns.h"

// CE of the neutrons in the carbon blocks, split into the beam neutron (CECreator "primary", TrackID 1) and the secondary neutrons made by hadronic interactions.
// 3_CEpolY_distribution.pdf: CEPolY of each group, CEPol is the spin the sim carries on the CE track (track->GetPolarization(), unit vectors, the beam
// PPValue is the fraction along +y; Geant4 gives secondaries none). 2_CE_count_by_creator.pdf: CE count per CECreator, the process that made the CE neutron.
// Needs ChargeExchange_values with CECreator and CEPolX/Y/Z (sim of 27 Oct 2026 or later) in output0_t*.root and settings.csv.
void SecondaryCE(){
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    gStyle->SetOptStat(0);
    // the carbon-block CE (CarbonCeCut) come from cuts/cut_CE.root (CutOpen.h), so every row passes; the ZDC protons of their children and the first tracker protons likewise
    // (children of a CE track: only those match a groupOf key below)
    TTree* cePtr = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
    if (!cePtr){Error("SecondaryCE", "no cuts/cut_CE.root (run CutEvents.cpp)"); return;}
    TTree& ce = *cePtr;
    gROOT->cd();
    if (!ce.GetBranch("CECreator") || !ce.GetBranch("CEPolX")){Error("SecondaryCE", "CECreator/CEPolX missing in ChargeExchange_values (rerun the sim)"); return;}
    const std::string outDir = "plots/7_secondaryCE/";
    for (const char* sub : {"1_counts", "2_polarization", "3_energy", "4_position", "5_phi_cosfit"}) gSystem->mkdir((outDir + sub).c_str(), true); // one folder per kind of plot

    // carbon-block CE only, like the other CE plots (CarbonCeCut in FitCore.h)
    // ce is read once, multithreaded, into one vector per column (Columns, TreeColumns.h; CECreator is a char leaf, read as a std::string), every row of the cache in entry order
    const bool hasVz = ce.GetBranch("CEVertexZ");
    std::vector<std::string> ceExprs = {"fEvent", "TrackID", "CEPolY", "LastAP", "LastPhiLab", "fX", "fY"};
    if (hasVz) ceExprs.push_back("CEVertexZ");
    const auto ceCol = Columns(&ce, ceExprs);
    const auto& ceEvent = ceCol[0], &ceTrack = ceCol[1], &cePolY = ceCol[2], &ceAP = ceCol[3], &cePhiLab = ceCol[4], &ceX = ceCol[5], &ceYv = ceCol[6];
    const auto ceCreator = Columns<std::string>(&ce, {"std::string(CECreator.data())"})[0];
    const Long64_t nCE = ceEvent.size();
    const auto creatorGroup = [](const std::string& c){return c == "primary" ? 0 : c == "neutronInelastic" ? 1 : c == "protonInelastic" ? 2 : 3;};
    std::map<std::pair<int, int>, int> groupOf; // (fEvent, CE neutron TrackID) -> 0 beam neutron, 1 neutronInelastic, 2 protonInelastic, 3 every other creator
    double sA[2] = {0, 0}, sCos[2] = {0, 0}; // model A and cos(lab phi) sums of the beam neutron CE (0) and the secondary CE (1), for the P estimate
    TH1D hYprim("hYprim", "", 42, -1.05, 1.05), hYsec("hYsec", "", 42, -1.05, 1.05); // CEPolY of the beam neutron and of the secondary neutrons
    std::map<std::string, long long> byCreator;
    long long nPrim = 0, nSec = 0;
    for (Long64_t k = 0; k < nCE; k++){
        const int fEvent = (int)ceEvent[k], trackID = (int)ceTrack[k];
        const int g = trackID == 1 ? 0 : 1;
        (g == 0 ? hYprim : hYsec).Fill(cePolY[k]);
        sA[g] += ceAP[k]; sCos[g] += std::cos(cePhiLab[k]);
        if (g == 0) nPrim++; else nSec++;
        byCreator[ceCreator[k]]++;
        groupOf[{fEvent, trackID}] = creatorGroup(ceCreator[k]);
    }
    const double fracSec = 100.0 * nSec / std::max<long long>(1, nPrim + nSec);

    // same look as the other plots: 1500 x 1200 canvas, right margin 0.30 for the legend and stats, line histograms, PlotTitle, stats.csv
    gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetNdivisions(505, "X");
    TCanvas canvas("cSecCE", "Plots", 1500, 1200);
    canvas.SetRightMargin(0.30);
    std::set<std::string> numbersStarted;
    const std::string creatorFile = "1_counts/2_CE_count_by_creator.pdf";

    // y component of the CE neutron polarization, every CE: beam neutron (black) and secondary neutrons (red), log y
    const std::string polYFile = "2_polarization/3_CEpolY_distribution.pdf";
    canvas.SetLogy(1);
    hYprim.SetLineColor(kBlack); hYsec.SetLineColor(kRed);
    hYprim.SetLineWidth(2); hYsec.SetLineWidth(2); hYprim.SetStats(0);
    hYprim.SetMinimum(0.5);
    hYprim.SetMaximum(2.0 * std::max(hYprim.GetMaximum(), hYsec.GetMaximum()));
    hYprim.SetTitle(PlotTitle(polYFile, "y component of the polarization of the neutron doing CE in the carbon blocks;CEPolY;Count").c_str());
    hYprim.Draw("hist"); hYsec.Draw("hist same");
    TLegend legY(0.745, 0.75, 0.99, 0.88); legY.SetTextSize(0.016);
    legY.AddEntry(&hYprim, Form("beam neutron (%lld), mean %.3f", nPrim, hYprim.GetMean()), "l");
    legY.AddEntry(&hYsec, Form("secondary (%lld), mean %.3f", nSec, hYsec.GetMean()), "l");
    legY.Draw();
    WritePlotNumbers(&hYprim, outDir + polYFile, numbersStarted);
    canvas.SaveAs((outDir + polYFile).c_str());
    canvas.SetLogy(0);

    const int nCreators = byCreator.size();
    TH1D hCreator("hCreator", "", nCreators, 0, nCreators);
    int bin = 1;
    // tallest bar on the left (the beam neutron, "primary", is the tallest), then down by count
    std::vector<std::pair<std::string, long long>> ordered(byCreator.begin(), byCreator.end());
    std::sort(ordered.begin(), ordered.end(), [](const auto& l, const auto& r){return l.second > r.second;});
    for (const auto& [name, n] : ordered){hCreator.SetBinContent(bin, n); hCreator.GetXaxis()->SetBinLabel(bin, name.c_str()); bin++;}
    hCreator.SetLineColor(kBlack); hCreator.SetLineWidth(2); hCreator.SetStats(0);
    hCreator.SetMinimum(0.5);
    hCreator.GetXaxis()->SetLabelSize(0.025); hCreator.GetXaxis()->SetTitleOffset(1.6); // creator names are long, keep them clear of the axis title
    canvas.SetBottomMargin(0.15);
    hCreator.SetTitle(PlotTitle(creatorFile, Form("CE in the carbon blocks by the process that made the CE neutron, %.1f%% from secondary neutrons;Creator process;Count", fracSec)).c_str());
    canvas.Clear(); canvas.SetLogy(1);
    hCreator.Draw("hist");
    WritePlotNumbers(&hCreator, outDir + creatorFile, numbersStarted);
    canvas.SaveAs((outDir + creatorFile).c_str());

    // 4_CE_count_primary_vs_secondary.pdf: the beam neutron on the left, every other neutron source added up on the right
    const std::string splitFile = "1_counts/4_CE_count_primary_vs_secondary.pdf";
    TH1D hSplit("hSplit", "", 2, 0, 2);
    hSplit.SetBinContent(1, nPrim); hSplit.GetXaxis()->SetBinLabel(1, "primary");
    hSplit.SetBinContent(2, nSec); hSplit.GetXaxis()->SetBinLabel(2, "all other neutron sources");
    hSplit.SetLineColor(kBlack); hSplit.SetLineWidth(2); hSplit.SetStats(0);
    hSplit.SetMinimum(0.5);
    hSplit.GetXaxis()->SetLabelSize(0.03); hSplit.GetXaxis()->SetTitleOffset(1.6);
    hSplit.SetTitle(PlotTitle(splitFile, Form("CE in the carbon blocks, beam neutron and all other neutron sources, %.1f%% from the other sources;Neutron source;Count", fracSec)).c_str());
    canvas.Clear(); canvas.SetLogy(1);
    hSplit.Draw("hist");
    WritePlotNumbers(&hSplit, outDir + splitFile, numbersStarted);
    canvas.SaveAs((outDir + splitFile).c_str());

    // 5_P_fit_if_inelastic_fixed.pdf: fitted P of the CE protons now, and what it would be if the neutrons made by neutronInelastic (and protonInelastic) were polarized like the beam neutron.
    // P = 2 <cos phi> / <A> (the unweighted cos fit with fixed <A>, same as 3_ReconCalcVals/1_phiCE_cosfit.pdf). The two "fixed" bars are estimates: the CE of those creators
    // count as polarized like the beam-neutron CE, P_est = P(beam neutron CE) * (N_beam + N_fixed) / N_all.
    const std::string pFile = "2_polarization/5_P_fit_if_inelastic_fixed.pdf";
    const double nAll = nPrim + nSec, pPrim = std::clamp(2 * sCos[0] / sA[0], -1.0, 1.0), pAll = std::clamp(2 * (sCos[0] + sCos[1]) / (sA[0] + sA[1]), -1.0, 1.0); // P is a polarization
    const double errPrim = std::sqrt(2.0 / nPrim) / std::fabs(sA[0] / nPrim), errAll = std::sqrt(2.0 / nAll) / std::fabs((sA[0] + sA[1]) / nAll);
    const long long nNeutronInel = byCreator.count("neutronInelastic") ? byCreator["neutronInelastic"] : 0, nProtonInel = byCreator.count("protonInelastic") ? byCreator["protonInelastic"] : 0;
    TH1D hP("hPfix", "", 5, 0, 5);
    const char* pName[5] = {"input PPValue", "now, all CE", "neutronInelastic fixed (est.)", "neutron + protonInelastic fixed (est.)", "beam neutron CE only"};
    const double pVal[5] = {ReadSettings().at("PPValue"), pAll, pPrim * (nPrim + nNeutronInel) / nAll, pPrim * (nPrim + nNeutronInel + nProtonInel) / nAll, pPrim};
    const double pErr[5] = {0, errAll, errPrim * (nPrim + nNeutronInel) / nAll, errPrim * (nPrim + nNeutronInel + nProtonInel) / nAll, errPrim};
    for (int i = 0; i < 5; i++){hP.SetBinContent(i + 1, pVal[i]); hP.SetBinError(i + 1, pErr[i]); hP.GetXaxis()->SetBinLabel(i + 1, pName[i]);}
    hP.SetLineColor(kBlack); hP.SetLineWidth(2); hP.SetStats(0);
    hP.SetMinimum(0); hP.SetMaximum(1.0);
    hP.GetXaxis()->SetLabelSize(0.022); hP.GetXaxis()->SetTitleOffset(1.6);
    hP.SetTitle(PlotTitle(pFile, "Fitted polarization P of the CE protons, now and if the neutronInelastic and protonInelastic neutrons were polarized like the beam neutron;;P").c_str());
    canvas.Clear(); canvas.SetLogy(0);
    hP.Draw("hist");
    TLine pRef(0, ReadSettings().at("PPValue"), 5, ReadSettings().at("PPValue")); pRef.SetLineStyle(3); pRef.SetLineWidth(2); pRef.SetLineColor(kRed); pRef.Draw(); // input PPValue
    WritePlotNumbers(&hP, outDir + pFile, numbersStarted);
    canvas.SaveAs((outDir + pFile).c_str());
    // 6_CE_proton_KE_by_creator.pdf and 7_CE_proton_smeared_KE_by_creator.pdf: KE at the ZDC of the protons made by the CE, true (fKEnergy) and with the ZDC energy smearing
    // (fDEnergy, /zdc/resStochastic and /zdc/resConstant in the sim, equal to fKEnergy when the smearing is off), one curve per process that made the CE neutron.
    // The proton's ParentID is the CE neutron's TrackID (like PlotDiffPhiCuts.cpp), matched on (fEvent, ParentID).
    const char* gName[4] = {"beam neutron", "neutronInelastic", "protonInelastic", "all other creators"};
    const int gColor[4] = {kBlack, kRed, kBlue, kGreen + 2};
    TH1D hKE[2][4]; // [0] true KE, [1] smeared KE
    for (int q = 0; q < 2; q++) for (int g = 0; g < 4; g++) hKE[q][g] = TH1D(Form("hKE%d_%d", q, g), "", 100, 0, 30000);
    TTree* piPtr = OpenCut("cuts/cut_ParentZdc.root", "ParticleInfo");
    if (!piPtr){Error("SecondaryCE", "no cuts/cut_ParentZdc.root (run CutEvents.cpp)"); return;}
    TTree& pi = *piPtr;
    gROOT->cd();
    // pi is read once, multithreaded, every row of the cache in entry order: the ZDC protons of the CE children
    const auto piCol = Columns(&pi, {"fEvent", "ParentID", "TrackID", "fKEnergy", "fDEnergy"});
    const Long64_t nPi = piCol[0].size();
    const double cutKE = 0.80 * kBeamKE; // MeV, the 80% KE cut of the other plots (initialKE * 0.80, Plot3D.cpp), on the smeared energy fDEnergy
    long long nAbove[4] = {0, 0, 0, 0}; // CE protons of each creator above the cut
    for (Long64_t k = 0; k < nPi; k++){
        const int piEvent = (int)piCol[0][k], piParent = (int)piCol[1][k]; const double piKE = piCol[3][k], piDE = piCol[4][k];
        const auto it = groupOf.find({piEvent, piParent});
        if (it != groupOf.end()){hKE[0][it->second].Fill(piKE); hKE[1][it->second].Fill(piDE); if (piDE > cutKE) nAbove[it->second]++;}
    }
    const std::string keFile[2] = {"3_energy/6_CE_proton_KE_by_creator.pdf", "3_energy/7_CE_proton_smeared_KE_by_creator.pdf"};
    const char* keTitle[2] = {"True Kinetic Energy at ZDC of the CE protons, by the process that made the CE neutron;Kinetic Energy (MeV);Count",
                              "Smeared Energy at ZDC of the CE protons, by the process that made the CE neutron;Smeared Energy fDEnergy (MeV);Count"};
    for (int q = 0; q < 2; q++){
        TLegend legKE(0.745, 0.65, 0.99, 0.88); legKE.SetTextSize(0.016);
        double keMax = 0;
        for (int g = 0; g < 4; g++){hKE[q][g].SetLineColor(gColor[g]); hKE[q][g].SetLineWidth(2); hKE[q][g].SetStats(0); keMax = std::max(keMax, hKE[q][g].GetMaximum());}
        hKE[q][0].SetMinimum(0.5); hKE[q][0].SetMaximum(2.0 * keMax);
        hKE[q][0].SetTitle(PlotTitle(keFile[q], keTitle[q]).c_str());
        canvas.Clear(); canvas.SetLogy(1);
        for (int g = 0; g < 4; g++){hKE[q][g].Draw(g ? "hist same" : "hist"); legKE.AddEntry(&hKE[q][g], Form("%s (%.0f)", gName[g], hKE[q][g].GetEntries()), "l");}
        legKE.Draw();
        WritePlotNumbers(&hKE[q][0], outDir + keFile[q], numbersStarted);
        canvas.SaveAs((outDir + keFile[q]).c_str());
        canvas.SetLogy(0);
    }
    // 8_CE_proton_smeared_KE_cut80_by_creator.pdf: the smeared energy plot with the 80% KE cut (dashed), how many CE protons of each creator stay above it
    const std::string cutFile = "3_energy/8_CE_proton_smeared_KE_cut80_by_creator.pdf";
    canvas.Clear(); canvas.SetLogy(1);
    hKE[1][0].SetTitle(PlotTitle(cutFile, Form("Smeared Energy at ZDC of the CE protons, by the process that made the CE neutron, 80%% KE cut = %.0f MeV;Smeared Energy fDEnergy (MeV);Count", cutKE)).c_str());
    TLegend legCut(0.745, 0.55, 0.99, 0.88); legCut.SetTextSize(0.016);
    for (int g = 0; g < 4; g++){
        hKE[1][g].Draw(g ? "hist same" : "hist");
        legCut.AddEntry(&hKE[1][g], Form("%s %lld/%.0f (%.0f%%)", gName[g], nAbove[g], hKE[1][g].GetEntries(), 100.0 * nAbove[g] / std::max(1.0, hKE[1][g].GetEntries())), "l");
    }
    TLine cutLine(cutKE, 0.5, cutKE, hKE[1][0].GetMaximum()); cutLine.SetLineStyle(2); cutLine.SetLineWidth(2); cutLine.Draw();
    legCut.Draw();
    WritePlotNumbers(&hKE[1][0], outDir + cutFile, numbersStarted);
    canvas.SaveAs((outDir + cutFile).c_str());
    canvas.SetLogy(0);
    // 9_CE_x_position_by_creator.pdf: x of the neutron at the CE (fX, mm), one curve per process that made the CE neutron, same groups as the KE plots
    const std::string xFile = "4_position/9_CE_x_position_by_creator.pdf";
    std::vector<double> xs[4];
    double xLo = 1e30, xHi = -1e30;
    std::map<std::pair<int, int>, double> xOf; // (fEvent, CE neutron TrackID) -> x at the CE
    for (Long64_t k = 0; k < nCE; k++){
        const int g = creatorGroup(ceCreator[k]);
        const double fX = ceX[k];
        xs[g].push_back(fX); xLo = std::min(xLo, fX); xHi = std::max(xHi, fX);
        xOf[{(int)ceEvent[k], (int)ceTrack[k]}] = fX;
    }
    const double xPad = 0.02 * (xHi - xLo) + 1e-6;
    TH1D hX[4];
    double xMax = 0;
    for (int g = 0; g < 4; g++){
        hX[g] = TH1D(Form("hX%d", g), "", 100, xLo - xPad, xHi + xPad);
        for (double x : xs[g]) hX[g].Fill(x);
        hX[g].SetLineColor(gColor[g]); hX[g].SetLineWidth(2); hX[g].SetStats(0); xMax = std::max(xMax, hX[g].GetMaximum());
    }
    hX[0].SetMinimum(0.5); hX[0].SetMaximum(2.0 * xMax);
    hX[0].SetTitle(PlotTitle(xFile, "x position of the neutron at the CE in the carbon blocks, by the process that made the CE neutron;x (mm);Count").c_str());
    canvas.Clear(); canvas.SetLogy(1);
    TLegend legX(0.745, 0.65, 0.99, 0.88); legX.SetTextSize(0.016);
    for (int g = 0; g < 4; g++){hX[g].Draw(g ? "hist same" : "hist"); legX.AddEntry(&hX[g], Form("%s (%.0f)", gName[g], hX[g].GetEntries()), "l");}
    legX.Draw();
    WritePlotNumbers(&hX[0], outDir + xFile, numbersStarted);
    canvas.SaveAs((outDir + xFile).c_str());
    // 10_CE_x_position_cut80_by_creator.pdf: the same x plot for the CE whose proton ends above the 80% KE cut on the smeared ZDC energy (plot 8)
    const std::string xCutFile = "4_position/10_CE_x_position_cut80_by_creator.pdf";
    TH1D hXc[4];
    double xcMax = 0;
    for (int g = 0; g < 4; g++) hXc[g] = TH1D(Form("hXc%d", g), "", 100, xLo - xPad, xHi + xPad);
    for (Long64_t k = 0; k < nPi; k++){
        const int piEvent = (int)piCol[0][k], piParent = (int)piCol[1][k]; const double piDE = piCol[4][k];
        const auto it = groupOf.find({piEvent, piParent});
        if (it != groupOf.end() && piDE > cutKE) hXc[it->second].Fill(xOf.at({piEvent, piParent}));
    }
    for (int g = 0; g < 4; g++){hXc[g].SetLineColor(gColor[g]); hXc[g].SetLineWidth(2); hXc[g].SetStats(0); xcMax = std::max(xcMax, hXc[g].GetMaximum());}
    hXc[0].SetMinimum(0.5); hXc[0].SetMaximum(2.0 * xcMax);
    hXc[0].SetTitle(PlotTitle(xCutFile, Form("x position of the neutron at the CE, by the process that made the CE neutron, 80%% KE cut = %.0f MeV;x (mm);Count", cutKE)).c_str());
    canvas.Clear();
    TLegend legXc(0.745, 0.65, 0.99, 0.88); legXc.SetTextSize(0.016);
    for (int g = 0; g < 4; g++){hXc[g].Draw(g ? "hist same" : "hist"); legXc.AddEntry(&hXc[g], Form("%s %.0f/%.0f (%.0f%%)", gName[g], hXc[g].GetEntries(), hX[g].GetEntries(), 100.0 * hXc[g].GetEntries() / std::max(1.0, hX[g].GetEntries())), "l");}
    legXc.Draw();
    WritePlotNumbers(&hXc[0], outDir + xCutFile, numbersStarted);
    canvas.SaveAs((outDir + xCutFile).c_str());
    // 18_phiCE_cosfit_beam_vs_all.pdf: the cos phi fit of 3_ReconCalcVals/1_phiCE_cosfit.pdf (N(phi) = N0 (1 + P <A> cos phi), N0 = entries / bins and <A> fixed, only P fitted, 50 bins of the CE
    // proton lab phi) for the CE of the beam neutron only (TrackID 1, red) and for every carbon-block CE (black)
    // keep(k) picks the CE rows (it reads the vectors of the ce entry k), what is added to the title; 19_ is the same plot for the CE protons that hit the first tracker inside a square
    const auto phiFitPlot = [&](const std::string& phiFile, const std::string& what, const std::function<bool(Long64_t)>& keep){
        TH1D hPhi[2] = {TH1D("hPhiAll", "", 50, 0, 360), TH1D("hPhiBeam", "", 50, 0, 360)};
        double sAng[2] = {0, 0}; long long nPhi[2] = {0, 0};
        for (Long64_t k = 0; k < nCE; k++){
            if (!keep(k)) continue;
            const double phiRad = cePhiLab[k], lastAP = ceAP[k]; const int trackID = (int)ceTrack[k];
            const double phiDeg = (phiRad < 0 ? phiRad + TMath::TwoPi() : phiRad) * 180 / TMath::Pi();
            for (int g = 0; g < 2; g++){
                if (g == 1 && trackID != 1) continue;
                hPhi[g].Fill(phiDeg); sAng[g] += lastAP; nPhi[g]++;
            }
        }
        const char* phiName[2] = {"all CE", "beam neutron CE (TrackID 1)"};
        const int phiColor[2] = {kBlack, kRed};
        TF1 phiFit[2];
        double pFit[2], pFitErr[2];
        double phiMax = 0;
        for (int g = 0; g < 2; g++){
            phiFit[g] = TF1(Form("fPhiCos%d", g), "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
            phiFit[g].FixParameter(0, hPhi[g].GetEntries() / 50); phiFit[g].FixParameter(2, sAng[g] / std::max<long long>(1, nPhi[g])); phiFit[g].SetParameter(1, 0.5); phiFit[g].SetParLimits(1, -1, 1);
            hPhi[g].Fit(&phiFit[g], "Q0");
            pFit[g] = phiFit[g].GetParameter(1); pFitErr[g] = phiFit[g].GetParError(1);
            hPhi[g].SetLineColor(phiColor[g]); hPhi[g].SetLineWidth(2); hPhi[g].SetStats(0);
            phiFit[g].SetLineColor(phiColor[g]); phiFit[g].SetLineStyle(2); phiFit[g].SetLineWidth(2);
            phiMax = std::max(phiMax, hPhi[g].GetMaximum());
        }
        hPhi[0].SetMinimum(0); hPhi[0].SetMaximum(1.5 * phiMax); // headroom for the fit text
        hPhi[0].SetTitle(PlotTitle(phiFile, "CE proton lab phi, cos fit N_{0} (1 + P #LTA#GT cos#phi) with N_{0} and #LTA#GT fixed, beam neutron CE vs all CE" + what + ";Angle (degrees);Count").c_str());
        canvas.Clear(); canvas.SetLogy(0);
        hPhi[0].Draw("hist"); hPhi[1].Draw("hist same");
        phiFit[0].DrawCopy("same"); phiFit[1].DrawCopy("same");
        TLegend legPhi(0.745, 0.45, 0.99, 0.55); legPhi.SetTextSize(0.016);
        for (int g = 0; g < 2; g++) legPhi.AddEntry(&hPhi[g], Form("%s (%lld)", g ? "beam CE" : "all CE", nPhi[g]), "l");
        legPhi.Draw();
        // the fit results go in the plot area, the legend box is too narrow for them
        TLatex txPhi; txPhi.SetNDC(); txPhi.SetTextSize(0.022);
        for (int g = 0; g < 2; g++){ // two lines per curve, kept inside the plot area (the legend box sits lower in the right margin)
            txPhi.SetTextColor(phiColor[g]);
            txPhi.DrawLatex(0.13, 0.875 - 0.10 * g, Form("%s,   P = %.3f #pm %.3f", phiName[g], pFit[g], pFitErr[g]));
            txPhi.DrawLatex(0.13, 0.845 - 0.10 * g, Form("#LTA#GT = %.4f,   #chi^{2}/ndf = %.1f/%d", sAng[g] / std::max<long long>(1, nPhi[g]), phiFit[g].GetChisquare(), phiFit[g].GetNDF()));        }
        WritePlotNumbers(&hPhi[0], outDir + phiFile, numbersStarted);
        canvas.SaveAs((outDir + phiFile).c_str());
    };
    phiFitPlot("5_phi_cosfit/18_phiCE_cosfit_beam_vs_all.pdf", "", [](Long64_t){return true;});
    // 19_phiCE_cosfit_beam_vs_all_vertex1mm.pdf: plot 18 for the CE whose vertex (fX, fY, the neutron position at the CE) is inside the 1 mm square, a neutron that scattered before the CE has moved off the beam axis
    phiFitPlot("5_phi_cosfit/19_phiCE_cosfit_beam_vs_all_vertex1mm.pdf", ", CE vertex within |x|,|y| #leq 1 mm", [&](Long64_t k){return std::fabs(ceX[k]) <= 1 && std::fabs(ceYv[k]) <= 1;});
    // 16_CE_neutron_vertex_z_by_creator.pdf and 17_..._cut80_by_creator.pdf: z where the CE neutron was made (CEVertexZ, mm), by the process that made it; 17 keeps the CE whose proton is above the 80% KE cut.
    // Needs CEVertexZ in ChargeExchange_values (sim of 28 Oct 2026 or later).
    if (hasVz){
        std::map<std::pair<int, int>, double> vzOf;
        double vLo = 1e30, vHi = -1e30;
        for (Long64_t k = 0; k < nCE; k++){const double vz = ceCol[7][k]; vzOf[{(int)ceEvent[k], (int)ceTrack[k]}] = vz; vLo = std::min(vLo, vz); vHi = std::max(vHi, vz);}
        const double vPad = 0.02 * (vHi - vLo) + 1e-6;
        const std::string vFile[2] = {"4_position/16_CE_neutron_vertex_z_by_creator.pdf", "4_position/17_CE_neutron_vertex_z_cut80_by_creator.pdf"};
        TH1D hV[2][4];
        for (int q = 0; q < 2; q++) for (int g = 0; g < 4; g++) hV[q][g] = TH1D(Form("hV%d_%d", q, g), "", 100, vLo - vPad, vHi + vPad);
        for (const auto& [key, g] : groupOf) hV[0][g].Fill(vzOf.at(key));
        for (Long64_t k = 0; k < nPi; k++){
            const int piEvent = (int)piCol[0][k], piParent = (int)piCol[1][k]; const double piDE = piCol[4][k];
            const auto it = groupOf.find({piEvent, piParent});
            if (it != groupOf.end() && piDE > cutKE) hV[1][it->second].Fill(vzOf.at({piEvent, piParent}));
        }
        for (int q = 0; q < 2; q++){
            double vMax = 0;
            for (int g = 0; g < 4; g++){hV[q][g].SetLineColor(gColor[g]); hV[q][g].SetLineWidth(2); hV[q][g].SetStats(0); vMax = std::max(vMax, hV[q][g].GetMaximum());}
            hV[q][0].SetMinimum(0.5); hV[q][0].SetMaximum(2.0 * vMax);
            const std::string cutNote = q ? std::string(Form(", 80%% KE cut = %.0f MeV", cutKE)) : "";
            hV[q][0].SetTitle(PlotTitle(vFile[q], "z position where the CE neutron was made, by the process that made it" + cutNote + ";z (mm);Count").c_str());
            canvas.Clear(); canvas.SetLogy(1);
            TLegend legV(0.745, 0.65, 0.99, 0.88); legV.SetTextSize(0.016);
            for (int g = 0; g < 4; g++){hV[q][g].Draw(g ? "hist same" : "hist"); legV.AddEntry(&hV[q][g], Form("%s (%.0f)", gName[g], hV[q][g].GetEntries()), "l");}
            legV.Draw();
            WritePlotNumbers(&hV[q][0], outDir + vFile[q], numbersStarted);
            canvas.SaveAs((outDir + vFile[q]).c_str());
        }
        canvas.SetLogy(0);
    }
    // 14_first_tracker_x_by_creator.pdf and 15_..._cut80_by_creator.pdf: x where the CE proton first enters a tracker (ParticleInfo boundaryType 0, the smallest fZ of the proton), by the process that
    // made the CE neutron; 15 keeps the protons above the 80% KE cut on their smeared ZDC energy. Protons are matched to their CE neutron on (fEvent, ParentID) like plot 6.
    std::map<std::pair<int, int>, double> zdcDE; // (fEvent, proton TrackID) -> smeared energy at the ZDC
    for (Long64_t k = 0; k < nPi; k++) zdcDE[{(int)piCol[0][k], (int)piCol[2][k]}] = piCol[4][k];
    struct FirstHit {double z, x, y; int group;};
    std::map<std::pair<int, int>, FirstHit> firstHit;
    TTree* ptPtr = OpenCut("cuts/cut_ParentTracker.root", "ParticleInfo"); // first tracker protons of the CE children
    if (!ptPtr){Error("SecondaryCE", "no cuts/cut_ParentTracker.root (run CutEvents.cpp)"); return;}
    TTree& pt = *ptPtr;
    gROOT->cd();
    const auto ptCol = Columns(&pt, {"fEvent", "ParentID", "TrackID", "fX", "fY", "fZ"}); // read once, multithreaded, every row in entry order
    for (size_t k = 0; k < ptCol[0].size(); k++){
        const int piEvent = (int)ptCol[0][k], piParent = (int)ptCol[1][k], piTrack = (int)ptCol[2][k]; const double piX = ptCol[3][k], piY = ptCol[4][k], piZ = ptCol[5][k];
        const auto it = groupOf.find({piEvent, piParent});
        if (it == groupOf.end()) continue;
        const auto [hit, isNew] = firstHit.try_emplace({piEvent, piTrack}, FirstHit{piZ, piX, piY, it->second});
        if (!isNew && piZ < hit->second.z) hit->second = {piZ, piX, piY, it->second};
    }
    const std::string trkFile[2] = {"4_position/14_first_tracker_x_by_creator.pdf", "4_position/15_first_tracker_x_cut80_by_creator.pdf"};
    for (int q = 0; q < 2; q++){
        TH1D hT[4];
        double tMax = 0;
        for (int g = 0; g < 4; g++) hT[g] = TH1D(Form("hT%d_%d", q, g), "", 100, -5, 5);
        std::vector<double> tx[4];
        double tLo = 1e30, tHi = -1e30;
        for (const auto& [key, h] : firstHit){
            if (q && !(zdcDE.count(key) && zdcDE[key] > cutKE)) continue;
            tx[h.group].push_back(h.x); tLo = std::min(tLo, h.x); tHi = std::max(tHi, h.x);
        }
        const double tPad = 0.02 * (tHi - tLo) + 1e-6;
        for (int g = 0; g < 4; g++){
            hT[g] = TH1D(Form("hT%d_%d", q, g), "", 100, tLo - tPad, tHi + tPad);
            for (double x : tx[g]) hT[g].Fill(x);
            hT[g].SetLineColor(gColor[g]); hT[g].SetLineWidth(2); hT[g].SetStats(0); tMax = std::max(tMax, hT[g].GetMaximum());
        }
        const std::string cutNote = q ? std::string(Form(", 80%% KE cut = %.0f MeV", cutKE)) : "";
        hT[0].SetMinimum(0.5); hT[0].SetMaximum(2.0 * tMax);
        hT[0].SetTitle(PlotTitle(trkFile[q], "x position where the CE proton first enters a tracker, by the process that made the CE neutron" + cutNote + ";x (mm);Count").c_str());
        canvas.Clear(); canvas.SetLogy(1);
        TLegend legT(0.745, 0.65, 0.99, 0.88); legT.SetTextSize(0.016);
        for (int g = 0; g < 4; g++){hT[g].Draw(g ? "hist same" : "hist"); legT.AddEntry(&hT[g], Form("%s (%.0f)", gName[g], hT[g].GetEntries()), "l");}
        legT.Draw();
        WritePlotNumbers(&hT[0], outDir + trkFile[q], numbersStarted);
        canvas.SaveAs((outDir + trkFile[q]).c_str());
    }
    canvas.SetLogy(0);
    // 12_CE_protons_per_event_by_creator.pdf and 13_..._cut80_by_creator.pdf: the same count split by the process that made the CE neutron, each curve over the events with at least one CE of that creator
    std::map<int, int> nEvG[2][4]; // [0] every CE, [1] above the cut
    for (Long64_t k = 0; k < nPi; k++){
        const int piEvent = (int)piCol[0][k], piParent = (int)piCol[1][k]; const double piDE = piCol[4][k];
        const auto it = groupOf.find({piEvent, piParent});
        if (it == groupOf.end()) continue;
        nEvG[0][it->second][piEvent]++;
        if (piDE > cutKE) nEvG[1][it->second][piEvent]++;
    }
    const std::string nGFile[2] = {"1_counts/12_CE_protons_per_event_by_creator.pdf", "1_counts/13_CE_protons_per_event_cut80_by_creator.pdf"};
    for (int q = 0; q < 2; q++){
        TH1D hNG[4];
        double nGMax = 0;
        for (int g = 0; g < 4; g++){
            hNG[g] = TH1D(Form("hNG%d_%d", q, g), "", 10, 0.5, 10.5);
            for (const auto& [ev, n] : nEvG[q][g]) hNG[g].Fill(n);
            hNG[g].SetLineColor(gColor[g]); hNG[g].SetLineWidth(2); hNG[g].SetStats(0); nGMax = std::max(nGMax, hNG[g].GetMaximum());
        }
        hNG[0].SetMinimum(0.5); hNG[0].SetMaximum(2.0 * nGMax);
        const std::string cutNote = q ? std::string(Form(", 80%% KE cut = %.0f MeV", cutKE)) : "";
        hNG[0].SetTitle(PlotTitle(nGFile[q], "Reconstructed CE protons per event by the process that made the CE neutron" + cutNote + ";Reconstructed CE protons of that creator in the event;Events").c_str());
        canvas.Clear(); canvas.SetLogy(1);
        TLegend legNG(0.745, 0.65, 0.99, 0.88); legNG.SetTextSize(0.016);
        for (int g = 0; g < 4; g++){hNG[g].Draw(g ? "hist same" : "hist"); legNG.AddEntry(&hNG[g], Form("%s (%.0f events), mean %.3f", gName[g], hNG[g].GetEntries(), hNG[g].GetMean()), "l");}
        legNG.Draw();
        WritePlotNumbers(&hNG[0], outDir + nGFile[q], numbersStarted);
        canvas.SaveAs((outDir + nGFile[q]).c_str());
    }
    canvas.SetLogy(0);
    printf("P: all %.3f +- %.3f, beam only %.3f, est fix nInel %.3f, est fix n+p inel %.3f\n", pAll, errAll, pPrim, pVal[2], pVal[3]);
    printf("SecondaryCE: %lld beam neutron CE, %lld secondary CE (%.1f%%)\n", nPrim, nSec, fracSec);
}
