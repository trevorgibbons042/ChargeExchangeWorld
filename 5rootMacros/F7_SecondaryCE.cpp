#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"

#include "CutFiles.h"
#include "FitCore.h"
#include "TreeColumns.h"

// CE of the neutrons in the carbon blocks, split into the beam neutron (CECreator "primary", TrackID 1) and the secondary neutrons made by hadronic interactions.
// CEPol is the spin the sim carries on the CE track (track->GetPolarization(), unit vectors, the beam PPValue is the fraction along +y; Geant4 gives secondaries none).
// Needs ChargeExchange_values with CECreator, CEPolX/Y/Z and CEVertexZ (cuts/cut_CE.root), AnglesRaw (cuts/cut_angles_all.root) and cuts/cut_ParentTracker.root, and settings.csv.
// Writes plots/7_secondaryCE/<folder>/<n>_<name>.pdf, by creator group (beam neutron, neutronInelastic, protonInelastic, all other creators) unless noted:
//   1_counts: 2_CE_count_by_creator, 4_CE_count_primary_vs_secondary, 12_CE_protons_per_event_by_creator, 13_CE_protons_per_event_cut80_by_creator
//   2_polarization: 3_CEpolY_distribution (beam neutron vs secondary), 5_P_fit_if_inelastic_fixed
//   3_energy: 8_CE_proton_smeared_KE_cut80_by_creator
//   4_position: 9_CE_x_position_by_creator, 10_CE_x_position_cut80_by_creator, 14_first_tracker_x_by_creator, 15_first_tracker_x_cut80_by_creator,
//               16_CE_neutron_vertex_z_by_creator, 17_CE_neutron_vertex_z_cut80_by_creator
//   5_phi_cosfit: 18_phiRe_cosfit_beam_vs_all, 19_phiRe_cosfit_beam_vs_all_vertex1mm (beam neutron CE vs all CE)
void F7_SecondaryCE(){
    FastJit();
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    gStyle->SetOptStat(0);
    // the carbon-block CE (CarbonCeCut) come from cuts/cut_CE.root (CutFiles.h), so every row passes; the first tracker protons of their children likewise
    // (children of a CE track: only those match a groupOf key below)
    TTree* cePtr = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
    if (!cePtr){Error("SecondaryCE", "no cuts/cut_CE.root (run CutPlots.cpp (events))"); return;}
    TTree& ce = *cePtr;
    gROOT->cd();
    const std::string outDir = "plots/7_secondaryCE/";
    for (const char* sub : {"1_counts", "2_polarization", "3_energy", "4_position", "5_phi_cosfit"}) gSystem->mkdir((outDir + sub).c_str(), true); // one folder per kind of plot

    // carbon-block CE only, like the other CE plots (CarbonCeCut in FitCore.h)
    // ce is read once, multithreaded, into one vector per column (Columns, TreeColumns.h; CECreator is a char leaf, read as a std::string), every row of the cache in entry order
    std::vector<std::string> ceExprs = {"fEvent", "TrackID", "CEPolY", "CEPhiNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ)", "fX", "fY", "CEThetaNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "CEVertexZ"};
    const auto ceCol = Columns(&ce, ceExprs);
    const auto& ceEvent = ceCol[0], &ceTrack = ceCol[1], &cePolY = ceCol[2], &cePhiLab = ceCol[3], &ceX = ceCol[4], &ceYv = ceCol[5], &ceTheta = ceCol[6]; // ceTheta: CETHETA_NEUFRAME (degrees), the A_eff of the phi fits comes from its 0.05 degree bins (AEffOfTheta)
    // the CE recon events in the estimated neutron frame, for 5_phi_cosfit: fEvent -> RECONPHI_ESTNEU, RECONTHETA_ESTNEU (deg; thetaRelN is NAN without a tube, so those events are left out)
    std::map<int, std::array<double, 5>> reconEst; // phi, theta, then the est neutron range: theta min, theta max (deg), phi half width (rad)
    TTree* raw = OpenCut("cuts/cut_angles_all.root", "AnglesRaw"); // also the ZDC KE columns below
    if (!raw){Error("SecondaryCE", "no cuts/cut_angles_all.root (run ReconAngleAll.cpp)"); return;}
    {
        const auto rr = Columns(raw, {"fEvent", "phiRelN * 180 / TMath::Pi()", "thetaRelN * 180 / TMath::Pi()", "thetaRelNMin * 180 / TMath::Pi()", "thetaRelNMax * 180 / TMath::Pi()", "phiRelNHalf"}, "isFromCE == 1 && thetaRelN >= 0");
        for (size_t i = 0; i < rr[0].size(); i++) reconEst[(int)rr[0][i]] = {rr[1][i], rr[2][i], rr[3][i], rr[4][i], rr[5][i]};
        gROOT->cd();
    }
    const double apScaleCE = ReadSettings().at("APScale"), ppValue = ReadSettings().at("PPValue");
    const auto ceCreator = Columns<std::string>(&ce, {"std::string(CECreator.data())"})[0];
    const Long64_t nCE = ceEvent.size();
    const auto creatorGroup = [](const std::string& c){return c == "primary" ? 0 : c == "neutronInelastic" ? 1 : c == "protonInelastic" ? 2 : 3;};
    std::map<std::pair<int, int>, int> groupOf; // (fEvent, CE neutron TrackID) -> 0 beam neutron, 1 neutronInelastic, 2 protonInelastic, 3 every other creator
    double sA[2] = {0, 0}, sCos[2] = {0, 0}; // model A (from CETHETA_NEUFRAME) and cos(CEPHI_NEUFRAME) sums of the beam neutron CE (0) and the secondary CE (1), for the P estimate
    TH1D hYprim("hYprim", "", 42, -1.05, 1.05), hYsec("hYsec", "", 42, -1.05, 1.05); // CEPolY of the beam neutron and of the secondary neutrons
    std::map<std::string, long long> byCreator;
    long long nPrim = 0, nSec = 0;
    std::vector<double> xs[4]; // x of the neutron at the CE by creator group, for the x plots below
    double xLo = 1e30, xHi = -1e30;
    std::map<std::pair<int, int>, double> xOf; // (fEvent, CE neutron TrackID) -> x at the CE
    for (Long64_t k = 0; k < nCE; k++){
        const int fEvent = (int)ceEvent[k], trackID = (int)ceTrack[k];
        const int g = trackID == 1 ? 0 : 1;
        (g == 0 ? hYprim : hYsec).Fill(cePolY[k]);
        sA[g] += AEffOfTheta(ceTheta[k], apScaleCE); sCos[g] += std::cos(cePhiLab[k]);
        if (g == 0) nPrim++; else nSec++;
        byCreator[ceCreator[k]]++;
        const int creator = creatorGroup(ceCreator[k]);
        groupOf[{fEvent, trackID}] = creator;
        xs[creator].push_back(ceX[k]); xLo = std::min(xLo, ceX[k]); xHi = std::max(xHi, ceX[k]);
        xOf[{fEvent, trackID}] = ceX[k];
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
    hYprim.SetLineColor(kBlack); hYsec.SetLineColor(CbCE());
    hYprim.SetLineWidth(2); hYsec.SetLineWidth(2); hYprim.SetStats(0);
    hYprim.SetMinimum(0.5);
    hYprim.SetMaximum(2.0 * std::max(hYprim.GetMaximum(), hYsec.GetMaximum()));
    hYprim.SetTitle(PlotTitle(polYFile, "y component of the polarization of the neutron doing CE in the carbon blocks, no cuts;CEPolY;Count").c_str());
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
    hCreator.SetTitle(PlotTitle(creatorFile, Form("CE in the carbon blocks by the process that made the CE neutron, %.1f%% from secondary neutrons, no cuts;Creator process;Count", fracSec)).c_str());
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
    hSplit.SetTitle(PlotTitle(splitFile, Form("CE in the carbon blocks, beam neutron and all other neutron sources, %.1f%% from the other sources, no cuts;Neutron source;Count", fracSec)).c_str());
    canvas.Clear(); canvas.SetLogy(1);
    hSplit.Draw("hist");
    WritePlotNumbers(&hSplit, outDir + splitFile, numbersStarted);
    canvas.SaveAs((outDir + splitFile).c_str());

    // 5_P_fit_if_inelastic_fixed.pdf: fitted P of the CE protons now, and what it would be if the neutrons made by neutronInelastic (and protonInelastic) were polarized like the beam neutron.
    // P = 2 <cos phi> / <A> (the unweighted cos fit with fixed <A>). The two "fixed" bars are estimates: the CE of those creators
    // count as polarized like the beam-neutron CE, P_est = P(beam neutron CE) * (N_beam + N_fixed) / N_all.
    const std::string pFile = "2_polarization/5_P_fit_if_inelastic_fixed.pdf";
    const double nAll = nPrim + nSec, pPrim = std::clamp(2 * sCos[0] / sA[0], -1.0, 1.0), pAll = std::clamp(2 * (sCos[0] + sCos[1]) / (sA[0] + sA[1]), -1.0, 1.0); // P is a polarization
    const double errPrim = std::sqrt(2.0 / nPrim) / std::fabs(sA[0] / nPrim), errAll = std::sqrt(2.0 / nAll) / std::fabs((sA[0] + sA[1]) / nAll);
    const long long nNeutronInel = byCreator.count("neutronInelastic") ? byCreator["neutronInelastic"] : 0, nProtonInel = byCreator.count("protonInelastic") ? byCreator["protonInelastic"] : 0;
    TH1D hP("hPfix", "", 5, 0, 5);
    const char* pName[5] = {"input PPValue", "now, all CE", "neutronInelastic fixed (est.)", "neutron + protonInelastic fixed (est.)", "beam neutron CE only"};
    const double pVal[5] = {ppValue, pAll, pPrim * (nPrim + nNeutronInel) / nAll, pPrim * (nPrim + nNeutronInel + nProtonInel) / nAll, pPrim};
    const double pErr[5] = {0, errAll, errPrim * (nPrim + nNeutronInel) / nAll, errPrim * (nPrim + nNeutronInel + nProtonInel) / nAll, errPrim};
    for (int i = 0; i < 5; i++){hP.SetBinContent(i + 1, pVal[i]); hP.SetBinError(i + 1, pErr[i]); hP.GetXaxis()->SetBinLabel(i + 1, pName[i]);}
    hP.SetLineColor(kBlack); hP.SetLineWidth(2); hP.SetStats(0);
    hP.SetMinimum(0); hP.SetMaximum(1.0);
    hP.GetXaxis()->SetLabelSize(0.022); hP.GetXaxis()->SetTitleOffset(1.6);
    hP.SetTitle(PlotTitle(pFile, "Fitted polarization P of the CE protons, CEPHI_NEUFRAME (true neutron frame), now and if the neutronInelastic and protonInelastic neutrons were polarized like the beam neutron, no cuts;;P").c_str());
    canvas.Clear(); canvas.SetLogy(0);
    hP.Draw("hist");
    TLine pRef(0, ppValue, 5, ppValue); pRef.SetLineStyle(3); pRef.SetLineWidth(2); pRef.SetLineColor(CbCE()); pRef.Draw(); // input PPValue
    WritePlotNumbers(&hP, outDir + pFile, numbersStarted);
    canvas.SaveAs((outDir + pFile).c_str());
    // ZDC KE of the event of each CE, smeared (zdcKEAllD, /zdc/resStochastic and /zdc/resConstant in the sim, equal to the true KE when the smearing is off), filled for plot 8 and the 80% KE cut of the plots below,
    // one curve per process that made the CE neutron. The KE is matched on fEvent only (the CE neutron of the event, TrackID = ParentID of its proton, is the key of groupOf).
    const char* gName[4] = {"beam neutron", "neutronInelastic", "protonInelastic", "all other creators"};
    const int gColor[4] = {kBlack, CbCE(), CbProton(), CbAll()};
    // line colour and width, no stats box and the log y range (min 0.5, max twice the tallest curve, set on h[0] that is drawn first) of the 4 creator group histograms h[0..3]
    const auto styleGroups = [&](TH1D* h){
        double hMax = 0;
        for (int g = 0; g < 4; g++){h[g].SetLineColor(gColor[g]); h[g].SetLineWidth(2); h[g].SetStats(0); hMax = std::max(hMax, h[g].GetMaximum());}
        h[0].SetMinimum(0.5); h[0].SetMaximum(2.0 * hMax);
    };
    // the 4 creator curves of h[0..3] overlaid (h[0] first), the legend with label(g) of each curve (its lower edge legLo), extra (a line) on top, then the numbers csv and the pdf (the canvas is cleared and its log y set by the caller)
    const auto saveGroups = [&](TH1D* h, const std::string& file, double legLo, auto&& label, TObject* extra = nullptr){
        TLegend leg(0.745, legLo, 0.99, 0.88); leg.SetTextSize(0.016);
        for (int g = 0; g < 4; g++){h[g].Draw(g ? "hist same" : "hist"); leg.AddEntry(&h[g], label(g).c_str(), "l");}
        if (extra) extra->Draw();
        leg.Draw();
        WritePlotNumbers(&h[0], outDir + file, numbersStarted);
        canvas.SaveAs((outDir + file).c_str());
    };
    TH1D hKE[4]; // smeared ZDC KE
    for (int g = 0; g < 4; g++) hKE[g] = TH1D(Form("hKE%d", g), "", 100, 0, KeAxisMax()); // 30000 MeV at 25 GeV/c, FitCore.h
    // The ZDC KE of an event is the sum over all particles, charged and neutral, that hit the ZDC (zdcKEAll true, zdcKEAllD smeared, -1 / NAN without the sim column, 0 without a hit): the columns of the recon events in AnglesRaw.
    // an event column, the same on every row of the event: read once per event from the Events tree (nRecoEv > 0 are the events that have rows in AnglesRaw), in any order
    if (!EventsTree()){Error("SecondaryCE", "no cuts/cut_events.root (run ReconAngleAll.cpp)"); return;}
    const auto keCol = ColumnsUnordered(EventsTree(), {"fEvent", "zdcKEAllD"}, "nRecoEv > 0");
    std::map<int, double> evKE; // fEvent -> smeared ZDC KE, only the events that have a ZDC KE
    for (size_t k = 0; k < keCol[0].size(); k++) if (std::isfinite(keCol[1][k]) && keCol[1][k] > 0) evKE[(int)keCol[0][k]] = keCol[1][k];
    const double cutKE = 0.80 * kBeamKE; // MeV, 80% of the beam KE, the cut on the smeared ZDC KE zdcKEAllD of the event
    long long nAbove[4] = {0, 0, 0, 0}; // CE of each creator whose event is above the cut
    // f(group, event, CE neutron TrackID, smeared KE) for every CE with a creator group whose event has a ZDC KE: the KE is the one of the event, shared by its CE
    const auto forPi = [&](auto&& f){
        for (const auto& [key, g] : groupOf){
            const auto it = evKE.find(key.first);
            if (it != evKE.end()) f(g, key.first, key.second, it->second);
        }
    };
    forPi([&](int g, int, int, double piDE){hKE[g].Fill(piDE); if (piDE > cutKE) nAbove[g]++;});
    const std::string kZdcNote = "event has ZDC KE"; // the cut of every plot built from evKE, named in its title
    styleGroups(hKE);
    // 8_CE_proton_smeared_KE_cut80_by_creator.pdf: the smeared energy plot with the 80% KE cut (dashed), how many CE protons of each creator stay above it
    const std::string cutFile = "3_energy/8_CE_proton_smeared_KE_cut80_by_creator.pdf";
    canvas.Clear(); canvas.SetLogy(1);
    hKE[0].SetTitle(PlotTitle(cutFile, Form("Smeared ZDC KE sum of the event, one entry per CE, by the process that made the CE neutron, %s, 80%% KE cut = %.0f MeV;Smeared ZDC KE zdcKEAllD (MeV);Count", kZdcNote.c_str(), cutKE)).c_str());
    TLine cutLine(cutKE, 0.5, cutKE, hKE[0].GetMaximum()); cutLine.SetLineStyle(2); cutLine.SetLineWidth(2);
    saveGroups(hKE, cutFile, 0.55, [&](int g){return std::string(Form("%s %lld/%.0f (%.0f%%)", gName[g], nAbove[g], hKE[g].GetEntries(), 100.0 * nAbove[g] / std::max(1.0, hKE[g].GetEntries())));}, &cutLine);
    canvas.SetLogy(0);
    // 9_CE_x_position_by_creator.pdf: x of the neutron at the CE (fX, mm), one curve per process that made the CE neutron, same groups as the KE plots
    const std::string xFile = "4_position/9_CE_x_position_by_creator.pdf";
    const double xPad = 0.02 * (xHi - xLo) + 1e-6;
    TH1D hX[4];
    for (int g = 0; g < 4; g++){
        hX[g] = TH1D(Form("hX%d", g), "", 100, xLo - xPad, xHi + xPad);
        for (double x : xs[g]) hX[g].Fill(x);
    }
    styleGroups(hX);
    hX[0].SetTitle(PlotTitle(xFile, "x position of the neutron at the CE in the carbon blocks, by the process that made the CE neutron, no cuts;x (mm);Count").c_str());
    canvas.Clear(); canvas.SetLogy(1);
    saveGroups(hX, xFile, 0.65, [&](int g){return std::string(Form("%s (%.0f)", gName[g], hX[g].GetEntries()));});
    // 10_CE_x_position_cut80_by_creator.pdf: the same x plot for the CE whose event is above the 80% KE cut on the smeared ZDC KE (plot 8)
    const std::string xCutFile = "4_position/10_CE_x_position_cut80_by_creator.pdf";
    TH1D hXc[4];
    for (int g = 0; g < 4; g++) hXc[g] = TH1D(Form("hXc%d", g), "", 100, xLo - xPad, xHi + xPad);
    forPi([&](int g, int piEvent, int piParent, double piDE){if (piDE > cutKE) hXc[g].Fill(xOf.at({piEvent, piParent}));});
    styleGroups(hXc);
    hXc[0].SetTitle(PlotTitle(xCutFile, Form("x position of the neutron at the CE, by the process that made the CE neutron, %s, 80%% KE cut = %.0f MeV;x (mm);Count", kZdcNote.c_str(), cutKE)).c_str());
    canvas.Clear();
    saveGroups(hXc, xCutFile, 0.65, [&](int g){return std::string(Form("%s %.0f/%.0f (%.0f%%)", gName[g], hXc[g].GetEntries(), hX[g].GetEntries(), 100.0 * hXc[g].GetEntries() / std::max(1.0, hX[g].GetEntries())));});
    // 18_phiRe_cosfit_beam_vs_all.pdf: the cos phi fit (N(phi) = N0 (1 + P <A> cos phi), N0 = entries / bins and <A> fixed, only P fitted, 50 bins of the CE
    // proton RECONPHI_ESTNEU, A_eff from RECONTHETA_ESTNEU) for the CE of the beam neutron only (TrackID 1, red) and for every carbon-block CE (black)
    // keep(k) picks the CE rows (it reads the vectors of the ce entry k), what is added to the title; 19_ is the same plot for the CE whose vertex (fX, fY) is within 1 mm of the beam axis in x and y
    const auto phiFitPlot = [&](const std::string& phiFile, const std::string& what, const std::function<bool(Long64_t)>& keep){
        TH1D hPhi[2] = {TH1D("hPhiAll", "", 50, 0, 360), TH1D("hPhiBeam", "", 50, 0, 360)};
        PhiRange phiRng[2]; AMean sAng[2]; long long nPhi[2] = {0, 0}; // sAng: A_eff of each set (AEffOfTheta at the own theta of each proton), its mean and error are the fixed A of the fit
        std::vector<double> phiV[2]; // phi of every CE proton of each set: the unbinned fit below uses them, the histogram only shows them
        for (Long64_t k = 0; k < nCE; k++){
            if (!keep(k)) continue;
            const auto est = reconEst.find((int)ceEvent[k]); // the recon fit of the event of this CE, in the estimated neutron frame
            if (est == reconEst.end()) continue;
            const double phiDeg = est->second[0], aEff = AEffOfTheta(est->second[1], apScaleCE), aLo = AEffOfTheta(est->second[2], apScaleCE), aHi = AEffOfTheta(est->second[3], apScaleCE); const int trackID = (int)ceTrack[k];
            for (int g = 0; g < 2; g++){
                if (g == 1 && trackID != 1) continue;
                hPhi[g].Fill(phiDeg); phiV[g].push_back(phiDeg); nPhi[g]++;
                if (std::isfinite(aLo) && std::isfinite(aHi)){sAng[g].Add(aEff, aLo, aHi); if (std::isfinite(est->second[4])) phiRng[g].Add(est->second[4]);}
                else sAng[g].Add(aEff);
            }
        }
        const char* phiName[2] = {"all CE", "beam neutron CE (TrackID 1)"};
        const int phiColor[2] = {kBlack, CbCE()};
        TF1 phiFit[2];
        PFit ufBox[2];
        double phiMax = 0;
        for (int g = 0; g < 2; g++){
            ufBox[g] = UnbinnedP(phiV[g], sAng[g], 0, 360, nullptr, true, &phiRng[g]); // unbinned fit of every proton's phi (FitCore.h), A fixed to the mean
            phiFit[g] = *AddUnbinnedFit(&hPhi[g], Form("fPhiCos%d", g), ufBox[g], sAng[g].Mean(), sAng[g].Err(), false); // in the function list for fit_values.csv, a copy here for the dashed curve
            phiFit[g].ResetBit(TF1::kNotDraw);
            hPhi[g].SetLineColor(phiColor[g]); hPhi[g].SetLineWidth(2); hPhi[g].SetStats(0);
            phiFit[g].SetLineColor(phiColor[g]); phiFit[g].SetLineStyle(2); phiFit[g].SetLineWidth(2);
            phiMax = std::max(phiMax, hPhi[g].GetMaximum());
        }
        hPhi[0].SetMinimum(0); hPhi[0].SetMaximum(1.5 * phiMax); // headroom for the fit text
        hPhi[0].SetTitle(PlotTitle(phiFile, "CE recon protons, RECONPHI_ESTNEU (estimated neutron frame), unbinned cos fit N(#phi) #propto 1 + P A_{eff} cos#phi with A_{eff} fixed, beam neutron CE vs all CE" + what + ";Angle (degrees);Count").c_str());
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
            txPhi.DrawLatex(0.13, 0.875 - 0.10 * g, Form("%s,   P = %.3f #pm %.3f", phiName[g], ufBox[g].p, ufBox[g].err));
            txPhi.DrawLatex(0.13, 0.845 - 0.10 * g, Form("%s = %.4f #pm %.4f,   %lld CE protons", AEffEq(), sAng[g].Mean(), sAng[g].Err(), nPhi[g]));        }
        WritePlotNumbers(&hPhi[0], outDir + phiFile, numbersStarted);
        canvas.SaveAs((outDir + phiFile).c_str());
    };
    phiFitPlot("5_phi_cosfit/18_phiRe_cosfit_beam_vs_all.pdf", ", no cuts", [](Long64_t){return true;});
    // 19_phiRe_cosfit_beam_vs_all_vertex1mm.pdf: plot 18 for the CE whose vertex (fX, fY, the neutron position at the CE) is within 1 mm of the beam axis in x and y, a neutron that scattered before the CE has moved off the beam axis
    phiFitPlot("5_phi_cosfit/19_phiRe_cosfit_beam_vs_all_vertex1mm.pdf", ", CE vertex within |x|,|y| #leq 1 mm", [&](Long64_t k){return std::fabs(ceX[k]) <= 1 && std::fabs(ceYv[k]) <= 1;});
    // 16_CE_neutron_vertex_z_by_creator.pdf and 17_..._cut80_by_creator.pdf: z where the CE neutron was made (CEVertexZ, mm), by the process that made it; 17 keeps the CE whose proton is above the 80% KE cut.
    std::map<std::pair<int, int>, double> vzOf;
    double vLo = 1e30, vHi = -1e30;
    for (Long64_t k = 0; k < nCE; k++){const double vz = ceCol[7][k]; vzOf[{(int)ceEvent[k], (int)ceTrack[k]}] = vz; vLo = std::min(vLo, vz); vHi = std::max(vHi, vz);}
    const double vPad = 0.02 * (vHi - vLo) + 1e-6;
    const std::string vFile[2] = {"4_position/16_CE_neutron_vertex_z_by_creator.pdf", "4_position/17_CE_neutron_vertex_z_cut80_by_creator.pdf"};
    TH1D hV[2][4];
    for (int q = 0; q < 2; q++) for (int g = 0; g < 4; g++) hV[q][g] = TH1D(Form("hV%d_%d", q, g), "", 100, vLo - vPad, vHi + vPad);
    for (const auto& [key, g] : groupOf) hV[0][g].Fill(vzOf.at(key));
    forPi([&](int g, int piEvent, int piParent, double piDE){if (piDE > cutKE) hV[1][g].Fill(vzOf.at({piEvent, piParent}));});
    for (int q = 0; q < 2; q++){
        styleGroups(hV[q]);
        const std::string cutNote = q ? std::string(Form(", %s, 80%% KE cut = %.0f MeV", kZdcNote.c_str(), cutKE)) : ", no cuts";
        hV[q][0].SetTitle(PlotTitle(vFile[q], "z position where the CE neutron was made, by the process that made it" + cutNote + ";z (mm);Count").c_str());
        canvas.Clear(); canvas.SetLogy(1);
        saveGroups(hV[q], vFile[q], 0.65, [&](int g){return std::string(Form("%s (%.0f)", gName[g], hV[q][g].GetEntries()));});
    }
    canvas.SetLogy(0);
    // 14_first_tracker_x_by_creator.pdf and 15_..._cut80_by_creator.pdf: x where the CE proton first enters a tracker (ParticleInfo boundaryType 0, the smallest fZ of the proton), by the process that
    // made the CE neutron; 15 keeps the protons above the 80% KE cut on their smeared ZDC energy. Protons are matched to their CE neutron on (fEvent, ParentID).
    struct FirstHit {double z, x, y; int group;};
    std::map<std::pair<int, int>, FirstHit> firstHit;
    TTree* ptPtr = OpenCut("cuts/cut_ParentTracker.root", "ParticleInfo"); // first tracker protons of the CE children
    if (!ptPtr){Error("SecondaryCE", "no cuts/cut_ParentTracker.root (run CutPlots.cpp (events))"); return;}
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
        std::vector<double> tx[4];
        double tLo = 1e30, tHi = -1e30;
        for (const auto& [key, h] : firstHit){
            if (q && !(evKE.count(key.first) && evKE.at(key.first) > cutKE)) continue;
            tx[h.group].push_back(h.x); tLo = std::min(tLo, h.x); tHi = std::max(tHi, h.x);
        }
        const double tPad = 0.02 * (tHi - tLo) + 1e-6;
        for (int g = 0; g < 4; g++){
            hT[g] = TH1D(Form("hT%d_%d", q, g), "", 100, tLo - tPad, tHi + tPad);
            for (double x : tx[g]) hT[g].Fill(x);
        }
        styleGroups(hT);
        const std::string cutNote = q ? std::string(Form(", hit tracker and %s, 80%% KE cut = %.0f MeV", kZdcNote.c_str(), cutKE)) : ", hit tracker (cut_ParentTracker)";
        hT[0].SetTitle(PlotTitle(trkFile[q], "x position where the CE proton first enters a tracker, by the process that made the CE neutron" + cutNote + ";x (mm);Count").c_str());
        canvas.Clear(); canvas.SetLogy(1);
        saveGroups(hT, trkFile[q], 0.65, [&](int g){return std::string(Form("%s (%.0f)", gName[g], hT[g].GetEntries()));});
    }
    canvas.SetLogy(0);
    // 12_CE_protons_per_event_by_creator.pdf and 13_..._cut80_by_creator.pdf: the same count split by the process that made the CE neutron, each curve over the events with at least one CE of that creator
    std::map<int, int> nEvG[2][4]; // [0] every CE, [1] above the cut
    forPi([&](int g, int piEvent, int, double piDE){nEvG[0][g][piEvent]++; if (piDE > cutKE) nEvG[1][g][piEvent]++;});
    const std::string nGFile[2] = {"1_counts/12_CE_protons_per_event_by_creator.pdf", "1_counts/13_CE_protons_per_event_cut80_by_creator.pdf"};
    for (int q = 0; q < 2; q++){
        TH1D hNG[4];
        for (int g = 0; g < 4; g++){
            hNG[g] = TH1D(Form("hNG%d_%d", q, g), "", 10, 0.5, 10.5);
            for (const auto& [ev, n] : nEvG[q][g]) hNG[g].Fill(n);
        }
        styleGroups(hNG);
        const std::string cutNote = q ? std::string(Form(", %s, 80%% KE cut = %.0f MeV", kZdcNote.c_str(), cutKE)) : ", " + kZdcNote + ", no KE cut";
        hNG[0].SetTitle(PlotTitle(nGFile[q], "CE per event by the process that made the CE neutron" + cutNote + ";CE of that creator in the event;Events").c_str());
        canvas.Clear(); canvas.SetLogy(1);
        saveGroups(hNG, nGFile[q], 0.65, [&](int g){return std::string(Form("%s (%.0f events), mean %.3f", gName[g], hNG[g].GetEntries(), hNG[g].GetMean()));});
    }
    canvas.SetLogy(0);
    printf("P: all %.3f +- %.3f, beam only %.3f, est fix nInel %.3f, est fix n+p inel %.3f\n", pAll, errAll, pPrim, pVal[2], pVal[3]);
    printf("SecondaryCE: %lld beam neutron CE, %lld secondary CE (%.1f%%)\n", nPrim, nSec, fracSec);
}
