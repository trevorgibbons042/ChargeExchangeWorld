#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

#include "FitCore.h"
#include "CutOpen.h"
#include "TreeColumns.h"

// Response matrices of the recon block angles (thetaReD, phiReD, CalcValsReD in cut_APandT.root) against the truth of the CE proton, for the random beam angle of the sim (/randGun/randPAngle: theta uniform
// 0 to 5 mrad = 0.29 deg, phi uniform). Two truths per CE proton: the lab truth (thetaCE, phiCE in CE: angles about the z axis, the frame the recon block angles are measured in) and the truth relative to the
// incoming neutron (the lab proton direction rotated so the neutron direction CENDirX/Y/Z of ChargeExchange_values, the last row of the track as in CalcAPandT.cpp, points along z).
// The lab truth shows the detector smearing alone, the neutron truth the detector smearing and the beam angle. Only the MC knows the neutron direction: the matrix carries the beam
// angle spread, real data needs only the matrix.
// The detector matrices R_det_* (block and _recon) are the detector smearing as its own matrix: the recon angle against the lab truth (thetaCE, phiCE) only, with their truth T_det_*. UnfoldPhi.cpp unfolds with R_det_*, then with R_beam_thetaphi.
// Two steps that multiply, in the true neutron frame (CE) and the estimated neutron frame (recon): physics (R_*_recon: CE true neutron frame to recon estimated neutron frame, truth T_*) and block (R_block_*: recon to recon block, both estimated neutron frame, truth T_rec_*, the flat one is M_thetaphi_recon). The product of the two flattened probability matrices is R_detector_physics_thetaphi (CE true neutron frame to recon block estimated neutron frame), checked against the direct R_thetaphi. The lab frame matrices (R_beam_*, R_det_*) stay in the root file for UnfoldPhi.cpp.
// Writes cuts/smear_matrix.root (R_theta, R_phi_all, R_phi_<k> per true theta slice, R_thetaphi flattened, their truth and recon histograms, so the matrices can go into RooUnfold) and
// plots/8_smear/1_matrix holds only the 6 pages of the M_recB_Fswitch matrices (recon block, true neutron frame vs estimated neutron frame: matrices and counts; UnfoldPhi.cpp adds 2_phi_cosfit and 3_P_vs_iteration); SavePlot2D skips every other path and the residual block is off (kResidualTypes), the matrices stay in the root file. What the other plots were: 1_Fswitch (M_CE_Fswitch, M_rec_Fswitch, M_recB_Fswitch: CE, recon and recon block angles, neutron frame vs lab frame; M_CE_lab_to_est, M_CE_neu_to_est, M_rec_lab_to_est, M_recB_lab_to_est: CE, recon and recon block lab angles, and the CE true neutron frame, in the estimated neutron frame), 2_M_physics (CE true neutron frame vs recon estimated neutron frame), 3_M_detector (recon vs recon block, estimated neutron frame), 4_chain (their product vs the direct matrix), 5_direct (CE true neutron frame or CE lab vs every recon type), 6_residuals (residuals, dilution, smear_dilution.csv). The SiCE protons are in the root file only, not plotted.
// SiCE_theta, SiCE_thetaphi (the bins of M_thetaphi) and SiCE_theta_phi are the SiCE protons: the protons of AnglesRaw (cut_angles_all.root) that are not children of a carbon CE track (isFromCE == 0; silicon CE in a noPhys run), recon block
// angles only, no truth. For RooUnfoldResponse::Fake: the measured sample of a cut holds them next to the CE protons.
// Rows: every CE proton of CalcValsReD with a neutron direction. Cut: true (neutron frame) and recon block theta both 0 to kThetaMax for the matrices and residuals, recon block theta
// 0 to kThetaMax for the dilution; no other cut.

static constexpr double kStep = 0.1;     // deg, theta bin of R_theta and slice of the dilution
static constexpr int kPhiFine = 36;      // phi bins of R_phi_all and R_phi_<k>
static constexpr int kPhiFlat = 50;      // phi bins of R_thetaphi, same as the PhiCosFit plots of PlotAngle.cpp
static constexpr int kSlices = 10;      // true theta slices of R_phi_<k> and theta bins of R_thetaphi

static double ThetaMaxDeg(){return 5.0 * ThetaScale();} // same range as PlotPolarization.cpp
// seconds since the last call (the first call starts the clock): the ">>> SmearMatrix ... took" lines of the log
static double Secs(){static auto last = std::chrono::steady_clock::now(); const auto now = std::chrono::steady_clock::now(); const double d = std::chrono::duration<double>(now - last).count(); last = now; return d;}
static double Wrap360(double a){a = std::fmod(a, 360.0); return a < 0 ? a + 360.0 : a;}
static double Wrap180(double a){return std::remainder(a, 360.0);}
// The flattened (theta slice, phi bin) index of R_thetaphi for theta in 0 to tMax and phi in 0 to 360 deg.
static int FlatBin(double theta, double phi, double tMax){
    return std::min(kSlices - 1, (int)(theta / (tMax / kSlices))) * kPhiFlat + std::min(kPhiFlat - 1, (int)(phi / (360.0 / kPhiFlat)));
}
struct NDir {double x, y, z, thLab, phLab;}; // neutron direction, and the lab theta and phi (deg) of the CE proton (LastTheta, LastPhiLab)
struct LabRecon {double thB, phB, th, ph, ex, ey, ez;}; // the recon block and recon angles in the neutron estimate frame (deg), outDTheta, outDPhi, outTheta, outPhi of AnglesRaw; ex, ey, ez: the estimated neutron direction, the tube centre at the block middle (mean of tubeX0/tubeX1, tubeY0/tubeY1, tubeZ0/tubeZ1), NaN without a tube
static ULong64_t EvKey(int ev, int tr){return ((ULong64_t)(unsigned)ev << 32) | (unsigned)tr;}

// theta and phi (deg) of the lab direction (thetaDeg, phiDeg) in the frame where the neutron direction d points along z (minimal rotation, axis d x z).
static void RelToNeutron(double thetaDeg, double phiDeg, const NDir& d, double& theta, double& phi){
    TVector3 p;
    p.SetMagThetaPhi(1.0, thetaDeg * TMath::DegToRad(), phiDeg * TMath::DegToRad());
    TVector3 n(d.x, d.y, d.z);
    n = n.Unit();
    const TVector3 axis = n.Cross(TVector3(0, 0, 1));
    if (axis.Mag() > 1e-12) p.Rotate(std::acos(std::max(-1.0, std::min(1.0, n.Z()))), axis);
    theta = p.Theta() * TMath::RadToDeg();
    phi = Wrap360(p.Phi() * TMath::RadToDeg());
}

// Divides each true column (x bin) of R by the truth histogram of that bin, so a column sums to the efficiency of that bin and R reads as P(recon | true).
static void NormColumns(TH2D* R, const TH1D* truth){
    for (int i = 1; i <= R->GetNbinsX(); i++){
        const double t = truth->GetBinContent(i);
        for (int j = 1; j <= R->GetNbinsY(); j++){
            R->SetBinContent(i, j, t > 0 ? R->GetBinContent(i, j) / t : 0.0);
            R->SetBinError(i, j, 0.0);
        }
    }
}

// Same page as the PlotAngle.cpp and PlotPolarization.cpp plots: 1500 x 1200 canvas, the plot left of a wide right margin that holds the legend or palette.
static void PageMargins(TCanvas& c){c.SetRightMargin(0.30); c.SetLeftMargin(0.14); c.SetBottomMargin(0.13);}

// axMax > 0: both axes shown from 0 to axMax only (the theta axes: ThetaAxisMax, like the other theta plots)
// only the recon block, true neutron frame vs estimated neutron frame, pages are plotted (the rest of 8_smear was dropped); the matrices stay in cuts/smear_matrix.root
static constexpr const char* kPlotDir = "plots/8_smear/1_Fswitch/3_M_recB_Fswitch/";
static constexpr const char* kPlotDir2 = "plots/8_smear/1_Fswitch/8_M_recB_est_to_lab/"; // recon block, estimated neutron frame (x) to recon block lab (y)
static constexpr int kResidualTypes = 0; // 6_residuals (residuals, dilution, csv) dropped; 4 brings them back
static void SavePlot2D(TH2D* h, const std::string& pathIn, const std::string& title, double axMax = 0, const char* zTitle = "P(recon | true)"){
    std::string path; // the matrix pages sit in plots/8_smear/1_matrix (true neutron frame to estimated neutron frame) and 2_matrix2 (estimated neutron frame to lab)
    if (pathIn.rfind(kPlotDir, 0) == 0) path = std::string("plots/8_smear/1_matrix/") + pathIn.substr(std::string(kPlotDir).size());
    else if (pathIn.rfind(kPlotDir2, 0) == 0) path = std::string("plots/8_smear/2_matrix2/") + pathIn.substr(std::string(kPlotDir2).size());
    else if (pathIn.rfind("plots/8_smear/3_matrix3/", 0) == 0) path = pathIn; // true neutron theta vs estimated neutron theta, lab frame
    else return;
    TCanvas c("cSmear2D", "Plots", 1500, 1200);
    PageMargins(c);
    h->SetTitle(PlotTitle(path, title).c_str());
    h->SetStats(0);
    h->GetZaxis()->SetTitle(zTitle);
    h->GetZaxis()->SetTitleOffset(1.6); // the palette labels of a small maximum (0.25, 0.3) are wider than the default gap and run into the title
    if (axMax > 0){h->GetXaxis()->SetRangeUser(0, axMax); h->GetYaxis()->SetRangeUser(0, axMax);}
    const bool diverging = std::string(zTitle) == "chain - direct"; // difference pages: white at 0, the viridis end colours on each side
    if (diverging){
        Double_t stops[3] = {0.0, 0.5, 1.0}, red[3] = {0.267, 1.0, 0.993}, green[3] = {0.005, 1.0, 0.906}, blue[3] = {0.329, 1.0, 0.144};
        TColor::CreateGradientColorTable(3, stops, red, green, blue, 255);
        gStyle->SetNumberContours(255);
    }
    if (!diverging){ // log colour scale from the smallest filled bin, empty bins stay white
        double lo = 0;
        for (int i = 1; i <= h->GetNbinsX(); i++) for (int j = 1; j <= h->GetNbinsY(); j++){const double v = h->GetBinContent(i, j); if (v > 0 && (lo == 0 || v < lo)) lo = v;}
        if (lo > 0) h->SetMinimum(lo);
        c.SetLogz();
    }
    h->Draw("COLZ");
    c.SaveAs(path.c_str());
    if (diverging){gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255);}
}

// The matrix as proton counts, beside its probability plot. Drawn from a clone: SavePlot2D narrows the axis range of the histogram it draws, and the file keeps the full one.
static void SaveCounts(const TH2D& h, const std::string& path, const std::string& title, double axMax = 0){
    TH2D* c = (TH2D*)h.Clone(Form("%s_plot", h.GetName()));
    SavePlot2D(c, path, title, axMax, "Protons");
    delete c;
}

// Product of two probability matrices of the same square binning: c(i,l) = sum_j a(i,j) b(j,l), a clone of a with that content.
static TH2D* MatMul(const TH2D& a, const TH2D& b, const char* name){
    TH2D* c = (TH2D*)a.Clone(name);
    c->Reset();
    for (int i = 1; i <= a.GetNbinsX(); i++) for (int l = 1; l <= b.GetNbinsY(); l++){
        double v = 0;
        for (int j = 1; j <= a.GetNbinsY(); j++) v += a.GetBinContent(i, j) * b.GetBinContent(j, l);
        c->SetBinContent(i, l, v);
    }
    return c;
}

static void SaveResidual(TH1D* hNeutron, TH1D* hLab, const char* angle, const std::string& path, const std::string& title){ // angle: "theta" or "phi"
    TCanvas c("cSmearRes", "Plots", 1500, 1200);
    PageMargins(c);
    hNeutron->SetTitle(PlotTitle(path, title).c_str());
    hNeutron->SetStats(0);
    hNeutron->SetLineColor(kRed + 1); hNeutron->SetLineWidth(2);
    hLab->SetLineColor(kBlue + 1); hLab->SetLineWidth(2);
    hNeutron->SetMaximum(1.15 * std::max(hNeutron->GetMaximum(), hLab->GetMaximum()));
    hNeutron->Draw("HIST");
    hLab->Draw("HIST SAME");
    TLegend leg(0.715, 0.68, 0.99, 0.88); // in the right margin, where the stats box of the other histograms sits
    leg.SetTextSize(0.017); // fixed, so the long entries do not shrink to fit the box
    leg.AddEntry(hNeutron, Form("#splitline{#splitline{recon block %s (neutron estimate frame)}{- CE %s (real neutron frame)}}{RMS %.3f}", angle, angle, hNeutron->GetRMS()), "l");
    leg.AddEntry(hLab, Form("#splitline{#splitline{recon block %s (neutron estimate frame)}{- CE %s (neutron estimate frame)}}{RMS %.3f}", angle, angle, hLab->GetRMS()), "l");
    leg.Draw();
    c.SaveAs(path.c_str());
}

void SmearMatrix(){
    Secs();
    ROOT::EnableImplicitMT(NThreads()); // Columns / ColumnsUnordered (TreeColumns.h) read on all the threads only after this (NThreads: ANA_NTHREADS, FitCore.h); without it the two AnglesRaw reads ran on one core, 74 s and 83 s on LCRC
    FastJit();
    TH1::AddDirectory(kFALSE); // the histograms are written to the file by hand at the end, so closing it does not delete the stack ones
    const double tMax = ThetaMaxDeg();
    const int nTheta = (int)std::lround(tMax / kStep);
    double maxTilt = 0; // deg, the largest neutron tilt from z over the CE rows (0: no random beam angle)

    TTree* ce = OpenCut("cuts/cut_APandT.root", "CE");
    TTree* red = OpenCut("cuts/cut_APandT.root", "CalcValsReD");
    TTree* rec = OpenCut("cuts/cut_APandT.root", "CalcValsRe"); // recon theta and phi (the track fit before the block positions), next to the recon block ones
    TTree* cev = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
    if (!ce || !red || !rec || !cev){printf("SmearMatrix: cut_APandT.root (CE, CalcValsReD, CalcValsRe) or cut_CE.root (ChargeExchange_values) missing\n"); return;}
    if (!cev->GetBranch("CENDirX")){printf("SmearMatrix: ChargeExchange_values has no CENDirX (sim older than the CENDir columns)\n"); return;}
    if (ce->GetEntries() != red->GetEntries() || ce->GetEntries() != rec->GetEntries()){printf("SmearMatrix: CE, CalcValsReD and CalcValsRe of cut_APandT.root differ in length\n"); return;}

    // neutron direction of the last CE row of every (event, track), file order, like the CE row CalcAPandT.cpp takes
    std::unordered_map<ULong64_t, NDir> dirs;
    dirs.reserve((size_t)cev->GetEntries());
    {
        // read on all threads (Columns, TreeColumns.h), rows in entry order so the last row of a (event, track) wins as before
        const auto c = Columns(cev, {"fEvent", "TrackID", "CENDirX", "CENDirY", "CENDirZ", "LastTheta", "LastPhiLab"});
        for (size_t i = 0; i < c[0].size(); i++){
            dirs[EvKey((int)c[0][i], (int)c[1][i])] = {c[2][i], c[3][i], c[4][i], c[5][i] * TMath::RadToDeg(), Wrap360(c[6][i] * TMath::RadToDeg())};
            maxTilt = std::max(maxTilt, std::acos(std::min(1.0, c[4][i])) * TMath::RadToDeg());
        }
        printf(">>> SmearMatrix ChargeExchange_values read took %.0f s\n", Secs()); fflush(stdout);
    }

    gSystem->mkdir("plots/8_smear/1_matrix", kTRUE); gSystem->mkdir("plots/8_smear/2_matrix2", kTRUE); gSystem->mkdir("plots/8_smear/3_matrix3", kTRUE);
    TFile out("cuts/smear_matrix.root", "RECREATE");
    const int nFlatT = kSlices * kPhiFlat;
    TH2D hRTheta("R_theta", "", nTheta, 0, tMax, nTheta, 0, tMax);
    TH1D hTTheta("T_theta", "", nTheta, 0, tMax), hMTheta("M_theta", "", nTheta, 0, tMax); // true theta of every row in range; recon theta of every row in range
    TH2D hRPhi("R_phi_all", "", kPhiFine, 0, 360, kPhiFine, 0, 360);
    TH1D hTPhi("T_phi_all", "", kPhiFine, 0, 360);
    // beam angle only, no detector: the true theta and phi of the proton about the neutron direction (x) against about the lab z axis (y, thetaCE and phiCE)
    TH2D hRBeamTheta("R_beam_theta", "", nTheta, 0, tMax, nTheta, 0, tMax), hRBeamPhi("R_beam_phi_all", "", kPhiFine, 0, 360, kPhiFine, 0, 360);
    // for UnfoldPhi.cpp: the measured histogram of the recon (not block) angles, and the mean A_eff of a true theta slice (AEffOfTheta of the neutron theta, as the P fits of PlotAngle.cpp take it)
    TH1D hMFlatRe("M_thetaphi_recon", "", nFlatT, 0, nFlatT);
    TProfile hAPTheta("AP_theta", "", kSlices, 0, tMax);
    const double apScale = ReadSettings().at("APScale");
    TH2D hRBeamFlat("R_beam_thetaphi", "", nFlatT, 0, nFlatT, nFlatT, 0, nFlatT); // the same two angles on the flattened (theta slice, phi bin) numbers of R_thetaphi
    TH2D hRFlat("R_thetaphi", "", nFlatT, 0, nFlatT, nFlatT, 0, nFlatT);
    // CE neutron frame to recon (neutron estimate frame): the same true bins (T_theta, T_phi_all, T_thetaphi) against the recon theta and phi
    TH2D hRThetaRe("R_theta_recon", "", nTheta, 0, tMax, nTheta, 0, tMax), hRPhiRe("R_phi_all_recon", "", kPhiFine, 0, 360, kPhiFine, 0, 360), hRFlatRe("R_thetaphi_recon", "", nFlatT, 0, nFlatT, nFlatT, 0, nFlatT);
    // detector only, no beam angle: the recon (block) angle against the true lab angle of the CE proton (thetaCE, phiCE); truth T_det_*, recon block R_det_*, recon R_det_*_recon
    TH1D hTDetTheta("T_det_theta", "", nTheta, 0, tMax), hTDetPhi("T_det_phi_all", "", kPhiFine, 0, 360), hTDetFlat("T_det_thetaphi", "", nFlatT, 0, nFlatT);
    TH2D hRDetTheta("R_det_theta", "", nTheta, 0, tMax, nTheta, 0, tMax), hRDetPhi("R_det_phi_all", "", kPhiFine, 0, 360, kPhiFine, 0, 360), hRDetFlat("R_det_thetaphi", "", nFlatT, 0, nFlatT, nFlatT, 0, nFlatT);
    // block step: the recon lab angle against the recon block lab angle of the same proton (truth T_rec_*, the recon lab angles of every row in range)
    TH1D hTRecTheta("T_rec_theta", "", nTheta, 0, tMax), hTRecPhi("T_rec_phi_all", "", kPhiFine, 0, 360);
    TH2D hRBlockTheta("R_block_theta", "", nTheta, 0, tMax, nTheta, 0, tMax), hRBlockPhi("R_block_phi_all", "", kPhiFine, 0, 360, kPhiFine, 0, 360), hRBlockFlat("R_block_thetaphi", "", nFlatT, 0, nFlatT, nFlatT, 0, nFlatT);
    // direct matrices with the recon angles in the neutron frame (the other four are R_*, R_det_*, R_*_recon, R_det_*_recon): index = truth (CE lab, CE neutron frame) x 2 + reco (recon block, recon); each [theta, phi, flat]
    TH2D* hNew[4][3];
    // Fswitch: the frame switches. states 0 recon (estimated neutron frame, thetaRe), 1 recon in the true neutron frame, 2 recon block (estimated neutron frame), 3 recon block in the true neutron frame,
    // 4 CE lab, 5 CE estimated neutron frame, 6 recon lab, 7 recon block lab, 8 CE true neutron frame; pair k goes from exFrom[k] to exTo[k]
    const int kEx = 7, kSt = 9;
    const int exFrom[kEx] = {1, 3, 4, 8, 6, 7, 2}, exTo[kEx] = {0, 2, 5, 5, 0, 2, 7};
    const char* exName[kEx] = {"M_rec_Fswitch", "M_recB_Fswitch", "M_CE_lab_to_est", "M_CE_neu_to_est", "M_rec_lab_to_est", "M_recB_lab_to_est", "M_recB_est_to_lab"};
    const char* stName[kSt] = {"recon_labF", "recon_neuF", "reconB_labF", "reconB_neuF", "CE_labF", "CE_estF", "recon_rawLab", "reconB_rawLab", "CE_neuF"};
    double maxChkTh = 0, maxChkPh = 0; // largest difference of the recon lab angles rotated by the rebuilt estimated neutron direction from thetaRe, phiRe
    TH1D* tEx[kSt][3]; // truth of the "from" state: theta, phi, flat of every row whose theta in that state is in range
    TProfile* apSt[kSt]; // A_eff (AEffOfTheta at the own theta of each proton, in that state's frame) per theta slice, for the cos fits of UnfoldPhi.cpp
    for (int st = 0; st < kSt; st++) apSt[st] = new TProfile(Form("AP_%s", stName[st]), "", kSlices, 0, tMax);
    TH2D* hEx[kEx][3];
    for (int st = 0; st < kSt; st++){
        tEx[st][0] = new TH1D(Form("T_%s_theta", stName[st]), "", nTheta, 0, tMax); tEx[st][1] = new TH1D(Form("T_%s_phi_all", stName[st]), "", kPhiFine, 0, 360); tEx[st][2] = new TH1D(Form("T_%s_thetaphi", stName[st]), "", nFlatT, 0, nFlatT);
    }
    for (int k = 0; k < kEx; k++){
        const std::string tag = Form("%s_to_%s", stName[exFrom[k]], stName[exTo[k]]);
        hEx[k][0] = new TH2D(Form("R_%s_theta", tag.c_str()), "", nTheta, 0, tMax, nTheta, 0, tMax);
        hEx[k][1] = new TH2D(Form("R_%s_phi_all", tag.c_str()), "", kPhiFine, 0, 360, kPhiFine, 0, 360);
        hEx[k][2] = new TH2D(Form("R_%s_thetaphi", tag.c_str()), "", nFlatT, 0, nFlatT, nFlatT, 0, nFlatT);
    }
    {
        const char* tags[4] = {"CElabF_reconBneuF", "CElabF_reconneuF", "CEneuF_reconBneuF", "CEneuF_reconneuF"};
        for (int k = 0; k < 4; k++){
            hNew[k][0] = new TH2D(Form("R_%s_theta", tags[k]), "", nTheta, 0, tMax, nTheta, 0, tMax);
            hNew[k][1] = new TH2D(Form("R_%s_phi_all", tags[k]), "", kPhiFine, 0, 360, kPhiFine, 0, 360);
            hNew[k][2] = new TH2D(Form("R_%s_thetaphi", tags[k]), "", nFlatT, 0, nFlatT, nFlatT, 0, nFlatT);
        }
    }
    TH2D hRDetThetaRe("R_det_theta_recon", "", nTheta, 0, tMax, nTheta, 0, tMax), hRDetPhiRe("R_det_phi_all_recon", "", kPhiFine, 0, 360, kPhiFine, 0, 360), hRDetFlatRe("R_det_thetaphi_recon", "", nFlatT, 0, nFlatT, nFlatT, 0, nFlatT);
    TH1D hTFlat("T_thetaphi", "", nFlatT, 0, nFlatT), hMFlat("M_thetaphi", "", nFlatT, 0, nFlatT);
    std::vector<TH2D*> hRSlice;
    std::vector<TH1D*> hTSlice;
    for (int k = 0; k < kSlices; k++){
        hRSlice.push_back(new TH2D(Form("R_phi_%d", k), Form("true theta %.2f to %.2f deg", k * tMax / kSlices, (k + 1) * tMax / kSlices), kPhiFine, 0, 360, kPhiFine, 0, 360));
        hTSlice.push_back(new TH1D(Form("T_phi_%d", k), "", kPhiFine, 0, 360));
    }
    const double rt = 0.5 * ThetaScale(); // residual range, deg
    // per reco type (recon block lab, recon block neutron frame, recon lab, recon neutron frame): residuals against the CE neutron frame truth (N) and the CE lab truth (L), and the dilution sums
    const char* reDir[4] = {"2_reconB_labF", "4_reconB_neuF", "1_recon_labF", "3_recon_neuF"};
    const char* reText[4] = {"Recon block", "Recon block", "Recon", "Recon"};
    const char* reFrame[4] = {"neutron estimate frame", "real neutron frame", "neutron estimate frame", "real neutron frame"};
    struct ResSet {TH1D *dThN, *dThL, *dPhN, *dPhL; std::vector<double> sN, s2N, sL, s2L, nS, pN, pSN, pS2N, pL, pSL, pS2L;}; // p*: theta residual per truth phi bin (36 bins of 10 deg), count, sum, sum of squares, for the real neutron frame (N) and lab (L) truth
    ResSet rs[4];
    for (int r = 0; r < 4; r++){
        rs[r].dThN = new TH1D(Form("dTheta_neutron_%s", reDir[r]), "", 100, -rt, rt); rs[r].dThL = new TH1D(Form("dTheta_lab_%s", reDir[r]), "", 100, -rt, rt);
        rs[r].dPhN = new TH1D(Form("dPhi_neutron_%s", reDir[r]), "", 90, -180, 180); rs[r].dPhL = new TH1D(Form("dPhi_lab_%s", reDir[r]), "", 90, -180, 180);
        rs[r].sN.assign(nTheta, 0); rs[r].s2N.assign(nTheta, 0); rs[r].sL.assign(nTheta, 0); rs[r].s2L.assign(nTheta, 0); rs[r].nS.assign(nTheta, 0);
        for (std::vector<double>* v : {&rs[r].pN, &rs[r].pSN, &rs[r].pS2N, &rs[r].pL, &rs[r].pSL, &rs[r].pS2L}) v->assign(36, 0);
    }

    // dilution per reco theta slice: sums of cos(dphi), cos(dphi)^2, count, for the neutron truth and the lab truth (in ResSet)
    long nRows = 0, nNoDir = 0, nPair = 0, nNeuOut = 0;
    // lab theta (deg) of the true neutron direction (CENDirX/Y/Z, x) against the estimated neutron direction (beam origin to the tube centre at the block middle, y), one entry per CE proton row
    TH2D hNeuTheta("N_true_vs_est_theta", "", 50, 0, 1.0, 50, 0, 1.0);

    // thetaCE, phiCE and thetaRe, thetaReD, phiRe, phiReD of cut_APandT.root are the CE in the TRUE neutron frame (CETHETA_NEUFRAME) and the recon angles in the ESTIMATED neutron frame (RECONTHETA_ESTNEU, RECONTHETABL_ESTNEU).
    // The lab CE truth comes from ChargeExchange_values (dirs); the recon angles rotated by the TRUE neutron (the "recon neutron frame" states) are made from the lab recon angles of AnglesRaw (outTheta ...), matched on (fEvent, ParentID).
    std::unordered_map<ULong64_t, LabRecon> labRecon;
    if (TTree* raw = OpenCut("cuts/cut_angles_all.root", "AnglesRaw")){
        // the CE protons of AnglesRaw only, read on all threads with the filter in the dataframe (Columns, TreeColumns.h), entry order: the last row of a key wins as before
        const auto c = Columns(raw, {"fEvent", "ParentID", "outDTheta", "outDPhi", "outTheta", "outPhi", "tubeX0", "tubeY0", "tubeX1", "tubeY1", "tubeZ0", "tubeZ1"}, "PDG == 2212 && isFromCE");
        labRecon.reserve(c[0].size());
        for (size_t i = 0; i < c[0].size(); i++)
            labRecon[EvKey((int)c[0][i], (int)c[1][i])] = {c[2][i] * TMath::RadToDeg(), Wrap360(c[3][i] * TMath::RadToDeg()), c[4][i] * TMath::RadToDeg(), Wrap360(c[5][i] * TMath::RadToDeg()), 0.5 * (c[6][i] + c[8][i]), 0.5 * (c[7][i] + c[9][i]), 0.5 * (c[10][i] + c[11][i])};
        printf(">>> SmearMatrix AnglesRaw CE protons (%zu) read took %.0f s\n", c[0].size(), Secs()); fflush(stdout);
    } else {printf("SmearMatrix: cut_angles_all.root AnglesRaw missing\n"); return;}
    TTreeReader rc(ce), rr(red), rre(rec);
    TTreeReaderValue<int> ev(rc, "fEvent"), tr(rc, "TrackID"), ev2(rr, "fEvent"), tr2(rr, "TrackID");
    TTreeReaderValue<double> thCE(rc, "thetaCE"), phCE(rc, "phiCE"), thR(rr, "thetaReD"), phR(rr, "phiReD"), thRe(rre, "thetaRe"), phRe(rre, "phiRe");
    while (rc.Next() && rr.Next() && rre.Next()){
        nRows++;
        if (*ev != *ev2 || *tr != *tr2){printf("SmearMatrix: CE and CalcValsReD rows do not line up at row %ld\n", nRows); return;}
        auto it = dirs.find(EvKey(*ev, *tr));
        if (it == dirs.end()){nNoDir++; continue;}
        const double thT = *thCE, phT = Wrap360(*phCE), thCEv = it->second.thLab, phLab = it->second.phLab; // CE truth in the true real neutron frame (thetaCE, phiCE) and in the lab
        const double thRec = *thR, phRec = Wrap360(*phR);
        const bool tIn = thT >= 0 && thT < tMax, rIn = thRec >= 0 && thRec < tMax;
        const double thRc = *thRe, phRc = Wrap360(*phRe); // recon theta and phi
        if (thRc >= 0 && thRc < tMax) hMFlatRe.Fill(FlatBin(thRc, phRc, tMax) + 0.5);
        if (tIn && thCEv >= 0 && thCEv < tMax){hRBeamTheta.Fill(thT, thCEv); hRBeamPhi.Fill(phT, phLab); hRBeamFlat.Fill(FlatBin(thT, phT, tMax) + 0.5, FlatBin(thCEv, phLab, tMax) + 0.5);}
        if (tIn && thRc >= 0 && thRc < tMax){hRThetaRe.Fill(thT, thRc); hRPhiRe.Fill(phT, phRc); hRFlatRe.Fill(FlatBin(thT, phT, tMax) + 0.5, FlatBin(thRc, phRc, tMax) + 0.5);}
        double thRecN, phRecN, thRcN, phRcN; // recon block and recon angles in the real neutron frame
        const auto lr = labRecon.find(EvKey(*ev, *tr));
        if (lr != labRecon.end()){RelToNeutron(lr->second.thB, lr->second.phB, it->second, thRecN, phRecN); RelToNeutron(lr->second.th, lr->second.ph, it->second, thRcN, phRcN);}
        else thRecN = phRecN = thRcN = phRcN = NAN;
        {   // direct matrices of the recon angles in the real neutron frame, against the CE lab and CE real neutron frame truth
            const double trTh[2] = {thCEv, thT}, trPh[2] = {phLab, phT}, reTh[2] = {thRecN, thRcN}, rePh[2] = {phRecN, phRcN};
            for (int tr = 0; tr < 2; tr++) for (int re = 0; re < 2; re++){
                if (!(trTh[tr] >= 0 && trTh[tr] < tMax && reTh[re] >= 0 && reTh[re] < tMax)) continue;
                TH2D** h = hNew[tr * 2 + re];
                h[0]->Fill(trTh[tr], reTh[re]); h[1]->Fill(trPh[tr], rePh[re]);
                h[2]->Fill(FlatBin(trTh[tr], trPh[tr], tMax) + 0.5, FlatBin(reTh[re], rePh[re], tMax) + 0.5);
            }
        }
        {   // Fswitch: truth of each state and the state to state matrices
            double thCEe = NAN, phCEe = NAN, thRcLab = NAN, phRcLab = NAN, thRecLab = NAN, phRecLab = NAN; // CE in the estimated real neutron frame, recon and recon block lab angles
            if (lr != labRecon.end()){
                thRcLab = lr->second.th; phRcLab = lr->second.ph; thRecLab = lr->second.thB; phRecLab = lr->second.phB;
                if (std::isfinite(lr->second.ez)){
                    const NDir est = {lr->second.ex, lr->second.ey, lr->second.ez, 0, 0};
                    const double thNeuT = std::acos(std::max(-1.0, std::min(1.0, it->second.z / std::sqrt(it->second.x * it->second.x + it->second.y * it->second.y + it->second.z * it->second.z)))) * TMath::RadToDeg();
                    const double thNeuE = std::acos(std::max(-1.0, std::min(1.0, est.z / std::sqrt(est.x * est.x + est.y * est.y + est.z * est.z)))) * TMath::RadToDeg();
                    if (thNeuT < 1.0 && thNeuE < 1.0) hNeuTheta.Fill(thNeuT, thNeuE); else nNeuOut++;
                    RelToNeutron(thCEv, phLab, est, thCEe, phCEe);
                    double thChk, phChk; RelToNeutron(thRcLab, phRcLab, est, thChk, phChk);
                    if (std::isfinite(thRc) && thRc > 0.05){maxChkTh = std::max(maxChkTh, std::fabs(thChk - thRc)); maxChkPh = std::max(maxChkPh, std::fabs(Wrap180(phChk - phRc)));}
                }
            }
            const double sTh[kSt] = {thRc, thRcN, thRec, thRecN, thCEv, thCEe, thRcLab, thRecLab, thT}, sPh[kSt] = {phRc, phRcN, phRec, phRecN, phLab, phCEe, phRcLab, phRecLab, phT};
            for (int st = 0; st < kSt; st++) if (sTh[st] >= 0 && sTh[st] < tMax){tEx[st][0]->Fill(sTh[st]); apSt[st]->Fill(sTh[st], AEffOfTheta(sTh[st], apScale)); tEx[st][1]->Fill(sPh[st]); tEx[st][2]->Fill(FlatBin(sTh[st], sPh[st], tMax) + 0.5);}
            for (int k = 0; k < kEx; k++){
                const int f = exFrom[k], t = exTo[k];
                if (!(sTh[f] >= 0 && sTh[f] < tMax && sTh[t] >= 0 && sTh[t] < tMax)) continue;
                hEx[k][0]->Fill(sTh[f], sTh[t]); hEx[k][1]->Fill(sPh[f], sPh[t]);
                hEx[k][2]->Fill(FlatBin(sTh[f], sPh[f], tMax) + 0.5, FlatBin(sTh[t], sPh[t], tMax) + 0.5);
            }
        }
        {   // residuals and dilution of each reco type against the CE real neutron frame truth and the CE lab truth
            const double reThs[4] = {thRec, thRecN, thRc, thRcN}, rePhs[4] = {phRec, phRecN, phRc, phRcN};
            for (int r = 0; r < 4; r++){
                if (!(reThs[r] >= 0 && reThs[r] < tMax)) continue;
                ResSet& q = rs[r];
                const double dN = Wrap180(rePhs[r] - phT), dL = Wrap180(rePhs[r] - phLab);
                const int b = std::min(nTheta - 1, (int)(reThs[r] / kStep));
                q.nS[b] += 1; q.sN[b] += std::cos(dN * TMath::DegToRad()); q.s2N[b] += std::pow(std::cos(dN * TMath::DegToRad()), 2);
                q.sL[b] += std::cos(dL * TMath::DegToRad()); q.s2L[b] += std::pow(std::cos(dL * TMath::DegToRad()), 2);
                if (tIn){
                    q.dThN->Fill(reThs[r] - thT); q.dThL->Fill(reThs[r] - thCEv); q.dPhN->Fill(dN); q.dPhL->Fill(dL);
                    const double rN = reThs[r] - thT, rL = reThs[r] - thCEv;
                    const int pbN = std::min(35, (int)(phT / 10)), pbL = std::min(35, (int)(phLab / 10));
                    q.pN[pbN] += 1; q.pSN[pbN] += rN; q.pS2N[pbN] += rN * rN;
                    q.pL[pbL] += 1; q.pSL[pbL] += rL; q.pS2L[pbL] += rL * rL;
                }
            }
        }
        if (thCEv >= 0 && thCEv < tMax){
            const int flatL = FlatBin(thCEv, phLab, tMax);
            hTDetTheta.Fill(thCEv); hTDetPhi.Fill(phLab); hTDetFlat.Fill(flatL + 0.5);
            if (rIn){hRDetTheta.Fill(thCEv, thRec); hRDetPhi.Fill(phLab, phRec); hRDetFlat.Fill(flatL + 0.5, FlatBin(thRec, phRec, tMax) + 0.5);}
            if (thRc >= 0 && thRc < tMax){hRDetThetaRe.Fill(thCEv, thRc); hRDetPhiRe.Fill(phLab, phRc); hRDetFlatRe.Fill(flatL + 0.5, FlatBin(thRc, phRc, tMax) + 0.5);}
        }
        if (thRc >= 0 && thRc < tMax){
            hTRecTheta.Fill(thRc); hTRecPhi.Fill(phRc);
            if (rIn){hRBlockTheta.Fill(thRc, thRec); hRBlockPhi.Fill(phRc, phRec); hRBlockFlat.Fill(FlatBin(thRc, phRc, tMax) + 0.5, FlatBin(thRec, phRec, tMax) + 0.5);}
        }
        const int kT = std::min(kSlices - 1, (int)(thT / (tMax / kSlices)));
        const int flatT = FlatBin(thT, phT, tMax), flatR = FlatBin(thRec, phRec, tMax);
        if (tIn){
            hTTheta.Fill(thT); hTPhi.Fill(phT); hTSlice[kT]->Fill(phT); hTFlat.Fill(flatT + 0.5); hAPTheta.Fill(thT, AEffOfTheta(thT, apScale));
        }
        if (rIn){
            hMTheta.Fill(thRec);
            hMFlat.Fill(flatR + 0.5);
            if (tIn){
                nPair++;
                hRTheta.Fill(thT, thRec); hRPhi.Fill(phT, phRec); hRSlice[kT]->Fill(phT, phRec); hRFlat.Fill(flatT + 0.5, flatR + 0.5);
            }
        }
    }
    printf("SmearMatrix: %ld CE rows, %ld without a neutron direction, %ld in the matrix\n", nRows, nNoDir, nPair);

    printf(">>> SmearMatrix main loop over the CE, block and recon rows (%ld rows) took %.0f s\n", nRows, Secs()); fflush(stdout);
    // SiCE protons: the recon block protons that are not children of a carbon CE track (silicon CE in a noPhys run), recon theta and phi only
    TH1D hSiCETheta("SiCE_theta", "", nTheta, 0, tMax), hSiCEFlat("SiCE_thetaphi", "", nFlatT, 0, nFlatT), hSiCEFlatRe("SiCE_thetaphi_recon", "", nFlatT, 0, nFlatT);
    TH2D hSiCEThPh("SiCE_theta_phi", "", nTheta, 0, tMax, kPhiFine, 0, 360);
    long nSiCE = 0;
    if (TTree* raw = OpenCut("cuts/cut_angles_all.root", "AnglesRaw")){
        // the non-CE protons with a theta in range (block or recon), read on all threads; the histograms do not depend on the row order (ColumnsUnordered, TreeColumns.h); NaN (no tube) fails the comparisons
        const std::string inRange = Form("(thetaRelND * %.17g >= 0 && thetaRelND * %.17g < %.17g) || (thetaRelN * %.17g >= 0 && thetaRelN * %.17g < %.17g)", TMath::RadToDeg(), TMath::RadToDeg(), tMax, TMath::RadToDeg(), TMath::RadToDeg(), tMax);
        const auto sc = ColumnsUnordered(raw, {"thetaRelND", "phiRelND", "thetaRelN", "phiRelN"}, "PDG == 2212 && isFromCE == 0 && (" + inRange + ")");
        printf(">>> SmearMatrix AnglesRaw SiCE protons (%zu rows) read took %.0f s\n", sc[0].size(), Secs()); fflush(stdout);
        for (size_t i = 0; i < sc[0].size(); i++){
            const double th = sc[0][i] * TMath::RadToDeg(), ph = Wrap360(sc[1][i] * TMath::RadToDeg());
            const double thr = sc[2][i] * TMath::RadToDeg(), phr = Wrap360(sc[3][i] * TMath::RadToDeg()); // recon theta and phi (the recon matrix of UnfoldPhi.cpp)
            if (thr >= 0 && thr < tMax) hSiCEFlatRe.Fill(FlatBin(thr, phr, tMax) + 0.5);
            if (!(th >= 0 && th < tMax)) continue; // NaN (no tube) fails
            hSiCETheta.Fill(th); hSiCEFlat.Fill(FlatBin(th, ph, tMax) + 0.5); hSiCEThPh.Fill(th, ph);
            nSiCE++;
        }
    } else printf("SmearMatrix: cut_angles_all.root AnglesRaw missing, no SiCE protons\n");
    printf("SmearMatrix: %ld SiCE protons (silicon CE, not from a carbon CE track) with recon block theta 0 to %.1f deg\n", nSiCE, tMax);

    // dilution csv: D = mean cos(reco phi - true phi) per reco theta slice, one csv and one pair of graphs per reco type
    TGraphErrors gN[4], gL[4], pgN[4], pgL[4]; // pg*: mean theta residual (point) and RMS (error bar) per truth phi bin
    for (int r = 0; r < kResidualTypes; r++){
        FILE* csv =fopen(Form("plots/8_smear/6_residuals/%s/smear_dilution.csv", reDir[r]), "w");
        fprintf(csv, "theta_lo,theta_hi,n,D_neutron,D_neutron_err,D_lab,D_lab_err\n");
        const ResSet& q = rs[r];
        for (int b = 0; b < nTheta; b++){
            if (q.nS[b] < 2) continue;
            const double n = q.nS[b];
            const double dN = q.sN[b] / n, dL = q.sL[b] / n;
            const double eN = std::sqrt(std::max(0.0, q.s2N[b] / n - dN * dN) / n), eL = std::sqrt(std::max(0.0, q.s2L[b] / n - dL * dL) / n);
            fprintf(csv, "%.4f,%.4f,%.0f,%.6f,%.6f,%.6f,%.6f\n", b * kStep, (b + 1) * kStep, n, dN, eN, dL, eL);
            const double x = (b + 0.5) * kStep;
            gN[r].SetPoint(gN[r].GetN(), x, dN); gN[r].SetPointError(gN[r].GetN() - 1, 0.5 * kStep, eN);
            gL[r].SetPoint(gL[r].GetN(), x, dL); gL[r].SetPointError(gL[r].GetN() - 1, 0.5 * kStep, eL);
        }
        fclose(csv);
        FILE* csvP = fopen(Form("plots/8_smear/6_residuals/%s/theta_residual_vs_phi.csv", reDir[r]), "w");
        fprintf(csvP, "phi_lo,phi_hi,n_neutron,mean_neutron,rms_neutron,n_lab,mean_lab,rms_lab\n");
        for (int b = 0; b < 36; b++){
            if (q.pN[b] < 2 || q.pL[b] < 2) continue;
            const double mN = q.pSN[b] / q.pN[b], mL = q.pSL[b] / q.pL[b];
            const double rN = std::sqrt(std::max(0.0, q.pS2N[b] / q.pN[b] - mN * mN)), rL = std::sqrt(std::max(0.0, q.pS2L[b] / q.pL[b] - mL * mL));
            fprintf(csvP, "%d,%d,%.0f,%.6f,%.6f,%.0f,%.6f,%.6f\n", b * 10, (b + 1) * 10, q.pN[b], mN, rN, q.pL[b], mL, rL);
            pgN[r].SetPoint(pgN[r].GetN(), b * 10 + 5.0, mN); pgN[r].SetPointError(pgN[r].GetN() - 1, 5.0, rN);
            pgL[r].SetPoint(pgL[r].GetN(), b * 10 + 5.0, mL); pgL[r].SetPointError(pgL[r].GetN() - 1, 5.0, rL);
        }
        fclose(csvP);
    }

    // column-normalised copies are what is plotted, the root file keeps the counts and the truth/recon histograms
    TH2D* nTheta2 = (TH2D*)hRTheta.Clone("R_theta_norm");
    TH2D* nFlat2 = (TH2D*)hRFlat.Clone("R_thetaphi_norm");
    TH2D* nPhi2 = (TH2D*)hRPhi.Clone("R_phi_all_norm");
    NormColumns(nTheta2, &hTTheta);
    NormColumns(nFlat2, &hTFlat);
    NormColumns(nPhi2, &hTPhi);
    TH2D* nBTh = (TH2D*)hRBeamTheta.Clone("R_beam_theta_norm");
    TH2D* nBPh = (TH2D*)hRBeamPhi.Clone("R_beam_phi_all_norm");
    NormColumns(nBTh, &hTTheta);
    NormColumns(nBPh, &hTPhi);
    TH2D* nBFl = (TH2D*)hRBeamFlat.Clone("R_beam_thetaphi_norm");
    NormColumns(nBFl, &hTFlat);

    gStyle->SetOptStat(0);
    PlainFonts(); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X"); // as PlotAngle.cpp
    const std::string sd = Form("beam tilt max %.2f deg", maxTilt);
    const std::string cutPair = Form("CE recon protons, both thetas 0 to %.1f deg, no other cuts", tMax);
    // the same four matrices as proton counts per bin (x_2_ next to x_), no column left blank
    // beam angle only (no detector): the lab angle the proton has about z against its angle about the neutron
    const std::string cutBeam = Form("CE protons, both thetas 0 to %.1f deg, no other cuts", tMax);
    SaveCounts(hNeuTheta, "plots/8_smear/3_matrix3/1_neutron_theta_true_vs_est.pdf", Form("Neutron theta, true (CENDir) vs estimated (beam origin to tube centre), lab frame, %s, CE protons with a tube, both thetas 0 to 1 deg (%ld outside), no other cuts;True neutron theta, lab (deg);Estimated neutron theta, lab (deg)", sd.c_str(), nNeuOut));
    SavePlot2D(nBTh, "plots/8_smear/1_Fswitch/1_M_CE_Fswitch/1_M_CE_Fswitch_theta.pdf", Form("M_CE_Fswitch: Response of the beam angle, CE theta (neutron estimate frame), %s, %s;CE theta, real neutron frame (deg);CE theta, neutron estimate frame (deg)", sd.c_str(), cutBeam.c_str()), ThetaAxisMax());
    SaveCounts(hRBeamTheta, "plots/8_smear/1_Fswitch/1_M_CE_Fswitch/4_M_CE_Fswitch_theta_counts.pdf", Form("M_CE_Fswitch: Protons per bin, beam angle, CE theta (neutron estimate frame), %s, %s;CE theta, real neutron frame (deg);CE theta, neutron estimate frame (deg)", sd.c_str(), cutBeam.c_str()), ThetaAxisMax());
    SavePlot2D(nBPh, "plots/8_smear/1_Fswitch/1_M_CE_Fswitch/2_M_CE_Fswitch_phi.pdf", Form("M_CE_Fswitch: Response of the beam angle, CE phi (neutron estimate frame), %s, %s;CE phi, real neutron frame (deg);CE phi, neutron estimate frame (deg)", sd.c_str(), cutBeam.c_str()));
    SaveCounts(hRBeamPhi, "plots/8_smear/1_Fswitch/1_M_CE_Fswitch/5_M_CE_Fswitch_phi_counts.pdf", Form("M_CE_Fswitch: Protons per bin, beam angle, CE phi (neutron estimate frame), %s, %s;CE phi, real neutron frame (deg);CE phi, neutron estimate frame (deg)", sd.c_str(), cutBeam.c_str()));
    const std::string flatDef = Form("bin number = theta slice x %d + phi bin (%d theta slices of %.2f deg, %d phi bins of %.1f deg)", kPhiFlat, kSlices, tMax / kSlices, kPhiFlat, 360.0 / kPhiFlat);
    SavePlot2D(nBFl, "plots/8_smear/1_Fswitch/1_M_CE_Fswitch/3_M_CE_Fswitch_thetaphi.pdf", Form("M_CE_Fswitch: Response of the beam angle, CE (theta, phi) (neutron estimate frame), %s, %s, %s;CE (theta, phi) bin number, real neutron frame;CE (theta, phi) bin number, neutron estimate frame", sd.c_str(), cutBeam.c_str(), flatDef.c_str()));
    SaveCounts(hRBeamFlat, "plots/8_smear/1_Fswitch/1_M_CE_Fswitch/6_M_CE_Fswitch_thetaphi_counts.pdf", Form("M_CE_Fswitch: Protons per bin, beam angle, CE (theta, phi) (neutron estimate frame), %s, %s, %s;CE (theta, phi) bin number, real neutron frame;CE (theta, phi) bin number, neutron estimate frame", sd.c_str(), cutBeam.c_str(), flatDef.c_str()));
    // detector only (no beam angle): the recon angle against the true lab angle of the CE proton, the first step of UnfoldPhi.cpp
    TH2D* nDTh = (TH2D*)hRDetTheta.Clone("R_det_theta_norm"); TH2D* nDPh = (TH2D*)hRDetPhi.Clone("R_det_phi_all_norm"); TH2D* nDFl = (TH2D*)hRDetFlat.Clone("R_det_thetaphi_norm");
    TH2D* nDThRe = (TH2D*)hRDetThetaRe.Clone("R_det_theta_recon_norm"); TH2D* nDPhRe = (TH2D*)hRDetPhiRe.Clone("R_det_phi_all_recon_norm"); TH2D* nDFlRe = (TH2D*)hRDetFlatRe.Clone("R_det_thetaphi_recon_norm");
    NormColumns(nDTh, &hTDetTheta); NormColumns(nDPh, &hTDetPhi); NormColumns(nDFl, &hTDetFlat);
    NormColumns(nDThRe, &hTDetTheta); NormColumns(nDPhRe, &hTDetPhi); NormColumns(nDFlRe, &hTDetFlat);
    const std::string cutDet = Form("CE protons, both thetas 0 to %.1f deg, no other cuts", tMax);
    // physics step, CE true neutron frame to recon estimated neutron frame (the truth is the CE theta, phi of the true neutron frame, T_theta, T_phi_all, T_thetaphi)
    TH2D* nNRTh = (TH2D*)hRThetaRe.Clone("R_theta_recon_norm"); TH2D* nNRPh = (TH2D*)hRPhiRe.Clone("R_phi_all_recon_norm"); TH2D* nNRFl = (TH2D*)hRFlatRe.Clone("R_thetaphi_recon_norm");
    NormColumns(nNRTh, &hTTheta); NormColumns(nNRPh, &hTPhi); NormColumns(nNRFl, &hTFlat);
    SavePlot2D(nNRTh, "plots/8_smear/2_M_physics/1_M_physics_theta.pdf", Form("M_physics: Physics response of recon theta (CE true real neutron frame vs recon estimated real neutron frame), %s;CE theta, true real neutron frame (deg);Recon theta, estimated real neutron frame (deg)", cutDet.c_str()), ThetaAxisMax());
    SaveCounts(hRThetaRe, "plots/8_smear/2_M_physics/4_M_physics_theta_counts.pdf", Form("M_physics: Protons per bin, physics response of recon theta (CE true real neutron frame vs recon estimated real neutron frame), %s;CE theta, true real neutron frame (deg);Recon theta, estimated real neutron frame (deg)", cutDet.c_str()), ThetaAxisMax());
    SavePlot2D(nNRPh, "plots/8_smear/2_M_physics/2_M_physics_phi.pdf", Form("M_physics: Physics response of recon phi (CE true real neutron frame vs recon estimated real neutron frame), %s;CE phi, true real neutron frame (deg);Recon phi, estimated real neutron frame (deg)", cutDet.c_str()));
    SaveCounts(hRPhiRe, "plots/8_smear/2_M_physics/5_M_physics_phi_counts.pdf", Form("M_physics: Protons per bin, physics response of recon phi (CE true real neutron frame vs recon estimated real neutron frame), %s;CE phi, true real neutron frame (deg);Recon phi, estimated real neutron frame (deg)", cutDet.c_str()));
    SavePlot2D(nNRFl, "plots/8_smear/2_M_physics/3_M_physics_thetaphi.pdf", Form("M_physics: Physics response of recon (theta, phi) (CE true real neutron frame vs recon estimated real neutron frame), %s, %s;CE (theta, phi) bin number, true real neutron frame;Recon (theta, phi) bin number, estimated real neutron frame", cutDet.c_str(), flatDef.c_str()));
    SaveCounts(hRFlatRe, "plots/8_smear/2_M_physics/6_M_physics_thetaphi_counts.pdf", Form("M_physics: Protons per bin, physics response of recon (theta, phi) (CE true real neutron frame vs recon estimated real neutron frame), %s, %s;CE (theta, phi) bin number, true real neutron frame;Recon (theta, phi) bin number, estimated real neutron frame", cutDet.c_str(), flatDef.c_str()));
    {   // all 8 direct matrices: truth (CE neutron estimate frame, CE real neutron frame) x reco (recon block lab, recon block neutron, recon lab, recon neutron), probability pages 1-3 then counts pages 4-6
        const char* trName[2] = {"CE_labF", "CE_neuF"}, *trText[2] = {"CE (neutron estimate frame)", "CE (real neutron frame)"}, *trFrame[2] = {"neutron estimate frame", "real neutron frame"};
        const char* reName[4] = {"reconB_labF", "reconB_neuF", "recon_labF", "recon_neuF"}, *reText[4] = {"recon block", "recon block", "recon", "recon"}, *reFrame[4] = {"neutron estimate frame", "real neutron frame", "neutron estimate frame", "real neutron frame"};
        const TH1D* truthH[2][3] = {{&hTDetTheta, &hTDetPhi, &hTDetFlat}, {&hTTheta, &hTPhi, &hTFlat}};
        const TH2D* cnt[2][4][3] = {
            {{&hRDetTheta, &hRDetPhi, &hRDetFlat}, {hNew[0][0], hNew[0][1], hNew[0][2]}, {&hRDetThetaRe, &hRDetPhiRe, &hRDetFlatRe}, {hNew[1][0], hNew[1][1], hNew[1][2]}},
            {{&hRTheta, &hRPhi, &hRFlat}, {hNew[2][0], hNew[2][1], hNew[2][2]}, {&hRThetaRe, &hRPhiRe, &hRFlatRe}, {hNew[3][0], hNew[3][1], hNew[3][2]}}};
        const char* angName[3] = {"theta", "phi", "thetaphi"};
        for (int tr = 0; tr < 2; tr++) for (int re = 0; re < 4; re++){
            const int reOrder[4] = {2, 4, 1, 3}; // folder number of each reco type: neutron estimate frame ones first
            const std::string dir = Form("plots/8_smear/5_direct/%d_%s/%d_%s", tr + 1, trName[tr], reOrder[re], reName[re]);
            for (int a = 0; a < 3; a++){
                TH2D* pr = (TH2D*)cnt[tr][re][a]->Clone(Form("R_direct_%d_%d_%d_norm", tr, re, a));
                NormColumns(pr, truthH[tr][a]);
                const std::string stem = Form("%s_to_%s_%s", trName[tr], reName[re], angName[a]);
                std::string title, axes;
                if (a == 2){
                    title = Form("Direct response of %s (theta, phi), %s to %s (%s), %s, %s", reText[re], trText[tr], reText[re], reFrame[re], cutDet.c_str(), flatDef.c_str());
                    axes = Form(";CE (theta, phi) bin number, %s;Recon (theta, phi) bin number, %s", trFrame[tr], reFrame[re]);
                } else {
                    const char* an = a == 0 ? "theta" : "phi";
                    title = Form("Direct response of %s %s, %s to %s (%s), %s", reText[re], an, trText[tr], reText[re], reFrame[re], cutDet.c_str());
                    axes = Form(";CE %s, %s (deg);Recon %s, %s (deg)", an, trFrame[tr], an, reFrame[re]);
                }
                const double axMax = a == 0 ? ThetaAxisMax() : 0;
                SavePlot2D(pr, Form("%s/%d_%s.pdf", dir.c_str(), a + 1, stem.c_str()), title + axes, axMax);
                SaveCounts(*cnt[tr][re][a], Form("%s/%d_%s_counts.pdf", dir.c_str(), a + 4, stem.c_str()), "Protons per bin, " + title + axes, axMax);
                delete pr;
            }
        }
    }
    printf("SmearMatrix: recon lab angles rotated by the rebuilt estimated neutron direction minus thetaRe, phiRe: largest theta difference %.2e deg, largest phi difference %.2e deg\n", maxChkTh, maxChkPh);
    {   // Fswitch: probability pages 1-3, counts pages 4-6, of the state to state matrices
        const char* stText[kSt] = {"recon", "recon", "recon block", "recon block", "CE", "CE", "recon", "recon block", "CE"}, *stFrame[kSt] = {"neutron estimate frame", "real neutron frame", "neutron estimate frame", "real neutron frame", "lab frame", "estimated real neutron frame", "lab frame", "lab frame", "true real neutron frame"};
        const char* angName[3] = {"theta", "phi", "thetaphi"};
        for (int k = 0; k < kEx; k++){
            const int f = exFrom[k], t = exTo[k];
            const std::string dir = Form("plots/8_smear/1_Fswitch/%d_%s", k + 2, exName[k]);
            for (int a = 0; a < 3; a++){
                TH2D* pr = (TH2D*)hEx[k][a]->Clone(Form("R_extra_%d_%d_norm", k, a));
                NormColumns(pr, tEx[f][a]);
                const std::string stem = Form("%s_%s", exName[k], angName[a]);
                std::string title, axes;
                if (a == 2){
                    title = Form("Direct response of (theta, phi), %s (%s) to %s (%s), %s, %s", stText[f], stFrame[f], stText[t], stFrame[t], cutDet.c_str(), flatDef.c_str());
                    axes = Form(";%s (theta, phi) bin number, %s;%s (theta, phi) bin number, %s", stText[f], stFrame[f], stText[t], stFrame[t]);
                } else {
                    const char* an = a == 0 ? "theta" : "phi";
                    title = Form("Direct response of %s, %s (%s) to %s (%s), %s", an, stText[f], stFrame[f], stText[t], stFrame[t], cutDet.c_str());
                    axes = Form(";%s %s, %s (deg);%s %s, %s (deg)", stText[f], an, stFrame[f], stText[t], an, stFrame[t]);
                }
                const double axMax = a == 0 ? ThetaAxisMax() : 0;
                SavePlot2D(pr, Form("%s/%d_%s.pdf", dir.c_str(), a + 1, stem.c_str()), title + axes, axMax);
                SaveCounts(*hEx[k][a], Form("%s/%d_%s_counts.pdf", dir.c_str(), a + 4, stem.c_str()), "Protons per bin, " + title + axes, axMax);
                delete pr;
            }
        }
    }
    // block step: recon lab against recon block lab
    TH2D* nBlTh = (TH2D*)hRBlockTheta.Clone("R_block_theta_norm"); TH2D* nBlPh = (TH2D*)hRBlockPhi.Clone("R_block_phi_all_norm"); TH2D* nBlFl = (TH2D*)hRBlockFlat.Clone("R_block_thetaphi_norm");
    NormColumns(nBlTh, &hTRecTheta); NormColumns(nBlPh, &hTRecPhi); NormColumns(nBlFl, &hMFlatRe);
    SavePlot2D(nBlTh, "plots/8_smear/3_M_detector/1_M_detector_theta.pdf", Form("M_detector: Block response of theta, recon (neutron estimate frame) vs recon block (neutron estimate frame), %s;Recon theta, neutron estimate frame (deg);Recon block theta, neutron estimate frame (deg)", cutDet.c_str()), ThetaAxisMax());
    SaveCounts(hRBlockTheta, "plots/8_smear/3_M_detector/4_M_detector_theta_counts.pdf", Form("M_detector: Protons per bin, block response of theta, recon (neutron estimate frame) vs recon block (neutron estimate frame), %s;Recon theta, neutron estimate frame (deg);Recon block theta, neutron estimate frame (deg)", cutDet.c_str()), ThetaAxisMax());
    SavePlot2D(nBlPh, "plots/8_smear/3_M_detector/2_M_detector_phi.pdf", Form("M_detector: Block response of phi, recon (neutron estimate frame) vs recon block (neutron estimate frame), %s;Recon phi, neutron estimate frame (deg);Recon block phi, neutron estimate frame (deg)", cutDet.c_str()));
    SaveCounts(hRBlockPhi, "plots/8_smear/3_M_detector/5_M_detector_phi_counts.pdf", Form("M_detector: Protons per bin, block response of phi, recon (neutron estimate frame) vs recon block (neutron estimate frame), %s;Recon phi, neutron estimate frame (deg);Recon block phi, neutron estimate frame (deg)", cutDet.c_str()));
    SavePlot2D(nBlFl, "plots/8_smear/3_M_detector/3_M_detector_thetaphi.pdf", Form("M_detector: Block response of (theta, phi), recon (neutron estimate frame) vs recon block (neutron estimate frame), %s, %s;Recon (theta, phi) bin number, neutron estimate frame;Recon block (theta, phi) bin number, neutron estimate frame", cutDet.c_str(), flatDef.c_str()));
    SaveCounts(hRBlockFlat, "plots/8_smear/3_M_detector/6_M_detector_thetaphi_counts.pdf", Form("M_detector: Protons per bin, block response of (theta, phi), recon (neutron estimate frame) vs recon block (neutron estimate frame), %s, %s;Recon (theta, phi) bin number, neutron estimate frame;Recon block (theta, phi) bin number, neutron estimate frame", cutDet.c_str(), flatDef.c_str()));
    // the two steps multiplied: physics (CE true neutron frame to recon, estimated neutron frame) and block (recon to recon block, estimated neutron frame), against the direct matrix R_thetaphi
    TH2D* nChain = (TH2D*)hRFlat.Clone("R_detector_physics_thetaphi"); nChain->Reset();
    TH2D* nDiff = (TH2D*)hRFlat.Clone("R_detector_physics_minus_direct_thetaphi"); nDiff->Reset();
    double maxDiff = 0, maxEffDiff = 0;
    {
        const int n = nFlatT;
        for (int i = 0; i < n; i++){
            double effChain = 0, effDirect = 0;
            std::vector<double> row(n, 0.0);
            for (int k = 0; k < n; k++){
                const double q = nNRFl->GetBinContent(i + 1, k + 1);
                if (q == 0) continue;
                for (int l = 0; l < n; l++) row[l] += q * nBlFl->GetBinContent(k + 1, l + 1);
            }
            for (int l = 0; l < n; l++){
                const double c = row[l];
                const double d = nFlat2->GetBinContent(i + 1, l + 1);
                nChain->SetBinContent(i + 1, l + 1, c); nDiff->SetBinContent(i + 1, l + 1, c - d);
                maxDiff = std::max(maxDiff, std::abs(c - d));
                effChain += c; effDirect += d;
            }
            if (hTFlat.GetBinContent(i + 1) > 0) maxEffDiff = std::max(maxEffDiff, std::abs(effChain - effDirect));
        }
    }
    printf("SmearMatrix: M_detector x M_physics minus direct R_thetaphi, largest bin difference %.4f, largest efficiency difference of a column %.4f\n", maxDiff, maxEffDiff);
    SavePlot2D(nChain, "plots/8_smear/4_chain/2_M_detector_x_M_physics/3_M_detector_x_M_physics_thetaphi.pdf", Form("M_detector x M_physics response of (theta, phi), CE (true real neutron frame) to recon block (estimated real neutron frame), %s, %s;CE (theta, phi) bin number, true real neutron frame;Recon block (theta, phi) bin number, estimated real neutron frame", cutDet.c_str(), flatDef.c_str()));
    nDiff->SetMinimum(-std::max(maxDiff, 1e-6)); nDiff->SetMaximum(std::max(maxDiff, 1e-6));
    SavePlot2D(nDiff, "plots/8_smear/4_chain/2_M_detector_x_M_physics/6_M_detector_x_M_physics_minus_direct_thetaphi.pdf", Form("M_detector x M_physics minus direct response R_thetaphi, largest difference %.4f, %s, %s;CE (theta, phi) bin number, true real neutron frame;Recon block (theta, phi) bin number, estimated real neutron frame", maxDiff, cutDet.c_str(), flatDef.c_str()), 0, "chain - direct");
    // the same product for the theta and the phi matrices alone (phi depends on theta, so these are approximate), each full product against its direct matrix
    std::vector<TH1*> chain1D;
    auto chain1y = [&](const char* ang, int i0, const TH2D* ph, const TH2D* bl, const TH2D* direct, double axMax){
        TH2D* full = MatMul(*ph, *bl, Form("R_detector_physics_%s", ang));
        TH2D* diff = (TH2D*)full->Clone(Form("R_detector_physics_minus_direct_%s", ang)); diff->Add(direct, -1);
        double mx = 0;
        for (int i = 1; i <= diff->GetNbinsX(); i++) for (int j = 1; j <= diff->GetNbinsY(); j++) mx = std::max(mx, std::abs(diff->GetBinContent(i, j)));
        printf("SmearMatrix: M_detector x M_physics of %s minus direct, largest bin difference %.4f\n", ang, mx);
        SavePlot2D(full, Form("plots/8_smear/4_chain/2_M_detector_x_M_physics/%d_M_detector_x_M_physics_%s.pdf", i0, ang), Form("M_detector x M_physics response of %s, CE (true real neutron frame) to recon block (estimated real neutron frame), %s;CE %s, true real neutron frame (deg);Recon block %s, estimated real neutron frame (deg)", ang, cutDet.c_str(), ang, ang), axMax);
        diff->SetMinimum(-std::max(mx, 1e-6)); diff->SetMaximum(std::max(mx, 1e-6));
        SavePlot2D(diff, Form("plots/8_smear/4_chain/2_M_detector_x_M_physics/%d_M_detector_x_M_physics_minus_direct_%s.pdf", i0 + 3, ang), Form("M_detector x M_physics minus direct response of %s, largest difference %.4f, %s;CE %s, true real neutron frame (deg);Recon block %s, estimated real neutron frame (deg)", ang, mx, cutDet.c_str(), ang, ang), axMax, "chain - direct");
        chain1D.push_back(full); chain1D.push_back(diff);
    };
    chain1y("theta", 1, nNRTh, nBlTh, nTheta2, ThetaAxisMax());
    chain1y("phi", 2, nNRPh, nBlPh, nPhi2, 0);
    for (int r = 0; r < kResidualTypes; r++){
        const std::string dirR =Form("plots/8_smear/6_residuals/%s/", reDir[r]);
        SaveResidual(rs[r].dThN, rs[r].dThL, "theta", dirR + "1_residual_theta.pdf", Form("%s theta (%s) - CE theta (real neutron frame or neutron estimate frame), %s, %s;#Delta#theta (deg);Count", reText[r], reFrame[r], sd.c_str(), cutPair.c_str()));
        SaveResidual(rs[r].dPhN, rs[r].dPhL, "phi", dirR + "2_residual_phi.pdf", Form("%s phi (%s) - CE phi (real neutron frame or neutron estimate frame), %s, %s;#Delta#phi (deg);Count", reText[r], reFrame[r], sd.c_str(), cutPair.c_str()));
        const std::string path = dirR + "3_dilution_vs_theta.pdf";
        TCanvas c("cSmearD", "Plots", 1500, 1200);
        PageMargins(c);
        TMultiGraph mg;
        gN[r].SetMarkerStyle(20); gN[r].SetMarkerSize(1.2); gN[r].SetMarkerColor(kRed + 1); gN[r].SetLineColor(kRed + 1);
        gL[r].SetMarkerStyle(20); gL[r].SetMarkerSize(1.2); gL[r].SetMarkerColor(kBlue + 1); gL[r].SetLineColor(kBlue + 1);
        mg.Add(&gN[r], "PE"); mg.Add(&gL[r], "PE");
        mg.SetTitle(PlotTitle(path, Form("Dilution D = <cos(#phi_{%s, %s} - #phi_{CE, real neutron or neutron estimate frame})>, %s, CE recon protons, %s theta 0 to %.1f deg, no other cuts;%s theta, %s (deg);D", reText[r], reFrame[r], sd.c_str(), reText[r], tMax, reText[r], reFrame[r])).c_str());
        mg.Draw("A");
        mg.GetYaxis()->SetRangeUser(-0.1, 1.1);
        mg.GetXaxis()->SetLimits(0, ThetaAxisMax());
        TLegend leg(0.715, 0.76, 0.99, 0.88);
        leg.SetTextSize(0.017);
        leg.AddEntry(&gN[r], "CE phi (real neutron frame)", "pe");
        leg.AddEntry(&gL[r], "CE phi (neutron estimate frame)", "pe");
        leg.Draw();
        c.SaveAs(path.c_str());
        {   // theta residual against the truth phi: mean (point) and RMS (error bar) per 10 deg phi bin
            const std::string pathP = dirR + "4_theta_residual_vs_phi.pdf";
            TCanvas cp("cSmearP", "Plots", 1500, 1200);
            PageMargins(cp);
            TMultiGraph mgP;
            pgN[r].SetMarkerStyle(20); pgN[r].SetMarkerSize(1.2); pgN[r].SetMarkerColor(kRed + 1); pgN[r].SetLineColor(kRed + 1);
            pgL[r].SetMarkerStyle(20); pgL[r].SetMarkerSize(1.2); pgL[r].SetMarkerColor(kBlue + 1); pgL[r].SetLineColor(kBlue + 1);
            mgP.Add(&pgN[r], "PE"); mgP.Add(&pgL[r], "PE");
            mgP.SetTitle(PlotTitle(pathP, Form("%s theta (%s) - CE theta, mean (point) and RMS (bar) per truth phi bin of 10 deg, %s, %s;CE phi (deg);#Delta#theta (deg)", reText[r], reFrame[r], sd.c_str(), cutPair.c_str())).c_str());
            mgP.Draw("A");
            double lo = 0, hi = 0;
            for (TGraphErrors* g : {&pgN[r], &pgL[r]}) for (int i = 0; i < g->GetN(); i++){lo = std::min(lo, g->GetY()[i] - g->GetEY()[i]); hi = std::max(hi, g->GetY()[i] + g->GetEY()[i]);}
            const double pad = 0.1 * std::max(hi - lo, 1e-3);
            mgP.GetYaxis()->SetRangeUser(lo - pad, hi + pad);
            mgP.GetXaxis()->SetLimits(0, 360);
            TLine zero(0, 0, 360, 0); zero.SetLineStyle(2); zero.Draw();
            TLegend legP(0.715, 0.76, 0.99, 0.88);
            legP.SetTextSize(0.017);
            legP.AddEntry(&pgN[r], "CE theta (real neutron frame)", "pe");
            legP.AddEntry(&pgL[r], "CE theta (neutron estimate frame)", "pe");
            legP.Draw();
            cp.SaveAs(pathP.c_str());
        }
    }

    out.cd();
    for (TH1* h : std::vector<TH1*>{&hRTheta, &hTTheta, &hMTheta, &hRPhi, &hTPhi, &hRFlat, &hTFlat, &hMFlat, nTheta2, nFlat2, nPhi2, &hSiCETheta, &hSiCEFlat, &hSiCEThPh, &hRBeamTheta, &hRBeamPhi, nBTh, nBPh, &hRBeamFlat, nBFl, &hMFlatRe, &hSiCEFlatRe, &hAPTheta, &hTDetTheta, &hTDetPhi, &hTDetFlat, &hRDetTheta, &hRDetPhi, &hRDetFlat, &hRDetThetaRe, &hRDetPhiRe, &hRDetFlatRe, nDTh, nDPh, nDFl, nDThRe, nDPhRe, nDFlRe, &hTRecTheta, &hTRecPhi, &hRBlockTheta, &hRBlockPhi, &hRBlockFlat, nBlTh, nBlPh, nBlFl, nChain, &hRThetaRe, &hRPhiRe, &hRFlatRe, nNRTh, nNRPh, nNRFl}) h->Write();
    for (TH1* h : hRSlice) h->Write();
    for (TH1* h : chain1D) h->Write();
    for (int st = 0; st < kSt; st++){for (int a = 0; a < 3; a++) tEx[st][a]->Write(); apSt[st]->Write();}
    for (int k = 0; k < kEx; k++) for (int a = 0; a < 3; a++) hEx[k][a]->Write();
    for (int r = 0; r < 4; r++){rs[r].dThN->Write(); rs[r].dThL->Write(); rs[r].dPhN->Write(); rs[r].dPhL->Write();}
    for (int k = 0; k < 4; k++) for (int a = 0; a < 3; a++) hNew[k][a]->Write();
    for (TH1* h : hTSlice) h->Write();
    out.Close();
    printf(">>> SmearMatrix matrices, plots and csv took %.0f s\n", Secs()); fflush(stdout);
}
