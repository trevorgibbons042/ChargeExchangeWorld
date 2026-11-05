#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
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
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

#include "FitCore.h"
#include "CutOpen.h"

// Response matrices of the recon block angles (thetaReD, phiReD, CalcValsReD in cut_APandT.root) against the truth of the CE proton, for the random beam angle of the sim (/randGun/randPAngle: theta uniform
// 0 to 5 mrad = 0.29 deg, phi uniform). Two truths per CE proton: the lab truth (thetaCE, phiCE in CE: angles about the z axis, the frame the recon block angles are measured in) and the truth relative to the
// incoming neutron (the lab proton direction rotated so the neutron direction CENDirX/Y/Z of ChargeExchange_values, the last row of the track as in CalcAPandT.cpp, points along z).
// The lab truth shows the detector smearing alone, the neutron truth the detector smearing and the beam angle. Only the MC knows the neutron direction: the matrix carries the beam
// angle spread, real data needs only the matrix.
// Writes cuts/smear_matrix.root (R_theta, R_phi_all, R_phi_<k> per true theta slice, R_thetaphi flattened, their truth and recon histograms, so the matrices can go into RooUnfold) and
// plots/8_smear/: 1_R_theta, 2_R_phi, 3_R_thetaphi, 4_residual_theta, 5_residual_phi, 6_dilution_vs_theta .pdf and smear_dilution.csv. The fakes are in the root file only, not plotted.
// F_theta, F_thetaphi (the bins of M_thetaphi) and F_theta_phi are the fakes: the protons of AnglesRaw (cut_angles_all.root) that are not children of a carbon CE track (isFromCE == 0; silicon CE in a noPhys run), recon block
// angles only, no truth. For RooUnfoldResponse::Fake: the measured sample of a cut holds them next to the CE protons.
// Rows: every CE proton of CalcValsReD with a neutron direction. Cut: true (neutron frame) and recon block theta both 0 to kThetaMax for the matrices and residuals, recon block theta
// 0 to kThetaMax for the dilution; no other cut.

static constexpr double kStep = 0.1;     // deg, theta bin of R_theta and slice of the dilution
static constexpr int kPhiFine = 36;      // phi bins of R_phi_all and R_phi_<k>
static constexpr int kPhiFlat = 50;      // phi bins of R_thetaphi, same as the PhiCosFit plots of PlotAngle.cpp
static constexpr int kSlices = 10;      // true theta slices of R_phi_<k> and theta bins of R_thetaphi

static double ThetaMaxDeg(){return 5.0 * ThetaScale();} // same range as PlotPolarization.cpp
static double Wrap360(double a){a = std::fmod(a, 360.0); return a < 0 ? a + 360.0 : a;}
static double Wrap180(double a){return std::remainder(a, 360.0);}
// The flattened (theta slice, phi bin) index of R_thetaphi for theta in 0 to tMax and phi in 0 to 360 deg.
static int FlatBin(double theta, double phi, double tMax){
    return std::min(kSlices - 1, (int)(theta / (tMax / kSlices))) * kPhiFlat + std::min(kPhiFlat - 1, (int)(phi / (360.0 / kPhiFlat)));
}
struct NDir {double x, y, z;};
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
static void SavePlot2D(TH2D* h, const std::string& path, const std::string& title, double axMax = 0, const char* zTitle = "P(recon | true)"){
    TCanvas c("cSmear2D", "Plots", 1500, 1200);
    PageMargins(c);
    h->SetTitle(PlotTitle(path, title).c_str());
    h->SetStats(0);
    h->GetZaxis()->SetTitle(zTitle);
    h->GetZaxis()->SetTitleOffset(1.6); // the palette labels of a small maximum (0.25, 0.3) are wider than the default gap and run into the title
    if (axMax > 0){h->GetXaxis()->SetRangeUser(0, axMax); h->GetYaxis()->SetRangeUser(0, axMax);}
    h->Draw("COLZ");
    c.SaveAs(path.c_str());
}

// The matrix as proton counts, beside its probability plot. Drawn from a clone: SavePlot2D narrows the axis range of the histogram it draws, and the file keeps the full one.
static void SaveCounts(const TH2D& h, const std::string& path, const std::string& title, double axMax = 0){
    TH2D* c = (TH2D*)h.Clone(Form("%s_plot", h.GetName()));
    SavePlot2D(c, path, title, axMax, "Protons");
    delete c;
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
    leg.AddEntry(hNeutron, Form("#splitline{#splitline{recon block %s (lab frame)}{- CE %s (neutron frame)}}{RMS %.3f}", angle, angle, hNeutron->GetRMS()), "l");
    leg.AddEntry(hLab, Form("#splitline{#splitline{recon block %s (lab frame)}{- CE %s (lab frame)}}{RMS %.3f}", angle, angle, hLab->GetRMS()), "l");
    leg.Draw();
    c.SaveAs(path.c_str());
}

void SmearMatrix(){
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
        TTreeReader r(cev);
        TTreeReaderValue<int> ev(r, "fEvent"), tr(r, "TrackID");
        TTreeReaderValue<double> dx(r, "CENDirX"), dy(r, "CENDirY"), dz(r, "CENDirZ");
        while (r.Next()){
            dirs[EvKey(*ev, *tr)] = {*dx, *dy, *dz};
            maxTilt = std::max(maxTilt, std::acos(std::min(1.0, *dz)) * TMath::RadToDeg());
        }
    }

    for (const char* d : {"1_reconBlock", "2_recon", "3_beamAngle", "4_residuals"}) gSystem->mkdir(Form("plots/8_smear/%s", d), kTRUE);
    TFile out("cuts/smear_matrix.root", "RECREATE");
    const int nFlatT = kSlices * kPhiFlat;
    TH2D hRTheta("R_theta", "", nTheta, 0, tMax, nTheta, 0, tMax);
    TH1D hTTheta("T_theta", "", nTheta, 0, tMax), hMTheta("M_theta", "", nTheta, 0, tMax); // true theta of every row in range; recon theta of every row in range
    TH2D hRPhi("R_phi_all", "", kPhiFine, 0, 360, kPhiFine, 0, 360);
    TH1D hTPhi("T_phi_all", "", kPhiFine, 0, 360);
    TH2D hRThetaRe("R_theta_recon", "", nTheta, 0, tMax, nTheta, 0, tMax); // the same true bins (T_theta, T_phi_all) against the recon theta and phi
    TH2D hRPhiRe("R_phi_all_recon", "", kPhiFine, 0, 360, kPhiFine, 0, 360);
    // beam angle only, no detector: the true theta and phi of the proton about the neutron direction (x) against about the lab z axis (y, thetaCE and phiCE)
    TH2D hRBeamTheta("R_beam_theta", "", nTheta, 0, tMax, nTheta, 0, tMax), hRBeamPhi("R_beam_phi_all", "", kPhiFine, 0, 360, kPhiFine, 0, 360);
    // for UnfoldPhi.cpp: the flattened recon (not block) matrix and measured histogram, and the mean A_eff of a true theta slice (AEffOfTheta of the neutron theta, as the P fits of PlotAngle.cpp take it)
    TH2D hRFlatRe("R_thetaphi_recon", "", nFlatT, 0, nFlatT, nFlatT, 0, nFlatT);
    TH1D hMFlatRe("M_thetaphi_recon", "", nFlatT, 0, nFlatT);
    TProfile hAPTheta("AP_theta", "", kSlices, 0, tMax);
    const double apScale = ReadSettings().at("APScale");
    TH2D hRBeamFlat("R_beam_thetaphi", "", nFlatT, 0, nFlatT, nFlatT, 0, nFlatT); // the same two angles on the flattened (theta slice, phi bin) numbers of R_thetaphi
    TH2D hRFlat("R_thetaphi", "", nFlatT, 0, nFlatT, nFlatT, 0, nFlatT);
    TH1D hTFlat("T_thetaphi", "", nFlatT, 0, nFlatT), hMFlat("M_thetaphi", "", nFlatT, 0, nFlatT);
    std::vector<TH2D*> hRSlice;
    std::vector<TH1D*> hTSlice;
    for (int k = 0; k < kSlices; k++){
        hRSlice.push_back(new TH2D(Form("R_phi_%d", k), Form("true theta %.2f to %.2f deg", k * tMax / kSlices, (k + 1) * tMax / kSlices), kPhiFine, 0, 360, kPhiFine, 0, 360));
        hTSlice.push_back(new TH1D(Form("T_phi_%d", k), "", kPhiFine, 0, 360));
    }
    const double rt = 0.5 * ThetaScale(); // residual range, deg
    TH1D hDThN("dTheta_neutron", "", 100, -rt, rt), hDThL("dTheta_lab", "", 100, -rt, rt);
    TH1D hDPhN("dPhi_neutron", "", 90, -180, 180), hDPhL("dPhi_lab", "", 90, -180, 180);

    // dilution per recon block theta slice: sums of cos(dphi), cos(dphi)^2, count, for the neutron truth and the lab truth
    std::vector<double> sN(nTheta, 0), s2N(nTheta, 0), sL(nTheta, 0), s2L(nTheta, 0), nS(nTheta, 0);
    long nRows = 0, nNoDir = 0, nPair = 0;

    TTreeReader rc(ce), rr(red), rre(rec);
    TTreeReaderValue<int> ev(rc, "fEvent"), tr(rc, "TrackID"), ev2(rr, "fEvent"), tr2(rr, "TrackID");
    TTreeReaderValue<double> thCE(rc, "thetaCE"), phCE(rc, "phiCE"), thR(rr, "thetaReD"), phR(rr, "phiReD"), thRe(rre, "thetaRe"), phRe(rre, "phiRe");
    while (rc.Next() && rr.Next() && rre.Next()){
        nRows++;
        if (*ev != *ev2 || *tr != *tr2){printf("SmearMatrix: CE and CalcValsReD rows do not line up at row %ld\n", nRows); return;}
        auto it = dirs.find(EvKey(*ev, *tr));
        if (it == dirs.end()){nNoDir++; continue;}
        double thT, phT;
        RelToNeutron(*thCE, *phCE, it->second, thT, phT);
        const double thRec = *thR, phRec = Wrap360(*phR), phLab = Wrap360(*phCE);
        const bool tIn = thT >= 0 && thT < tMax, rIn = thRec >= 0 && thRec < tMax;
        const double thRc = *thRe, phRc = Wrap360(*phRe); // recon theta and phi
        if (thRc >= 0 && thRc < tMax) hMFlatRe.Fill(FlatBin(thRc, phRc, tMax) + 0.5);
        if (tIn && thRc >= 0 && thRc < tMax){hRThetaRe.Fill(thT, thRc); hRPhiRe.Fill(phT, phRc); hRFlatRe.Fill(FlatBin(thT, phT, tMax) + 0.5, FlatBin(thRc, phRc, tMax) + 0.5);}
        if (tIn && *thCE >= 0 && *thCE < tMax){hRBeamTheta.Fill(thT, *thCE); hRBeamPhi.Fill(phT, phLab); hRBeamFlat.Fill(FlatBin(thT, phT, tMax) + 0.5, FlatBin(*thCE, phLab, tMax) + 0.5);}
        const int kT = std::min(kSlices - 1, (int)(thT / (tMax / kSlices)));
        const int flatT = FlatBin(thT, phT, tMax), flatR = FlatBin(thRec, phRec, tMax);
        if (tIn){
            hTTheta.Fill(thT); hTPhi.Fill(phT); hTSlice[kT]->Fill(phT); hTFlat.Fill(flatT + 0.5); hAPTheta.Fill(thT, AEffOfTheta(thT, apScale));
        }
        if (rIn){
            hMTheta.Fill(thRec);
            hMFlat.Fill(flatR + 0.5);
            const double dN = Wrap180(phRec - phT), dL = Wrap180(phRec - phLab);
            const int b = std::min(nTheta - 1, (int)(thRec / kStep));
            nS[b] += 1; sN[b] += std::cos(dN * TMath::DegToRad()); s2N[b] += std::pow(std::cos(dN * TMath::DegToRad()), 2);
            sL[b] += std::cos(dL * TMath::DegToRad()); s2L[b] += std::pow(std::cos(dL * TMath::DegToRad()), 2);
            if (tIn){
                nPair++;
                hRTheta.Fill(thT, thRec); hRPhi.Fill(phT, phRec); hRSlice[kT]->Fill(phT, phRec); hRFlat.Fill(flatT + 0.5, flatR + 0.5);
                hDThN.Fill(thRec - thT); hDThL.Fill(thRec - *thCE);
                hDPhN.Fill(dN); hDPhL.Fill(dL);
            }
        }
    }
    printf("SmearMatrix: %ld CE rows, %ld without a neutron direction, %ld in the matrix\n", nRows, nNoDir, nPair);

    // fakes: the recon block protons that are not children of a carbon CE track (silicon CE in a noPhys run), recon theta and phi only
    TH1D hFTheta("F_theta", "", nTheta, 0, tMax), hFFlat("F_thetaphi", "", nFlatT, 0, nFlatT), hFFlatRe("F_thetaphi_recon", "", nFlatT, 0, nFlatT);
    TH2D hFThPh("F_theta_phi", "", nTheta, 0, tMax, kPhiFine, 0, 360);
    long nFake = 0;
    if (TTree* raw = OpenCut("cuts/cut_angles_all.root", "AnglesRaw")){
        TTreeReader rf(raw);
        TTreeReaderValue<int> pdg(rf, "PDG"), fromCE(rf, "isFromCE");
        TTreeReaderValue<double> thF(rf, "outDTheta"), phF(rf, "outDPhi"), thF2(rf, "outTheta"), phF2(rf, "outPhi");
        while (rf.Next()){
            if (*pdg != 2212 || *fromCE) continue;
            const double th = *thF * TMath::RadToDeg(), ph = Wrap360(*phF * TMath::RadToDeg());
            const double thr = *thF2 * TMath::RadToDeg(), phr = Wrap360(*phF2 * TMath::RadToDeg()); // recon theta and phi (the recon matrix of UnfoldPhi.cpp)
            if (thr >= 0 && thr < tMax) hFFlatRe.Fill(FlatBin(thr, phr, tMax) + 0.5);
            if (th < 0 || th >= tMax) continue;
            hFTheta.Fill(th); hFFlat.Fill(FlatBin(th, ph, tMax) + 0.5); hFThPh.Fill(th, ph);
            nFake++;
        }
    } else printf("SmearMatrix: cut_angles_all.root AnglesRaw missing, no fakes\n");
    printf("SmearMatrix: %ld fake protons (silicon CE, not from a carbon CE track) with recon block theta 0 to %.1f deg\n", nFake, tMax);

    // dilution csv: D = mean cos(recon block phi - true phi) per recon block theta slice
    FILE* csv = fopen("plots/8_smear/4_residuals/smear_dilution.csv", "w");
    fprintf(csv, "theta_lo,theta_hi,n,D_neutron,D_neutron_err,D_lab,D_lab_err\n");
    TGraphErrors gN, gL;
    for (int b = 0; b < nTheta; b++){
        if (nS[b] < 2) continue;
        const double n = nS[b];
        const double dN = sN[b] / n, dL = sL[b] / n;
        const double eN = std::sqrt(std::max(0.0, s2N[b] / n - dN * dN) / n), eL = std::sqrt(std::max(0.0, s2L[b] / n - dL * dL) / n);
        fprintf(csv, "%.4f,%.4f,%.0f,%.6f,%.6f,%.6f,%.6f\n", b * kStep, (b + 1) * kStep, n, dN, eN, dL, eL);
        const double x = (b + 0.5) * kStep;
        gN.SetPoint(gN.GetN(), x, dN); gN.SetPointError(gN.GetN() - 1, 0.5 * kStep, eN);
        gL.SetPoint(gL.GetN(), x, dL); gL.SetPointError(gL.GetN() - 1, 0.5 * kStep, eL);
    }
    fclose(csv);

    // column-normalised copies are what is plotted, the root file keeps the counts and the truth/recon histograms
    TH2D* nTheta2 = (TH2D*)hRTheta.Clone("R_theta_norm");
    TH2D* nFlat2 = (TH2D*)hRFlat.Clone("R_thetaphi_norm");
    TH2D* nPhi2 = (TH2D*)hRPhi.Clone("R_phi_all_norm");
    NormColumns(nTheta2, &hTTheta);
    NormColumns(nFlat2, &hTFlat);
    NormColumns(nPhi2, &hTPhi);
    TH2D* nThRe = (TH2D*)hRThetaRe.Clone("R_theta_recon_norm");
    TH2D* nPhRe = (TH2D*)hRPhiRe.Clone("R_phi_all_recon_norm");
    NormColumns(nThRe, &hTTheta);
    NormColumns(nPhRe, &hTPhi);
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
    SavePlot2D(nTheta2, "plots/8_smear/1_reconBlock/1_R_theta.pdf", Form("Response of recon block theta (lab frame), %s, %s;CE theta, neutron frame (deg);Recon block theta, lab frame (deg)", sd.c_str(), cutPair.c_str()), ThetaAxisMax());
    SavePlot2D(nPhi2, "plots/8_smear/1_reconBlock/2_R_phi.pdf", Form("Response of recon block phi (lab frame), %s, %s;CE phi, neutron frame (deg);Recon block phi, lab frame (deg)", sd.c_str(), cutPair.c_str()));
    SavePlot2D(nFlat2, "plots/8_smear/1_reconBlock/3_R_thetaphi.pdf", Form("Response of recon block (theta, phi) (lab frame), %s, %s, bin number = theta slice x %d + phi bin (%d theta slices of %.2f deg, %d phi bins of %.1f deg);CE (theta, phi) bin number, neutron frame;Recon block (theta, phi) bin number, lab frame", sd.c_str(), cutPair.c_str(), kPhiFlat, kSlices, tMax / kSlices, kPhiFlat, 360.0 / kPhiFlat));
    const std::string cutPairRe = Form("CE recon protons, both thetas 0 to %.1f deg, no other cuts", tMax);
    SavePlot2D(nThRe, "plots/8_smear/2_recon/8_R_theta_recon.pdf", Form("Response of recon theta (lab frame), %s, %s;CE theta, neutron frame (deg);Recon theta, lab frame (deg)", sd.c_str(), cutPairRe.c_str()), ThetaAxisMax());
    SavePlot2D(nPhRe, "plots/8_smear/2_recon/9_R_phi_recon.pdf", Form("Response of recon phi (lab frame), %s, %s;CE phi, neutron frame (deg);Recon phi, lab frame (deg)", sd.c_str(), cutPairRe.c_str()));
    // the same four matrices as proton counts per bin (x_2_ next to x_), no column left blank
    SaveCounts(hRTheta, "plots/8_smear/1_reconBlock/1_2_R_theta_counts.pdf", Form("Protons per bin, recon block theta (lab frame), %s, %s;CE theta, neutron frame (deg);Recon block theta, lab frame (deg)", sd.c_str(), cutPair.c_str()), ThetaAxisMax());
    SaveCounts(hRPhi, "plots/8_smear/1_reconBlock/2_2_R_phi_counts.pdf", Form("Protons per bin, recon block phi (lab frame), %s, %s;CE phi, neutron frame (deg);Recon block phi, lab frame (deg)", sd.c_str(), cutPair.c_str()));
    SaveCounts(hRThetaRe, "plots/8_smear/2_recon/8_2_R_theta_recon_counts.pdf", Form("Protons per bin, recon theta (lab frame), %s, %s;CE theta, neutron frame (deg);Recon theta, lab frame (deg)", sd.c_str(), cutPairRe.c_str()), ThetaAxisMax());
    SaveCounts(hRPhiRe, "plots/8_smear/2_recon/9_2_R_phi_recon_counts.pdf", Form("Protons per bin, recon phi (lab frame), %s, %s;CE phi, neutron frame (deg);Recon phi, lab frame (deg)", sd.c_str(), cutPairRe.c_str()));
    // beam angle only (no detector): the lab angle the proton has about z against its angle about the neutron
    const std::string cutBeam = Form("CE protons, both thetas 0 to %.1f deg, no other cuts", tMax);
    SavePlot2D(nBTh, "plots/8_smear/3_beamAngle/10_R_beam_theta.pdf", Form("Response of the beam angle, CE theta (lab frame), %s, %s;CE theta, neutron frame (deg);CE theta, lab frame (deg)", sd.c_str(), cutBeam.c_str()), ThetaAxisMax());
    SaveCounts(hRBeamTheta, "plots/8_smear/3_beamAngle/10_2_R_beam_theta_counts.pdf", Form("Protons per bin, beam angle, CE theta (lab frame), %s, %s;CE theta, neutron frame (deg);CE theta, lab frame (deg)", sd.c_str(), cutBeam.c_str()), ThetaAxisMax());
    SavePlot2D(nBPh, "plots/8_smear/3_beamAngle/11_R_beam_phi.pdf", Form("Response of the beam angle, CE phi (lab frame), %s, %s;CE phi, neutron frame (deg);CE phi, lab frame (deg)", sd.c_str(), cutBeam.c_str()));
    SaveCounts(hRBeamPhi, "plots/8_smear/3_beamAngle/11_2_R_beam_phi_counts.pdf", Form("Protons per bin, beam angle, CE phi (lab frame), %s, %s;CE phi, neutron frame (deg);CE phi, lab frame (deg)", sd.c_str(), cutBeam.c_str()));
    const std::string flatDef = Form("bin number = theta slice x %d + phi bin (%d theta slices of %.2f deg, %d phi bins of %.1f deg)", kPhiFlat, kSlices, tMax / kSlices, kPhiFlat, 360.0 / kPhiFlat);
    SavePlot2D(nBFl, "plots/8_smear/3_beamAngle/12_R_beam_thetaphi.pdf", Form("Response of the beam angle, CE (theta, phi) (lab frame), %s, %s, %s;CE (theta, phi) bin number, neutron frame;CE (theta, phi) bin number, lab frame", sd.c_str(), cutBeam.c_str(), flatDef.c_str()));
    SaveCounts(hRBeamFlat, "plots/8_smear/3_beamAngle/12_2_R_beam_thetaphi_counts.pdf", Form("Protons per bin, beam angle, CE (theta, phi) (lab frame), %s, %s, %s;CE (theta, phi) bin number, neutron frame;CE (theta, phi) bin number, lab frame", sd.c_str(), cutBeam.c_str(), flatDef.c_str()));
    SaveResidual(&hDThN, &hDThL, "theta", "plots/8_smear/4_residuals/4_residual_theta.pdf", Form("Recon block theta (lab frame) - CE theta (neutron frame or lab frame), %s, %s;#Delta#theta (deg);Count", sd.c_str(), cutPair.c_str()));
    SaveResidual(&hDPhN, &hDPhL, "phi", "plots/8_smear/4_residuals/5_residual_phi.pdf", Form("Recon block phi (lab frame) - CE phi (neutron frame or lab frame), %s, %s;#Delta#phi (deg);Count", sd.c_str(), cutPair.c_str()));

    {
        const std::string path = "plots/8_smear/4_residuals/6_dilution_vs_theta.pdf";
        TCanvas c("cSmearD", "Plots", 1500, 1200);
        PageMargins(c);
        TMultiGraph mg;
        gN.SetMarkerStyle(20); gN.SetMarkerSize(1.2); gN.SetMarkerColor(kRed + 1); gN.SetLineColor(kRed + 1);
        gL.SetMarkerStyle(20); gL.SetMarkerSize(1.2); gL.SetMarkerColor(kBlue + 1); gL.SetLineColor(kBlue + 1);
        mg.Add(&gN, "PE"); mg.Add(&gL, "PE");
        mg.SetTitle(PlotTitle(path, Form("Dilution D = <cos(#phi_{recon block, lab frame} - #phi_{CE, neutron or lab frame})>, %s, CE recon protons, recon block theta 0 to %.1f deg, no other cuts;Recon block theta, lab frame (deg);D", sd.c_str(), tMax)).c_str());
        mg.Draw("A");
        mg.GetYaxis()->SetRangeUser(-0.1, 1.1);
        mg.GetXaxis()->SetLimits(0, ThetaAxisMax());
        TLegend leg(0.715, 0.76, 0.99, 0.88);
        leg.SetTextSize(0.017);
    leg.AddEntry(&gN, "#splitline{CE phi (neutron frame)}{detector + beam angle}", "pe");
        leg.AddEntry(&gL, "#splitline{CE phi (lab frame)}{detector only}", "pe");
        leg.Draw();
        c.SaveAs(path.c_str());
    }

    out.cd();
    for (TH1* h : std::vector<TH1*>{&hRTheta, &hTTheta, &hMTheta, &hRPhi, &hTPhi, &hRFlat, &hTFlat, &hMFlat, nTheta2, nFlat2, nPhi2, &hDThN, &hDThL, &hDPhN, &hDPhL, &hFTheta, &hFFlat, &hFThPh, &hRThetaRe, &hRPhiRe, nThRe, nPhRe, &hRBeamTheta, &hRBeamPhi, nBTh, nBPh, &hRBeamFlat, nBFl, &hRFlatRe, &hMFlatRe, &hFFlatRe, &hAPTheta}) h->Write();
    for (TH1* h : hRSlice) h->Write();
    for (TH1* h : hTSlice) h->Write();
    out.Close();
}
