#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TF1.h>
#include <TNamed.h>
#include <ROOT/TThreadExecutor.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "FitCore.h"
#include "CutOpen.h"
#include "TreeColumns.h"

// Iterative Bayesian (D'Agostini) unfolding of the recon block phi of the CE protons from the estimated neutron frame back to the true neutron frame with the frame switch matrix M_recB_Fswitch
// (R_reconB_neuF_to_reconB_labF_thetaphi, truth T_reconB_neuF_thetaphi; the matrix names say "labF" for the estimated frame), flat prior, kIterDraw iterations, then the P fit of PlotAngle.cpp on the phi:
// N(phi) = N0 (1 + P A_eff cos(phi)), A fixed. Reads cuts/smear_matrix.root (SmearMatrix.cpp), run SmearMatrix.cpp first, from a run folder.
// Writes plots/8_smear/3_phi_cosfit/1_M_recB_phi_cosfit_unfolded.pdf (the estimated frame phi unfolded to the true neutron frame) and 2_M_recB_phi_cosfit_trueNeutron.pdf (the true neutron frame truth, the answer the unfolding aims at), same y axis.
// Each response is P(recon bin j | true bin i) = R_ij / T_i with T the truth of every row (so a column sums to the efficiency, the misses are the rest). The bins are the flattened
// (theta slice x phi bin) numbers of SmearMatrix.cpp. The error of the unfolded bins and of P comes from kToys Poisson toys of the measured histogram (the statistical error of the
// matrix is not in it). A_eff of a curve is sum_i AP(theta_i) / N over its protons: the mean AP of each theta slice in the frame of the curve (AP_reconB_neuF for the truth and the unfolded curve, AP_reconB_labF for the measured estimated frame curve, SmearMatrix.cpp), weighted by the counts of the slice.

static constexpr int kPhiBins = 50;  // phi bins of the flattened matrix (SmearMatrix.cpp kPhiFlat)
static constexpr int kSlices = 10;   // theta slices of the flattened matrix (SmearMatrix.cpp kSlices)
static constexpr int kBins = kSlices * kPhiBins;
static constexpr int kIterMax = 8;   // iterations kept for the P vs iteration plot
static constexpr int kIterDraw = 4;  // iterations of the drawn unfolded histogram
static constexpr int kToys = 200;

static void PageMargins(TCanvas& c){c.SetRightMargin(0.30); c.SetLeftMargin(0.14); c.SetBottomMargin(0.13);}

struct PResult {double p = 0, err = 0, a = 0; bool ok = false;};

// Weighted least squares of the bin counts n[b] = a + b c_b, c_b the mean of cos(phi) over the bin (equal to the fit N0 (1 + P A cos(phi)) with b = a P A), so P = b / (a A_eff).
// sig: the error of each bin. The error of P is the delta method with the covariance of a and b.
static PResult FitP(const double* n, const double* sig, double aEff){
    PResult r;
    double s = 0, sc = 0, scc = 0, sn = 0, snc = 0;
    const double w = 2 * TMath::Pi() / kPhiBins;
    for (int b = 0; b < kPhiBins; b++){
        const double lo = b * w, hi = lo + w, c = (std::sin(hi) - std::sin(lo)) / w, wt = 1.0 / (sig[b] * sig[b]);
        s += wt; sc += wt * c; scc += wt * c * c; sn += wt * n[b]; snc += wt * n[b] * c;
    }
    const double det = s * scc - sc * sc;
    if (!(det > 0) || aEff == 0) return r;
    const double a = (scc * sn - sc * snc) / det, bb = (s * snc - sc * sn) / det;
    const double vaa = scc / det, vbb = s / det, vab = -sc / det;
    if (!(a > 0)) return r;
    const double dPb = 1.0 / (a * aEff), dPa = -bb / (a * a * aEff);
    r.p = bb / (a * aEff); r.a = a;
    r.err = std::sqrt(std::max(0.0, dPb * dPb * vbb + dPa * dPa * vaa + 2 * dPb * dPa * vab));
    r.ok = true;
    return r;
}

// the flattened bin index is theta slice x kPhiBins + phi bin: the phi histogram sums the slices, sliceN is the count of each slice
static void PhiOf(const std::vector<double>& v, double* phi, double* sliceN){
    for (int b = 0; b < kPhiBins; b++) phi[b] = 0;
    for (int k = 0; k < kSlices; k++) sliceN[k] = 0;
    for (int i = 0; i < kBins; i++){phi[i % kPhiBins] += v[i]; sliceN[i / kPhiBins] += v[i];}
}

static double AEffOf(const double* sliceN, const TProfile& ap){
    double sum = 0, tot = 0;
    for (int k = 0; k < kSlices; k++){sum += sliceN[k] * ap.GetBinContent(k + 1); tot += sliceN[k];}
    return tot > 0 ? sum / tot : 0.0;
}

struct Response {std::vector<std::vector<double>> p; std::vector<double> eff;};

static Response MakeResponse(const TH2D& R, const TH1D& T){
    Response r;
    r.p.assign(kBins, std::vector<double>(kBins, 0.0));
    r.eff.assign(kBins, 0.0);
    for (int i = 0; i < kBins; i++){
        const double t = T.GetBinContent(i + 1);
        if (t <= 0) continue;
        for (int j = 0; j < kBins; j++){r.p[i][j] = R.GetBinContent(i + 1, j + 1) / t; r.eff[i] += r.p[i][j];}
    }
    return r;
}

// D'Agostini: t_i <- t_i / eff_i * sum_j p_ij m_j / f_j with f_j = sum_k p_kj t_k, from a flat prior over the bins that have a response. Returns the truth after each of nIter iterations.
static std::vector<std::vector<double>> Bayes(const Response& r, const std::vector<double>& m, int nIter){
    std::vector<double> t(kBins);
    for (int i = 0; i < kBins; i++) t[i] = r.eff[i] > 0 ? 1.0 : 0.0;
    std::vector<std::vector<double>> out;
    for (int it = 0; it < nIter; it++){
        std::vector<double> f(kBins, 0.0), next(kBins, 0.0);
        for (int k = 0; k < kBins; k++) if (t[k] > 0) for (int j = 0; j < kBins; j++) f[j] += r.p[k][j] * t[k];
        for (int i = 0; i < kBins; i++){
            if (!(r.eff[i] > 0) || !(t[i] > 0)) continue;
            double s = 0;
            for (int j = 0; j < kBins; j++) if (f[j] > 0) s += r.p[i][j] * m[j] / f[j];
            next[i] = t[i] / r.eff[i] * s;
        }
        t = next;
        out.push_back(t);
    }
    return out;
}

static std::vector<double> ContentsOf(const TH1D& h){
    std::vector<double> v(kBins);
    for (int i = 0; i < kBins; i++) v[i] = h.GetBinContent(i + 1);
    return v;
}

static TH1D* PhiHist(const char* name, const double* v, const double* sig, int color, int width, int marker){
    TH1D* h = new TH1D(name, "", kPhiBins, 0, 360);
    for (int b = 0; b < kPhiBins; b++){h->SetBinContent(b + 1, v[b]); h->SetBinError(b + 1, sig ? sig[b] : std::sqrt(std::max(v[b], 1.0)));}
    h->SetLineColor(color); h->SetLineWidth(width); h->SetMarkerColor(color);
    if (marker){h->SetMarkerStyle(marker); h->SetMarkerSize(1.2);}
    return h;
}

struct Curve {std::string file, name; const double* v; const double* sig; PResult r; double a; int color; bool points;};

// one curve page: the phi histogram with its fit and P, the look of the PhiCosFit plots of PlotAngle.cpp; id makes the ROOT object names unique
static void DrawCurve(const Curve& cv, const std::string& path, const std::string& id, const std::string& cuts, const std::string& extra, double top, double tMax){
    TCanvas c("cUnfold", "Plots", 1500, 1200); // the look of the PhiCosFit plots of PlotAngle.cpp
    c.SetLeftMargin(0.14); c.SetBottomMargin(0.13); c.SetRightMargin(0.25); c.SetTopMargin(0.17);
    TH1D* h = PhiHist(("h" + id).c_str(), cv.v, cv.sig, kRed, 1, cv.points ? 20 : 0);
    double nSum = 0; for (int b = 0; b < kPhiBins; b++) nSum += cv.v[b];
    h->SetEntries(nSum); h->SetMarkerSize(0.8);
    h->SetTitle(PlotTitle(path, Form("Phi, %s, %s%s, %d phi bins of %.1f deg;Angle (degrees);Count", cv.name.c_str(), cuts.c_str(), extra.c_str(), kPhiBins, 360.0 / kPhiBins)).c_str());
    h->SetMinimum(0); h->SetMaximum(1.1 * top);
    TF1* f = new TF1(("f" + id).c_str(), "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
    f->SetParNames("N0 (fit)", "P", "A"); f->SetParameters(cv.r.a, cv.r.p, cv.a); f->SetParError(1, cv.r.err);
    f->SetLineColor(kBlack); f->SetNDF(0); f->SetChisquare(0);
    h->GetListOfFunctions()->Add(f);
    h->SetName("Fit results"); // the first line of the stats box is the histogram name
    h->Draw(cv.points ? "E1" : "HIST"); // plain line like PhiCosFit, only the unfolded curve keeps its error bars
    f->Draw("SAME");
    c.Update();
    if (auto* st = (TPaveStats*)h->FindObject("stats")){st->SetX1NDC(0.76); st->SetX2NDC(0.995); st->SetY1NDC(0.55); st->SetY2NDC(0.83);}
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
    tx.DrawLatex(0.13, 0.90, "Binned fit  N(#phi) #propto N_{0} (1 + P #upoint A_{eff} cos#phi) (weighted least squares on the phi bins)");
    tx.DrawLatex(0.13, 0.865, Form("%s = %.4f (fixed, %.2f deg theta slices),   P = %.2f #pm %.2f", AEffEq(), cv.a, tMax / kSlices, cv.r.p, cv.r.err));
    c.Modified(); c.Update();
    c.SaveAs(path.c_str());
}

// Unfolds the measured counts meas (flattened (theta slice, phi bin), estimated neutron frame) for the first nIter iterations: the phi histogram of each iteration, its error (the spread of kToys Poisson toys of
// meas; the statistical error of the matrix is not in it), the A_eff of that curve (apTrue, the AP profile of the true neutron frame) and its P fit (the error of P from the toys). The kToys toys are
// independent: they run on all the threads (ANA_NTHREADS, NThreads()), each toy with its own seeded generator, so the result does not depend on the thread count.
static void UnfoldIterations(const Response& resp, const std::vector<double>& meas, const TProfile& apTrue, int nIter, double (*phiU)[kPhiBins], double (*sigU)[kPhiBins], double* aU, PResult* pU){
    const auto nominal = Bayes(resp, meas, nIter);
    std::vector<double> toyPhi((size_t)nIter * kToys * kPhiBins);
    ROOT::TThreadExecutor pool(NThreads());
    pool.Foreach([&](int t){
        TRandom3 rng(12345 + t);
        std::vector<double> mt(kBins);
        for (int j = 0; j < kBins; j++) mt[j] = rng.Poisson(meas[j]);
        const auto u = Bayes(resp, mt, nIter);
        double sl[kSlices];
        for (int it = 0; it < nIter; it++) PhiOf(u[it], &toyPhi[((size_t)it * kToys + t) * kPhiBins], sl);
    }, ROOT::TSeqI(kToys));
    double sliceN[kSlices];
    for (int it = 0; it < nIter; it++){
        PhiOf(nominal[it], phiU[it], sliceN);
        aU[it] = AEffOf(sliceN, apTrue);
        for (int b = 0; b < kPhiBins; b++){
            double s = 0, s2 = 0;
            for (int t = 0; t < kToys; t++){const double x = toyPhi[((size_t)it * kToys + t) * kPhiBins + b]; s += x; s2 += x * x;}
            const double mean = s / kToys;
            sigU[it][b] = std::max(std::sqrt(std::max(0.0, s2 / kToys - mean * mean)), 1e-6);
        }
        pU[it] = FitP(phiU[it], sigU[it], aU[it]);
        double sp = 0, sp2 = 0; int ok = 0;
        for (int t = 0; t < kToys; t++){
            const PResult r = FitP(&toyPhi[((size_t)it * kToys + t) * kPhiBins], sigU[it], aU[it]);
            if (!r.ok) continue;
            sp += r.p; sp2 += r.p * r.p; ok++;
        }
        if (ok > 1){const double mean = sp / ok; pU[it].err = std::sqrt(std::max(0.0, sp2 / ok - mean * mean));}
    }
}

// Frame switch unfolding: the recon block phi counts of the estimated neutron frame (TDet) back to the true neutron frame (T) with the frame switch matrix R, then the two pages, the unfolded curve and the truth curve
static void UnfoldFswitch(const TH2D& RBeam, const TH1D& T, const TH1D& TDet, const TProfile& apTrue, const TProfile& apEst, double tMax, const std::string& pathUnfolded, const std::string& pathTruth){
    const Response resp = MakeResponse(RBeam, T);
    const std::vector<double> lab = ContentsOf(TDet), truthNeutron = ContentsOf(T);
    double phiU[kIterMax][kPhiBins], sigU[kIterMax][kPhiBins], aU[kIterMax];
    PResult pU[kIterMax];
    UnfoldIterations(resp, lab, apTrue, kIterMax, phiU, sigU, aU, pU);
    // the y axis of the two pages is the same
    double phiL[kPhiBins], phiN[kPhiBins], sigN[kPhiBins], slL[kSlices], slN[kSlices], top = 0;
    PhiOf(lab, phiL, slL); PhiOf(truthNeutron, phiN, slN);
    for (int b = 0; b < kPhiBins; b++) sigN[b] = std::sqrt(std::max(phiN[b], 1.0));
    const double aN = AEffOf(slN, apTrue);
    const PResult pN = FitP(phiN, sigN, aN);
    double sigL[kPhiBins]; // the estimated neutron frame phi as measured, no unfolding
    for (int b = 0; b < kPhiBins; b++) sigL[b] = std::sqrt(std::max(phiL[b], 1.0));
    const double aL = AEffOf(slL, apEst);
    const PResult pL = FitP(phiL, sigL, aL);
    const int d = kIterDraw - 1;
    for (int b = 0; b < kPhiBins; b++) top = std::max({top, phiN[b] + sigN[b], phiL[b] + sigL[b], phiU[d][b] + sigU[d][b]});
    printf("UnfoldPhi recon block phi: P unfolded from the estimated neutron frame to the true neutron frame (%d it) %.3f +- %.3f, true neutron frame truth %.3f +- %.3f\n", kIterDraw, pU[d].p, pU[d].err, pN.p, pN.err);
    const std::string cuts = Form("CE protons, recon block theta 0 to %.1f deg, no other cuts", tMax);
    const Curve cvU = {"unfolded", "recon block phi, estimated neutron frame unfolded to the true neutron frame", phiU[d], sigU[d], pU[d], aU[d], kBlack, true};
    DrawCurve(cvU, pathUnfolded, "recBUnfolded", cuts, Form(", unfolded with M_recB_Fswitch only, %d iterations, flat prior", kIterDraw), top, tMax);
    const Curve cvT = {"truth", "recon block phi, true neutron frame", phiN, sigN, pN, aN, kGreen + 2, false};
    DrawCurve(cvT, pathTruth, "recBTruth", cuts, "", top, tMax);
    const Curve cvE = {"measured", "recon block phi, estimated neutron frame (no unfolding)", phiL, sigL, pL, aL, kBlue + 1, false};
    DrawCurve(cvE, "plots/8_smear/3_phi_cosfit/3_M_recB_phi_cosfit_estNeutron.pdf", "recBEst", cuts, "", top, tMax);
    {   // page 9: P of the unfolded phi against the number of times the matrix was applied, with the true neutron frame truth P as a line
        const std::string path = "plots/8_smear/3_phi_cosfit/4_M_recB_P_vs_iteration.pdf";
        TCanvas c("cUnfoldIt", "Plots", 1500, 1200);
        PageMargins(c);
        TGraphErrors g;
        for (int it = 0; it < kIterMax; it++){g.SetPoint(it, it + 1, pU[it].p); g.SetPointError(it, 0, pU[it].err);}
        g.SetMarkerStyle(20); g.SetMarkerSize(1.2);
        g.SetTitle(PlotTitle(path, Form("P of the unfolded recon block phi (estimated neutron frame to true neutron frame, M_recB_Fswitch) against the times the matrix was applied, flat prior, %s, %d phi bins of %.1f deg;Times the matrix was applied (iterations);P", cuts.c_str(), kPhiBins, 360.0 / kPhiBins)).c_str());
        g.Draw("AP");
        g.GetXaxis()->SetLimits(0.5, kIterMax + 0.5);
        const double pSim = ReadSettings().at("PPValue"); // the polarization the sim was run with (settings.csv)
        TLine truth(0.5, pSim, kIterMax + 0.5, pSim); truth.SetLineColor(kGreen + 2); truth.SetLineWidth(2); truth.SetLineStyle(2); truth.Draw();
        TLegend leg(0.715, 0.76, 0.99, 0.88);
        leg.SetTextSize(0.017);
        leg.AddEntry(&g, "unfolded phi", "pe");
        leg.AddEntry(&truth, Form("true polarization (sim input), P = %.2f", pSim), "l");
        leg.Draw();
        c.SaveAs(path.c_str());
    }
}

// The measured counts of a sample: the recon block angles in the estimated neutron frame (thetaRelND, phiRelND of AnglesRaw: RECONTHETABL_ESTNEU, RECONPHIBL_ESTNEU) in the flattened (theta slice, phi bin)
// numbers of the matrix (SmearMatrix.cpp FlatBin), the rows that pass filter. Read on all threads, the order does not matter.
static std::vector<double> MeasuredOf(TTree* t, const std::string& filter, double tMax){
    std::vector<double> m(kBins, 0.0);
    const auto c = ColumnsUnordered(t, {"thetaRelND * 57.29577951308232", "phiRelND * 57.29577951308232"}, filter);
    for (size_t i = 0; i < c[0].size(); i++){
        const double th = c[0][i], ph = std::fmod(std::fmod(c[1][i], 360.0) + 360.0, 360.0);
        if (!(th >= 0 && th < tMax)) continue; // NaN (no tube) fails
        m[std::min(kSlices - 1, (int)(th / (tMax / kSlices))) * kPhiBins + std::min(kPhiBins - 1, (int)(ph / (360.0 / kPhiBins)))] += 1;
    }
    return m;
}

// The unfolding of 1_phiReD_cosfit.pdf (all recon protons) and 1_phiReD_cosfit_CE.pdf (the CE recon protons) of 6_cuts/9_megaCut3/1_1particlecut, each into its own pdf beside them: the recon block angles in the
// estimated neutron frame of cut_megacut3_1.root unfolded to the true neutron frame with the matrix of that cut (M_recB_Fswitch in cuts/smear_matrix_9_megaCut3_1_1particlecut.root: SmearMatrix.cpp with SMEAR_CUT,
// smear_cuts.sh), kIterDraw iterations, flat prior. The block phi spike at 0 stays in the measured counts (it is in the matrix too). The matrix describes CE protons only, so the all recon protons page
// says the non-CE protons are not in it.
static void UnfoldMegaCut3(double tMax){
    const char* matrixFile = "cuts/smear_matrix_9_megaCut3_1_1particlecut.root";
    TFile* file = TFile::Open(matrixFile, "READ");
    if (!file || file->IsZombie()){printf("UnfoldPhi: %s missing, no megaCut3 unfolding (run smear_cuts.sh: SmearMatrix.cpp with SMEAR_CUT=9_megaCut3/1_1particlecut)\n", matrixFile); return;}
    const TH2D* r = (TH2D*)file->Get("R_reconB_neuF_to_reconB_labF_thetaphi"); const TH1D* tTrue = (TH1D*)file->Get("T_reconB_neuF_thetaphi"); const TProfile* apTrue = (TProfile*)file->Get("AP_reconB_neuF");
    const TNamed* cutText = (TNamed*)file->Get("cut_text");
    if (!r || !tTrue || !apTrue || tTrue->GetNbinsX() != kBins){printf("UnfoldPhi: a matrix is missing or has not %d bins in %s (rerun SmearMatrix.cpp with SMEAR_CUT)\n", kBins, matrixFile); return;}
    TTree* raw = OpenCut("cuts/cut_megacut3_1.root", "AnglesRaw");
    if (!raw) return;
    const std::string dir = "plots/6_cuts/9_megaCut3/1_1particlecut/2_phi_cosfit/";
    gSystem->mkdir(dir.c_str(), kTRUE);
    const Response resp = MakeResponse(*r, *tTrue);
    const std::string cuts = cutText ? cutText->GetTitle() : "9_megaCut3/1_1particlecut";
    struct Sample {const char *filter, *pdf, *name, *note;};
    const Sample samples[2] = {{"PDG == 2212", "1_phiReD_cosfit_unfolded.pdf", "all recon protons", ", the matrix describes CE protons only: the non-CE protons are not in it"},
                               {"PDG == 2212 && isFromCE", "1_phiReD_cosfit_CE_unfolded.pdf", "CE recon protons", ""}};
    for (int s = 0; s < 2; s++){
        const std::vector<double> meas = MeasuredOf(raw, samples[s].filter, tMax);
        double n = 0; for (double x : meas) n += x;
        if (!(n > 0)){printf("UnfoldPhi: no %s in cut_megacut3_1.root with a recon block theta 0 to %.1f deg\n", samples[s].name, tMax); continue;}
        double phiU[kIterMax][kPhiBins], sigU[kIterMax][kPhiBins], aU[kIterMax];
        PResult pU[kIterMax];
        UnfoldIterations(resp, meas, *apTrue, kIterDraw, phiU, sigU, aU, pU); // only the iterations that are drawn
        const int d = kIterDraw - 1;
        double top = 0;
        for (int b = 0; b < kPhiBins; b++) top = std::max(top, phiU[d][b] + sigU[d][b]);
        printf("UnfoldPhi megaCut3 %s (%.0f protons): P unfolded to the true neutron frame (%d it) %.3f +- %.3f\n", samples[s].name, n, kIterDraw, pU[d].p, pU[d].err);
        const Curve cv = {"unfolded", std::string("recon block phi, ") + samples[s].name + ", estimated neutron frame unfolded to the true neutron frame", phiU[d], sigU[d], pU[d], aU[d], kBlack, true};
        DrawCurve(cv, dir + samples[s].pdf, Form("mega%d", s), cuts, Form(", unfolded with M_recB_Fswitch of this cut, %d iterations, flat prior, block phi spike at 0 kept%s", kIterDraw, samples[s].note), top, tMax);
    }
}

void UnfoldPhi(){
    FastJit();
    TH1::AddDirectory(kFALSE);
    PlainFonts(); gStyle->SetOptStat(1111); gStyle->SetOptFit(11);
    UnfoldMegaCut3(5.0 * ThetaScale()); // the two megaCut3 phi cosfits unfolded with the matrix of smear_cuts.sh; first, it does not need the full sample matrix below
    TFile* file = TFile::Open("cuts/smear_matrix.root", "READ");
    if (!file || file->IsZombie()){printf("UnfoldPhi: cuts/smear_matrix.root missing (run SmearMatrix.cpp)\n"); return;}
    const TH2D* r = (TH2D*)file->Get("R_reconB_neuF_to_reconB_labF_thetaphi"); const TH1D* tTrue = (TH1D*)file->Get("T_reconB_neuF_thetaphi"); const TH1D* tEst = (TH1D*)file->Get("T_reconB_labF_thetaphi");
    const TProfile* apTrue = (TProfile*)file->Get("AP_reconB_neuF"); const TProfile* apEst = (TProfile*)file->Get("AP_reconB_labF");
    if (!r || !tTrue || !tEst || !apTrue || !apEst){printf("UnfoldPhi: a matrix is missing in smear_matrix.root (rerun SmearMatrix.cpp)\n"); return;}
    if (tTrue->GetNbinsX() != kBins){printf("UnfoldPhi: the matrix has %d bins, the macro expects %d\n", tTrue->GetNbinsX(), kBins); return;}
    gSystem->mkdir("plots/8_smear/3_phi_cosfit", kTRUE);
    UnfoldFswitch(*r, *tTrue, *tEst, *apTrue, *apEst, 5.0 * ThetaScale(), "plots/8_smear/3_phi_cosfit/1_M_recB_phi_cosfit_unfolded.pdf", "plots/8_smear/3_phi_cosfit/2_M_recB_phi_cosfit_trueNeutron.pdf");
}

