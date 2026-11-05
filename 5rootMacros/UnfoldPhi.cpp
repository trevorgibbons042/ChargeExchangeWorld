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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "FitCore.h"

// Iterative Bayesian (D'Agostini) unfolding of the (theta, phi) distribution of the CE protons from the recon (block) angles back to the neutron frame, then the P fit of PlotAngle.cpp on
// the unfolded phi: N(phi) = N0 (1 + P A_eff cos(phi)), A fixed. Reads cuts/smear_matrix.root (SmearMatrix.cpp), run SmearMatrix.cpp first, from a run folder.
// Two sets: recon block (R_thetaphi, M_thetaphi, F_thetaphi) and recon (the _recon ones). Per set, four curves in plots/8_smear/5_unfold/:
//   CE phi (lab)        true angles about the lab z axis (the CE phi of 3_ReconCalcVals/1_phiCE_cosfit), no detector
//   CE phi (neutron frame)             true angles about the neutron direction, the answer the unfolding aims at
//   recon (block) phi   recon (block) angles of the CE protons plus the fakes (the protons that are not from a carbon CE track)
//   unfolded recon (block) phi   the recon (block) phi minus the fakes, unfolded with the response R / T and a flat prior, kIterDraw iterations
// The response is P(recon bin j | true bin i) = R_ij / T_i with T the truth of every row (so a column sums to the efficiency, the misses are the rest). The bins are the flattened
// (theta slice x phi bin) numbers of SmearMatrix.cpp. The error of the unfolded bins and of P comes from kToys Poisson toys of the measured histogram (the statistical error of the
// matrix is not in it). A_eff of a curve is the mean of AP_theta (A_eff of the neutron theta per true slice) over the theta slices of that curve.
// Writes per set (5 files each, block 1 to 5, recon 6 to 10) one graph per curve (CE_phi_lab, CE_neutronframe, CE_reconblock or CE_recon, unfolded_recon_phi, same y axis), the P vs iteration graph, and unfold_P.csv.

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

// one set: the unfolding, its toys, the fits, the two pages and the csv rows
static void UnfoldSet(int setIdx, const std::string& tag, const std::string& label, const std::string& phiName, const TH2D& R, const TH1D& T, const TH1D& M, const TH1D& F, const TH2D& RBeam, const TProfile& ap, double tMax, FILE* csv){
    const std::string dir = setIdx == 0 ? "2_reconBlock" : "3_recon"; // folder of this set in plots/8_smear/5_unfold/
    const Response resp = MakeResponse(R, T);
    const std::vector<double> truthNeutron = ContentsOf(T), fakes = ContentsOf(F);
    std::vector<double> measured(kBins), m0(kBins), truthLab(kBins, 0.0);
    for (int j = 0; j < kBins; j++){measured[j] = M.GetBinContent(j + 1) + fakes[j]; m0[j] = std::max(measured[j] - fakes[j], 0.0);}
    for (int j = 0; j < kBins; j++) for (int i = 0; i < kBins; i++) truthLab[j] += RBeam.GetBinContent(i + 1, j + 1);

    const auto nominal = Bayes(resp, m0, kIterMax);
    // toys of the measured histogram, unfolded: the phi histogram of each toy at each iteration
    TRandom3 rng(12345);
    std::vector<double> toyPhi((size_t)kIterMax * kToys * kPhiBins);
    for (int t = 0; t < kToys; t++){
        std::vector<double> mt(kBins);
        for (int j = 0; j < kBins; j++) mt[j] = std::max((double)rng.Poisson(measured[j]) - fakes[j], 0.0);
        const auto u = Bayes(resp, mt, kIterMax);
        double sl[kSlices];
        for (int it = 0; it < kIterMax; it++) PhiOf(u[it], &toyPhi[((size_t)it * kToys + t) * kPhiBins], sl);
    }

    double phiU[kIterMax][kPhiBins], sigU[kIterMax][kPhiBins], aU[kIterMax];
    PResult pU[kIterMax];
    double sliceN[kSlices];
    for (int it = 0; it < kIterMax; it++){
        PhiOf(nominal[it], phiU[it], sliceN);
        aU[it] = AEffOf(sliceN, ap);
        for (int b = 0; b < kPhiBins; b++){
            double s = 0, s2 = 0;
            for (int t = 0; t < kToys; t++){const double x = toyPhi[((size_t)it * kToys + t) * kPhiBins + b]; s += x; s2 += x * x;}
            const double mean = s / kToys;
            sigU[it][b] = std::max(std::sqrt(std::max(0.0, s2 / kToys - mean * mean)), 1e-6);
        }
        pU[it] = FitP(phiU[it], sigU[it], aU[it]);
        // the error of P is the spread of the P of the toys, fitted with the same bin errors (this keeps the correlation of the bins)
        double sp = 0, sp2 = 0; int ok = 0;
        for (int t = 0; t < kToys; t++){
            const PResult r = FitP(&toyPhi[((size_t)it * kToys + t) * kPhiBins], sigU[it], aU[it]);
            if (!r.ok) continue;
            sp += r.p; sp2 += r.p * r.p; ok++;
        }
        if (ok > 1){const double mean = sp / ok; pU[it].err = std::sqrt(std::max(0.0, sp2 / ok - mean * mean));}
    }

    // the three curves without unfolding
    double phiTN[kPhiBins], phiTL[kPhiBins], phiM[kPhiBins], sigTN[kPhiBins], sigTL[kPhiBins], sigM[kPhiBins], slTN[kSlices], slTL[kSlices], slM[kSlices];
    PhiOf(truthNeutron, phiTN, slTN); PhiOf(truthLab, phiTL, slTL); PhiOf(measured, phiM, slM);
    for (int b = 0; b < kPhiBins; b++){sigTN[b] = std::sqrt(std::max(phiTN[b], 1.0)); sigTL[b] = std::sqrt(std::max(phiTL[b], 1.0)); sigM[b] = std::sqrt(std::max(phiM[b], 1.0));}
    const double aTN = AEffOf(slTN, ap), aTL = AEffOf(slTL, ap), aM = AEffOf(slM, ap);
    const PResult pTN = FitP(phiTN, sigTN, aTN), pTL = FitP(phiTL, sigTL, aTL), pM = FitP(phiM, sigM, aM);

    fprintf(csv, "%s,CE_phi_lab,0,%.5f,%.5f\n%s,CE_neutronframe,0,%.5f,%.5f\n%s,%s,0,%.5f,%.5f\n", tag.c_str(), pTL.p, pTL.err, tag.c_str(), pTN.p, pTN.err, tag.c_str(), phiName.c_str(), pM.p, pM.err);
    for (int it = 0; it < kIterMax; it++) fprintf(csv, "%s,unfolded %s,%d,%.5f,%.5f\n", tag.c_str(), phiName.c_str(), it + 1, pU[it].p, pU[it].err);
    printf("UnfoldPhi %s  P CE phi (lab) %.3f +- %.3f  CE phi (neutron frame) %.3f +- %.3f  CE %s (lab frame) %.3f +- %.3f  unfolded %s (%d it) %.3f +- %.3f\n", tag.c_str(), pTL.p, pTL.err, pTN.p, pTN.err, phiName.c_str(), pM.p, pM.err, phiName.c_str(), kIterDraw, pU[kIterDraw - 1].p, pU[kIterDraw - 1].err);

    const std::string cuts = Form("CE protons, theta 0 to %.1f deg, no other cuts", tMax);
    {   // four graphs, one per curve, all with the same y axis: the phi histogram with its fit and P
        const int d = kIterDraw - 1;
        struct Curve {std::string file, name; const double* v; const double* sig; PResult r; double a; int color; bool points;};
        const std::vector<Curve> curves = {
            {"CE_phi_lab", "CE phi (lab)", phiTL, sigTL, pTL, aTL, kBlue + 1, false},
            {"CE_neutronframe", "CE phi (neutron frame)", phiTN, sigTN, pTN, aTN, kGreen + 2, false},
            {setIdx == 0 ? "CE_reconblock" : "CE_recon", "CE " + phiName + " (lab frame)", phiM, sigM, pM, aM, kRed + 1, false},
            {"unfolded_recon_phi", "unfolded " + phiName, phiU[d], sigU[d], pU[d], aU[d], kBlack, true}};
        double top = 0;
        for (const Curve& cv : curves) for (int b = 0; b < kPhiBins; b++) top = std::max(top, cv.v[b] + cv.sig[b]);
        for (size_t k = 0; k < curves.size(); k++){
            const Curve& cv = curves[k];
            if (k < 2 && setIdx != 0) continue; // the two CE truth curves (lab, neutron frame) are the same in both sets: one copy, from set 0
            const std::string path = k < 2 ? Form("plots/8_smear/5_unfold/1_CE_physicsangles/%d_%s.pdf", (int)k + 1, cv.file.c_str()) : Form("plots/8_smear/5_unfold/%s/%d_%s.pdf", dir.c_str(), (int)k - 1, cv.file.c_str());
            const std::string extra = cv.points ? Form(", unfolded to the neutron frame, %d iterations, flat prior", kIterDraw) : "";
            TCanvas c("cUnfold", "Plots", 1500, 1200); // the look of the PhiCosFit plots of PlotAngle.cpp
            c.SetLeftMargin(0.14); c.SetBottomMargin(0.13); c.SetRightMargin(0.25); c.SetTopMargin(0.17);
            TH1D* h = PhiHist(Form("h%d_%s", (int)k, tag.c_str()), cv.v, cv.sig, kRed, 1, cv.points ? 20 : 0);
            double nSum = 0; for (int b = 0; b < kPhiBins; b++) nSum += cv.v[b];
            h->SetEntries(nSum); h->SetMarkerSize(0.8);
            h->SetTitle(PlotTitle(path, Form("Phi, %s, %s%s, %d phi bins of %.1f deg;Angle (degrees);Count", cv.name.c_str(), cuts.c_str(), extra.c_str(), kPhiBins, 360.0 / kPhiBins)).c_str());
            h->SetMinimum(0); h->SetMaximum(1.1 * top);
            TF1* f = new TF1(Form("f%d_%s", (int)k, tag.c_str()), "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
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
            tx.DrawLatex(0.13, 0.865, Form("A_{eff} = #Sigma N_{i} AP(#theta_{i}) / N_{tot} = %.4f (fixed, %.2f deg theta slices),   P = %.2f #pm %.2f", cv.a, tMax / kSlices, cv.r.p, cv.r.err));
            c.Modified(); c.Update();
            c.SaveAs(path.c_str());
        }
    }
    {   // page 2: P of the unfolded phi against the number of iterations, with the three other P as lines
        const std::string path = Form("plots/8_smear/5_unfold/%s/3_P_vs_iteration.pdf", dir.c_str());
        TCanvas c("cUnfoldIt", "Plots", 1500, 1200);
        PageMargins(c);
        TGraphErrors g;
        for (int it = 0; it < kIterMax; it++){g.SetPoint(it, it + 1, pU[it].p); g.SetPointError(it, 0, pU[it].err);}
        g.SetMarkerStyle(20); g.SetMarkerSize(1.2);
        g.SetTitle(PlotTitle(path, Form("P of the unfolded phi of %s against the iterations, flat prior, %s, %d phi bins of %.1f deg;Iterations;P", label.c_str(), cuts.c_str(), kPhiBins, 360.0 / kPhiBins)).c_str());
        g.Draw("APE");
        g.GetXaxis()->SetLimits(0.5, kIterMax + 0.5);
        double lo = std::min({pM.p, pTN.p, pTL.p}), hi = std::max({pM.p, pTN.p, pTL.p});
        for (int it = 0; it < kIterMax; it++){lo = std::min(lo, pU[it].p - pU[it].err); hi = std::max(hi, pU[it].p + pU[it].err);}
        const double pad = 0.1 * std::max(hi - lo, 0.05);
        g.GetYaxis()->SetRangeUser(lo - pad, hi + pad);
        TLine lTN(0.5, pTN.p, kIterMax + 0.5, pTN.p), lTL(0.5, pTL.p, kIterMax + 0.5, pTL.p), lM(0.5, pM.p, kIterMax + 0.5, pM.p);
        lTN.SetLineColor(kGreen + 2); lTL.SetLineColor(kBlue + 1); lM.SetLineColor(kRed + 1);
        for (TLine* l : {&lTN, &lTL, &lM}){l->SetLineWidth(2); l->SetLineStyle(2); l->Draw();}
        TLegend leg(0.715, 0.70, 0.99, 0.88);
        leg.AddEntry(&g, Form("unfolded %s", phiName.c_str()), "pe"); leg.AddEntry(&lTN, "CE phi (neutron frame)", "l"); leg.AddEntry(&lTL, "CE phi (lab)", "l"); leg.AddEntry(&lM, ("CE " + phiName + " (lab frame)").c_str(), "l");
        leg.Draw();
        c.SaveAs(path.c_str());
    }
}

void UnfoldPhi(){
    FastJit();
    TH1::AddDirectory(kFALSE);
    TFile* file = TFile::Open("cuts/smear_matrix.root", "READ");
    if (!file || file->IsZombie()){printf("UnfoldPhi: cuts/smear_matrix.root missing (run SmearMatrix.cpp)\n"); return;}
    auto get = [&](const char* name){TObject* o = file->Get(name); if (!o) printf("UnfoldPhi: %s missing in smear_matrix.root (rerun SmearMatrix.cpp)\n", name); return o;};
    const TH2D* rBlock = (TH2D*)get("R_thetaphi"); const TH2D* rRecon = (TH2D*)get("R_thetaphi_recon"); const TH2D* rBeam = (TH2D*)get("R_beam_thetaphi");
    const TH1D* t = (TH1D*)get("T_thetaphi"); const TH1D* mBlock = (TH1D*)get("M_thetaphi"); const TH1D* mRecon = (TH1D*)get("M_thetaphi_recon");
    const TH1D* fBlock = (TH1D*)get("F_thetaphi"); const TH1D* fRecon = (TH1D*)get("F_thetaphi_recon"); const TProfile* ap = (TProfile*)get("AP_theta");
    if (!rBlock || !rRecon || !rBeam || !t || !mBlock || !mRecon || !fBlock || !fRecon || !ap) return;
    if (t->GetNbinsX() != kBins){printf("UnfoldPhi: the matrix has %d bins, the macro expects %d\n", t->GetNbinsX(), kBins); return;}

    for (const char* d : {"1_CE_physicsangles", "2_reconBlock", "3_recon"}) gSystem->mkdir(Form("plots/8_smear/5_unfold/%s", d), kTRUE);
    PlainFonts(); gStyle->SetOptStat(1111); gStyle->SetOptFit(11);
    const double tMax = 5.0 * ThetaScale();
    FILE* csv = fopen("plots/8_smear/5_unfold/unfold_P.csv", "w");
    fprintf(csv, "set,curve,iteration,P,P_err\n");
    UnfoldSet(0, "block", "CE recon block protons", "recon block phi", *rBlock, *t, *mBlock, *fBlock, *rBeam, *ap, tMax, csv);
    UnfoldSet(1, "recon", "CE recon protons", "recon phi", *rRecon, *t, *mRecon, *fRecon, *rBeam, *ap, tMax, csv);
    fclose(csv);
}
