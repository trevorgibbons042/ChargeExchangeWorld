#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TMath.h>
#include <TStyle.h>
#include <ROOT/RDataFrame.hxx>

#include <cmath>
#include <cstdio>
#include <vector>

#include "FitCore.h"

// Unbinned fit of P to the recon block phi of the CE protons, the sample of 3_ReconCalcVals/5_phiReD_cosfit.pdf (PlotAngle.cpp): CalcValsReD of cuts/cut_APandT.root, phiReD from kBlockPhiLo to
// kBlockPhiHi, <A> fixed to the mean apCE of the CE tree. Each proton enters the likelihood with its own phi, no bins: density (1 + P <A> cos phi) / N(P) on [lo, hi],
// N(P) = (hi - lo) + P <A> (sin hi - sin lo) (angles in rad). P is found with Newton's method and its error is 1 / sqrt(-d2lnL/dP2). No limit on P (the binned fit limits it to -1 to 1).
// Writes UnbinnedPhiReD.pdf with the 50 bin histogram of the binned plot, the binned fit (dashed black) and the unbinned fit (blue). Run in a run folder: root -l -b -q UnbinnedPhiReD.cpp
void UnbinnedPhiReD(){
    ROOT::EnableImplicitMT(NThreads());
    ROOT::RDataFrame dfCE("CE", "cuts/cut_APandT.root"), dfPhi("CalcValsReD", "cuts/cut_APandT.root");
    const double aMean = *dfCE.Mean("apCE");
    auto phiTake = dfPhi.Filter(Form("phiReD >= %g && phiReD <= %g", kBlockPhiLo, kBlockPhiHi)).Take<double>("phiReD");
    const std::vector<double> phi = *phiTake;
    const double n = (double)phi.size(), lo = kBlockPhiLo * TMath::Pi() / 180, hi = kBlockPhiHi * TMath::Pi() / 180;
    const double sd = std::sin(hi) - std::sin(lo);

    const PFit uf = UnbinnedP(phi, aMean, kBlockPhiLo, kBlockPhiHi); // FitCore.h, the same fit as every plot of the plot chain
    const double p = uf.p, err = uf.err;

    TH1D h("hUnb", "", 50, 0, 360);
    for (double x : phi) h.Fill(x);
    TF1 fb("fBinned", "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
    fb.FixParameter(0, h.GetEntries() / 50); fb.FixParameter(2, aMean); fb.SetParameter(1, 0.5); fb.SetParLimits(1, -1, 1);
    h.Fit(&fb, "Q0");
    const double w = 360.0 / 50, k = n * (TMath::Pi() / 180) * w / ((hi - lo) + p * aMean * sd); // expected counts per bin of the unbinned fit
    TF1 fu("fUnbinned", Form("%g * (1 + %g * %g * cos(x * TMath::Pi() / 180))", k, p, aMean), 0.5, 359.5);

    TCanvas cv("cvUnb", "", 1000, 800);
    cv.SetLeftMargin(0.14); cv.SetTopMargin(0.17); cv.SetRightMargin(0.04);
    h.SetStats(0); h.SetMinimum(0); h.SetMaximum(1.25 * h.GetMaximum());
    h.SetLineColor(kRed);
    h.SetTitle("Reconstructed CE protons, recon block phi Distribution, fit phi 0.5 to 359.5 deg (spike at 0 left out), no other cuts;Angle (degrees);Count");
    h.Draw("hist");
    fb.SetLineColor(kBlack); fb.SetLineStyle(2); fb.SetLineWidth(3); fb.Draw("same");
    fu.SetLineColor(kBlue); fu.SetLineWidth(3); fu.Draw("same");
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032);
    tx.DrawLatex(0.16, 0.86, "N(#phi) #propto 1 + P #upoint #LTA#GT cos#phi, #LTA#GT fixed");
    tx.SetTextColor(kBlue); tx.DrawLatex(0.16, 0.81, Form("unbinned fit (blue): P = %.3f #pm %.3f,  %.0f protons, #LTA#GT = %.4f", p, err, n, aMean));
    tx.SetTextColor(kBlack); tx.DrawLatex(0.16, 0.76, Form("binned #chi^{2} fit (dashed): P = %.3f #pm %.3f,  #chi^{2}/ndf = %.0f / %d", fb.GetParameter(1), fb.GetParError(1), h.Chisquare(&fb), 49));
    cv.SaveAs("UnbinnedPhiReD.pdf");
    printf("unbinned P = %.4f +- %.4f (n = %.0f, <A> = %.5f), binned P = %.4f +- %.4f\n", p, err, n, aMean, fb.GetParameter(1), fb.GetParError(1));
}
