#pragma once
#include <TCanvas.h>
#include "Palette.h"
#include <TF1.h>
#include <TH1.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TTree.h>

#include <string>
#include <vector>

#include "FitCore.h"
#include "TreeColumns.h"

// The name of the phi on the x axis, from the branch or expression that is fitted: CEPHI_NEUFRAME (true CE phi about the true neutron), RECONPHI_ESTNEU, RECONPHIBL_ESTNEU (block), or the lab phi (no neutron frame).
inline std::string PhiAxisName(const std::string& phiBranch){
    if (phiBranch.rfind("CEPhiNeu", 0) == 0) return "CEPHI_NEUFRAME";
    if (phiBranch.rfind("phiReDLab", 0) == 0) return "recon block phi about the lab z axis";
    if (phiBranch.rfind("phiReLab", 0) == 0) return "recon phi about the lab z axis";
    if (phiBranch.rfind("phiReD", 0) == 0 || phiBranch.rfind("phiRelND", 0) == 0) return "RECONPHIBL_ESTNEU"; // phiRelND: the same angle in AnglesRaw (rad), used by PlotReconNew.cpp
    if (phiBranch.rfind("phiRe", 0) == 0) return "RECONPHI_ESTNEU";
    if (phiBranch.find("LastPhiLab") != std::string::npos) return "CE phi about the lab z axis";
    return "Angle";
}
// phi histogram (50 bins) of t with the polarization fit N(phi) = N0 (1 + P A_eff cos(phi)): N0 = entries / bins (the protons that pass cut) and A_eff = sum_i N_i AP(theta_i) / N0 are fixed, only P is fitted.
// theta_i is the own theta of each of those protons, no theta bins (AEffOfTheta, FitCore.h); thetaBranch is the theta of the plot in degrees, the CETHETA_NEUFRAME for a CEPHI_NEUFRAME and the RECONTHETA_ESTNEU for a RECONPHI_ESTNEU.
// Used by PlotAngle.cpp (3_ReconCalcVals).
// notes: extra text lines under the fit formula (which frame the angles are in, which neutron), none for the 3_ReconCalcVals plots.
// range: {thetaMin, thetaMax, phiHalf} branches (degrees, degrees, rad) of the est neutron range of each proton (CalcAPandT.cpp), none for a fit with no est neutron frame. With them A_eff gets the
// theta range error (A at the two ends, AMean) and P the phi range error (PhiRange), and the side box splits the errors of A_eff and P by source.
// range2: the same three branches for the block grid range alone (the slope error of the block position fit, a range of its own in quadrature), only on the block fits.
inline void PhiCosFit(TCanvas* canvas, TTree* t, const std::string& phiBranch, const std::string& thetaBranch, const std::string& cut, const std::string& name, const std::string& title, const std::string& path, const std::vector<std::string>& notes = {}, const std::string& aBranch = "", const std::vector<std::string>& range = {}, const std::vector<std::string>& range2 = {}, int color = -1){ // color: the histogram line, -1 the CE color
    const int nBins = 50;
    canvas->Clear(); canvas->SetRightMargin(0.25); canvas->SetTopMargin(0.17);
    TH1* h = MakeHist(name + "(" + std::to_string(nBins) + ",0,360)");
    std::vector<double> phiAll;
    AMean aMean; // A_eff of the protons of the fit, with its error
    PhiRange phiRange, phiRange2; // the est neutron phi range of those protons, and the block grid one
    if (h){
        std::vector<std::string> exprs = {phiBranch, thetaBranch};
        if (!aBranch.empty()) exprs.push_back(aBranch); // the A the model weighted each event with (LastAP), used as it is
        const size_t iRange = exprs.size();
        if (range.size() == 3) for (const auto& r : range) exprs.push_back(r); // est neutron range: thetaMin, thetaMax (degrees), phi half width (rad)
        const size_t iRange2 = exprs.size();
        if (range2.size() == 3) for (const auto& r : range2) exprs.push_back(r); // block grid range: thetaMin, thetaMax (degrees), phi half width (rad)
        const auto cols = Columns(t, exprs, cut); // phi, theta (degrees), A (if aBranch), range (if range), range2 (if range2)
        phiAll = cols[0];
        const double apScale = ReadSettings().at("APScale"); // settings.csv from the sim, read once
        for (size_t k = 0; k < phiAll.size(); k++){
            h->Fill(phiAll[k]);
            const double a = aBranch.empty() ? AEffOfTheta(cols[1][k], apScale) : cols[2][k];
            const bool r1 = range.size() == 3 && std::isfinite(cols[iRange][k]) && std::isfinite(cols[iRange + 1][k]);
            const bool r2 = range2.size() == 3 && std::isfinite(cols[iRange2][k]) && std::isfinite(cols[iRange2 + 1][k]);
            aMean.Add(a, r1 ? AEffOfTheta(cols[iRange][k], apScale) : a, r1 ? AEffOfTheta(cols[iRange + 1][k], apScale) : a, r2 ? AEffOfTheta(cols[iRange2][k], apScale) : a, r2 ? AEffOfTheta(cols[iRange2 + 1][k], apScale) : a);
            if (r1 && std::isfinite(cols[iRange + 2][k])) phiRange.Add(cols[iRange + 2][k]);
            if (r2 && std::isfinite(cols[iRange2 + 2][k])) phiRange2.Add(cols[iRange2 + 2][k]);
        }
        h->Draw();
    }
    if (!h || h->GetEntries() < 10){canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1); return;}
    const double nFit = h->GetEntries();
    h->SetTitle((title + Form(", %d phi bins of %.1f deg", nBins, 360.0 / nBins) + ";" + PhiAxisName(phiBranch) + " (degrees);Count").c_str());
    h->SetLineColor(color >= 0 ? color : CbCE()); // CE protons are vermillion
    const bool block = phiBranch.rfind("phiReD", 0) == 0 || phiBranch.rfind("phiRelND", 0) == 0; // the block phi fit leaves out the spike at 0 (kBlockPhiLo to kBlockPhiHi), the likelihood is normalised over that range
    const PFit uf = UnbinnedP(phiAll, aMean, block ? kBlockPhiLo : 0, block ? kBlockPhiHi : 360, nullptr, true, range.size() == 3 ? &phiRange : nullptr, range2.size() == 3 ? &phiRange2 : nullptr); // unbinned fit of every proton's phi (FitCore.h), A fixed to A_eff
    TF1& fit = *AddUnbinnedFit(h, "fPhiCos", uf, aMean.Mean(), aMean.Err(), true);
    fit.SetParameter(0, nFit / nBins);
    fit.SetLineColor(kBlack);
    h->SetMinimum(0); h->SetMaximum(1.1 * h->GetMaximum());
    h->SetName("Fit results"); // the first line of the stats box is the histogram name
    canvas->Update();
    if (auto* st = (TPaveStats*)h->FindObject("stats")){st->SetX1NDC(0.76); st->SetX2NDC(0.995); st->SetY1NDC(0.55); st->SetY2NDC(0.83);}
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
    tx.DrawLatex(0.13, 0.90, "Unbinned fit  N(#phi) #propto N_{0} (1 + P #upoint A_{eff} cos#phi) (every proton, no bins), N_{0} = protons / bins, fixed");
    tx.DrawLatex(0.13, 0.865, Form(aBranch.empty() ? "%s = %.4f #pm %.4f (fixed, per-proton #theta),   P = %.2f #pm %.2f" : "%s = %.4f #pm %.4f (fixed, AP = LastAP of each CE),   P = %.2f #pm %.2f", AEffEq(), aMean.Mean(), aMean.Err(), uf.p, uf.err));
    tx.SetTextSize(0.026);
    DrawErrorBySource(aMean, uf, range.size() == 3); // side box (under the stats box): the errors of A_eff and P by source
    for (size_t i = 0; i < notes.size(); i++) tx.DrawLatex(0.13, 0.83 - 0.034 * i, notes[i].c_str());
    canvas->Modified(); canvas->Update();
    canvas->SaveAs(path.c_str());
    canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
}
