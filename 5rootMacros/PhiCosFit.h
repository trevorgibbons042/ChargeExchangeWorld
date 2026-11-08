#pragma once
#include <TCanvas.h>
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

// phi histogram (50 bins) of t with the polarization fit N(phi) = N0 (1 + P A_eff cos(phi)): N0 = entries / bins (the protons that pass cut) and A_eff = sum_i N_i AP(theta_i) / N0 are fixed, only P is fitted.
// theta_i is the own theta of each of those protons, no theta bins (AEffOfTheta, FitCore.h); thetaBranch is the theta of the plot in degrees, the CETHETA_NEUFRAME for a CEPHI_NEUFRAME and the RECONTHETA_ESTNEU for a RECONPHI_ESTNEU.
// Used by PlotAngle.cpp (3_ReconCalcVals).
// notes: extra text lines under the fit formula (which frame the angles are in, which neutron), none for the 3_ReconCalcVals plots.
// range: {thetaMin, thetaMax, phiHalf} branches (degrees, degrees, rad) of the est neutron range of each proton (CalcAPandT.cpp), none for a fit with no est neutron frame. With them A_eff gets the
// theta range error (A at the two ends, AMean) and P the phi range error (PhiRange), and the side box splits the errors of A_eff and P by source.
inline void PhiCosFit(TCanvas* canvas, TTree* t, const std::string& phiBranch, const std::string& thetaBranch, const std::string& cut, const std::string& name, const std::string& title, const std::string& path, const std::vector<std::string>& notes = {}, const std::string& aBranch = "", const std::vector<std::string>& range = {}){
    const int nBins = 50;
    canvas->Clear(); canvas->SetRightMargin(0.25); canvas->SetTopMargin(0.17);
    TH1* h = MakeHist(name + "(" + std::to_string(nBins) + ",0,360)");
    std::vector<double> phiAll;
    AMean aMean; // A_eff of the protons of the fit, with its error
    PhiRange phiRange; // the est neutron phi range of those protons
    if (h){
        std::vector<std::string> exprs = {phiBranch, thetaBranch};
        if (!aBranch.empty()) exprs.push_back(aBranch); // the A the model weighted each event with (LastAP), used as it is
        const size_t iRange = exprs.size();
        if (range.size() == 3) for (const auto& r : range) exprs.push_back(r); // est neutron range: thetaMin, thetaMax (degrees), phi half width (rad)
        const auto cols = Columns(t, exprs, cut); // phi, theta (degrees), A (if aBranch), range (if range)
        phiAll = cols[0];
        const double apScale = ReadSettings().at("APScale"); // settings.csv from the sim, read once
        for (size_t k = 0; k < phiAll.size(); k++){
            h->Fill(phiAll[k]);
            const double a = aBranch.empty() ? AEffOfTheta(cols[1][k], apScale) : cols[2][k];
            if (range.size() == 3 && std::isfinite(cols[iRange][k]) && std::isfinite(cols[iRange + 1][k])){
                aMean.Add(a, AEffOfTheta(cols[iRange][k], apScale), AEffOfTheta(cols[iRange + 1][k], apScale));
                if (std::isfinite(cols[iRange + 2][k])) phiRange.Add(cols[iRange + 2][k]);
            } else aMean.Add(a);
        }
        h->Draw();
    }
    if (!h || h->GetEntries() < 10){canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1); return;}
    const double nFit = h->GetEntries();
    h->SetTitle((title + Form(", %d phi bins of %.1f deg", nBins, 360.0 / nBins) + ";Angle (degrees);Count").c_str());
    h->SetLineColor(kRed); // CE protons are red
    const bool block = phiBranch.rfind("phiReD", 0) == 0; // the block phi fit leaves out the spike at 0 (kBlockPhiLo to kBlockPhiHi), the likelihood is normalised over that range
    const PFit uf = UnbinnedP(phiAll, aMean, block ? kBlockPhiLo : 0, block ? kBlockPhiHi : 360, nullptr, true, range.size() == 3 ? &phiRange : nullptr); // unbinned fit of every proton's phi (FitCore.h), A fixed to A_eff
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
    // side box (under the stats box): the errors of A_eff and P by source, each in quadrature in the total on the title line
    TPaveText* eb = new TPaveText(0.755, 0.06, 0.995, 0.52, "NDC");
    eb->SetTextAlign(12); eb->SetTextSize(0.024); eb->SetFillColor(kWhite); eb->SetBorderSize(1);
    eb->AddText("#bf{Error by source}");
    eb->AddText(Form("#bf{%s} #pm %.4f", AEffEq(), aMean.Err()));
    eb->AddText(Form(" est n #theta range  %.4f", aMean.RangeErr()));
    eb->AddText(Form(" statistics (#theta spread)  %.4f", aMean.StatErr()));
    eb->AddText(Form(" energy range (beam)  %.4f", aMean.BeamErr()));
    eb->AddText(Form("#bf{P} #pm %.3f", uf.err));
    eb->AddText(Form(" est n #phi range  %.3f", uf.errPhi));
    eb->AddText(Form(" est n #theta (via A)  %.3f", uf.errARange));
    eb->AddText(Form(" #theta spread (via A)  %.3f", uf.errAStat));
    eb->AddText(Form(" energy range (via A)  %.3f", uf.errABeam));
    eb->AddText(Form(" statistics (likelihood)  %.3f", uf.errFit));
    if (range.size() != 3) eb->AddText("no est neutron range here");
    eb->Draw();
    for (size_t i = 0; i < notes.size(); i++) tx.DrawLatex(0.13, 0.83 - 0.034 * i, notes[i].c_str());
    canvas->Modified(); canvas->Update();
    canvas->SaveAs(path.c_str());
    canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
}
