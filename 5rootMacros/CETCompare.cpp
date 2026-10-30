#include <TCanvas.h>
#include <TChain.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include <string>

#include <TROOT.h>

#include "CutCache.h"
#include "DrawPlan.h"
#include "FitCore.h"

// t of the charge exchange (LastT) for the CE in the carbon blocks (CarbonCeCut) and for the CE in the silicon tracker planes (every other CE),
// ChargeExchange_values of the sim. Normalized to unit area in a second plot to compare the shapes. plots/4_CEblock/1_CarbonVsSi/t_carbon_vs_tracker*.pdf
void CETCompare(){
    const std::string carbon = CarbonCeCut(), tracker = "didChargeExchange>0 && !(" + carbon + ")";
    // the three selections come from the cut cache (CutCache.h): every row of the tree passes, so the draws need no selection
    const CacheTree cAll = GetCutTree("output0_t*.root", "ChargeExchange_values", {"didChargeExchange>0"}); // all CE
    const CacheTree cCarbon = GetCutTree("output0_t*.root", "ChargeExchange_values", {carbon});             // CE in the carbon blocks (seed cut_CE.root)
    const CacheTree cTracker = GetCutTree("output0_t*.root", "ChargeExchange_values", {"didChargeExchange>0", "!(" + carbon + ")"}); // CE in the silicon trackers, from the all-CE cache
    TTree *tAll = cAll.tree, *tCarbon = cCarbon.tree, *tTracker = cTracker.tree;
    if (!tAll || !tCarbon || !tTracker){Error("CETCompare", "cut caches missing (run CutEvents.cpp first)"); return;}
    gROOT->cd();
    ROOT::EnableImplicitMT(NThreads());
    TH1D hC("hTCarbon", "", 100, 0, 0.05), hT("hTTracker", "", 100, 0, 0.05);
    DrawPlan plan; // the five fills (z all/carbon/tracker, t carbon/tracker) are collected in one multithreaded pass per tree, then filled here in entry order
    gSystem->mkdir("plots/4_CEblock/1_CarbonVsSi", true);
    gStyle->SetOptStat(0);
    TCanvas c("cT", "", 1500, 1200); c.SetLeftMargin(0.14); c.SetBottomMargin(0.13); c.SetRightMargin(0.05);
    // z where the CE happens, all CE of the sim: carbon blocks, silicon trackers, and both
    {
        // the stack of the run (2construction.cc): nBlocks carbon blocks and (nBlocks + 1) * nBlocksTracking tracker planes, dz apart, centred at 4 m
        const auto st = ReadSettings();
        const int nElements = (int)std::lround((st.at("nBlocks") + 1) * st.at("nBlocksTracking") + st.at("nBlocks"));
        const double halfSpan = (nElements - 1) * st.at("dz") * 1000.0 / 2.0;
        const double z0 = 4000 - halfSpan - 200, z1 = 4000 + halfSpan + 200; const int nz = (int)((z1 - z0) / 5.0);
        TH1D hAll("hZAll", "", nz, z0, z1), hZC("hZCarbon", "", nz, z0, z1), hZT("hZTracker", "", nz, z0, z1);
        plan.H(tAll, "fZ>>hZAll", "", "", &hAll); plan.H(tCarbon, "fZ>>hZCarbon", "", "", &hZC); plan.H(tTracker, "fZ>>hZTracker", "", "", &hZT);
        plan.H(tCarbon, "LastT>>hTCarbon", "", "", &hC); plan.H(tTracker, "LastT>>hTTracker", "", "", &hT);
        plan.Run();
        const long long nC = tCarbon->GetEntries(), nT = tTracker->GetEntries(), nA = tAll->GetEntries();
        hAll.SetLineColor(kBlack); hZC.SetLineColor(kBlue); hZT.SetLineColor(kRed);
        for (TH1D* h : {&hAll, &hZC, &hZT}) h->SetLineWidth(2);
        hAll.SetTitle("Z position where Charge Exchange happens, all CE protons (sim);Z Position (mm);Count");
        hAll.SetMaximum(1.15 * hAll.GetMaximum());
        hAll.Draw("hist"); hZC.Draw("hist same"); hZT.Draw("hist same");
        TLegend lz(0.62, 0.74, 0.95, 0.90); lz.SetTextSize(0.024);
        lz.AddEntry(&hAll, Form("all CE, %lld", nA), "l");
        lz.AddEntry(&hZC, Form("CE in carbon blocks, %lld", nC), "l");
        lz.AddEntry(&hZT, Form("CE in silicon trackers, %lld", nT), "l");
        lz.Draw();
        c.SaveAs("plots/4_CEblock/1_CarbonVsSi/z_CE_blocks_vs_trackers.pdf");
    }
    printf("carbon CE: N=%.0f mean t=%.5g | tracker CE: N=%.0f mean t=%.5g\n", hC.GetEntries(), hC.GetMean(), hT.GetEntries(), hT.GetMean());
    for (const bool norm : {false, true}){
        TH1D* a = (TH1D*)hC.Clone(), *b = (TH1D*)hT.Clone();
        if (norm){a->Scale(1.0 / a->Integral()); b->Scale(1.0 / b->Integral());}
        a->SetLineColor(kBlue); b->SetLineColor(kRed); a->SetLineWidth(2); b->SetLineWidth(2);
        a->SetTitle(Form("t of the charge exchange, all CE protons (sim), carbon blocks vs silicon trackers%s;t;%s", norm ? ", unit area" : "", norm ? "Fraction" : "Count"));
        a->SetMaximum(1.15 * std::max(a->GetMaximum(), b->GetMaximum()));
        a->Draw("hist"); b->Draw("hist same");
        TLegend leg(0.40, 0.78, 0.94, 0.90); leg.SetTextSize(0.025);
        leg.AddEntry(a, Form("CE in carbon blocks (%.0f), mean t = %.4g", hC.GetEntries(), hC.GetMean()), "l");
        leg.AddEntry(b, Form("CE in silicon trackers (%.0f), mean t = %.4g", hT.GetEntries(), hT.GetMean()), "l");
        leg.Draw();
        c.SaveAs(norm ? "plots/4_CEblock/1_CarbonVsSi/t_carbon_vs_tracker_norm.pdf" : "plots/4_CEblock/1_CarbonVsSi/t_carbon_vs_tracker.pdf");
    }
}
