#include <TBox.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <TROOT.h>
#include "FitCore.h"
#include "QualityCut.h"
#include "TreeColumns.h"

// Picks the thresholds in QualityCut.h. Reads AnglesRaw (every fitted proton, before the quality cut, ReconAngleAll.cpp) and compares the
// CE protons (isFromCE == 1) with all recon protons for the four quality variables (residual and impact, true and block position).
// The four cuts are set at one common percentile of the CE distribution, chosen so that the four together keep 95% of the CE protons.
// Writes plots/6_cuts/3_ResCuts/: one CE vs recon protons plot per variable with its threshold, and thresholds.csv. Prints the numbers to paste into QualityCut.h.
static double Quantile(std::vector<double> v, double q){
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, (size_t)(q * v.size()))];
}

void QualityStudy(double keepCE = 0.95){
    FastJit();
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    gStyle->SetPalette(kViridis);
    TFile* f = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree* t = f ? (TTree*)f->Get("AnglesRaw") : nullptr;
    if (!t){printf("QualityStudy: AnglesRaw missing in cut_angles_all.root (run ReconAngleAll.cpp first)\n"); return;}

    const std::vector<std::string> names = {"resRms", "resDRms", "impact", "impactD"};
    const std::vector<std::string> labels = {"residual, true position (mm)", "residual, block position (mm)", "impact to beam axis, true position (mm)", "impact to beam axis, block position (mm)"};
    std::vector<std::vector<double>> ce(4), recon(4); // [variable][track]; recon = every fitted proton, CE ones included
    // one multithreaded read of every branch the loops below use (TreeColumns.h), the rows in entry order; the optional branches only when the tree has them
    const bool hasR2 = t->GetBranch("r2"), hasCand = t->GetBranch("nProtonsEv");
    std::vector<std::string> want = {"resRms", "resDRms", "impact", "impactD", "isFromCE"}; // columns 0-3 the four variables, 4 isFromCE
    int iR2 = -1, iCand = -1;
    if (hasR2){iR2 = (int)want.size(); want.push_back("r2"); want.push_back("r2D");}
    if (hasCand){iCand = (int)want.size(); for (const char* b : {"nProtonsEv", "nHits", "firstXD", "firstYD", "resHalfRatioD", "outTheta", "outDTheta"}) want.push_back(b);}
    const std::vector<std::vector<double>> col = Columns(t, want);
    const size_t nRows = col[0].size();
    auto isCE = [&](size_t e){return col[4][e] != 0;};
    auto inBase = [&](size_t e){ // baseline of 6_cuts: recon theta from 0 to 1 degree (NaN fails)
        if (!hasCand) return true;
        const double deg = col[iCand + 5][e] * TMath::RadToDeg();
        return deg >= 0 && deg <= 1;
    };
    for (size_t e = 0; e < nRows; e++){
        if (!inBase(e)) continue;
        for (int i = 0; i < 4; i++){
            recon[i].push_back(col[i][e]);
            if (isCE(e)) ce[i].push_back(col[i][e]);
        }
    }
    const size_t nCE = ce[0].size(), nRecon = recon[0].size();
    if (nCE == 0){printf("QualityStudy: no CE protons in AnglesRaw\n"); return;}

    // joint pass: a track passes when all four variables are below their thresholds (the same percentile q of the CE values for each variable)
    auto thresholds = [&](double q){std::vector<double> th(4); for (int i = 0; i < 4; i++) th[i] = Quantile(ce[i], q); return th;};
    auto jointKept = [&](const std::vector<std::vector<double>>& s, const std::vector<double>& th){
        size_t n = s[0].size(), kept = 0;
        for (size_t k = 0; k < n; k++){
            bool ok = true;
            for (int i = 0; i < 4; i++) ok = ok && s[i][k] <= th[i];
            kept += ok;
        }
        return n ? double(kept) / n : 0.0;
    };
    double lo = keepCE, hi = 1.0; // bisection on q, kept fraction grows with q
    for (int it = 0; it < 40; it++){
        const double mid = 0.5 * (lo + hi);
        (jointKept(ce, thresholds(mid)) >= keepCE ? hi : lo) = mid;
    }
    const std::vector<double> th = thresholds(hi);
    const double ceKept = jointKept(ce, th), reconKept = jointKept(recon, th);

    const std::string dir = "plots/6_cuts/3_ResCuts/";
    gSystem->mkdir(dir.c_str(), true);
    FILE* out = fopen((dir + "thresholds.csv").c_str(), "w");
    if (out) fprintf(out, "variable,threshold,CE_kept_alone_pct,recon_kept_alone_pct\n");
    printf("QualityStudy: %zu CE protons, %zu recon protons\n", nCE, nRecon);

    TCanvas c("cQuality", "", 1500, 1200);
    c.SetLeftMargin(0.14); c.SetBottomMargin(0.13); c.SetRightMargin(0.30); c.SetLogy(); // legend and equations in the right margin, outside the graph
    for (int i = 1; i < 2; i++){ // only the block residual is plotted: the true position residual (0) and the impact variables (2, 3) enter the joint cut above but get no plot
        size_t ceBelow = 0, reconBelow = 0;
        for (double x : ce[i]) ceBelow += x <= th[i];
        for (double x : recon[i]) reconBelow += x <= th[i];
        const double cePct = 100.0 * ceBelow / nCE, reconPct = 100.0 * reconBelow / nRecon;
        printf("  %-8s <= %.4g mm   alone keeps %.1f%% CE, %.1f%% recon\n", names[i].c_str(), th[i], cePct, reconPct);
        if (out) fprintf(out, "%s,%.6g,%.3f,%.3f\n", names[i].c_str(), th[i], cePct, reconPct);

        const double xMax = std::max(1e-6, 1.2 * std::max(th[i], Quantile(recon[i], 0.99)));
        TH1D hCE("hQCE", "", 100, 0, xMax), hRecon("hQRecon", "", 100, 0, xMax);
        for (double x : ce[i]) hCE.Fill(x);
        for (double x : recon[i]) hRecon.Fill(x);
        hCE.SetLineColor(kBlue); hRecon.SetLineColor(kRed); hCE.SetLineWidth(2); hRecon.SetLineWidth(2);
        hRecon.SetTitle((labels[i] + ", CE recon protons vs all recon protons, recon theta 0-1 deg, no KE cut;" + labels[i] + ";Count").c_str());
        hRecon.SetStats(0); hCE.SetStats(0);
        hRecon.SetMinimum(0.5);
        hRecon.SetMaximum(2.0 * std::max(hRecon.GetMaximum(), hCE.GetMaximum()));
        hRecon.DrawCopy(); hCE.DrawCopy("same");
        TLegend leg(0.725, 0.80, 0.99, 0.90); leg.SetTextSize(0.014); // no cut line: the thresholds are in thresholds.csv
        leg.AddEntry(&hCE, Form("CE recon protons (%zu)", nCE), "l");
        leg.AddEntry(&hRecon, Form("all recon protons (%zu)", nRecon), "l");
        leg.Draw();
        // definition of the variable (ReconAngleAll.cpp), fitted line x = x_{0} + b_{x} z, y = y_{0} + b_{y} z
        TPaveText eq(0.725, 0.48, 0.99, 0.72, "NDC");
        eq.SetBorderSize(1); eq.SetFillColor(0); eq.SetTextAlign(12); eq.SetTextSize(0.016);
        eq.AddText(i % 2 ? "line fit to block hits:" : "line fit to hits:");
        eq.AddText("#hat{x} = x_{0} + b_{x}z, #hat{y} = y_{0} + b_{y}z");
        if (i < 2){
            eq.AddText("SS_{res} = #Sigma_{i} (x_{i} - #hat{x}_{i})^{2} + (y_{i} - #hat{y}_{i})^{2}");
            eq.AddText("residual = #sqrt{SS_{res} / n}");
        } else {
            eq.AddText("impact = |x_{0} b_{y} - y_{0} b_{x}|");
            eq.AddText("         / #sqrt{b_{x}^{2} + b_{y}^{2}}");
        }
        eq.Draw();
        c.SaveAs((dir + std::to_string(i + 1) + "_" + names[i] + ".pdf").c_str());
    }

    // residual_cut_efficiency.csv: the block residual cut (the only residual cut, QualityCut.h) at the res_scan.sh thresholds and the default,
    // Totals before and after the cut.
    if (FILE* rc = fopen((dir + "residual_cut_efficiency.csv").c_str(), "w")){
        fprintf(rc, "cut,block_residual_max_mm,CE_before,CE_after,CE_kept_pct,recon_before,recon_after,recon_left_out_pct\n");
        auto row = [&](const char* name, double resMax){
            auto kept = [&](const std::vector<std::vector<double>>& s){
                size_t k = 0;
                for (size_t i = 0; i < s[1].size(); i++)
                    k += s[1][i] <= resMax;
                return k;
            };
            const size_t ceAfter = kept(ce), reconAfter = kept(recon);
            fprintf(rc, "%s,%g,%zu,%zu,%.4f,%zu,%zu,%.4f\n", name, resMax, nCE, ceAfter, 100.0 * ceAfter / nCE, nRecon, reconAfter, 100.0 * (nRecon - reconAfter) / nRecon);
        };
        for (double r : {2.0, 3.1, 5.0}) row(Form("block residual < %g mm", r), r);
        row("block residual default", kMaxResidualDMm);
        fclose(rc);
    }

    // Cut candidates, plotted only (nothing is cut): CE vs recon protons, one folder each, from all of AnglesRaw
    struct Cand {std::string branch, folder, label; double lo, hi; int nBins;};
    const std::vector<Cand> cands = {
        {"nProtonsEv", "4_nProtons", "protons with a track in the event", -0.5, 19.5, 20},
        {"nHits", "8_nHits", "tracker hits", -0.5, 19.5, 20}};
    {
        if (hasCand){
            // col[iCand + k]: 0 nProtonsEv, 1 nHits, 2 firstXD, 3 firstYD, 4 resHalfRatioD, 5 outTheta, 6 outDTheta
            std::vector<std::vector<double>> candCE(cands.size()), candRecon(cands.size());
            for (size_t e = 0; e < nRows; e++){
                if (!inBase(e)) continue;
                const double vals[2] = {col[iCand][e], col[iCand + 1][e]}; // in the order of cands
                for (size_t i = 0; i < cands.size(); i++){
                    if (!std::isfinite(vals[i])) continue;
                    candRecon[i].push_back(vals[i]);
                    if (isCE(e)) candCE[i].push_back(vals[i]);
                }
            }
            for (size_t i = 0; i < cands.size(); i++){
                const Cand& cd = cands[i];
                const std::string cdir = "plots/6_cuts/" + cd.folder + "/";
                gSystem->mkdir(cdir.c_str(), true);
                TH1D hCE("hCandCE", "", cd.nBins, cd.lo, cd.hi), hRecon("hCandRecon", "", cd.nBins, cd.lo, cd.hi);
                for (double x : candCE[i]) hCE.Fill(x);
                for (double x : candRecon[i]) hRecon.Fill(x);
                hCE.SetLineColor(kBlue); hRecon.SetLineColor(kRed); hCE.SetLineWidth(2); hRecon.SetLineWidth(2);
                hRecon.SetTitle((cd.label + ", CE recon protons vs all recon protons, recon theta 0-1 deg, no KE cut;" + cd.label + ";Count").c_str());
                hRecon.SetStats(0); hCE.SetStats(0);
                if (cd.branch == "nHits") for (int bin = 1; bin <= cd.nBins; bin++) hRecon.GetXaxis()->SetBinLabel(bin, bin % 4 == 1 && bin > 1 ? Form("%d", bin - 1) : ""); // labels at 4, 8, 12, 16 hits only
                hRecon.SetMinimum(0.5);
                hRecon.SetMaximum(2.0 * std::max(hRecon.GetMaximum(), hCE.GetMaximum()));
                hRecon.DrawCopy(); hCE.DrawCopy("same");
                TLegend leg(0.725, 0.80, 0.99, 0.90); leg.SetTextSize(0.014);
                leg.AddEntry(&hCE, Form("CE recon protons (%zu)", candCE[i].size()), "l");
                leg.AddEntry(&hRecon, Form("all recon protons (%zu)", candRecon[i].size()), "l");
                leg.Draw();
                if (cd.branch == "nProtonsEv"){ // the 1 and 2 proton cuts (nProtonsEv <= n): lines between bins, kept counts in csv/cut_efficiency.csv
                    const double ymax = hRecon.GetMaximum();
                    TLine l1(1.5, 0.5, 1.5, ymax), l2(2.5, 0.5, 2.5, ymax);
                    for (TLine* l : {&l1, &l2}){l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();}
                    gSystem->mkdir((cdir + "csv").c_str(), true);
                    FILE* nf = fopen((cdir + "csv/cut_efficiency.csv").c_str(), "w");
                    fprintf(nf, "cut,CE_before,CE_after,CE_kept_pct,recon_before,recon_after,recon_left_out_pct\n");
                    for (int nMax : {1, 2}){
                        const size_t ceAfter = std::count_if(candCE[i].begin(), candCE[i].end(), [&](double x){return x <= nMax;});
                        const size_t reAfter = std::count_if(candRecon[i].begin(), candRecon[i].end(), [&](double x){return x <= nMax;});
                        fprintf(nf, "protons in event <= %d,%zu,%zu,%.4f,%zu,%zu,%.4f\n", nMax, candCE[i].size(), ceAfter, 100.0 * ceAfter / candCE[i].size(),
                                candRecon[i].size(), reAfter, 100.0 * (1.0 - double(reAfter) / candRecon[i].size()));
                    }
                    fclose(nf);
                }
                if (cd.branch == "nHits"){ // the cut keeping 4, 8, 12 or 16 tracker hits
                    gSystem->mkdir((cdir + "csv").c_str(), true);
                    FILE* hf = fopen((cdir + "csv/cut_efficiency.csv").c_str(), "w");
                    fprintf(hf, "cut,CE_before,CE_after,CE_kept_pct,recon_before,recon_after,recon_left_out_pct\n");
                    auto inRange = [](double x){return x == 4 || x == 8 || x == 12 || x == 16;};
                    const size_t ceAfter = std::count_if(candCE[i].begin(), candCE[i].end(), inRange);
                    const size_t reAfter = std::count_if(candRecon[i].begin(), candRecon[i].end(), inRange);
                    fprintf(hf, "4, 8, 12 or 16 tracker hits,%zu,%zu,%.4f,%zu,%zu,%.4f\n", candCE[i].size(), ceAfter, 100.0 * ceAfter / candCE[i].size(),
                            candRecon[i].size(), reAfter, 100.0 * (1.0 - double(reAfter) / candRecon[i].size()));
                    fclose(hf);
                }
                c.SaveAs((cdir + (cd.branch == "thetaBlock" ? "2_" : "1_") + cd.branch + ".pdf").c_str());
            }
        }
    }
    if (out){
        fprintf(out, "all four together,,%.3f,%.3f\n", 100.0 * ceKept, 100.0 * reconKept);
        fclose(out);
    }
    printf("  all four together keep %.1f%% CE, %.1f%% recon (%.1f%% of recon removed)\n", 100.0 * ceKept, 100.0 * reconKept, 100.0 * (1.0 - reconKept));
    printf("Paste into QualityCut.h:\n  kMaxResidualMm = %.6g, kMaxResidualDMm = %.6g, kMaxImpactMm = %.6g, kMaxImpactDMm = %.6g\n", th[0], th[1], th[2], th[3]);
    BoxAndCropPlots();
}
