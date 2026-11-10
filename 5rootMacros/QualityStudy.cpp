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
// Writes plots/6_cuts/3_ResCuts/thresholds.csv (no pdf of the residual). Prints the numbers to paste into QualityCut.h.
static double Quantile(std::vector<double> v, double q){
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, (size_t)(q * v.size()))];
}

void QualityStudy(double keepCE = 0.95){
    FastJit();
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    gStyle->SetPalette(kViridis);
    // cuts/cut_base.root (CutPlots.cpp samples) holds only the rows of the theta baseline, which is all this macro uses; used when it is newer than cut_angles_all.root, else every fitted proton is read
    Long_t idB, flagsB, mtBase = 0, mtAll = 1; Long64_t sizeB;
    const bool haveBase = !gSystem->GetPathInfo("cuts/cut_base.root", &idB, &sizeB, &flagsB, &mtBase) && !gSystem->GetPathInfo("cuts/cut_angles_all.root", &idB, &sizeB, &flagsB, &mtAll) && mtBase >= mtAll;
    TFile* f = TFile::Open(haveBase ? "cuts/cut_base.root" : "cuts/cut_angles_all.root", "READ");
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
    if (hasCand){iCand = (int)want.size(); for (const char* b : {"nProtonsEv", "nHits", "firstXD", "firstYD", "resHalfRatioD", "thetaRelN", "thetaRelND"}) want.push_back(b);}
    const std::vector<std::vector<double>> col = ColumnsUnordered(t, want); // row order does not matter here, and AnglesRaw (many clusters) is then read on all threads
    const size_t nRows = col[0].size();
    auto isCE = [&](size_t e){return col[4][e] != 0;};
    auto inBase = [&](size_t e){ // baseline of 6_cuts: RECONTHETA_ESTNEU and RECONTHETABL_ESTNEU from 0 to 2 degrees (NaN fails)
        if (!hasCand) return true;
        const double deg = col[iCand + 5][e] * TMath::RadToDeg(), degBlock = col[iCand + 6][e] * TMath::RadToDeg();
        return deg >= 0 && deg <= 2 && degBlock >= 0 && degBlock <= 2;
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
    // sorted copies, made once: the bisection below asks for the thresholds 40 times and Quantile would sort a copy of every variable each time.
    // ce itself stays in row order, jointKept needs the four values of one track together.
    std::vector<std::vector<double>> ceSorted = ce;
    for (auto& v : ceSorted) std::sort(v.begin(), v.end());
    auto thresholds = [&](double q){
        std::vector<double> th(4);
        for (int i = 0; i < 4; i++) th[i] = ceSorted[i].empty() ? NAN : ceSorted[i][std::min(ceSorted[i].size() - 1, (size_t)(q * ceSorted[i].size()))]; // same index as Quantile
        return th;
    };
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

    }

    // residual_cut_efficiency.csv: study only, no plot applies a residual cut. CE and recon protons kept at three block residual thresholds,
    // totals before and after.
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
        fclose(rc);
    }

    if (out){
        fprintf(out, "all four together,,%.3f,%.3f\n", 100.0 * ceKept, 100.0 * reconKept);
        fclose(out);
    }
    printf("  all four together keep %.1f%% CE, %.1f%% recon (%.1f%% of recon removed)\n", 100.0 * ceKept, 100.0 * reconKept, 100.0 * (1.0 - reconKept));
    BoxAndCropPlots();
}
