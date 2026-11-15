#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1F.h>
#include "FitCore.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

// Block scan comparison, one pdf per quantity with the number of blocks on x. Run it on the Mac, outside a run folder:
//   root -l -b -q 'PlotBlockScan.cpp("build/root/<name>_<total>m/25GeV", "6_megaCut")'   (the default scanDir "25GeV" is for running inside build/root/<name>_<total>m; 6jobs/ScanBlocks.sh makes build/root/<name>_<total>m/<E>GeV/<n>block, one run folder per block count)
// The scan folder (one <E>GeV folder) holds <n>block/6_ReconNew_1/<folder>/scan_stats.csv (the file F9_PlotReconNew.cpp writes in every folder of plots/6_ReconNew_1), or a copy of it as <n>block/<folder>/scan_stats.csv, n = 1 to 40.
// <folder> is one of 1_noCuts, 2_KECUT80, 3_nPlanesUsedEven, 4_Theta0_1, 5_avgTrackerHitsPerPlane1, 6_megaCut (the default); the cuts of it are in every title.
// Every number is of the event recon (cuts/cut_angles_all.root). Signal = events with a CE hit (CE + CE mix), background = non-CE events. S, B are the counts that pass the cuts of the folder, S_all, B_all the counts with no cut.
// The polarization has one pdf per recon sample (CE + CE mix, all recon particles). Points only, error bars only on the curves that have an error.
// Output: <scan folder>/block_scan_plots/<folder>/<n>_<name>.pdf. The pdfs already in that folder are removed first, the numbers can shift.
namespace {
struct ScanRow {
    int nBlocks = 0;
    double neutrons = 0, pTrue = 0, nCETrue = 0, sAll = 0, bAll = 0, s = 0, b = 0, meanA = 0, pMix = 0, eMix = 0, pAll = 0, eAll = 0, pCE = NAN, eCE = NAN;
};
std::string CutText(const std::string& folder){
    if (folder == "1_noCuts") return "no cuts";
    if (folder == "2_KECUT80") return "KE above 80% of the beam KE";
    if (folder == "3_nPlanesUsedEven") return "planes in the fit even";
    if (folder == "4_Theta0_1") return "theta 0 to 1 degree";
    if (folder == "5_avgTrackerHitsPerPlane1") return "average 1 hit per plane";
    return "KE above 80% of the beam KE, planes in the fit even, theta 0 to 1 degree, average 1 hit per plane";
}

bool ReadRow(const std::string& path, ScanRow& r){
    std::ifstream in(path);
    std::string head, line;
    if (!in || !std::getline(in, head) || !std::getline(in, line)) return false;
    std::vector<double> v;
    std::stringstream ss(line); std::string f;
    while (std::getline(ss, f, ',')) v.push_back(std::stod(f));
    if (v.size() < 15) return false;
    r.nBlocks = (int)std::lround(v[0]); r.neutrons = v[1]; r.pTrue = v[2]; r.nCETrue = v[3]; r.sAll = v[4]; r.bAll = v[5]; r.s = v[6]; r.b = v[7]; r.meanA = v[8];
    r.pMix = v[10]; r.eMix = v[11]; r.pAll = v[13]; r.eAll = v[14]; // v[9] and v[12] are the fit sizes nFitCEMix and nFitAll, not plotted
    if (v.size() >= 18){r.pCE = v[16]; r.eCE = v[17];} // nFitCE, pCE, errCE: scan_stats.csv of F9_PlotReconNew.cpp from the CE columns on
    return true;
}
double Binom(double k, double n){return n > 0 ? std::sqrt(std::max(0.0, (k / n) * (1.0 - k / n) / n)) : 0.0;}

struct Curve {
    const char* label; int color; int marker;
    std::vector<double> x, y, ey;
};

void Draw(const std::string& dir, int number, const char* name, const char* title, const char* yTitle, const std::string& cuts, const std::vector<Curve>& curves, double refLine, int nMax, double yLo = NAN, double yHi = NAN){ // yLo, yHi: a fixed y axis range, so two plots can share one
    TCanvas c("cBlockScan", "", 900, 650);
    c.SetLeftMargin(0.13); c.SetRightMargin(0.04); c.SetTopMargin(0.14);
    double lo = 1e300, hi = -1e300;
    for (const Curve& k : curves) for (size_t i = 0; i < k.y.size(); i++){
        if (!std::isfinite(k.y[i])) continue;
        lo = std::min(lo, k.y[i] - k.ey[i]); hi = std::max(hi, k.y[i] + k.ey[i]);
    }
    if (refLine == refLine){lo = std::min(lo, refLine); hi = std::max(hi, refLine);}
    if (!(hi > lo)){lo -= 0.5; hi += 0.5;}
    const double pad = 0.12 * (hi - lo);
    TH1F frame("frameBlockScan", "", 1, 0.5, nMax + 0.5);
    frame.SetMinimum(yLo == yLo ? yLo : lo - pad); frame.SetMaximum(yHi == yHi ? yHi : hi + 2.0 * pad);
    frame.SetStats(0);
    frame.SetTitle(Form(";Number of blocks;%s", yTitle));
    frame.GetXaxis()->SetNdivisions(nMax, "N");
    frame.Draw("axis");
    // the "axis" draw option does not paint the histogram title, so the title is written by hand above the frame
    TLatex text; text.SetNDC(); text.SetTextFont(42);
    text.SetTextSize(0.036); text.DrawLatex(0.13, 0.955, title);
    text.SetTextSize(0.026); text.DrawLatex(0.13, 0.915, Form("Selection: %s", cuts.c_str()));
    if (refLine == refLine){TLine l(0.5, refLine, nMax + 0.5, refLine); l.SetLineStyle(2); l.SetLineColor(kGray + 2); l.DrawClone();}
    TLegend leg(0.50, 0.74, 0.95, 0.84);
    leg.SetBorderSize(0); leg.SetFillStyle(0);
    std::vector<TGraph*> gs;
    for (const Curve& k : curves){
        // error bars only where the curve has an error, points only, no connecting line
        const bool hasErr = std::any_of(k.ey.begin(), k.ey.end(), [](double e){return e > 0;});
        TGraph* g = hasErr ? new TGraphErrors((int)k.x.size()) : new TGraph((int)k.x.size());
        for (size_t i = 0; i < k.x.size(); i++){
            g->SetPoint((int)i, k.x[i], k.y[i]);
            if (hasErr) ((TGraphErrors*)g)->SetPointError((int)i, 0, k.ey[i]);
        }
        g->SetMarkerStyle(k.marker); g->SetMarkerColor(k.color); g->SetLineColor(k.color); g->SetMarkerSize(1.0);
        g->Draw("P same");
        leg.AddEntry(g, k.label, hasErr ? "ep" : "p");
        gs.push_back(g);
    }
    if (refLine == refLine) leg.AddEntry((TObject*)nullptr, Form("dashed line, true P = %.2f", refLine), "");
    if (curves.size() > 1 || refLine == refLine) leg.Draw();
    c.SaveAs(Form("%s/%d_%s.pdf", dir.c_str(), number, name));
    for (TGraph* g : gs) delete g;
}
}

void PlotBlockScan(const char* scanDir = "25GeV", const char* folder = "6_megaCut"){
    const std::string base = scanDir, fold = folder;
    std::vector<ScanRow> rows;
    for (int n = 1; n <= 40; n++){
        ScanRow r;
        // either the copied layout <n>block/<folder>/ or the full plots layout <n>block/6_ReconNew_1/<folder>/
        if (ReadRow(Form("%s/%dblock/%s/scan_stats.csv", base.c_str(), n, fold.c_str()), r) ||
            ReadRow(Form("%s/%dblock/6_ReconNew_1/%s/scan_stats.csv", base.c_str(), n, fold.c_str()), r)) rows.push_back(r);
    }
    if (rows.empty()){Error("PlotBlockScan", "no %s/<n>block/%s/scan_stats.csv", base.c_str(), fold.c_str()); return;}
    const std::string out = base + "/block_scan_plots/" + fold, cuts = CutText(fold);
    gSystem->mkdir(out.c_str(), true);
    // the numbers of the pdfs can shift between versions, so the old pdfs go first
    if (void* d = gSystem->OpenDirectory(out.c_str())){
        while (const char* e = gSystem->GetDirEntry(d)){
            const std::string f = e;
            if (f.size() > 4 && f.compare(f.size() - 4, 4, ".pdf") == 0) gSystem->Unlink((out + "/" + f).c_str());
        }
        gSystem->FreeDirectory(d);
    }
    gStyle->SetOptStat(0);
    int nMax = 0;
    for (const ScanRow& r : rows) nMax = std::max(nMax, r.nBlocks);

    // one curve of the per-run values, y and its error come from the two functions
    const auto make = [&](const char* label, int color, int marker, const std::function<double(const ScanRow&)>& y, const std::function<double(const ScanRow&)>& ey){
        Curve k; k.label = label; k.color = color; k.marker = marker;
        for (const ScanRow& r : rows){k.x.push_back(r.nBlocks); k.y.push_back(y(r)); k.ey.push_back(ey(r));}
        return k;
    };
    const auto zero = [](const ScanRow&){return 0.0;};
    int number = 0;
    const double pTrue = rows.front().pTrue;
    // CE only fit (the actual CE angle): the pCE columns of scan_stats.csv, or for runs without them <scan folder>/block_scan_plots/ce_only_from_pdf.csv (cut,nBlocks,nFitCE,pCE,errCE, read from 1_phiBlock_cosfit_CE.pdf)
    Curve ceCurve; ceCurve.label = "CE recon protons"; ceCurve.color = CbCE(); ceCurve.marker = 22;
    {
        std::ifstream in(base + "/block_scan_plots/ce_only_from_pdf.csv");
        std::string line;
        std::getline(in, line);
        while (std::getline(in, line)){
            std::vector<std::string> f; std::stringstream ss(line); std::string w;
            while (std::getline(ss, w, ',')) f.push_back(w);
            if (f.size() < 5 || f[0] != fold) continue;
            ceCurve.x.push_back(std::stod(f[1])); ceCurve.y.push_back(std::stod(f[3])); ceCurve.ey.push_back(std::stod(f[4]));
        }
        if (ceCurve.x.empty()) for (const ScanRow& r : rows) if (std::isfinite(r.pCE)){ceCurve.x.push_back(r.nBlocks); ceCurve.y.push_back(r.pCE); ceCurve.ey.push_back(r.eCE);}
    }
    // the polarization plots share one y range: the span of the samples with their errors and the true P, padded like Draw
    double pLo = pTrue, pHi = pTrue;
    for (const ScanRow& r : rows){
        for (const double v : {r.pMix - r.eMix, r.pMix + r.eMix, r.pAll - r.eAll, r.pAll + r.eAll}) if (std::isfinite(v)){pLo = std::min(pLo, v); pHi = std::max(pHi, v);}
    }
    for (size_t i = 0; i < ceCurve.x.size(); i++) for (const double v : {ceCurve.y[i] - ceCurve.ey[i], ceCurve.y[i] + ceCurve.ey[i]}) if (std::isfinite(v)){pLo = std::min(pLo, v); pHi = std::max(pHi, v);}
    const double pPad = 0.12 * (pHi - pLo);
    const auto drawP = [&](const char* name, const char* title, const Curve& curve){Draw(out, ++number, name, title, "Fitted polarization P", cuts, {curve}, pTrue, nMax, pLo - pPad, pHi + 2.0 * pPad);};
    drawP("polarization_CE_mix", "Fitted polarization vs. number of blocks, CE + CE mix recon",
          make("CE + CE mix recon", CbProton(), 21, [](const ScanRow& r){return r.pMix;}, [](const ScanRow& r){return r.eMix;}));
    drawP("polarization_all_recon", "Fitted polarization vs. number of blocks, all recon particles",
          make("all recon particles", CbAll(), 20, [](const ScanRow& r){return r.pAll;}, [](const ScanRow& r){return r.eAll;}));
    if (!ceCurve.x.empty()) drawP("polarization_CE", "Fitted polarization vs. number of blocks, CE recon protons at the actual CE angle", ceCurve);
    // FOM = (signal in the cuts / neutrons) x mean(A_eff)^2, the mean taken over the signal events in the cuts
    Draw(out, ++number, "FOM", "Figure of merit vs. number of blocks", "FOM = (S / N_{neutrons}) #times <A_{eff}>^{2}", cuts, {
        make("signal in the cuts", CbCE(), 20, [](const ScanRow& r){return r.s / r.neutrons * r.meanA * r.meanA;}, [](const ScanRow& r){return std::sqrt(r.s) / r.neutrons * r.meanA * r.meanA;})}, NAN, nMax);
    Draw(out, ++number, "CE_over_neutrons", "CE yield per incident neutron vs. number of blocks", "N_{CE} / N_{neutrons}", "none (every CE in the carbon blocks)", {
        make("every CE in the carbon blocks", CbCE(), 20, [](const ScanRow& r){return r.nCETrue / r.neutrons;}, zero)}, NAN, nMax);
    Draw(out, ++number, "CErecon_in_cuts_over_neutrons", "Reconstructed CE yield per incident neutron vs. number of blocks", "N_{CE recon, in the cuts} / N_{neutrons}", cuts, {
        make("CE + CE mix recon events in the cuts", CbProton(), 20, [](const ScanRow& r){return r.s / r.neutrons;}, [](const ScanRow& r){return std::sqrt(r.s) / r.neutrons;})}, NAN, nMax);
    Draw(out, ++number, "purity", "Signal purity vs. number of blocks", "Purity = S / (S + B)", cuts, {
        make("signal in the cuts", CbCE(), 20, [](const ScanRow& r){return r.s / (r.s + r.b);}, [](const ScanRow& r){return Binom(r.s, r.s + r.b);})}, NAN, nMax);
    Draw(out, ++number, "contamination", "Background contamination vs. number of blocks", "Contamination = B / (S + B)", cuts, {
        make("non-CE events in the cuts", CbAll(), 20, [](const ScanRow& r){return r.b / (r.s + r.b);}, [](const ScanRow& r){return Binom(r.b, r.s + r.b);})}, NAN, nMax);
    Draw(out, ++number, "signal_efficiency", "Signal efficiency vs. number of blocks", "Signal efficiency = S / S_{all}", cuts, {
        make("CE + CE mix recon events", CbCE(), 20, [](const ScanRow& r){return r.s / r.sAll;}, [](const ScanRow& r){return Binom(r.s, r.sAll);})}, NAN, nMax);
    Draw(out, ++number, "signal_significance", "Signal significance vs. number of blocks", "Signal significance = S / #sqrt{S + B}", cuts, {
        make("signal in the cuts", CbCE(), 20, [](const ScanRow& r){return r.s / std::sqrt(r.s + r.b);}, zero)}, NAN, nMax);
    Draw(out, ++number, "background_efficiency", "Background efficiency vs. number of blocks", "Background efficiency = B / B_{all}", cuts, {
        make("non-CE events", CbAll(), 20, [](const ScanRow& r){return r.b / r.bAll;}, [](const ScanRow& r){return Binom(r.b, r.bAll);})}, NAN, nMax);
    Draw(out, ++number, "background_rejection", "Background rejection vs. number of blocks", "Background rejection = 1 - B / B_{all}", cuts, {
        make("non-CE events", CbAll(), 20, [](const ScanRow& r){return 1.0 - r.b / r.bAll;}, [](const ScanRow& r){return Binom(r.b, r.bAll);})}, NAN, nMax);
    printf(">>> PlotBlockScan: %d plots in %s (%zu runs)\n", number, out.c_str(), rows.size());
}
