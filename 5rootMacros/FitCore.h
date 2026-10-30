#pragma once
#include <memory>
#include <TF1.h>
#include <TH1.h>

#include <TLegend.h>
#include <TList.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <set>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>
#include <thread>
#include <cstdlib>
#include <sys/stat.h>
#include <TStyle.h>

// Every plotting macro includes this header: plain (not bold) Helvetica for every title, label, stats box, legend and text box.
static const bool kPlainFonts = [](){
    gStyle->SetTextFont(42); gStyle->SetStatFont(42); gStyle->SetLegendFont(42);
    gStyle->SetTitleFont(42, "XYZ"); gStyle->SetLabelFont(42, "XYZ"); gStyle->SetTitleFont(42, ""); // "": the pad title
    return true;
}();

// Beam: the generator (5generator.cc) fires neutrons with 25 GeV/c momentum, so the beam KE is sqrt(p^2 + m^2) - m = 24078 MeV, not 25 GeV.
static const double kNeutronMass = 939.565;   // MeV
static const double kBeamMomentum = 25000.0;  // MeV/c
static const double kBeamKE = std::sqrt(kBeamMomentum * kBeamMomentum + kNeutronMass * kNeutronMass) - kNeutronMass; // MeV

// Transverse-only t in GeV^2: (p_n sin(theta))^2, theta = proton lab polar angle (rad), p_n = lab beam momentum
// (MeV). The full invariant also has an energy term, which the measured final KE (after energy loss) overwhelms.
inline double TFromAngle(double pn, double theta){
    const double pT = pn * std::sin(theta) / 1000.0;
    return pT * pT;
}

// Same formula as the model (G4ChargeExchange.cc, per GeV there), including its minus sign. pbeam in MeV, t in GeV^2.
// AP times apScale (the sim's /CE/APScale), clamped to -1.0 at the low end like the model.
inline double APFromT(double pbeam, double t, double apScale){
    static const double kAPA = 0.7079, kAPk = 0.3860;
    return std::max(-1.0, apScale * -kAPA * (1 - std::exp(-kAPk * pbeam / 1000.0)) * std::sqrt(std::max(t, 0.0)));
}

// Run settings the sim wrote to settings.csv (a row of names, a row of values), as name -> value (.at("APScale"), ...).
// Throws if the file is missing or the two rows do not match.
inline std::map<std::string, double> ReadSettings(const std::string& path = "settings.csv"){
    // memoised by path and mtime+size (called per plot): returns a copy, so callers see no difference
    static std::map<std::string, std::pair<std::pair<long long, long long>, std::map<std::string, double>>> memo;
    struct stat st;
    const bool haveStat = (::stat(path.c_str(), &st) == 0);
    const std::pair<long long, long long> stamp = haveStat ? std::make_pair((long long)st.st_mtime, (long long)st.st_size) : std::make_pair(-1LL, -1LL);
    if (haveStat){
        auto it = memo.find(path);
        if (it != memo.end() && it->second.first == stamp) return it->second.second;
    }
    std::ifstream in(path);
    std::string names, values;
    if (!std::getline(in, names) || !std::getline(in, values))
        throw std::runtime_error("cannot read " + path + " (copy it from the sim's root/<folder>/ next to the ROOT files)");
    std::map<std::string, double> settings;
    std::istringstream nameStream(names), valueStream(values);
    std::string name, value;
    while (std::getline(nameStream, name, ',') && std::getline(valueStream, value, ',')){
        if (!name.empty() && name[0] == '#') name.erase(0, 1); // # marks a default value, the name is the same
        settings[name] = std::stod(value);
    }
    if (settings.empty()) throw std::runtime_error(path + " has no settings");
    if (haveStat) memo[path] = {stamp, settings};
    return settings;
}

// Threads for ROOT::EnableImplicitMT: env ANA_NTHREADS, else the hardware thread count (plot.sh splits the cores between its lanes).
inline unsigned NThreads(){
    if (const char* e = std::getenv("ANA_NTHREADS")){
        const long v = std::atol(e);
        if (v > 0) return (unsigned)v;
    }
    const unsigned h = std::thread::hardware_concurrency();
    return h ? h : 1;
}

// The charge exchange selection of every stage: a CE with its vertex inside one of the carbon blocks (CE in the tracker silicon is left out).
// The block z positions mirror MyDetectorConstruction (2construction.cc): centre of the stack at 4 m, a block after each group of nBlocksTracking tracker planes, dz apart.
inline std::string CarbonCeCut(){
    static std::string memo; static std::map<std::string, double> memoSettings;
    const auto s = ReadSettings();
    if (!memo.empty() && s == memoSettings) return memo;
    memoSettings = s;
    const int nBlocks = (int)std::lround(s.at("nBlocks")), nTracking = (int)std::lround(s.at("nBlocksTracking"));
    const double dzMm = s.at("dz") * 1000.0, halfMm = s.at("CarbonThickness") * 1000.0 / 2.0;
    const int nElements = (nBlocks + 1) * nTracking + nBlocks;
    const double zStart = 4000.0 - (nElements - 1) * dzMm / 2.0;
    std::string windows;
    for (int i = 0; i < nBlocks; i++){
        const double zc = zStart + (i * (nTracking + 1) + nTracking) * dzMm;
        char w[96]; std::snprintf(w, sizeof w, "%s(fZ>%.3f && fZ<%.3f)", windows.empty() ? "" : " || ", zc - halfMm, zc + halfMm);
        windows += w;
    }
    memo = "didChargeExchange>0 && fZ<7000 && (" + windows + ")";
    return memo;
}

// Plot title in the shared "N. Quantity, Label" format. title is "title;xaxis;yaxis". N is the leading digit groups of the
// file name joined by dots (1_1_xy_all.pdf -> 1.1), omitted if the file has none. Spaces around every ';' segment are trimmed.
inline std::string PlotTitle(const std::string& path, const std::string& title){
    std::string file = path.substr(path.find_last_of('/') + 1);
    std::string num;
    for (size_t i = 0; i < file.size() && std::isdigit((unsigned char)file[i]); ){
        size_t j = i;
        while (j < file.size() && std::isdigit((unsigned char)file[j])) j++;
        num += (num.empty() ? "" : ".") + file.substr(i, j - i);
        i = (j + 1 < file.size() && file[j] == '_') ? j + 1 : file.size();
    }
    std::string out = num.empty() ? "" : num + ". ";
    size_t start = 0;
    for (bool first = true; start <= title.size(); first = false){
        size_t end = title.find(';', start);
        if (end == std::string::npos) end = title.size();
        size_t a = title.find_first_not_of(' ', start), b = title.find_last_not_of(' ', end ? end - 1 : 0);
        std::string seg = (a == std::string::npos || a >= end) ? "" : title.substr(a, b - a + 1);
        out += (first ? "" : ";") + seg;
        start = end + 1;
    }
    return out;
}

// Seed from median and IQR (robust to sparse bins and far-away outliers),
// then refit on a window around the fitted mean until it settles.
// A fit drawn over a sub-range stops at the range edge. Adds dotted continuations of f from each edge out to the
// histogram axis limits, so it is visible where the fit curve would keep going. Appended after the fit function, so the
// stats box still reports the real fit.
inline void DrawFitTails(TH1* h, TF1* f){
    double lo = f->GetXmin(), hi = f->GetXmax();
    auto tail = [&](double a, double b, const char* tag){
        if (b <= a) return;
        TF1* t = (TF1*)f->Clone(Form("%s_%s", f->GetName(), tag));
        t->SetRange(a, b);
        t->SetLineStyle(kDotted);
        t->SetLineColor(f->GetLineColor());
        h->GetListOfFunctions()->Add(t);
    };
    tail(h->GetXaxis()->GetXmin(), lo, "tailLo");
    tail(hi, h->GetXaxis()->GetXmax(), "tailHi");
}

// The fit equation in a box in the right margin under the stats box (the fit plots use a 0.30 right margin), added to h's
// function list so it is drawn with h. Parameter names match the stats box.
static const std::string kGausEq = "Constant #times exp(-(x - Mean)^{2} / 2 Sigma^{2})";
// On every plot of all recon protons: the CE protons among them under the same cuts, bottom right below the stats and fit boxes.
// nCEAll: the CE recon protons before any cut, for the share of them the plot's cuts lose (left out when not given), in its own box (DrawCELost).
// nRecon < 0: a CE protons only plot, only the CE lost box.
// nReconAll: every recon proton before any cut, for the share of them in the plot's cut region (left out when not given or nRecon < 0).
// Call it after the legend (and any stats-like info box) is drawn on gPad: the CE lost box is placed under it.
inline void DrawCELost(double nCE, double nCEAll);
inline void DrawCEPercent(double nCE, double nRecon, double nCEAll = -1, double nReconAll = -1){
    const bool keptLine = nRecon >= 0 && nReconAll > 0;
    const int lines = (nRecon >= 0) + keptLine;
    if (nCEAll > 0) DrawCELost(nCE, nCEAll);
    if (lines == 0) return;
    TPaveText* pt = new TPaveText(0.745, 0.02, 0.99, 0.02 + 0.04 * lines, "NDC"); // right of the graph and its tick labels on canvases with a 0.28 or 0.30 right margin
    pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetTextAlign(12); pt->SetTextSize(0.016);
    if (nRecon >= 0) pt->AddText(Form("CE / recon = %.0f / %.0f = %.2f%%", nCE, nRecon, nRecon > 0 ? 100.0 * nCE / nRecon : 0.0));
    if (keptLine) pt->AddText(Form("recon in cut / all = %.0f / %.0f = %.2f%%", nRecon, nReconAll, 100.0 * nRecon / nReconAll));
    pt->Draw();
}

// "CE lost by cuts" in its own one line box right under the legend on gPad, same width and text size; without a legend where the legend would sit:
// under the stats box (the histogram's, or a stats-like text box in the right margin), else under where gStyle puts the stats box.
inline void DrawCELost(double nCE, double nCEAll){
    TLegend* leg = nullptr; TH1* h = nullptr; TPave* info = nullptr;
    for (TObject* o : *gPad->GetListOfPrimitives()){
        if (auto* l = dynamic_cast<TLegend*>(o)) leg = l; // the last legend drawn
        else if (auto* hh = dynamic_cast<TH1*>(o)){if (!h) h = hh;}
        else if (auto* p = dynamic_cast<TPave*>(o)) if (p->GetX1NDC() >= 0.7 && p->GetY1NDC() >= 0.4 && strcmp(p->GetName(), "fitEquation")) info = p; // a fit info box in the right margin (PlotPolarization)
    }
    double x1 = 0.745, x2 = 0.99, top = gStyle->GetStatY() - gStyle->GetStatH(), size = 0.016;
    if (leg){x1 = leg->GetX1NDC(); x2 = leg->GetX2NDC(); top = leg->GetY1NDC(); if (leg->GetTextSize() > 0) size = leg->GetTextSize();}
    else {
        TPave* st = h && !h->TestBit(TH1::kNoStats) ? (TPave*)h->FindObject("stats") : nullptr;
        if (h && !h->TestBit(TH1::kNoStats) && !st){gPad->Modified(); gPad->Update(); st = (TPave*)h->FindObject("stats");} // ROOT makes the stats box when it paints
        if (st) top = st->GetY1NDC();
        else if (info) top = info->GetY1NDC();
    }
    TPaveText* pt = new TPaveText(x1, top - 0.05, x2, top - 0.01, "NDC"); // 0.01 gap: cut_side_box.py tells the frames apart
    pt->SetName("ceLost");
    pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetTextAlign(12); pt->SetTextSize(size);
    pt->AddText(Form("CE lost by cuts = %.0f / %.0f = %.1f%%", nCEAll - nCE, nCEAll, 100.0 * (nCEAll - nCE) / nCEAll));
    pt->Draw();
}

// The fit's stats box goes in the right margin too: ROOT widens a stats box with fit values into the graph, but keeps the place of one that
// already exists, so it is made here at the margin's width. The equation box sits above the CE box of DrawCEPercent (at most 0.02 to 0.10).
inline void AddFitEquation(TH1* h, const std::vector<std::string>& lines){
    if (!h->GetListOfFunctions()->FindObject("stats")){
        TPaveStats* st = new TPaveStats(0.745, 0.58, 0.99, 0.935, "brNDC");
        st->SetName("stats"); st->SetParent(h);
        st->SetOptStat(gStyle->GetOptStat()); st->SetOptFit(gStyle->GetOptFit());
        st->SetBorderSize(1); st->SetFillColor(0); st->SetTextFont(42);
        h->GetListOfFunctions()->Add(st);
    }
    TPaveText* eq = new TPaveText(0.745, 0.17, 0.99, 0.17 + 0.04 * (lines.size() + 1), "NDC");
    eq->SetName("fitEquation");
    eq->SetBorderSize(1); eq->SetFillColor(0); eq->SetTextAlign(12); eq->SetTextSize(0.017);
    eq->AddText("Fit function:");
    for (const auto& l : lines) eq->AddText(l.c_str());
    h->GetListOfFunctions()->Add(eq);
}

// Phi distribution (degrees, 0 to 360) fit with N0 (1 + AP P cos(phi)): AP fixed to ap (the mean AP of the CE protons under the same
// cuts), N0 and the polarization P (limited to -1 to 1) free.
inline void FitCosPhi(TH1* h, double ap){
    if (h->GetEntries() < 10 || std::fabs(ap) < 1e-12) return;
    TF1 f("fCosPhi", "[0]*(1+[1]*[2]*cos(x*TMath::DegToRad()))", 0, 360);
    f.SetParNames("N_{0}", "AP", "P");
    f.SetParameters(h->GetEntries() / h->GetNbinsX(), ap, 0.0);
    f.FixParameter(1, ap);
    f.SetParLimits(2, -1, 1);
    h->Fit(&f, "QL");
    h->SetMinimum(0); h->SetMaximum(1.1 * h->GetMaximum());
    if (h->GetFunction("fCosPhi")) AddFitEquation(h, {"N_{0} (1 + AP #times P #times cos#phi)", "AP fixed: mean AP of the CE protons"});
}

// Block positions (fDx) are discrete cell centers, so h needs one bin per block (bin width = block pitch). A plain
// point-wise gaus fit on such a histogram is meaningless: the "I" option fits the gaussian integrated over each
// block's width, seeded from the histogram mean and RMS. No windowed refit, since a block cannot be split.
inline void FitBlocks(TH1* h){
    if (h->GetEntries() < 3) return;
    double bw = h->GetBinWidth(1);
    TF1 f("fBlocks", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    f.SetParameters(h->GetMaximum(), h->GetMean(), std::max(h->GetRMS(), 0.5 * bw));
    h->Fit(&f, "QLI");
    if (h->GetFunction("fBlocks")) AddFitEquation(h, {kGausEq, "integrated over each block width"});
}

inline void FitCore(TH1* h){
    if (h->GetEntries() < 3) return;

    // Seed from the peak bin and its FWHM, not the quantiles: a flat background (e.g. phi difference, wrapped at
    // +-180) inflates the IQR until the window covers the whole histogram and the fit goes broad.
    double bw = h->GetBinWidth(1);
    // The seed (peak bin, background, FWHM) is read from a lightly smoothed copy: with few entries one noisy bin
    // looked like a one-bin-wide peak and the fit locked onto it. The fit itself still runs on h.
    std::unique_ptr<TH1> hs(static_cast<TH1*>(h->Clone("fitCoreSeed")));
    hs->SetDirectory(nullptr);
    hs->Smooth(2);
    int pk = hs->GetMaximumBin();
    double mean = hs->GetBinCenter(pk);

    // Flat background level = median bin content (robust to the peak). Only when it is a real fraction of the peak
    // (e.g. phi difference wrapped at +-180) is the background fitted; otherwise the plain gaussian is kept.
    std::vector<double> contents;
    for (int b = 1; b <= h->GetNbinsX(); b++) contents.push_back(hs->GetBinContent(b));
    std::nth_element(contents.begin(), contents.begin() + contents.size() / 2, contents.end());
    double bkg = contents[contents.size() / 2];
    bool useBkg = bkg > 0.05 * hs->GetBinContent(pk);
    if (!useBkg) bkg = 0;

    double half = bkg + 0.5 * (hs->GetBinContent(pk) - bkg);
    int lo = pk, hi = pk;
    while (lo > 1 && hs->GetBinContent(lo - 1) > half) lo--;
    while (hi < h->GetNbinsX() && hs->GetBinContent(hi + 1) > half) hi++;
    double sigma = std::max((h->GetBinCenter(hi) - h->GetBinCenter(lo) + bw) / 2.355, bw);

    if (useBkg){
        TF1 f("fCoreBkg", "gaus(0) + pol0(3)", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        f.SetParNames("Constant", "Mean", "Sigma", "Background");
        for (int i = 0; i < 4; i++){
            double halfWin = std::max(3.0 * sigma, 1.5 * bw);
            f.SetParameters(hs->GetBinContent(pk) - bkg, mean, sigma, bkg);
            f.SetParLimits(3, 0, hs->GetBinContent(pk));
            h->Fit(&f, "QL", "", mean - halfWin, mean + halfWin);
            TF1* g = h->GetFunction("fCoreBkg");
            if (!g) return;
            double m = g->GetParameter(1), s = std::fabs(g->GetParameter(2));
            if (s < bw || !std::isfinite(m) || !std::isfinite(s)) return;
            bool settled = std::fabs(m - mean) < 0.05 * sigma && std::fabs(s - sigma) < 0.05 * sigma;
            mean = m; sigma = s;
            if (settled) break;
        }
        TF1* g = h->GetFunction("fCoreBkg");
        if (!g) return;
        // A flat background cannot take up shoulders around the core, so the window fit pushes it up to the shoulder level.
        // When the fit is poor, refit the whole range with a second, wider gaussian for the shoulders: the flat part is then
        // set by the far wings, and the narrow gaussian stays the core (Constant, Mean, Sigma).
        if (g->GetNDF() > 0 && g->GetChisquare() / g->GetNDF() > 3){
            const double c = g->GetParameter(0), m = g->GetParameter(1), s = std::fabs(g->GetParameter(2)), b0 = g->GetParameter(3);
            const double xlo = h->GetXaxis()->GetXmin(), xhi = h->GetXaxis()->GetXmax();
            TF1 f2("fCore2Bkg", "gaus(0) + gaus(3) + pol0(6)", xlo, xhi);
            f2.SetParNames("Constant", "Mean", "Sigma", "Constant2", "Mean2", "Sigma2", "Background");
            f2.SetParameters(c, m, s, 0.3 * c, m, 4 * s, 0.3 * b0);
            f2.SetParLimits(2, bw, 2 * s);
            f2.SetParLimits(4, m - 4 * s, m + 4 * s);
            f2.SetParLimits(5, 1.5 * s, 0.5 * (xhi - xlo));
            f2.SetParLimits(6, 0, b0);
            h->Fit(&f2, "QL");
            if (TF1* g2 = h->GetFunction("fCore2Bkg")) g = g2;
        }
        if (std::string(g->GetName()) == "fCore2Bkg")
            AddFitEquation(h, {kGausEq, "+ Constant2 #times exp(-(x - Mean2)^{2} / 2 Sigma2^{2})", "+ Background"});
        else
            AddFitEquation(h, {kGausEq, "+ Background"});
        DrawFitTails(h, g);
        return;
    }

    for (int i = 0; i < 4; i++){
        double halfWin = std::max(2.5 * sigma, 1.5 * bw);
        h->Fit("gaus", "QL", "", mean - halfWin, mean + halfWin);
        TF1* f = h->GetFunction("gaus");
        if (!f) return;
        double m = f->GetParameter(1), s = std::fabs(f->GetParameter(2));
        if (s < bw || !std::isfinite(m) || !std::isfinite(s)) return;
        bool settled = std::fabs(m - mean) < 0.05 * sigma && std::fabs(s - sigma) < 0.05 * sigma;
        mean = m; sigma = s;
        if (settled) break;
    }
    if (TF1* f = h->GetFunction("gaus")){
        AddFitEquation(h, {kGausEq});
        DrawFitTails(h, f);
    }
}

// The first line of a stats box is the histogram's name: "Fit results" when h has a fit function, else what (or "Statistics").
// Set it with h->SetName right before the canvas is painted and saved, and put the old name back after, so Draw ">>name" and gDirectory->Get still find h.
inline std::string StatsTitle(TH1* h, const std::string& what = ""){
    for (TObject* o : *h->GetListOfFunctions()) if (dynamic_cast<TF1*>(o)) return "Fit results";
    return what.empty() ? "Statistics" : what;
}

// Appends the fit of h (first fit function on it) to fit_values.csv in the folder of plotPath. The first write to a
// file in a run truncates it, so reruns do not pile up.
inline void WriteFitValues(TH1* h, const std::string& plotPath, std::set<std::string>& started){
    TF1* f = nullptr;
    for (TObject* o : *h->GetListOfFunctions()) {
        if ((f = dynamic_cast<TF1*>(o))) break;
    }
    if (!f) return;
    const std::string csv = std::string(gSystem->GetDirName(plotPath.c_str())) + "/fit_values.csv";
    FILE* outCsv = fopen(csv.c_str(), started.insert(csv).second ? "w" : "a");
    if (!outCsv) return;
    if (started.insert(csv + "#header").second) fprintf(outCsv, "plot,function,entries,chi2,ndf,parameter,value,error\n");
    for (int i = 0; i < f->GetNpar(); i++) {
        fprintf(outCsv, "%s,%s,%.0f,%.6g,%d,%s,%.8g,%.4g\n", gSystem->BaseName(plotPath.c_str()), f->GetName(), h->GetEntries(), f->GetChisquare(), f->GetNDF(), f->GetParName(i), f->GetParameter(i), f->GetParError(i));
    }
    fclose(outCsv);
}

// Appends entries, mean and standard deviation of h (x, and y too for 2D histograms) to stats.csv in the folder of
// plotPath. Same truncate-on-first-write rule as WriteFitValues.
inline void WriteHistStats(TH1* h, const std::string& plotPath, std::set<std::string>& started){
    const std::string file = std::string(gSystem->GetDirName(plotPath.c_str())) + "/stats.csv";
    FILE* out = fopen(file.c_str(), started.insert(file).second ? "w" : "a");
    if (!out) return;
    if (started.insert(file + "#header").second) fprintf(out, "plot,entries,mean,sd,mean_y,sd_y\n");
    fprintf(out, "%s,%.0f,%.8g,%.8g", gSystem->BaseName(plotPath.c_str()), h->GetEntries(), h->GetMean(1), h->GetStdDev(1));
    if (h->GetDimension() == 2) fprintf(out, ",%.8g,%.8g\n", h->GetMean(2), h->GetStdDev(2));
    else fprintf(out, ",,\n");
    fclose(out);
}

// Appends one row (the numbers drawn on a graph, not a histogram) to stats.csv in the folder of plotPath, after a header
// row written once per file. Same truncate-on-first-write rule as WriteHistStats; header and row are comma-joined strings.
// name: the csv in the plot's folder, for a second table with its own header next to stats.csv
inline void WriteStatsRow(const std::string& plotPath, const std::string& header, const std::string& row, std::set<std::string>& started, const std::string& name = "stats.csv"){
    const std::string file = std::string(gSystem->GetDirName(plotPath.c_str())) + "/" + name;
    FILE* out = fopen(file.c_str(), started.insert(file).second ? "w" : "a");
    if (!out) return;
    if (started.insert(file + "#header").second) fprintf(out, "%s\n", header.c_str());
    fprintf(out, "%s\n", row.c_str());
    fclose(out);
}

// Every number worth keeping from one saved plot: stats.csv always, fit_values.csv when it has a fit.
inline void WritePlotNumbers(TH1* h, const std::string& plotPath, std::set<std::string>& started){
    WriteHistStats(h, plotPath, started);
    WriteFitValues(h, plotPath, started);
}

// The side box with the cut numbers on every 6_cuts pdf and the crop of every pdf to its content (cut_side_box.py, crop_pdfs.py in the run folder, plot.sh runs them last
// too): called at the end of the macros that write 6_cuts plots, so a single macro run leaves finished plots.
inline void BoxAndCropPlots(){
    if (std::getenv("PLOT_PIPELINE")) return; // plot.sh / Plot.pbs run both scripts after all lanes finished; running them here races with other lanes
    for (const char* script : {"cut_side_box.py", "crop_pdfs.py"})
        if (!gSystem->AccessPathName(script)) gSystem->Exec((std::string("python3 -I ") + script + " plots").c_str());
}
