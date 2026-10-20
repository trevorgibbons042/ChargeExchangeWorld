#pragma once
#include <memory>
#include <TF1.h>
#include <TH1.h>

#include <TList.h>
#include <TSystem.h>

#include <algorithm>
#include <cstdio>
#include <set>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>

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

// Block positions (fDx) are discrete cell centers, so h needs one bin per block (bin width = block pitch). A plain
// point-wise gaus fit on such a histogram is meaningless: the "I" option fits the gaussian integrated over each
// block's width, seeded from the histogram mean and RMS. No windowed refit, since a block cannot be split.
inline void FitBlocks(TH1* h){
    if (h->GetEntries() < 3) return;
    double bw = h->GetBinWidth(1);
    TF1 f("fBlocks", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    f.SetParameters(h->GetMaximum(), h->GetMean(), std::max(h->GetRMS(), 0.5 * bw));
    h->Fit(&f, "QLI");
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
        if (TF1* g = h->GetFunction("fCoreBkg")) DrawFitTails(h, g);
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
    if (TF1* f = h->GetFunction("gaus")) DrawFitTails(h, f);
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

// Every number worth keeping from one saved plot: stats.csv always, fit_values.csv when it has a fit.
inline void WritePlotNumbers(TH1* h, const std::string& plotPath, std::set<std::string>& started){
    WriteHistStats(h, plotPath, started);
    WriteFitValues(h, plotPath, started);
}
