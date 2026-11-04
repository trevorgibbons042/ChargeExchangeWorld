// KE cut overlays of the megaCut folders, in plots/6_cuts/<6_megaCut|7_megaCut2|9_megaCut3>/1_1particlecut/2_phi_cosfit/<1_Recon|2_ReconD>/KE_overlay/.
// One pdf per curve of Plot3D.cpp part 4 (P vs start theta, its percent error, the ratios 5 to 13), each with the ZDC KE cuts of the folder on top of each other:
//   7_megaCut2, 9_megaCut3: 80% (the folder itself) and 85, 90, 95% (its _KE85, _KE90, _KE95 folders, cut_megacut2_p_ke85.root ... from CutPlots.cpp)
//   6_megaCut: its own sample without a ZDC KE cut, then the 7_megaCut2 curves (6_megaCut plus ZDC KE > x is the megacut2 sample, so those folders hold its 80 to 95% curves)
// Reads the csv files the part 4 plots write next to their pdfs (P_vs_start_theta.csv: plot,sample,start_theta_deg,value,error,pct_kept; ratios_vs_start_theta.csv: plot,start_theta_deg,value),
// so no fit is redone. All recon protons, all recon particles and the CE recon protons of the P plots, and the proton and particle ratios, are separate pdfs because they are separate plots there.
// Run from the run folder after every part 4 and before tidy_plots.py, which moves the csv files (plot.sh, Plot.pbs pdfs step).
#include "TCanvas.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TStyle.h"
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>

namespace {

struct KePt {double x, y, ey;};
using KeCurves = std::map<std::string, std::vector<KePt>>; // "plot name|sample" (sample empty for the ratio csv) -> points in file order

struct KeLevel {std::string label, folder; int color;}; // folder: under plots/6_cuts

std::vector<std::string> KeSplit(const std::string& line){
    std::vector<std::string> f;
    std::stringstream ss(line);
    std::string c;
    while (std::getline(ss, c, ',')) f.push_back(c);
    return f;
}

bool KeExists(const std::string& path){return !gSystem->AccessPathName(path.c_str());}

// the csv of one folder: P_vs_start_theta.csv (with the sample column) or ratios_vs_start_theta.csv
KeCurves KeRead(const std::string& path, bool hasSample){
    KeCurves out;
    std::ifstream in(path);
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line)){
        const auto f = KeSplit(line);
        if (hasSample && f.size() >= 5) out[f[0] + "|" + f[1]].push_back({atof(f[2].c_str()), atof(f[3].c_str()), atof(f[4].c_str())});
        else if (!hasSample && f.size() >= 3) out[f[0] + "|"].push_back({atof(f[1].c_str()), atof(f[2].c_str()), 0});
    }
    return out;
}

std::string KeReplace(std::string s, char from, char to){for (auto& c : s) if (c == from) c = to; return s;}

// "5_phiReD_purity_vs_start_theta.pdf" -> "purity"; "2_phiReD_P_vs_start_theta.pdf" -> "P"
std::string KeQuantity(const std::string& plot){
    std::string q = plot.substr(0, plot.rfind(".pdf"));
    q = q.substr(q.find('_') + 1);                       // the plot number
    q = q.substr(q.find('_') + 1);                       // phiRe or phiReD
    if (const size_t k = q.find("_vs_start_theta"); k != std::string::npos) q = q.substr(0, k) + q.substr(k + 15);
    for (const char* suffix : {"_CE", "_proton", "_particle"}){ // the P scans of one sample each (Plot3D.cpp part 4)
        const size_t n = strlen(suffix);
        if (q.size() > n && q.compare(q.size() - n, n, suffix) == 0){q = q.substr(0, q.size() - n); break;}
    }
    return q;
}

// the true P of the run (PPValue of settings.csv), the reference of the percent error plots
double KePTrue(){
    std::ifstream in("settings.csv");
    std::string names, values, a, b;
    std::getline(in, names); std::getline(in, values);
    std::stringstream sn(names), sv(values);
    while (std::getline(sn, a, ',') && std::getline(sv, b, ',')) if (a == "PPValue") return std::atof(b.c_str());
    return 1.0;
}

void KeOverlayFolder(const std::string& family, const std::string& name, const std::string& sub){
    const std::string base = "plots/6_cuts/";
    const bool ke6 = family == "6_megaCut";
    const std::string fam = ke6 ? "7_megaCut2" : family; // the folders that hold the KE levels
    std::vector<KeLevel> levels;
    if (ke6) levels.push_back({"no ZDC KE cut", family + "/" + name, kBlack});
    levels.push_back({"ZDC KE > 80%", fam + "/" + name, kRed});
    const int colors[3] = {kBlue, kGreen + 2, kMagenta + 1}, pcts[3] = {85, 90, 95};
    for (int i = 0; i < 3; i++) levels.push_back({Form("ZDC KE > %d%%", pcts[i]), fam + "/" + name + Form("_KE%d", pcts[i]), colors[i]});

    const std::string csvDir = "/2_phi_cosfit/" + sub + "/csv/";
    const std::string outDir = sub == "2_ReconD" ? base + family + "/KE_overlay/" : base + family + "/" + name + "/2_phi_cosfit/" + sub + "/KE_overlay/"; // the recon block overlays sit in the megaCut folder itself
    const bool recon = sub == "1_Recon";
    // description of the cuts every curve applies, for the plot (every cut a plot applies is written on it)
    std::vector<std::string> cutLines = {Form("%s/%s:", family.c_str(), name.c_str()),
        name[2] == '1' ? "at most 1 recon particle in the event" : "at most 2 recon particles in the event",
        "2, 4, 6 or 8 tracker hits", "recon theta, block theta 0 to 1 deg"};
    if (family != "6_megaCut") cutLines.push_back("a recon ZDC proton, ZDC KE cut of the legend");
    if (family == "9_megaCut3"){cutLines.push_back("2, 4, 6 or 8 charged tracker hits"); cutLines.push_back("exactly 1 charged ZDC particle");}
    if (ke6) cutLines.push_back("the other curves also have the ZDC KE cut");

    for (const bool isP : {true, false}){
        std::vector<KeCurves> data;
        std::vector<const KeLevel*> have;
        for (const auto& lv : levels){
            const std::string path = base + lv.folder + csvDir + (isP ? "P_vs_start_theta.csv" : "ratios_vs_start_theta.csv");
            if (!KeExists(path)) continue;
            data.push_back(KeRead(path, isP));
            have.push_back(&lv);
        }
        if (data.size() < 2) continue; // an overlay needs at least two curves
        std::map<std::string, int> keys; // every curve the first level has
        for (const auto& [k, pts] : data[0]) keys[k] = int(pts.size());
        for (const auto& [key, n] : keys){
            (void)n;
            const std::string plot = key.substr(0, key.find('|')), sample = key.substr(key.find('|') + 1);
            const bool pct = plot.find("pct_error") != std::string::npos, ce = plot.size() > 7 && plot.compare(plot.size() - 7, 7, "_CE.pdf") == 0;
            const std::string quantity = KeQuantity(plot);
            auto* mg = new TMultiGraph();
            auto* leg = new TLegend(0.705, 0.62, 0.99, 0.90);
            leg->SetTextSize(0.016); leg->SetBorderSize(0);
            int drawn = 0;
            for (size_t i = 0; i < data.size(); i++){
                const auto it = data[i].find(key);
                if (it == data[i].end() || it->second.empty()) continue;
                auto* g = new TGraphAsymmErrors();
                for (const auto& p : it->second){ // the error bars of P stop at -1 and 1, the physical range, as in the single plots of Plot3D.cpp
                    g->SetPoint(g->GetN(), p.x, p.y);
                    const bool bars = isP && !pct;
                    g->SetPointError(g->GetN() - 1, 0, 0, bars ? p.y - std::max(p.y - p.ey, -1.0) : 0, bars ? std::min(p.y + p.ey, 1.0) - p.y : 0);
                }
                if (isP && !pct) g->SetLineColorAlpha(have[i]->color, 0.3); else g->SetLineColor(have[i]->color); // the error bars of the four curves overlap: see-through
                g->SetMarkerColor(have[i]->color); g->SetMarkerStyle(20); g->SetMarkerSize(isP && !pct ? 1.0 : 0.5); g->SetLineWidth(1);
                mg->Add(g, isP && !pct ? "P" : "PL");
                leg->AddEntry(g, have[i]->label.c_str(), "lp");
                drawn++;
            }
            if (drawn < 2){delete mg; delete leg; continue;}
            const std::string what = (isP ? std::string(pct ? "P error (%)" : "P") + (ce ? ", CE recon protons" : ", " + sample) : quantity);
            auto* canvas = new TCanvas("cKEov", "", 800, 640);
            canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
            mg->SetTitle((what + " vs start theta, ZDC KE cut overlay, " + family + "/" + name + ", " + (recon ? "recon" : "recon block") + " theta from the start value to 1 deg;Start theta (deg);" + what).c_str());
            if (isP){ // P is between -1 and 1 (the percent error plots show the same range)
                const double pTrue = KePTrue();
                mg->SetMinimum(pct ? 100 * (-1.05 - pTrue) / pTrue : -1.05); mg->SetMaximum(pct ? 100 * (1.05 - pTrue) / pTrue : 1.05);
            }
            mg->Draw("A");
            mg->GetXaxis()->SetLimits(0, 1.0);
            leg->Draw();
            auto* box = new TPaveText(0.705, 0.20, 0.99, 0.60, "NDC");
            box->SetTextSize(0.014); box->SetTextAlign(12); box->SetBorderSize(0); box->SetFillStyle(0);
            for (const auto& l : cutLines) box->AddText(l.c_str());
            box->Draw();
            gSystem->mkdir(outDir.c_str(), true);
            const std::string file = plot.substr(0, plot.rfind(".pdf")); // the P plots name their sample (_CE, _proton, _particle)
            canvas->SaveAs((outDir + file + ".pdf").c_str());
            delete canvas;
        }
    }
}

} // namespace

void PlotKEOverlay(){
    gStyle->SetOptStat(0);
    for (const char* family : {"7_megaCut2", "9_megaCut3"}) // 6_megaCut has no ZDC KE cut, so no overlay
        for (const char* name : {"1_1particlecut"})
            for (const char* sub : {"2_ReconD"}) // recon block only, no overlay for 1_Recon
                KeOverlayFolder(family, name, sub);
}
