#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <set>
#include <vector>

#include "FitCore.h"
#include "FomCommon.h"

// Fraction of CE protons among all protons per theta, N_CE / N, one curve per KE cut. Reads the 1_csv/ke_*.csv files written by PlotAngle.cpp
// and writes the 1_* ce_fraction_vs_theta plots and the 2_* cumulative_ce_percent plots (log y axis only; the numbers 1, 3 of each group are the unused linear plots) into the 2_plots folder beside the csv folder. The csv bins are added into 0.1 deg bins (full range) or 0.05 deg bins (0.3 to 2 deg), each plotted as a point at its center.
void PlotCEFraction(){
    gStyle->SetOptStat(0);
    const auto &energies = kEnergies, &sets = kSets, &kes = kKes;
    const auto &colors = kColors;
    const auto &energyName = kEnergyName, &setName = kSetName;
    const auto keName = KeName;

    const double thetaMax = ReadThetaMax();

    // full theta range in 0.1 deg bins, and the 0.3 to 2 deg integrating region in 0.05 deg bins
    struct Variant {double binWidth, xMin, xMax; std::string tag, label; double minN;}; // minN: bins with fewer protons are left out (too noisy)
    std::vector<Variant> variants = {{0.1, 0.0, thetaMax, "", "0.1", 20}, {0.05, 0.3, 2.0, "_0p3to2", "0.05", 1}};

    // variants with a minN: the x range ends at the last bin that holds at least minN protons for any plot, so every plot of the variant shares it
    auto readCsv = [](const std::string& path){
        std::vector<std::array<double, 4>> rows; // theta_lo, theta_hi, N, N_CE
        std::ifstream in(path);
        std::string line; std::getline(in, line);
        while (std::getline(in, line)){
            std::vector<double> x; std::stringstream ss(line); std::string f;
            while (std::getline(ss, f, ',')) x.push_back(std::stod(f));
            if (x.size() < 4 || x[1] - x[0] > 0.5) continue; // skip the total row
            rows.push_back({x[0], x[1], x[2], x[3]});
        }
        return rows;
    };
    for (auto& v : variants){
        if (v.minN <= 1) continue;
        const int nAll = (int)std::lround((v.xMax - v.xMin) / v.binWidth);
        int lastBin = 0;
        for (const auto& en : energies) for (const auto& set : sets) for (const auto& ke : kes){
            std::vector<double> sumN(nAll, 0.0);
            for (const auto& r : readCsv("plots/6_fomCalc/" + en + "/" + set + "/1_csv/" + ke + ".csv")){
                const int bin = (int)std::floor((0.5 * (r[0] + r[1]) - v.xMin) / v.binWidth);
                if (bin >= 0 && bin < nAll) sumN[bin] += r[2];
            }
            for (int i = 0; i < nAll; i++) if (sumN[i] >= v.minN) lastBin = std::max(lastBin, i);
        }
        v.xMax = v.xMin + (lastBin + 1) * v.binWidth;
    }

    // numbers behind each plot go to 2_plots/stats.csv, shared with PlotFomCumulative.cpp: this macro replaces only its own rows (plot names starting "1_" or "2_")
    std::set<std::string> statsStarted;
    const std::string statsHeader = "plot,ke_cut,N,N_CE,pct_CE,fom";
    for (const auto& en : energies) for (const auto& set : sets){
        const std::string file = "plots/6_fomCalc/" + en + "/" + set + "/2_plots/stats.csv";
        std::ifstream in(file);
        std::string line, kept;
        while (std::getline(in, line)) if (!line.empty() && line[0] >= '3' && line != statsHeader) kept += line + "\n";
        in.close();
        if (kept.empty()) continue;
        gSystem->mkdir(gSystem->GetDirName(file.c_str()), true);
        std::ofstream out(file);
        out << statsHeader << "\n" << kept;
        statsStarted.insert(file); statsStarted.insert(file + "#header");
    }
    auto statsRow = [&](const std::string& plotPath, const std::string& ke, double n, double nCE){
        WriteStatsRow(plotPath, statsHeader, Form("%s,%s,%.0f,%.0f,%.6g,", gSystem->BaseName(plotPath.c_str()), ke.c_str(), n, nCE, n > 0 ? 100.0 * nCE / n : 0.0), statsStarted);
    };

    TCanvas c("c", "", 1200, 900);
    c.SetLeftMargin(0.14); c.SetRightMargin(0.2); // room for the legend outside the graph
    for (size_t vi = 0; vi < variants.size(); vi++) for (const auto& en : energies) for (const auto& set : sets){
        const auto& v = variants[vi];
        const std::string dir = "plots/6_fomCalc/" + en + "/" + set + "/";
        gSystem->mkdir((dir + "2_plots").c_str(), true);
        TMultiGraph mg;
        TLegend leg(0.81, 0.74, 0.99, 0.88);
        leg.SetBorderSize(0); leg.SetTextSize(0.022);
        const int nBins = (int)std::lround((v.xMax - v.xMin) / v.binWidth);
        for (size_t k = 0; k < kes.size(); k++){
            std::ifstream in(dir + "1_csv/" + kes[k] + ".csv");
            if (!in) continue;
            std::string line; std::getline(in, line);
            // csv bins have unequal widths: add each one to the fixed width bin holding its center
            std::vector<double> sumN(nBins, 0.0), sumCE(nBins, 0.0);
            while (std::getline(in, line)){
                std::vector<double> x; std::stringstream ss(line); std::string f;
                while (std::getline(ss, f, ',')) x.push_back(std::stod(f));
                if (x.size() < 4 || x[1] - x[0] > 0.5) continue; // skip the total row
                const int bin = (int)std::floor((0.5 * (x[0] + x[1]) - v.xMin) / v.binWidth);
                if (bin < 0 || bin >= nBins) continue;
                sumN[bin] += x[2]; sumCE[bin] += x[3];
            }
            TGraph* g = new TGraph();
            double plottedN = 0, plottedCE = 0; // protons in the bins that made it onto the plot
            for (int i = 0; i < nBins; i++) if (sumN[i] >= v.minN){
                g->SetPoint(g->GetN(), v.xMin + (i + 0.5) * v.binWidth, sumCE[i] / sumN[i]);
                plottedN += sumN[i]; plottedCE += sumCE[i];
            }
            statsRow(dir + "2_plots/1_" + std::to_string(2 * vi + 2) + "_ce_fraction_vs_theta" + v.tag + "_log.pdf", kes[k], plottedN, plottedCE);
            g->SetLineColor(colors[k]); g->SetMarkerColor(colors[k]); g->SetMarkerStyle(20); g->SetMarkerSize(0.4); g->SetLineWidth(1);
            mg.Add(g, "LP");
            leg.AddEntry(g, keName(kes[k]).c_str(), "lp");
        }
        if (!mg.GetListOfGraphs() || mg.GetListOfGraphs()->IsEmpty()) continue;
        c.Clear();
        mg.SetTitle(("CE protons / all protons vs #theta, " + energyName.at(en) + ", " + setName.at(set) +
                     (v.minN > 1 ? Form(" (bins with < %g protons left out)", v.minN) : "") + ";recon #theta (deg), " + v.label + " deg bins;N_{CE} / N_{all protons}").c_str());
        mg.Draw("A");
        mg.GetXaxis()->SetLimits(v.xMin, v.xMax);
        leg.DrawClone();
        // file numbers: 1_ CE fraction per bin (1 full range, 2 its log, 3 window, 4 its log), 2_ cumulative CE percent (same order); PlotFomCumulative.cpp uses 3_ to 5_
        const std::string fracStem = std::string("_ce_fraction_vs_theta") + v.tag;
        SaveLogCopy(c, mg, dir + "2_plots/1_" + std::to_string(2 * vi + 2) + fracStem + "_log.pdf");
    }

    // 2_*: running percent of CE protons, sum of N_CE / sum of N over the bins up to each theta; full range, and only the 0.3 to 2 deg integrating region
    const std::vector<std::array<double, 2>> cumRanges = {{0.0, thetaMax}, {0.3, 2.0}};
    const std::vector<std::string> cumTags = {"", "_0p3to2"};
    const double minCumN = 50; // the running percent starts once this many protons are summed
    for (size_t ci = 0; ci < cumRanges.size(); ci++) for (const auto& en : energies) for (const auto& set : sets){
        const double lo = cumRanges[ci][0], hi = cumRanges[ci][1];
        const std::string dir = "plots/6_fomCalc/" + en + "/" + set + "/";
        TMultiGraph mg;
        TLegend leg(0.81, 0.74, 0.99, 0.88);
        leg.SetBorderSize(0); leg.SetTextSize(0.022);
        for (size_t k = 0; k < kes.size(); k++){
            const auto rows = readCsv(dir + "1_csv/" + kes[k] + ".csv");
            if (rows.empty()) continue;
            TGraph* g = new TGraph();
            double sumN = 0, sumCE = 0;
            for (const auto& r : rows){
                const double center = 0.5 * (r[0] + r[1]);
                if (center < lo || center >= hi) continue;
                sumN += r[2]; sumCE += r[3];
                if (sumN >= minCumN) g->SetPoint(g->GetN(), r[1], 100.0 * sumCE / sumN); // the first few bins hold too few protons for a steady ratio
            }
            statsRow(dir + "2_plots/2_" + std::to_string(2 * ci + 2) + "_cumulative_ce_percent" + cumTags[ci] + "_log.pdf", kes[k], sumN, sumCE); // final point of the curve
            g->SetLineColor(colors[k]); g->SetLineWidth(1);
            mg.Add(g, "L");
            leg.AddEntry(g, keName(kes[k]).c_str(), "l");
        }
        if (!mg.GetListOfGraphs() || mg.GetListOfGraphs()->IsEmpty()) continue;
        c.Clear();
        const std::string from = ci == 1 ? " from 0.3 to #theta" : " up to #theta";
        mg.SetTitle(("Cumulative percent of CE protons" + from + ", " + energyName.at(en) + ", " + setName.at(set) + Form(" (starts at N #geq %g protons)", minCumN) +
                     ";recon #theta (deg);cumulative 100 #times N_{CE} / N_{all protons} (%)").c_str());
        mg.SetMinimum(0); mg.SetMaximum(100); // a percent cannot pass 100
        mg.Draw("A");
        mg.GetXaxis()->SetLimits(lo, hi);
        leg.DrawClone();
        const std::string cumStem = "_cumulative_ce_percent" + cumTags[ci];
        SaveLogCopy(c, mg, dir + "2_plots/2_" + std::to_string(2 * ci + 2) + cumStem + "_log.pdf");
    }
}
