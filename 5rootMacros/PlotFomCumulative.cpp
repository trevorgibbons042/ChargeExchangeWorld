#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "FitCore.h"
#include "FomCommon.h"

// Running sum over theta of each fom.csv column (eff, eff_x_AP2), one curve per KE cut,
// written next to the csv files as cumulative_<column>.pdf. Reads plots/6_fomCalc written by PlotAngle.cpp.
static int colIndex(const std::string& header, const std::string& name){
    std::stringstream ss(header); std::string f; int i = 0;
    while (std::getline(ss, f, ',')){ if (f == name) return i; i++; }
    return 0;
}

void PlotFomCumulative(){
    gStyle->SetOptStat(0);
    const auto &energies = kEnergies, &sets = kSets, &kes = kKes;
    const std::vector<std::string> colNames = {"eff", "eff_x_AP2", "eff_x_AP2_x_CE"}; // fom.csv header names; the last is eff_x_AP2 times the running N_CE / N (not a csv column)
    // plot labels: csv column -> axis quantity, directory name -> readable name
    const std::map<std::string, std::string> quantity = {
        {"eff", "#varepsilon = N_{protons} / N_{neutrons}"}, {"eff_x_AP2", "#varepsilon A_{y}^{2}"},
        {"eff_x_AP2_x_CE", "#varepsilon A_{y}^{2} #times N_{CE} / N_{all protons}"}};
    const auto &energyName = kEnergyName, &setName = kSetName;
    const auto keName = KeName;
    const auto &colors = kColors;

    const double thetaMax = ReadThetaMax();

    std::set<std::string> statsStarted;
    TCanvas c("c", "", 1200, 900);
    c.SetLeftMargin(0.14); c.SetRightMargin(0.2); // room for the legend outside the graph
    for (const auto& en : energies) for (const auto& set : sets){
        const std::string dir = "plots/6_fomCalc/" + en + "/" + set + "/";
        gSystem->mkdir((dir + "2_plots").c_str(), true);
        // window: the sum only runs over the integrating region (bin centers from windowLo to windowHi deg), for the eff_x_AP2 plots only
        const double windowLo = 0.3, windowHi = 2.0;
        for (const bool window : {false, true}) for (const auto& colName : colNames){
            if (window && colName == "eff") continue;
            const auto col = std::make_pair(0, colName);
            TMultiGraph mg;
            // eff_x_AP2 plots (the FOM) list each curve's final value, so they need a wider legend outside the graph
            const bool weightCE = colName == "eff_x_AP2_x_CE"; // each point: running eff_x_AP2 times running N_CE / N
            const bool showFom = colName != "eff";
            // file numbers: 1_ and 2_ belong to PlotCEFraction.cpp (CE fraction per bin, cumulative CE percent). Here one group per column:
            // 3_ eff, 4_ eff_x_AP2, 5_ eff_x_AP2_x_CE; inside a group 1_ full range, 2_ its log, 3_ window, 4_ its log
            const size_t group = std::find(colNames.begin(), colNames.end(), colName) - colNames.begin() + 3;
            const std::string numOf = std::to_string(group) + "_";
            const std::string stem = "_cumulative_" + colName + (window ? "_0p3to2" : "");
            const std::string plotFile = dir + "2_plots/" + numOf + std::to_string(window ? 3 : 1) + stem + ".pdf";
            const std::string logFile = dir + "2_plots/" + numOf + std::to_string(window ? 4 : 2) + stem + "_log.pdf";
            c.SetRightMargin(showFom ? 0.3 : 0.2);
            TLegend leg(showFom ? 0.72 : 0.81, 0.74, 0.99, 0.88);
            TPaveText counts(0.70, 0.58, 0.99, 0.72, "NDC"); // protons summed into each curve, under the legend
            counts.SetBorderSize(0); counts.SetFillStyle(0); counts.SetTextAlign(12); counts.SetTextSize(0.02);
            leg.SetTextSize(0.022);
            leg.SetBorderSize(0);
            // the numbers shown on the eff_x_AP2 plots, one row per KE cut: 2_plots/summary_full_range.csv or summary_0p3to2.csv (_x_CE suffix for the CE-weighted plots)
            std::ofstream summary;
            if (showFom){
                summary.open(dir + "2_plots/summary_" + (window ? "0p3to2" : "full_range") + (weightCE ? "_x_CE" : "") + ".csv");
                summary << "ke_cut,N,N_CE,pct_CE,FOM\n";
            }
            for (size_t k = 0; k < kes.size(); k++){
                std::ifstream in(dir + "1_csv/" + kes[k] + ".csv");
                if (!in) continue;
                TGraph* g = new TGraph();
                std::string line; std::getline(in, line);
                const int iCol = colIndex(line, weightCE ? "eff_x_AP2" : colName), iN = colIndex(line, "N"), iCE = colIndex(line, "N_CE");
                double run = 0, nAll = 0, nCE = 0, y = 0; // nAll, nCE: protons in the summed bins; y: the plotted value (run, times nCE / nAll for the CE plot)
                while (std::getline(in, line)){
                    std::vector<double> v; std::stringstream ss(line); std::string f;
                    while (std::getline(ss, f, ',')) v.push_back(std::stod(f));
                    if (v.size() <= (size_t)std::max({iCol, iN, iCE}) || v[1] - v[0] > 0.5) continue; // skip the total row
                    if (window && (0.5 * (v[0] + v[1]) < windowLo || 0.5 * (v[0] + v[1]) >= windowHi)) continue;
                    run += v[iCol]; nAll += v[iN]; nCE += v[iCE];
                    y = weightCE ? (nAll > 0 ? run * nCE / nAll : 0.0) : run;
                    g->SetPoint(g->GetN(), v[1], y);
                }
                g->SetLineColor(colors[k]); g->SetLineWidth(1);
                mg.Add(g, "L");
                leg.AddEntry(g, (keName(kes[k]) + (showFom ? Form(" - FOM %.3e", y) : "")).c_str(), "l");
                // 2_plots/stats.csv: the final point of each curve (PlotCEFraction.cpp adds its own rows to the same file); "fom" is that value
                WriteStatsRow(logFile, "plot,ke_cut,N,N_CE,pct_CE,fom", Form("%s,%s,%.0f,%.0f,%.6g,%.6g", gSystem->BaseName(logFile.c_str()), kes[k].c_str(), nAll, nCE, nAll > 0 ? 100.0 * nCE / nAll : 0.0, y), statsStarted);
                if (showFom) summary << kes[k] << "," << nAll << "," << nCE << "," << (nAll > 0 ? 100.0 * nCE / nAll : 0.0) << "," << y << "\n";
                if (showFom) counts.AddText(Form("%s: N = %.0f (CE %.0f, %.1f%%)", keName(kes[k]).c_str(), nAll, nCE, nAll > 0 ? 100.0 * nCE / nAll : 0.0));
            }
            if (!mg.GetListOfGraphs() || mg.GetListOfGraphs()->IsEmpty()) continue;
            c.Clear();
            const std::string q = quantity.at(col.second);
            const std::string from = window ? " from 0.3 to #theta" : " up to #theta";
            mg.SetTitle(("Cumulative " + q + from + ", " + energyName.at(en) + ", " + setName.at(set) + ";" +
                         "recon #theta (deg)" + ";" +
                         "cumulative " + q + from).c_str());
            mg.Draw("A");
            mg.GetXaxis()->SetLimits(window ? windowLo : 0, window ? windowHi : thetaMax);
            leg.DrawClone();
            if (showFom) counts.DrawClone();
            SaveLogCopy(c, mg, logFile); // log y axis only, plotFile (the linear plot) is no longer saved
        }
    }
}
