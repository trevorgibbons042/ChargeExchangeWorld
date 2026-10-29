#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
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

// Running percent of CE protons among reconstructed protons up to each theta, 100 * sum N_CE / sum N over the 1_csv/ke_*.csv bins written by
// PlotAngle.cpp, one curve per KE cut: ke_all/1_cumulative_ce_percent_log.pdf (log y axis only), beside the csv folder.
void PlotCEFraction(){
    gStyle->SetOptStat(0);
    const auto &sets = kSets, &kes = kKes;
    const auto &colors = kColors;
    const auto &setName = kSetName;
    const auto keName = KeName;

    const double thetaMax = 15.0; // deg, fixed axis range for every plot (ReadThetaMax() ran to the largest recon theta)

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

    // numbers behind each plot go to ke_all/stats.csv, shared with PlotFomCumulative.cpp: this macro replaces only its own rows (plot names starting "1_")
    std::set<std::string> statsStarted;
    const std::string statsHeader = "plot,ke_cut,N,N_CE,pct_CE,fom";
    for (const auto& set : sets){
        const std::string file = "plots/5_fomCalc/" + set + "/ke_all/stats.csv";
        std::ifstream in(file);
        std::string line, kept;
        while (std::getline(in, line)) if (!line.empty() && line[0] >= '2' && line != statsHeader) kept += line + "\n";
        in.close();
        if (kept.empty()) continue;
        gSystem->mkdir(gSystem->GetDirName(file.c_str()), true);
        std::ofstream out(file);
        out << statsHeader << "\n" << kept;
        statsStarted.insert(file); statsStarted.insert(file + "#header");
    }

    TCanvas c("c", "", 1200, 900);
    c.SetLeftMargin(0.14); c.SetRightMargin(0.2); // room for the legend outside the graph
    for (const auto& set : sets){
        const std::string dir = "plots/5_fomCalc/" + set + "/";
        gSystem->mkdir((dir + "ke_all").c_str(), true);
        // Not cumulative: percent of the simulated neutrons whose CE proton lands in each theta bin, 100 * N_CE(bin) / N_neutrons (TotalNeutrons, settings.csv)
        {
            const double binW = 0.05, nNeutrons = ReadSettings().at("TotalNeutrons");
            const auto rows = readCsv(dir + "1_csv/ke_all.csv");
            const int nBins = (int)std::ceil(thetaMax / binW);
            std::vector<double> nCE(nBins, 0.0), nRecon(nBins, 0.0);
            for (const auto& r : rows){
                const int b = (int)(0.5 * (r[0] + r[1]) / binW);
                if (b >= 0 && b < nBins){nCE[b] += r[3]; nRecon[b] += r[2];}
            }
            TGraph g;
            for (int b = 0; b < nBins; b++) g.SetPoint(g.GetN(), (b + 0.5) * binW, 100.0 * nCE[b] / nNeutrons);
            c.Clear(); c.SetLogy(0);
            g.SetTitle(("100 #times N_{CE} / N_{neutrons} per #theta bin (0.05 deg), detector KE, N_{CE} = CE recon protons, " + setName.at(set) + ";" +
                        std::string(set == "2_reconBlockTheta" ? "recon block #theta" : "recon #theta") + " (deg);100 #times N_{CE} / N_{neutrons} (%)").c_str());
            g.SetMarkerStyle(20); g.SetMarkerSize(0.6); g.SetLineWidth(1);
            g.Draw("APL");
            g.GetXaxis()->SetLimits(0, thetaMax);
            c.SaveAs((dir + "ke_all/3_ce_percent_of_neutrons_per_bin.pdf").c_str());

            // The same, times A_y^2 at the bin center: t from theta (TFromAngle), A_y from t (APFromT, the model formula) with the run's APScale
            const double apScale = ReadSettings().at("APScale");
            TGraph gA2;
            for (int b = 0; b < nBins; b++){
                const double theta = (b + 0.5) * binW;
                const double ap = APFromT(kBeamMomentum, TFromAngle(kBeamMomentum, theta * TMath::Pi() / 180), apScale);
                gA2.SetPoint(gA2.GetN(), theta, ap * ap * 100.0 * nCE[b] / nNeutrons);
            }
            c.Clear();
            gA2.SetTitle(("100 #times N_{CE} / N_{neutrons} #times A_{y}^{2}(#theta) per #theta bin (0.05 deg), detector KE, N_{CE} = CE recon protons, " + setName.at(set) + ";" +
                          std::string(set == "2_reconBlockTheta" ? "recon block #theta" : "recon #theta") + " (deg);A_{y}^{2} #times 100 #times N_{CE} / N_{neutrons} (%)").c_str());
            gA2.SetMarkerStyle(20); gA2.SetMarkerSize(0.6); gA2.SetLineWidth(1);
            gA2.Draw("APL");
            gA2.GetXaxis()->SetLimits(0, thetaMax);
            c.SaveAs((dir + "ke_all/4_ce_percent_x_A2_per_bin.pdf").c_str());

            // Running integral up to each bin: sum of N_CE / N_neutrons * A_y^2(theta) * delta theta (0.05 deg)
            TGraph gInt;
            double integral = 0;
            for (int b = 0; b < nBins; b++){
                const double theta = (b + 0.5) * binW;
                const double ap = APFromT(kBeamMomentum, TFromAngle(kBeamMomentum, theta * TMath::Pi() / 180), apScale);
                integral += nCE[b] / nNeutrons * ap * ap * binW;
                gInt.SetPoint(gInt.GetN(), (b + 1) * binW, integral);
            }
            c.Clear();
            gInt.SetTitle(("#Sigma N_{CE} / N_{neutrons} #times A_{y}^{2}(#theta) #times #Delta#theta up to #theta (#Delta#theta = 0.05 deg), detector KE, N_{CE} = CE recon protons, " + setName.at(set) + ";" +
                           std::string(set == "2_reconBlockTheta" ? "recon block #theta" : "recon #theta") + " (deg);#Sigma N_{CE} / N_{neutrons} #times A_{y}^{2} #times #Delta#theta (deg)").c_str());
            gInt.SetMarkerStyle(20); gInt.SetMarkerSize(0.6); gInt.SetLineWidth(1);
            gInt.Draw("APL");
            gInt.GetXaxis()->SetLimits(0, thetaMax);
            TPaveText finalValue(0.81, 0.78, 0.99, 0.88, "NDC"); // outside the graph, in the right margin: the integral up to the last bin
            finalValue.SetFillColor(0); finalValue.SetBorderSize(0); finalValue.SetTextAlign(12); finalValue.SetTextSize(0.028);
            finalValue.AddText("final value");
            finalValue.AddText(Form("%.4g", integral));
            finalValue.Draw();
            c.SaveAs((dir + "ke_all/5_integral_ce_x_A2_dtheta.pdf").c_str());

            // Per theta bin, the share of the reconstructed protons that are CE (8_) and not CE (9_): 100 * N_CE / N_recon and 100 * (N_recon - N_CE) / N_recon
            for (const bool nonCE : {false, true}){
                TGraph gShare;
                for (int b = 0; b < nBins; b++){
                    if (nRecon[b] <= 0) continue; // an empty bin has no share
                    gShare.SetPoint(gShare.GetN(), (b + 0.5) * binW, 100.0 * (nonCE ? nRecon[b] - nCE[b] : nCE[b]) / nRecon[b]);
                }
                if (gShare.GetN() == 0) continue;
                c.Clear();
                gShare.SetTitle(std::string(nonCE ? "100 #times N_{non-CE} / N_{recon}" : "100 #times N_{CE} / N_{recon}").append(" per #theta bin (0.05 deg), N_{CE} = CE recon protons, N_{recon} = all recon protons, detector KE, ").append(setName.at(set)).append(";")
                                .append(set == "2_reconBlockTheta" ? "recon block #theta" : "recon #theta").append(nonCE ? " (deg);100 #times N_{non-CE} / N_{recon protons} (%)" : " (deg);100 #times N_{CE} / N_{recon protons} (%)").c_str());
                gShare.SetMarkerStyle(20); gShare.SetMarkerSize(0.6); gShare.SetLineWidth(1);
                gShare.Draw("APL");
                gShare.GetXaxis()->SetLimits(0, thetaMax);
                c.SaveAs((dir + (nonCE ? "ke_all/2_nonce_percent_of_recon_per_bin.pdf" : "ke_all/1_ce_percent_of_recon_per_bin.pdf")).c_str());
            }
        }
    }
}
