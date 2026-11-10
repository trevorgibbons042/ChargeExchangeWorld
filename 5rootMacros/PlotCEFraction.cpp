#include <TCanvas.h>
#include "Palette.h"
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
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <set>
#include <vector>

#include <TROOT.h>
#include "FitCore.h"
#include "FomCommon.h"

// Running percent of CE protons among reconstructed protons up to each theta, 100 * sum N_CE / sum N over the 1_csv/ke_*.csv bins written by
// PlotAngle.cpp, one curve per KE cut: ke_gt80/1_cumulative_ce_percent_log.pdf (log y axis only), beside the csv folder.
void PlotCEFraction(){
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    gStyle->SetOptStat(0);
    const auto &sets = Sets(), &kes = Kes();
    const int* colors = FomColors();
    const auto &setName = SetNames();
    auto cutNote = [](const std::string&){return std::string("no KE cut");};

    const double thetaMax = 15.0; // deg, fixed axis range for every plot (ReadThetaMax() ran to the largest RECONTHETA_ESTNEU)

    auto readCsv = [](const std::string& path){
        std::vector<std::array<double, 5>> rows; // theta_lo, theta_hi, N (recon protons), N_CE, N_particles (every fitted track; = N for a csv without that column)
        std::ifstream in(path);
        if (!in) Error("PlotCEFraction", "%s missing (PlotAngle.cpp writes it, and needs cuts/cut_APandT.root from CalcAPandT.cpp)", path.c_str());
        std::string line; std::getline(in, line);
        while (std::getline(in, line)){
            std::vector<double> x; std::stringstream ss(line); std::string f;
            while (std::getline(ss, f, ',')) x.push_back(std::stod(f));
            if (x.size() < 4 || x[1] - x[0] > 0.5) continue; // skip the total row
            rows.push_back({x[0], x[1], x[2], x[3], x.size() > 8 ? x[8] : x[2]});
        }
        return rows;
    };

    // numbers behind each plot go to ke_gt80/stats.csv, shared with PlotFomCumulative.cpp: this macro replaces only its own rows (plot names starting "1_")
    std::set<std::string> statsStarted;
    const std::string statsHeader = "plot,ke_cut,N,N_CE,pct_CE,fom";
    for (const auto& set : sets){
        const std::string file = "plots/5_fomCalc/" + set + "/ke_gt80/stats.csv";
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
        gSystem->mkdir((dir + "ke_gt80").c_str(), true);
        const std::string xTitle = set == "2_reconBlockTheta_estneu" ? "RECONTHETABL_ESTNEU" : "RECONTHETA_ESTNEU";
        // clears the canvas, titles the graph (title; x axis; yTitle), styles it in color and draws it over 0 to thetaMax; the caller saves it
        const auto drawGraph = [&](TGraph& g, const std::string& title, const std::string& yTitle, int color){
            c.Clear();
            g.SetTitle((title + ";" + xTitle + " (deg);" + yTitle).c_str());
            DrawPointsLine(g, color, thetaMax);
        };
        // Not cumulative: percent of the simulated neutrons whose CE proton lands in each theta bin, 100 * N_CE(bin) / N_neutrons (TotalNeutrons, settings.csv)
        {
            const double binW = 0.05, nNeutrons = ReadSettings().at("TotalNeutrons");
            const auto rows = readCsv(dir + "1_csv/ke_gt80.csv");
            const int nBins = (int)std::ceil(thetaMax / binW);
            std::vector<double> nCE(nBins, 0.0), nRecon(nBins, 0.0), nPart(nBins, 0.0);
            for (const auto& r : rows){
                const int b = (int)(0.5 * (r[0] + r[1]) / binW);
                if (b >= 0 && b < nBins){nCE[b] += r[3]; nRecon[b] += r[2]; nPart[b] += r[4];}
            }
            TGraph g;
            for (int b = 0; b < nBins; b++) g.SetPoint(g.GetN(), (b + 0.5) * binW, 100.0 * nCE[b] / nNeutrons);
            c.SetLogy(0);
            drawGraph(g, "100 #times N_{CE} / N_{neutrons} per #theta bin (0.05 deg), " + cutNote(set) + ", N_{CE} = CE recon protons, " + setName.at(set), "100 #times N_{CE} / N_{neutrons} (%)", CbAll()); // 3_: green
            c.SaveAs((dir + "ke_gt80/3_ce_percent_of_neutrons_per_bin.pdf").c_str());

            // The same, times A_y^2 at the bin center: t from theta (TFromAngle), A_y from t (APFromT, the model formula) with the run's APScale
            const double apScale = ReadSettings().at("APScale");
            TGraph gA2;
            for (int b = 0; b < nBins; b++){
                const double theta = (b + 0.5) * binW;
                const double ap = APFromT(kBeamMomentum, TFromAngle(kBeamMomentum, theta * TMath::Pi() / 180), apScale);
                gA2.SetPoint(gA2.GetN(), theta, ap * ap * 100.0 * nCE[b] / nNeutrons);
            }
            drawGraph(gA2, "100 #times N_{CE} / N_{neutrons} #times A_{y}^{2}(#theta) per #theta bin (0.05 deg), " + cutNote(set) + ", N_{CE} = CE recon protons, " + setName.at(set), "A_{y}^{2} #times 100 #times N_{CE} / N_{neutrons} (%)", CbAll()); // 4_: green
            c.SaveAs((dir + "ke_gt80/4_FOM_per_bin.pdf").c_str());

            // Running integral up to each bin: sum of N_CE / N_neutrons * A_y^2(theta) * delta theta (0.05 deg)
            TGraph gInt;
            double integral = 0;
            for (int b = 0; b < nBins; b++){
                const double theta = (b + 0.5) * binW;
                const double ap = APFromT(kBeamMomentum, TFromAngle(kBeamMomentum, theta * TMath::Pi() / 180), apScale);
                integral += nCE[b] / nNeutrons * ap * ap * binW;
                gInt.SetPoint(gInt.GetN(), (b + 1) * binW, integral);
            }
            drawGraph(gInt, "#Sigma N_{CE} / N_{neutrons} #times A_{y}^{2}(#theta) #times #Delta#theta up to #theta (#Delta#theta = 0.05 deg), " + cutNote(set) + ", N_{CE} = CE recon protons, " + setName.at(set), "#Sigma N_{CE} / N_{neutrons} #times A_{y}^{2} #times #Delta#theta (deg)", CbCE()); // CE protons are red
            TPaveText finalValue(0.81, 0.78, 0.99, 0.88, "NDC"); // outside the graph, in the right margin: the integral up to the last bin
            finalValue.SetFillColor(0); finalValue.SetBorderSize(0); finalValue.SetTextAlign(12); finalValue.SetTextSize(0.028);
            finalValue.AddText("final value");
            finalValue.AddText(Form("%.4g", integral));
            finalValue.Draw();
            c.SaveAs((dir + "ke_gt80/5_FOM_cumulative.pdf").c_str());

            // Per theta bin, the share of the reconstructed protons (1_, 2_) and of every reconstructed charged particle (6_, 7_) that are CE (100 * N_CE / N) and not CE (100 * (N - N_CE) / N)
            for (const bool part : {true}) for (const bool nonCE : {false, true}){ // all recon charged particles only, no all recon protons share
                const auto& nDen = part ? nPart : nRecon;
                const std::string den = part ? "N_{recon,particle}" : "N_{recon,proton}", num = nonCE ? "N_{non-CE,recon}" : "N_{CE,recon}";
                TGraph gShare;
                for (int b = 0; b < nBins; b++){
                    if (nDen[b] <= 0) continue; // an empty bin has no share
                    gShare.SetPoint(gShare.GetN(), (b + 0.5) * binW, 100.0 * (nonCE ? nDen[b] - nCE[b] : nCE[b]) / nDen[b]);
                }
                if (gShare.GetN() == 0) continue;
                // 1_, 2_ (all recon protons as the denominator) blue, 6_, 7_ (all recon charged particles) red
                drawGraph(gShare, std::string(nonCE ? "Contamination, " : "Purity, ") + (part ? "all recon charged particles: " : "all recon protons: ") + "100 #times " + num + " / " + den + " per #theta bin (0.05 deg), N_{CE,recon} = CE recon protons, " + den + " = all recon " + (part ? "charged particles (any pdg)" : "protons (pdg 2212)") + ", " + cutNote(set) + ", " + setName.at(set),
                          "100 #times " + num + " / " + den + " (%)", part ? CbCE() : CbProton());
                c.SaveAs((dir + "ke_gt80/" + (part ? (nonCE ? "7_contamination_particle_percent_per_bin.pdf" : "6_purity_particle_percent_per_bin.pdf") : (nonCE ? "2_contamination_proton_percent_per_bin.pdf" : "1_purity_proton_percent_per_bin.pdf"))).c_str());
            }

            // 8 to 15: the ratios of the cut folders per theta bin (the bin is the cut), CE = CE recon protons, prot = recon protons (pdg 2212), part = every recon charged particle; out / before are against the whole theta range
            {
                double totCE = 0, totProt = 0, totPart = 0;
                for (int b = 0; b < nBins; b++){totCE += nCE[b]; totProt += nRecon[b]; totPart += nPart[b];}
                (void)totProt;
                const double nonPartBefore = totPart - totCE;
                struct Ratio {const char* file, *label; int color; std::function<double(double, double, double)> f;}; // f(CE, prot, part) in the bin
                const int fileNum[5] = {8, 9, 11, 13, 15}; // the numbers the particle ratios kept when the all recon protons ones (10, 12, 14) went
                const Ratio ratios[5] = {
                    {"purity_particle", "N_{CE,recon,in bin} / N_{recon,particle,in bin}", CbAll(), [&](double ce, double, double pa){return pa > 0 ? ce / pa : NAN;}},
                    {"contamination_particle", "N_{nonCE,recon,particle,in bin} / N_{recon,particle,in bin}", CbAll(), [&](double ce, double, double pa){return pa > 0 ? (pa - ce) / pa : NAN;}},
                    {"background_in_over_out_particle", "N_{nonCE,recon,particle,in bin} / N_{nonCE,recon,particle,out of bin}", CbAll(), [&](double ce, double, double pa){const double out = nonPartBefore - (pa - ce); return out > 0 ? (pa - ce) / out : NAN;}},
                    {"background_efficiency_particle", "N_{nonCE,recon,particle,in bin} / N_{nonCE,recon,particle,all bins}", CbAll(), [&](double ce, double, double pa){return nonPartBefore > 0 ? (pa - ce) / nonPartBefore : NAN;}},
                    {"signal_significance_particle", "N_{CE,recon,in bin} / sqrt(N_{nonCE,recon,particle,in bin})", CbAll(), [&](double ce, double, double pa){return pa - ce > 0 ? ce / std::sqrt(pa - ce) : NAN;}}};
                for (int k = 0; k < 5; k++){
                    TGraph gR;
                    for (int b = 0; b < nBins; b++){const double v = ratios[k].f(nCE[b], nRecon[b], nPart[b]); if (std::isfinite(v)) gR.SetPoint(gR.GetN(), (b + 0.5) * binW, v);}
                    if (gR.GetN() < 2) continue;
                    std::string quantity = ratios[k].file; // the file name without underscores leads the title: "purity particle", "signal significance proton"
                    for (char& ch : quantity) if (ch == '_') ch = ' ';
                    drawGraph(gR, quantity + ": " + ratios[k].label + " per #theta bin (0.05 deg), N_{CE,recon} = CE recon protons, " + cutNote(set) + ", " + setName.at(set), ratios[k].label, ratios[k].color);
                    c.SaveAs((dir + "ke_gt80/" + std::to_string(fileNum[k]) + "_" + ratios[k].file + "_per_bin.pdf").c_str());
                }
            }
        }
    }
}
