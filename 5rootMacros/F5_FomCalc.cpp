#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "FitCore.h"
#include "CutFiles.h"
#include "TreeColumns.h"

// plots/5_fomCalc: the figure of merit and the CE purity per theta bin, from the event recon (AnglesRaw of cuts/cut_angles_all.root, ReconAngleAll.cpp), no cuts.
// Per bin of 0.05 deg of the est neutron theta (true position: thetaRelN, block position: thetaRelND):
//   N_CE   = fitted protons of the events with a CE hit (evClass 0 and 1, the signal of F9_PlotReconNew.cpp)
//   N_part = every fitted proton (evClass 0 to 2), so N_part - N_CE is the background
//   FOM per bin = N_CE / N_neutrons * A_y^2(theta of the bin center), N_neutrons = TotalNeutrons of settings.csv, A_y from AEffOfTheta (FitCore.h) with the APScale of the run.
// Numbers behind the plots: <set>/1_noCuts/bins.csv.
void F5_FomCalc(){
    FastJit();
    ROOT::EnableImplicitMT(NThreads());
    gStyle->SetOptStat(0);

    TTree* t = OpenCut("cuts/cut_angles_all.root", kAngleTree); // with the Events tree as a friend
    if (!t){Error("F5_FomCalc", "no AnglesRaw in cuts/cut_angles_all.root (run ReconAngleAll.cpp)"); return;}

    const auto settings = ReadSettings();
    const double nNeutrons = settings.at("TotalNeutrons"), apScale = settings.at("APScale");
    const double binW = 0.05, thetaMax = 15.0; // deg
    const int nBins = (int)std::lround(thetaMax / binW);

    // theta of the est neutron frame in degrees for both positions and the event class: one pass over the tree on all threads, each slot counts into its own bins (no column vectors of every row in memory)
    ROOT::RDataFrame df0(*t);
    ROOT::RDF::RNode df = df0.Define("thTrue", "thetaRelN*180/TMath::Pi()").Define("thBlock", "thetaRelND*180/TMath::Pi()").Define("cls", "(int)evClass");
    struct Counts {std::vector<double> ce, part;};
    std::vector<std::array<Counts, 2>> slotCounts(df0.GetNSlots()); // [slot][set]: set 0 true position, set 1 block position
    for (auto& sc : slotCounts) for (auto& cn : sc){cn.ce.assign(nBins, 0.0); cn.part.assign(nBins, 0.0);}
    auto nRowsRes = df0.Count(); // booked with the ForeachSlot below: the same pass counts the rows
    df.ForeachSlot([&](unsigned slot, double thTrue, double thBlock, int cls){
        if (cls < 0 || cls > 2) return;
        const double th2[2] = {thTrue, thBlock};
        for (int k = 0; k < 2; k++){
            if (!std::isfinite(th2[k]) || th2[k] < 0) continue;
            const int b = (int)(th2[k] / binW);
            if (b >= nBins) continue;
            slotCounts[slot][k].part[b]++;
            if (cls <= 1) slotCounts[slot][k].ce[b]++;
        }
    }, {"thTrue", "thBlock", "cls"});
    const Long64_t nRows = (Long64_t)*nRowsRes;

    struct Set {const char* dir; const char* axis; const char* name; int col;};
    const Set sets[2] = {{"1_reconTheta_estneu", "RECONTHETA_ESTNEU", "true position", 0}, {"2_reconBlockTheta_estneu", "RECONTHETABL_ESTNEU", "block position", 1}};

    TCanvas c("c", "", 1200, 900);
    c.SetLeftMargin(0.14); c.SetRightMargin(0.2);
    int made = 0;
    for (const Set& set : sets){
        // counts per bin, the slots added (exact: whole numbers)
        std::vector<double> nCE(nBins, 0.0), nPart(nBins, 0.0);
        for (const auto& sc : slotCounts) for (int b = 0; b < nBins; b++){nCE[b] += sc[set.col].ce[b]; nPart[b] += sc[set.col].part[b];}
        double totCE = 0, totPart = 0;
        for (int b = 0; b < nBins; b++){totCE += nCE[b]; totPart += nPart[b];}
        const double nonCEAll = totPart - totCE;

        const std::string dir = std::string("plots/5_fomCalc/") + set.dir + "/1_noCuts";
        gSystem->mkdir(dir.c_str(), true);
        const std::string xTitle = set.axis;
        const std::string note = std::string("no cuts, N_{CE} = CE recon protons, ") + set.name;
        auto ay = [&](int b){return AEffOfTheta((b + 0.5) * binW, apScale);};

        { // the table behind the plots
            FILE* out = fopen((dir + "/bins.csv").c_str(), "w");
            if (out){
                fprintf(out, "theta_lo,theta_hi,N_CE,N_nonCE,N_particles,A_y,N_CE_over_neutrons,FOM\n");
                for (int b = 0; b < nBins; b++)
                    fprintf(out, "%g,%g,%.0f,%.0f,%.0f,%.8g,%.8g,%.8g\n", b * binW, (b + 1) * binW, nCE[b], nPart[b] - nCE[b], nPart[b], ay(b), nCE[b] / nNeutrons, nCE[b] / nNeutrons * ay(b) * ay(b));
                fclose(out);
            }
        }

        // round markers (size 0.6) joined by a width 1 line in one color, x from 0 to thetaMax
        auto save = [&](TGraph& g, const std::string& file, const std::string& title, const std::string& yTitle, int color, double finalValue = NAN){
            if (g.GetN() < 2) return;
            const std::string path = dir + "/" + file;
            c.Clear();
            g.SetTitle(PlotTitle(path, title + ";" + xTitle + " (deg);" + yTitle).c_str());
            g.SetMarkerStyle(20); g.SetMarkerSize(0.6); g.SetLineWidth(1); g.SetLineColor(color); g.SetMarkerColor(color);
            g.Draw("APL");
            g.GetXaxis()->SetLimits(0, thetaMax);
            TPaveText pave(0.81, 0.78, 0.99, 0.88, "NDC"); // right margin: the integral up to the last bin
            if (std::isfinite(finalValue)){
                pave.SetFillColor(0); pave.SetBorderSize(0); pave.SetTextAlign(12); pave.SetTextSize(0.028);
                pave.AddText("final value");
                pave.AddText(Form("%.4g", finalValue));
                pave.Draw();
            }
            c.SaveAs(path.c_str());
            made++;
        };

        TGraph gPct, gFom, gInt;
        double integral = 0;
        for (int b = 0; b < nBins; b++){
            const double th = (b + 0.5) * binW, share = 100.0 * nCE[b] / nNeutrons;
            gPct.SetPoint(gPct.GetN(), th, share);
            gFom.SetPoint(gFom.GetN(), th, ay(b) * ay(b) * share);
            integral += nCE[b] / nNeutrons * ay(b) * ay(b) * binW;
            gInt.SetPoint(gInt.GetN(), (b + 1) * binW, integral);
        }
        save(gPct, "3_ce_percent_of_neutrons_per_bin.pdf", "100 #times N_{CE} / N_{neutrons} per #theta bin (0.05 deg), " + note, "100 #times N_{CE} / N_{neutrons} (%)", CbAll());
        save(gFom, "4_FOM_per_bin.pdf", "100 #times N_{CE} / N_{neutrons} #times A_{y}^{2}(#theta) per #theta bin (0.05 deg), " + note, "A_{y}^{2} #times 100 #times N_{CE} / N_{neutrons} (%)", CbAll());
        // running integral up to each bin: sum of N_CE / N_neutrons * A_y^2 * delta theta
        save(gInt, "5_FOM_cumulative.pdf", "#Sigma N_{CE} / N_{neutrons} #times A_{y}^{2}(#theta) #times #Delta#theta up to #theta (#Delta#theta = 0.05 deg), " + note, "#Sigma N_{CE} / N_{neutrons} #times A_{y}^{2} #times #Delta#theta", CbAll(), integral);

        // purity and contamination of the recon particles, in percent
        for (const bool contamination : {false, true}){
            TGraph g;
            for (int b = 0; b < nBins; b++){
                if (nPart[b] <= 0) continue;
                g.SetPoint(g.GetN(), (b + 0.5) * binW, 100.0 * (contamination ? nPart[b] - nCE[b] : nCE[b]) / nPart[b]);
            }
            const std::string num = contamination ? "N_{non-CE,recon}" : "N_{CE,recon}";
            save(g, contamination ? "7_contamination_particle_percent_per_bin.pdf" : "6_purity_particle_percent_per_bin.pdf",
                 std::string(contamination ? "Contamination, " : "Purity, ") + "all recon particles: 100 #times " + num + " / N_{recon,particle} per #theta bin (0.05 deg), " + note,
                 "100 #times " + num + " / N_{recon,particle} (%)", CbCE());
        }

        // 11, 13 and 15: ratios per theta bin (the bin is the cut), against the whole theta range
        struct Ratio {const char* file; int number; const char* label; double (*f)(double ce, double part, double nonCEAll);};
        const Ratio ratios[3] = {
            {"background_in_over_out_particle", 11, "N_{nonCE,recon,particle,in bin} / N_{nonCE,recon,particle,out of bin}", [](double ce, double pa, double all){const double out = all - (pa - ce); return out > 0 ? (pa - ce) / out : NAN;}},
            {"background_efficiency_particle", 13, "N_{nonCE,recon,particle,in bin} / N_{nonCE,recon,particle,all bins}", [](double ce, double pa, double all){return all > 0 ? (pa - ce) / all : NAN;}},
            {"signal_significance_particle", 15, "N_{CE,recon,in bin} / sqrt(N_{nonCE,recon,particle,in bin})", [](double ce, double pa, double){return pa - ce > 0 ? ce / std::sqrt(pa - ce) : NAN;}}};
        for (const Ratio& r : ratios){
            TGraph g;
            for (int b = 0; b < nBins; b++){const double v = r.f(nCE[b], nPart[b], nonCEAll); if (std::isfinite(v)) g.SetPoint(g.GetN(), (b + 0.5) * binW, v);}
            std::string quantity = r.file;
            for (char& ch : quantity) if (ch == '_') ch = ' ';
            save(g, std::to_string(r.number) + "_" + r.file + "_per_bin.pdf", quantity + ": " + r.label + " per #theta bin (0.05 deg), " + note, r.label, CbAll());
        }
    }
    printf(">>> F5_FomCalc: %d plots in plots/5_fomCalc (%lld recon rows)\n", made, nRows);
}
