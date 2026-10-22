#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "FitCore.h"

// Polarization from the phi modulation in thin theta slices. Per slice (kThetaStep deg wide, 0 to kThetaMax) the phi counts are fit with
// N0 * (1 + AP * P * cos(phi)), AP fixed to the mean AP of the CE protons in the slice and P free. Done for the CE protons and for all protons
// (their phi is fit with the same AP, so the unpolarized non-CE protons pull P down), for true and block position and every KE cut.
// Writes plots/6_fomCalc/<energy>/3_polarization/<theta set: recon or recon block>/<slice width>/<ap source>/<protons>/: csv/<ke>.csv, 1_P_vs_theta.pdf, 2_pct_error_vs_theta.pdf (vs PPValue
// in settings.csv) and 3_phi_fits_ke_all.pdf (one page per fitted slice, no KE cut). Reads cut_APandT.root and cut_angles_all.root.
static double kThetaStep = 0.1; // deg, set per pass over the slice widths
static const double kThetaMax = 5.0; // deg
static std::string gApName; // which AP the current pass uses, named in the fit box
static const int kPhiBins = 100;
static const size_t kMinN = 10; // protons needed in a slice for its fit

struct KECut {std::string dir, cut, label;};
struct Slice {double lo, hi, nPhi, ap, p, pErr;};

// Draws expr (nCols columns separated by ':') once under cut and copies the columns out. The same (tree, expr, cut) comes back from a cache:
// the slice width and AP-source loops ask for the same columns many times.
static void Gather(TTree* t, const std::string& expr, const std::string& cut, int nCols, std::vector<std::vector<double>>& cols){
    static std::map<std::string, std::vector<std::vector<double>>> cache;
    const std::string key = std::string(t->GetName()) + "|" + expr + "|" + cut;
    auto it = cache.find(key);
    if (it == cache.end()){
        t->SetEstimate(t->GetEntries() + 1);
        const Long64_t n = t->Draw(expr.c_str(), cut.c_str(), "goff");
        std::vector<std::vector<double>> drawn(nCols);
        for (int i = 0; i < nCols; i++) drawn[i].assign(t->GetVal(i), t->GetVal(i) + n);
        it = cache.emplace(key, std::move(drawn)).first;
    }
    cols = it->second;
}

static int SliceOf(double theta){return (theta >= 0 && theta < kThetaMax) ? (int)(theta / kThetaStep) : -1;}

// Fits one slice's phi values (deg). With draw, the histogram and fit are drawn in the current pad.
static bool FitSlice(const std::vector<double>& phis, Slice& s, bool draw){
    if (phis.size() < kMinN || std::fabs(s.ap) < 1e-12) return false;
    TH1D h("hPolPhi", "", kPhiBins, 0, 360);
    h.SetDirectory(nullptr);
    for (double phi : phis) h.Fill(phi);
    TF1 f("fPol", "[0]*(1+[1]*[2]*cos(x*TMath::DegToRad()))", 0, 360);
    f.SetParameters(h.GetEntries() / kPhiBins, 0.5, s.ap);
    f.SetParName(0, "N0"); f.SetParName(1, "P"); f.SetParName(2, "AP");
    f.FixParameter(2, s.ap);
    f.SetParLimits(1, -1, 1); // P is a polarization
    if (h.Fit(&f, draw ? "LQ" : "LQ0") != 0) return false; // drawn pages keep the fit in the histogram, like the other plots' fit boxes
    s.nPhi = h.GetEntries(); s.p = f.GetParameter(1); s.pErr = f.GetParError(1);
    if (draw){
        h.SetTitle(PlotTitle("3_phi_fits.pdf", Form("Phi, theta %.2f to %.2f deg, AP = %.4f", s.lo, s.hi, s.ap) + std::string(";Angle (degrees);Count")).c_str());
        h.GetXaxis()->SetNdivisions(505);
        h.SetStats(0);
        h.DrawCopy();
        // fit info box in the right margin, outside the graph
        TPaveText* box = new TPaveText(0.76, 0.58, 0.995, 0.88, "NDC");
        box->SetBorderSize(1); box->SetFillColor(0); box->SetTextAlign(12); box->SetTextSize(0.018);
        box->AddText("Fit function:");
        box->AddText("N_{0} (1 + AP #times P #times cos#phi)");
        box->AddText(gApName.c_str());
        box->AddText(Form("AP_{avg} = %.4f  (fixed)", s.ap));
        box->AddText(Form("P = %.3f #pm %.3f", s.p, s.pErr));
        box->AddText(Form("N_{0} = %.2f #pm %.2f", f.GetParameter(0), f.GetParError(0)));
        box->AddText(Form("entries = %.0f", s.nPhi));
        box->Draw();
    }
    return true;
}

void PlotPolarization(std::string energy = "fKEnergy"){ // energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the input trees and plot folder
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    const std::string energyDir = (energy == "fDEnergy") ? "2_KEDetector" : "1_KEtrue";
    const std::string keLabel = (energy == "fDEnergy") ? "ZDC KE" : "true KE";
    const double pRef = ReadSettings().at("PPValue");
    std::set<std::string> statsStarted; // stats.csv files already truncated in this run

    TFile* ceFile = TFile::Open("cut_APandT.root", "READ");
    TFile* allFile = TFile::Open("cut_angles_all.root", "READ");
    TTree* all = allFile ? (TTree*)allFile->Get(("Angles" + suffix).c_str()) : nullptr;
    if (!ceFile || !all){printf("PlotPolarization: cut_APandT.root or cut_angles_all.root missing\n"); return;}

    const double initialKE = 25000; // MeV
    std::vector<KECut> kes = {{"ke_all", "", "no KE cut"}};
    for (int pct : {50, 60, 70, 80, 90}) kes.push_back({"ke_gt" + std::to_string(pct), "finalKE > " + std::to_string((int)(initialKE * pct / 100.0)), "KE cut " + std::to_string(pct) + "%"});

    struct PosSet {std::string dir, name, ceTree, theta, phi, ap, allTheta, allPhi;};
    const std::vector<PosSet> sets = {
        {"1_reconTheta", "recon theta", "CalcValsRe", "thetaRe", "phiRe", "apRe", "outTheta", "outPhi"},
        {"2_reconBlockTheta", "recon block theta", "CalcValsReD", "thetaReD", "phiReD", "apReD", "outDTheta", "outDPhi"}};
    const std::vector<std::string> typeDirs = {"1_ce_protons", "2_all_protons"}, typeNames = {"CE protons", "all protons"};
    const int colors[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta, kOrange + 7};
    const std::string deg = "*180/TMath::Pi()";

    TCanvas pageC("cPolPage", "", 1500, 1200); // phi pages and graphs, same canvas and margins as the PlotAngle.cpp histograms
    pageC.SetRightMargin(0.28); pageC.SetLeftMargin(0.14); pageC.SetBottomMargin(0.13);

    struct ApSrc {std::string dir, name, tree, branch;}; // where each proton's AP comes from; theta and phi come from the theta set (recon or recon block), rows line up between the trees
    const std::vector<ApSrc> apSrcs = {
        {"1_apCE", "AP from CE", "CE" + suffix, "apCE"},
        {"2_apRecon", "AP from recon", "CalcValsRe" + suffix, "apRe"},
        {"3_apReconBlock", "AP from recon block", "CalcValsReD" + suffix, "apReD"}};

    for (const double width : {0.1, 0.05, 0.01}) for (const auto& apSrc : apSrcs){
    kThetaStep = width; gApName = apSrc.name;
    const int nSlices = (int)std::round(kThetaMax / kThetaStep);
    const std::string widthDir = Form("%gdeg_slices/", width);
    if (!ceFile->Get(apSrc.tree.c_str())){printf("PlotPolarization: tree %s missing\n", apSrc.tree.c_str()); continue;}
    for (const auto& set : sets){
        TTree* ce = (TTree*)ceFile->Get((set.ceTree + suffix).c_str()); // theta and phi of the protons, binned in this theta
        if (!ce){printf("PlotPolarization: tree %s missing\n", (set.ceTree + suffix).c_str()); continue;}
        const bool apInOwnTree = apSrc.tree == set.ceTree + suffix;
        if (!apInOwnTree && !ce->GetFriend(apSrc.tree.c_str())) ce->AddFriend(apSrc.tree.c_str());

        std::vector<std::vector<Slice>> fits[2]; // [protons type][ke cut] -> fitted slices
        for (int type = 0; type < 2; type++){
            const std::string base = "plots/6_fomCalc/" + energyDir + "/3_polarization/" + set.dir + "/" + widthDir + apSrc.dir + "/" + typeDirs[type] + "/";
            gSystem->mkdir((base + "csv").c_str(), true);

            fits[type].resize(kes.size());
            for (size_t k = 0; k < kes.size(); k++){
                const std::string cut = kes[k].cut.empty() ? "1" : kes[k].cut;
                gSystem->mkdir((base + kes[k].dir).c_str(), true);
                const std::string pagesPdf = base + kes[k].dir + "/3_phi_fits.pdf"; // one page per fitted slice
                bool pagesOpen = false;
                std::vector<std::vector<double>> ceCols, allCols;
                Gather(ce, set.theta + ":" + set.phi + ":" + (apInOwnTree ? "" : apSrc.tree + ".") + apSrc.branch, cut, 3, ceCols);
                Gather(all, set.allTheta + deg + ":" + set.allPhi + deg, cut, 2, allCols);

                std::vector<double> apSum(nSlices, 0.0), apN(nSlices, 0.0);
                std::vector<std::vector<double>> cePhi(nSlices), allPhi(nSlices);
                for (size_t i = 0; i < ceCols[0].size(); i++){
                    const int b = SliceOf(ceCols[0][i]);
                    if (b < 0) continue;
                    cePhi[b].push_back(ceCols[1][i]); apSum[b] += ceCols[2][i]; apN[b]++;
                }
                for (size_t i = 0; i < allCols[0].size(); i++){
                    const int b = SliceOf(allCols[0][i]);
                    if (b >= 0) allPhi[b].push_back(allCols[1][i]);
                }

                FILE* out = fopen((base + "csv/" + kes[k].dir + ".csv").c_str(), "w");
                if (out) fprintf(out, "theta_lo,theta_hi,N,AP_mean,P,P_err,pct_err\n");
                for (int b = 0; b < nSlices; b++){
                    if (apN[b] == 0) continue; // no AP without a CE proton in the slice
                    Slice s{b * kThetaStep, (b + 1) * kThetaStep, 0, apSum[b] / apN[b], 0, 0};
                    pageC.cd(); pageC.Clear();
                    if (!FitSlice(type == 0 ? cePhi[b] : allPhi[b], s, true)) continue;
                    if (!pagesOpen){pageC.Print((pagesPdf + "[").c_str()); pagesOpen = true;}
                    pageC.Print(pagesPdf.c_str());
                    fits[type][k].push_back(s);
                    if (out) fprintf(out, "%g,%g,%.0f,%.8g,%.8g,%.8g,%.8g\n", s.lo, s.hi, s.nPhi, s.ap, s.p, s.pErr, 100.0 * (s.p - pRef) / pRef);
                }
                if (out) fclose(out);
                if (pagesOpen) pageC.Print((pagesPdf + "]").c_str());

                // One fit over a whole theta range instead of slices: 0.5 to 2 deg and 0 to kThetaMax deg, AP = mean AP of the CE protons in the range.
                // Independent of the slice width, so only written during the first width pass: 3_polarization/total_fits/<range>/<ap>/<protons>/<ke>/, one folder per range.
                if (width == 0.1){
                    for (const auto& range : std::vector<std::pair<double, double>>{{0.5, 2.0}, {0.0, kThetaMax}}){
                        const std::string totDir = "plots/6_fomCalc/" + energyDir + "/3_polarization/" + set.dir + "/total_fits/" + (range.first > 0 ? "1_0p5to2deg" : "2_0to5deg") + "/" + apSrc.dir + "/" + typeDirs[type] + "/" + kes[k].dir + "/";
                        gSystem->mkdir(totDir.c_str(), true);
                        const std::string totHeader = "plot,theta_lo,theta_hi,N,AP_mean,P,P_err,pct_err";
                        double apSum2 = 0, apN2 = 0;
                        for (size_t i = 0; i < ceCols[0].size(); i++){
                            if (ceCols[0][i] < range.first || ceCols[0][i] >= range.second) continue;
                            apSum2 += ceCols[2][i]; apN2++;
                        }
                        if (apN2 == 0) continue;
                        std::vector<double> phis;
                        if (type == 0){
                            for (size_t i = 0; i < ceCols[0].size(); i++)
                                if (ceCols[0][i] >= range.first && ceCols[0][i] < range.second) phis.push_back(ceCols[1][i]);
                        } else {
                            for (size_t i = 0; i < allCols[0].size(); i++)
                                if (allCols[0][i] >= range.first && allCols[0][i] < range.second) phis.push_back(allCols[1][i]);
                        }
                        Slice t{range.first, range.second, 0, apSum2 / apN2, 0, 0};
                        pageC.cd(); pageC.Clear();
                        if (!FitSlice(phis, t, true)) continue;
                        const std::string totPdf = totDir + "1_phi_fit_total.pdf";
                        pageC.SaveAs(totPdf.c_str());
                        WriteStatsRow(totPdf, totHeader, Form("%s,%g,%g,%.0f,%.8g,%.8g,%.8g,%.8g", gSystem->BaseName(totPdf.c_str()), t.lo, t.hi, t.nPhi, t.ap, t.p, t.pErr, 100.0 * (t.p - pRef) / pRef), statsStarted);
                    }
                }
            }

            // P vs theta (pctPlot false) and its percent error vs PPValue (pctPlot true), one plot per KE cut: <base>/<ke>/1_P_vs_theta.pdf, 2_pct_error_vs_theta.pdf.
            // stepPlot: histogram-like steps, each fitted slice a horizontal line from its lower to its upper edge. Otherwise (_points files)
            // one marker per slice at its center, without the slice width.
            for (size_t k = 0; k < kes.size(); k++){
                if (fits[type][k].empty()) continue;
                // <base>/<ke>/stats.csv: one row per plot in the folder, from the fitted slices that plot shows (the phi pages hold the same slices)
                {
                    double nTot = 0, pSum = 0, pctSum = 0, pctMin = 1e30, pctMax = -1e30;
                    for (const auto& s : fits[type][k]){
                        const double pct = 100.0 * (s.p - pRef) / pRef;
                        nTot += s.nPhi; pSum += s.p; pctSum += pct; pctMin = std::min(pctMin, pct); pctMax = std::max(pctMax, pct);
                    }
                    const double nFit = fits[type][k].size();
                    const std::string row = Form("%.0f,%.0f,%.8g,%.8g,%.8g,%.8g", nFit, nTot, pSum / nFit, pctSum / nFit, pctMin, pctMax);
                    for (const char* f : {"3_phi_fits.pdf", "1_P_vs_theta.pdf", "1_P_vs_theta_points.pdf", "2_pct_error_vs_theta.pdf", "2_pct_error_vs_theta_points.pdf", "2_pct_error_vs_theta_log.pdf", "2_pct_error_vs_theta_points_log.pdf"})
                        WriteStatsRow(base + kes[k].dir + "/" + f, "plot,slices_fitted,N_in_slices,P_mean,pct_err_mean,pct_err_min,pct_err_max", std::string(f) + "," + row, statsStarted);
                }
                for (const bool stepPlot : {true, false}) for (const bool pctPlot : {false, true}) for (const bool logPlot : {false, true}){
                    if (logPlot && !pctPlot) continue; // the log plots show the absolute percent error (log can't show negative values)
                    const std::string file = std::string(pctPlot ? "2_pct_error_vs_theta" : "1_P_vs_theta") + (stepPlot ? "" : "_points") + (logPlot ? "_log" : "") + ".pdf";
                    std::vector<TGraph> steps; // one polyline per run of touching slices
                    TGraph points;
                    const double xPad = 0.02 * (fits[type][k].back().hi - fits[type][k].front().lo), xLo = fits[type][k].front().lo - xPad, xHi = fits[type][k].back().hi + xPad; // x axis fits the fitted slices
                    double prevHi = -1, yMin = logPlot ? 1e30 : (pctPlot ? 0.0 : pRef), yMax = logPlot ? 0.0 : yMin; // the reference line stays inside the plot
                    for (const auto& s : fits[type][k]){
                        const double pct = 100.0 * (s.p - pRef) / pRef;
                        const double y = pctPlot ? (logPlot ? std::max(std::fabs(pct), 1e-6) : pct) : s.p;
                        if (stepPlot){
                            if (steps.empty() || s.lo > prevHi + 1e-9) steps.emplace_back();
                            steps.back().AddPoint(s.lo, y); steps.back().AddPoint(s.hi, y);
                        } else {
                            points.AddPoint(0.5 * (s.lo + s.hi), y);
                        }
                        prevHi = s.hi;
                        yMin = std::min(yMin, y); yMax = std::max(yMax, y);
                    }
                    pageC.cd(); pageC.Clear();
                    const double pad = 0.1 * (yMax - yMin > 0 ? yMax - yMin : 1.0);
                    const std::string title = std::string(pctPlot ? (logPlot ? "Absolute percent error of fitted polarization vs PPValue, log scale" : "Percent error of fitted polarization vs PPValue") : "Fitted polarization P") + " vs #theta, " +
                                typeNames[type] + ", " + apSrc.name + ", " + set.name + ", " + keLabel + ", " + kes[k].label + (stepPlot ? Form(", %g deg #theta slices", kThetaStep) : "") +
                                ";recon #theta (deg);" + (pctPlot ? (logPlot ? "|100 #times (P_{fit} - P_{true}) / P_{true}| (%)" : "100 #times (P_{fit} - P_{true}) / P_{true} (%)") : "P_{fit}");
                    TH1F frame("framePol", PlotTitle(file, title).c_str(), 1, xLo, xHi); // axes and title only, the data is drawn on top
                    frame.SetStats(0);
                    frame.SetMinimum(logPlot ? yMin / 2.0 : yMin - pad); frame.SetMaximum(logPlot ? yMax * 2.0 : yMax + pad);
                    frame.GetXaxis()->SetNdivisions(505);
                    frame.SetBinContent(1, logPlot ? yMin * 1e-3 : yMin - 100.0 * (yMax - yMin + 1.0)); // far below the range so its line stays off the plot ("AXIS" would not draw the title)
                    pageC.SetLogy(logPlot);
                    frame.Draw("HIST");
                    for (auto& g : steps){g.SetLineWidth(1); g.DrawClone("L");}
                    if (!stepPlot){points.SetMarkerStyle(20); points.SetMarkerSize(1.2); points.DrawClone("P");}
                    if (!logPlot){ // the reference (true P, 0% error) is not on a log axis
                        TLine ref(xLo, pctPlot ? 0.0 : pRef, xHi, pctPlot ? 0.0 : pRef);
                        ref.SetLineStyle(2);
                        ref.DrawClone();
                    }
                    pageC.SaveAs((base + kes[k].dir + "/" + file).c_str());
                    pageC.SetLogy(0);
                }
            }
        }
    }
    }
}
