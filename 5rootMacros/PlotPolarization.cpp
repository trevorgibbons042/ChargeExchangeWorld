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
#include <TPad.h>
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

#include <TROOT.h>
#include "FitCore.h"
#include "CutFiles.h"
#include "TreeColumns.h"

// Polarization from the phi modulation in thin theta slices. Per slice (kThetaStep deg wide, 0 to kThetaMax) the phi counts are fit with
// N0 * (1 + AP * P * cos(phi)), AP fixed to the mean AP of the CE protons in the slice and P free. Done for the CE protons and for all protons
// (their phi is fit with the same AP, so the unpolarized non-CE protons pull P down), for true and block position and every KE cut.
// Writes plots/5_fomCalc/<theta set: recon or recon block>/<ke>/: slices_<tag>.csv, 1_P_vs_theta_<tag>.pdf, 2_pct_error_vs_theta_<tag>.pdf (vs PPValue
// in settings.csv), tag = <CE or recon> in <ke>/1_apCE/ and <ke>/2_apRecon/ (recon theta), apReconBlock_<CE or recon> directly in <ke>/ (recon block theta). Recon theta uses AP from CE and from recon; recon block theta only AP from recon block, with the block residual cut.
// The fit of every slice is in <ke>/slices_<tag>.csv, not plotted. Reads cut_APandT.root and cut_angles_all.root.
static double kThetaStep = 0.1; // deg, set per pass over the slice widths
static const double kThetaMax = 5.0; // deg
static std::string gApName; // which AP the current pass uses, named in the fit box
static const int kPhiBins = 100;
static const size_t kMinN = 10; // protons needed in a slice for its fit

struct Slice {double lo, hi, nPhi, ap, p, pErr, nCE = 0, nCEAll = -1, nReconAll = -1;}; // nCEAll, nReconAll: CE and recon protons of the theta set before any cut, for the CE lost and the recon share // nCE: CE protons in the slice (= nPhi for the CE fits)

// Reads expr (nCols columns separated by ':') once, multithreaded (Columns, TreeColumns.h), under the KE cut keDir (block: and the block residual cut, CutFiles.h) and copies the columns
// out. The same (tree, expr, cut) comes back from a cache: the slice width and AP-source loops ask for the same columns many times.
// listName: the cut list of a view of t (Angles + suffix of AnglesRaw, CutFiles.h); empty for a tree whose list is named like it.
static void Gather(TTree* t, const std::string& expr, const std::string& keDir, bool block, int nCols, std::vector<std::vector<double>>& cols, const std::string& listName = ""){
    static std::map<std::string, std::vector<std::vector<double>>> cache;
    const std::string key = (listName.empty() ? std::string(t->GetName()) : listName) + "|" + expr + "|" + keDir + (block ? "|block" : "");
    auto it = cache.find(key);
    if (it == cache.end()){
        ApplyCut(t, keDir, block, listName);
        std::vector<std::vector<double>> drawn = Columns(t, SplitDraw(expr).exprs);
        t->SetEntryList(nullptr);
        drawn.resize(nCols);
        it = cache.emplace(key, std::move(drawn)).first;
    }
    cols = it->second;
}

static int SliceOf(double theta){return (theta >= 0 && theta < kThetaMax) ? (int)(theta / kThetaStep) : -1;}

// Fills h (kPhiBins over 0 to 360) with phis (deg) and fits it with f = N0 (1 + AP P cos(phi)), AP fixed to ap. False when there are too few
// protons or the fit fails. opt "LQ" keeps the fit in h for drawing, "LQ0" does not.
static bool FitPhi(const std::vector<double>& phis, double ap, TH1D& h, TF1& f, const char* opt){
    if (phis.size() < kMinN || std::fabs(ap) < 1e-12) return false;
    h.SetDirectory(nullptr);
    for (double phi : phis) h.Fill(phi);
    f.SetParameters(h.GetEntries() / kPhiBins, 0.5, ap);
    f.SetParName(0, "N0"); f.SetParName(1, "P"); f.SetParName(2, "AP");
    f.FixParameter(2, ap);
    f.SetParLimits(1, -1, 1); // P is a polarization
    return h.Fit(&f, opt) == 0;
}

// Fits one slice's phi values (deg). With draw, the histogram and fit are drawn in the current pad.
static bool FitSlice(const std::vector<double>& phis, Slice& s, bool draw, const std::string& file = "", const std::string& cutNote = ""){ // file: the drawn page, numbers the title; cutNote: the cuts applied, named in the title
    TH1D h("hPolPhi", "", kPhiBins, 0, 360);
    TF1 f("fPol", "[0]*(1+[1]*[2]*cos(x*TMath::DegToRad()))", 0, 360);
    if (!FitPhi(phis, s.ap, h, f, draw ? "LQ" : "LQ0")) return false; // drawn pages keep the fit in the histogram, like the other plots' fit boxes
    s.nPhi = h.GetEntries(); s.p = f.GetParameter(1); s.pErr = f.GetParError(1);
    if (draw){
        h.SetTitle(PlotTitle(file, Form("Phi, theta %.2f to %.2f deg, AP = %.4f", s.lo, s.hi, s.ap) + (cutNote.empty() ? "" : ", " + cutNote) + std::string(";Angle (degrees);Count")).c_str());
        h.GetXaxis()->SetNdivisions(505);
        h.SetStats(0);
        h.SetMinimum(0); h.SetMaximum(1.1 * h.GetMaximum());
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
        box->Draw(); // before DrawCEPercent: the CE lost box goes under it
        DrawCEPercent(s.nCE, s.nCE != s.nPhi ? s.nPhi : -1, s.nCEAll, s.nReconAll); // recon protons fit: CE share and CE lost, CE fit: CE lost
    }
    return true;
}

// 9 of the fitted slices (evenly spaced over them) on one page, 3x3: each slice's phi histogram with its N0 (1 + AP P cos(phi)) fit, P and AP in
// the pad title. phis: the phi values (deg) of every slice, by slice index.
static void DrawSliceGrid(const std::vector<std::vector<double>>& phis, const std::vector<Slice>& fitted, const std::string& path, const std::string& what){
    if (fitted.empty()) return;
    TCanvas c("cSliceGrid", "", 1500, 1200);
    TPad top("top", "", 0, 0.94, 1, 1), body("body", "", 0, 0, 1, 0.94);
    top.Draw(); body.Draw();
    top.cd();
    TPaveText head(0, 0, 1, 1, "NDC");
    head.SetBorderSize(0); head.SetFillColor(0); head.SetTextSize(0.32);
    head.AddText(PlotTitle(path, "Phi fits N_{0} (1 + AP #times P #times cos#phi), AP fixed, per #theta slice, " + what).c_str());
    head.Draw();
    body.Divide(3, 3, 0.004, 0.004);
    const float titleSize = gStyle->GetTitleFontSize();
    gStyle->SetTitleFontSize(0.065); // readable pad titles at a third of the page
    const size_t nShow = std::min<size_t>(9, fitted.size());
    for (size_t i = 0; i < nShow; i++){
        const Slice& s = fitted[nShow == 1 ? 0 : i * (fitted.size() - 1) / (nShow - 1)];
        body.cd(i + 1);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.12);
        TH1D h(Form("hGrid%zu", i), "", kPhiBins, 0, 360);
        TF1 f("fGrid", "[0]*(1+[1]*[2]*cos(x*TMath::DegToRad()))", 0, 360);
        if (!FitPhi(phis[(int)std::lround(s.lo / kThetaStep)], s.ap, h, f, "LQ")) continue;
        h.SetTitle(Form("#theta %.2f to %.2f deg: P = %.3f #pm %.3f, AP = %.4f, N = %.0f;#phi (deg);Count", s.lo, s.hi, f.GetParameter(1), f.GetParError(1), s.ap, h.GetEntries()));
        h.SetStats(0);
        h.GetXaxis()->SetNdivisions(505);
        for (TAxis* a : {h.GetXaxis(), h.GetYaxis()}){a->SetTitleSize(0.055); a->SetLabelSize(0.05);}
        h.SetMinimum(0); h.SetMaximum(1.1 * h.GetMaximum());
        h.DrawCopy();
    }
    c.SaveAs(path.c_str());
    gStyle->SetTitleFontSize(titleSize);
}

void PlotPolarization(std::string energy = "fKEnergy"){ // energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the input trees and plot folder
    FastJit();
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    const std::string keLabel = (energy == "fDEnergy") ? "ZDC KE" : "true KE";
    const double pRef = ReadSettings().at("PPValue");
    std::set<std::string> statsStarted; // stats.csv files already truncated in this run

    TFile* ceFile = TFile::Open("cuts/cut_APandT.root", "READ");
    TFile* allFile = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree* all = allFile ? (TTree*)allFile->Get(kAngleTree.c_str()) : nullptr; // every fitted proton
    if (!ceFile || !all){printf("PlotPolarization: cut_APandT.root or AnglesRaw in cut_angles_all.root missing\n"); return;}
    const double nReconAll = (double)Columns(all, {"fEvent"})[0].size(); // every recon proton (AnglesRaw row), before any cut

    const std::vector<KeCut> kes = {KeCuts(keLabel).front()}; // CutFiles.h, the one cut ke_gt80: no KE and no residual cut folders

    // block: the block position set, which only takes the block AP and the block residual cut (kBlockCut, QualityCut.h); the recon set takes AP from CE and from recon
    struct PosSet {std::string dir, name, ceTree, theta, phi, ap, allTheta, allPhi; bool block;};
    const std::vector<PosSet> sets = {
        {"1_reconTheta", "recon theta", "CalcValsRe", "thetaRe", "phiRe", "apRe", "outTheta", "outPhi", false},
        {"2_reconBlockTheta", "recon block theta", "CalcValsReD", "thetaReD", "phiReD", "apReD", "outDTheta", "outDPhi", true}};
    const std::vector<std::string> typeTags = {"CE", "recon"}, typeNames = {"CE recon protons", "all recon protons"}; // CE and recon files share a folder, the tag ends the file name
    const int colors[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta, kOrange + 7, kCyan + 2, kGray + 2};
    const std::string deg = "*180/TMath::Pi()";

    TCanvas pageC("cPolPage", "", 1500, 1200); // phi pages and graphs, same canvas and margins as the PlotAngle.cpp histograms
    pageC.SetRightMargin(0.28); pageC.SetLeftMargin(0.14); pageC.SetBottomMargin(0.13);

    struct ApSrc {std::string tag, name, tree, branch; bool block; std::string sub;}; // sub: its subfolder of the KE cut folder ("" for the only AP of the block set). Where each proton's AP comes from; theta and phi come from the theta set (recon or recon block), rows line up between the trees
    const std::vector<ApSrc> apSrcs = {
        {"apReconBlock", "AP from recon block", "CalcValsReD" + suffix, "apReD", true, ""}};

    // Every slice fit of this run in one file (the same rows as the <ke>/slices_<tag>.csv files), with the folder names as columns
    gSystem->mkdir("plots/5_fomCalc", true);
    FILE* mega = fopen("plots/5_fomCalc/all_slice_fits.csv", "w");
    if (mega) fprintf(mega, "theta_set,slice_width_deg,ap_source,protons,ke_cut,theta_lo,theta_hi,N,AP_mean,P,P_err,pct_err\n");

    for (const double width : {kThetaMax / 100}) for (const auto& apSrc : apSrcs){ // 100 slices over 0 to kThetaMax
    kThetaStep = width; gApName = apSrc.name;
    const int nSlices = (int)std::round(kThetaMax / kThetaStep);
    if (!ceFile->Get(apSrc.tree.c_str())){printf("PlotPolarization: tree %s missing\n", apSrc.tree.c_str()); continue;}
    for (const auto& set : sets){
        if (apSrc.block != set.block) continue;
        TTree* ce = (TTree*)ceFile->Get((set.ceTree + suffix).c_str()); // theta and phi of the protons, binned in this theta
        if (!ce){printf("PlotPolarization: tree %s missing\n", (set.ceTree + suffix).c_str()); continue;}
        const bool apInOwnTree = apSrc.tree == set.ceTree + suffix;
        if (!apInOwnTree && !ce->GetFriend(apSrc.tree.c_str())) ce->AddFriend(apSrc.tree.c_str());

        std::vector<std::vector<Slice>> fits[2]; // [protons type][ke cut] -> fitted slices
        for (int type = 0; type < 2; type++){
            // one number per file in a folder: plot kind k (1 P, 2 pct error, 3, 4 total fits; PlotFomPolarization.cpp continues; 7 the slice fit grid), CE file 2k-1, recon file 2k
            auto fileNum = [type](int kind){return std::to_string(2 * (kind - 1) + 1 + type) + "_";};
            // every AP source and protons type of a theta set shares the folders, the tag <ap>_<CE or recon> ends the file name
            // one folder per KE cut holds everything of that cut: slice plots, total fits, fit range maps (PlotFomPolarization.cpp)
            const std::string base = "plots/5_fomCalc/" + set.dir + "/", tag = (apSrc.sub.empty() ? apSrc.tag + "_" : "") + typeTags[type]; // in an AP subfolder the AP is not repeated in the name

            fits[type].resize(kes.size());
            for (size_t k = 0; k < kes.size(); k++){
                gSystem->mkdir((base + kes[k].dir + "/" + apSrc.sub).c_str(), true);
                std::vector<std::vector<double>> ceCols, allCols;
                Gather(ce, set.theta + ":" + set.phi + ":" + (apInOwnTree ? "" : apSrc.tree + ".") + apSrc.branch, kes[k].dir, set.block, 3, ceCols);
                Gather(all, set.allTheta + deg + ":" + set.allPhi + deg, kes[k].dir, set.block, 2, allCols, "Angles" + suffix);

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

                FILE* out = fopen((base + kes[k].dir + "/" + apSrc.sub + "slices_" + tag + ".csv").c_str(), "w");
                if (out) fprintf(out, "theta_lo,theta_hi,N,AP_mean,P,P_err,pct_err\n");
                for (int b = 0; b < nSlices; b++){
                    if (apN[b] == 0) continue; // no AP without a CE proton in the slice
                    Slice s{b * kThetaStep, (b + 1) * kThetaStep, 0, apSum[b] / apN[b], 0, 0, apN[b]};
                    if (!FitSlice(type == 0 ? cePhi[b] : allPhi[b], s, false)) continue; // slice fits are not plotted, only the csv
                    fits[type][k].push_back(s);
                    if (out) fprintf(out, "%g,%g,%.0f,%.8g,%.8g,%.8g,%.8g\n", s.lo, s.hi, s.nPhi, s.ap, s.p, s.pErr, 100.0 * (s.p - pRef) / pRef);
                    if (mega) fprintf(mega, "%s,%g,%s,%s,%s,%g,%g,%.0f,%.8g,%.8g,%.8g,%.8g\n", set.dir.c_str(), width, apSrc.tag.c_str(), typeTags[type].c_str(), kes[k].dir.c_str(), s.lo, s.hi, s.nPhi, s.ap, s.p, s.pErr, 100.0 * (s.p - pRef) / pRef);
                }
                if (out) fclose(out);
                // 13, 14_phi_slice_fits_<tag>.pdf: some of the slice fits above on one page
                if (!set.block) DrawSliceGrid(type == 0 ? cePhi : allPhi, fits[type][k], base + kes[k].dir + "/" + apSrc.sub + fileNum(7) + "phi_slice_fits_" + tag + ".pdf",
                              typeNames[type] + ", " + apSrc.name + ", " + set.name + ", " + kes[k].label);

                // One fit over a whole theta range instead of slices: 0.5 to 2 deg and 0 to kThetaMax deg, AP = mean AP of the CE protons in the range.
                // Independent of the slice width: 5_fomCalc/<theta set>/<ke>/<ap>/5, 6_phi_fit_0p5to2deg_<tag>.pdf and 7, 8_phi_fit_0to5deg_<tag>.pdf, numbers in total_fits.csv.
                // Recon theta only, the block set has no total fits.
                if (!set.block){
                    for (const auto& range : std::vector<std::pair<double, double>>{{0.5, 2.0}, {0.0, kThetaMax}}){
                        const std::string totDir = base + kes[k].dir + "/" + apSrc.sub;
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
                        Slice t{range.first, range.second, 0, apSum2 / apN2, 0, 0, apN2, (double)ce->GetEntries(), nReconAll};
                        pageC.cd(); pageC.Clear();
                        const std::string totPdf = totDir + (range.first > 0 ? fileNum(3) + "phi_fit_0p5to2deg_" : fileNum(4) + "phi_fit_0to5deg_") + tag + ".pdf";
                        if (!FitSlice(phis, t, true, gSystem->BaseName(totPdf.c_str()), kes[k].label + (set.block ? ", block residual < 2 mm" : ""))) continue;
                        pageC.SaveAs(totPdf.c_str());
                        WriteStatsRow(totPdf, totHeader, Form("%s,%g,%g,%.0f,%.8g,%.8g,%.8g,%.8g", gSystem->BaseName(totPdf.c_str()), t.lo, t.hi, t.nPhi, t.ap, t.p, t.pErr, 100.0 * (t.p - pRef) / pRef), statsStarted, "total_fits.csv");
                    }
                }
            }

            // P vs theta (1_P_vs_theta) and the absolute percent error vs PPValue on a log axis (2_pct_error_vs_theta_log), one plot per KE cut,
            // one marker per fitted slice at its center.
            for (size_t k = 0; k < kes.size(); k++){
                if (fits[type][k].empty() || set.block) continue; // the block set has no P vs theta plots
                // <base>/<ke>/stats.csv: one row per plot in the folder, from the fitted slices that plot shows
                {
                    double nTot = 0, pSum = 0, pctSum = 0, pctMin = 1e30, pctMax = -1e30;
                    for (const auto& s : fits[type][k]){
                        const double pct = 100.0 * (s.p - pRef) / pRef;
                        nTot += s.nPhi; pSum += s.p; pctSum += pct; pctMin = std::min(pctMin, pct); pctMax = std::max(pctMax, pct);
                    }
                    const double nFit = fits[type][k].size();
                    const std::string row = Form("%.0f,%.0f,%.8g,%.8g,%.8g,%.8g", nFit, nTot, pSum / nFit, pctSum / nFit, pctMin, pctMax);
                    for (int kind : {1, 2}){
                        const std::string f = fileNum(kind) + (kind == 1 ? "P_vs_theta_" : "pct_error_vs_theta_log_") + tag + ".pdf";
                        WriteStatsRow(base + kes[k].dir + "/" + apSrc.sub + f, "plot,slices_fitted,N_in_slices,P_mean,pct_err_mean,pct_err_min,pct_err_max", f + "," + row, statsStarted);
                    }
                }
                for (const bool pctPlot : {false, true}){ // the error plot is log only, so it shows the absolute percent error (log can't show negative values)
                    const std::string file = (pctPlot ? fileNum(2) + "pct_error_vs_theta_log_" : fileNum(1) + "P_vs_theta_") + tag + ".pdf";
                    TGraph points;
                    // fixed axes, the same on every P plot (CE and recon, every KE cut and AP source) so they compare directly:
                    // theta 0 to kThetaMax, P -1 to 1 (the fit limits), |percent error| 0.01 to 1000 % on the log axis
                    const double xLo = 0, xHi = kThetaMax, yMin = pctPlot ? 1e-2 : -1.0, yMax = pctPlot ? 1e3 : 1.0;
                    for (const auto& s : fits[type][k]){
                        const double pct = 100.0 * (s.p - pRef) / pRef;
                        points.AddPoint(0.5 * (s.lo + s.hi), pctPlot ? std::max(std::fabs(pct), yMin) : s.p); // below the axis floor drawn on it
                    }
                    pageC.cd(); pageC.Clear();
                    const std::string title = std::string(pctPlot ? "Absolute percent error of fitted polarization vs PPValue, log scale" : "Fitted polarization P") + " vs #theta, " +
                                typeNames[type] + ", " + apSrc.name + ", " + set.name + ", " + kes[k].label + (set.block ? ", block residual < 2 mm" : "") + ", theta < 5 deg" +
                                ";" + std::string(set.block ? "recon block #theta" : "recon #theta") + " (deg);" + (pctPlot ? "|100 #times (P_{fit} - P_{true}) / P_{true}| (%)" : "P_{fit}");
                    TH1F frame("framePol", PlotTitle(file, title).c_str(), 1, xLo, xHi); // axes and title only, the data is drawn on top
                    frame.SetStats(0);
                    frame.SetMinimum(yMin); frame.SetMaximum(yMax);
                    frame.GetXaxis()->SetNdivisions(505);
                    frame.SetBinContent(1, pctPlot ? yMin * 1e-3 : yMin - 100.0); // far below the range so its line stays off the plot ("AXIS" would not draw the title)
                    pageC.SetLogy(pctPlot);
                    frame.Draw("HIST");
                    points.SetMarkerStyle(20); points.SetMarkerSize(1.2); points.DrawClone("P");
                    { // recon protons: the CE ones in the plotted slices; both: the CE lost vs every CE proton of the theta set, no cuts
                        double nCE = 0, nRec = 0;
                        for (const auto& s : fits[type][k]){nCE += s.nCE; nRec += s.nPhi;}
                        DrawCEPercent(nCE, type == 1 ? nRec : -1, ce->GetEntries(), nReconAll);
                    }
                    if (!pctPlot){ // the reference (true P), not on the log axis
                        TLine ref(xLo, pRef, xHi, pRef);
                        ref.SetLineStyle(2);
                        ref.DrawClone();
                    }
                    pageC.SaveAs((base + kes[k].dir + "/" + apSrc.sub + file).c_str());
                    pageC.SetLogy(0);
                }
            }
        }
    }
    }
    if (mega) fclose(mega);
}
