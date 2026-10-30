#include <tuple>
#include <map>
#include <ROOT/RDataFrame.hxx>
#include <set>
#include <TFile.h>
#include <TTree.h>

#include <TSystem.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDatime.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TMath.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "FitCore.h"
#include <TROOT.h>
#include "CutCache.h"
#include "CutFiles.h"


// One output file: tree, draw expression, title text, axis titles, whether to fit the core, file name.
struct PlotSpec {
    TTree* tree;
    std::string expr, title, axes;
    bool fit;
    std::string file;
    std::string opt = ""; // TTree::Draw option, "colz" for the 2D plots
    bool block = false; // block position quantity: the block residual cut too (the <tree>_block list of the KE cut, CutFiles.h)
    // phi distribution: N0 (1 + AP P cos(phi)) fit (FitCosPhi, FitCore.h), AP fixed to the mean apVar of the CE protons in apTree under the
    // same KE and theta cut (apTheta is that tree's theta in degrees)
    TTree* apTree = nullptr;
    std::string apVar = "", apTheta = "";
};

struct ThetaCut {std::string dir, op, label;}; // op is appended to the theta expression

// Draws every plot in specs once per KE cut without a theta cut, in <groupDir>/<ke.dir>/<spec file> (a spec file may carry its own subfolder,
// 1_theta/1_thetaCE.pdf), and once per theta cut where cutPath(spec file, theta cut folder) gives a path; it returns "" for specs that get no
// theta cuts (only the phi plots do). thetaExpr is the recon theta in degrees for the tree being drawn.
using CutPath = std::function<std::string(const std::string& file, const std::string& thDir)>;
static void DrawGroup(TCanvas* canvas, const std::string& groupDir, const std::vector<KeCut>& kes,
                      const std::vector<PlotSpec>& specs, const std::string& label, const std::string& baseCut,
                      const std::string& extraCut, const std::function<std::string(TTree*)>& thetaExpr, const CutPath& cutPath)
{
    const std::vector<ThetaCut> thetas = {
        {"th_all", "", ""},
        {"th_gt0p1", " > 0.1", "theta > 0.1 deg"},
        {"th_gt0p3", " > 0.3", "theta > 0.3 deg"}};

    std::set<std::string> fitFilesStarted; // fit_values.txt files already truncated in this call
    for (const auto& ke : kes) {
        for (const auto& th : thetas) {
            for (const auto& s : specs) {
                const std::string rel = th.op.empty() ? s.file : cutPath(s.file, th.dir);
                if (rel.empty()) continue;
                const std::string path = groupDir + "/" + ke.dir + "/" + rel;
                gSystem->mkdir(gSystem->GetDirName(path.c_str()), true);

                std::vector<std::string> parts;
                for (const std::string& c : {baseCut, extraCut}) {
                    if (!c.empty()) parts.push_back("(" + c + ")");
                }
                if (!th.op.empty()) parts.push_back("(" + thetaExpr(s.tree) + th.op + ")");
                std::string cut;
                for (size_t i = 0; i < parts.size(); i++) cut += (i ? " && " : "") + parts[i];

                std::string title = s.title + ", " + label;
                if (!ke.label.empty()) title += ", " + ke.label;
                if (!th.label.empty()) title += ", " + th.label;

                // CE protons before any cut (the entry list is not set yet), for the share the cuts lose
                const bool isAngles = s.tree->GetBranch("isFromCE"); // Angles tree = all recon protons, cut_APandT.root trees = CE protons only
                const double nCEAll = isAngles ? s.tree->GetEntries("isFromCE == 1") : s.tree->GetEntries();
                const bool anyCut = !cut.empty() || ke.dir != "ke_all" || s.block;
                double ap = 0;
                if (s.apTree) {
                    ApplyCut(s.apTree, ke.dir, s.block);
                    const Long64_t nAP = s.apTree->Draw(s.apVar.c_str(), th.op.empty() ? "1" : (s.apTheta + th.op).c_str(), "goff");
                    if (nAP > 0) ap = TMath::Mean(nAP, s.apTree->GetV1());
                    s.apTree->SetEntryList(nullptr);
                }
                ApplyCut(s.tree, ke.dir, s.block); // the KE (and block residual) cut, CutTrees.cpp
                if (cut.empty()) cut = "1";
                canvas->Clear();
                s.tree->Draw(s.expr.c_str(), cut.c_str(), s.opt.c_str());
                TH1* h = s.tree->GetHistogram();
                if (!h){s.tree->SetEntryList(nullptr); continue;}
                h->SetTitle(PlotTitle(s.file, title + ";" + s.axes).c_str());
                if (s.fit && h->GetEntries() >= 10) FitCore(h);
                if (s.apTree) FitCosPhi(h, ap);
                if (h->GetDimension() == 1 && s.expr.find("phi") != std::string::npos){h->SetMinimum(0); h->SetMaximum(1.1 * h->GetMaximum());} // every 1D phi plot, 0 to max
                if (isAngles) // the CE ones among the recon protons under the same cut
                    DrawCEPercent(s.tree->GetEntries(("(" + cut + ") && isFromCE == 1").c_str()), h->GetEntries(), nCEAll, s.tree->GetEntries()); // GetEntries(): every entry, the list ignored
                else if (anyCut)
                    DrawCEPercent(s.tree->GetEntries(cut.c_str()), -1, nCEAll);
                s.tree->SetEntryList(nullptr);
                WritePlotNumbers(h, path, fitFilesStarted);
                canvas->Update();
                canvas->SaveAs(path.c_str());
            }
        }
    }
}

const double totalNeutrons = ReadSettings().at("TotalNeutrons"); // settings.csv from the sim
// Theta bin edges (deg) for the FOM table, from the data: every pair of block (cell center) hit points of a reconstructed track gives an angle
// atan(hypot(dx, dy) / dz). The unique angles up to maxTheta are ordered and each gets one bin: edges sit midway between consecutive angles,
// the first edge is 0 and the last is maxTheta. keys: (fEvent, TrackID) of the reconstructed tracks (the Angles tree).
// FOM theta bins: fixed width from 0 to maxTheta (degrees), the last bin ends exactly at maxTheta
static std::vector<double> ThetaEdges(double maxTheta, double binWidth)
{
    std::vector<double> edges;
    for (int i = 0; i * binWidth < maxTheta; i++) edges.push_back(i * binWidth);
    edges.push_back(maxTheta);
    return edges;
}

// FOM table, csv only, one file per KE cut: <groupDir>/1_csv/<ke.dir>.csv. t is the tree with the recon values (sfx: "Re" true position, "ReD" block position);
// all is the tree of every proton from every event (cut_angles_all.root), allThetaExpr its recon theta in degrees for the same position type.
// Per theta bin: N = reconstructed protons in the bin (the Angles tree only holds tracks with a fit), N_CE = CE protons in the bin (the ones with AP), eff = N / total neutrons, mean AP of the CE protons,
// eff*|<AP>| (the FOM) and eff*<AP>^2 (the F^2 term). The last row is the sum over all bins.
// block: the block residual cut on both trees too (block position)
static void WriteFomCsv(const std::string& groupDir, const std::vector<KeCut>& kes, TTree* t, const std::string& sfx, TTree* all, const std::string& allThetaExpr, const std::vector<double>& thetaEdges,
                        bool block = false)
{
    const size_t nb = thetaEdges.size() - 1;
    auto binOf = [&](double theta){ // bin = [lo, hi), -1 outside
        auto it = std::upper_bound(thetaEdges.begin(), thetaEdges.end(), theta);
        return (it == thetaEdges.begin() || it == thetaEdges.end()) ? -1 : int(it - thetaEdges.begin() - 1);
    };
    for (const auto& ke : kes) {
        t->SetEstimate(t->GetEntries() + 1);
        ApplyCut(t, ke.dir, block);
        const Long64_t n = t->Draw(("theta" + sfx + ":ap" + sfx).c_str(), "", "goff");
        t->SetEntryList(nullptr);
        std::vector<double> countCE(nb, 0.0), countAll(nb, 0.0), sumAP(nb, 0.0);
        for (Long64_t i = 0; i < n; i++) {
            const int b = binOf(t->GetV1()[i]);
            if (b < 0) continue;
            countCE[b]++;
            sumAP[b] += t->GetV2()[i];
        }
        all->SetEstimate(all->GetEntries() + 1);
        ApplyCut(all, ke.dir, block);
        const Long64_t nAll = all->Draw(allThetaExpr.c_str(), "", "goff");
        all->SetEntryList(nullptr);
        for (Long64_t i = 0; i < nAll; i++) {
            const int b = binOf(all->GetV1()[i]);
            if (b >= 0) countAll[b]++;
        }
        gSystem->mkdir((groupDir + "/1_csv").c_str(), true);
        FILE* out = fopen((groupDir + "/1_csv/" + ke.dir + ".csv").c_str(), "w");
        if (!out) continue;
        fprintf(out, "theta_lo,theta_hi,N,N_CE,eff,mean_AP,eff_x_AP,eff_x_AP2\n");
        double totN = 0, totCE = 0, totAP = 0, totFom = 0, totFom2 = 0;
        for (size_t b = 0; b < nb; b++) {
            const double mean = countCE[b] > 0 ? sumAP[b] / countCE[b] : 0.0, eff = countAll[b] / totalNeutrons;
            fprintf(out, "%g,%g,%.0f,%.0f,%.8g,%.8g,%.8g,%.8g\n", thetaEdges[b], thetaEdges[b + 1], countAll[b], countCE[b], eff, mean, eff * std::fabs(mean), eff * mean * mean);
            totN += countAll[b]; totCE += countCE[b]; totAP += sumAP[b]; totFom += eff * std::fabs(mean); totFom2 += eff * mean * mean;
        }
        fprintf(out, "%g,%g,%.0f,%.0f,%.8g,%.8g,%.8g,%.8g\n", thetaEdges.front(), thetaEdges.back(), totN, totCE, totN / totalNeutrons, totCE > 0 ? totAP / totCE : 0.0, totFom, totFom2);
        fclose(out);
    }
}

// phi histogram (50 bins) of t with the polarization fit N(phi) = N0 (1 + P <A> cos(phi)): N0 = entries / bins and <A> (mean sim analyzing power apCE of the CE protons) are fixed, only P is fitted
// aExpr (the per-event analyzing power, e.g. "LastAP", all negative) switches to the A-weighted version: phi is filled with weight -A, so a bin holds sum(-A) (1 + P A_eff cos(phi)) with A_eff = sum(A^2) / sum(A)
// (N0 = sum(-A) / bins and A_eff fixed), and the line P = 2 sum(A cos(phi)) / sum(A^2) +- sqrt(2 / sum(A^2)) is printed.
static void PhiCosFit(TCanvas* canvas, TTree* t, const std::string& phiBranch, const std::string& cut, const std::string& name, const std::string& title, double aMean, const std::string& path, const std::string& aExpr = ""){
    const int nBins = 50;
    const bool weighted = !aExpr.empty();
    canvas->Clear(); canvas->SetRightMargin(0.25); canvas->SetTopMargin(weighted ? 0.2 : 0.17); // weighted: a third text line above the frame
    double aFit = aMean, nFit = 0, pSum = 0, pSumErr = 0;
    TH1* h = nullptr;
    if (weighted){
        t->SetEstimate(t->GetEntries() + 1);
        const Long64_t n = t->Draw((aExpr + ":" + phiBranch).c_str(), cut.c_str(), "goff");
        double sA = 0, sA2 = 0, sAc = 0;
        for (Long64_t k = 0; k < n; k++){const double a = t->GetV1()[k]; sA += a; sA2 += a * a; sAc += a * std::cos(t->GetV2()[k] * TMath::Pi() / 180);}
        if (n < 10 || sA == 0){canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1); return;}
        aFit = sA2 / sA; pSum = 2 * sAc / sA2; pSumErr = std::sqrt(2 / sA2);
        h = new TH1D(name.c_str(), "", nBins, 0, 360);
        h->Sumw2();
        t->Draw((phiBranch + ">>" + name).c_str(), ("-(" + aExpr + ")" + (cut.empty() ? "" : "*(" + cut + ")")).c_str(), "hist"); // "hist": a histogram with Sumw2 is drawn as error bars by default
        nFit = h->GetSumOfWeights();
    } else {
        t->Draw((phiBranch + ">>" + name + "(" + std::to_string(nBins) + ",0,360)").c_str(), cut.c_str());
        h = t->GetHistogram();
        if (h) nFit = h->GetEntries();
    }
    if (!h || h->GetEntries() < 10){canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1); return;}
    h->SetTitle((title + ";Angle (degrees);" + (weighted ? "Weighted count (weight -A)" : "Count")).c_str());
    TF1 fit("fPhiCos", "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
    fit.SetParNames("N0", "P", "A");
    fit.FixParameter(0, nFit / nBins); fit.FixParameter(2, aFit); fit.SetParameter(1, 0.5); fit.SetParLimits(1, -1, 1);
    h->Fit(&fit, "Q");
    h->SetMinimum(0); h->SetMaximum(1.1 * h->GetMaximum());
    h->SetName("Fit results"); // the first line of the stats box is the histogram name
    if (weighted) fit.DrawCopy("same"); // the "hist" draw option leaves out the fitted function
    canvas->Update();
    if (auto* st = (TPaveStats*)h->FindObject("stats")){st->SetX1NDC(0.76); st->SetX2NDC(0.995); st->SetY1NDC(0.55); st->SetY2NDC(0.83);}
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
    if (weighted){
        tx.DrawLatex(0.13, 0.90, "Fit  N(#phi) = N_{0} (1 + P #upoint A_{eff} cos#phi),   weight -A per event,   N_{0} = #Sigma(-A) / bins (fixed)");
        tx.DrawLatex(0.13, 0.865, Form("A_{eff} = #SigmaA^{2} / #SigmaA = %.4f (fixed),   P = %.2f #pm %.2f", aFit, fit.GetParameter(1), fit.GetParError(1)));
        tx.DrawLatex(0.13, 0.83, Form("P = 2 #SigmaA cos#phi / #SigmaA^{2} = %.2f #pm %.2f", pSum, pSumErr));
    } else {
        tx.DrawLatex(0.13, 0.90, "Fit  N(#phi) = N_{0} (1 + P #upoint #LTA#GT cos#phi),   N_{0} = Entries / bins (fixed)");
        tx.DrawLatex(0.13, 0.865, Form("#LTA#GT = %.4f (fixed, mean sim A (LastAP)),   P = %.2f #pm %.2f", aMean, fit.GetParameter(1), fit.GetParError(1)));
    }
    canvas->Modified(); canvas->Update();
    canvas->SaveAs(path.c_str());
    canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
}

void PlotAngle(std::string energy = "fKEnergy"){ // energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the input files and plot folder
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    gStyle->SetOptFit(1111); // probability, chi2/ndf, fit parameters with errors
    std::string foldermain = "plots/3_ReconCalcVals/";
    gSystem->mkdir(foldermain.c_str(), true);


    TFile* File = TFile::Open("cuts/cut_APandT.root", "READ");

    TTree *events = (TTree*)File->Get(("CE" + suffix).c_str());
    TTree *events2 = (TTree*)File->Get(("CalcValsRe" + suffix).c_str());
    TTree *events3 = (TTree*)File->Get(("CalcValsReD" + suffix).c_str());
    TTree *events4 = (TTree*)File->Get(("CalcValsReDIFF" + suffix).c_str());
    TTree *events5 = (TTree*)File->Get(("CalcValsReDDIFF" + suffix).c_str());

    // All five trees are filled in the same loop (CalcAPandT.cpp), so entries line up and the recon theta can
    // be used to cut every tree.
    for (TTree* t : {events, events3, events4, events5}) t->AddFriend(events2);
    // the diff trees also need the CE value (diff = recon - CE) for the percent error plots
    for (TTree* t : {events4, events5}) t->AddFriend(events);

    TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
    gPad->SetRightMargin(0.30);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); // room for the rotated Y title and bottom tick labels
    gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X"); // stats/fit box to the right, fewer x labels

    // Same x range on every theta distribution so plots compare directly: 0.023 deg bins, 0 to 2.3 deg (theta acceptance 0.04 rad)
    auto thetaHist = [](const std::string& name){return ">>" + name + "(100,0,2.3)";};
    const std::string ang = "Angle (degrees);Count";
    const std::string diff = "Angle Difference (degrees);Count";
    const std::string pct = "Percent Error (%);Count";
    const std::string thetaAP = "Theta (degrees);AP";
    const std::string phiAP = "Phi (degrees);AP";
    auto pctErr = [&](const std::string& diffBranch, const std::string& ceBranch, const std::string& hist){
        return "100*" + diffBranch + "/CE" + suffix + "." + ceBranch + ">>" + hist + "(100,-100,100)";
    };

    //CE protons, cut_APandT.root. Theta cuts use the reconstructed (tracker) theta, in degrees.
    // No cut, plus a gt (>) cut at 50-99% of the beam KE, each also with the 2 mm block residual cut (KeCuts, CutFiles.h).
    const std::string keLabel = (energy == "fDEnergy") ? "ZDC KE" : "true KE"; // which energy the cuts are on, shown in the titles
    const std::vector<KeCut> keCE = KeCuts(keLabel);
    std::vector<PlotSpec> specsCE = {
        {events,  "thetaCE" + thetaHist("hThetaCE"),  "Reconstructed CE protons, CE theta Distribution", ang, false, "1_theta/1_thetaCE.pdf"},
        {events2, "thetaRe" + thetaHist("hThetaRe"),  "Reconstructed CE protons, recon theta Distribution", ang, false, "1_theta/2_thetaRe.pdf"},
        {events4, "thetaReD>>hThetaDiff(100,-1,1)", "Reconstructed CE protons, difference between CE and Recon Theta", diff, true, "1_theta/4_thetaDiff.pdf"},
        {events3, "thetaReD" + thetaHist("hThetaBlock"), "Reconstructed CE protons, recon block theta Distribution", ang, false, "1_theta/3_thetaBlock.pdf"},
        {events5, "thetaReD>>hThetaBlockDiff(100,-1,1)", "Reconstructed CE protons, difference between CE and Recon Block Theta", diff, true, "1_theta/5_thetaBlockDiff.pdf"},
        {events,  "phiCE>>hPhiCE(100,0,360)",    "Reconstructed CE protons, CE phi Distribution", ang, false, "2_phi/1_phiCE.pdf"},
        {events2, "phiRe>>hPhiRe(100,0,360)",    "Reconstructed CE protons, recon phi Distribution", ang, false, "2_phi/2_phiRe.pdf"},
        {events4, "phiReD>>hPhiDiff(100,-180,180)",   "Reconstructed CE protons, difference between CE and Recon Phi", diff, true, "2_phi/4_phiDiff.pdf"},
        {events3, "phiReD>>hPhiBlock(100,0,360)",   "Reconstructed CE protons, recon block phi Distribution", ang, false, "2_phi/3_phiBlock.pdf"},
        {events5, "phiReD>>hPhiBlockDiff(100,-180,180)",   "Reconstructed CE protons, difference between CE and Recon Block Phi", diff, true, "2_phi/5_phiBlockDiff.pdf"},
        {events,  "tCE>>hTCE(100,0,0.05)", "Reconstructed CE protons, CE t Distribution", "t;Count", false, "3_t/1_tCE.pdf"},
        {events2, "tRe>>hTRe(100,0,0.05)", "Reconstructed CE protons, recon t Distribution", "t;Count", false, "3_t/2_tRe.pdf"},
        {events4, "tReD>>hTDiff(100,-0.05,0.05)", "Reconstructed CE protons, difference between CE and Recon t", "t Difference;Count", true, "3_t/4_tDiff.pdf"},
        {events3, "tReD>>hTBlock(100,0,0.05)", "Recon block t Distribution", "t;Count", false, "3_t/3_tBlock.pdf"},
        {events5, "tReD>>hTBlockDiff(100,-0.05,0.05)", "Difference between CE and Recon Block t", "t Difference;Count", true, "3_t/5_tBlockDiff.pdf"},
        {events,  "apCE>>hAPCE(100,-1,0)",     "Reconstructed CE protons, CE AP Distribution", "AP;Count", false, "4_ap/1_apCE.pdf"},
        {events2, "apRe>>hAPRe(100,-1,0)",     "Recon AP Distribution", "AP;Count", false, "4_ap/2_apRe.pdf"},
        {events4, "apReD>>hAPDiff(100,-0.15,0.1)",    "Difference between CE and Recon AP", "AP Difference;Count", true, "4_ap/4_apDiff.pdf"},
        {events3, "apReD>>hAPBlock(100,-1,0)",    "Recon block AP Distribution", "AP;Count", false, "4_ap/3_apBlock.pdf"},
        {events5, "apReD>>hAPBlockDiff(100,-0.15,0.1)",    "Difference between CE and Recon Block AP", "AP Difference;Count", true, "4_ap/5_apBlockDiff.pdf"},
        // percent error = 100 * (recon - CE) / CE, CE value from the friend CE tree; clipped to +-100%
        {events4, pctErr("thetaReD", "thetaCE", "hPctTheta"), "Percent Error of Recon Theta vs CE", pct, true, "1_theta/6_thetaPct.pdf"},
        {events5, pctErr("thetaReD", "thetaCE", "hPctThetaB"), "Percent Error of Recon Block Theta vs CE", pct, true, "1_theta/6_thetaBlockPct.pdf"},
        {events4, pctErr("phiReD", "phiCE", "hPctPhi"), "Percent Error of Recon Phi vs CE", pct, true, "2_phi/6_phiPct.pdf"},
        {events5, pctErr("phiReD", "phiCE", "hPctPhiB"), "Percent Error of Recon Block Phi vs CE", pct, true, "2_phi/7_phiBlockPct.pdf"},
        {events4, pctErr("tReD", "tCE", "hPctT"), "Percent Error of Recon t vs CE", pct, true, "3_t/6_tPct.pdf"},
        {events5, pctErr("tReD", "tCE", "hPctTB"), "Percent Error of Recon Block t vs CE", pct, true, "3_t/7_tBlockPct.pdf"},
        {events4, pctErr("apReD", "apCE", "hPctAP"), "Percent Error of Recon AP vs CE", pct, true, "4_ap/6_apPct.pdf"},
        {events5, pctErr("apReD", "apCE", "hPctAPB"), "Percent Error of Recon Block AP vs CE", pct, true, "4_ap/7_apBlockPct.pdf"},
        // angle (degrees) vs AP; x is AP, y is the angle
        {events,  "apCE:thetaCE>>hThetaAPCE(100,0,2.3,100,-1,0)",    "Reconstructed CE protons, CE theta vs AP",          thetaAP, false, "5_angle_vs_ap/1_thetaVsAP/1_thetaVsApCE.pdf", "colz"},
        {events2, "apRe:thetaRe>>hThetaAPRe(100,0,2.3,100,-1,0)",    "Recon theta vs AP",       thetaAP, false, "5_angle_vs_ap/1_thetaVsAP/2_thetaVsApRe.pdf", "colz"},
        {events3, "apReD:thetaReD>>hThetaAPBlock(100,0,2.3,100,-1,0)", "Recon block theta vs AP", thetaAP, false, "5_angle_vs_ap/1_thetaVsAP/3_thetaVsApBlock.pdf", "colz"},
        {events,  "apCE:phiCE>>hPhiAPCE(100,0,360,100,-1,0)",        "Reconstructed CE protons, CE phi vs AP",            phiAP, false, "5_angle_vs_ap/2_phiVsAP/1_phiVsApCE.pdf", "colz"},
        {events2, "apRe:phiRe>>hPhiAPRe(100,0,360,100,-1,0)",        "Recon phi vs AP",         phiAP, false, "5_angle_vs_ap/2_phiVsAP/2_phiVsApRe.pdf", "colz"},
        {events3, "apReD:phiReD>>hPhiAPBlock(100,0,360,100,-1,0)",   "Recon block phi vs AP",   phiAP, false, "5_angle_vs_ap/2_phiVsAP/3_phiVsApBlock.pdf", "colz"}};

    // AP difference (recon - CE) and its percent error against the CE angle, for the true position recon (events4) and the block recon (events5)
    {
        struct DiffSet {TTree* tree; std::string ap, name, what;};
        const std::vector<DiffSet> dsets = {{events4, "apReD", "Re", "Recon"}, {events5, "apReD", "Block", "Recon block"}};
        struct AngleAxis {std::string ce, dir, name, label, xrange;};
        const std::vector<AngleAxis> axesSets = {{"thetaCE", "1_thetaVsAP", "theta", "Theta (degrees)", "0,2.3"}, {"phiCE", "2_phiVsAP", "phi", "Phi (degrees)", "0,360"}};
        for (const auto& ax : axesSets) {
            // after the three 2D plots already in the folder: AP differences first (4-5), then percent errors (6-7)
            for (size_t i = 0; i < dsets.size(); i++) {
                const auto& d = dsets[i];
                const std::string x = "CE" + suffix + "." + ax.ce;
                const std::string apCE = "CE" + suffix + ".apCE";
                const std::string base = "5_angle_vs_ap/" + ax.dir + "/";
                specsCE.push_back({d.tree, d.ap + ":" + x + ">>h" + ax.name + "ApDiff" + d.name + "(100," + ax.xrange + ",100,-0.15,0.1)",
                                   "AP difference (" + d.what + " - CE) vs CE " + ax.name, ax.label + ";AP Difference", false,
                                   base + std::to_string(4 + i) + "_" + ax.name + "ApDiff" + d.name + ".pdf", "colz"});
                specsCE.push_back({d.tree, "100*" + d.ap + "/" + apCE + ":" + x + ">>h" + ax.name + "ApPct" + d.name + "(100," + ax.xrange + ",100,-100,100)",
                                   "AP percent error (" + d.what + " vs CE) vs CE " + ax.name, ax.label + ";AP Percent Error (%)", false,
                                   base + std::to_string(6 + i) + "_" + ax.name + "ApPct" + d.name + ".pdf", "colz"});
            }
        }
    }
    // block residual cut on every block position plot (CalcValsReD, CalcValsReDDIFF), only there (QualityCut.h)
    for (auto& s : specsCE)
        if (s.tree == events3 || s.tree == events5) s.block = true;

    // FOM from the true position recon (events2, "Re") and from the block position recon (events3, "ReD")
    const std::string fomDir = "plots/5_fomCalc/";

    //All protons and CE protons, theta/phi and block theta/phi
    TFile* FileAll = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree *all = (TTree*)FileAll->Get(("Angles" + suffix).c_str());
    TTree *allCE = (TTree*)FileAll->Get(("AnglesCE" + suffix).c_str());

    printf("[PlotAngle] FOM setup %s\n", TDatime().AsSQLString()); fflush(stdout);
    // up to the largest recon theta of any reconstructed proton (true or block), so every one lands in a bin
    const double maxTheta = std::max(all->GetMaximum("outTheta"), all->GetMaximum("outDTheta")) * 180.0 / TMath::Pi() + 1e-6;
    const std::vector<double> thetaEdges = ThetaEdges(maxTheta, 0.001);

    printf("[PlotAngle] %zu theta bins, writing FOM csv %s\n", thetaEdges.size() - 1, TDatime().AsSQLString()); fflush(stdout);
    // FOM tables: reconstructed protons from all events (N) and the CE protons with AP (N_CE) per theta bin
    const std::vector<KeCut> keFom = {keCE.front()}; // the FOM tables have no KE and no residual cuts
    WriteFomCsv(fomDir + "1_reconTheta", keFom, events2, "Re", all, "outTheta*180/TMath::Pi()", thetaEdges);
    WriteFomCsv(fomDir + "2_reconBlockTheta", keFom, events3, "ReD", all, "outDTheta*180/TMath::Pi()", thetaEdges, true);

    printf("[PlotAngle] FOM csv done %s\n", TDatime().AsSQLString()); fflush(stdout);
    // every recon proton, no theta acceptance cut, so the theta plots run up to the largest recon theta
    auto thetaHistAll = [&](const std::string& name){return Form(">>%s(100,0,%g)", name.c_str(), std::ceil(maxTheta));};
    const std::string deg = "*180/TMath::Pi()";

    std::vector<KeCut> keAll; // the KE cuts without the residual cut twins, for 1_CEvalues
    for (const auto& ke : keCE) if (ke.dir.size() < std::string(kRes2Dir).size() || ke.dir.compare(ke.dir.size() - std::string(kRes2Dir).size(), std::string::npos, kRes2Dir) != 0) keAll.push_back(ke);

    // 3_ReconCalcVals only holds detector energy plots, so the true KE pass stops after its FOM tables
    if (energy != "fDEnergy") return;
    // theta cuts only on the phi plots: <ke>/2_phi/<th>/<file>
    DrawGroup(canvas, foldermain.substr(0, foldermain.size() - 1), {keAll.front()}, specsCE, "CE recon protons", "", "",
              [&](TTree* t){return std::string(t == events2 ? "thetaRe" : "CalcValsRe" + suffix + ".thetaRe");},
              [](const std::string& f, const std::string& th){return f.rfind("2_phi/", 0) == 0 ? "2_phi/" + th + "/" + f.substr(6) : std::string();});
    // phi fits of the CE protons (CE angles) and of the same events at the recon angles
    {
        const std::string phiDir = foldermain.substr(0, foldermain.size() - 1) + "/" + keAll.front().dir + "/";
        gSystem->mkdir(phiDir.c_str(), true);
        // CE protons: every carbon-block CE of the sim (ChargeExchange_values, lab phi of the proton at the reaction), not only the ones with a recon track
        // every carbon-block CE of the sim: the cut cache of CarbonCeCut (seed cut_CE.root), so every row passes and the cut is empty
        const CacheTree cacheCE = GetCutTree("output0_t*.root", "ChargeExchange_values", {CarbonCeCut()});
        TTree* rawCE = cacheCE.tree;
        gROOT->cd();
        const std::string ceCut = "";
        const Long64_t nRaw = rawCE->Draw("LastAP", ceCut.c_str(), "goff");
        double aRaw = 0; for (Long64_t k = 0; k < nRaw; k++) aRaw += rawCE->GetV1()[k] / nRaw;
        PhiCosFit(canvas, rawCE, "(LastPhiLab < 0 ? LastPhiLab + 2 * TMath::Pi() : LastPhiLab) * 180 / TMath::Pi()", ceCut, "hPhiCEfit", "All CE protons (every carbon-block CE), CE phi Distribution", aRaw, phiDir + "1_phiCE_cosfit.pdf");
        // the same events at the recon angles (CE protons with a recon track), A is the mean sim A of those events
        const Long64_t nA = events->Draw("apCE", "", "goff");
        double aMean = 0; for (Long64_t k = 0; k < nA; k++) aMean += events->GetV1()[k] / nA;
        PhiCosFit(canvas, events, "phiRe", "", "hPhiRefit", "Reconstructed CE protons, recon phi Distribution", aMean, phiDir + "2_phiRe_cosfit.pdf");
        // the same fit at the block position recon phi (phiReD), A is the mean recon block A of those events
        const Long64_t nAD = events3->Draw("apReD", "", "goff");
        double aMeanD = 0; for (Long64_t k = 0; k < nAD; k++) aMeanD += events3->GetV1()[k] / nAD;
        PhiCosFit(canvas, events3, "phiReD", "", "hPhiReDfit", "Reconstructed CE protons, recon block phi Distribution", aMeanD, phiDir + "5_phiReD_cosfit.pdf");
    }
    printf("[PlotAngle] CE plots done %s\n", TDatime().AsSQLString()); fflush(stdout);
}
