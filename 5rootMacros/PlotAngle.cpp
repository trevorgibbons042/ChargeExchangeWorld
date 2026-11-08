#include <tuple>
#include <map>
#include <ROOT/RDataFrame.hxx>
#include <set>
#include <TFile.h>
#include <TTree.h>

#include <TSystem.h>
#include <TStyle.h>
#include <TColor.h>
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
#include "CutOpen.h"
#include "CutFiles.h"
#include "TreeColumns.h"
#include "PhiCosFit.h"


// One output file: tree, draw expression, title text, axis titles, whether to fit the core, file name.
struct PlotSpec {
    TTree* tree;
    std::string expr, title, axes;
    bool fit;
    std::string file;
    std::string opt = ""; // draw option, "colz" for the 2D plots
};

struct ThetaCut {std::string dir, op, label; double min;}; // op: the theta cut as text (empty: no cut), min: its lower limit in degrees

// Draws every plot in specs once per KE cut without a theta cut, in <groupDir>/<ke.dir>/<spec file> (a spec file may carry its own subfolder,
// 1_theta/1_thetaCE.pdf), and once per theta cut where cutPath(spec file, theta cut folder) gives a path; it returns "" for specs that get no
// theta cuts (only the phi plots do). thetaExpr is the RECONTHETA_ESTNEU in degrees for the tree being drawn.
// Per KE cut and tree one Columns call (TreeColumns.h) reads the draw expressions of all the specs of that tree and its theta; the histograms are filled from those columns
// in entry order, so a plot is the one TTree::Draw made, the theta cut is a mask over the rows.
using CutPath = std::function<std::string(const std::string& file, const std::string& thDir)>;
// CE protons are red: a white-to-red palette for the 2D plots (the other plots keep viridis, set in PlotAngle())
static void RedPalette(){
    double stops[2] = {0.0, 1.0}, red[2] = {1.0, 0.8}, green[2] = {1.0, 0.0}, blue[2] = {1.0, 0.0};
    TColor::CreateGradientColorTable(2, stops, red, green, blue, 255);
    gStyle->SetNumberContours(255);
}
static void ViridisPalette(){gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255);}

static void DrawGroup(TCanvas* canvas, const std::string& groupDir, const std::vector<KeCut>& kes,
                      const std::vector<PlotSpec>& specs, const std::string& label, const std::string& baseCut,
                      const std::string& extraCut, const std::function<std::string(TTree*)>& thetaExpr, const CutPath& cutPath)
{
    const std::vector<ThetaCut> thetas = {
        {"th_all", "", "", 0},
        {"th_gt0p1", " > 0.1", "theta > 0.1 deg", 0.1},
        {"th_gt0p3", " > 0.3", "theta > 0.3 deg", 0.3}};

    std::string filter; // the cut that is the same for every plot
    for (const std::string& c : {baseCut, extraCut}) {
        if (!c.empty()) filter += (filter.empty() ? "" : " && ") + std::string("(") + c + ")";
    }

    // the columns to read per tree: the draw expressions of its specs, and the theta when a theta cut is drawn of that tree
    using Key = TTree*;
    std::map<Key, std::vector<std::string>> exprsOf;
    auto slot = [&](const Key& k, const std::string& e){
        auto& v = exprsOf[k];
        const size_t at = std::find(v.begin(), v.end(), e) - v.begin();
        if (at == v.size()) v.push_back(e);
        return (int)at;
    };
    std::vector<DrawSpec> draws;
    std::vector<std::vector<int>> specCol(specs.size());
    std::vector<int> specTheta(specs.size(), -1);
    for (size_t i = 0; i < specs.size(); i++) {
        const Key k = specs[i].tree;
        draws.push_back(SplitDraw(specs[i].expr));
        for (const std::string& e : draws[i].exprs) specCol[i].push_back(slot(k, e));
        for (const auto& th : thetas)
            if (!th.op.empty() && !cutPath(specs[i].file, th.dir).empty()) specTheta[i] = slot(k, thetaExpr(specs[i].tree));
    }

    std::set<std::string> fitFilesStarted; // fit_values.txt files already truncated in this call
    for (const auto& ke : kes) {
        std::map<Key, std::vector<std::vector<double>>> cols; // one vector per expression, one entry per row under the KE cut
        std::map<TTree*, double> nCEAll; // CE protons before any cut (the entry list is not set yet), for the share the cuts lose
        for (const auto& [key, exprs] : exprsOf) {
            nCEAll[key] = key->GetEntries(); // GetEntries(): every entry, the list ignored
            if (ke.dir.empty()) { // no KE cut: no entry list
                cols[key] = Columns(key, exprs, filter);
                continue;
            }
            ApplyCut(key, ke.dir); // the KE cut, CutPlots.cpp (trees)
            cols[key] = Columns(key, exprs, filter);
            key->SetEntryList(nullptr);
        }
        for (const auto& th : thetas) {
            for (size_t i = 0; i < specs.size(); i++) {
                const auto& s = specs[i];
                const std::string rel = th.op.empty() ? s.file : cutPath(s.file, th.dir);
                if (rel.empty()) continue;
                const std::string path = groupDir + "/" + (ke.dir.empty() ? "" : ke.dir + "/") + rel; // an empty ke.dir (no KE cut) puts the plots straight in groupDir
                gSystem->mkdir(gSystem->GetDirName(path.c_str()), true);

                std::string title = s.title + ", " + label;
                if (!ke.label.empty()) title += ", " + ke.label;
                if (!th.label.empty()) title += ", " + th.label;

                const auto& c = cols.at(s.tree);
                canvas->Clear();
                TH1* h = MakeHist(draws[i].target); // after the Clear: a histogram of the same name from the last theta cut leaves the pad first
                if (!h) continue;
                const std::vector<double>& theta = c[th.op.empty() ? 0 : specTheta[i]];
                const std::vector<double>& x = c[specCol[i].back()];
                const std::vector<double>* y = specCol[i].size() > 1 ? &c[specCol[i].front()] : nullptr; // "y:x": y is the first expression
                Long64_t nPass = 0; // rows that pass the theta cut
                std::vector<double> xFit; // the values the unbinned fit (FitCore) takes: the rows that pass, no bins
                for (size_t r = 0; r < x.size(); r++) {
                    if (!th.op.empty() && !(theta[r] > th.min)) continue;
                    nPass++;
                    if (y) ((TH2*)h)->Fill(x[r], (*y)[r]);
                    else {h->Fill(x[r]); if (s.fit) xFit.push_back(x[r]);}
                }
                const bool is1D = h->GetDimension() == 1;
                if (is1D) h->SetLineColor(kRed); else RedPalette(); // CE protons are red
                h->Draw(s.opt.c_str()); // the histogram is the pad's current object, as after TTree::Draw
                h->SetTitle(PlotTitle(s.file, title + ";" + s.axes).c_str());
                if (s.fit && h->GetEntries() >= 10) FitCore(h, xFit);
                if (is1D) for (TObject* o : *h->GetListOfFunctions()) if (TF1* f = dynamic_cast<TF1*>(o)) f->SetLineColor(kBlack); // a red histogram would hide a red fit
                if (h->GetDimension() == 1 && s.expr.find("phi") != std::string::npos){h->SetMinimum(0); h->SetMaximum(1.1 * h->GetMaximum());} // every 1D phi plot, 0 to max
                if (h->GetDimension() == 1 && h->GetXaxis()->GetXmin() == -180 && h->GetXaxis()->GetXmax() == 180) h->GetXaxis()->SetNdivisions(508, kFALSE); // phi differences, -180 to 180 every 45 deg
                DrawCEPercent(nPass, -1, nCEAll.at(s.tree)); // the cuts of the plot (KE cut of ke, theta) are all counted as lost
                WritePlotNumbers(h, path, fitFilesStarted);
                canvas->Update();
                canvas->SaveAs(path.c_str());
                if (!is1D) ViridisPalette();
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

// N_particles: every fitted track of any pdg in the bin (N is the pdg 2212 ones).
// FOM table, csv only, one file per KE cut: <groupDir>/1_csv/<ke.dir>.csv. t is the tree with the recon values (sfx: "Re" true position, "ReD" block position);
// all is the tree of every proton from every event (AnglesRaw of cut_angles_all.root; the cut list "Angles" + allSfx keeps those that pass the KE cut), allThetaExpr its RECONTHETA_ESTNEU in degrees for the same position type.
// Per theta bin: N = reconstructed protons in the bin (the Angles tree only holds tracks with a fit), N_CE = CE protons in the bin (the ones with AP), eff = N / total neutrons, mean AP of the CE protons,
// eff*|<AP>| (the FOM) and eff*<AP>^2 (the F^2 term). The last row is the sum over all bins.
static void WriteFomCsv(const std::string& groupDir, const std::vector<KeCut>& kes, TTree* t, const std::string& sfx, TTree* all, const std::string& allSfx, const std::string& allThetaExpr, const std::vector<double>& thetaEdges)
{
    const size_t nb = thetaEdges.size() - 1;
    auto binOf = [&](double theta){ // bin = [lo, hi), -1 outside
        auto it = std::upper_bound(thetaEdges.begin(), thetaEdges.end(), theta);
        return (it == thetaEdges.begin() || it == thetaEdges.end()) ? -1 : int(it - thetaEdges.begin() - 1);
    };
    for (const auto& ke : kes) {
        ApplyCut(t, ke.dir);
        const auto ce = ColumnsUnordered(t, {"theta" + sfx, "ap" + sfx}); // theta, AP of every CE proton under the cut; the bin counts do not depend on the row order
        t->SetEntryList(nullptr);
        std::vector<double> countCE(nb, 0.0), countAll(nb, 0.0), countPart(nb, 0.0), sumAP(nb, 0.0); // countAll: recon protons (PDG 2212), countPart: every fitted track (any pdg)
        // the binning loops run on all threads (MapRows, FitCore.h): each row range counts into its own arrays, added together after. The counts are exact; sumAP is added in a
        // different order than a serial loop, which changes it by about 1e-16 relative.
        struct BinCounts {std::vector<double> a, b;};
        const auto ceParts = MapRows((Long64_t)ce[0].size(), [&](Long64_t lo, Long64_t hi){
            BinCounts p; p.a.assign(nb, 0.0); p.b.assign(nb, 0.0); // a: CE count, b: sum of AP
            for (Long64_t i = lo; i < hi; i++) {
                const int bin = binOf(ce[0][i]);
                if (bin < 0) continue;
                p.a[bin]++;
                p.b[bin] += ce[1][i];
            }
            return p;
        });
        for (const auto& p : ceParts) for (size_t b = 0; b < nb; b++) {countCE[b] += p.a[b]; sumAP[b] += p.b[b];}
        ApplyCut(all, ke.dir, "Angles" + allSfx);
        const auto allCols = ColumnsUnordered(all, {allThetaExpr, "PDG"}); const auto& allTheta = allCols[0]; // AnglesRaw has many clusters: read on all threads
        all->SetEntryList(nullptr);
        const auto allParts = MapRows((Long64_t)allTheta.size(), [&](Long64_t lo, Long64_t hi){
            BinCounts p; p.a.assign(nb, 0.0); p.b.assign(nb, 0.0); // a: every fitted track, b: the protons (PDG 2212)
            for (Long64_t i = lo; i < hi; i++) {
                const int bin = binOf(allTheta[i]);
                if (bin < 0) continue;
                p.a[bin]++;
                if (allCols[1][i] == 2212) p.b[bin]++;
            }
            return p;
        });
        for (const auto& p : allParts) for (size_t b = 0; b < nb; b++) {countPart[b] += p.a[b]; countAll[b] += p.b[b];}
        gSystem->mkdir((groupDir + "/1_csv").c_str(), true);
        FILE* out = fopen((groupDir + "/1_csv/" + ke.dir + ".csv").c_str(), "w");
        if (!out) continue;
        fprintf(out, "theta_lo,theta_hi,N,N_CE,eff,mean_AP,eff_x_AP,eff_x_AP2,N_particles\n");
        double totPart = 0, totN = 0, totCE = 0, totAP = 0, totFom = 0, totFom2 = 0;
        for (size_t b = 0; b < nb; b++) {
            const double mean = countCE[b] > 0 ? sumAP[b] / countCE[b] : 0.0, eff = countAll[b] / totalNeutrons;
            fprintf(out, "%g,%g,%.0f,%.0f,%.8g,%.8g,%.8g,%.8g,%.0f\n", thetaEdges[b], thetaEdges[b + 1], countAll[b], countCE[b], eff, mean, eff * std::fabs(mean), eff * mean * mean, countPart[b]);
            totPart += countPart[b]; totN += countAll[b]; totCE += countCE[b]; totAP += sumAP[b]; totFom += eff * std::fabs(mean); totFom2 += eff * mean * mean;
        }
        fprintf(out, "%g,%g,%.0f,%.0f,%.8g,%.8g,%.8g,%.8g,%.0f\n", thetaEdges.front(), thetaEdges.back(), totN, totCE, totN / totalNeutrons, totCE > 0 ? totAP / totCE : 0.0, totFom, totFom2, totPart);
        fclose(out);
    }
}

void PlotAngle(std::string energy = "fKEnergy"){ // energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the input files and plot folder
    FastJit();
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    gStyle->SetOptFit(11); // fit parameters with errors (the fits are unbinned: no chi2 or probability)
    std::string foldermain = "plots/3_ReconCalcVals/";
    gSystem->mkdir(foldermain.c_str(), true);


    TFile* File = TFile::Open("cuts/cut_APandT.root", "READ");

    TTree *events = (TTree*)File->Get(("CE" + suffix).c_str());
    TTree *events2 = (TTree*)File->Get(("CalcValsRe" + suffix).c_str());
    TTree *events3 = (TTree*)File->Get(("CalcValsReD" + suffix).c_str());
    TTree *events4 = (TTree*)File->Get(("CalcValsReDIFF" + suffix).c_str());
    TTree *events5 = (TTree*)File->Get(("CalcValsReDDIFF" + suffix).c_str());

    // All five trees are filled in the same loop (CalcAPandT.cpp), so entries line up and the RECONTHETA_ESTNEU can
    // be used to cut every tree.
    for (TTree* t : {events, events3, events4, events5}) t->AddFriend(events2);
    // the diff trees also need the CE value (diff = recon - CE) for the percent error plots
    for (TTree* t : {events4, events5}) t->AddFriend(events);
    events3->AddFriend(events); // the block recon vs CETHETA_NEUFRAME 2D plot (5_thetaBlockDiff.pdf) needs the CE value next to the block value

    TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
    gPad->SetRightMargin(0.30);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); // room for the rotated Y title and bottom tick labels
    gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X"); // stats/fit box to the right, fewer x labels

    // Same x range on every theta distribution so plots compare directly: 100 bins, 0 to ThetaAxisMax deg (2.3 deg at 25 GeV/c, scaled with 1 / p, FitCore.h)
    auto thetaHist = [](const std::string& name){return std::string(Form(">>%s(100,0,%g)", name.c_str(), ThetaAxisMax()));};
    const std::string ang = "Angle (degrees);Count";
    const std::string diff = "Angle Difference (degrees);Count";
    const std::string pct = "Percent Error (%);Count";
    auto pctErr = [&](const std::string& diffBranch, const std::string& ceBranch, const std::string& hist){
        // 100 * (a / b), bracketed: TTree::Draw divides first, and a * 100 / b would put values at exactly -100% a rounding step off
        return "100*(" + diffBranch + "/CE" + suffix + "." + ceBranch + ")>>" + hist + "(100,-100,100)";
    };

    //CE protons, cut_APandT.root. Theta cuts use the reconstructed (tracker) theta, in degrees.
    // No cut, plus a gt (>) cut at 50-99% of the beam KE (KeCuts, CutFiles.h). No residual cut.
    const std::string keLabel = (energy == "fDEnergy") ? "ZDC KE" : "true KE"; // which energy the cuts are on, shown in the titles
    const std::vector<KeCut> keCE = KeCuts(keLabel);
    auto reDiff = [&](const std::string& re, const std::string& block){return "(CalcValsRe" + suffix + "." + re + "-" + block + ")";}; // recon minus recon block
    auto wrap360 = [](const std::string& d){return "(" + d + "-360*(" + d + ">180)+360*(" + d + "<-180))";}; // phi difference into -180..180
    std::vector<PlotSpec> specsCE = {
        {events,  "thetaCE" + thetaHist("hThetaCE"),  "Reconstructed CE protons, CETHETA_NEUFRAME Distribution", ang, false, "1_theta/1_thetaCE.pdf"},
        {events2, "thetaRe" + thetaHist("hThetaRe"),  "Reconstructed CE protons, RECONTHETA_ESTNEU Distribution", ang, false, "1_theta/2_thetaRe.pdf"},
        // theta matrices, x is the truth (CETHETA_NEUFRAME) and y the reconstructed value; the file names keep the old Diff names
        {events4, Form("CalcValsRe%s.thetaRe:CE%s.thetaCE>>hThetaDiff(100,0,%g,100,0,%g)", suffix.c_str(), suffix.c_str(), ThetaAxisMax(), ThetaAxisMax()), "Reconstructed CE protons, RECONTHETA_ESTNEU vs CETHETA_NEUFRAME", "CETHETA_NEUFRAME (degrees);RECONTHETA_ESTNEU (degrees)", false, "1_theta/4_thetaDiff.pdf", "colz"},
        {events3, "thetaReD" + thetaHist("hThetaBlock"), "Reconstructed CE protons, RECONTHETABL_ESTNEU Distribution", ang, false, "1_theta/3_thetaBlock.pdf"},
        {events3, Form("thetaReD:CE%s.thetaCE>>hThetaBlockDiff(100,0,%g,100,0,%g)", suffix.c_str(), ThetaAxisMax(), ThetaAxisMax()), "Reconstructed CE protons, RECONTHETABL_ESTNEU vs CETHETA_NEUFRAME", "CETHETA_NEUFRAME (degrees);RECONTHETABL_ESTNEU (degrees)", false, "1_theta/5_thetaBlockDiff.pdf", "colz"},
        {events,  "phiCE>>hPhiCE(100,0,360)",    "Reconstructed CE protons, CEPHI_NEUFRAME Distribution", ang, false, "2_phi/1_phiCE.pdf"},
        {events2, "phiRe>>hPhiRe(100,0,360)",    "Reconstructed CE protons, RECONPHI_ESTNEU Distribution", ang, false, "2_phi/2_phiRe.pdf"},
        {events4, "phiReD>>hPhiDiff(100,-180,180)",   "Reconstructed CE protons, difference between CEPHI_NEUFRAME and RECONPHI_ESTNEU", diff, true, "2_phi/4_phiDiff.pdf"},
        {events3, "phiReD>>hPhiBlock(100,0,360)",   "Reconstructed CE protons, RECONPHIBL_ESTNEU Distribution", ang, false, "2_phi/3_phiBlock.pdf"},
        {events5, "phiReD>>hPhiBlockDiff(100,-180,180)",   "Reconstructed CE protons, difference between CEPHI_NEUFRAME and RECONPHIBL_ESTNEU", diff, true, "2_phi/5_phiBlockDiff.pdf"},
        {events,  "tCE>>hTCE(100,0,0.05)", "Reconstructed CE protons, t from CETHETA_NEUFRAME Distribution", "t;Count", false, "3_t/1_tCE.pdf"},
        {events2, "tRe>>hTRe(100,0,0.05)", "Reconstructed CE protons, t from RECONTHETA_ESTNEU Distribution", "t;Count", false, "3_t/2_tRe.pdf"},
        {events4, "tReD>>hTDiff(100,-0.05,0.05)", "Reconstructed CE protons, difference between t from CETHETA_NEUFRAME and t from RECONTHETA_ESTNEU", "t Difference;Count", true, "3_t/4_tDiff.pdf"},
        {events3, "tReD>>hTBlock(100,0,0.05)", "t from RECONTHETABL_ESTNEU Distribution", "t;Count", false, "3_t/3_tBlock.pdf"},
        {events5, "tReD>>hTBlockDiff(100,-0.05,0.05)", "Difference between t from CETHETA_NEUFRAME and t from RECONTHETABL_ESTNEU", "t Difference;Count", true, "3_t/5_tBlockDiff.pdf"},
        {events,  "apCE>>hAPCE(100,-1,0)",     "Reconstructed CE protons, AP from CETHETA_NEUFRAME Distribution", "AP;Count", false, "4_ap/1_apCE.pdf"},
        {events2, "apRe>>hAPRe(100,-1,0)",     "AP from RECONTHETA_ESTNEU Distribution", "AP;Count", false, "4_ap/2_apRe.pdf"},
        {events4, "apReD>>hAPDiff(100,-0.15,0.1)",    "Difference between AP from CETHETA_NEUFRAME and AP from RECONTHETA_ESTNEU", "AP Difference;Count", true, "4_ap/4_apDiff.pdf"},
        {events3, "apReD>>hAPBlock(100,-1,0)",    "AP from RECONTHETABL_ESTNEU Distribution", "AP;Count", false, "4_ap/3_apBlock.pdf"},
        {events5, "apReD>>hAPBlockDiff(100,-0.15,0.1)",    "Difference between AP from CETHETA_NEUFRAME and AP from RECONTHETABL_ESTNEU", "AP Difference;Count", true, "4_ap/5_apBlockDiff.pdf"},
        // recon minus recon block (events3 has the recon tree events2 as friend, same entries)
        {events3, Form("thetaReD:CalcValsRe%s.thetaRe>>hThetaReBlockDiff(100,0,%g,100,0,%g)", suffix.c_str(), ThetaAxisMax(), ThetaAxisMax()), "Reconstructed CE protons, RECONTHETABL_ESTNEU vs RECONTHETA_ESTNEU", "RECONTHETA_ESTNEU (degrees);RECONTHETABL_ESTNEU (degrees)", false, "1_theta/7_thetaReBlockDiff.pdf", "colz"},
        {events3, wrap360(reDiff("phiRe", "phiReD")) + ">>hPhiReBlockDiff(72,-180,180)", "Reconstructed CE protons, difference between RECONPHI_ESTNEU and RECONPHIBL_ESTNEU", diff, true, "2_phi/8_phiReBlockDiff.pdf"},
        {events3, reDiff("tRe", "tReD") + ">>hTReBlockDiff(100,-0.05,0.05)", "Difference between t from RECONTHETA_ESTNEU and t from RECONTHETABL_ESTNEU", "t Difference;Count", true, "3_t/8_tReBlockDiff.pdf"},
        {events3, reDiff("apRe", "apReD") + ">>hAPReBlockDiff(100,-0.15,0.1)", "Difference between AP from RECONTHETA_ESTNEU and AP from RECONTHETABL_ESTNEU", "AP Difference;Count", true, "4_ap/8_apReBlockDiff.pdf"},
        // percent error = 100 * (recon - CE) / CE, CE value from the friend CE tree; clipped to +-100%
        {events4, pctErr("thetaReD", "thetaCE", "hPctTheta"), "Percent Error of RECONTHETA_ESTNEU vs CETHETA_NEUFRAME", pct, true, "1_theta/6_thetaPct.pdf"},
        {events5, pctErr("thetaReD", "thetaCE", "hPctThetaB"), "Percent Error of RECONTHETABL_ESTNEU vs CETHETA_NEUFRAME", pct, true, "1_theta/6_thetaBlockPct.pdf"},
        {events4, pctErr("phiReD", "phiCE", "hPctPhi"), "Percent Error of RECONPHI_ESTNEU vs CEPHI_NEUFRAME", pct, true, "2_phi/6_phiPct.pdf"},
        {events5, pctErr("phiReD", "phiCE", "hPctPhiB"), "Percent Error of RECONPHIBL_ESTNEU vs CEPHI_NEUFRAME", pct, true, "2_phi/7_phiBlockPct.pdf"},
        {events4, pctErr("tReD", "tCE", "hPctT"), "Percent Error of t from RECONTHETA_ESTNEU vs CE t from CETHETA_NEUFRAME", pct, true, "3_t/6_tPct.pdf"},
        {events5, pctErr("tReD", "tCE", "hPctTB"), "Percent Error of t from RECONTHETABL_ESTNEU vs CE t from CETHETA_NEUFRAME", pct, true, "3_t/7_tBlockPct.pdf"},
        {events4, pctErr("apReD", "apCE", "hPctAP"), "Percent Error of AP from RECONTHETA_ESTNEU vs CE AP from CETHETA_NEUFRAME", pct, true, "4_ap/6_apPct.pdf"},
        {events5, pctErr("apReD", "apCE", "hPctAPB"), "Percent Error of AP from RECONTHETABL_ESTNEU vs CE AP from CETHETA_NEUFRAME", pct, true, "4_ap/7_apBlockPct.pdf"},
        // angle (degrees) vs AP; x is AP, y is the angle
        {events,  Form("apCE:thetaCE>>hThetaAPCE(100,0,%g,100,-1,0)", ThetaAxisMax()),    "Reconstructed CE protons, CETHETA_NEUFRAME vs AP",          std::string("CETHETA_NEUFRAME (degrees);AP"), false, "5_angle_vs_ap/1_thetaVsAP/1_thetaVsApCE.pdf", "colz"},
        {events2, Form("apRe:thetaRe>>hThetaAPRe(100,0,%g,100,-1,0)", ThetaAxisMax()),    "RECONTHETA_ESTNEU vs AP",       std::string("RECONTHETA_ESTNEU (degrees);AP"), false, "5_angle_vs_ap/1_thetaVsAP/2_thetaVsApRe.pdf", "colz"},
        {events3, Form("apReD:thetaReD>>hThetaAPBlock(100,0,%g,100,-1,0)", ThetaAxisMax()), "RECONTHETABL_ESTNEU vs AP", std::string("RECONTHETABL_ESTNEU (degrees);AP"), false, "5_angle_vs_ap/1_thetaVsAP/3_thetaVsApBlock.pdf", "colz"},
        {events,  "apCE:phiCE>>hPhiAPCE(100,0,360,100,-1,0)",        "Reconstructed CE protons, CEPHI_NEUFRAME vs AP",            std::string("CEPHI_NEUFRAME (degrees);AP"), false, "5_angle_vs_ap/2_phiVsAP/1_phiVsApCE.pdf", "colz"},
        {events2, "apRe:phiRe>>hPhiAPRe(100,0,360,100,-1,0)",        "RECONPHI_ESTNEU vs AP",         std::string("RECONPHI_ESTNEU (degrees);AP"), false, "5_angle_vs_ap/2_phiVsAP/2_phiVsApRe.pdf", "colz"},
        {events3, "apReD:phiReD>>hPhiAPBlock(100,0,360,100,-1,0)",   "RECONPHIBL_ESTNEU vs AP",   std::string("RECONPHIBL_ESTNEU (degrees);AP"), false, "5_angle_vs_ap/2_phiVsAP/3_phiVsApBlock.pdf", "colz"}};

    // AP difference (recon - CE) and its percent error against the CE angle, for the true position recon (events4) and the block recon (events5)
    {
        struct DiffSet {TTree* tree; std::string ap, name, what;};
        const std::vector<DiffSet> dsets = {{events4, "apReD", "Re", "AP from RECONTHETA_ESTNEU"}, {events5, "apReD", "Block", "AP from RECONTHETABL_ESTNEU"}};
        struct AngleAxis {std::string ce, dir, name, label, xrange;};
        const std::vector<AngleAxis> axesSets = {{"thetaCE", "1_thetaVsAP", "theta", "CETHETA_NEUFRAME (degrees)", Form("0,%g", ThetaAxisMax())}, {"phiCE", "2_phiVsAP", "phi", "CEPHI_NEUFRAME (degrees)", "0,360"}};
        for (const auto& ax : axesSets) {
            // after the three 2D plots already in the folder: AP differences first (4-5), then percent errors (6-7)
            for (size_t i = 0; i < dsets.size(); i++) {
                const auto& d = dsets[i];
                const std::string x = "CE" + suffix + "." + ax.ce;
                const std::string apCE = "CE" + suffix + ".apCE";
                const std::string base = "5_angle_vs_ap/" + ax.dir + "/";
                specsCE.push_back({d.tree, d.ap + ":" + x + ">>h" + ax.name + "ApDiff" + d.name + "(100," + ax.xrange + ",100,-0.15,0.1)",
                                   "AP difference (" + d.what + " - CE AP) vs " + (ax.name == "theta" ? "CETHETA_NEUFRAME" : "CEPHI_NEUFRAME"), ax.label + ";AP Difference", false,
                                   base + std::to_string(4 + i) + "_" + ax.name + "ApDiff" + d.name + ".pdf", "colz"});
                specsCE.push_back({d.tree, "100*(" + d.ap + "/" + apCE + "):" + x + ">>h" + ax.name + "ApPct" + d.name + "(100," + ax.xrange + ",100,-100,100)",
                                   "AP percent error (" + d.what + " vs CE AP) vs " + (ax.name == "theta" ? "CETHETA_NEUFRAME" : "CEPHI_NEUFRAME"), ax.label + ";AP Percent Error (%)", false,
                                   base + std::to_string(6 + i) + "_" + ax.name + "ApPct" + d.name + ".pdf", "colz"});
            }
        }
    }
    // FOM from the true position recon (events2, "Re") and from the block position recon (events3, "ReD")
    const std::string fomDir = "plots/5_fomCalc/";

    //All protons and CE protons, theta/phi and block theta/phi
    TFile* FileAll = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree *all = (TTree*)FileAll->Get(kAngleTree); // every fitted proton
    if (!all){Error("PlotAngle", "AnglesRaw missing in cuts/cut_angles_all.root (rerun ReconAngleAll.cpp)"); return;}

    printf("[PlotAngle] FOM setup %s\n", TDatime().AsSQLString()); fflush(stdout);
    // up to the largest RECONTHETA_ESTNEU of any reconstructed proton (true or block), so every one lands in a bin
    ROOT::RDataFrame dfAll(*all);
    auto maxT = dfAll.Max<double>("thetaRelN"), maxDT = dfAll.Max<double>("thetaRelND");
    const double maxTheta = std::max(*maxT, *maxDT) * 180.0 / TMath::Pi() + 1e-6;
    const std::vector<double> thetaEdges = ThetaEdges(maxTheta, 0.001);

    printf("[PlotAngle] %zu theta bins, writing FOM csv %s\n", thetaEdges.size() - 1, TDatime().AsSQLString()); fflush(stdout);
    // FOM tables: reconstructed protons from all events (N) and the CE protons with AP (N_CE) per theta bin
    const std::vector<KeCut> keFom = {keCE.front()}; // the FOM tables have no KE and no residual cuts
    WriteFomCsv(fomDir + "1_reconTheta_estneu", keFom, events2, "Re", all, suffix, "thetaRelN*180/TMath::Pi()", thetaEdges);
    WriteFomCsv(fomDir + "2_reconBlockTheta_estneu", keFom, events3, "ReD", all, suffix, "thetaRelND*180/TMath::Pi()", thetaEdges);

    printf("[PlotAngle] FOM csv done %s\n", TDatime().AsSQLString()); fflush(stdout);
    // every recon proton, no theta acceptance cut, so the theta plots run up to the largest RECONTHETA_ESTNEU
    const std::string deg = "*180/TMath::Pi()";

    // 3_ReconCalcVals only holds detector energy plots, so the true KE pass stops after its FOM tables
    if (energy != "fDEnergy") return;
    // theta cuts only on the phi plots: <ke>/2_phi/<th>/<file>
    // no KE cut on these plots (empty dir): the tracks without a ZDC hit (KE NAN) are in them too
    DrawGroup(canvas, foldermain.substr(0, foldermain.size() - 1), {KeCut{"", "", "no KE cut", 0}}, specsCE, "CE recon protons", "", "",
              [&](TTree* t){return std::string(t == events2 ? "thetaRe" : "CalcValsRe" + suffix + ".thetaRe");},
              [](const std::string& f, const std::string& th){return f.rfind("2_phi/", 0) == 0 ? "2_phi/" + th + "/" + f.substr(6) : std::string();});
    // phi fits of the CE protons (CE angles) and of the same events at the recon angles
    {
        const std::string phiDir = foldermain; // these fits have no KE cut either
        gSystem->mkdir(phiDir.c_str(), true);
        // CE protons: every carbon-block CE of the sim (ChargeExchange_values, lab phi of the proton at the reaction), not only the ones with a recon track
        // every carbon-block CE of the sim: cuts/cut_CE.root (CarbonCeCut, CutOpen.h), so every row passes and the cut is empty
        TTree* rawCE = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
        gROOT->cd();
        const std::string ceCut = "";
        // A_eff of each fit is built inside PhiCosFit from the own theta of each proton of that fit (no theta bins): the CETHETA_NEUFRAME (LastTheta, rad) for the CEPHI_NEUFRAME, the RECONTHETA_ESTNEU for the RECONPHI_ESTNEU. Only 1_phiCE_cosfit passes "LastAP": it takes the mean of the A the model weighted each CE with, as a closure check
        PhiCosFit(canvas, rawCE, "CEPhiNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "CEThetaNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", ceCut, "hPhiCEfit", "All CE protons (every carbon-block CE), CEPHI_NEUFRAME Distribution, no cuts", phiDir + "1_phiCE_cosfit.pdf", {}, "LastAP");
        // the same events at the recon angles (CE protons with a recon track)
        PhiCosFit(canvas, events, "phiRe", "thetaRe", "", "hPhiRefit", "Reconstructed CE protons, RECONPHI_ESTNEU Distribution, no cuts", phiDir + "2_phiRe_cosfit.pdf", {}, "", {"thetaReMin", "thetaReMax", "phiReHalf"});
        // the same fit at the block position RECONPHI_ESTNEU (phiReD) and block theta (thetaReD)
        PhiCosFit(canvas, events3, "phiReD", "thetaReD", Form("phiReD >= %g && phiReD <= %g", kBlockPhiLo, kBlockPhiHi), "hPhiReDfit", std::string("Reconstructed CE protons, RECONPHIBL_ESTNEU Distribution, fit phi ") + "0.5 to 359.5 deg (spike at 0 left out), no other cuts", phiDir + "5_phiReD_cosfit.pdf", {}, "", {"thetaReDMin", "thetaReDMax", "phiReDHalf"});
        // the same three fits with the angles about the lab z axis, no neutron frame and no neutron estimate (CE: LastPhiLab, LastTheta; recon: outPhi/outTheta, outDPhi/outDTheta of ReconAngleAll)
        PhiCosFit(canvas, rawCE, "std::fmod(LastPhiLab * 180 / TMath::Pi() + 360.0, 360.0)", "LastTheta * 180 / TMath::Pi()", ceCut, "hPhiCELabfit", "All CE protons (every carbon-block CE), CE phi about the lab z axis (no neutron frame), no cuts", phiDir + "1_phiCELab_cosfit.pdf");
        PhiCosFit(canvas, events2, "phiReLab", "thetaReLab", "", "hPhiReLabfit", "Reconstructed CE protons, recon phi about the lab z axis (no neutron estimate), no cuts", phiDir + "2_phiReLab_cosfit.pdf");
        PhiCosFit(canvas, events3, "phiReDLab", "thetaReDLab", Form("phiReDLab >= %g && phiReDLab <= %g", kBlockPhiLo, kBlockPhiHi), "hPhiReDLabfit", std::string("Reconstructed CE protons, block recon phi about the lab z axis (no neutron estimate), fit phi ") + "0.5 to 359.5 deg (spike at 0 left out), no other cuts", phiDir + "5_phiReDLab_cosfit.pdf");
    }
    printf("[PlotAngle] CE plots done %s\n", TDatime().AsSQLString()); fflush(stdout);
}
