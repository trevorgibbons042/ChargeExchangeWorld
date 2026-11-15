#include <map>
#include <ROOT/RDataFrame.hxx>
#include <set>
#include <TTree.h>

#include <TSystem.h>
#include <TStyle.h>
#include <TColor.h>
#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TMath.h>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

#include "FitCore.h"
#include <TROOT.h>
#include "CutFiles.h"
#include "TreeColumns.h"


// One output file: tree, draw expression, title text, axis titles, whether to fit the core, file name.

struct PlotSpec {
    TTree* tree;
    std::string expr, title, axes;
    bool fit;
    std::string file;
    std::string opt = ""; // draw option, "colz" for the 2D plots
};


// CE protons are red: a white-to-red palette for the 2D plots (the other plots keep viridis)
static void RedPalette(){
    double stops[2] = {0.0, 1.0}, red[2] = {1.0, 0.8}, green[2] = {1.0, 0.0}, blue[2] = {1.0, 0.0};
    TColor::CreateGradientColorTable(2, stops, red, green, blue, 255);
    gStyle->SetNumberContours(255);
}
static void ViridisPalette(){gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255);}

// Draws every plot of specs, no cuts, in <groupDir>/<spec file> (a spec file may carry its own subfolder, 1_theta/1_thetaCE.pdf).
// One Columns call (TreeColumns.h) per tree reads the draw expressions of all the specs of that tree; the histograms are filled from those columns in entry order.
static void DrawGroup(TCanvas* canvas, const std::string& groupDir, const std::vector<PlotSpec>& specs, const std::string& label)
{
    std::map<TTree*, std::vector<std::string>> exprsOf; // the columns to read per tree
    auto slot = [&](TTree* k, const std::string& e){
        auto& v = exprsOf[k];
        const size_t at = std::find(v.begin(), v.end(), e) - v.begin();
        if (at == v.size()) v.push_back(e);
        return (int)at;
    };
    std::vector<DrawSpec> draws;
    std::vector<std::vector<int>> specCol(specs.size());
    for (size_t i = 0; i < specs.size(); i++) {
        draws.push_back(SplitDraw(specs[i].expr));
        for (const std::string& e : draws[i].exprs) specCol[i].push_back(slot(specs[i].tree, e));
    }
    std::map<TTree*, std::vector<std::vector<double>>> cols; // one vector per expression, one entry per row
    for (const auto& [key, exprs] : exprsOf) cols[key] = Columns(key, exprs, "");

    std::set<std::string> fitFilesStarted; // fit_values.txt files already truncated in this call
    for (size_t i = 0; i < specs.size(); i++) {
        const auto& s = specs[i];
        const std::string path = groupDir + "/" + s.file;
        gSystem->mkdir(gSystem->GetDirName(path.c_str()), true);
        const std::string title = s.title + ", " + label + ", no cuts";
        const auto& c = cols.at(s.tree);
        canvas->Clear();
        TH1* h = MakeHist(draws[i].target); // after the Clear: a histogram of the same name from the last plot leaves the pad first
        if (!h) continue;
        const std::vector<double>& x = c[specCol[i].back()];
        const std::vector<double>* y = specCol[i].size() > 1 ? &c[specCol[i].front()] : nullptr; // "y:x": y is the first expression
        std::vector<double> xFit; // the values the unbinned fit (FitCore) takes: every row, no bins
        for (size_t r = 0; r < x.size(); r++) {
            if (y) ((TH2*)h)->Fill(x[r], (*y)[r]);
            else {h->Fill(x[r]); if (s.fit) xFit.push_back(x[r]);}
        }
        const bool is1D = h->GetDimension() == 1;
        if (is1D) h->SetLineColor(CbCE()); else RedPalette(); // CE protons are red
        h->Draw(s.opt.c_str()); // the histogram is the pad's current object, as after TTree::Draw
        h->SetTitle(PlotTitle(s.file, title + ";" + s.axes).c_str());
        if (s.fit && h->GetEntries() >= 10) FitCore(h, xFit);
        if (is1D) for (TObject* o : *h->GetListOfFunctions()) if (TF1* f = dynamic_cast<TF1*>(o)) f->SetLineColor(kBlack); // a red histogram would hide a red fit
        if (is1D && s.expr.find("phi") != std::string::npos){h->SetMinimum(0); h->SetMaximum(1.1 * h->GetMaximum());} // every 1D phi plot, 0 to max
        if (is1D && h->GetXaxis()->GetXmin() == -180 && h->GetXaxis()->GetXmax() == 180) h->GetXaxis()->SetNdivisions(508, kFALSE); // phi differences, -180 to 180 every 45 deg
        DrawCEPercent((Long64_t)x.size(), -1, s.tree->GetEntries());
        WritePlotNumbers(h, path, fitFilesStarted);
        canvas->Update();
        canvas->SaveAs(path.c_str());
        if (!is1D) ViridisPalette();
    }
}

void F3_ReconCalcVals(std::string energy = "fDEnergy"){ // energy: only "fDEnergy" (detector energy), the one CalcAPandT.cpp writes (the _DE trees of cut_APandT.root)
    FastJit();
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    if (energy != "fDEnergy"){Error("F3_ReconCalcVals", "energy must be fDEnergy: cut_APandT.root has only the detector KE trees (_DE)"); return;}
    const std::string suffix = "_DE";
    gStyle->SetOptFit(11); // fit parameters with errors (the fits are unbinned: no chi2 or probability)
    std::string foldermain = "plots/3_ReconCalcVals/1_noCuts/"; // no cuts: the same folder name as plots/6_ReconNew_1
    gSystem->mkdir(foldermain.c_str(), true);


    TTree *events = OpenCut("cuts/cut_APandT.root", "CE" + suffix);
    TTree *events2 = OpenCut("cuts/cut_APandT.root", "CalcValsRe" + suffix);
    TTree *events3 = OpenCut("cuts/cut_APandT.root", "CalcValsReD" + suffix);
    TTree *events4 = OpenCut("cuts/cut_APandT.root", "CalcValsReDIFF" + suffix);
    TTree *events5 = OpenCut("cuts/cut_APandT.root", "CalcValsReDDIFF" + suffix);
    if (!events || !events2 || !events3 || !events4 || !events5) return; // OpenCut printed the missing tree (run CalcAPandT.cpp)
    // All five trees are filled in the same loop (CalcAPandT.cpp), so their entries line up. The friends give a tree the columns of another that its plots use (all plots are no cuts).
    for (TTree* t : {events, events3, events4}) t->AddFriend(events2);
    // the diff trees also need the CE value (diff = recon - CE) for the percent error plots
    for (TTree* t : {events4, events5}) t->AddFriend(events);
    events3->AddFriend(events); // the block recon vs CETHETA_NEUFRAME 2D plot (5_thetaBlockDiff.pdf) needs the CE value next to the block value

    TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
    gPad->SetRightMargin(0.30);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); // room for the rotated Y title and bottom tick labels
    gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X"); // stats/fit box to the right, fewer x labels

    // Same x range on every theta distribution so plots compare directly: 100 bins, 0 to ThetaAxisMax deg (2.3 deg at 25 GeV/c, scaled with 1 / p, FitCore.h)
    auto thetaHist = [](const std::string& name){return std::string(Form(">>%s(100,0,%g)", name.c_str(), ThetaAxisMax()));};
    const std::string diff = "Angle Difference (degrees);Count";
    const std::string pct = "Percent Error (%);Count";
    auto pctErr = [&](const std::string& diffBranch, const std::string& ceBranch, const std::string& hist){
        // 100 * (a / b), bracketed: TTree::Draw divides first, and a * 100 / b would put values at exactly -100% a rounding step off
        return "100*(" + diffBranch + "/CE" + suffix + "." + ceBranch + ")>>" + hist + "(100,-100,100)";
    };

    //CE protons, cut_APandT.root. Theta cuts use the reconstructed (tracker) theta, in degrees.
    auto reDiff = [&](const std::string& re, const std::string& block){return "(CalcValsRe" + suffix + "." + re + "-" + block + ")";}; // recon minus recon block
    auto wrap360 = [](const std::string& d){return "(" + d + "-360*(" + d + ">180)+360*(" + d + "<-180))";}; // phi difference into -180..180
    std::vector<PlotSpec> specsCE = {
        {events,  "thetaCE" + thetaHist("hThetaCE"),  "Reconstructed CE protons, CETHETA_NEUFRAME Distribution", "CETHETA_NEUFRAME (degrees);Count", false, "1_theta/1_thetaCE.pdf"},
        {events2, "thetaRe" + thetaHist("hThetaRe"),  "Reconstructed CE protons, RECONTHETA_ESTNEU Distribution", "RECONTHETA_ESTNEU (degrees);Count", false, "1_theta/2_thetaRe.pdf"},
        // theta matrices, x is the truth (CETHETA_NEUFRAME) and y the reconstructed value; the file names keep the old Diff names
        {events4, Form("CalcValsRe%s.thetaRe:CE%s.thetaCE>>hThetaDiff(100,0,%g,100,0,%g)", suffix.c_str(), suffix.c_str(), ThetaAxisMax(), ThetaAxisMax()), "Reconstructed CE protons, RECONTHETA_ESTNEU vs CETHETA_NEUFRAME", "CETHETA_NEUFRAME (degrees);RECONTHETA_ESTNEU (degrees)", false, "1_theta/4_thetaDiff.pdf", "colz"},
        {events3, "thetaReD" + thetaHist("hThetaBlock"), "Reconstructed CE protons, RECONTHETABL_ESTNEU Distribution", "RECONTHETABL_ESTNEU (degrees);Count", false, "1_theta/3_thetaBlock.pdf"},
        {events3, Form("thetaReD:CE%s.thetaCE>>hThetaBlockDiff(100,0,%g,100,0,%g)", suffix.c_str(), ThetaAxisMax(), ThetaAxisMax()), "Reconstructed CE protons, RECONTHETABL_ESTNEU vs CETHETA_NEUFRAME", "CETHETA_NEUFRAME (degrees);RECONTHETABL_ESTNEU (degrees)", false, "1_theta/5_thetaBlockDiff.pdf", "colz"},
        {events,  "phiCE>>hPhiCE(100,0,360)",    "Reconstructed CE protons, CEPHI_NEUFRAME Distribution", "CEPHI_NEUFRAME (degrees);Count", false, "2_phi/1_phiCE.pdf"},
        {events2, "phiRe>>hPhiRe(100,0,360)",    "Reconstructed CE protons, RECONPHI_ESTNEU Distribution", "RECONPHI_ESTNEU (degrees);Count", false, "2_phi/2_phiRe.pdf"},
        {events4, "phiReD>>hPhiDiff(100,-180,180)",   "Reconstructed CE protons, difference between CEPHI_NEUFRAME and RECONPHI_ESTNEU", diff, true, "2_phi/4_phiDiff.pdf"},
        {events3, "phiReD>>hPhiBlock(100,0,360)",   "Reconstructed CE protons, RECONPHIBL_ESTNEU Distribution", "RECONPHIBL_ESTNEU (degrees);Count", false, "2_phi/3_phiBlock.pdf"},
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
        struct DiffSet {TTree* tree; std::string name, what;};
        const std::vector<DiffSet> dsets = {{events4, "Re", "AP from RECONTHETA_ESTNEU"}, {events5, "Block", "AP from RECONTHETABL_ESTNEU"}};
        struct AngleAxis {std::string ce, dir, name, label, xrange;};
        const std::vector<AngleAxis> axesSets = {{"thetaCE", "1_thetaVsAP", "theta", "CETHETA_NEUFRAME (degrees)", Form("0,%g", ThetaAxisMax())}, {"phiCE", "2_phiVsAP", "phi", "CEPHI_NEUFRAME (degrees)", "0,360"}};
        for (const auto& ax : axesSets) {
            // after the three 2D plots already in the folder: AP differences first (4-5), then percent errors (6-7)
            for (size_t i = 0; i < dsets.size(); i++) {
                const auto& d = dsets[i];
                const std::string x = "CE" + suffix + "." + ax.ce;
                const std::string apCE = "CE" + suffix + ".apCE";
                const std::string what = ax.name == "theta" ? "CETHETA_NEUFRAME" : "CEPHI_NEUFRAME";
                const std::string base = "5_angle_vs_ap/" + ax.dir + "/";
                specsCE.push_back({d.tree, "apReD:" + x + ">>h" + ax.name + "ApDiff" + d.name + "(100," + ax.xrange + ",100,-0.15,0.1)",
                                   "AP difference (" + d.what + " - CE AP) vs " + what, ax.label + ";AP Difference", false,
                                   base + std::to_string(4 + i) + "_" + ax.name + "ApDiff" + d.name + ".pdf", "colz"});
                specsCE.push_back({d.tree, "100*(apReD/" + apCE + "):" + x + ">>h" + ax.name + "ApPct" + d.name + "(100," + ax.xrange + ",100,-100,100)",
                                   "AP percent error (" + d.what + " vs CE AP) vs " + what, ax.label + ";AP Percent Error (%)", false,
                                   base + std::to_string(6 + i) + "_" + ax.name + "ApPct" + d.name + ".pdf", "colz"});
            }
        }
    }

    DrawGroup(canvas, foldermain.substr(0, foldermain.size() - 1), specsCE, "CE recon protons");

    // phi fits of the CE protons (CE angles) and of the same events at the recon angles
    {
        const std::string phiDir = foldermain; // these fits have no cuts either
        // CE protons: every carbon-block CE of the sim (ChargeExchange_values, lab phi of the proton at the reaction), not only the ones with a recon track
        // cuts/cut_CE.root holds only the carbon-block CE (CarbonCeCut, FitCore.h), so no cut is needed here
        TTree* rawCE = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
        gROOT->cd();
        // A_eff of each fit is built inside PhiCosFit from the own theta of each proton of that fit (no theta bins): the CETHETA_NEUFRAME (LastTheta, rad) for the CEPHI_NEUFRAME, the RECONTHETA_ESTNEU for the RECONPHI_ESTNEU. Only 1_phiCE_cosfit passes "LastAP": it takes the mean of the A the model weighted each CE with, as a closure check
        PhiCosFit(canvas, rawCE, "CEPhiNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "CEThetaNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "", "hPhiCEfit", "All CE protons (every carbon-block CE), CEPHI_NEUFRAME Distribution, no cuts", phiDir + "1_phiCE_cosfit.pdf", {}, "LastAP");
        // the same events at the recon angles (CE protons with a recon track)
        PhiCosFit(canvas, events, "phiRe", "thetaRe", "", "hPhiRefit", "Reconstructed CE protons, RECONPHI_ESTNEU Distribution, no cuts", phiDir + "2_phiRe_cosfit.pdf", {}, "", {"thetaReMin", "thetaReMax", "phiReHalf"});
        // the same fit at the block position RECONPHI_ESTNEU (phiReD) and block theta (thetaReD)
        PhiCosFit(canvas, events3, "phiReD", "thetaReD", Form("phiReD > %g && phiReD < %g", kBlockPhiLo, kBlockPhiHi), "hPhiReDfit", std::string("Reconstructed CE protons, RECONPHIBL_ESTNEU Distribution, fit phi ") + "above 0 and below 360 deg (spike at 0 left out), no other cuts", phiDir + "5_phiReD_cosfit.pdf", {}, "", {"thetaReDMin", "thetaReDMax", "phiReDHalf"}, {"thetaReDBMin", "thetaReDBMax", "phiReDBHalf"});
        // the same three fits with the angles about the lab z axis, no neutron frame and no neutron estimate (CE: LastPhiLab, LastTheta; recon: outPhi/outTheta, outDPhi/outDTheta of ReconAngleAll)
        PhiCosFit(canvas, rawCE, "std::fmod(LastPhiLab * 180 / TMath::Pi() + 360.0, 360.0)", "LastTheta * 180 / TMath::Pi()", "", "hPhiCELabfit", "All CE protons (every carbon-block CE), CE phi about the lab z axis (no neutron frame), no cuts", phiDir + "1_phiCELab_cosfit.pdf");
        PhiCosFit(canvas, events2, "phiReLab", "thetaReLab", "", "hPhiReLabfit", "Reconstructed CE protons, recon phi about the lab z axis (no neutron estimate), no cuts", phiDir + "2_phiReLab_cosfit.pdf");
        PhiCosFit(canvas, events3, "phiReDLab", "thetaReDLab", Form("phiReDLab > %g && phiReDLab < %g", kBlockPhiLo, kBlockPhiHi), "hPhiReDLabfit", std::string("Reconstructed CE protons, block recon phi about the lab z axis (no neutron estimate), fit phi ") + "above 0 and below 360 deg (spike at 0 left out), no other cuts", phiDir + "5_phiReDLab_cosfit.pdf", {}, "", {}, {"thetaReDLabMin", "thetaReDLabMax", "phiReDLabHalf"});
    }
    printf("[F3_ReconCalcVals] CE plots done %s\n", TDatime().AsSQLString()); fflush(stdout);
}
