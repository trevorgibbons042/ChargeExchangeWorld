#include <set>
#include <TFile.h>
#include <TTree.h>

#include <TSystem.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TGraph.h>
#include <TMath.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "FitCore.h"


// One output file: tree, draw expression, title text, axis titles, whether to fit the core, file name.
struct PlotSpec {
    TTree* tree;
    std::string expr, title, axes;
    bool fit;
    std::string file;
    std::string opt = ""; // TTree::Draw option, "colz" for the 2D plots
    std::string cut = ""; // extra cut for just this plot (AP slices)
};

struct KECut {std::string dir, cut, label;};
struct ThetaCut {std::string dir, op, label;}; // op is appended to the theta expression

// Draws every plot in specs once per (KE cut, theta cut). Plots go in <groupDir>/<ke.dir>/<th.dir>/<spec file>, where a spec file may carry its own subfolder (1_theta/1_thetaCE.pdf). thetaExpr is the recon theta in degrees for the tree being drawn.
static void DrawGroup(TCanvas* canvas, const std::string& groupDir, const std::vector<KECut>& kes,
                      const std::vector<PlotSpec>& specs, const std::string& label, const std::string& baseCut,
                      const std::string& extraCut, const std::function<std::string(TTree*)>& thetaExpr)
{
    const std::vector<ThetaCut> thetas = {
        {"th_all", "", ""},
        {"th_gt0p1", " > 0.1", "theta > 0.1 deg"},
        {"th_gt0p3", " > 0.3", "theta > 0.3 deg"}};

    std::set<std::string> fitFilesStarted; // fit_values.txt files already truncated in this call
    for (const auto& ke : kes) {
        for (const auto& th : thetas) {
            for (const auto& s : specs) {
                const std::string path = groupDir + "/" + ke.dir + "/" + th.dir + "/" + s.file;
                gSystem->mkdir(gSystem->GetDirName(path.c_str()), true);

                std::vector<std::string> parts;
                for (const std::string& c : {baseCut, ke.cut, extraCut, s.cut}) {
                    if (!c.empty()) parts.push_back("(" + c + ")");
                }
                if (!th.op.empty()) parts.push_back("(" + thetaExpr(s.tree) + th.op + ")");
                std::string cut;
                for (size_t i = 0; i < parts.size(); i++) cut += (i ? " && " : "") + parts[i];

                std::string title = s.title + ", " + label;
                if (!ke.label.empty()) title += ", " + ke.label;
                if (!th.label.empty()) title += ", " + th.label;

                canvas->Clear();
                s.tree->Draw(s.expr.c_str(), cut.c_str(), s.opt.c_str());
                TH1* h = s.tree->GetHistogram();
                if (!h) continue;
                h->SetTitle(PlotTitle(s.file, title + ";" + s.axes).c_str());
                if (s.fit && h->GetEntries() >= 10) FitCore(h);
                WritePlotNumbers(h, path, fitFilesStarted);
                canvas->Update();
                canvas->SaveAs(path.c_str());
            }
        }
    }
}

const double totalNeutrons = 100e6; // neutrons simulated, the FOM denominator
const double beamPol = 0.7;         // beam polarization P in dN/dphi ~ 1 + P*AP*cos(phi)

// Pulls (AP, phi) of every event passing cut out of tree t.
static void GetApPhi(TTree* t, const std::string& sfx, const std::string& cut, std::vector<double>& ap, std::vector<double>& phi){
    t->SetEstimate(t->GetEntries() + 1);
    Long64_t n = t->Draw(("ap" + sfx + ":phi" + sfx).c_str(), cut.c_str(), "goff");
    ap.assign(t->GetV1(), t->GetV1() + n);
    phi.assign(t->GetV2(), t->GetV2() + n);
}

// For every KE folder and theta bin: the recon phi distribution (position set by sfx: "Re" true, "ReD" block) fitted with N*(1 + P*<AP>*cos(phi)), P = beamPol
// and <AP> = mean AP of the events in the bin, both fixed, so only the normalization N floats.
static void DrawPhiBins(TCanvas* canvas, const std::string& groupDir, const std::vector<KECut>& kes, TTree* t,
                        const std::string& sfx, const std::string& what, const std::vector<double>& thetaEdges, bool showPhiBins)
{
    std::set<std::string> numberFiles; // stats and fit files already truncated in this call
    gSystem->mkdir(groupDir.c_str(), true);
    for (const auto& ke : kes) {
        std::string keCut = ke.cut.empty() ? "1" : ke.cut;
        std::vector<double> ap, phi;

        std::vector<double> fom(thetaEdges.size() - 1, 0.0);
        struct BinNumbers {double lo, hi, n, meanAP, fitN, fitNErr, chi2ndf;};
        std::vector<BinNumbers> rows; // written to numbers.txt next to the FOM plot
        for (size_t b = 0; b + 1 < thetaEdges.size(); b++) {
            double lo = thetaEdges[b], hi = thetaEdges[b + 1];
            std::string name = Form("th_%g_%g", lo, hi);
            std::replace(name.begin(), name.end(), '.', 'p');

            GetApPhi(t, sfx, Form("(%s) && theta%s >= %g && theta%s < %g", keCut.c_str(), sfx.c_str(), lo, sfx.c_str(), hi), ap, phi);
            if (ap.empty()) {rows.push_back({lo, hi, 0, 0, NAN, NAN, NAN}); continue;}
            double meanAP = 0;
            for (double x : ap) meanAP += x;
            meanAP /= ap.size();
            fom[b] = ap.size() * std::fabs(meanAP) / totalNeutrons;
            rows.push_back({lo, hi, (double)ap.size(), meanAP, NAN, NAN, NAN});

            if (!showPhiBins) continue; // phi-bin fit plots are optional, the FOM above does not need them
            gSystem->mkdir((groupDir + "/" + ke.dir + "/phi_fit").c_str(), true);

            TH1D h("hPhiBin", "", 100, 0, 360);
            for (double x : phi) h.Fill(x < 0 ? x + 360 : x);
            h.SetStats(1); // fit box (gStyle OptFit) lists the fit parameters
            h.SetTitle(PlotTitle("phi_fit.pdf", Form("Recon %s phi, CE events CE protons, %s, theta %g to %g deg, <AP> = %.4f, P = %.2f;Phi (degrees);Count",
                            what.c_str(), ke.label.c_str(), lo, hi, meanAP, beamPol)).c_str());
            TF1 f("fPhiBin", "[0]*(1 + [1]*cos(x*TMath::DegToRad()))", 0, 360);
            f.SetParNames("N", "P*<AP>");
            f.SetParameters(h.GetEntries() * h.GetBinWidth(1) / 360.0, beamPol * meanAP);
            f.FixParameter(1, beamPol * meanAP);
            h.Fit(&f, "QR");
            rows.back().fitN = f.GetParameter(0);
            rows.back().fitNErr = f.GetParError(0);
            rows.back().chi2ndf = f.GetNDF() ? f.GetChisquare() / f.GetNDF() : NAN;

            canvas->Clear();
            const int oldOptStat = gStyle->GetOptStat();
            gStyle->SetOptStat(10); // entries only, plus the fit box
            h.SetLineWidth(2);
            h.Draw("HIST");
            f.SetLineColor(kRed); f.SetLineWidth(3);
            f.Draw("SAME");
            canvas->Update();
            canvas->SaveAs((groupDir + "/" + ke.dir + "/phi_fit/" + name + ".pdf").c_str());
            WritePlotNumbers(&h, groupDir + "/" + ke.dir + "/phi_fit/" + name + ".pdf", numberFiles);
            gStyle->SetOptStat(oldOptStat);
        }

        // FOM = (#CE events in the bin * |mean AP|) / total neutrons, one bar per theta bin
        gSystem->mkdir((groupDir + "/" + ke.dir).c_str(), true);
        if (FILE* out = fopen((groupDir + "/" + ke.dir + "/numbers.txt").c_str(), "w")) {
            fprintf(out, "FOM and phi fit, Recon %s, CE events CE protons, %s, N_neutrons = %.3g, P = %.2f\n", what.c_str(), ke.label.empty() ? "no KE cut" : ke.label.c_str(), totalNeutrons, beamPol);
            fprintf(out, "%-14s %8s %10s %12s %10s %10s %14s %9s\n", "theta (deg)", "N", "<AP>", "FOM", "fit N", "fit N err", "fit P*<AP> fixed", "chi2/ndf");
            for (size_t b = 0; b < rows.size(); b++) {
                const auto& r = rows[b];
                fprintf(out, "%5g - %-6g %8.0f %10.4f %12.4e %10.4g %10.4g %14.4g %9.3g\n", r.lo, r.hi, r.n, r.meanAP, fom[b], r.fitN, r.fitNErr, beamPol * r.meanAP, r.chi2ndf);
            }
            fclose(out);
        }
        TGraph gFom;
        for (size_t b = 0; b < fom.size(); b++) gFom.SetPoint(b, 0.5 * (thetaEdges[b] + thetaEdges[b + 1]), fom[b]);
        gFom.SetMarkerStyle(20);
        gFom.SetTitle(PlotTitle("FOM_vs_theta.pdf", Form("FOM = N_{CE} #times |#LTAP#GT| / N_{neutrons}, N_{neutrons} = %.2g, CE events CE protons, %s;Recon %s theta (degrees, bin center);FOM", totalNeutrons, ke.label.c_str(), what.c_str())).c_str());
        canvas->Clear();
        gFom.Draw("APL");
        canvas->Update();
        gSystem->mkdir((groupDir + "/" + ke.dir).c_str(), true);
        canvas->SaveAs((groupDir + "/" + ke.dir + "/FOM_vs_theta.pdf").c_str());
    }
}

void PlotAngle(std::string energy = "fKEnergy", bool showPhiBins = true){ // energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the input files and plot folder; showPhiBins: also save the per-theta-bin phi fit plots
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    gStyle->SetOptFit(1111); // probability, chi2/ndf, fit parameters with errors
    const std::string energyDir = (energy == "fDEnergy") ? "2_KEDetector" : "1_KEtrue";
    std::string foldermain = "plots/3_angle_recon/";
    gSystem->mkdir(foldermain.c_str(), true);


    TFile* File = TFile::Open("cut_APandT.root", "READ");

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
    gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetNdivisions(505, "X"); // stats/fit box to the right, fewer x labels

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
    // No cut, plus a gt (>) cut at 50-90% of the 25 GeV beam energy.
    const double initialKE = 25000; // MeV
    const std::string keLabel = (energy == "fDEnergy") ? "ZDC KE" : "true KE"; // which energy the cuts are on, shown in the titles
    auto makeKECuts = [&](){
        std::vector<KECut> kes = {{"ke_all", "", "no " + keLabel + " cut"}};
        for (int pct : {50, 60, 70, 80, 90}){
            std::string cutKE = std::to_string((int)(initialKE * pct / 100.0));
            std::string pctStr = std::to_string(pct);
            kes.push_back({"ke_gt" + pctStr, "finalKE > " + cutKE, keLabel + " > " + cutKE + " (" + pctStr + "%)"});
        }
        return kes;
    };
    const std::vector<KECut> keCE = makeKECuts();
    std::vector<PlotSpec> specsCE = {
        {events,  "thetaCE" + thetaHist("hThetaCE"),  "CE theta Distribution", ang, false, "1_theta/1_thetaCE.pdf"},
        {events2, "thetaRe" + thetaHist("hThetaRe"),  "Recon theta Distribution", ang, false, "1_theta/2_thetaRe.pdf"},
        {events4, "thetaReD", "Difference between CE and Recon Theta", diff, true, "1_theta/4_thetaDiff.pdf"},
        {events3, "thetaReD" + thetaHist("hThetaBlock"), "Recon block theta Distribution", ang, false, "1_theta/3_thetaBlock.pdf"},
        {events5, "thetaReD", "Difference between CE and Recon Block Theta", diff, true, "1_theta/5_thetaBlockDiff.pdf"},
        {events,  "phiCE>>hPhiCE(100,0,360)",    "CE phi Distribution", ang, false, "2_phi/1_phiCE.pdf"},
        {events2, "phiRe>>hPhiRe(100,0,360)",    "Recon phi Distribution", ang, false, "2_phi/2_phiRe.pdf"},
        {events4, "phiReD>>hPhiDiff(100,-180,180)",   "Difference between CE and Recon Phi", diff, true, "2_phi/4_phiDiff.pdf"},
        {events3, "phiReD>>hPhiBlock(100,0,360)",   "Recon block phi Distribution", ang, false, "2_phi/3_phiBlock.pdf"},
        {events5, "phiReD>>hPhiBlockDiff(100,-180,180)",   "Difference between CE and Recon Block Phi", diff, true, "2_phi/5_phiBlockDiff.pdf"},
        {events,  "tCE>>hTCE(100,0,0.05)", "CE t Distribution", "t;Count", false, "3_t/1_tCE.pdf"},
        {events2, "tRe>>hTRe(100,0,0.05)", "Recon t Distribution", "t;Count", false, "3_t/2_tRe.pdf"},
        {events4, "tReD>>hTDiff(100,-0.05,0)", "Difference between CE and Recon t", "t Difference;Count", true, "3_t/4_tDiff.pdf"},
        {events3, "tReD>>hTBlock(100,0,0.05)", "Recon block t Distribution", "t;Count", false, "3_t/3_tBlock.pdf"},
        {events5, "tReD>>hTBlockDiff(100,-0.05,0)", "Difference between CE and Recon Block t", "t Difference;Count", true, "3_t/5_tBlockDiff.pdf"},
        {events,  "apCE>>hAPCE(100,-1,0)",     "CE AP Distribution", "AP;Count", false, "4_ap/1_apCE.pdf"},
        {events2, "apRe>>hAPRe(100,-1,0)",     "Recon AP Distribution", "AP;Count", false, "4_ap/2_apRe.pdf"},
        {events4, "apReD>>hAPDiff(100,-0.15,0.1)",    "Difference between CE and Recon AP", "AP Difference;Count", true, "4_ap/4_apDiff.pdf"},
        {events3, "apReD>>hAPBlock(100,-1,0)",    "Recon block AP Distribution", "AP;Count", false, "4_ap/3_apBlock.pdf"},
        {events5, "apReD>>hAPBlockDiff(100,-0.15,0.1)",    "Difference between CE and Recon Block AP", "AP Difference;Count", true, "4_ap/5_apBlockDiff.pdf"},
        // percent error = 100 * (recon - CE) / CE, CE value from the friend CE tree; clipped to +-100%
        {events4, pctErr("thetaReD", "thetaCE", "hPctTheta"), "Percent Error of Recon Theta vs CE", pct, true, "1_theta/6_thetaPct.pdf"},
        {events5, pctErr("thetaReD", "thetaCE", "hPctThetaB"), "Percent Error of Recon Block Theta vs CE", pct, true, "1_theta/7_thetaBlockPct.pdf"},
        {events4, pctErr("phiReD", "phiCE", "hPctPhi"), "Percent Error of Recon Phi vs CE", pct, true, "2_phi/6_phiPct.pdf"},
        {events5, pctErr("phiReD", "phiCE", "hPctPhiB"), "Percent Error of Recon Block Phi vs CE", pct, true, "2_phi/7_phiBlockPct.pdf"},
        {events4, pctErr("tReD", "tCE", "hPctT"), "Percent Error of Recon t vs CE", pct, true, "3_t/6_tPct.pdf"},
        {events5, pctErr("tReD", "tCE", "hPctTB"), "Percent Error of Recon Block t vs CE", pct, true, "3_t/7_tBlockPct.pdf"},
        {events4, pctErr("apReD", "apCE", "hPctAP"), "Percent Error of Recon AP vs CE", pct, true, "4_ap/6_apPct.pdf"},
        {events5, pctErr("apReD", "apCE", "hPctAPB"), "Percent Error of Recon Block AP vs CE", pct, true, "4_ap/7_apBlockPct.pdf"},
        // angle (degrees) vs AP; x is AP, y is the angle
        {events,  "apCE:thetaCE>>hThetaAPCE(100,0,2.3,100,-1,0)",    "CE theta vs AP",          thetaAP, false, "5_angle_vs_ap/1_thetaVsAP/1_thetaVsApCE.pdf", "colz"},
        {events2, "apRe:thetaRe>>hThetaAPRe(100,0,2.3,100,-1,0)",    "Recon theta vs AP",       thetaAP, false, "5_angle_vs_ap/1_thetaVsAP/2_thetaVsApRe.pdf", "colz"},
        {events3, "apReD:thetaReD>>hThetaAPBlock(100,0,2.3,100,-1,0)", "Recon block theta vs AP", thetaAP, false, "5_angle_vs_ap/1_thetaVsAP/3_thetaVsApBlock.pdf", "colz"},
        {events,  "apCE:phiCE>>hPhiAPCE(100,0,360,100,-1,0)",        "CE phi vs AP",            phiAP, false, "5_angle_vs_ap/2_phiVsAP/1_phiVsApCE.pdf", "colz"},
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
    std::vector<double> thetaEdges; // 0.1 deg bins, 0 to 2.3 deg (theta acceptance 0.04 rad)
    for (int i = 0; i <= 23; i++) thetaEdges.push_back(i * 0.1);

    const std::string reconDir = foldermain + "1_CEvalues/" + energyDir + "/";
    DrawGroup(canvas, reconDir + "1_CEprotons", keCE, specsCE, "CE events CE protons", "", "",
              [&](TTree* t){return std::string(t == events2 ? "thetaRe" : "CalcValsRe" + suffix + ".thetaRe");});
    // FOM from the true position recon (events2, "Re") and from the block position recon (events3, "ReD")
    DrawPhiBins(canvas, reconDir + "2_fomCalc/1_true", keCE, events2, "Re", "true position", thetaEdges, showPhiBins); // theta bin edges, degrees
    DrawPhiBins(canvas, reconDir + "2_fomCalc/2_block", keCE, events3, "ReD", "block position", thetaEdges, showPhiBins);

    //All protons and CE protons, theta/phi and block theta/phi
    TFile* FileAll = TFile::Open("cut_angles_all.root", "READ");
    TTree *all = (TTree*)FileAll->Get(("Angles" + suffix).c_str());
    TTree *allCE = (TTree*)FileAll->Get(("AnglesCE" + suffix).c_str());

    // Detector acceptance is square: above this theta (rad) phi is only accepted toward the corners.
    const double thetaMax = 0.04;
    const std::string thetaMaxCut = Form("outTheta < %g", thetaMax);
    const std::string deg = "*180/TMath::Pi()";

    const std::vector<KECut> keAll = makeKECuts();
    auto angleSpecs = [&](TTree* t){
        return std::vector<PlotSpec>{
            {t, "outTheta" + deg + thetaHist("hAngTheta"),  "Recon theta Distribution", ang, false, "1_theta.pdf"},
            {t, "outPhi" + deg + ">>hAngPhi(100,0,360)",    "Recon phi Distribution", ang, false, "2_phi.pdf"},
            {t, "outDTheta" + deg + thetaHist("hAngBlockTheta"), "Recon block theta Distribution", ang, false, "3_blocktheta.pdf"},
            {t, "outDPhi" + deg + ">>hAngBlockPhi(100,0,360)",   "Recon block phi Distribution", ang, false, "4_blockphi.pdf"}};
    };
    auto thetaDeg = [&](TTree*){return std::string("outTheta") + deg;};

    // These trees have no calc values (t, AP), only recon angles, so they sit outside 1_CEvalues
    const std::string nonCEDir = foldermain + "2_nonCEangles/" + energyDir + "/";
    DrawGroup(canvas, nonCEDir + "2_allprotons", keAll, angleSpecs(all),
              "all events all protons", "", thetaMaxCut, thetaDeg);
    //Protons in CE events, may or may not come from the CE process
    DrawGroup(canvas, nonCEDir + "1_ce_event_protons", keAll, angleSpecs(allCE),
              "CE events all protons", "", thetaMaxCut, thetaDeg);

    // KE cut efficiency on protons in CE events (same theta acceptance as the plots above). Real = from the CE
    // (isFromCE == 1), fake = in a CE event but not from the CE. A good cut keeps most real and few fake.
    const Long64_t nReal = allCE->GetEntries((thetaMaxCut + " && isFromCE == 1").c_str());
    const Long64_t nFake = allCE->GetEntries((thetaMaxCut + " && isFromCE == 0").c_str());
    auto percent = [](Long64_t n, Long64_t total){return total ? 100.0 * n / total : 0.0;};
    FILE* eff = fopen((nonCEDir + "ke_cut_efficiency.txt").c_str(), "w");
    if (eff) {
        fprintf(eff, "%s cut efficiency (finalKE), CE events all protons, %s < %g. real = %lld, fake = %lld\n",
                keLabel.c_str(), "outTheta", thetaMax, (long long)nReal, (long long)nFake);
        fprintf(eff, "%-22s %12s %14s\n", "cut", "fake kept %", "real left out %");
        for (const auto& ke : keAll) {
            if (ke.cut.empty()) continue;
            const Long64_t realKept = allCE->GetEntries((thetaMaxCut + " && isFromCE == 1 && " + ke.cut).c_str());
            const Long64_t fakeKept = allCE->GetEntries((thetaMaxCut + " && isFromCE == 0 && " + ke.cut).c_str());
            fprintf(eff, "%-22s %12.2f %14.2f\n", ke.cut.c_str(), percent(fakeKept, nFake), 100.0 - percent(realKept, nReal));
        }
        fclose(eff);
    }
}
