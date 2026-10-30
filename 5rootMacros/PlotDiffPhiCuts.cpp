#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <TTree.h>

#include <cmath>
#include <map>
#include <set>
#include <string>
#include <utility>

#include "CutCache.h"
#include "DrawPlan.h"
#include "FitCore.h"
#include "CutFiles.h"

// Protons whose recon block phi is off from the CE phi by phiLo to phiHi degrees (CalcValsReDDIFF.phiReD in cut_APandT.root).
// Their ZDC positions come from cut_ParentID.root, matched on (fEvent, ParentID); their KE is finalKE from cut_APandT.root.
// energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the input trees like PlotAngle.cpp. Next to each window plot the same for every recon
// proton, and the block fit residual, R^2 and impact of every recon proton and the CE ones.
void PlotDiffPhiCuts(double phiLo = -180, double phiHi = -50, std::string energy = "fKEnergy"){
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    const std::string folder = std::string("plots/4_CEblock/3_diffphiBlock/");
    gSystem->mkdir(folder.c_str(), true);
    gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X");

    TFile* fAP = TFile::Open("cuts/cut_APandT.root", "READ");
    TTree* diff = (TTree*)fAP->Get(("CalcValsReDDIFF" + suffix).c_str());
    // ZDC protons of the CE children, from the cut cache (CutCache.h, built from the seed cut_ParentID.root): every row passes the ZDC proton cut
    const CacheTree cachePos = GetCutTree("output0_t*.root", "ParticleInfo", {"@isParentCE", "boundaryType == 2 && pdg == 2212"});
    TTree* pos = cachePos.tree;
    if (!diff || !pos){Error("PlotDiffPhiCuts", "missing tree in cut_APandT.root or cut_ParentID.root"); return;}

    // Every column pull of this macro is booked here and collected in one multithreaded pass per tree (DrawPlan.h), the render code below reads the columns in entry order.
    ROOT::EnableImplicitMT(NThreads());
    DrawPlan plan;
    // events in the phi diff window: block residual cut as the list string (CutTrees.cpp), the list itself is not set on the tree
    auto cDiff = plan.C(diff, {"fEvent", "TrackID", "finalKE"}, Form("phiReD >= %g && phiReD < %g", phiLo, phiHi), ListCut("ke_all", true));
    // Every recon proton (Angles tree of cut_angles_all.root, all events): (fEvent, TrackID) for the positions, finalKE, isFromCE and the
    // block line fit quality (ReconAngleAll.cpp)
    TFile* fAng = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree* ang = fAng ? (TTree*)fAng->Get(("Angles" + suffix).c_str()) : nullptr;
    if (!ang || !ang->GetBranch("r2D")){Error("PlotDiffPhiCuts", "Angles tree with r2D missing in cut_angles_all.root (rerun ReconAngleAll.cpp)"); return;}
    auto cAng = plan.C(ang, {"fEvent", "TrackID", "finalKE", "isFromCE", "resDRms", "r2D", "impactD", "outPhi", "outDPhi"});
    DrawPlan::Cols cEx;
    if (ang->GetBranch("nProtonsEv")) cEx = plan.C(ang, {"resRms", "impact", "r2", "nProtonsEv", "zOrderOK", "resHalfRatioD", "outTheta*180/TMath::Pi()", "outDTheta*180/TMath::Pi()", "nHits"});
    // ZDC protons of the CE children (pos) and of every event (allEv), the cut caches: every row passes the ZDC proton cut
    const CacheTree cacheAllEv = GetCutTree("output0_t*.root", "ParticleInfo", {"boundaryType == 2 && pdg == 2212"});
    TTree* allEv = cacheAllEv.tree;
    if (!allEv){Error("PlotDiffPhiCuts", "ParticleInfo cut caches missing (run CutEvents.cpp)"); return;}
    auto cPos = plan.C(pos, {"fEvent", "ParentID", "fDx", "fDy", "fX", "fY", "TrackID"});
    auto cAllEv = plan.C(allEv, {"fEvent", "TrackID", "fDx", "fDy", "fX", "fY"});
    plan.Run();
    gROOT->cd();
    const Long64_t nSel = cDiff.N();
    std::set<std::pair<long long, long long>> keys;
    TH1D hKE("hKE", "", 100, 0, 30000);
    for (Long64_t i = 0; i < nSel; i++){
        keys.insert({llround(cDiff.V(0)[i]), llround(cDiff.V(1)[i])});
        hKE.Fill(cDiff.V(2)[i]);
    }

    const Long64_t nAng = cAng.N();
    std::vector<const double*> angCols(9); // views of the DrawPlan columns, no copy
    for (int c = 0; c < 9; c++) angCols[c] = cAng.V(c);
    // every recon proton (all events, CE ones included) in the same block phi diff window: the CE protons keep their block phi - CE phi
    // window above, the non-CE ones have no CE phi so their diff is recon block phi - recon tracker phi (outPhi), wrapped to -180 to 180
    // like CalcAPandT.cpp, with the same block residual cut
    std::set<std::pair<long long, long long>> reconKeys;
    for (Long64_t i = 0; i < nAng; i++){
        const double dPhi = std::remainder((angCols[8][i] - angCols[7][i]) * 180 / TMath::Pi(), 360.0);
        if (angCols[3][i] < 0.5 && dPhi >= phiLo && dPhi < phiHi && PassesBlockResidual(angCols[4][i])) reconKeys.insert({llround(angCols[0][i]), llround(angCols[1][i])});
    }

    // positions at the ZDC, block (fDx, fDy) and true (fX, fY): the window protons from cut_ParentID.root (cut_APandT.root TrackID is the
    // CE neutron's track, so it matches the proton's ParentID there), every recon proton from cut_AllEvProtons.root (matched on TrackID)
    TH2D hBlock("hBlock", "", 100, -200, 200, 100, -200, 200), hTrue("hTrue", "", 100, -200, 200, 100, -200, 200);
    TH2D hBlockAll("hBlockAll", "", 100, -200, 200, 100, -200, 200), hTrueAll("hTrueAll", "", 100, -200, 200, 100, -200, 200);
    const Long64_t nPos = cPos.N();
    std::set<std::pair<long long, long long>> ceWinKeys; // the window CE protons as (fEvent, proton TrackID), the Angles tree key
    for (Long64_t i = 0; i < nPos; i++){
        if (!keys.count({llround(cPos.V(0)[i]), llround(cPos.V(1)[i])})) continue;
        ceWinKeys.insert({llround(cPos.V(0)[i]), llround(cPos.V(6)[i])});
        hBlock.Fill(cPos.V(2)[i], cPos.V(3)[i]);
        hTrue.Fill(cPos.V(4)[i], cPos.V(5)[i]);
    }
    const Long64_t nAllEv = cAllEv.N();
    for (Long64_t i = 0; i < nAllEv; i++){
        if (!reconKeys.count({llround(cAllEv.V(0)[i]), llround(cAllEv.V(1)[i])})) continue;
        hBlockAll.Fill(cAllEv.V(2)[i], cAllEv.V(3)[i]);
        hTrueAll.Fill(cAllEv.V(4)[i], cAllEv.V(5)[i]);
    }
    hBlockAll.Add(&hBlock);
    hTrueAll.Add(&hTrue);

    TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
    gPad->SetRightMargin(0.30);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);

    const std::string cutLabel = Form("CE recon protons, block phi diff %g to %g deg, N = %lld", phiLo, phiHi, nSel);
    std::set<std::string> numbersStarted;
    auto save = [&](TH1* h, const std::string& title, const std::string& file, const char* opt){
        h->SetTitle(PlotTitle(file, title).c_str());
        canvas->Clear();
        h->Draw(opt);
        WritePlotNumbers(h, folder + file, numbersStarted);
        canvas->Update();
        canvas->SaveAs((folder + file).c_str());
    };
    const std::string xy = ";X Position (mm);Y Position (mm)";
    save(&hBlock, "Position at ZDC BLOCK, " + cutLabel + xy, "1_xy_block.pdf", "colz");
    const std::string reconLabel = Form("all recon protons, block phi diff %g to %g deg (CE: block - CE, non-CE: block - tracker), N = %lld", phiLo, phiHi,
                                            nSel + (Long64_t)reconKeys.size());
    save(&hBlockAll, "Position at ZDC BLOCK, " + reconLabel + xy, "2_xy_block_recon.pdf", "colz");
    save(&hTrue, "Position at ZDC, " + cutLabel + xy, "3_xy.pdf", "colz");
    save(&hTrueAll, "Position at ZDC, " + reconLabel + xy, "4_xy_recon.pdf", "colz");

    // Recon proton histograms on log y: all recon protons (black), the CE ones among them (blue), and for the KE the phi diff window (red)
    auto drawOverlay = [&](TH1D& hRecon, TH1D& hCE, TH1D* hWin, const std::string& label, const std::string& what, const std::string& file,
                           const std::vector<std::string>& eqLines){
        hRecon.SetLineColor(kBlack); hCE.SetLineColor(kBlue);
        for (TH1D* h : {&hRecon, &hCE}) h->SetLineWidth(2);
        if (hWin){hWin->SetLineColor(kRed); hWin->SetLineWidth(2);}
        hRecon.SetStats(0);
        hRecon.SetMinimum(0.5); // log y: the recon protons outnumber the CE ones by far
        hRecon.SetMaximum(2.0 * std::max({hRecon.GetMaximum(), hCE.GetMaximum(), hWin ? hWin->GetMaximum() : 0.0}));
        hRecon.SetTitle(PlotTitle(file, label + ", " + what + ";" + label + ";Count").c_str());
        canvas->Clear();
        canvas->SetLogy(1);
        hRecon.Draw(); hCE.Draw("same"); if (hWin) hWin->Draw("same");
        TLegend leg(0.745, 0.75, 0.99, 0.88); leg.SetTextSize(0.016);
        leg.AddEntry(&hRecon, Form("all recon protons (%.0f)", hRecon.GetEntries()), "l");
        leg.AddEntry(&hCE, Form("CE recon protons in phi window (%.0f, %.2f%%)", hCE.GetEntries(), hRecon.GetEntries() > 0 ? 100.0 * hCE.GetEntries() / hRecon.GetEntries() : 0.0), "l");
        if (hWin) leg.AddEntry(hWin, Form("phi diff window (%.0f)", hWin->GetEntries()), "l");
        leg.Draw();
        TPaveText eq(0.745, 0.70 - 0.05 * (eqLines.size() + 1), 0.99, 0.70, "NDC");
        eq.SetBorderSize(1); eq.SetFillColor(0); eq.SetTextAlign(12); eq.SetTextSize(0.016);
        for (const auto& l : eqLines) eq.AddText(l.c_str());
        if (!eqLines.empty()) eq.Draw();
        WritePlotNumbers(&hRecon, folder + file, numbersStarted);
        canvas->Update();
        canvas->SaveAs((folder + file).c_str());
        canvas->SetLogy(0);
    };

    // KE: every recon proton and the CE ones in the phi diff window (finalKE of the Angles tree)
    {
        TH1D hKERecon("hKERecon", "", 100, 0, 30000), hKECE("hKECE", "", 100, 0, 30000);
        for (Long64_t i = 0; i < nAng; i++){
            hKERecon.Fill(angCols[2][i]);
            if (ceWinKeys.count({llround(angCols[0][i]), llround(angCols[1][i])})) hKECE.Fill(angCols[2][i]);
        }
        drawOverlay(hKERecon, hKECE, nullptr, "Kinetic Energy at ZDC (MeV)", Form("all recon protons and CE recon protons in phi diff %g to %g deg", phiLo, phiHi), "5_ke.pdf", {});
    }

    // Block line fit quality of every recon proton and the CE ones in the phi diff window: 6_resDRms.pdf, 7_r2D.pdf, 8_impactD.pdf
    struct QualVar {std::string label, file; double lo, hi; std::vector<std::string> eq;}; // eq: short lines, they sit in the right margin
    const std::string line = "fit to block hits: #hat{x} = x_{0} + b_{x}z, #hat{y} = y_{0} + b_{y}z";
    const std::string ssRes = "SS_{res} = #Sigma_{i} (x_{i} - #hat{x}_{i})^{2} + (y_{i} - #hat{y}_{i})^{2}";
    const std::vector<QualVar> vars = { // same order as columns 4, 5, 6 of angCols
        {"block residual (mm)", "6_resDRms.pdf", 0, 5, {line, "residual = #sqrt{SS_{res} / n}", ssRes}},
        {"block R^{2}", "7_r2D.pdf", 0, 1, {line, "R^{2} = 1 - SS_{res} / SS_{tot}", ssRes, "SS_{tot} = #Sigma_{i} (x_{i} - #bar{x})^{2} + (y_{i} - #bar{y})^{2}"}},
        {"block impact to beam axis (mm)", "8_impactD.pdf", 0, 20, {line, "impact = |x_{0} b_{y} - y_{0} b_{x}| / #sqrt{b_{x}^{2} + b_{y}^{2}}"}}};
    for (size_t v = 0; v < vars.size(); v++){
        TH1D hRecon(Form("hRecon%zu", v), "", 100, vars[v].lo, vars[v].hi), hCE(Form("hCE%zu", v), "", 100, vars[v].lo, vars[v].hi);
        for (Long64_t i = 0; i < nAng; i++){
            const double x = angCols[4 + v][i];
            if (!std::isfinite(x)) continue;
            hRecon.Fill(x);
            if (ceWinKeys.count({llround(angCols[0][i]), llround(angCols[1][i])})) hCE.Fill(x);
        }
        drawOverlay(hRecon, hCE, nullptr, vars[v].label, Form("all recon protons and CE recon protons in phi diff %g to %g deg", phiLo, phiHi), vars[v].file, vars[v].eq);
    }

    // The other cut candidates (6_cuts/3_ResCuts to 8_nHits), shown only, nothing is cut on them: the true position residual, impact and R^2, protons per event,
    // z order, residual of the 2nd half over the 1st half of the tracker hits, recon theta, recon block theta and tracker hits. 9_ to 17_
    // The Angles tree already dropped the tracks that fail the z order, so 13_zOrder.pdf only holds protons that pass it.
    if (ang->GetBranch("nProtonsEv")){
        const Long64_t nEx = cEx.N();
        if (nEx != nAng){Error("PlotDiffPhiCuts", "extra columns have %lld rows, not %lld", nEx, nAng); return;}
        struct ExVar {std::string label, file; double lo, hi; int nBins; std::vector<std::string> eq;};
        const std::string lineT = "fit to hits: #hat{x} = x_{0} + b_{x}z, #hat{y} = y_{0} + b_{y}z";
        const std::vector<ExVar> exVars = {
            {"residual, true position (mm)", "9_resRms.pdf", 0, 5, 100, {lineT, "residual = #sqrt{SS_{res} / n}"}},
            {"impact to beam axis, true position (mm)", "10_impact.pdf", 0, 20, 100, {lineT, "impact = |x_{0} b_{y} - y_{0} b_{x}| / #sqrt{b_{x}^{2} + b_{y}^{2}}"}},
            {"R^{2} of the line fit, true position", "11_r2.pdf", 0, 1, 100, {lineT, "R^{2} = 1 - SS_{res} / SS_{tot}"}},
            {"protons with a track in the event", "12_nProtons.pdf", -0.5, 19.5, 20, {}},
            {"planes hit in increasing z (1 yes, 0 no)", "13_zOrder.pdf", -0.5, 1.5, 2, {}},
            {"block residual, hits after first 2 / first 2 hits (>2 hits)", "14_resHalf.pdf", 0, 10, 100, {}},
            {"recon theta (degrees)", "15_theta.pdf", 0, 2.3, 100, {}},
            {"recon block theta (degrees)", "16_thetaBlock.pdf", 0, 2.3, 100, {}},
            {"tracker hits", "17_nHits.pdf", -0.5, 19.5, 20, {}}};
        std::vector<const double*> exCols(exVars.size());
        for (size_t c = 0; c < exVars.size(); c++) exCols[c] = cEx.V(c);
        for (size_t v = 0; v < exVars.size(); v++){
            TH1D hRecon(Form("hExRecon%zu", v), "", exVars[v].nBins, exVars[v].lo, exVars[v].hi), hCE(Form("hExCE%zu", v), "", exVars[v].nBins, exVars[v].lo, exVars[v].hi);
            for (Long64_t i = 0; i < nAng; i++){
                const double x = exCols[v][i];
                if (!std::isfinite(x)) continue;
                hRecon.Fill(x);
                if (ceWinKeys.count({llround(angCols[0][i]), llround(angCols[1][i])})) hCE.Fill(x);
            }
            drawOverlay(hRecon, hCE, nullptr, exVars[v].label, Form("all recon protons and CE recon protons in phi diff %g to %g deg", phiLo, phiHi), exVars[v].file, exVars[v].eq);
        }
    }
    BoxAndCropPlots();
}
