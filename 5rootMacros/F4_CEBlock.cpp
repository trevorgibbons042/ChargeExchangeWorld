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
#include <memory>
#include <set>
#include <string>
#include "BlockSquares.h"
#include "DrawPlan.h"
#include "FitCore.h"
#include "CutFiles.h"
#include <TKey.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TList.h>
#include <algorithm>
#include <vector>
#include <TROOT.h>
#include "TreeColumns.h"
#include <TH1D.h>

// ===== PlotDiffPhiCuts =====
// Events whose RECONPHIBL_ESTNEU is off from the CEPHI_NEUFRAME by absLo to absHi degrees either way, -absHi to -absLo or absLo to absHi (CalcValsReDDIFF.phiReD in cut_APandT.root).
// The block phi diff window is the only cut: no KE cut, no block residual cut. The xy plots show the position of the last tracker hit of the recon fit (lastXD, lastYD of AnglesRaw), matched on fEvent.
// energy: "fDEnergy" (detector energy), the only input trees cut_APandT.root has. Next to each window plot the same for every recon
// event, and the block fit residual of every recon event and the CE ones.
void PlotDiffPhiCuts(double absLo = 50, double absHi = 180, std::string energy = "fDEnergy"){
    FastJit();
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    const std::string folder = "plots/4_CEblock/3_diffphiBlock/";
    gSystem->mkdir(folder.c_str(), true);
    gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X");

    TFile* fAP = TFile::Open("cuts/cut_APandT.root", "READ");
    TTree* diff = (TTree*)fAP->Get(("CalcValsReDDIFF" + suffix).c_str());
    if (!diff){Error("PlotDiffPhiCuts", "missing tree in cut_APandT.root"); return;}

    // Every column pull of this macro is booked here and collected in one multithreaded pass per tree (DrawPlan.h), the render code below reads the columns in entry order.
    ROOT::EnableImplicitMT(NThreads());
    DrawPlan plan;
    plan.unordered = true; // nothing below depends on the row order, so AnglesRaw (many clusters) is read on all threads
    // CE events in the phi diff window, either sign: the only cut
    auto cDiff = plan.C(diff, {"fEvent"}, Form("fabs(phiReD) >= %g && fabs(phiReD) <= %g", absLo, absHi));
    // Every recon event (the rows of AnglesRaw in cut_angles_all.root): fEvent for the positions,
    // the KE of the energy (zdcKEAll / zdcKEAllD), isFromCE and the block line fit quality (ReconAngleAll.cpp)
    TTree* ang = OpenCut("cuts/cut_angles_all.root", kAngleTree); // with the Events tree as a friend: zdcKEAll, zdcKEAllD, inCE, nChargedParticle by name
    if (!ang){Error("PlotDiffPhiCuts", "no AnglesRaw in cuts/cut_angles_all.root (rerun ReconAngleAll.cpp)"); return;}
    const std::string angSel = AngleSelection("Angles"); // the same rows for both column sets, so they line up
    auto cAng = plan.C(ang, {"fEvent", KeVar(suffix), "isFromCE", "resDRms", "phiRelN", "phiRelND"}, angSel);
    auto cEx = plan.C(ang, {"resRms", "nChargedParticle", "resHalfRatioD", "thetaRelN*180/TMath::Pi()", "thetaRelND*180/TMath::Pi()", "nHits"}, angSel);
    // the position of the last tracker hit of every recon event, block (lastXD, lastYD): the window events and every recon event pick their rows by fEvent below
    auto cPos = plan.C(ang, {"fEvent", "lastXD", "lastYD"}, angSel);
    plan.Run();
    gROOT->cd();
    const Long64_t nSel = cDiff.N();
    const time_t tKeys = time(nullptr);
    // Key sets: sorted vectors of fEvent keys, not std::set, which is slow and large for every recon event of the run.
    struct KeySet {
        std::vector<uint32_t> k;
        void Finish(){std::sort(k.begin(), k.end()); k.erase(std::unique(k.begin(), k.end()), k.end());}
        size_t count(long long ev) const {return std::binary_search(k.begin(), k.end(), uint32_t(ev));}
        size_t size() const {return k.size();}
    };
    KeySet keys;
    for (Long64_t i = 0; i < nSel; i++){
        keys.k.push_back(uint32_t(llround(cDiff.V(0)[i])));
    }
    keys.Finish();

    const Long64_t nAng = cAng.N();
    std::vector<const double*> angCols(6); // views of the DrawPlan columns, no copy
    for (int c = 0; c < 6; c++) angCols[c] = cAng.V(c);
    // every recon event (CE ones included) in the same block phi diff window: the CE events keep their block phi - CEPHI_NEUFRAME
    // window above, the non-CE ones have no CEPHI_NEUFRAME so their diff is RECONPHIBL_ESTNEU - recon tracker phi (phiRelN), wrapped to -180 to 180
    // like CalcAPandT.cpp. ceWinKeys: the window CE events (fEvent), the AnglesRaw rows of an event in keys
    KeySet reconKeys, ceWinKeys;
    struct AngKeys {std::vector<uint32_t> ce, recon;};
    const auto angParts = MapRows(nAng, [&](Long64_t lo, Long64_t hi){ // row ranges on all threads (MapRows, FitCore.h), each collects its own keys
        AngKeys part;
        for (Long64_t i = lo; i < hi; i++){
            if (angCols[2][i] > 0.5){
                if (keys.count(llround(angCols[0][i]))) part.ce.push_back(uint32_t(llround(angCols[0][i])));
                continue;
            }
            const double dPhi = std::remainder((angCols[5][i] - angCols[4][i]) * 180 / TMath::Pi(), 360.0);
            if (fabs(dPhi) >= absLo && fabs(dPhi) <= absHi) part.recon.push_back(uint32_t(llround(angCols[0][i])));
        }
        return part;
    });
    for (const auto& part : angParts){
        ceWinKeys.k.insert(ceWinKeys.k.end(), part.ce.begin(), part.ce.end());
        reconKeys.k.insert(reconKeys.k.end(), part.recon.begin(), part.recon.end());
    }
    ceWinKeys.Finish(); reconKeys.Finish();

    // positions at the last tracker, block (lastXD, lastYD): the window events and every recon event, rows of AnglesRaw matched on fEvent
    // the block maps have one bin per tracker cell (nCols_blocks x nRows_blocks over 400 mm), drawn as one square per cell (BlockSquares.h)
    const int nColsD = (int)ReadSettings().at("nCols_blocks"), nRowsD = (int)ReadSettings().at("nRows_blocks");
    TH2D hBlock("hBlock", "", nColsD, -200, 200, nRowsD, -200, 200);
    TH2D hBlockAll("hBlockAll", "", nColsD, -200, 200, nRowsD, -200, 200);
    // the rows of a column set whose fEvent is in a key set, found on all threads; the histograms are then filled here in row order, so the stats sums are those of a serial fill
    auto matching = [&](const DrawPlan::Cols& c, const KeySet& set){
        const auto parts = MapRows(c.N(), [&](Long64_t lo, Long64_t hi){
            std::vector<Long64_t> rows;
            for (Long64_t i = lo; i < hi; i++) if (set.count(llround(c.V(0)[i]))) rows.push_back(i);
            return rows;
        });
        std::vector<Long64_t> rows;
        for (const auto& p : parts) rows.insert(rows.end(), p.begin(), p.end());
        return rows;
    };
    for (Long64_t i : matching(cPos, keys)){
        hBlock.Fill(cPos.V(1)[i], cPos.V(2)[i]);
    }
    for (Long64_t i : matching(cPos, reconKeys)){
        hBlockAll.Fill(cPos.V(1)[i], cPos.V(2)[i]);
    }
    printf(">>> PlotDiffPhiCuts key sets and position fills took %.0f s\n", difftime(time(nullptr), tKeys)); fflush(stdout);
    hBlockAll.Add(&hBlock);

    TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
    gPad->SetRightMargin(0.30);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);

    const std::string cutLabel = Form("CE recon events, RECONPHIBL_ESTNEU - CEPHI_NEUFRAME -%g to -%g and %g to %g deg, no other cut, N = %lld", absHi, absLo, absLo, absHi, nSel);
    const std::string overlayWhat = Form("all recon events and CE recon events, RECONPHIBL_ESTNEU - CEPHI_NEUFRAME -%g to -%g and %g to %g deg", absHi, absLo, absLo, absHi); // the title of the log y overlays
    std::set<std::string> numbersStarted;
    auto save = [&](TH1* h, const std::string& title, const std::string& file, const char* opt, bool squares = false){
        h->SetTitle(PlotTitle(file, title).c_str());
        canvas->Clear();
        h->Draw(opt);
        if (squares) DrawBlockSquares(h, 0, false);
        WritePlotNumbers(h, folder + file, numbersStarted);
        canvas->Update();
        canvas->SaveAs((folder + file).c_str());
    };
    const std::string xy = ";X Position (mm);Y Position (mm)";
    save(&hBlock, "Position at last tracker BLOCK, " + cutLabel + ", 1 square per tracker cell" + xy, "1_xy_block.pdf", "colz", true);
    const std::string reconLabel = Form("all recon events, RECONPHIBL_ESTNEU - CEPHI_NEUFRAME -%g to -%g and %g to %g deg (CE: RECONPHIBL_ESTNEU - CEPHI_NEUFRAME, non-CE: RECONPHIBL_ESTNEU - RECONPHI_ESTNEU), no other cut, N = %lld", absHi, absLo, absLo, absHi,
                                            nSel + (Long64_t)reconKeys.size());
    save(&hBlockAll, "Position at last tracker BLOCK, " + reconLabel + ", 1 square per tracker cell" + xy, "2_xy_block_recon.pdf", "colz", true);

    // Recon histograms on log y: all recon events (green, every row of AnglesRaw) and the CE recon events in the phi diff window (red)
    auto drawOverlay = [&](TH1D& hRecon, TH1D& hCE, const std::string& label, const std::string& what, const std::string& file,
                           const std::vector<std::string>& eqLines){
        hRecon.SetLineColor(CbAll()); hCE.SetLineColor(CbCE());
        for (TH1D* h : {&hRecon, &hCE}) h->SetLineWidth(2);
        hRecon.SetStats(0);
        hRecon.SetMinimum(0.5); // log y: the recon events outnumber the CE ones by far
        hRecon.SetMaximum(2.0 * std::max({hRecon.GetMaximum(), hCE.GetMaximum()}));
        hRecon.SetTitle(PlotTitle(file, label + ", " + what + ";" + label + ";Count").c_str());
        canvas->Clear();
        canvas->SetLogy(1);
        hRecon.Draw(); hCE.Draw("same");
        TLegend leg(0.745, 0.72, 0.99, 0.88); leg.SetTextSize(0.014);
        leg.AddEntry(&hRecon, Form("all recon events, no cuts (%.0f)", hRecon.GetEntries()), "l");
        leg.AddEntry(&hCE, Form("CE recon events in window (%.0f, %.2f%%)", hCE.GetEntries(), hRecon.GetEntries() > 0 ? 100.0 * hCE.GetEntries() / hRecon.GetEntries() : 0.0), "l");
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

    // Per row, once for all the overlays below: 1 = a CE event in the window (ceWinKeys). Row ranges on all threads, each writes its own part.
    std::vector<unsigned char> rowFlag(nAng);
    MapRows(nAng, [&](Long64_t lo, Long64_t hi){
        for (Long64_t i = lo; i < hi; i++)
            rowFlag[i] = ceWinKeys.count(llround(angCols[0][i])) ? 1 : 0;
        return 0;
    });
    // Fills hRecon (every row with a finite x) and hCE (the CE ones in the window) from the column x. Each row range fills its own copies on all threads (MapRows), added
    // into the two histograms in range order after: the bin contents are exact, the mean and RMS sums can differ from a serial fill in the last digits.
    struct Hist2 {std::shared_ptr<TH1D> recon, ce;};
    auto fill2 = [&](const double* x, TH1D& hRecon, TH1D& hCE){
        const bool addDir = TH1::AddDirectoryStatus(); TH1::AddDirectory(false); // the copies of the ranges are not registered in the current directory
        const auto parts = MapRows(nAng, [&](Long64_t lo, Long64_t hi){
            Hist2 p{std::make_shared<TH1D>(), std::make_shared<TH1D>()};
            hRecon.Copy(*p.recon); hCE.Copy(*p.ce); // same binning, still empty
            for (Long64_t i = lo; i < hi; i++){
                if (!std::isfinite(x[i])) continue;
                p.recon->Fill(x[i]);
                if (rowFlag[i]) p.ce->Fill(x[i]);
            }
            return p;
        });
        TH1::AddDirectory(addDir);
        for (const Hist2& p : parts){hRecon.Add(p.recon.get()); hCE.Add(p.ce.get());}
    };

    // KE: every recon event and the CE ones in the phi diff window (zdcKEAll of the Angles tree)
    {
        const double keMax = KeAxisMax(); // 30000 MeV at 25 GeV/c, FitCore.h
        TH1D hKERecon("hKERecon", "", 100, 0, keMax), hKECE("hKECE", "", 100, 0, keMax);
        fill2(angCols[1], hKERecon, hKECE); // no ZDC hit: no KE, not filled
        drawOverlay(hKERecon, hKECE, "Kinetic Energy at ZDC (MeV)", overlayWhat, "5_ke.pdf", {});
    }

    // Block line fit quality of every recon event and the CE ones in the phi diff window: 6_resDRms.pdf
    struct QualVar {std::string label, file; double lo, hi; std::vector<std::string> eq;}; // eq: short lines, they sit in the right margin
    const std::string line = "fit to block hits: #hat{x} = x_{0} + b_{x}z, #hat{y} = y_{0} + b_{y}z";
    const std::string ssRes = "SS_{res} = #Sigma_{i} (x_{i} - #hat{x}_{i})^{2} + (y_{i} - #hat{y}_{i})^{2}";
    const std::vector<QualVar> vars = { // same order as column 3 of angCols
        {"block residual (mm)", "6_resDRms.pdf", 0, 5, {line, "residual = #sqrt{SS_{res} / n}", ssRes}}};
    for (size_t v = 0; v < vars.size(); v++){
        TH1D hRecon(Form("hRecon%zu", v), "", 100, vars[v].lo, vars[v].hi), hCE(Form("hCE%zu", v), "", 100, vars[v].lo, vars[v].hi);
        fill2(angCols[3 + v], hRecon, hCE);
        drawOverlay(hRecon, hCE, vars[v].label, overlayWhat, vars[v].file, vars[v].eq);
    }

    // The other cut candidates, shown only, nothing is cut on them: the true position residual, protons per event,
    // residual of the 2nd half over the 1st half of the tracker hits, RECONTHETA_ESTNEU, RECONTHETABL_ESTNEU and tracker hits. 9_ to 17_
    const Long64_t nEx = cEx.N();
    if (nEx != nAng){Error("PlotDiffPhiCuts", "extra columns have %lld rows, not %lld", nEx, nAng); return;}
    struct ExVar {std::string label, file; double lo, hi; int nBins; std::vector<std::string> eq;};
    const std::string lineT = "fit to hits: #hat{x} = x_{0} + b_{x}z, #hat{y} = y_{0} + b_{y}z";
    const std::vector<ExVar> exVars = {
        {"residual, true position (mm)", "9_resRms.pdf", 0, 5, 100, {lineT, "residual = #sqrt{SS_{res} / n}"}},
        {"tracks of the event with a tracker hit", "12_nParticles.pdf", -0.5, 19.5, 20, {}},
        {"block residual, hits after first 2 / first 2 hits (>2 hits)", "14_resHalf.pdf", 0, 10, 100, {}},
        {"RECONTHETA_ESTNEU (degrees)", "15_theta.pdf", 0, ThetaAxisMax(), 100, {}}, // 2.3 deg at 25 GeV/c, scaled with 1 / p (FitCore.h)
        {"RECONTHETABL_ESTNEU (degrees)", "16_thetaBlock.pdf", 0, ThetaAxisMax(), 100, {}},
        {"tracker hits", "17_nHits.pdf", -0.5, MaxTrackerHits() + 2.5, MaxTrackerHits() + 3, {}}}; // one bin per count from 0 to 2 * nBlocks + 2
    std::vector<const double*> exCols(exVars.size());
    for (size_t c = 0; c < exVars.size(); c++) exCols[c] = cEx.V(c);
    for (size_t v = 0; v < exVars.size(); v++){
        TH1D hRecon(Form("hExRecon%zu", v), "", exVars[v].nBins, exVars[v].lo, exVars[v].hi), hCE(Form("hExCE%zu", v), "", exVars[v].nBins, exVars[v].lo, exVars[v].hi);
        fill2(exCols[v], hRecon, hCE);
        drawOverlay(hRecon, hCE, exVars[v].label, overlayWhat, exVars[v].file, exVars[v].eq);
    }
    BoxAndCropPlots();
}

// ===== GraphBlocks =====
void GraphBlocks(){
    FastJit();
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    const TString filename = "cuts/cut_blocks.root";
    const TString outDir = "plots/4_CEblock";
    const int nBins = 100;
    const char *ceBase = "all CE protons (sim, reconstructed or not)";
    const TString ceLabel = TString(ceBase) + ", no cuts";

    TFile *file = TFile::Open(filename, "READ");

    std::vector<TTree*> blocks;
    for (auto key : *file->GetListOfKeys()) {
        TTree *t = dynamic_cast<TTree*>(((TKey*)key)->ReadObj());
        if (t && TString(t->GetName()).BeginsWith("Block")) blocks.push_back(t);
    }
    if (blocks.empty()) return;

    gSystem->mkdir(outDir, true);
    TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
    gPad->SetRightMargin(0.30);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); // room for the rotated Y title and bottom tick labels
    gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X"); // stats/fit box to the right, fewer x labels

    // 2x2 style grid (more blocks make it grow)
    int nCols = std::ceil(std::sqrt((double)blocks.size()));
    int nRows = std::ceil(blocks.size() / (double)nCols);
    const int gridW = std::max(1500, nCols * 750), gridH = std::max(1200, nRows * 600); // 750 x 600 px per pad: the 2 x 2 grid of 4 blocks stays 1500 x 1200, more blocks get a larger page
    TCanvas* gridCanvas = new TCanvas("gridCanvas", "CE per block", gridW, gridH);
    gridCanvas->Divide(nCols, nRows);

    // every branch the plots and counts below use, read once per block tree by one multithreaded Columns call (TreeColumns.h, rows in entry order)
    const std::vector<std::string> colNames = {"fZ", "fX", "fY", "trackerX", "trackerY", "detectorX", "detectorY", "trackerDX", "trackerDY", "detectorDX", "detectorDY", "reachedDetector", "reachedNextTracker", "fEvent"};
    auto colIdx = [&](const std::string& n){return int(std::find(colNames.begin(), colNames.end(), n) - colNames.begin());};
    const int iX = colIdx("fX"), iY = colIdx("fY"), iZ = colIdx("fZ"), iReachedDet = colIdx("reachedDetector"), iReachedTrk = colIdx("reachedNextTracker"), iEvent = colIdx("fEvent");
    std::vector<std::vector<std::vector<double>>> blockCols(blocks.size()); // [block][column][row], kept for the recon grids
    std::vector<double> blockZLo(blocks.size(), 0.0), blockZHi(blocks.size(), 0.0);

    // TTree::Draw gave its histogram the line, fill and marker attributes of the tree (black line, not the style's blue one): do the same for the filled histograms
    auto drawStyle = [](TH1* h, TTree* t){
        h->SetLineColor(t->GetLineColor()); h->SetLineStyle(t->GetLineStyle()); h->SetLineWidth(t->GetLineWidth());
        h->SetFillColor(t->GetFillColor()); h->SetFillStyle(t->GetFillStyle());
        h->SetMarkerColor(t->GetMarkerColor()); h->SetMarkerStyle(t->GetMarkerStyle()); h->SetMarkerSize(t->GetMarkerSize());
    };

    std::set<std::string> numberFiles; // stats csv files already truncated in this run
    for (size_t b = 0; b < blocks.size(); b++) {
        TTree *t = blocks[b];
        TString name = t->GetName();
        TString dir = Form("%s/%s", outDir.Data(), name.Data());
        gSystem->mkdir(dir, true);

        blockCols[b] = Columns(t, colNames);
        const std::vector<std::vector<double>>& cols = blockCols[b];
        double zLo = 0, zHi = 0;
        if (!cols[iZ].empty()){zLo = *std::min_element(cols[iZ].begin(), cols[iZ].end()); zHi = *std::max_element(cols[iZ].begin(), cols[iZ].end());}
        blockZLo[b] = zLo; blockZHi[b] = zHi;
        double pad = 0.05 * (zHi - zLo) + 1.0;
        Long64_t nCE = cols[iZ].size();
        Long64_t nZDC = std::count(cols[iReachedDet].begin(), cols[iReachedDet].end(), 1.0);
        TString rate = Form("ZDC/CE = %lld/%lld = %.4f", nZDC, nCE, nCE ? (double)nZDC / nCE : 0.0);
        Long64_t nTracker = std::count(cols[iReachedTrk].begin(), cols[iReachedTrk].end(), 1.0);
        TString trackerRate = Form("Tracker/CE = %lld/%lld = %.4f", nTracker, nCE, nCE ? (double)nTracker / nCE : 0.0);

        TString zHist = Form("fZ>>h(%d,%f,%f)", nBins, zLo - pad, zHi + pad);

        struct ZPlot { const char *file; const char *label; };
        ZPlot zPlots[] = {
            {"1_CE_z", "Z Position of Charge Exchange"},
        };
        for (auto &p : zPlots) {
            canvas->cd();
            canvas->Clear();
            TH1* hz = MakeHist(SplitDraw(zHist.Data()).target); // no cut in zPlots: every row fills
            drawStyle(hz, t);
            hz->SetLineColor(CbProton()); // blue, the colorblind-safe blue of CbProton()
            for (double z : cols[iZ]) hz->Fill(z);
            hz->Draw("");
            TString zTitle = Form("%s, %s, %s: %lld", p.label, ceLabel.Data(), name.Data(), nCE);
            hz->SetTitle(PlotTitle(p.file, (zTitle + ";Z Position (mm);Counts").Data()).c_str());
            canvas->Update();
            WritePlotNumbers(hz, Form("%s/%s.pdf", dir.Data(), p.file), numberFiles);
            canvas->SaveAs(Form("%s/%s.pdf", dir.Data(), p.file));
        }

        // same CE z plot into the grid
        gridCanvas->cd(b + 1);
        gPad->SetRightMargin(0.05);
        TH1* hGrid = MakeHist(Form("hGrid%zu(%d,%f,%f)", b + 1, nBins, zLo - pad, zHi + pad));
        drawStyle(hGrid, t);
        for (double z : cols[iZ]) hGrid->Fill(z);
        hGrid->Draw("");
        hGrid->SetTitle(Form("Z Position of Charge Exchange, %s, %s: %lld;Z Position (mm);Counts", ceLabel.Data(), name.Data(), nCE));
        hGrid->SetStats(0); // small grid pads have no margin for a stats box, it would sit on the graph

        const double pitchT = 400.0 / ReadSettings().at("nCols_blocks"), pitchD = 400.0 / ReadSettings().at("nCols_detector"); // cell pitch, mm
        struct XYPlot { const char *file; TString expr; const char *cut; const char *label; };
        XYPlot xyPlots[] = {
            {"2_CE_xy",       "fY:fX>>h2(100,-50,50,100,-50,50)",                 "",                  "Position of Charge Exchange"},
            {"3_Tracker_xy",  "trackerY:trackerX>>h2(100,-50,50,100,-50,50)",    "reachedNextTracker==1", "Position at Next Tracker"},
            {"4_ZDC_xy", "detectorY:detectorX>>h2(100,-50,50,100,-50,50)",   "reachedDetector==1", "Position at ZDC"},
            // the detector's version of 3_ and 4_: block cell centers over the range of 25 cells centered on the beam axis, one bin per cell (25x25), drawn as one square per cell (BlockSquares.h)
            {"5_Tracker_xy_detector", Form("trackerDY:trackerDX>>h2(25,%g,%g,25,%g,%g)", -12.5 * pitchT, 12.5 * pitchT, -12.5 * pitchT, 12.5 * pitchT), "reachedNextTracker==1", "Position at Next Tracker, detector block position"},
            {"6_ZDC_xy_detector", Form("detectorDY:detectorDX>>h2(25,%g,%g,25,%g,%g)", -12.5 * pitchD, 12.5 * pitchD, -12.5 * pitchD, 12.5 * pitchD), "reachedDetector==1", "Position at ZDC, detector block position"},
        };
        for (auto &p : xyPlots) {
            canvas->cd();
            canvas->Clear();
            const DrawSpec ds = SplitDraw(p.expr.Data()); // exprs {y, x}; the cut is "<branch>==1" or empty
            TH1* hxy = MakeHist(ds.target);
            drawStyle(hxy, t);
            const std::vector<double>& xs = cols[colIdx(ds.exprs[1])]; const std::vector<double>& ys = cols[colIdx(ds.exprs[0])];
            const std::string cutBranch = TString(p.cut).Contains("==") ? std::string(p.cut).substr(0, std::string(p.cut).find("==")) : "";
            const std::vector<double>* mask = cutBranch.empty() ? nullptr : &cols[colIdx(cutBranch)];
            for (size_t k = 0; k < xs.size(); k++) if (!mask || (*mask)[k] == 1) ((TH2*)hxy)->Fill(xs[k], ys[k]);
            hxy->Draw("colz");
            const bool squares = TString(p.file).EndsWith("_detector"); // block positions: one square per cell
            if (squares) DrawBlockSquares(hxy, 0, false);
            const bool atTracker = TString(p.file).Contains("Tracker"), atZdc = TString(p.file).Contains("ZDC");
            const TString xyCuts = atTracker ? TString(ceBase) + ", no cuts besides reaching the tracker" : atZdc ? TString(ceBase) + ", no cuts besides reaching the ZDC" : ceLabel;
            TString xyTitle = Form("%s, %s, %s", p.label, xyCuts.Data(), name.Data());
            if (TString(p.file) == "4_ZDC_xy" || TString(p.file) == "6_ZDC_xy_detector") xyTitle += ", " + rate;
            if (TString(p.file) == "3_Tracker_xy" || TString(p.file) == "5_Tracker_xy_detector") xyTitle += ", " + trackerRate;
            if (squares) xyTitle += ", 1 square per cell";
            hxy->SetTitle(PlotTitle(p.file, (xyTitle + ";X Position (mm);Y Position (mm)").Data()).c_str());
            canvas->Update();
            WritePlotNumbers(hxy, Form("%s/%s.pdf", dir.Data(), p.file), numberFiles);
            canvas->SaveAs(Form("%s/%s.pdf", dir.Data(), p.file));

            // the same rows as a histogram of the radius r = sqrt(x^2 + y^2) from the beam axis (2_CE_xy -> 2_CE_r, 5_Tracker_xy_detector -> 5_Tracker_r_detector)
            std::string rFile = p.file; rFile.replace(rFile.find("_xy"), 3, "_r");
            std::vector<double> rs;
            for (size_t k = 0; k < xs.size(); k++) if (!mask || (*mask)[k] == 1) rs.push_back(std::hypot(xs[k], ys[k]));
            const double rMax = rs.empty() ? 1.0 : 1.05 * *std::max_element(rs.begin(), rs.end()) + 1e-6;
            canvas->Clear();
            TH1D* hr = new TH1D(Form("hR_%s_%zu", rFile.c_str(), b), "", nBins, 0, rMax);
            drawStyle(hr, t);
            for (double r : rs) hr->Fill(r);
            hr->Draw("");
            hr->SetTitle(PlotTitle(rFile, (std::string(p.label) + ", radius, " + std::string(xyCuts.Data()) + ", " + name.Data() + ";r = #sqrt{X^{2} + Y^{2}} (mm);Counts").c_str()).c_str());
            canvas->Update();
            WritePlotNumbers(hr, Form("%s/%s.pdf", dir.Data(), rFile.c_str()), numberFiles);
            canvas->SaveAs(Form("%s/%s.pdf", dir.Data(), rFile.c_str()));
        }
        for (size_t k = 0; k < colNames.size(); k++) // the recon grids need only these columns
            if ((int)k != iX && (int)k != iY && (int)k != iZ && (int)k != iReachedDet && (int)k != iEvent) std::vector<double>().swap(blockCols[b][k]);
    }

    gridCanvas->Update();
    gridCanvas->SaveAs(Form("%s/CE_per_block.pdf", outDir.Data()));
    std::set<std::string> gridStats; // one row per block: the CE count in each panel's title
    for (TTree* t : blocks) WriteStatsRow(Form("%s/CE_per_block.pdf", outDir.Data()), "plot,block,ce_events", Form("CE_per_block.pdf,%s,%lld", t->GetName(), (Long64_t)t->GetEntries()), gridStats);

    // Reconstructed CE: the CE tracks of the events with a recon fit in AnglesCE (cut_angles_all.root, no cuts). Same z grid as above with only these CE,
    // then only those that also reached the ZDC.
    std::set<int> reconKeys; // the recon events (fEvent) with a carbon CE
    {
        TTree *a = OpenCut("cuts/cut_angles_all.root", kAngleTree); // opening it attaches the Events tree (EventsTree())
        if (!a || !EventsTree() || !EventsTree()->GetBranch("inCE")){Error("GraphBlocks", "AnglesRaw or the Events tree with inCE missing in cuts/ (rerun ReconAngleAll.cpp)"); return;}
        // the AnglesCE view is the events with inCE set that have a recon fit (nRecoEv > 0, the events of AnglesRaw): read from the Events tree, one row per event, not from every AnglesRaw row
        const auto keys = ColumnsUnordered(EventsTree(), {"fEvent"}, "inCE != 0 && nRecoEv > 0"); // multithreaded read in any order (TreeColumns.h), the keys go into a set
        for (size_t i = 0; i < keys[0].size(); i++) reconKeys.insert((int)keys[0][i]);
    }
    struct ReconGrid {const char *file, *label; bool zdcOnly;};
    for (const ReconGrid &g : {ReconGrid{"CE_per_block_recon.pdf", "reconstructed CE (inCE, no other cuts)", false}, ReconGrid{"CE_per_block_recon_zdc.pdf", "reconstructed CE (inCE) that hit ZDC, no other cuts", true}}){
        TCanvas* reconCanvas = new TCanvas(Form("canvas_%s", g.file), g.label, gridW, gridH);
        reconCanvas->Divide(nCols, nRows);
        const TString path = Form("%s/%s", outDir.Data(), g.file);
        for (size_t b = 0; b < blocks.size(); b++){
            TTree *t = blocks[b];
            const std::vector<std::vector<double>>& cols = blockCols[b];
            std::vector<char> keep(cols[iZ].size(), 0); // the rows of the CE tracks of the events with a recon fit (and, zdcOnly, that reached the ZDC)
            Long64_t nKept = 0;
            for (size_t i = 0; i < keep.size(); i++)
                if (reconKeys.count((int)cols[iEvent][i]) && (!g.zdcOnly || cols[iReachedDet][i] == 1)){keep[i] = 1; nKept++;}

            double zLo = blockZLo[b], zHi = blockZHi[b];
            double pad = 0.05 * (zHi - zLo) + 1.0;
            reconCanvas->cd(b + 1);
            gPad->SetRightMargin(0.05);
            TH1* hRecon = MakeHist(Form("hRecon%d_%zu(%d,%f,%f)", (int)g.zdcOnly, b + 1, nBins, zLo - pad, zHi + pad));
            drawStyle(hRecon, t);
            for (size_t i = 0; i < keep.size(); i++) if (keep[i]) hRecon->Fill(cols[iZ][i]);
            hRecon->Draw("");
            hRecon->SetTitle(Form("Z Position of Charge Exchange, %s, %s: %lld;Z Position (mm);Counts", g.label, t->GetName(), nKept));
            hRecon->SetStats(0);
            WriteStatsRow(path.Data(), "plot,block,ce_events", Form("%s,%s,%lld", g.file, t->GetName(), nKept), gridStats);
        }
        reconCanvas->Update();
        reconCanvas->SaveAs(path);
    }

    // The same three selections as 2D color plots: Z of the CE (x) against its radius r = sqrt(fX^2 + fY^2) from the beam axis (y), one pad per block
    struct ZRGrid {const char *file, *label; int mode;}; // mode 0 every CE, 1 the CE with a reconstructed proton (inCE), 2 those that also reached the ZDC
    for (const ZRGrid &g : {ZRGrid{"CE_per_block_zr.pdf", ceLabel.Data(), 0}, ZRGrid{"CE_per_block_recon_zr.pdf", "reconstructed CE (inCE, no other cuts)", 1}, ZRGrid{"CE_per_block_recon_zdc_zr.pdf", "reconstructed CE (inCE) that hit ZDC, no other cuts", 2}}){
        TCanvas* zrCanvas = new TCanvas(Form("canvas_%s", g.file), g.label, gridW, gridH);
        zrCanvas->Divide(nCols, nRows);
        const TString path = Form("%s/%s", outDir.Data(), g.file);
        for (size_t b = 0; b < blocks.size(); b++){
            const std::vector<std::vector<double>>& cols = blockCols[b];
            std::vector<double> zs, rs; // the kept rows
            for (size_t i = 0; i < cols[iZ].size(); i++){
                if (g.mode >= 1 && !reconKeys.count((int)cols[iEvent][i])) continue;
                if (g.mode == 2 && cols[iReachedDet][i] != 1) continue;
                zs.push_back(cols[iZ][i]); rs.push_back(std::hypot(cols[iX][i], cols[iY][i]));
            }
            const double zLo = blockZLo[b], zHi = blockZHi[b], pad = 0.05 * (zHi - zLo) + 1.0;
            const double rHi = rs.empty() ? 1.0 : 1.05 * *std::max_element(rs.begin(), rs.end()) + 1e-6;
            zrCanvas->cd(b + 1);
            gPad->SetRightMargin(0.15);
            TH2D* hZR = new TH2D(Form("hZR%d_%zu", g.mode, b + 1), "", nBins, zLo - pad, zHi + pad, 50, 0, rHi);
            for (size_t i = 0; i < zs.size(); i++) hZR->Fill(zs[i], rs[i]);
            hZR->SetStats(0);
            hZR->SetTitle(Form("Z and radius of Charge Exchange, %s, %s: %zu;Z Position (mm);r = #sqrt{X^{2} + Y^{2}} (mm)", g.label, blocks[b]->GetName(), zs.size()));
            hZR->Draw("colz");
            WriteStatsRow(path.Data(), "plot,block,ce_events", Form("%s,%s,%zu", g.file, blocks[b]->GetName(), zs.size()), gridStats);
        }
        zrCanvas->Update();
        zrCanvas->SaveAs(path);
    }

    // The CE of each block over all the neutrons the sim shot (TotalNeutrons of settings.csv): every CE of the block, and the CE with a reconstructed proton (inCE, no other cuts). One point per block, binomial error.
    {
        const double nNeutrons = ReadSettings().at("TotalNeutrons");
        std::vector<double> nCEb(blocks.size(), 0.0), nReconB(blocks.size(), 0.0);
        for (size_t b = 0; b < blocks.size(); b++){
            nCEb[b] = (double)blocks[b]->GetEntries();
            for (size_t i = 0; i < blockCols[b][iZ].size(); i++) if (reconKeys.count((int)blockCols[b][iEvent][i])) nReconB[b] += 1;
        }
        struct RatioPage {const char *file, *label; const std::vector<double>* n;};
        for (const RatioPage& g : {RatioPage{"CE_total_over_neutrons.pdf", "N_{CE,total}", &nCEb}, RatioPage{"CE_recon_over_neutrons.pdf", "N_{CE,recon}", &nReconB}}){
            TCanvas c("cRatio", "Plots", 1500, 1200);
            c.SetRightMargin(0.05); c.SetLeftMargin(0.14); c.SetBottomMargin(0.13);
            const int nb = (int)blocks.size();
            TGraphErrors gr(nb + 1); // the blocks, then the total of all of them at the far right
            double ymax = 0, sumN = 0;
            for (int b = 0; b < nb; b++){
                const double p = (*g.n)[b] / nNeutrons, e = std::sqrt(std::max(0.0, p * (1 - p) / nNeutrons));
                gr.SetPoint(b, b + 1, p); gr.SetPointError(b, 0, e);
                ymax = std::max(ymax, p + e); sumN += (*g.n)[b];
            }
            {
                const double p = sumN / nNeutrons, e = std::sqrt(std::max(0.0, p * (1 - p) / nNeutrons));
                gr.SetPoint(nb, nb + 1, p); gr.SetPointError(nb, 0, e);
                ymax = std::max(ymax, p + e);
            }
            TH1D frame("hRatioFrame", "", nb + 1, 0.5, nb + 1.5);
            for (int b = 0; b < nb; b++) frame.GetXaxis()->SetBinLabel(b + 1, blocks[b]->GetName());
            frame.GetXaxis()->SetBinLabel(nb + 1, "Total");
            frame.SetStats(0); frame.SetMinimum(0); frame.SetMaximum(1.2 * ymax);
            const std::string path = std::string(outDir.Data()) + "/" + g.file;
            frame.SetTitle(PlotTitle(g.file, Form("%s / N_{neutrons} per block, N_{neutrons} = %.0f (TotalNeutrons)%s;Block;%s / N_{neutrons}", g.label, nNeutrons, g.n == &nReconB ? ", reconstructed CE (inCE, no other cuts)" : ", every CE of the block", g.label)).c_str());
            frame.Draw();
            gr.SetMarkerStyle(20); gr.SetMarkerSize(1.5); gr.SetMarkerColor(CbCE()); gr.SetLineColor(CbCE());
            gr.Draw("P same");
            c.SaveAs(path.c_str());
            for (int b = 0; b < nb; b++) WriteStatsRow(path, "plot,block,count,over_neutrons", Form("%s,%s,%.0f,%.6g", g.file, blocks[b]->GetName(), (*g.n)[b], (*g.n)[b] / nNeutrons), gridStats);
            WriteStatsRow(path, "plot,block,count,over_neutrons", Form("%s,Total,%.0f,%.6g", g.file, sumN, sumN / nNeutrons), gridStats);
        }
        { // both on one page: every CE of the block (blue) and the CE with a reconstructed proton (red, inCE, no other cuts)
            TCanvas c("cRatioOv", "Plots", 1500, 1200);
            c.SetRightMargin(0.05); c.SetLeftMargin(0.14); c.SetBottomMargin(0.13);
            const int nb = (int)blocks.size();
            TGraphErrors grAll(nb + 1), grRecon(nb + 1); // the blocks, then the total of all of them at the far right
            double ymax = 0, sumA = 0, sumR = 0;
            for (int b = 0; b <= nb; b++){
                if (b < nb){sumA += nCEb[b]; sumR += nReconB[b];}
                const double nA = b < nb ? nCEb[b] : sumA, nR = b < nb ? nReconB[b] : sumR;
                const double pA = nA / nNeutrons, pR = nR / nNeutrons;
                const double eA = std::sqrt(std::max(0.0, pA * (1 - pA) / nNeutrons)), eR = std::sqrt(std::max(0.0, pR * (1 - pR) / nNeutrons));
                grAll.SetPoint(b, b + 0.95, pA); grAll.SetPointError(b, 0, eA);
                grRecon.SetPoint(b, b + 1.05, pR); grRecon.SetPointError(b, 0, eR);
                ymax = std::max(ymax, pA + eA);
            }
            TH1D frame("hRatioOvFrame", "", nb + 1, 0.5, nb + 1.5);
            for (int b = 0; b < nb; b++) frame.GetXaxis()->SetBinLabel(b + 1, blocks[b]->GetName());
            frame.GetXaxis()->SetBinLabel(nb + 1, "Total");
            frame.SetStats(0); frame.SetMinimum(0); frame.SetMaximum(1.2 * ymax);
            const std::string path = std::string(outDir.Data()) + "/CE_over_neutrons_overlay.pdf";
            frame.SetTitle(PlotTitle("CE_over_neutrons_overlay.pdf", Form("N_{CE} / N_{neutrons} per block, N_{neutrons} = %.0f (TotalNeutrons), every CE of the block and reconstructed CE (inCE, no other cuts);Block;N_{CE} / N_{neutrons}", nNeutrons)).c_str());
            frame.Draw();
            grAll.SetMarkerStyle(20); grAll.SetMarkerSize(1.5); grAll.SetMarkerColor(CbProton()); grAll.SetLineColor(CbProton()); grAll.Draw("P same");
            grRecon.SetMarkerStyle(21); grRecon.SetMarkerSize(1.5); grRecon.SetMarkerColor(CbCE()); grRecon.SetLineColor(CbCE()); grRecon.Draw("P same");
            TLegend leg(0.60, 0.78, 0.94, 0.88); leg.SetTextSize(0.025);
            leg.AddEntry(&grAll, "all CE, N_{CE,total} / N_{neutrons}", "lp"); leg.AddEntry(&grRecon, "reconstructed CE, N_{CE,recon} / N_{neutrons}", "lp"); leg.Draw();
            c.SaveAs(path.c_str());
        }
    }
}

// ===== CETCompare =====
// t of the charge exchange (LastT) for the CE in the carbon blocks (CarbonCeCut) and for the CE in the silicon tracker planes (every other CE),
// ChargeExchange_values of the sim. Normalized to unit area in a second plot to compare the shapes. plots/4_CEblock/1_CarbonVsSi/t_carbon_vs_tracker.pdf, t_carbon_vs_tracker_norm.pdf and z_CE_blocks_vs_trackers.pdf (z of the CE, all, carbon and tracker)
void CETCompare(){
    FastJit();
    // all CE and the CE in the carbon blocks are cut files (CutFiles.h): every row of the tree passes, so the draws need no selection; the tracker CE is all CE without the carbon ones
    const std::string noCarbon = "!(" + CarbonCeCut() + ")";
    TTree *tAll = OpenCut("cuts/cut_CEall.root", "ChargeExchange_values");              // all CE
    TTree *tCarbon = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");              // CE in the carbon blocks
    if (!tAll || !tCarbon){Error("CETCompare", "cut files missing (run CutPlots.cpp (events) first)"); return;}
    gROOT->cd();
    ROOT::EnableImplicitMT(NThreads());
    TH1D hC("hTCarbon", "", 100, 0, 0.05), hT("hTTracker", "", 100, 0, 0.05);
    DrawPlan plan; // the five fills (z all/carbon/tracker, t carbon/tracker) are collected in one multithreaded pass per tree, then filled here in entry order
    gSystem->mkdir("plots/4_CEblock/1_CarbonVsSi", true);
    gStyle->SetOptStat(0);
    TCanvas c("cT", "", 1500, 1200); c.SetLeftMargin(0.14); c.SetBottomMargin(0.13); c.SetRightMargin(0.05);
    // z where the CE happens, all CE of the sim: carbon blocks, silicon trackers, and both
    {
        // the stack of the run (2construction.cc): nBlocks carbon blocks and the tracker planes (StackZ, FitCore.h), blocks evenly spread from z = 0 to the ZDC face, trackers either side, dz the gap between block and tracker edges (an old run: centred at 4 m, dz centre to centre)
        const auto st = ReadSettings();
        const StackLayout stack = StackZ(st);
        const double zFirst = std::min(stack.blockZ.front(), stack.planeZ.empty() ? 1e9 : stack.planeZ.front()), zLast = std::max(stack.blockZ.back(), stack.planeZ.empty() ? -1e9 : stack.planeZ.back());
        // the range also covers every CE row (the tracker CE can sit past the stack, e.g. at the ZDC end), so none falls in the overflow
        ROOT::RDataFrame dfZ(*tAll); // min and max of fZ in one multithreaded pass (GetMinimum and GetMaximum each read the whole branch on one thread)
        auto zMin = dfZ.Min<double>("fZ"), zMax = dfZ.Max<double>("fZ");
        const double z0 = std::min(zFirst - 200, *zMin - 50), z1 = std::max(zLast + 200, *zMax + 50); const int nz = 100;
        TH1D hAll("hZAll", "", nz, z0, z1), hZC("hZCarbon", "", nz, z0, z1), hZT("hZTracker", "", nz, z0, z1);
        plan.H(tAll, "fZ>>hZAll", "", &hAll); plan.H(tCarbon, "fZ>>hZCarbon", "", &hZC); plan.H(tAll, "fZ>>hZTracker", noCarbon, &hZT);
        plan.H(tCarbon, "LastT>>hTCarbon", "", &hC); plan.H(tAll, "LastT>>hTTracker", noCarbon, &hT);
        plan.Run();
        const long long nC = tCarbon->GetEntries(), nT = (long long)hZT.GetEntries(), nA = tAll->GetEntries();
        hAll.SetLineColor(kGray + 1); hZC.SetLineColor(CbCE()); hZT.SetLineColor(kBlack); // carbon CE red (CE protons are red), tracker CE black, the total gray so it stays apart from the tracker curve
        for (TH1D* h : {&hAll, &hZC, &hZT}) h->SetLineWidth(2);
        hAll.SetTitle("Z position where Charge Exchange happens, all CE protons (sim);Z Position (mm);Count");
        hAll.SetMaximum(1.15 * hAll.GetMaximum());
        hAll.Draw("hist"); hZC.Draw("hist same"); hZT.Draw("hist same");
        TLegend lz(0.62, 0.74, 0.95, 0.90); lz.SetTextSize(0.024);
        lz.AddEntry(&hAll, Form("all CE, %lld", nA), "l");
        lz.AddEntry(&hZC, Form("CE in carbon blocks, %lld", nC), "l");
        lz.AddEntry(&hZT, Form("CE in silicon trackers, %lld", nT), "l");
        lz.Draw();
        c.SaveAs("plots/4_CEblock/1_CarbonVsSi/z_CE_blocks_vs_trackers.pdf");
    }
    printf("carbon CE: N=%.0f mean t=%.5g | tracker CE: N=%.0f mean t=%.5g\n", hC.GetEntries(), hC.GetMean(), hT.GetEntries(), hT.GetMean());
    for (const bool norm : {false, true}){
        TH1D* a = (TH1D*)hC.Clone(), *b = (TH1D*)hT.Clone();
        if (norm){a->Scale(1.0 / a->Integral()); b->Scale(1.0 / b->Integral());}
        a->SetLineColor(CbCE()); b->SetLineColor(kBlack); a->SetLineWidth(2); b->SetLineWidth(2);
        a->SetTitle(Form("t of the charge exchange, all CE protons (sim), carbon blocks vs silicon trackers%s;t;%s", norm ? ", unit area" : "", norm ? "Fraction" : "Count"));
        a->SetMaximum(1.15 * std::max(a->GetMaximum(), b->GetMaximum()));
        a->Draw("hist"); b->Draw("hist same");
        TLegend leg(0.40, 0.78, 0.94, 0.90); leg.SetTextSize(0.025);
        leg.AddEntry(a, Form("CE in carbon blocks (%.0f), mean t = %.4g", hC.GetEntries(), hC.GetMean()), "l");
        leg.AddEntry(b, Form("CE in silicon trackers (%.0f), mean t = %.4g", hT.GetEntries(), hT.GetMean()), "l");
        leg.Draw();
        c.SaveAs(norm ? "plots/4_CEblock/1_CarbonVsSi/t_carbon_vs_tracker_norm.pdf" : "plots/4_CEblock/1_CarbonVsSi/t_carbon_vs_tracker.pdf");
    }
}

// 4_CEblock: the three steps of 4_CEblock, in order. absLo, absHi, energy go to PlotDiffPhiCuts.
void F4_CEBlock(double absLo = 50, double absHi = 180, std::string energy = "fDEnergy"){
    PlotDiffPhiCuts(absLo, absHi, energy);
    GraphBlocks();
    CETCompare();
}
