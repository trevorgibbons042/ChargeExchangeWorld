#include <TStyle.h>
#include "Palette.h"
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TSystem.h>
#include <TList.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <TROOT.h>
#include "BlockSquares.h"
#include "FitCore.h"
#include "CutFiles.h"
#include "TreeColumns.h"

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
    const std::vector<std::string> colNames = {"fZ", "fX", "fY", "trackerX", "trackerY", "detectorX", "detectorY", "trackerDX", "trackerDY", "detectorDX", "detectorDY", "reachedDetector", "reachedNextTracker", "fEvent", "TrackID"};
    auto colIdx = [&](const std::string& n){return int(std::find(colNames.begin(), colNames.end(), n) - colNames.begin());};
    const int iX = colIdx("fX"), iY = colIdx("fY"), iZ = colIdx("fZ"), iReachedDet = colIdx("reachedDetector"), iReachedTrk = colIdx("reachedNextTracker"), iEvent = colIdx("fEvent"), iTrack = colIdx("TrackID");
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

        if (FILE* out = fopen((dir + "/numbers.txt").Data(), "w")) {
            fprintf(out, "%s, %s\n", ceBase, name.Data());
            fprintf(out, "CE events: %lld\n", nCE);
            fprintf(out, "%s\n%s\n", trackerRate.Data(), rate.Data());
            fclose(out);
        }

        TString zHist = Form("fZ>>h(%d,%f,%f)", nBins, zLo - pad, zHi + pad);

        struct ZPlot { const char *file; const char *cut; const char *label; };
        ZPlot zPlots[] = {
            {"1_CE_z",        "",                     "Z Position of Charge Exchange"},
        };
        for (auto &p : zPlots) {
            canvas->cd();
            canvas->Clear();
            TH1* hz = MakeHist(SplitDraw(zHist.Data()).target); // no cut in zPlots: every row fills
            drawStyle(hz, t);
            hz->SetLineColor(CbCE()); // CE protons are red in the other plots
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
            if ((int)k != iX && (int)k != iY && (int)k != iZ && (int)k != iReachedDet && (int)k != iEvent && (int)k != iTrack) std::vector<double>().swap(blockCols[b][k]);
    }

    gridCanvas->Update();
    gridCanvas->SaveAs(Form("%s/CE_per_block.pdf", outDir.Data()));
    std::set<std::string> gridStats; // one row per block: the CE count in each panel's title
    for (TTree* t : blocks) WriteStatsRow(Form("%s/CE_per_block.pdf", outDir.Data()), "plot,block,ce_events", Form("CE_per_block.pdf,%s,%lld", t->GetName(), (Long64_t)t->GetEntries()), gridStats);

    // Reconstructed CE: the CE tracks (event, track id) with a reconstructed proton in AnglesCE (cut_angles_all.root, no cuts). Same z grid as above with only these CE,
    // then only those that also reached the ZDC.
    std::set<std::pair<int,int>> reconKeys;
    {
        TFile *fAng = TFile::Open("cuts/cut_angles_all.root", "READ");
        TTree *a = fAng ? (TTree*)fAng->Get(kAngleTree) : nullptr; // the AnglesCE view: AnglesRaw rows with inCE set
        if (!a || !a->GetBranch("inCE")){Error("GraphBlocks", "AnglesRaw with inCE missing in cuts/cut_angles_all.root (rerun ReconAngleAll.cpp)"); return;}
        const auto keys = ColumnsUnordered(a, {"fEvent", "ParentID"}, "inCE != 0"); // multithreaded read in any order (TreeColumns.h), the keys go into a set
        for (size_t i = 0; i < keys[0].size(); i++) reconKeys.insert({(int)keys[0][i], (int)keys[1][i]});
    }
    struct ReconGrid {const char *file, *label; bool zdcOnly;};
    for (const ReconGrid &g : {ReconGrid{"CE_per_block_recon.pdf", "reconstructed CE (inCE, no other cuts)", false}, ReconGrid{"CE_per_block_recon_zdc.pdf", "reconstructed CE (inCE) that hit ZDC, no other cuts", true}}){
        TCanvas* reconCanvas = new TCanvas(Form("canvas_%s", g.file), g.label, gridW, gridH);
        reconCanvas->Divide(nCols, nRows);
        const TString path = Form("%s/%s", outDir.Data(), g.file);
        for (size_t b = 0; b < blocks.size(); b++){
            TTree *t = blocks[b];
            const std::vector<std::vector<double>>& cols = blockCols[b];
            std::vector<char> keep(cols[iZ].size(), 0); // the rows of the CE tracks with a reconstructed proton (and, zdcOnly, that reached the ZDC)
            Long64_t nKept = 0;
            for (size_t i = 0; i < keep.size(); i++)
                if (reconKeys.count({(int)cols[iEvent][i], (int)cols[iTrack][i]}) && (!g.zdcOnly || cols[iReachedDet][i] == 1)){keep[i] = 1; nKept++;}

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
                if (g.mode >= 1 && !reconKeys.count({(int)cols[iEvent][i], (int)cols[iTrack][i]})) continue;
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
            for (size_t i = 0; i < blockCols[b][iZ].size(); i++) if (reconKeys.count({(int)blockCols[b][iEvent][i], (int)blockCols[b][iTrack][i]})) nReconB[b] += 1;
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
