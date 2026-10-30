#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TString.h>
#include <TSystem.h>
#include <TList.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <TROOT.h>
#include "FitCore.h"
#include "CutFiles.h"
#include "TreeColumns.h"

void GraphBlocks(){
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    const TString filename = "cuts/cut_blocks.root";
    const TString outDir = "plots/4_CEblock";
    const int nBins = 100;
    const char *ceLabel = "all CE protons (sim, reconstructed or not)";

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
    TCanvas* gridCanvas = new TCanvas("gridCanvas", "CE per block", 1500, 1200);
    gridCanvas->Divide(nCols, nRows);

    // every branch the plots and counts below use, read once per block tree by one multithreaded Columns call (TreeColumns.h, rows in entry order)
    const std::vector<std::string> colNames = {"fZ", "fX", "fY", "trackerX", "trackerY", "detectorX", "detectorY", "trackerDX", "trackerDY", "detectorDX", "detectorDY", "reachedDetector", "reachedNextTracker", "fEvent", "TrackID"};
    auto colIdx = [&](const std::string& n){return int(std::find(colNames.begin(), colNames.end(), n) - colNames.begin());};
    const int iZ = colIdx("fZ"), iReachedDet = colIdx("reachedDetector"), iReachedTrk = colIdx("reachedNextTracker"), iEvent = colIdx("fEvent"), iTrack = colIdx("TrackID");
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
            fprintf(out, "%s, %s\n", ceLabel, name.Data());
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
            for (double z : cols[iZ]) hz->Fill(z);
            hz->Draw("");
            TString zTitle = Form("%s, %s, %s: %lld", p.label, ceLabel, name.Data(), nCE);
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
        hGrid->SetTitle(Form("Z Position of Charge Exchange, %s, %s: %lld;Z Position (mm);Counts", ceLabel, name.Data(), nCE));
        hGrid->SetStats(0); // small grid pads have no margin for a stats box, it would sit on the graph

        const double pitchT = 400.0 / ReadSettings().at("nCols_blocks"), pitchD = 400.0 / ReadSettings().at("nCols_detector"); // cell pitch, mm
        struct XYPlot { const char *file; TString expr; const char *cut; const char *label; };
        XYPlot xyPlots[] = {
            {"2_CE_xy",       "fY:fX>>h2(100,-50,50,100,-50,50)",                 "",                  "Position of Charge Exchange"},
            {"3_Tracker_xy",  "trackerY:trackerX>>h2(100,-50,50,100,-50,50)",    "reachedNextTracker==1", "Position at Next Tracker"},
            {"4_ZDC_xy", "detectorY:detectorX>>h2(100,-50,50,100,-50,50)",   "reachedDetector==1", "Position at ZDC"},
            // the detector's version of 3_ and 4_: block cell centers over the range of 25 cells centered on the beam axis, 100x100 bins
            {"5_Tracker_xy_detector", Form("trackerDY:trackerDX>>h2(100,%g,%g,100,%g,%g)", -12.5 * pitchT, 12.5 * pitchT, -12.5 * pitchT, 12.5 * pitchT), "reachedNextTracker==1", "Position at Next Tracker, detector block position"},
            {"6_ZDC_xy_detector", Form("detectorDY:detectorDX>>h2(100,%g,%g,100,%g,%g)", -12.5 * pitchD, 12.5 * pitchD, -12.5 * pitchD, 12.5 * pitchD), "reachedDetector==1", "Position at ZDC, detector block position"},
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
            TString xyTitle = Form("%s, %s, %s", p.label, ceLabel, name.Data());
            if (TString(p.file) == "4_ZDC_xy" || TString(p.file) == "6_ZDC_xy_detector") xyTitle += ", " + rate;
            if (TString(p.file) == "3_Tracker_xy" || TString(p.file) == "5_Tracker_xy_detector") xyTitle += ", " + trackerRate;
            hxy->SetTitle(PlotTitle(p.file, (xyTitle + ";X Position (mm);Y Position (mm)").Data()).c_str());
            canvas->Update();
            WritePlotNumbers(hxy, Form("%s/%s.pdf", dir.Data(), p.file), numberFiles);
            canvas->SaveAs(Form("%s/%s.pdf", dir.Data(), p.file));
        }
        for (size_t k = 0; k < colNames.size(); k++) // the recon grids need only these columns
            if ((int)k != iZ && (int)k != iReachedDet && (int)k != iEvent && (int)k != iTrack) std::vector<double>().swap(blockCols[b][k]);
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
        TTree *a = fAng ? (TTree*)fAng->Get(kAngleTree.c_str()) : nullptr; // the AnglesCE view: AnglesRaw rows with passQuality and inCE set
        if (!a || !a->GetBranch("inCE")){Error("GraphBlocks", "AnglesRaw with inCE missing in cuts/cut_angles_all.root (rerun ReconAngleAll.cpp)"); return;}
        const auto keys = Columns(a, {"fEvent", "ParentID"}, "passQuality != 0 && inCE != 0"); // multithreaded read (TreeColumns.h)
        for (size_t i = 0; i < keys[0].size(); i++) reconKeys.insert({(int)keys[0][i], (int)keys[1][i]});
    }
    struct ReconGrid {const char *file, *label; bool zdcOnly;};
    for (const ReconGrid &g : {ReconGrid{"CE_per_block_recon.pdf", "reconstructed CE", false}, ReconGrid{"CE_per_block_recon_zdc.pdf", "reconstructed CE that hit ZDC", true}}){
        TCanvas* reconCanvas = new TCanvas(Form("canvas_%s", g.file), g.label, 1500, 1200);
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
}
