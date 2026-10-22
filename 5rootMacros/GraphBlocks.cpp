#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TString.h>
#include <TSystem.h>
#include <TList.h>

#include <cmath>
#include <vector>

#include "FitCore.h"

void GraphBlocks(){
    const TString filename = "cut_blocks.root";
    const TString outDir = "plots/5_CEblock";
    const int nBins = 100;
    const char *ceLabel = "CE events CE protons";

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

    std::set<std::string> numberFiles; // stats csv files already truncated in this run
    for (size_t b = 0; b < blocks.size(); b++) {
        TTree *t = blocks[b];
        TString name = t->GetName();
        TString dir = Form("%s/%s", outDir.Data(), name.Data());
        gSystem->mkdir(dir, true);

        double zLo = t->GetMinimum("fZ");
        double zHi = t->GetMaximum("fZ");
        double pad = 0.05 * (zHi - zLo) + 1.0;
        Long64_t nCE = t->GetEntries();
        Long64_t nZDC = t->GetEntries("reachedDetector==1");
        TString rate = Form("ZDC/CE = %lld/%lld = %.4f", nZDC, nCE, nCE ? (double)nZDC / nCE : 0.0);
        Long64_t nTracker = t->GetEntries("reachedNextTracker==1");
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
            t->Draw(zHist, p.cut, "");
            TString zTitle = Form("%s, %s, %s: %lld", p.label, ceLabel, name.Data(), t->GetSelectedRows());
            t->GetHistogram()->SetTitle(PlotTitle(p.file, (zTitle + ";Z Position (mm);Counts").Data()).c_str());
            canvas->Update();
            WritePlotNumbers(t->GetHistogram(), Form("%s/%s.pdf", dir.Data(), p.file), numberFiles);
            canvas->SaveAs(Form("%s/%s.pdf", dir.Data(), p.file));
        }

        // same CE z plot into the grid
        gridCanvas->cd(b + 1);
        gPad->SetRightMargin(0.05);
        t->Draw(Form("fZ>>hGrid%zu(%d,%f,%f)", b + 1, nBins, zLo - pad, zHi + pad), "", "");
        t->GetHistogram()->SetTitle(Form("Z Position of Charge Exchange, %s, %s: %lld;Z Position (mm);Counts", ceLabel, name.Data(), t->GetSelectedRows()));

        struct XYPlot { const char *file; const char *expr; const char *cut; const char *label; };
        XYPlot xyPlots[] = {
            {"2_CE_xy",       "fY:fX>>h2(100,-50,50,100,-50,50)",                 "",                  "Position of Charge Exchange"},
            {"3_Tracker_xy",  "trackerY:trackerX>>h2(100,-50,50,100,-50,50)",    "reachedNextTracker==1", "Position at Next Tracker"},
            {"4_ZDC_xy", "detectorY:detectorX>>h2(100,-50,50,100,-50,50)",   "reachedDetector==1", "Position at ZDC"},
        };
        for (auto &p : xyPlots) {
            canvas->cd();
            canvas->Clear();
            t->Draw(p.expr, p.cut, "colz");
            TString xyTitle = Form("%s, %s, %s", p.label, ceLabel, name.Data());
            if (TString(p.file) == "4_ZDC_xy") xyTitle += ", " + rate;
            if (TString(p.file) == "3_Tracker_xy") xyTitle += ", " + trackerRate;
            t->GetHistogram()->SetTitle(PlotTitle(p.file, (xyTitle + ";X Position (mm);Y Position (mm)").Data()).c_str());
            canvas->Update();
            WritePlotNumbers(t->GetHistogram(), Form("%s/%s.pdf", dir.Data(), p.file), numberFiles);
            canvas->SaveAs(Form("%s/%s.pdf", dir.Data(), p.file));
        }
    }

    gridCanvas->Update();
    gridCanvas->SaveAs(Form("%s/CE_per_block.pdf", outDir.Data()));
    std::set<std::string> gridStats; // one row per block: the CE count in each panel's title
    for (TTree* t : blocks) WriteStatsRow(Form("%s/CE_per_block.pdf", outDir.Data()), "plot,block,ce_events", Form("CE_per_block.pdf,%s,%lld", t->GetName(), (Long64_t)t->GetEntries()), gridStats);
}
