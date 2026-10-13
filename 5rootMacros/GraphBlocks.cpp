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

void GraphBlocks(){
    const TString filename = "cut_blocks.root";
    const TString outDir = "pngs/Blocks";
    const int nBins = 30;

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

    // 2x2 style grid (more blocks make it grow)
    int nCols = std::ceil(std::sqrt((double)blocks.size()));
    int nRows = std::ceil(blocks.size() / (double)nCols);
    TCanvas* gridCanvas = new TCanvas("gridCanvas", "CE per block", 1500, 1200);
    gridCanvas->Divide(nCols, nRows);

    for (size_t b = 0; b < blocks.size(); b++) {
        TTree *t = blocks[b];
        TString name = t->GetName();
        TString dir = Form("%s/%s", outDir.Data(), name.Data());
        gSystem->mkdir(dir, true);

        double zLo = t->GetMinimum("fZ");
        double zHi = t->GetMaximum("fZ");
        double pad = 0.05 * (zHi - zLo) + 1.0;
        TString zHist = Form("fZ>>h(%d,%f,%f)", nBins, zLo - pad, zHi + pad);

        struct ZPlot { const char *file; const char *cut; const char *label; };
        ZPlot zPlots[] = {
            {"CE_z",          "",                     "CE"},
            {"NextTracker_z", "reachedNextTracker==1", "CE Protons Reaching Next Tracker"},
            {"Detector_z",    "reachedDetector==1",    "CE Protons Entering Detector"},
        };
        for (auto &p : zPlots) {
            canvas->cd();
            canvas->Clear();
            t->Draw(zHist, p.cut, "");
            t->GetHistogram()->SetTitle(Form("%s, %s: %lld;Z Position (mm);Counts", name.Data(), p.label, t->GetSelectedRows()));
            canvas->Update();
            canvas->SaveAs(Form("%s/%s.png", dir.Data(), p.file));
        }

        // same CE z plot into the grid
        gridCanvas->cd(b + 1);
        gPad->SetRightMargin(0.05);
        t->Draw(Form("fZ>>hGrid%zu(%d,%f,%f)", b + 1, nBins, zLo - pad, zHi + pad), "", "");
        t->GetHistogram()->SetTitle(Form("%s, CE: %lld;Z Position (mm);Counts", name.Data(), t->GetSelectedRows()));

        struct XYPlot { const char *file; const char *expr; const char *cut; const char *label; };
        XYPlot xyPlots[] = {
            {"CE_xy",       "fY:fX>>h2(100,-50,50,100,-50,50)",                 "",                  "Position of Charge Exchange"},
            {"Detector_xy", "detectorY:detectorX>>h2(100,-50,50,100,-50,50)",   "reachedDetector==1", "Position at Detector, CE Protons"},
        };
        for (auto &p : xyPlots) {
            canvas->cd();
            canvas->Clear();
            t->Draw(p.expr, p.cut, "colz");
            t->GetHistogram()->SetTitle(Form("%s, %s;X Position (mm);Y Position (mm)", name.Data(), p.label));
            canvas->Update();
            canvas->SaveAs(Form("%s/%s.png", dir.Data(), p.file));
        }
    }

    gridCanvas->Update();
    gridCanvas->SaveAs(Form("%s/CE_per_block.png", outDir.Data()));
    gSystem->Unlink("pngs/CE_per_block.png");
}
