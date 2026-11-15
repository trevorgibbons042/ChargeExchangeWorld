#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TLegend.h>
#include <TTree.h>
#include <TROOT.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <functional>
#include <set>
#include <TSystem.h>
#include <TStyle.h>

#include "BlockSquares.h"
#include "CutFiles.h"
#include "FitCore.h"
#include "Plot3DCols.h"
#include "TreeColumns.h"

// The maps (gP3Keys, gP3Flags, gP3Vals), isReco/rawVal/rawOk and the columns (cut_Plot3DCE.root, cut_Plot3DAllEv.root) or aliases (AnglesRaw) that replace them in the Draw cuts (P3) are in Plot3DCols.h.

// Plot3D writes 1_position (1_true and 2_block) and 2_energy. P3_SUB=n runs one group only (0: 1_position/1_true, 1: 1_position/2_block, 2: 2_energy), so the three can run as three processes side by side. Without it all three run.
static bool RunSub(int i){const char* e = std::getenv("P3_SUB"); return !e || !*e || std::atoi(e) == i;}

void F1_2_Plot3D(){
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    FastJit();
    std::string folder = "plots";

    // step folders only; CE and recon proton plots share a folder, file names end in _CE or _recon

    gSystem->mkdir(folder.c_str(), true);

    // the maps are loaded by CutPlots.cpp (zdc) (Plot3DColumns), not here


    // The three sets are rows of AnglesRaw (one recon fit per track the finder finds) split by evClass, with the Plot3D columns (Plot3DCols.h, CutPlots.cpp (zdc)):
    // events2 CE recon protons (evClass 0, cut_Plot3DCE.root), eventsAllEv CE + CE mix recon (evClass 0 or 1, cut_Plot3DAllEv.root), eventsAllPart all recon (cut_Plot3DAllPart.root).
    // The position of a row is the last tracker hit (lastX, lastY, lastXD, lastYD), its KE the ZDC sum of the event (finalKE, finalKED); every Draw cut goes through P3() so isReco(fEvent) etc. read the columns
    TTree *events2 = OpenCut("cuts/cut_Plot3DCE.root", kAngleTree);

    // the position plots (1_position) have no ZDC requirement: the rows of cut_Plot3DAllEv / cut_Plot3DAllPart; the energy plots (2_energy) read the same files and need a ZDC KE through finalKE >= 0 in their cut strings
    TTree* posAllEv = OpenCut("cuts/cut_Plot3DAllEv.root", kAngleTree);
    TTree* posAllPart = OpenCut("cuts/cut_Plot3DAllPart.root", kAngleTree);
    TTree* eventsAllEv = posAllEv; // CE + CE mix recon and all recon: the same rows as the position plots, the ZDC KE comes from the cut strings
    TTree* eventsAllPart = posAllPart;
    if (!events2 || !posAllEv || !posAllPart){Error("Plot3D", "missing cuts/cut_Plot3DCE.root, cut_Plot3DAllEv.root or cut_Plot3DAllPart.root (run CutPlots.cpp first)"); return;}


    {

        TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
        gPad->SetRightMargin(0.30);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); // room for the rotated Y title and bottom tick labels
        gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X"); // stats/fit box to the right, fewer x labels

        // Every plot goes through this: draw expr from t with cut, set the title (title;xaxis;yaxis), let tweak adjust the
        // histogram (range, fit, ...), save to path. opt is the Draw option, e.g. "colz".
        std::set<std::string> fitFilesStarted; // stats and fit files already truncated in this run
        std::string zdcNote = " with ZDC KE"; // the legend phrase of the plots that need a ZDC KE (empty for the position plots)
        std::string titleSuffix; // added to the title text (before the axis titles), the KE cut of the position plots
        double nCERecoAll = -1, nCEZdcAll = -1; // CE protons at the ZDC before any cut (with a reco track, and all), counted on the first plot that shows them
        std::function<void()> countZdcAll = []{}; // replaced below, once the cut strings exist
        // The plots are made in 2 passes over the same code (plotAll below). Pass 1 (collect) only books the histograms and counts each plot reads. Then one event loop per tree fills the HistsMT / CountsMT memos (TreeColumns.h),
        // and pass 2 draws every plot from the memo. Without it each plot, overlay curve and percent box is an event loop of its own.
        bool collect = false;
        std::map<TTree*, std::vector<std::pair<std::string, std::string>>> bookedHists; // tree -> (drawExpr, cut) as HistMT gets them
        std::vector<std::string> bookedCounts; // cuts counted on events2 (the percent boxes and countZdcAll)
        std::map<TTree*, std::vector<std::pair<std::string, std::string>>> bookedFits; // tree -> (expr, cut) of the plots with an unbinned fit
        std::map<std::string, std::vector<double>> fitValues; // key expr + cut -> the values the fit reads, filled after pass 1 from one read per tree
        std::function<void()> bookZdcCounts = []{}; // replaced below: books the two cuts countZdcAll counts
        TH1D dummyH("dummyH", "", 100, -1000, 1000); dummyH.SetDirectory(nullptr); // handed to the tweaks in pass 1, which draw nothing then
        auto drawPlot = [&](TTree* t, const std::string& expr, const std::string& cut, const std::string& title,
                            const std::string& path, const char* opt = "", const std::function<void(TH1*)>& tweak = nullptr, int fitMode = kFitNone){ // fitMode: the unbinned gaussian fit of the plotted values (FitCore.h), before the tweak
            gSystem->mkdir(gSystem->GetDirName(path.c_str()), true);
            if (collect){ // pass 1: book what this plot reads (the histogram, the CE count of its cut, the ZDC counts), run the tweak on dummyH so its overlays book theirs, draw nothing
                bookedHists[t].emplace_back(P3(expr), P3(cut));
                bookedCounts.push_back(P3(cut)); bookZdcCounts();
                if (fitMode != kFitNone) bookedFits[t].emplace_back(P3(expr), P3(cut));
                if (tweak) tweak(&dummyH);
                return;
            }
            canvas->Clear();
            TH1* h = HistMT(t, P3(expr), P3(cut)); // all the threads (an answer from the memo after pass 1), then drawn like TTree::Draw did
            if (!h) return;
            h->Draw(opt);
            const size_t axes = title.find(';');
            h->SetTitle(PlotTitle(path, title.substr(0, axes) + titleSuffix + (axes == std::string::npos ? "" : title.substr(axes))).c_str());
            const bool ceOnly1D = t == events2 && h->GetDimension() == 1;
            const bool allPart1D = t == eventsAllPart && h->GetDimension() == 1;
            if (ceOnly1D) h->SetLineColor(CbCE()); // CE protons are red in every plot; the overlay tweaks below set their own colors
            if (fitMode != kFitNone){ // the same values and cut as the histogram, read once more for the likelihood (no bins)
                const auto booked = fitValues.find(P3(expr) + "\x1f" + P3(cut)); // read in pass 1's flush; a fit that was not booked reads its own
                const std::vector<double> xs = booked != fitValues.end() ? booked->second : Columns(t, SplitDraw(P3(expr)).exprs, P3(cut))[0];
                if (fitMode == kFitCore) FitCore(h, xs); else FitBlocks(h, xs);
            }
            if (allPart1D) h->SetLineColor(CbAll()); // all recon charged particles are green, as in the x overlays
            if (tweak) tweak(h);
            if (ceOnly1D || allPart1D) for (TObject* o : *h->GetListOfFunctions()) if (TF1* fit = dynamic_cast<TF1*>(o)) fit->SetLineColor(kBlack); // a red histogram would hide a red fit
            const std::string hName = h->GetName(); // the stats box shows the name: a readable one until the pdf is saved
            h->SetName(StatsTitle(h, title.substr(0, std::min(title.find(','), axes)) + (t == events2 ? ", CE recon protons" : t == eventsAllEv ? ", CE + CE mix recon" : t == eventsAllPart ? ", all recon" : "")).c_str()); // quantity: the title up to its first comma
            if (t == eventsAllEv || t == eventsAllPart || t == events2) countZdcAll(); // the ZDC totals the percent boxes use, counted once
            if (t == eventsAllEv || t == eventsAllPart){ // all recon protons (or particles): the CE ones (events2) under the same cut
                canvas->cd();
                DrawCEPercent(CountMT(events2, P3(cut)), h->GetEntries(), nCERecoAll); // lost: vs every reco CE proton at the ZDC
            } else if (t == events2 && nCERecoAll > 0){ // CE protons only: the share the cuts lose, vs every CE proton at the ZDC (reco or not, as the cut asks)
                const double nCE = CountMT(events2, P3(cut));
                const double nBase = cut.find("isReco") != std::string::npos ? nCERecoAll : nCEZdcAll;
                canvas->cd();
                if (nCE < nBase) DrawCEPercent(nCE, -1, nBase);
            }
            WritePlotNumbers(h, path, fitFilesStarted); // stats.csv, plus fit_values.csv if h has a fit, next to the plot
            canvas->Update();
            canvas->SaveAs(path.c_str());
            if (gSystem->AccessPathName(path.c_str())){ // the pdf is not there: TPDF::Open failed once on LCRC (2_block/2_xy_recon.pdf, 46 processes at the same time), so the folder is made again and the pdf saved a second time
                gSystem->mkdir(gSystem->GetDirName(path.c_str()), true);
                canvas->SaveAs(path.c_str());
            }
            h->SetName(hName.c_str());
        };
        // x:y position maps: linear color scale from 0 up to the fullest bin (empty bins stay blank)
        auto drawXY = [&](TTree* t, const std::string& expr, const std::string& cut, const std::string& title, const std::string& path){
            drawPlot(t, expr, cut, title, path, "colz", [](TH1* h){h->SetMinimum(0);});
        };
        // block position map, one bin per tracker cell (the histogram keeps the stats box, colorbar and percent boxes): each occupied cell is a square of the cell size (the pitch of the bin, edge cells may be cut by the frame) at its exact position, colored by
        // its count on the same log scale. The squares touch, no gaps.
        auto drawXYDots = [&](TTree* t, const std::string& expr, const std::string& cut, const std::string& title, const std::string& path, double halfRange = 0){ // halfRange > 0: hard frame of +-halfRange mm, the cells at the edge are clipped by it
            drawPlot(t, expr, cut, title, path, "colz", [halfRange, &collect](TH1* h){if (!collect) DrawBlockSquares(h, halfRange, false);}); // BlockSquares.h, linear color scale
        };

        // 1D plot on a hard frame of +-half mm: the bars of the cells past it are clipped by the frame (like 4_xy_zoom_recon.pdf); call in drawPlot's tweak, after the fit
        auto clipX = [&](TH1* h, double half){
            if (collect) return;
            canvas->Clear();
            TH1* fr = gPad->DrawFrame(-half, 0, half, 1.1 * h->GetMaximum());
            fr->SetTitle(h->GetTitle()); fr->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle()); fr->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
            h->Draw("sames");
        };

        // call in drawPlot's tweak: h (all recon protons, blue) with the CE protons (events2, red) and all recon charged particles (eventsAllPart, green) of the same quantity drawn over it, and a legend below the stats box.
        // logY: log y from 0.5, else a linear y 20% above the highest bin.
        auto overlayCEPart = [&](TH1* h, const std::string& exprCE, const std::string& cutCE, const std::string& exprPart, const std::string& cutPart, const char* ceLabel, bool logY){
            if (collect){bookedHists[events2].emplace_back(P3(exprCE), P3(cutCE)); bookedHists[eventsAllPart].emplace_back(P3(exprPart), P3(cutPart)); return;}
            TH1* hCE = HistMT(events2, P3(exprCE), P3(cutCE)); hCE->Draw("same");
            TH1* hPart = HistMT(eventsAllPart, P3(exprPart), P3(cutPart)); hPart->Draw("same");
            h->SetLineColor(CbProton()); hCE->SetLineColor(CbCE()); hPart->SetLineColor(CbAll()); hCE->SetStats(0); hPart->SetStats(0);
            const double top = std::max({h->GetMaximum(), hCE->GetMaximum(), hPart->GetMaximum()});
            if (logY){h->SetMinimum(0.5); h->SetMaximum(2.0 * top); canvas->SetLogy(1);} else h->SetMaximum(1.2 * top);
            TLegend* leg = new TLegend(0.745, 0.30, 0.99, 0.42); leg->SetTextSize(0.012); // below the stats box
            leg->AddEntry(hPart, ("all recon" + zdcNote).c_str(), "l");
            leg->AddEntry(h, ("CE + CE mix recon" + zdcNote).c_str(), "l");
            leg->AddEntry(hCE, ceLabel, "l");
            leg->Draw();
        };

        std::string zdcProtons = "finalKE >= 0"; // the event has a ZDC KE (NaN fails)
        // position plots: only protons with a reconstructed track
        std::string zdcReco = zdcProtons + " && isReco(fEvent)";
        // the same for every recorded particle (any pdg, eventsAllPart): ZDC hit and a recon track
        std::string zdcRecoPartAll = "finalKE >= 0 && isReco(fEvent)";
        countZdcAll = [&]{
            if (nCERecoAll >= 0) return;
            {const auto c = CountsMT(events2, {P3(zdcReco), P3(zdcProtons)}); nCERecoAll = c[0]; nCEZdcAll = c[1];} // one event loop for both
        };

        bookZdcCounts = [&]{bookedCounts.push_back(P3(zdcReco)); bookedCounts.push_back(P3(zdcProtons));};

        auto plotAll = [&]{ // every plot below, run twice (collect false / true, see bookedHists)
        // Position plots have no KE cut and no ZDC requirement: 1_position/1_true/ and 1_position/2_block/
        {
            // the CE particles at the last tracker, whether or not they reach the ZDC
            zdcProtons = "true"; zdcReco = zdcRecoPartAll = "isReco(fEvent) == 1"; zdcNote = "";
            nCERecoAll = nCEZdcAll = -1; // the percent boxes count against these cuts
            const std::string folder_pos = folder + "/1_position/1_true/", folder_posD = folder + "/1_position/2_block/";
            const std::string xySlice = zdcReco + " && fabs(lastY) < 2"; // the block position plots (D) have the same cut as the true position ones
            titleSuffix = ", has recon track, no ZDC or KE cut"; // the cut of the position plots: a recon track
            const std::string xySliceP = zdcRecoPartAll + " && fabs(lastY) < 2";
            if (RunSub(0)){ // P3_SUB 0: 1_true
            // ZDC position, true (lastX, lastY)
            drawXY(eventsAllEv, "lastY:lastX>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at last tracker, CE + CE mix recon;X Position (mm);Y Position (mm)", folder_pos + "2_xy_recon.pdf");
            drawXY(eventsAllEv, "lastY:lastX>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at last tracker, CE + CE mix recon;X Position (mm);Y Position (mm)", folder_pos + "4_xy_zoom_recon.pdf");
            drawXY(events2, "lastY:lastX>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at last tracker, CE recon protons;X Position (mm);Y Position (mm)", folder_pos + "1_xy_CE.pdf");
            drawXY(events2, "lastY:lastX>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at last tracker, CE recon protons;X Position (mm);Y Position (mm)", folder_pos + "3_xy_zoom_CE.pdf");

            gStyle->SetOptFit(11); // fit parameters with errors (the fits are unbinned: no chi2 or probability)
            drawPlot(events2, "lastX>>hSliceX200(100,-200,200)", xySlice, "X Position at last tracker, CE recon protons, |Y| < 2 mm;X Position (mm);Count", folder_pos + "5_x_fit_200_CE.pdf", "", nullptr, kFitCore);
            drawPlot(events2, "lastX>>hSliceX(100,-50,50)", xySlice, "X Position at last tracker, CE recon protons, |Y| < 2 mm;X Position (mm);Count", folder_pos + "7_x_fit_50_CE.pdf", "", nullptr, kFitCore);

            // the x slices as one overlay each, no fit: CE recon protons (red), all recon protons (blue), all recon charged particles (green)
            struct XOv {std::string bins, file, tag;};
            for (const auto& o : std::vector<XOv>{{"(100,-200,200)", "9_x_overlay_200.pdf", "200"}, {"(100,-50,50)", "10_x_overlay_50.pdf", "50"}}){
                drawPlot(eventsAllEv, "lastX>>hXOv" + o.tag + o.bins, xySlice, "X Position at last tracker, all recon, CE + CE mix recon and CE recon protons, |Y| < 2 mm;X Position (mm);Count", folder_pos + o.file, "", [&](TH1* h){
                    overlayCEPart(h, "lastX>>hXOvCE" + o.tag + o.bins, xySlice, "lastX>>hXOvPart" + o.tag + o.bins, xySliceP, "CE recon protons", false);
                });
            }
            // r = sqrt(lastX^2 + lastY^2) of the ZDC hit (the xy plots as a radius), one overlay per range like the xy plots: CE recon protons (red), all recon protons (blue), all recon charged particles (green)
            for (const auto& o : std::vector<std::array<std::string, 3>>{{"200", "11_r_overlay_200.pdf", "(100,0,200)"}, {"50", "12_r_overlay_50.pdf", "(100,0,50)"}}){
                const std::string rE = "sqrt(lastX*lastX+lastY*lastY)>>hROv";
                drawPlot(eventsAllEv, rE + o[0] + o[2], zdcReco, "Radius at last tracker, all recon, CE + CE mix recon and CE recon protons;r = #sqrt{X^{2} + Y^{2}} (mm);Count", folder_pos + o[1], "", [&](TH1* h){
                    overlayCEPart(h, "sqrt(lastX*lastX+lastY*lastY)>>hROvCE" + o[0] + o[2], zdcReco, "sqrt(lastX*lastX+lastY*lastY)>>hROvPart" + o[0] + o[2], zdcRecoPartAll, "CE recon protons", false);
                });
            }

            }
            if (RunSub(1)){ // P3_SUB 1: 2_block
            // ZDC position, block (lastXD, lastYD): the last tracker plane is nCols_blocks x nRows_blocks blocks over 400 mm. Every bin holds the same whole number of blocks (BlockBins, FitCore.h):
            // 100 bins over 100 mm on a pitch of 400 / 801 mm put an extra block in the bin that starts at x = 0 (a cross at x = 0 and y = 0), and in one bin of ten along the x slices (spikes).
            const int nColsD = (int)ReadSettings().at("nCols_blocks"), nRowsD = (int)ReadSettings().at("nRows_blocks");
            const std::string binsSlice200 = BlockBins(nColsD, 400.0 / nColsD, 200), binsSlice50 = BlockBins(nColsD, 400.0 / nColsD, 53.5); // x slices: one bin per tracker cell, the +-50 mm ones with whole cells past 50 (clipped by the frame, clipX)
            // the 4 xy plots: one marker per occupied tracker cell (nCols_blocks x nRows_blocks over 400 mm), at the exact cell position (drawXYDots)
            const std::string binsCell = BlockBins(nColsD, 400.0 / nColsD, 200), binsCellY = BlockBins(nRowsD, 400.0 / nRowsD, 200), binsCellZoom = BlockBins(nColsD, 400.0 / nColsD, 53.5), binsCellZoomY = BlockBins(nRowsD, 400.0 / nRowsD, 53.5); // whole cells past +-50: the frame is +-50, the edge cells are clipped
            const std::string cellNote = Form(", 1 square per tracker cell (%d x %d over 400 mm)", nColsD, nRowsD), zoomNote = cellNote + ", edge cells clipped at 50 mm";
            drawXYDots(eventsAllEv, "lastYD:lastXD>>h2(" + binsCell + "," + binsCellY + ")", zdcReco, "Position at last tracker BLOCK, CE + CE mix recon" + cellNote + ";X Position (mm);Y Position (mm)", folder_posD + "2_xy_recon.pdf");
            drawXYDots(eventsAllEv, "lastYD:lastXD>>h2(" + binsCellZoom + "," + binsCellZoomY + ")", zdcReco, "Position at last tracker BLOCK, CE + CE mix recon" + zoomNote + ";X Position (mm);Y Position (mm)", folder_posD + "4_xy_zoom_recon.pdf", 50);
            drawXYDots(events2, "lastYD:lastXD>>h2(" + binsCell + "," + binsCellY + ")", zdcReco, "Position at last tracker BLOCK, CE recon protons" + cellNote + ";X Position (mm);Y Position (mm)", folder_posD + "1_xy_CE.pdf");
            drawXYDots(events2, "lastYD:lastXD>>h2(" + binsCellZoom + "," + binsCellZoomY + ")", zdcReco, "Position at last tracker BLOCK, CE recon protons" + zoomNote + ";X Position (mm);Y Position (mm)", folder_posD + "3_xy_zoom_CE.pdf", 50);

            // the six x fits: ZDC hit and a recon track only
            const std::string xySliceB = zdcReco + " && fabs(lastYD) < 2", xySliceBP = zdcRecoPartAll + " && fabs(lastYD) < 2";
            drawPlot(events2, "lastXD>>hSliceDX200(" + binsSlice200 + ")", xySliceB, "X Position at last tracker BLOCK, CE recon protons, |Y| < 2 mm;X Position (mm);Count", folder_posD + "5_x_fit_200_CE.pdf", "", nullptr, kFitBlocks); //whole blocks per bin
            drawPlot(events2, "lastXD>>hSliceDX(" + binsSlice50 + ")", xySliceB, "X Position at last tracker BLOCK, CE recon protons, |Y| < 2 mm, edge cells clipped at 50 mm;X Position (mm);Count", folder_posD + "7_x_fit_50_CE.pdf", "", [&](TH1* h){clipX(h, 50);}, kFitBlocks);
            // the x slices as one overlay each, no fit: CE recon protons (red), all recon protons (blue), all recon charged particles (green)
            for (const auto& o : std::vector<std::array<std::string, 4>>{{"200", binsSlice200, "6_x_overlay_200.pdf", ""}, {"50", binsSlice50, "8_x_overlay_50.pdf", ", edge cells clipped at 50 mm"}}){
                const std::string tag = o[0], bins = "(" + o[1] + ")";
                drawPlot(eventsAllEv, "lastXD>>hXOvD" + tag + bins, xySliceB, "X Position at last tracker BLOCK, all recon, CE + CE mix recon and CE recon protons, |Y| < 2 mm" + o[3] + ";X Position (mm);Count", folder_posD + o[2], "", [&](TH1* h){
                    overlayCEPart(h, "lastXD>>hXOvDCE" + tag + bins, xySliceB, "lastXD>>hXOvDPart" + tag + bins, xySliceBP, "CE recon protons", false);
                    if (tag == "50") h->GetXaxis()->SetRangeUser(-50, 50);
                });
            }

            // r = sqrt(lastXD^2 + lastYD^2) of the last tracker block (the block xy plots as a radius), same overlay as 1_true
            for (const auto& o : std::vector<std::array<std::string, 3>>{{"200", "13_r_overlay_200.pdf", "(100,0,200)"}, {"50", "14_r_overlay_50.pdf", "(100,0,50)"}}){
                drawPlot(eventsAllEv, "sqrt(lastXD*lastXD+lastYD*lastYD)>>hROvD" + o[0] + o[2], zdcReco, "Radius at last tracker BLOCK, all recon, CE + CE mix recon and CE recon protons;r = #sqrt{X^{2} + Y^{2}} (mm);Count", folder_posD + o[1], "", [&](TH1* h){
                    overlayCEPart(h, "sqrt(lastXD*lastXD+lastYD*lastYD)>>hROvDCE" + o[0] + o[2], zdcReco, "sqrt(lastXD*lastXD+lastYD*lastYD)>>hROvDPart" + o[0] + o[2], zdcRecoPartAll, "CE recon protons", false);
                });
            }

            // RECONPHIBL_ESTNEU (rawVal 15, 0 to 360 deg) of the same |Y| < 2 mm events as the x fits above, no fit
            const std::string phiOk = " && rawOk(fEvent, 15)";
            // one version per x window, like the x fits: |X| < 200 mm (file 9) and |X| < 50 mm (file 10)
            struct PhiWin {int w; std::string nCE;};
            for (const auto& pw : std::vector<PhiWin>{{200, "9"}, {50, "10"}}){
                const std::string tag = std::to_string(pw.w), xWin = Form(" && fabs(lastXD) < %d", pw.w), winNote = Form(", |X| < %d mm", pw.w);
                drawPlot(events2, "rawVal(fEvent, 15)>>hPhiDY2CE" + tag + "(100,0,360)", xySliceB + phiOk + xWin, "RECONPHIBL_ESTNEU, CE recon protons, |Y| < 2 mm" + winNote + ";RECONPHIBL_ESTNEU (degrees);Count", folder_posD + pw.nCE + "_phiBlock_" + tag + "_CE.pdf");
            }
            // RECONTHETABL_ESTNEU (rawVal 11, deg) of the same events, the theta version of the phi plots above (files 11, 12), same |X| windows
            const std::string thOk = " && rawOk(fEvent, 11)", thBins = Form("(100,0,%g)", ThetaAxisMax()), thE = "rawVal(fEvent, 11)>>";
            struct ThWin {int w; std::string nCE;};
            for (const auto& tw : std::vector<ThWin>{{200, "11"}, {50, "12"}}){
                const std::string tag = std::to_string(tw.w), xWin = Form(" && fabs(lastXD) < %d", tw.w), winNote = Form(", |X| < %d mm", tw.w);
                drawPlot(events2, thE + "hThDY2CE" + tag + thBins, xySliceB + thOk + xWin, "RECONTHETABL_ESTNEU, CE recon protons, |Y| < 2 mm" + winNote + ";RECONTHETABL_ESTNEU (degrees);Count", folder_posD + tw.nCE + "_thetaBlock_" + tag + "_CE.pdf");
            }
            }
            titleSuffix.clear();
        }

        // KE x axes: 0 to KeAxisMax (FitCore.h, 30k MeV at 25 GeV/c) for the 2_energy plots; the tag names the range in k MeV
        struct KERange {double lo, hi; std::string tag;};
        const double keAxisMax = KeAxisMax();
        auto keTag = [](double lo, double hi){return std::string(Form("_%gto%gk", lo / 1000, hi / 1000));};
        const std::vector<KERange> energyRanges = {{0, keAxisMax, keTag(0, keAxisMax)}}; // 2_energy: 0 to 1.2 x the beam KE (KeAxisMax), full range only
        auto binSpec = [](const KERange& r){return std::string(Form("(100,%g,%g)", r.lo, r.hi));};

        // 2_energy/ holds the energy plots
        if (RunSub(2)) {
            // same rows as the position plots, the ZDC KE comes from the cut strings below
            zdcProtons = "finalKE >= 0"; zdcReco = zdcProtons + " && isReco(fEvent)"; zdcRecoPartAll = "finalKE >= 0 && isReco(fEvent)"; zdcNote = " with ZDC KE";
            nCERecoAll = nCEZdcAll = -1;
            const std::string folder_smear = folder + "/2_energy/";
            gSystem->mkdir(folder_smear.c_str(), true);
            titleSuffix = ", no KE cut";
            // Energy plots: 2_energy/<n>_<quantity><range>_<CE or recon>.pdf, numbered 1 (true KE), 2 (smeared KE), 3 (total KE, charged + neutral), 5, 7, 9 and 10 (smearing)
            // True KE (finalKE) and smeared KE (finalKED) at the ZDC.
            struct KEQty {std::string branch, file, axis, title;};
            std::vector<KEQty> keQtys = {{"finalKE", "true_ke", "ZDC KE sum of the event (MeV)", "True ZDC KE (sum of the charged particles of the event)"},
                                       {"finalKED", "smeared_ke", "Smeared ZDC KE sum of the event (MeV)", "Smeared ZDC KE (sum of the charged particles of the event)"}};
            // one log y plot per quantity: CE protons (red), all recon protons (blue) and all recon charged particles (green): 1_true_ke, 2_smeared_ke
            for (const auto& q : keQtys){
                for (const auto& r : energyRanges){
                    drawPlot(eventsAllEv, q.branch + ">>hKE" + r.tag.substr(1) + binSpec(r), zdcReco, q.title + ", all recon, CE + CE mix recon and CE recon protons with ZDC KE;" + q.axis + ";Count",
                             folder_smear + std::to_string(&q - &keQtys[0] + 1) + "_" + q.file + r.tag + ".pdf", "", [&](TH1* h){
                        overlayCEPart(h, q.branch + ">>hKECE" + r.tag.substr(1) + binSpec(r), zdcProtons, q.branch + ">>hKEPart" + r.tag.substr(1) + binSpec(r), zdcRecoPartAll, "CE recon protons with ZDC KE", true); // every charged particle: contains the protons
                    });
                    canvas->SetLogy(0);
                }
            }

            // 3_total_ke: the true KE of every particle at the ZDC, charged and neutral (zdcKEAll, EventSummary of the sim). The cut is zdcKEAll >= 0 (an event with only neutrons at the ZDC counts), so the percent boxes are counted again against it
            {
                const std::string keepProtons = zdcProtons, keepReco = zdcReco, keepPartAll = zdcRecoPartAll, keepNote = zdcNote;
                zdcProtons = "zdcKEAll >= 0"; zdcReco = zdcProtons + " && isReco(fEvent)"; zdcRecoPartAll = zdcReco; zdcNote = " with total ZDC KE";
                nCERecoAll = nCEZdcAll = -1;
                for (const auto& r : energyRanges){
                    drawPlot(eventsAllEv, "zdcKEAll>>hKETot" + r.tag.substr(1) + binSpec(r), zdcReco, "Total ZDC KE (sum of all particles of the event, charged and neutral), all recon, CE + CE mix recon and CE recon protons with total ZDC KE;Total ZDC KE sum of the event (MeV);Count",
                             folder_smear + "3_total_ke" + r.tag + ".pdf", "", [&](TH1* h){
                        overlayCEPart(h, "zdcKEAll>>hKETotCE" + r.tag.substr(1) + binSpec(r), zdcProtons, "zdcKEAll>>hKETotPart" + r.tag.substr(1) + binSpec(r), zdcRecoPartAll, "CE recon protons with total ZDC KE", true);
                    });
                    canvas->SetLogy(0);
                }
                zdcProtons = keepProtons; zdcReco = keepReco; zdcRecoPartAll = keepPartAll; zdcNote = keepNote;
                nCERecoAll = nCEZdcAll = -1;
            }

            // ZDC energy smearing (finalKED = Gaussian smear of finalKE, /zdc/resStochastic and /zdc/resConstant in the sim)
            {

                // finalKE = true KE, finalKED = smeared KE, dE = finalKED - finalKE
                // the smear sigma is resStochastic * sqrt(KE in GeV) GeV, so the dE axis is +-10000 MeV at the 25 GeV/c beam (24078 MeV KE) and scales with sqrt(KE)
                const double smearLim = 10000.0 * std::sqrt(kBeamKE / 24078.0);

                drawPlot(events2, Form("finalKED-finalKE>>hSmearDiffCE(100,-%g,%g)", smearLim, smearLim), zdcProtons, "Energy Difference at ZDC, CE recon protons with ZDC KE;finalKED - finalKE (MeV);Count", folder_smear + "5_diff_CE.pdf", "", nullptr, kFitCore);

                // percent error = 100 * (finalKED - finalKE) / finalKE
                const std::string pct = "100*((finalKED-finalKE)/finalKE)"; // bracketed: the C++ of HistMT would divide 100*(a-b) by c, TTree::Draw 100*(a-b/c)
                drawPlot(events2, pct + ">>hSmearPctCE(100,-100,100)", zdcProtons, "Percent Error at ZDC, CE recon protons with ZDC KE;Percent Error (%);Count", folder_smear + "7_pct_CE.pdf", "", nullptr, kFitCore);
                

                // the smearing quantities as one log y overlay each: CE protons (red), all recon protons (blue), all recon charged particles (green)
                struct SmearOv {std::string expr, bins, file, axis, title, tag;};
                const std::vector<SmearOv> smearOvs = {{"finalKED-finalKE", Form("(100,-%g,%g)", smearLim, smearLim), "9_diff_overlay.pdf", "finalKED - finalKE (MeV)", "Energy Difference at ZDC", "Diff"},
                                                       {pct, "(100,-100,100)", "10_pct_overlay.pdf", "Percent Error (%)", "Percent Error at ZDC", "Pct"}};
                for (const auto& o : smearOvs){
                    drawPlot(eventsAllEv, o.expr + ">>hOv" + o.tag + o.bins, zdcReco, o.title + ", all recon, CE + CE mix recon and CE recon protons with ZDC KE;" + o.axis + ";Count",
                             folder_smear + o.file, "", [&](TH1* h){
                        overlayCEPart(h, o.expr + ">>hOvCE" + o.tag + o.bins, zdcProtons, o.expr + ">>hOvPart" + o.tag + o.bins, zdcRecoPartAll, "CE recon protons with ZDC KE", true);
                    });
                    canvas->SetLogy(0);
                }
            }
            titleSuffix.clear();
        }

        };

        collect = true; plotAll(); collect = false; // pass 1: book
        for (auto& [tree, reqs] : bookedHists){ // one event loop per tree; a histogram name twice in one HistsMT call would replace the first (MakeHist), so such a request goes to a later round
            std::sort(reqs.begin(), reqs.end()); reqs.erase(std::unique(reqs.begin(), reqs.end()), reqs.end());
            std::vector<std::vector<std::pair<std::string, std::string>>> rounds;
            std::vector<std::set<std::string>> roundNames;
            for (const auto& rq : reqs){
                const std::string target = SplitDraw(rq.first).target, name = target.substr(0, target.find('('));
                size_t k = 0;
                while (k < rounds.size() && roundNames[k].count(name)) k++;
                if (k == rounds.size()){rounds.emplace_back(); roundNames.emplace_back();}
                rounds[k].push_back(rq); roundNames[k].insert(name);
            }
            for (const auto& round : rounds) for (TH1* h : HistsMT(tree, round)) delete h; // the memo keeps its own clone
        }
        std::sort(bookedCounts.begin(), bookedCounts.end()); bookedCounts.erase(std::unique(bookedCounts.begin(), bookedCounts.end()), bookedCounts.end());
        CountsMT(events2, bookedCounts);
        for (auto& [tree, reqs] : bookedFits){ // the fit values of every fit plot of a tree in one read: the value and a 0/1 flag of the cut per fit, the rows with flag 1 in entry order are what Columns(tree, expr, cut) gives
            std::sort(reqs.begin(), reqs.end()); reqs.erase(std::unique(reqs.begin(), reqs.end()), reqs.end());
            std::vector<std::string> exprs;
            for (const auto& rq : reqs){exprs.push_back(SplitDraw(rq.first).exprs[0]); exprs.push_back("((" + (rq.second.empty() ? std::string("1") : rq.second) + ") ? 1.0 : 0.0)");}
            const auto cols = Columns(tree, exprs);
            for (size_t j = 0; j < reqs.size(); j++){
                std::vector<double>& xs = fitValues[reqs[j].first + "\x1f" + reqs[j].second];
                for (size_t i = 0; i < cols[2 * j].size(); i++) if (cols[2 * j + 1][i] > 0.5) xs.push_back(cols[2 * j][i]);
            }
        }
        plotAll(); // pass 2: draw, every HistMT and CountMT answered from the memo
    }
    BoxAndCropPlots();
}
