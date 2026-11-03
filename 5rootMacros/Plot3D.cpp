#include <TColor.h>
#include <TH2D.h>
#include <TMarker.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraph2D.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TROOT.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <functional>
#include <memory>
#include <set>
#include <TSystem.h>
#include <sstream>
#include <chrono>
#include <unistd.h> // getpid
#include <TStyle.h>

#include "BlockSquares.h"
#include "CutOpen.h"
#include "FitCore.h"
#include "Plot3DCols.h"
#include "QualityCut.h"
#include "TreeColumns.h"

// The maps (gP3Keys, gP3Flags, gP3Vals), isReco/rawVal/rawOk and the columns (cut_ZdcCE.root, cut_ZdcAllEv.root) or aliases (AnglesRaw) that replace them in the Draw cuts (P3) are in Plot3DCols.h.

struct SampleCut {std::string dir, expr, label, onFolder, file; bool zdcKE80 = false;}; // zdcKE80: the cut file itself needs a ZDC hit with ZDC KE > 80% // see sampleCuts in Plot3D

// Plot3D has four independent parts that each write their own plots/ folders: 1 position, energy, 2_KECUT80 and 1_hitZDC; 2 the cut ratios and rank csv files of the 6_cuts folders;
// 3 their 1_additionalinfo plots; 4 their phi fits. P3_PART=n runs part n only, so the four can run as four processes side by side (Plot.pbs, plot.sh). Without it all four run.
// Part 2 is one process per cut with P3_CUT=i (0: the 1_hitZDC baseline, 1 to 11: sampleCuts in order), each writing its rank rows to cuts/rank_rows_<i>.tsv; part 5 reads those
// files and writes the rank_*.csv files. Parts 3 and 4 also take P3_CUT=i to do one cut folder only (0: 1_hitZDC and 2_KECUT80). Without P3_PART all of it runs in one process.
static bool RunPart(int n){const char* e = std::getenv("P3_PART"); return !e || !*e || std::atoi(e) == n;}
static bool RunCut(int i){const char* e = std::getenv("P3_CUT"); return !e || !*e || std::atoi(e) == i;}
// P3_SUB (0 to 3): in part 1 the four groups 1_position/1_true, 1_position/2_block, 2_energy and the 80% KE cut plots (1_hitZDC, 2_KECUT80, 0_cuts); in parts 3 and 4 of the four folders P3_CUT=0 stands for in parts 3 and 4 (1_hitZDC, 2_KECUT80, 0_cuts/1_noCuts, 0_cuts/2_theta0_1) only this one, so they can run as four processes. Without it all four.
// The six counts of AnglesRaw every part 2, 4 and 5 process needs, in one event loop: CE protons / pdg 2212 / any particle with the theta baseline, then CE protons / pdg 2212 / any particle without it.
// Kept in cuts/plot3d_totals.csv (first line the theta baseline, second the counts) and reused while the file is newer than cut_angles_all.root and the baseline is the same; else counted again (written to a
// temporary file and renamed, so processes that start together do not clash).
static std::array<double, 6> RawTotals(TTree* t, const std::string& thetaBase){
    const char* file = "cuts/plot3d_totals.csv";
    Long_t id, flags, mt = 0, mtRaw = 0; Long64_t size;
    std::array<double, 6> r{};
    if (!gSystem->GetPathInfo(file, &id, &size, &flags, &mt) && !gSystem->GetPathInfo("cuts/cut_angles_all.root", &id, &size, &flags, &mtRaw) && mt >= mtRaw){
        if (FILE* f = fopen(file, "r")){
            char line[4096]; bool ok = fgets(line, sizeof line, f) && std::string(line) == thetaBase + "\n";
            ok = ok && fscanf(f, "%lf,%lf,%lf,%lf,%lf,%lf", &r[0], &r[1], &r[2], &r[3], &r[4], &r[5]) == 6;
            fclose(f);
            if (ok) return r;
        }
    }
    const std::vector<std::string> cuts = {"isFromCE==1 && " + thetaBase, "PDG == 2212 && " + thetaBase, thetaBase, "isFromCE==1", "PDG == 2212", "1"};
    std::map<std::string, std::string> aliases;
    ROOT::RDataFrame df0(*t);
    ROOT::RDF::RNode df = WithEntryNumber(df0, df0.GetNSlots());
    std::vector<std::string> ex;
    for (const auto& c : cuts) ex.push_back(c == "1" ? "" : ColumnsExpr(t, P3(c), aliases));
    for (const auto& [alias, column] : aliases) df = df.Alias(alias, column);
    std::vector<ROOT::RDF::RResultPtr<ULong64_t>> n;
    for (const auto& e : ex) n.push_back(e.empty() ? df.Count() : df.Filter(e).Count());
    for (size_t i = 0; i < n.size(); i++) r[i] = (double)n[i].GetValue(); // the first value runs the one event loop for all six
    const std::string tmp = Form("%s.%d.tmp", file, (int)getpid());
    if (FILE* f = fopen(tmp.c_str(), "w")){
        fprintf(f, "%s\n%.0f,%.0f,%.0f,%.0f,%.0f,%.0f\n", thetaBase.c_str(), r[0], r[1], r[2], r[3], r[4], r[5]);
        fclose(f);
        std::rename(tmp.c_str(), file);
    }
    return r;
}
static bool RunSub(int i){const char* e = std::getenv("P3_SUB"); return !e || !*e || std::atoi(e) == i;}
// P2Log: one line per step with the seconds since the first call, flushed, so the job log shows where part 2 spends its time (grep "^>>> p2")
static void P2Log(const std::string& what){
    static const auto t0 = std::chrono::steady_clock::now();
    printf(">>> p2 %6.1f s  %s\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), what.c_str()); fflush(stdout);
}

void Plot3D(){
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    FastJit();
    std::string folder = "plots";

    // step folders only; CE and recon proton plots share a folder, file names end in _CE or _recon
    std::string folder_pos = folder + "/1_position/1_true/";

    std::string folder_KEcuts = folder + "/6_cuts/";


    gSystem->mkdir(folder.c_str(), true);
    gSystem->mkdir(folder_KEcuts.c_str(), true);

    double Graph = 1;

    // the maps are loaded by CutPlots.cpp (zdc) (Plot3DColumns), not here


    TFile* File = TFile::Open("cuts/cut_CE.root", "READ");

    TTree *events = (TTree*)File->Get("ParticleInfo");
    TTree *ChargeExchange = (TTree*)File->Get("ChargeExchange_values");

    // events2: ZDC protons of the CE children (cut_ZdcCE.root) and eventsAllEv: ZDC protons of every event that have a reconstructed track (cut_hitZDC_AllEv.root, CutPlots.cpp (samples): every cut on
    // it asks for isReco, so its rows are the same as cut_ZdcAllEv.root with isReco), with the Plot3D columns (Plot3DCols.h, CutPlots.cpp (zdc));
    // every Draw cut goes through P3() so isReco(fEvent, TrackID) etc. read those columns
    TTree *events2 = OpenCut("cuts/cut_ZdcCE.root", "ParticleInfo");

    // every event, not only the CE events kept in cut_CE.root: third "all protons in all events" set next to CE protons and all protons in CE events
    TTree* eventsAllEv = OpenCut("cuts/cut_hitZDC_AllEv.root", "ParticleInfo");
    // eventsAllPart: the same for every recorded particle (charged, any pdg) at the ZDC with a reconstructed track (cut_hitZDC_AllPart.root), the "all recon particles" plots
    TTree* eventsAllPart = OpenCut("cuts/cut_hitZDC_AllPart.root", "ParticleInfo");
    if (!events2 || !eventsAllEv || !eventsAllPart){Error("Plot3D", "missing cuts/cut_ZdcCE.root, cut_hitZDC_AllEv.root or cut_hitZDC_AllPart.root (run CutPlots.cpp first)"); return;}


    if (Graph == 1){

        TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
        gPad->SetRightMargin(0.30);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); // room for the rotated Y title and bottom tick labels
        gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X"); // stats/fit box to the right, fewer x labels

        // Every plot goes through this: draw expr from t with cut, set the title (title;xaxis;yaxis), let tweak adjust the
        // histogram (range, fit, ...), save to path. opt is the Draw option, e.g. "colz".
        std::set<std::string> fitFilesStarted; // stats and fit files already truncated in this run
        std::string titleSuffix; // added to the title text (before the axis titles), the KE cut of the position plots
        double nCERecoAll = -1, nCEZdcAll = -1, nRecoAll = -1, nRecoAllPart = -1; // CE protons at the ZDC before any cut (with a reco track, and all), and every reco proton (nRecoAllPart: every reco particle) there, counted on the first plot that shows them
        std::function<void()> countZdcAll = []{}; // replaced below, once the cut strings exist
        auto drawPlot = [&](TTree* t, const std::string& expr, const std::string& cut, const std::string& title,
                            const std::string& path, const char* opt = "", const std::function<void(TH1*)>& tweak = nullptr){
            gSystem->mkdir(gSystem->GetDirName(path.c_str()), true);
            canvas->Clear();
            TH1* h = HistMT(t, P3(expr), P3(cut)); // all the threads, then drawn like TTree::Draw did
            if (!h) return;
            h->Draw(opt);
            const size_t axes = title.find(';');
            h->SetTitle(PlotTitle(path, title.substr(0, axes) + titleSuffix + (axes == std::string::npos ? "" : title.substr(axes))).c_str());
            const bool ceOnly1D = t == events2 && h->GetDimension() == 1;
            const bool allPart1D = t == eventsAllPart && h->GetDimension() == 1;
            if (ceOnly1D) h->SetLineColor(kRed); // CE protons are red in every plot; the overlay tweaks below set their own colors
            if (allPart1D) h->SetLineColor(kGreen + 2); // all recon particles are green, as in the x overlays
            if (tweak) tweak(h);
            if (ceOnly1D || allPart1D) for (TObject* o : *h->GetListOfFunctions()) if (TF1* fit = dynamic_cast<TF1*>(o)) fit->SetLineColor(kBlack); // a red histogram would hide a red fit
            const std::string hName = h->GetName(); // the stats box shows the name: a readable one until the pdf is saved
            h->SetName(StatsTitle(h, title.substr(0, std::min(title.find(','), axes)) + (t == events2 ? ", CE protons" : t == eventsAllEv ? ", all recon protons" : t == eventsAllPart ? ", all recon particles" : "")).c_str()); // quantity: the title up to its first comma
            if (t == eventsAllEv || t == eventsAllPart || t == events2) countZdcAll(); // the ZDC totals the percent boxes use, counted once
            if (t == eventsAllEv || t == eventsAllPart){ // all recon protons (or particles): the CE ones (cut_ParentZdc.root) under the same cut
                canvas->cd();
                DrawCEPercent(CountMT(events2, P3(cut)), h->GetEntries(), nCERecoAll, t == eventsAllPart ? nRecoAllPart : nRecoAll); // lost: vs every reco CE proton at the ZDC
            } else if (t == events2 && nCERecoAll > 0){ // CE protons only: the share the cuts lose, vs every CE proton at the ZDC (reco or not, as the cut asks)
                const double nCE = CountMT(events2, P3(cut));
                const double nBase = cut.find("isReco") != std::string::npos ? nCERecoAll : nCEZdcAll;
                canvas->cd();
                if (nCE < nBase) DrawCEPercent(nCE, -1, nBase);
            }
            if (path.find("/6_cuts/") == std::string::npos) WritePlotNumbers(h, path, fitFilesStarted); // the 6_cuts plots keep only the pdf; stats.csv, plus fit_values.csv if h has a fit, next to the plot
            canvas->Update();
            canvas->SaveAs(path.c_str());
            h->SetName(hName.c_str());
        };
        // x:y position maps: linear color scale from 0 up to the fullest bin (empty bins stay blank)
        auto drawXY = [&](TTree* t, const std::string& expr, const std::string& cut, const std::string& title, const std::string& path){
            drawPlot(t, expr, cut, title, path, "colz", [](TH1* h){h->SetMinimum(0);});
        };
        // block position map, one bin per ZDC cell (the histogram keeps the stats box, colorbar and percent boxes): each occupied cell is a square of the cell size (the pitch of the bin, edge cells may be cut by the frame) at its exact position, colored by
        // its count on the same log scale. The squares touch, no gaps.
        auto drawXYDots = [&](TTree* t, const std::string& expr, const std::string& cut, const std::string& title, const std::string& path, double halfRange = 0){ // halfRange > 0: hard frame of +-halfRange mm, the cells at the edge are clipped by it
            drawPlot(t, expr, cut, title, path, "colz", [halfRange](TH1* h){DrawBlockSquares(h, halfRange, false);}); // BlockSquares.h, linear color scale
        };

        // 1D plot on a hard frame of +-half mm: the bars of the cells past it are clipped by the frame (like 4_xy_zoom_recon.pdf); call in drawPlot's tweak, after the fit
        auto clipX = [&](TH1* h, double half){
            canvas->Clear();
            TH1* fr = gPad->DrawFrame(-half, 0, half, 1.1 * h->GetMaximum());
            fr->SetTitle(h->GetTitle()); fr->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle()); fr->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
            h->Draw("sames");
        };

        const std::string zdcProtons = "boundaryType == 2 && pdg == 2212";
        // position plots: only protons with a reconstructed track
        const std::string zdcReco = zdcProtons + " && isReco(fEvent, TrackID)", zdcRecoAll = zdcReco;
        // the same for every recorded particle (any pdg, eventsAllPart): ZDC hit and a recon track, and for the block plots the block residual cut too
        const std::string zdcRecoPartAll = "boundaryType == 2 && isReco(fEvent, TrackID)", zdcRecoDPartAll = "boundaryType == 2 && isRecoD(fEvent, TrackID)";
        // baseline of every plot and count in 6_cuts: recon theta and recon block theta from 0 to 1 degree (NaN fails), on top of the ZDC hit and the cut
        const std::string thetaBase = "rawVal(fEvent, TrackID, 9) >= 0 && rawVal(fEvent, TrackID, 9) <= 1 && rawVal(fEvent, TrackID, 11) >= 0 && rawVal(fEvent, TrackID, 11) <= 1", zdcBase = zdcRecoAll + " && " + thetaBase, zdcBasePart = zdcRecoPartAll + " && " + thetaBase; // zdcBasePart: every recorded particle (eventsAllPart), not only the protons
        countZdcAll = [&]{
            if (nCERecoAll >= 0) return;
            nCERecoAll = CountMT(events2, P3(zdcRecoAll)); nCEZdcAll = CountMT(events2, P3(zdcProtons)); nRecoAll = CountMT(eventsAllEv, P3(zdcRecoAll)); nRecoAllPart = CountMT(eventsAllPart, P3(zdcRecoPartAll));
        };
        const std::string kNoCutNote = "no theta cut"; // 0_cuts: only the ZDC hit and a recon track
        const std::string kThetaNote = "recon theta and recon block theta 0-1 deg"; // the thetaBase cut, named in the titles of the 6_cuts plots
        const std::string zdcRecoDAll = zdcProtons + " && isRecoD(fEvent, TrackID)"; // block position plots: also the block residual cut

        // Position plots have no KE cut: 1_position/1_true/ and 1_position/2_block/
        std::vector<std::pair<std::string, std::string>> posKes = {{"", ""}}; // unused folder, extra cut
        std::vector<std::string> posKeLabels = {""};
        if (RunPart(1)) for (size_t pk = 0; pk < posKes.size(); pk++){
            const std::string folder_pos = folder + "/1_position/1_true/", folder_posD = folder + "/1_position/2_block/";
            const std::string zdcReco = zdcRecoAll + posKes[pk].second, zdcRecoD = zdcRecoDAll + posKes[pk].second;
            const std::string xySlice = zdcReco + " && fabs(fY) < 2", xySliceD = zdcRecoD + " && fabs(fDy) < 2";
            titleSuffix = ", has recon track, no KE cut" + posKeLabels[pk]; // the cut of the 1_true plots: hit ZDC and a recon track
            const std::string zdcRecoP = zdcRecoPartAll + posKes[pk].second, xySliceP = zdcRecoP + " && fabs(fY) < 2";
            if (RunSub(0)){ // P3_SUB 0: 1_true
            // ZDC position, true (fX, fY)
            drawXY(eventsAllEv, "fY:fX>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at ZDC, all recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_pos + "2_xy_recon.pdf");
            drawXY(eventsAllEv, "fY:fX>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at ZDC, all recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_pos + "4_xy_zoom_recon.pdf");
            drawXY(events2, "fY:fX>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at ZDC, CE recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_pos + "1_xy_CE.pdf");
            drawXY(events2, "fY:fX>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at ZDC, CE recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_pos + "3_xy_zoom_CE.pdf");

            gStyle->SetOptFit(11); // fit parameters with errors (the fits are unbinned: no chi2 or probability)
            drawPlot(events2, "fX>>hSliceX200(100,-200,200)", xySlice, "X Position at ZDC, CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_pos + "5_x_fit_200_CE.pdf", "", [](TH1* h){FitCore(h);});
            drawPlot(events2, "fX>>hSliceX(100,-50,50)", xySlice, "X Position at ZDC, CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_pos + "7_x_fit_50_CE.pdf", "", [](TH1* h){FitCore(h);});

            // the same plots for every recorded particle that hit the ZDC and has a recon track (any pdg), own pdf next to each all recon protons one: <name>_allParticles.pdf
            drawXY(eventsAllPart, "fY:fX>>h2(100,-200,200,100,-200,200)", zdcRecoP, "Position at ZDC, all recon particles that hit ZDC;X Position (mm);Y Position (mm)", folder_pos + "2_xy_recon_allParticles.pdf");
            drawXY(eventsAllPart, "fY:fX>>h2(100,-50,50,100,-50,50)", zdcRecoP, "Position at ZDC, all recon particles that hit ZDC;X Position (mm);Y Position (mm)", folder_pos + "4_xy_zoom_recon_allParticles.pdf");

            // the x slices as one overlay each, no fit: CE recon protons (red), all recon protons (blue), all recon particles (green)
            struct XOv {std::string bins, file, tag;};
            for (const auto& o : std::vector<XOv>{{"(100,-200,200)", "9_x_overlay_200.pdf", "200"}, {"(100,-50,50)", "10_x_overlay_50.pdf", "50"}}){
                drawPlot(eventsAllEv, "fX>>hXOv" + o.tag + o.bins, xySlice, "X Position at ZDC, all recon particles, all recon protons and CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_pos + o.file, "", [&](TH1* h){
                    TH1* hCE = HistMT(events2, P3("fX>>hXOvCE" + o.tag + o.bins), P3(xySlice)); hCE->Draw("same");
                    TH1* hPart = HistMT(eventsAllPart, P3("fX>>hXOvPart" + o.tag + o.bins), P3(xySliceP)); hPart->Draw("same");
                    h->SetLineColor(kBlue); hCE->SetLineColor(kRed); hPart->SetLineColor(kGreen + 2); hCE->SetStats(0); hPart->SetStats(0);
                    h->SetMaximum(1.2 * std::max({h->GetMaximum(), hCE->GetMaximum(), hPart->GetMaximum()}));
                    TLegend* leg = new TLegend(0.745, 0.30, 0.99, 0.42); leg->SetTextSize(0.012); // below the stats box
                    leg->AddEntry(hPart, "all recon particles that hit ZDC", "l");
                    leg->AddEntry(h, "all recon protons that hit ZDC", "l");
                    leg->AddEntry(hCE, "CE recon protons that hit ZDC", "l");
                    leg->Draw();
                });
            }


            }
            if (RunSub(1)){ // P3_SUB 1: 2_block
            titleSuffix = ", has recon track, block residual < 2 mm, no KE cut" + posKeLabels[pk]; // the cut of the 2_block plots: isRecoD
            // ZDC position, block (fDx, fDy): the ZDC is nCols_detector x nRows_detector blocks over 400 mm. Every bin holds the same whole number of blocks (BlockBins, FitCore.h):
            // 100 bins over 100 mm on a pitch of 400 / 801 mm put an extra block in the bin that starts at x = 0 (a cross at x = 0 and y = 0), and in one bin of ten along the x slices (spikes).
            const int nColsD = (int)ReadSettings().at("nCols_detector"), nRowsD = (int)ReadSettings().at("nRows_detector");
            const std::string binsFull = BlockBins(nColsD, 4.0, 200), binsFullY = BlockBins(nRowsD, 4.0, 200), binsZoom = BlockBins(nColsD, 1.0, 50), binsZoomY = BlockBins(nRowsD, 1.0, 50), binsSlice200 = BlockBins(nColsD, 400.0 / nColsD, 200), binsSlice50 = BlockBins(nColsD, 400.0 / nColsD, 53.5); // x slices: one bin per ZDC cell, the +-50 mm ones with whole cells past 50 (clipped by the frame, clipX)
            // the 4 xy plots: one marker per occupied ZDC cell (nCols_detector x nRows_detector over 400 mm), at the exact cell position (drawXYDots)
            const std::string binsCell = BlockBins(nColsD, 400.0 / nColsD, 200), binsCellY = BlockBins(nRowsD, 400.0 / nRowsD, 200), binsCellZoom = BlockBins(nColsD, 400.0 / nColsD, 53.5), binsCellZoomY = BlockBins(nRowsD, 400.0 / nRowsD, 53.5); // whole cells past +-50: the frame is +-50, the edge cells are clipped
            const std::string cellNote = Form(", 1 square per ZDC cell (%d x %d over 400 mm)", nColsD, nRowsD), zoomNote = cellNote + ", edge cells clipped at 50 mm";
            drawXYDots(eventsAllEv, "fDy:fDx>>h2(" + binsCell + "," + binsCellY + ")", zdcRecoD, "Position at ZDC BLOCK, all recon protons that hit ZDC" + cellNote + ";X Position (mm);Y Position (mm)", folder_posD + "2_xy_recon.pdf");
            drawXYDots(eventsAllEv, "fDy:fDx>>h2(" + binsCellZoom + "," + binsCellZoomY + ")", zdcRecoD, "Position at ZDC BLOCK, all recon protons that hit ZDC" + zoomNote + ";X Position (mm);Y Position (mm)", folder_posD + "4_xy_zoom_recon.pdf", 50);
            drawXYDots(events2, "fDy:fDx>>h2(" + binsCell + "," + binsCellY + ")", zdcRecoD, "Position at ZDC BLOCK, CE recon protons that hit ZDC" + cellNote + ";X Position (mm);Y Position (mm)", folder_posD + "1_xy_CE.pdf");
            drawXYDots(events2, "fDy:fDx>>h2(" + binsCellZoom + "," + binsCellZoomY + ")", zdcRecoD, "Position at ZDC BLOCK, CE recon protons that hit ZDC" + zoomNote + ";X Position (mm);Y Position (mm)", folder_posD + "3_xy_zoom_CE.pdf", 50);

            // the same block xy plots for every recorded particle that hit the ZDC and has a recon track passing the block residual cut (any pdg)
            const std::string zdcRecoDP = zdcRecoDPartAll + posKes[pk].second;
            drawXYDots(eventsAllPart, "fDy:fDx>>h2(" + binsCell + "," + binsCellY + ")", zdcRecoDP, "Position at ZDC BLOCK, all recon particles that hit ZDC" + cellNote + ";X Position (mm);Y Position (mm)", folder_posD + "2_xy_recon_allParticles.pdf");
            drawXYDots(eventsAllPart, "fDy:fDx>>h2(" + binsCellZoom + "," + binsCellZoomY + ")", zdcRecoDP, "Position at ZDC BLOCK, all recon particles that hit ZDC" + zoomNote + ";X Position (mm);Y Position (mm)", folder_posD + "4_xy_zoom_recon_allParticles.pdf", 50);

            // the six x fits have no block residual cut: ZDC hit and a recon track only
            titleSuffix = ", has recon track, no KE cut" + posKeLabels[pk];
            const std::string xySliceB = zdcReco + " && fabs(fDy) < 2", xySliceBP = zdcRecoP + " && fabs(fDy) < 2";
            drawPlot(events2, "fDx>>hSliceDX200(" + binsSlice200 + ")", xySliceB, "X Position at ZDC BLOCK, CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_posD + "5_x_fit_200_CE.pdf", "", [](TH1* h){FitBlocks(h);}); //whole blocks per bin
            drawPlot(events2, "fDx>>hSliceDX(" + binsSlice50 + ")", xySliceB, "X Position at ZDC BLOCK, CE recon protons that hit ZDC, |Y| < 2 mm, edge cells clipped at 50 mm;X Position (mm);Count", folder_posD + "7_x_fit_50_CE.pdf", "", [&](TH1* h){FitBlocks(h); clipX(h, 50);});
            drawPlot(eventsAllEv, "fDx>>hSliceDXAllEv200(" + binsSlice200 + ")", xySliceB, "X Position at ZDC BLOCK, all recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_posD + "6_x_fit_200_recon.pdf", "", [](TH1* h){FitBlocks(h);});
            drawPlot(eventsAllEv, "fDx>>hSliceDXAllEv(" + binsSlice50 + ")", xySliceB, "X Position at ZDC BLOCK, all recon protons that hit ZDC, |Y| < 2 mm, edge cells clipped at 50 mm;X Position (mm);Count", folder_posD + "8_x_fit_50_recon.pdf", "", [&](TH1* h){FitBlocks(h); clipX(h, 50);});
            drawPlot(eventsAllPart, "fDx>>hSliceDXAllPart200(" + binsSlice200 + ")", xySliceBP, "X Position at ZDC BLOCK, all recon particles that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_posD + "6_x_fit_200_recon_allParticles.pdf", "", [](TH1* h){FitBlocks(h);});
            drawPlot(eventsAllPart, "fDx>>hSliceDXAllPart(" + binsSlice50 + ")", xySliceBP, "X Position at ZDC BLOCK, all recon particles that hit ZDC, |Y| < 2 mm, edge cells clipped at 50 mm;X Position (mm);Count", folder_posD + "8_x_fit_50_recon_allParticles.pdf", "", [&](TH1* h){FitBlocks(h); clipX(h, 50);});

            // recon block phi (rawVal 15, 0 to 360 deg) of the same |Y| < 2 mm events as the x fits above, no fit
            const std::string phiOk = " && rawOk(fEvent, TrackID, 15)";
            // one version per x window, like the x fits: |X| < 200 mm (files 9, 10) and |X| < 50 mm (files 11, 12)
            struct PhiWin {int w; std::string nCE, nRecon;};
            for (const auto& pw : std::vector<PhiWin>{{200, "9", "10"}, {50, "11", "12"}}){
                const std::string tag = std::to_string(pw.w), xWin = Form(" && fabs(fDx) < %d", pw.w), winNote = Form(", |X| < %d mm", pw.w);
                drawPlot(events2, "rawVal(fEvent, TrackID, 15)>>hPhiDY2CE" + tag + "(100,0,360)", xySliceB + phiOk + xWin, "Recon block phi, CE recon protons that hit ZDC, |Y| < 2 mm" + winNote + ";Recon block phi (degrees);Count", folder_posD + pw.nCE + "_phiBlock_" + tag + "_CE.pdf");
                drawPlot(eventsAllEv, "rawVal(fEvent, TrackID, 15)>>hPhiDY2AllEv" + tag + "(100,0,360)", xySliceB + phiOk + xWin, "Recon block phi, all recon protons that hit ZDC, |Y| < 2 mm" + winNote + ";Recon block phi (degrees);Count", folder_posD + pw.nRecon + "_phiBlock_" + tag + "_recon.pdf");
                drawPlot(eventsAllPart, "rawVal(fEvent, TrackID, 15)>>hPhiDY2AllPart" + tag + "(100,0,360)", xySliceBP + phiOk + xWin, "Recon block phi, all recon particles that hit ZDC, |Y| < 2 mm" + winNote + ";Recon block phi (degrees);Count", folder_posD + pw.nRecon + "_phiBlock_" + tag + "_recon_allParticles.pdf");
            }
            }
            titleSuffix.clear();
        }

        // KE x axes: 0 to 30k MeV for the 3_KEcuts plots, 10k to 30k MeV for the 2_energy plots
        struct KERange {double lo, hi; std::string tag;};
        const std::vector<KERange> keRanges = {{0, 30000, "_0to30k"}}; // 3_KEcuts: full range only
        const std::vector<KERange> energyRanges = {{0, 30000, "_0to30k"}, {10000, 30000, "_10to30k"}}; // 2_energy: 0 to 30k and 10k to 30k
        auto binSpec = [](const KERange& r){return std::string(Form("(100,%g,%g)", r.lo, r.hi));};

        // 2_energy/ holds the energy plots, with the KE cut folders ke_cut_<pct>/ next to them
        if (RunPart(1) && RunSub(2)) {
            const std::string folder_smear = folder + "/2_energy/";
            gSystem->mkdir(folder_smear.c_str(), true);
            const std::string zdcProtonsE = zdcProtons, zdcRecoE = zdcReco;
            titleSuffix = ", no KE cut";
            // Energy plots: 2_energy/<n>_<quantity><range>_<CE or recon>.pdf, numbered 1, 2 (true KE CE, recon), 3, 4 (smeared KE), 5 to 8 (smearing)
            // True KE (fKEnergy) and smeared KE (fDEnergy) at the ZDC.
            struct KESet {TTree* tree; std::string tag, label, cut;}; // tag ends the file name, CE and recon files share a folder
            const std::vector<KESet> keSets = {{events2, "_CE", "CE protons that hit ZDC (reconstructed or not)", zdcProtonsE}, {eventsAllEv, "_recon", "all recon protons that hit ZDC", zdcRecoE}};
            struct KEQty {std::string branch, file, axis, title;};
            std::vector<KEQty> keQtys = {{"fKEnergy", "true_ke", "Kinetic Energy (MeV)", "True Kinetic Energy at ZDC"}};
            if (events->GetBranch("fDEnergy")) keQtys.push_back({"fDEnergy", "smeared_ke", "Smeared Energy fDEnergy (MeV)", "Smeared Energy at ZDC"}); // files from before fDEnergy was added don't have it
            // one log y plot per quantity: CE protons (red), all recon protons (blue) and all recon particles (green): 1_true_ke, 2_smeared_ke
            for (const auto& q : keQtys){
                for (const auto& r : energyRanges){
                    drawPlot(eventsAllEv, q.branch + ">>hKE" + r.tag.substr(1) + binSpec(r), zdcRecoE, q.title + ", all recon particles, all recon protons and CE protons that hit ZDC (reconstructed or not);" + q.axis + ";Count",
                             folder_smear + std::to_string(&q - &keQtys[0] + 1) + "_" + q.file + r.tag + ".pdf", "", [&](TH1* h){
                        TH1* hCE = HistMT(events2, P3(q.branch + ">>hKECE" + r.tag.substr(1) + binSpec(r)), P3(zdcProtonsE)); hCE->Draw("same");
                        TH1* hPart = HistMT(eventsAllPart, P3(q.branch + ">>hKEPart" + r.tag.substr(1) + binSpec(r)), P3(zdcRecoPartAll)); hPart->Draw("same"); // every charged particle: contains the protons
                        h->SetLineColor(kBlue); hCE->SetLineColor(kRed); hPart->SetLineColor(kGreen + 2); hCE->SetStats(0); hPart->SetStats(0);
                        h->SetMinimum(0.5); h->SetMaximum(2.0 * std::max({h->GetMaximum(), hCE->GetMaximum(), hPart->GetMaximum()})); canvas->SetLogy(1);
                        TLegend* leg = new TLegend(0.745, 0.30, 0.99, 0.42); leg->SetTextSize(0.012); // below the stats box
                        leg->AddEntry(hPart, "all recon particles that hit ZDC", "l");
                        leg->AddEntry(h, "all recon protons that hit ZDC", "l");
                        leg->AddEntry(hCE, "CE protons that hit ZDC (reco or not)", "l");
                        leg->Draw();
                    });
                    canvas->SetLogy(0);
                }
            }

            // ZDC energy smearing (fDEnergy = Gaussian smear of fKEnergy, /zdc/resStochastic and /zdc/resConstant in the sim)
            if (events->GetBranch("fDEnergy")){ // files from before fDEnergy was added don't have it
                gSystem->mkdir(folder_smear.c_str(), true);

                // fKEnergy = true KE, fDEnergy = smeared KE, dE = fDEnergy - fKEnergy

                drawPlot(events2, "fDEnergy-fKEnergy>>hSmearDiffCE(100,-10000,10000)", zdcProtonsE, "Energy Difference at ZDC, CE protons that hit ZDC (reconstructed or not);fDEnergy - fKEnergy (MeV);Count", folder_smear + "5_diff_CE.pdf", "", [](TH1* h){FitCore(h);});

                // percent error = 100 * (fDEnergy - fKEnergy) / fKEnergy
                const std::string pct = "100*((fDEnergy-fKEnergy)/fKEnergy)"; // bracketed: the C++ of HistMT would divide 100*(a-b) by c, TTree::Draw 100*(a-b/c)
                drawPlot(events2, pct + ">>hSmearPctCE(100,-100,100)", zdcProtonsE, "Percent Error at ZDC, CE protons that hit ZDC (reconstructed or not);Percent Error (%);Count", folder_smear + "7_pct_CE.pdf", "", [](TH1* h){FitCore(h);});
                

                // the smearing quantities as one log y overlay each: CE protons (red), all recon protons (blue), all recon particles (green)
                struct SmearOv {std::string expr, bins, file, axis, title, tag;};
                const std::vector<SmearOv> smearOvs = {{"fDEnergy-fKEnergy", "(100,-10000,10000)", "9_diff_overlay.pdf", "fDEnergy - fKEnergy (MeV)", "Energy Difference at ZDC", "Diff"},
                                                       {pct, "(100,-100,100)", "10_pct_overlay.pdf", "Percent Error (%)", "Percent Error at ZDC", "Pct"}};
                for (const auto& o : smearOvs){
                    drawPlot(eventsAllEv, o.expr + ">>hOv" + o.tag + o.bins, zdcRecoE, o.title + ", all recon particles, all recon protons and CE protons that hit ZDC (reconstructed or not);" + o.axis + ";Count",
                             folder_smear + o.file, "", [&](TH1* h){
                        TH1* hCE = HistMT(events2, P3(o.expr + ">>hOvCE" + o.tag + o.bins), P3(zdcProtonsE)); hCE->Draw("same");
                        TH1* hPart = HistMT(eventsAllPart, P3(o.expr + ">>hOvPart" + o.tag + o.bins), P3(zdcRecoPartAll)); hPart->Draw("same");
                        h->SetLineColor(kBlue); hCE->SetLineColor(kRed); hPart->SetLineColor(kGreen + 2); hCE->SetStats(0); hPart->SetStats(0);
                        h->SetMinimum(0.5); h->SetMaximum(2.0 * std::max({h->GetMaximum(), hCE->GetMaximum(), hPart->GetMaximum()})); canvas->SetLogy(1);
                        TLegend* leg = new TLegend(0.745, 0.30, 0.99, 0.42); leg->SetTextSize(0.012); // below the stats box
                        leg->AddEntry(hPart, "all recon particles that hit ZDC", "l");
                        leg->AddEntry(h, "all recon protons that hit ZDC", "l");
                        leg->AddEntry(hCE, "CE protons that hit ZDC (reco or not)", "l");
                        leg->Draw();
                    });
                    canvas->SetLogy(0);
                }
            }
            else{Warning("Plot3D", "no fDEnergy branch in cut_CE.root, skipping smearing plots");}
            titleSuffix.clear();
        }

        const double initialKE = kBeamKE; // MeV, beam KE (FitCore.h)
        // KE cuts on the ZDC-smeared energy (fDEnergy) only, in ke_cut_<pct>/. Files from before fDEnergy was added get no KE cuts.
        struct CutEnergy {std::string branch, dir, label;};
        std::vector<CutEnergy> cutEnergies;
        if (events->GetBranch("fDEnergy")) cutEnergies.push_back({"fDEnergy", "KEDetector", "ZDC KE"});
        if (RunPart(1) && RunSub(3)) for (int pct : {80}){
            const std::string zdcReco = zdcRecoAll;
            std::string cut = std::to_string(initialKE * pct / 100.0);
            std::string pctStr = std::to_string(pct) + "% of initial KE";
            for (const auto& e : cutEnergies){
                std::string folder_cut = folder_KEcuts + "2_KECUT" + std::to_string(pct) + "/";
                gSystem->mkdir(folder_cut.c_str(), true);
                const std::string axis = ";" + e.label + " (MeV);Count";
                const double cutVal = initialKE * pct / 100.0;
                const std::string realOutCut = zdcReco + " && " + e.branch + " < " + cut;
                const std::string fakeKeptCut = zdcReco + " && " + e.branch + " > " + cut;
                // counts over all KE (not only the plotted range): CE protons left out of the cut out of all CE protons, all recon protons kept in out of all recon protons
                const long long nRealOut = CountMT(events2, P3(realOutCut)), nRealTotal = CountMT(events2, P3(zdcReco));
                const long long nFakeKept = CountMT(eventsAllEv, P3(fakeKeptCut)), nFakeTotal = CountMT(eventsAllEv, P3(zdcReco));
                const double realOutPct = 100.0 * nRealOut / std::max<long long>(1, nRealTotal);
                const double fakeKeptPct = 100.0 * nFakeKept / std::max<long long>(1, nFakeTotal);
                for (const auto& r : keRanges){
                // real (reconstructed CE protons) and all reconstructed protons of all events on one plot, with the cut and what it does to each
                drawPlot(events2, e.branch + ">>hCutRealOv" + binSpec(r), zdcReco, e.label + " of Protons, CE recon protons vs all recon protons, both hit ZDC, " + pctStr + " cut" + axis, folder_cut + "3_real_vs_recon_" + e.dir + r.tag + ".pdf", "", [&](TH1* h){
                    h->SetLineColor(kRed);
                    TH1* hf = HistMT(eventsAllEv, P3(e.branch + ">>hCutFakeOv" + binSpec(r)), P3(zdcReco)); hf->Draw("same");
                    hf->SetLineColor(kBlue);
                    TH1* hp = HistMT(eventsAllPart, P3(e.branch + ">>hCutPartOv" + binSpec(r)), P3(zdcRecoPartAll)); hp->Draw("same"); // every recon particle (any pdg) that hit the ZDC
                    hp->SetLineColor(kGreen + 2);
                    const long long nPartKept = CountMT(eventsAllPart, P3(zdcRecoPartAll + " && " + e.branch + " > " + cut)), nPartTotal = CountMT(eventsAllPart, P3(zdcRecoPartAll));
                    h->SetMinimum(0.5); h->SetMaximum(2.0 * std::max(h->GetMaximum(), std::max(hf->GetMaximum(), hp->GetMaximum()))); canvas->SetLogy(1); // log y, full x range of the file name, the cut is the dashed line
                    TLine* cl = new TLine(cutVal, 0.5, cutVal, h->GetMaximum()); cl->SetLineStyle(2); cl->Draw();
                    TLegend* leg = new TLegend(0.745, 0.44, 0.99, 0.55); leg->SetTextSize(0.012); // below the stats box; both samples hit the ZDC (title)
                    leg->AddEntry(h, Form("CE recon protons, %lld/%lld = %.1f%% left out", nRealOut, nRealTotal, realOutPct), "l");
                    leg->AddEntry(hf, Form("all recon protons, %lld/%lld = %.1f%% kept in", nFakeKept, nFakeTotal, fakeKeptPct), "l");
                    leg->AddEntry(hp, Form("all recon particles, %lld/%lld = %.1f%% kept in", nPartKept, nPartTotal, 100.0 * nPartKept / std::max<long long>(1, nPartTotal)), "l");
                    leg->Draw();
                });
                canvas->SetLogy(0);
                }
            }
        }

        // Cuts on the reconstructed protons, no KE cut, each in its own folder under 6_cuts: dir, extra Draw cut, label, and the info variable folder the cut is on
        // (those plots sit in the cut folder itself, the other variables in its 1_additionalinfo). The same KE plot as 1_hitZDC with the cut applied goes in 1_additionalinfo.
        const std::vector<SampleCut> sampleCuts = [&]{
          std::vector<SampleCut> v = {
            {"3_ResCuts", "rawVal(fEvent, TrackID, 1) < " + std::to_string(kRes4Mm), "cut resDRms (block residual) < 4 mm", "3_ResCuts", "cut_ResCuts"},
            {"4_nParticles/1_particleCut", "rawVal(fEvent, TrackID, 6) <= 1", "cut at most 1 reconstructed particle in the event", "4_nParticles", "cut_nParticles1"},
            {"4_nParticles/2_particleCut", "rawVal(fEvent, TrackID, 6) <= 2", "cut at most 2 reconstructed particles in the event", "4_nParticles", "cut_nParticles2"},
            {"5_nHits", "(rawVal(fEvent, TrackID, 10) == 2 || rawVal(fEvent, TrackID, 10) == 4 || rawVal(fEvent, TrackID, 10) == 6 || rawVal(fEvent, TrackID, 10) == 8)", "cut 2, 4, 6 or 8 tracker hits", "5_nHits", "cut_nHits"},
          };
          // 6_megaCut (no ZDC requirement), 7_megaCut2 (ZDC KE > 80%) and 9_megaCut3 (megaCut2 and 2, 4, 6 or 8 charged tracker hits in the event, all particles added up), each with
          // 1_1particlecut: at most 1 particle, 2_2particlecut: at most 2
          static const char* const megaNames[2] = {"1_1particlecut", "2_2particlecut"}; // folder names of 6_megaCut, 7_megaCut2 and 9_megaCut3
          auto addMega = [&](int level){ // 0: megaCut, 1: megaCut2, 2: megaCut3
              for (int i = 1; i <= 2; i++){
                  const int maxP = i;
                  const bool zdc = level >= 1;
                  const std::string top = level == 0 ? "6_megaCut" : level == 1 ? "7_megaCut2" : "9_megaCut3";
                  std::string expr = "rawVal(fEvent, TrackID, 6) <= " + std::to_string(maxP) + " && (rawVal(fEvent, TrackID, 10) == 2 || rawVal(fEvent, TrackID, 10) == 4 || rawVal(fEvent, TrackID, 10) == 6 || rawVal(fEvent, TrackID, 10) == 8) && rawVal(fEvent, TrackID, 9) >= 0 && rawVal(fEvent, TrackID, 9) <= 1 && rawVal(fEvent, TrackID, 11) >= 0 && rawVal(fEvent, TrackID, 11) <= 1";
                  if (zdc) expr += " && fDEnergy > " + std::to_string(initialKE * 0.80);
                  if (level == 2) expr += " && (rawVal(fEvent, TrackID, 17) == 2 || rawVal(fEvent, TrackID, 17) == 4 || rawVal(fEvent, TrackID, 17) == 6 || rawVal(fEvent, TrackID, 17) == 8) && rawVal(fEvent, TrackID, 7) == 1";
                  const std::string label = std::string(level == 0 ? "MEGA CUT: at most " : level == 1 ? "MEGA CUT 2: ZDC KE > 80%, at most " : "MEGA CUT 3: ZDC KE > 80%, at most ") + std::to_string(maxP) + (maxP == 1 ? " reconstructed particle" : " reconstructed particles") + " in event AND number of times the reconstructed particle hits tracker [2,4,6,8], recon theta and recon block theta from 0 degrees to 1 degree" + (level == 2 ? " AND charged tracker plane hits in the event (all particles added up) [2,4,6,8] AND exactly 1 charged particle at the ZDC in the event" : "");
                  const char* file = level == 0 ? "cut_megacut_" : level == 1 ? "cut_megacut2_" : "cut_megacut3_";
                  v.push_back({top + "/" + megaNames[i - 1], expr, label, top, std::string(file) + std::to_string(i), zdc});
              }
          };
          addMega(0); addMega(1);
          v.push_back({"8_nZdcParticles/1_1zdcparticlecut", "rawVal(fEvent, TrackID, 7) <= 1", "cut at most 1 charged particle at the ZDC in the event", "8_nZdcParticles", "cut_nZdc1", false}); // P3_CUT 9
          addMega(2); // P3_CUT 10 and 11
          return v;
        }();
        struct RankRow {std::string cut, sample; double nCEin, nAllIn, nCEtotal, nCEreconBefore, fom[2], pol[2], polErr[2], nPartIn, nNonProtBefore, nNonPartBefore;}; // nAllIn: recon protons (pdg 2212) in the cut, nPartIn: every recon particle in the cut, nNon*Before: the non-CE protons / particles of the sample the cut is made on (the base of the row), so out of cut = before - in;  fom, pol (fitted P), polErr: [0] recon theta / phi, [1] recon block theta / phi
        double nCEreconTotal = 0; // every CE proton with a fitted track (AnglesRaw, isFromCE), no ZDC requirement: CE_before of the QualityStudy csv files, the CE_recon_before_cut of the csv files
        // AnglesRaw: every proton with a fitted track, ZDC hit or not. The cut folders count, fit and plot these (no ZDC requirement), except with the 80% ZDC KE cut.
        TTree* rawT = OpenCut("cuts/cut_angles_all.root", "AnglesRaw"); // every fitted proton, the Plot3D columns are aliases of its branches (Plot3DRawAliases)
        if (rawT) Plot3DRawAliases(rawT);
        gROOT->cd();
        // the rows of AnglesRaw that pass a cut of sampleCuts, from its cut file (CutPlots.cpp (samples)), with the same Plot3D aliases as rawT
        auto cutTree = [&](const SampleCut& sc){
            TTree* t = OpenCut("cuts/" + sc.file + ".root", "AnglesRaw");
            if (t && !t->GetAlias("raw_0")) Plot3DRawAliases(t);
            gROOT->cd();
            return t;
        };
        P2Log("start");
        double nNonProtTh = 0, nNonPartTh = 0, nNonProtAll = 0, nNonPartAll = 0; // non-CE fitted protons (pdg 2212) and particles before any cut: with the theta baseline (every folder but 0_cuts/1_noCuts) and without it
        if (rawT && (RunPart(0) || RunPart(2) || RunPart(4) || RunPart(5))){ // P3_PART=0 only makes cuts/plot3d_totals.csv (one pass of all the threads, Plot.pbs / plot.sh run it before the other parts)
            const std::array<double, 6> tot = RawTotals(rawT, thetaBase);
            nCEreconTotal = tot[0];
            nNonProtTh = tot[1] - tot[0]; nNonPartTh = tot[2] - tot[0];
            nNonProtAll = tot[4] - tot[3]; nNonPartAll = tot[5] - tot[3];
        }
        P2Log("totals done");
        const std::string rawCE = "isFromCE==1";
        std::vector<RankRow> rankRows; // every cut folder, with and without the 80% KE cut, for the rank_*.csv files below
        // FOM like 5_fomCalc/<theta set>/5_integral_ce_x_A2_dtheta.pdf (PlotCEFraction.cpp): sum over 0.05 deg theta bins of N_CE(bin) / N_neutrons * A_y^2(theta bin center) * 0.05 deg,
        // N_CE = CE recon protons that pass q (zdcRecoAll and the cut), t from theta and A_y from t with the model formula, theta from 0 to 15 deg
        const double nNeutrons = ReadSettings().at("TotalNeutrons"), apScaleFom = ReadSettings().at("APScale");
        auto fomStats = [&](TTree* t, const std::string& q, double* fom){
            // one pass over t for both thetas: raw value and rawOk flag of each, the rows without a finite value are skipped like the rawOk cut did
            const auto c = Columns(t, {P3("rawVal(fEvent, TrackID, 9)"), P3("rawOk(fEvent, TrackID, 9)"), P3("rawVal(fEvent, TrackID, 11)"), P3("rawOk(fEvent, TrackID, 11)")}, P3(q));
            for (int j = 0; j < 2; j++){
                double sum = 0;
                for (size_t i = 0; i < c[2 * j].size(); i++){
                    if (c[2 * j + 1][i] == 0) continue;
                    const double th = c[2 * j][i];
                    if (th < 0 || th >= 15) continue;
                    const double A = APFromT(kBeamMomentum, TFromAngle(kBeamMomentum, (std::floor(th / 0.05) + 0.5) * 0.05 * TMath::Pi() / 180), apScaleFom);
                    sum += A * A * 0.05 / nNeutrons;
                }
                fom[j] = sum;
            }
        };
        // phi cos fit of the CE recon protons that pass q (zdcRecoAll and the cut), like 3_ReconCalcVals/2_phiRe_cosfit.pdf, at their recon phi and recon theta (j 0) or recon block phi and block theta (j 1):
        // N(phi) = N0 (1 + P <A> cos(phi)), N0 = entries / 50 bins fixed, <A> = mean model A_y at the recon (block) theta of those protons fixed, only P fitted. Returns P, its error in
        // pol, polErr; the plot goes to path unless it is empty
        auto phiFit = [&](TTree* t, const std::string& qIn, int j, const std::string& title, const std::string& path, double* pol, double* polErr){
            const int kPhi = j ? 15 : 14, kTheta = j ? 11 : 9;
            std::string q = qIn + Form(" && rawOk(fEvent, TrackID, %d) && rawOk(fEvent, TrackID, %d)", kPhi, kTheta);
            if (j) q += Form(" && rawVal(fEvent, TrackID, %d) >= %g && rawVal(fEvent, TrackID, %d) <= %g", kPhi, kBlockPhiLo, kPhi, kBlockPhiHi); // block phi: no spike at 0
            *pol = 0; *polErr = 0;
            const auto cols = Columns(t, {P3(Form("rawVal(fEvent, TrackID, %d)", kTheta)), P3(Form("rawVal(fEvent, TrackID, %d)", kPhi))}, P3(q)); // theta, phi of the protons that pass q, in entry order
            const Long64_t n = cols[0].size();
            double aMean = 0;
            const double apScalePhi = ReadSettings().at("APScale"); // settings.csv is read from disk on every call, keep it out of the event loop
            for (Long64_t i = 0; i < n; i++)
                aMean += APFromT(kBeamMomentum, TFromAngle(kBeamMomentum, (std::floor(cols[0][i] / 0.05) + 0.5) * 0.05 * TMath::Pi() / 180), apScalePhi) / std::max<Long64_t>(1, n);
            canvas->Clear(); canvas->SetRightMargin(0.25); canvas->SetTopMargin(0.17);
            const std::string hn = j ? "hPhiReDfit" : "hPhiRefit";
            delete gDirectory->Get(hn.c_str());
            TH1* h = new TH1F(hn.c_str(), hn.c_str(), 50, 0, 360); // filled in entry order like the serial Draw (same mean and RMS)
            for (double ph : cols[1]) h->Fill(ph);
            t->TAttLine::Copy(*(TAttLine*)h); t->TAttFill::Copy(*(TAttFill*)h); t->TAttMarker::Copy(*(TAttMarker*)h);
            if (!path.empty()) h->Draw();
            if (h && h->GetEntries() >= 10){
                h->SetTitle((title + ";Angle (degrees);Count").c_str());
                const PFit uf = UnbinnedP(cols[1], aMean, j ? kBlockPhiLo : 0, j ? kBlockPhiHi : 360); // unbinned fit of the phi of every proton (FitCore.h), the histogram only shows it
                AddUnbinnedFit(h, "fPhiCosCut", uf, aMean, !path.empty());
                h->SetMinimum(0); h->SetMaximum(1.1 * h->GetMaximum());
                *pol = uf.p; *polErr = uf.err;
                if (!path.empty()){
                    h->SetName(StatsTitle(h).c_str()); // stats box title; h is deleted below
                    canvas->Update();
                    if (auto* st = (TPaveStats*)h->FindObject("stats")){st->SetX1NDC(0.76); st->SetX2NDC(0.995); st->SetY1NDC(0.55); st->SetY2NDC(0.83);}
                    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
                    tx.DrawLatex(0.13, 0.90, "Unbinned fit  N(#phi) #propto 1 + P #upoint #LTA#GT cos#phi (every proton, no bins)");
                    tx.DrawLatex(0.13, 0.865, Form("#LTA#GT = %.4f (fixed, model A at the %s theta),   P = %.2f #pm %.2f", aMean, j ? "recon block" : "recon", *pol, *polErr));
                    canvas->Modified(); canvas->Update();
                    canvas->SaveAs(path.c_str());
                    WritePlotNumbers(h, path, fitFilesStarted); // stats.csv and fit_values.csv (N0, P, A) next to the pdf
                }
            }
            if (h) delete h;
            canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
        };
        TTree* ceAllTree = OpenCut("cuts/cut_CEall.root", "ChargeExchange_values"); // every charge exchange the sim made (CutEvents.cpp), counted once
        if (!ceAllTree) return;
        const Long64_t nCEallZdcAll = ceAllTree->GetEntries();
        if (RunPart(2) || RunPart(5)){
        if (RunPart(2) && RunCut(0)) for (const auto& e : cutEnergies) for (const bool ke80 : {false, true}) for (const bool noTheta : {false, true}){ // the uncut samples (1_hitZDC with the theta window, 0_cuts without) are in every ranking as the baseline
            if (noTheta && ke80) continue; // 0_cuts has no KE cut
            const std::string q = (noTheta ? zdcRecoAll : zdcBase) + (ke80 ? " && " + e.branch + " > " + std::to_string(initialKE * 0.80) : "");
            double fom[2], pol[2], polErr[2]; fomStats(events2, q, fom);
            for (int j = 0; j < 2; j++) phiFit(events2, q, j, "", "", &pol[j], &polErr[j]);
            const std::string qPart = (noTheta ? zdcRecoPartAll : zdcBasePart) + (ke80 ? " && " + e.branch + " > " + std::to_string(initialKE * 0.80) : "");
            rankRows.push_back({noTheta ? "0_cuts (no cut)" : "1_hitZDC (no cut)", ke80 ? "cut and ZDC KE > 80%" : "cut", double(CountMT(events2, P3(q))), double(CountMT(eventsAllEv, P3(q))), double(nCEallZdcAll),
                                nCEreconTotal, {fom[0], fom[1]}, {pol[0], pol[1]}, {polErr[0], polErr[1]}, double(CountMT(eventsAllPart, P3(qPart))), noTheta ? nNonProtAll : nNonProtTh, noTheta ? nNonPartAll : nNonPartTh});
        }
        for (const auto& sc : sampleCuts) for (const auto& e : cutEnergies) for (const bool ke80 : {false, true}){ // 0_KE_hitZDC (the cut and hitting the ZDC) and 0_KE_KECUT80 (also the 80% KE cut)
            if (!RunPart(2) || !RunCut(int(&sc - &sampleCuts[0]) + 1)) continue;
            P2Log("cut " + sc.dir + (ke80 ? " ke80" : " nokecut") + " start");
            const std::string folder_res = folder_KEcuts + sc.dir + "/1_additionalinfo/";
            const double cutVal = initialKE * 0.80; // the 80% KE cut of 2_KECUT80
            const std::string qCut = zdcBase + " && " + sc.expr + (ke80 ? " && " + e.branch + " > " + std::to_string(cutVal) : "");
            const long long nCEq = CountMT(events2, P3(qCut)), nAllq = CountMT(eventsAllEv, P3(qCut));
            { // cut_ratios.csv in the cut folder, one row without and one with the 80% KE cut: CE recon protons in the cut (no ZDC requirement without the 80% KE cut) over all recon protons in the cut,
              // over all CE protons the sim made (no reconstruction or ZDC requirement), and over all CE recon protons that hit ZDC before the cut
                const long long nCEallZdc = nCEallZdcAll; // every charge exchange the sim made, one CE proton each
                gSystem->mkdir((folder_KEcuts + sc.dir).c_str(), true);
                if (!ke80){ // cuts.txt in the cut folder: states the cuts of this folder, straight from sampleCuts so it cannot drift from the cut
                    FILE* tf = fopen((folder_KEcuts + sc.dir + "/cuts.txt").c_str(), "w");
                    const char* zdcNote = sc.zdcKE80 ? "Sample: AnglesRaw (cuts/cut_angles_all.root), the fitted protons that are reconstructed ZDC protons with ZDC KE > 80% of initial KE (in the cut file)."
                                                     : "Sample: AnglesRaw (cuts/cut_angles_all.root), every proton with a fitted track. No ZDC hit required\n"
                                                       "(the \"cut and ZDC KE > 80%\" row of cut_ratios.csv and the 0_KE_KECUT80 plots also need a ZDC hit and ZDC KE > 80% of initial KE).";
                    fprintf(tf, "%s: %s\n\n%s\n\n"
                                "Cuts applied (all must pass):\n  baseline of every 6_cuts folder: recon theta and recon block theta 0 to 1 deg (%s)\n  this folder: %s\n\n"
                                "rawVal columns: 1 resDRms (mm), 6 nProtonsEv, 7 nZdcCharged, 9 recon theta (deg), 10 nHits, 17 nTrackerHitsEv, 18 nTrackerTracksEv\nCut file: cuts/%s.root (CutPlots.cpp (samples))\n",
                            sc.dir.c_str(), sc.label.c_str(), zdcNote, thetaBase.c_str(), sc.expr.c_str(), sc.file.c_str());
                    fclose(tf);
                }
                FILE* rf = fopen((folder_KEcuts + sc.dir + "/cut_ratios.csv").c_str(), ke80 ? "a" : "w");
                if (!ke80) fprintf(rf, "sample,CE_recon_proton_in_cut,recon_proton_in_cut,CE_recon_proton_in_cut/recon_proton_in_cut,CE_protons_total,CE_recon_before_cut,CE_recon_proton_in_cut/CE_recon_before_cut,nonCE_recon_proton_in_cut,nonCE_recon_proton_in_cut/recon_proton_in_cut,FOM_theta,FOM_thetaBlock,P_phiRe,P_phiRe_err,P_phiReD,P_phiReD_err,recon_particle_in_cut,CE_recon_proton_in_cut/recon_particle_in_cut,nonCE_recon_particle_in_cut/recon_particle_in_cut,nonCE_recon_proton_in_cut/nonCE_recon_proton_out_of_cut,nonCE_recon_particle_in_cut/nonCE_recon_particle_out_of_cut,nonCE_recon_proton_in_cut/nonCE_recon_proton_before_cut,nonCE_recon_particle_in_cut/nonCE_recon_particle_before_cut,CE_recon_proton_in_cut/sqrt(nonCE_recon_proton_in_cut),CE_recon_proton_in_cut/sqrt(nonCE_recon_particle_in_cut)\n");
                // without the 80% KE cut no ZDC hit is required (AnglesRaw); the ZDC KE cut needs the ZDC, so that row is CE recon protons that hit ZDC
                TTree* tS = ke80 ? events2 : cutTree(sc); // the cut file holds only the rows that pass the cut
                const std::string qS = ke80 ? qCut : "PDG == 2212", qSCE = ke80 ? qCut : rawCE; // qS: the recon protons of the cut
                const long long nCEc = CountMT(tS, P3(qSCE)), nAllc = ke80 ? nAllq : CountMT(tS, P3(qS));
                const long long nPartc = ke80 ? CountMT(eventsAllPart, P3(zdcBasePart + " && " + sc.expr + " && " + e.branch + " > " + std::to_string(cutVal))) : CountMT(tS, P3("1")); // every recon particle of the cut
                const double nonProtOut = nNonProtTh - double(nAllc - nCEc), nonPartOut = nNonPartTh - double(nPartc - nCEc); // the non-CE protons / particles the cut removes
                P2Log("  counts done");
                double fom[2], pol[2], polErr[2]; fomStats(tS, qSCE, fom);
                P2Log("  fomStats done");
                for (int j = 0; j < 2; j++) phiFit(tS, qSCE, j, "", "", &pol[j], &polErr[j]);
                P2Log("  phiFit done");
                fprintf(rf, "%s,%lld,%lld,%.6f,%lld,%.0f,%.6f,%lld,%.6f,%.6g,%.6g,%.5f,%.5f,%.5f,%.5f,%lld,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", ke80 ? "cut and ZDC KE > 80%" : "cut", nCEc, nAllc, double(nCEc) / std::max<long long>(1, nAllc),
                        nCEallZdc, nCEreconTotal, nCEc / std::max(1.0, nCEreconTotal),
                        nAllc - nCEc, double(nAllc - nCEc) / std::max<long long>(1, nAllc), fom[0], fom[1], pol[0], polErr[0], pol[1], polErr[1],
                        nPartc, double(nCEc) / std::max<long long>(1, nPartc), double(nPartc - nCEc) / std::max<long long>(1, nPartc), (nAllc - nCEc) / std::max(1.0, nonProtOut), (nPartc - nCEc) / std::max(1.0, nonPartOut),
                        (nAllc - nCEc) / std::max(1.0, nNonProtTh), (nPartc - nCEc) / std::max(1.0, nNonPartTh), nCEc / std::sqrt(std::max(1.0, double(nAllc - nCEc))), nCEc / std::sqrt(std::max(1.0, double(nPartc - nCEc))));
                fclose(rf);
                rankRows.push_back({sc.dir, ke80 ? "cut and ZDC KE > 80%" : "cut", double(nCEc), double(nAllc), double(nCEallZdc), nCEreconTotal, {fom[0], fom[1]}, {pol[0], pol[1]}, {polErr[0], polErr[1]}, double(nPartc), nNonProtTh, nNonPartTh});
            }
            P2Log("  rank row done");
            if (ke80) for (const auto& r : keRanges) for (const bool lin : {false, true}){ // the ZDC KE plot only with the KE cut: the cut alone needs no ZDC hit. lin: the same plot with a linear y axis, saved as ..._linear.pdf
                drawPlot(events2, e.branch + ">>hResCE" + binSpec(r), qCut, e.label + " of Protons that hit ZDC, " + sc.label + (ke80 ? ", ZDC KE > 80%" : "") + ", " + kThetaNote + ", CE recon protons vs all recon protons;" + e.label + " (MeV);Count", folder_res + "0_KE_KECUT80_" + e.dir + r.tag + (lin ? "_linear" : "") + ".pdf", "", [&](TH1* h){
                    h->SetLineColor(kRed);
                    TH1* hf = HistMT(eventsAllEv, P3(e.branch + ">>hResAll" + binSpec(r)), P3(qCut)); hf->Draw("same");
                    hf->SetLineColor(kBlue);
                    const std::string qCutPart = zdcBasePart + " && " + sc.expr + (ke80 ? " && " + e.branch + " > " + std::to_string(cutVal) : ""); // the same cut on every recon particle (any pdg)
                    TH1* hp = HistMT(eventsAllPart, P3(e.branch + ">>hResPart" + binSpec(r)), P3(qCutPart)); hp->Draw("same");
                    hp->SetLineColor(kGreen + 2);
                    if (lin){ // per bin counts and the CE share of the recon protons and of the recon particles in the cut, beside the linear pdf
                        FILE* bf = fopen((folder_res + "0_KE_KECUT80_" + e.dir + r.tag + "_perbin.csv").c_str(), "w");
                        fprintf(bf, "KE_lo_MeV,KE_hi_MeV,CE_recon_protons,all_recon_protons,all_recon_particles,CE_over_all_protons_percent,CE_over_all_particles_percent\n");
                        for (int b = 1; b <= h->GetNbinsX(); b++){
                            const double ce = h->GetBinContent(b), pr = hf->GetBinContent(b), pa = hp->GetBinContent(b);
                            fprintf(bf, "%g,%g,%g,%g,%g,%.3f,%.3f\n", h->GetBinLowEdge(b), h->GetBinLowEdge(b) + h->GetBinWidth(b), ce, pr, pa, pr > 0 ? 100 * ce / pr : 0., pa > 0 ? 100 * ce / pa : 0.);
                        }
                        fclose(bf);
                    }
                    h->SetMinimum(lin ? 0 : 0.5); h->SetMaximum((lin ? 1.1 : 2.0) * std::max(h->GetMaximum(), std::max(hf->GetMaximum(), hp->GetMaximum()))); canvas->SetLogy(lin ? 0 : 1);
                    TLegend* leg = new TLegend(0.745, 0.44, 0.99, 0.55); leg->SetTextSize(0.012);
                    leg->AddEntry(h, Form("CE recon protons that hit ZDC%s, %lld", ke80 ? ", KE > 80%" : "", nCEq), "l");
                    leg->AddEntry(hf, Form("all recon protons that hit ZDC%s, %lld", ke80 ? ", KE > 80%" : "", nAllq), "l");
                    leg->AddEntry(hp, Form("all recon particles that hit ZDC%s, %lld", ke80 ? ", KE > 80%" : "", CountMT(eventsAllPart, P3(qCutPart))), "l");
                    leg->Draw();
                });
                canvas->SetLogy(0);
            }
        }

        const char* partEnv = std::getenv("P3_PART");
        const bool splitRank = partEnv && std::atoi(partEnv) == 2 && std::getenv("P3_CUT") && *std::getenv("P3_CUT"); // one process of part 2: its rows go to a file
        const bool aggregate = partEnv && std::atoi(partEnv) == 5; // part 5: the rows of every part 2 process
        const auto rowFile = [](int i){return "cuts/rank_rows_" + std::to_string(i) + ".tsv";};
        if (splitRank){
            FILE* rf = fopen(rowFile(std::atoi(std::getenv("P3_CUT"))).c_str(), "w");
            for (const auto& r : rankRows)
                fprintf(rf, "%s\t%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", r.cut.c_str(), r.sample.c_str(), r.nCEin, r.nAllIn, r.nCEtotal, r.nCEreconBefore,
                        r.fom[0], r.fom[1], r.pol[0], r.pol[1], r.polErr[0], r.polErr[1], r.nPartIn, r.nNonProtBefore, r.nNonPartBefore);
            fclose(rf);
        }
        if (aggregate){
            rankRows.clear();
            for (int i = 0; i <= (int)sampleCuts.size(); i++){
                std::ifstream in(rowFile(i));
                if (!in){Error("Plot3D", "%s missing: run the part 2 process of cut %d first", rowFile(i).c_str(), i); return;}
                std::string line;
                while (std::getline(in, line)){
                    std::vector<std::string> f;
                    std::stringstream ss(line);
                    for (std::string x; std::getline(ss, x, '\t');) f.push_back(x);
                    if (f.size() != 15) continue;
                    RankRow r{f[0], f[1], std::stod(f[2]), std::stod(f[3]), std::stod(f[4]), std::stod(f[5]), {std::stod(f[6]), std::stod(f[7])}, {std::stod(f[8]), std::stod(f[9])}, {std::stod(f[10]), std::stod(f[11])}, std::stod(f[12]), std::stod(f[13]), std::stod(f[14])};
                    rankRows.push_back(r);
                }
            }
        }
        if (!splitRank){ // the cut folders ranked by each column of cut_ratios.csv (best first, nonCE / CE lowest first), in 6_cuts: rank_<column>.csv, every column in every file
            using Get = double (*)(const RankRow&);
            struct Col {const char* name; Get v; bool lowBest;};
            const Col cols[] = {
                {"CE_recon_proton_in_cut/recon_proton_in_cut", [](const RankRow& r){return r.nCEin / std::max(1.0, r.nAllIn);}, false},
                {"CE_recon_proton_in_cut/CE_recon_before_cut", [](const RankRow& r){return r.nCEin / std::max(1.0, r.nCEreconBefore);}, false},
                {"nonCE_recon_proton_in_cut/recon_proton_in_cut", [](const RankRow& r){return (r.nAllIn - r.nCEin) / std::max(1.0, r.nAllIn);}, true},
                {"FOM_theta", [](const RankRow& r){return r.fom[0];}, false}, {"FOM_thetaBlock", [](const RankRow& r){return r.fom[1];}, false},
                {"P_phiRe", [](const RankRow& r){return r.pol[0];}, false}, {"P_phiReD", [](const RankRow& r){return r.pol[1];}, false},
                {"CE_recon_proton_in_cut/recon_particle_in_cut", [](const RankRow& r){return r.nCEin / std::max(1.0, r.nPartIn);}, false},
                {"nonCE_recon_particle_in_cut/recon_particle_in_cut", [](const RankRow& r){return (r.nPartIn - r.nCEin) / std::max(1.0, r.nPartIn);}, true},
                {"nonCE_recon_proton_in_cut/nonCE_recon_proton_out_of_cut", [](const RankRow& r){return (r.nAllIn - r.nCEin) / std::max(1.0, r.nNonProtBefore - (r.nAllIn - r.nCEin));}, true},
                {"nonCE_recon_particle_in_cut/nonCE_recon_particle_out_of_cut", [](const RankRow& r){return (r.nPartIn - r.nCEin) / std::max(1.0, r.nNonPartBefore - (r.nPartIn - r.nCEin));}, true},
                {"nonCE_recon_proton_in_cut/nonCE_recon_proton_before_cut", [](const RankRow& r){return (r.nAllIn - r.nCEin) / std::max(1.0, r.nNonProtBefore);}, true},
                {"nonCE_recon_particle_in_cut/nonCE_recon_particle_before_cut", [](const RankRow& r){return (r.nPartIn - r.nCEin) / std::max(1.0, r.nNonPartBefore);}, true},
                {"CE_recon_proton_in_cut/sqrt(nonCE_recon_proton_in_cut)", [](const RankRow& r){return r.nCEin / std::sqrt(std::max(1.0, r.nAllIn - r.nCEin));}, false},
                {"CE_recon_proton_in_cut/sqrt(nonCE_recon_particle_in_cut)", [](const RankRow& r){return r.nCEin / std::sqrt(std::max(1.0, r.nPartIn - r.nCEin));}, false}};
            const char* files[] = {"rank_CE_recon_over_recon_proton_in_cut.csv", "rank_CE_recon_over_CE_recon_before_cut.csv", "rank_nonCE_recon_over_recon_proton_in_cut.csv",
                                   "rank_FOM_theta.csv", "rank_FOM_thetaBlock.csv", "rank_P_phiRe.csv", "rank_P_phiReD.csv",
                                   "rank_CE_recon_over_recon_particle_in_cut.csv", "rank_nonCE_recon_over_recon_particle_in_cut.csv", "rank_nonCE_recon_proton_in_over_out_of_cut.csv", "rank_nonCE_recon_particle_in_over_out_of_cut.csv",
                                   "rank_nonCE_recon_proton_kept.csv", "rank_nonCE_recon_particle_kept.csv", "rank_CE_over_sqrt_nonCE_recon_proton.csv", "rank_CE_over_sqrt_nonCE_recon_particle.csv"};
            for (size_t c = 0; c < sizeof(cols) / sizeof(cols[0]); c++){
                std::vector<RankRow> sorted = rankRows;
                std::stable_sort(sorted.begin(), sorted.end(), [&](const RankRow& a, const RankRow& b){return cols[c].lowBest ? cols[c].v(a) < cols[c].v(b) : cols[c].v(a) > cols[c].v(b);});
                FILE* rk = fopen((folder_KEcuts + files[c]).c_str(), "w");
                fprintf(rk, "rank,cut,sample");
                for (const auto& col : cols) fprintf(rk, ",%s", col.name);
                fprintf(rk, ",P_phiRe_err,P_phiReD_err,CE_recon_proton_in_cut,recon_proton_in_cut,CE_protons_total,CE_recon_before_cut,nonCE_recon_proton_in_cut,recon_particle_in_cut\n");
                for (size_t i = 0; i < sorted.size(); i++){
                    const RankRow& r = sorted[i];
                    fprintf(rk, "%zu,%s,%s", i + 1, r.cut.c_str(), r.sample.c_str());
                    for (const auto& col : cols) fprintf(rk, ",%.6g", col.v(r));
                    fprintf(rk, ",%.5f,%.5f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f\n", r.polErr[0], r.polErr[1], r.nCEin, r.nAllIn, r.nCEtotal, r.nCEreconBefore, r.nAllIn - r.nCEin, r.nPartIn);
                }
                fclose(rk);
            }
        }

        } // part 2

        // phi cos fits in every cut folder, like 3_ReconCalcVals/2_phiRe_cosfit.pdf, CE recon protons and all recon protons of the folder on the same graph, at their recon phi and recon block phi:
        // N(phi) = N0 (1 + P <A> cos(phi)), N0 = entries / 50 bins fixed, <A> = mean model A_y at the recon (block) theta of those protons fixed, only P fitted
        if (RunPart(4)) {
            struct PhiSet {std::string dir, expr; TTree* t; std::string base, basePart, note;}; // base: the sample the cut applies to (hit ZDC and reconstructed, or any fitted track for the cuts)
            std::vector<PhiSet> phiSets = {{"1_hitZDC", "", events2, zdcBase, zdcBasePart, kThetaNote}, {"2_KECUT80", " && fDEnergy > " + std::to_string(initialKE * 0.80), events2, zdcBase, zdcBasePart, kThetaNote}, {"0_cuts/1_noCuts", "", events2, zdcRecoAll, zdcRecoPartAll, kNoCutNote}, {"0_cuts/2_theta0_1", "", events2, zdcBase, zdcBasePart, kThetaNote}};
            for (const auto& sc : sampleCuts) phiSets.push_back({sc.dir, "", cutTree(sc), rawCE, "", ""});
            for (const int pct : {85, 90, 95}){ // KE scan of 9_megaCut3/1_1particlecut (CutPlots.cpp (samples)): P3_CUT 12 to 14
                TTree* t = OpenCut(Form("cuts/cut_megacut3_1_ke%d.root", pct), "AnglesRaw");
                if (t && !t->GetAlias("raw_0")) Plot3DRawAliases(t);
                gROOT->cd();
                if (t) phiSets.push_back({Form("9_megaCut3/1_1particlecut_KE%d", pct), "", t, rawCE, "", ""});
            }
            for (const auto& ps : phiSets){
                if (const int pi = int(&ps - &phiSets[0]); !RunCut(pi < 4 ? 0 : pi - 3) || (pi < 4 && !RunSub(pi))) continue; // P3_CUT: 0 the first four sets (P3_SUB picks one), i the sampleCuts[i - 1] set
                gSystem->mkdir((folder_KEcuts + ps.dir + "/2_phi_cosfit").c_str(), true);
                const bool isZdc = ps.t == events2; // the ZDC folders hold the protons that hit the ZDC, the cut folders the fitted protons that pass the cut
                struct Samp {TTree* t; std::string q, name; int color;}; // CE recon protons (red), all recon protons (pdg 2212, blue) and all recon particles (any pdg, green) of the folder, drawn together
                const Samp samps[3] = {{ps.t, ps.base + ps.expr, "CE recon protons", kRed},
                                       {isZdc ? eventsAllEv : ps.t, (isZdc ? ps.base : std::string("PDG == 2212")) + ps.expr, "all recon protons", kBlue},
                                       {isZdc ? eventsAllPart : ps.t, (isZdc ? ps.basePart : std::string("1")) + ps.expr, "all recon particles", kGreen + 2}};
                const std::string cutTxt = [&]() -> std::string { // every cut the folder applies, for the title
                    if (isZdc) return "hit ZDC, " + ps.note + (ps.expr.empty() ? ", no KE cut" : ", ZDC KE > 80%");
                    for (const auto& c : sampleCuts) if (c.dir == ps.dir) return c.label + (c.zdcKE80 ? ", hit ZDC" : ", no ZDC requirement");
                    if (const size_t k = ps.dir.rfind("_KE"); k != std::string::npos) return "9_megaCut3 1 particle cut, hit ZDC, ZDC KE > " + ps.dir.substr(k + 3) + "%";
                    return ps.dir;
                }();
                const std::string where = isZdc ? " that hit ZDC" : " in the cut", dirP = folder_KEcuts + ps.dir + "/2_phi_cosfit/";
                const double apScalePhi = ReadSettings().at("APScale"), pTrue = ReadSettings().at("PPValue");
                for (int j = 0; j < 2; j++){ // j 0: phiRe at thetaRe, 1: phiReD at thetaReD
                    const int kT = j ? 11 : 9, kPh = j ? 15 : 14;
                    const std::string rb = j ? "recon block" : "recon", stem = j ? "phiReD" : "phiRe";
                    const std::string dirJ = dirP + (j ? "2_ReconD/" : "1_Recon/"); // recon and recon block plots in a subfolder each
                    gSystem->mkdir(dirJ.c_str(), true);
                    std::vector<double> th[3], ph[3]; // theta and phi of every proton of the sample
                    StartScan scanned[3]; // per sample, the phi histograms and A sums for every start theta, from one pass over the protons (ScanStartTheta, Plot3DCols.h), and the mean A
                    for (int s = 0; s < 3; s++){
                        std::string qs = samps[s].q + Form(" && rawOk(fEvent, TrackID, %d) && rawOk(fEvent, TrackID, %d)", kT, kPh);
                        if (j) qs += std::string(Form(" && rawVal(fEvent, TrackID, %d) >= %g && rawVal(fEvent, TrackID, %d) <= %g", kPh, kBlockPhiLo, kPh, kBlockPhiHi)); // block phi: no spike at 0 (block theta 0 to 1 deg is in the baseline)
                        const std::string colT = Form("rawVal(fEvent, TrackID, %d)", kT), colPh = Form("rawVal(fEvent, TrackID, %d)", kPh);
                        // the rows only feed histograms and sums, so a multi-cluster tree (a cut file of AnglesRaw) is read on all threads in any order
                        auto c = MultiCluster(samps[s].t) ? ColumnsUnordered(samps[s].t, {P3(colT), P3(colPh)}, P3(qs)) : Columns(samps[s].t, {P3(colT), P3(colPh)}, P3(qs));
                        th[s] = std::move(c[0]); ph[s] = std::move(c[1]);
                        scanned[s] = ScanStartTheta(th[s], ph[s], apScalePhi, j ? kBlockPhiLo : 0, j ? kBlockPhiHi : 360);
                    }
                    // 1: phi distribution and cos fit of the three samples on one plot, N(phi) = N0 (1 + P <A> cos(phi)), N0 = entries / 50 bins fixed, <A> = mean model A of the sample fixed, only P fitted
                    const std::string allNames = "CE recon protons, all recon protons and all recon particles";
                    for (const bool ceOnly : {false, true}){ // the CE recon protons get a plot of their own (_CE) with their own y range; the other plot has the all recon protons and all recon particles
                        const int sLo = ceOnly ? 0 : 1, sHi = ceOnly ? 1 : 3;
                        TH1D* hP[3] = {nullptr, nullptr, nullptr};
                        double pol[3] = {0, 0, 0}, polErr[3] = {0, 0, 0}, aMean[3] = {0, 0, 0}, hMax = 0;
                        for (int s = sLo; s < sHi; s++){
                            delete gDirectory->Get(Form("hPhiOv%d", s));
                            hP[s] = new TH1D(Form("hPhiOv%d", s), "", 50, 0, 360);
                            for (double x : ph[s]) hP[s]->Fill(x);
                            aMean[s] = scanned[s].aAll;
                            hP[s]->SetLineColor(samps[s].color); hP[s]->SetLineWidth(2);
                            if (hP[s]->GetEntries() >= 10){
                                const PFit uf = UnbinnedP(ph[s], aMean[s], j ? kBlockPhiLo : 0, j ? kBlockPhiHi : 360); // unbinned fit of every proton's phi (FitCore.h)
                                AddUnbinnedFit(hP[s], Form("fPhiOv%d", s), uf, aMean[s], false); // the fit stays in the histogram for stats.csv and fit_values.csv, drawn below
                                pol[s] = uf.p; polErr[s] = uf.err;
                            }
                            hMax = std::max(hMax, hP[s]->GetMaximum());
                        }
                        canvas->Clear(); canvas->cd();
                        gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.36); gPad->SetBottomMargin(0.13); // the text lines sit above the frame
                        for (int s = sLo; s < sHi; s++){
                            hP[s]->SetStats(0); hP[s]->SetMinimum(0); hP[s]->SetMaximum(1.1 * hMax);
                            hP[s]->SetTitle(Form("%s%s, %s, %s phi distribution;Angle (degrees);Count", (ceOnly ? "CE recon protons" : "all recon protons and all recon particles"), where.c_str(), cutTxt.c_str(), rb.c_str()));
                            hP[s]->Draw(s > sLo ? "hist same" : "hist");
                        }
                        for (int s = sLo; s < sHi; s++) if (hP[s]->GetEntries() >= 10){ // the fit curve of each set in its color, a function of its own: a clone of the one in the histogram crashed when the canvas was cleared
                            TF1* fd = new TF1(Form("fPhiDraw%d", s), "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                            fd->SetParameters(hP[s]->GetEntries() / 50, pol[s], aMean[s]); fd->SetLineColor(samps[s].color); fd->SetLineStyle(2); fd->SetLineWidth(3); fd->Draw("same");
                        }
                        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
                        tx.DrawLatex(0.18, 0.90, "Unbinned fit  N(#phi) #propto 1 + P #upoint #LTA#GT cos#phi, dashed = fit, #LTA#GT fixed");
                        for (int s = sLo; s < sHi; s++){
                            tx.SetTextColor(samps[s].color);
                            tx.DrawLatex(0.18, 0.855 - 0.045 * (s - sLo), Form("%s: %.0f entries, #LTA#GT = %.4f (%s theta),  P = %.3g #pm %.2g", samps[s].name.c_str(), hP[s]->GetEntries(), aMean[s], rb.c_str(), pol[s], polErr[s]));
                        }
                        tx.SetTextColor(kBlack);
                        if (j) tx.DrawLatex(0.18, 0.72, BlockPhiNote());
                        const std::string path = dirJ + "1_" + stem + "_cosfit" + (ceOnly ? "_CE" : "") + ".pdf";
                        canvas->Modified(); canvas->Update();
                        canvas->SaveAs(path.c_str());
                        for (int s = sLo; s < sHi; s++){
                            hP[s]->SetName(StatsTitle(hP[s], samps[s].name + where).c_str());
                            WritePlotNumbers(hP[s], path, fitFilesStarted); // stats.csv and fit_values.csv (N0, P, A) next to the pdf, a row per sample
                            delete hP[s];
                        }
                        canvas->Clear();
                    }
                    // 2, 3: P of the same fit for the protons from a start theta up to 1 deg (recon theta, or recon block theta for j 1), start in 0.01 deg steps, and its percent error to the true P (PPValue, dotted line)
                    struct Scan {std::vector<double> x, y, ey, kept;}; // kept: % of all CE recon protons (CE_recon_before_cut), of the sample for all recon protons, that are in the sample with theta from the start value on
                    Scan scan[3];
                    struct Page {TH1D* h; double start, P, eP, A;}; // the phi distribution and fit behind the scan point every 0.1 deg, for 4_..._cosfit_by_start_theta
                    std::vector<Page> pages[3];
                    for (int s = 0; s < 3; s++){
                        const size_t n = th[s].size();
                        for (int b = 0; b < 100; b++){
                            const double start = b * 0.01;
                            TH1D hS("hPstart", "", 50, 0, 360); // the protons with start <= theta <= 1 deg, from the prepared bins (under and overflow included)
                            double nIn = 0;
                            for (int j = 0; j < kScanBins; j++){hS.SetBinContent(j, scanned[s].phi[b][j]); nIn += scanned[s].phi[b][j];}
                            hS.SetEntries(nIn);
                            const double aSum = scanned[s].aSum[b];
                            if (hS.GetEntries() < 10) break; // the minimum the cosfit plot fits with
                            const PFit uf = scanned[s].fit[b]; // unbinned P of the protons of this start (ScanStartTheta, Plot3DCols.h)
                            if (b % 10 == 0){TH1D* hc = (TH1D*)hS.Clone(Form("hPage%d_%d", s, b)); hc->SetDirectory(nullptr); pages[s].push_back({hc, start, uf.p, uf.err, aSum / hS.GetEntries()});}
                            scan[s].x.push_back(start); scan[s].y.push_back(uf.p); scan[s].ey.push_back(uf.err);
                            scan[s].kept.push_back(100.0 * hS.GetEntries() / std::max(1.0, s ? double(n) : nCEreconTotal));
                        }
                    }
                    { // 4: one pdf of many pages, a page per start theta every 0.1 deg, the three samples on each page
                        std::set<long> starts;
                        for (int s = 0; s < 3; s++) for (const auto& pg : pages[s]) starts.insert(std::lround(pg.start * 100));
                        const std::string path = dirJ + "4_" + stem + "_cosfit_by_start_theta.pdf";
                        size_t k = 0;
                        for (long key : starts){
                            const Page* pgOf[3] = {nullptr, nullptr, nullptr};
                            double hMax = 0;
                            for (int s = 0; s < 3; s++) for (const auto& pg : pages[s]) if (std::lround(pg.start * 100) == key){pgOf[s] = &pg; hMax = std::max(hMax, pg.h->GetMaximum());}
                            canvas->Clear(); canvas->SetLeftMargin(0.16); canvas->SetRightMargin(0.04); canvas->SetTopMargin(0.36); canvas->SetBottomMargin(0.13);
                            bool first = true;
                            for (int s = 0; s < 3; s++) if (pgOf[s]){
                                TH1D* h = pgOf[s]->h;
                                h->SetLineColor(samps[s].color); h->SetLineWidth(2); h->SetStats(0); h->SetMinimum(0); h->SetMaximum(1.1 * hMax);
                                h->SetTitle(Form("%s%s, %s, %s phi distribution, %s theta %.2f to 1 deg;Angle (degrees);Count", allNames.c_str(), where.c_str(), cutTxt.c_str(), rb.c_str(), rb.c_str(), key / 100.0));
                                h->Draw(first ? "hist" : "hist same"); first = false;
                            }
                            std::vector<TF1*> fds;
                            for (int s = 0; s < 3; s++) if (pgOf[s]){
                                TF1* fd = new TF1(Form("fPageDraw%d_%zu", s, k), "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                                fd->SetParameters(pgOf[s]->h->GetEntries() / 50, pgOf[s]->P, pgOf[s]->A); fd->SetLineColor(samps[s].color); fd->SetLineStyle(2); fd->SetLineWidth(3); fd->Draw("same");
                                fds.push_back(fd);
                            }
                            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.033);
                            tx.DrawLatex(0.18, 0.90, Form("start theta %.2f deg: unbinned fit N(#phi) #propto 1 + P #upoint #LTA#GT cos#phi, dashed = fit", key / 100.0));
                            for (int s = 0; s < 3; s++) if (pgOf[s]){
                                tx.SetTextColor(samps[s].color);
                                tx.DrawLatex(0.18, 0.855 - 0.045 * s, Form("%s: %.0f entries, #LTA#GT = %.4f,  P = %.3g #pm %.2g", samps[s].name.c_str(), pgOf[s]->h->GetEntries(), pgOf[s]->A, pgOf[s]->P, pgOf[s]->eP));
                            }
                            tx.SetTextColor(kBlack);
                            if (j) tx.DrawLatex(0.18, 0.72, BlockPhiNote());
                            canvas->Modified(); canvas->Update();
                            canvas->Print((path + (starts.size() > 1 ? (k == 0 ? "(" : k + 1 == starts.size() ? ")" : "") : "")).c_str());
                            for (TF1* fd : fds) delete fd;
                            k++;
                        }
                        for (int s = 0; s < 3; s++) for (auto& pg : pages[s]) delete pg.h;
                        canvas->Clear();
                    }
                    for (const int mode : {0, 1, 2, 3}){ // 2, 3: the P scans of the three samples on one axis, then the CE recon protons alone (_CE)
                        const bool pct = mode & 1, ceOnly = mode & 2;
                        bool any = false;
                        for (int s = ceOnly ? 0 : 1; s < (ceOnly ? 1 : 3); s++) if (scan[s].x.size() > 1) any = true;
                        if (!any) continue;
                        canvas->Clear(); canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
                        const double yRef = pct ? 0.0 : pTrue;
                        auto conv = [&](double P){return pct ? 100 * (P - pTrue) / pTrue : P;}; // P to the plotted value
                        double ymin = yRef, ymax = yRef; const double xMax = 1.0; // the scan always ends at 1 deg, the upper theta of the fits
                        for (int s = ceOnly ? 0 : 1; s < (ceOnly ? 1 : 3); s++) if (scan[s].x.size() > 1) for (size_t i = 0; i < scan[s].x.size(); i++){ // the error bars are cut at P = -1 and 1
                            ymin = std::min(ymin, conv(std::max(scan[s].y[i] - scan[s].ey[i], -1.0))); ymax = std::max(ymax, conv(std::min(scan[s].y[i] + scan[s].ey[i], 1.0)));
                        }
                        const double pad = 0.05 * (ymax - ymin);
                        TH2D axis("hPstartAxis", Form("Polarization fit P of %s%s, %s, %s phi, %s theta from the start value to 1 deg%s;Start theta (deg);%s", (ceOnly ? "CE recon protons" : "all recon protons and all recon particles"), where.c_str(), ps.dir.c_str(), rb.c_str(), rb.c_str(), j ? ", phi 0.5 to 359.5 deg" : "",
                                  pct ? "100 (P_{fit} - P_{true}) / P_{true} (%)" : "P"), 100, 0, xMax, 100, ymin - pad, ymax + pad);
                        axis.SetStats(0);
                        axis.Draw("axis");
                        TLine ref(0, yRef, xMax, yRef); ref.SetLineStyle(3); ref.SetLineWidth(2); ref.SetLineColor(kGray + 2); ref.Draw(); // true P, from the run's PPValue
                        TLegend leg(0.705, 0.76, 0.99, 0.90); leg.SetTextSize(0.014);
                        TGraphAsymmErrors* gs[3] = {nullptr, nullptr, nullptr};
                        const std::string path = dirJ + std::to_string(2 + pct) + "_" + stem + (pct ? "_P_vs_start_theta_pct_error" : "_P_vs_start_theta") + (ceOnly ? "_CE" : "") + ".pdf";
                        for (int s = ceOnly ? 0 : 1; s < (ceOnly ? 1 : 3); s++) if (scan[s].x.size() > 1){
                            gs[s] = new TGraphAsymmErrors();
                            for (size_t i = 0; i < scan[s].x.size(); i++){ // the error bar of P stops at -1 and 1, the physical range
                                const double y = scan[s].y[i], lo = std::max(y - scan[s].ey[i], -1.0), hi = std::min(y + scan[s].ey[i], 1.0);
                                gs[s]->SetPoint(i, scan[s].x[i], conv(y));
                                gs[s]->SetPointError(i, 0, 0, conv(y) - conv(lo), conv(hi) - conv(y));
                            }
                            gs[s]->SetLineColor(samps[s].color); gs[s]->SetMarkerColor(samps[s].color); gs[s]->SetMarkerStyle(20); gs[s]->SetMarkerSize(0.6);
                            if (gs[s]->GetN()) gs[s]->Draw("PZ same");
                            leg.AddEntry(gs[s], samps[s].name.c_str(), "lp");
                            for (size_t i = 0; i < scan[s].x.size(); i++) // the points drawn: start theta, value, error, % kept
                                WriteStatsRow(path, "plot,sample,start_theta_deg,value,error,pct_kept", Form("%s,%s,%.4g,%.6g,%.4g,%.4g", gSystem->BaseName(path.c_str()), samps[s].name.c_str(), scan[s].x[i], gs[s]->GetY()[i], pct ? 100 * scan[s].ey[i] / pTrue : scan[s].ey[i], scan[s].kept[i]), fitFilesStarted, "P_vs_start_theta.csv");
                        }
                        leg.Draw();
                        canvas->SaveAs(path.c_str());
                        for (auto* g : gs) delete g; // null for the sample not drawn
                    }
                    if (j == 1){ // 5 to 12: ratios of the counts with block theta from the start value to 1 deg ("in cut") against the start value, CE = the CE recon protons, prot = all recon protons, part = all recon particles of the folder
                        const double nonProtBefore = double(th[1].size()) - th[0].size(), nonPartBefore = double(th[2].size()) - th[0].size(); // before the start cut: the whole sample
                        struct Ratio {const char* file, *label; int color; std::function<double(double, double, double)> f;}; // f(CE, prot, part) in cut
                        const Ratio ratios[8] = {
                            {"CE_over_particle", "N_{CE,recon,in cut} / N_{recon,particle,in cut}", kGreen + 2, [&](double ce, double, double pa){return pa > 0 ? ce / pa : NAN;}},
                            {"nonCE_particle_over_particle", "N_{nonCE,recon,particle,in cut} / N_{recon,particle,in cut}", kGreen + 2, [&](double ce, double, double pa){return pa > 0 ? (pa - ce) / pa : NAN;}},
                            {"nonCE_proton_in_over_out", "N_{nonCE,recon,proton,in cut} / N_{nonCE,recon,proton,out of cut}", kBlue, [&](double ce, double pr, double){const double out = nonProtBefore - (pr - ce); return out > 0 ? (pr - ce) / out : NAN;}},
                            {"nonCE_particle_in_over_out", "N_{nonCE,recon,particle,in cut} / N_{nonCE,recon,particle,out of cut}", kGreen + 2, [&](double ce, double, double pa){const double out = nonPartBefore - (pa - ce); return out > 0 ? (pa - ce) / out : NAN;}},
                            {"nonCE_proton_kept", "N_{nonCE,recon,proton,in cut} / N_{nonCE,recon,proton,before cut}", kBlue, [&](double ce, double pr, double){return nonProtBefore > 0 ? (pr - ce) / nonProtBefore : NAN;}},
                            {"nonCE_particle_kept", "N_{nonCE,recon,particle,in cut} / N_{nonCE,recon,particle,before cut}", kGreen + 2, [&](double ce, double, double pa){return nonPartBefore > 0 ? (pa - ce) / nonPartBefore : NAN;}},
                            {"CE_over_sqrt_nonCE_proton", "N_{CE,recon,in cut} / sqrt(N_{nonCE,recon,proton,in cut})", kBlue, [&](double ce, double pr, double){return pr - ce > 0 ? ce / std::sqrt(pr - ce) : NAN;}},
                            {"CE_over_sqrt_nonCE_particle", "N_{CE,recon,in cut} / sqrt(N_{nonCE,recon,particle,in cut})", kGreen + 2, [&](double ce, double, double pa){return pa - ce > 0 ? ce / std::sqrt(pa - ce) : NAN;}}};
                        std::vector<double> nIn[3]; // counts with theta from the start value to 1 deg, per start value (0.01 deg steps)
                        for (int b = 0; b < 100; b++) for (int s = 0; s < 3; s++){
                            double n = 0; for (double t : th[s]) if (t >= b * 0.01 && t <= 1.0) n++;
                            nIn[s].push_back(n);
                        }
                        for (int k = 0; k < 8; k++){
                            TGraph g;
                            for (int b = 0; b < 100; b++){const double v = ratios[k].f(nIn[0][b], nIn[1][b], nIn[2][b]); if (std::isfinite(v)) g.SetPoint(g.GetN(), b * 0.01, v);}
                            if (g.GetN() < 2) continue;
                            const std::string path = dirJ + std::to_string(5 + k) + "_" + stem + "_" + ratios[k].file + "_vs_start_theta.pdf";
                            canvas->Clear(); canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
                            g.SetTitle((std::string(ratios[k].label) + " vs start theta, " + cutTxt + ", " + rb + " theta from the start value to 1 deg" + where + ";Start theta (deg);" + ratios[k].label).c_str());
                            g.SetMarkerStyle(20); g.SetMarkerSize(0.6); g.SetLineWidth(1); g.SetLineColor(ratios[k].color); g.SetMarkerColor(ratios[k].color);
                            g.Draw("APL"); g.GetXaxis()->SetLimits(0, 1.0);
                            TLegend leg(0.705, 0.82, 0.99, 0.90); leg.SetTextSize(0.014); leg.AddEntry(&g, ratios[k].file, "lp"); leg.Draw();
                            canvas->SaveAs(path.c_str());
                            for (int i = 0; i < g.GetN(); i++) WriteStatsRow(path, "plot,start_theta_deg,value", Form("%s,%.4g,%.6g", gSystem->BaseName(path.c_str()), g.GetX()[i], g.GetY()[i]), fitFilesStarted, "ratios_vs_start_theta.csv");
                        }
                        { // 13: the FOM of the CE recon protons with block theta from the start value to 1 deg, same sum as fomStats (A_y^2 at the 0.05 deg bin center * 0.05 deg / N_neutrons)
                            TGraph g;
                            for (int b = 0; b < 100; b++){
                                double sum = 0;
                                for (double t : th[0]) if (t >= b * 0.01 && t <= 1.0)
                                    sum += std::pow(APFromT(kBeamMomentum, TFromAngle(kBeamMomentum, (std::floor(t / 0.05) + 0.5) * 0.05 * TMath::Pi() / 180), apScaleFom), 2) * 0.05 / nNeutrons;
                                g.SetPoint(g.GetN(), b * 0.01, sum);
                            }
                            const std::string path = dirJ + "13_" + stem + "_FOM_vs_start_theta.pdf";
                            canvas->Clear(); canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
                            g.SetTitle(("FOM of CE recon protons vs start theta, " + cutTxt + ", " + rb + " theta from the start value to 1 deg" + where + ";Start theta (deg);FOM = #Sigma A_{y}^{2} #Delta#theta / N_{neutrons}").c_str());
                            g.SetMarkerStyle(20); g.SetMarkerSize(0.6); g.SetLineWidth(1); g.SetLineColor(kRed); g.SetMarkerColor(kRed);
                            g.Draw("APL"); g.GetXaxis()->SetLimits(0, 1.0);
                            canvas->SaveAs(path.c_str());
                            for (int i = 0; i < g.GetN(); i++) WriteStatsRow(path, "plot,start_theta_deg,value", Form("%s,%.4g,%.6g", gSystem->BaseName(path.c_str()), g.GetX()[i], g.GetY()[i]), fitFilesStarted, "ratios_vs_start_theta.csv");
                        }
                    }
                    canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
                }
            }
        }

        // recon protons that hit the ZDC, all KE (no cut): CE protons and every recon proton of all events
        if (RunPart(1)) for (const auto& e : cutEnergies) for (const auto& r : keRanges) for (const int var : {0, 1, 2}){ // 1_hitZDC, 0_cuts/1_noCuts (hit ZDC only), 0_cuts/2_theta0_1 (the same sample as 1_hitZDC)
            const bool noTheta = var == 1;
            const std::string hitDir = var == 0 ? "1_hitZDC/" : var == 1 ? "0_cuts/1_noCuts/" : "0_cuts/2_theta0_1/", hb = noTheta ? zdcRecoAll : zdcBase, hbp = noTheta ? zdcRecoPartAll : zdcBasePart, hnote = noTheta ? kNoCutNote : kThetaNote;
            std::string folder_hit = folder_KEcuts + hitDir;
            gSystem->mkdir(folder_hit.c_str(), true);
            drawPlot(events2, e.branch + ">>hHitCE" + binSpec(r), hb, e.label + " of Protons that hit the ZDC, CE recon protons that hit ZDC vs all events recon, " + hnote + ", no KE cut;" + e.label + " (MeV);Count", folder_hit + "1_real_vs_recon_" + e.dir + r.tag + ".pdf", "", [&](TH1* h){
                h->SetLineColor(kRed);
                TH1* hr = HistMT(eventsAllEv, P3(e.branch + ">>hHitRecon" + binSpec(r)), P3(hb)); hr->Draw("same");
                hr->SetLineColor(kBlue);
                TH1* hpr = HistMT(eventsAllPart, P3(e.branch + ">>hHitPart" + binSpec(r)), P3(hbp)); hpr->Draw("same"); // every recon particle (any pdg) that hit the ZDC
                hpr->SetLineColor(kGreen + 2);
                h->SetMinimum(0.5); h->SetMaximum(2.0 * std::max(h->GetMaximum(), std::max(hr->GetMaximum(), hpr->GetMaximum()))); canvas->SetLogy(1);
                TLegend* leg = new TLegend(0.745, 0.44, 0.99, 0.55); leg->SetTextSize(0.015);
                leg->AddEntry(h, Form("CE recon protons that hit ZDC, %lld", (long long)h->GetEntries()), "l");
                leg->AddEntry(hr, Form("all recon protons that hit ZDC, %lld", (long long)hr->GetEntries()), "l");
                leg->AddEntry(hpr, Form("all recon particles that hit ZDC, %lld", (long long)hpr->GetEntries()), "l");
                leg->Draw();
            });
            canvas->SetLogy(0);
        }

        // 1_additionalinfo in every KE cut folder: the cut candidates (6_cuts/3_ResCuts to 5_nHits) of the reconstructed protons that hit the ZDC and pass that KE cut
        // (1_hitZDC: no KE cut), reconstructed CE protons vs all reconstructed protons of all events. Nothing is cut on these variables.
        struct InfoVar {std::string name, folder, label; int k; double lo, hi; int nBins;}; // lo >= hi: range from the 1 to 99% quantiles of all fitted protons
        const std::vector<InfoVar> kInfoVars = {
            {"resDRms", "3_ResCuts", "residual, block position (mm)", 1, 0, 0, 100},
            {"nParticlesEv", "4_nParticles", "reconstructed particles in the event", 6, -0.5, 19.5, 20},
            {"nHits", "5_nHits", "tracker hits", 10, -0.5, 10.5, 11},
            {"nZdcCharged", "8_nZdcParticles", "charged particles at the ZDC in the event", 7, -0.5, 9.5, 10},
            {"nChargedTrackerHitsEv", "10_nChargedTrackerHits", "charged tracker plane hits in the event", 17, -0.5, 40.5, 41}, // all charged particles of the event added up
            {"thetaBlock", "7_theta", "recon block theta (degrees)", 11, 0, 2.3, 100},
            {"impactD", "13_impactD", "impact to the beam axis, block position (mm)", 3, 0, 0, 100}, // plots only, no cut on these two yet
            {"r2D", "14_r2D", "R^{2} of the block position line fit", 5, 0, 1.02, 102}};
        // 1st and 99th percentile of every fitted proton's value k (the raw_k column of the AnglesRaw cache where ok_k is 1), once per variable
        // kept in cuts/info_quantiles.tsv (k, 1%, 99%) while it is newer than cut_angles_all.root, so the passes over AnglesRaw run once per run folder, not in every Plot3D process
        std::map<int, std::pair<double, double>> quantileCache;
        bool quantilesRead = false;
        auto quantiles = [&](int k){
            const char* qFile = "cuts/info_quantiles.tsv";
            if (!quantilesRead){
                quantilesRead = true;
                Long_t id, size, flags, mtQ = 0, mtA = 1;
                if (!gSystem->GetPathInfo(qFile, &id, &size, &flags, &mtQ) && !gSystem->GetPathInfo("cuts/cut_angles_all.root", &id, &size, &flags, &mtA) && mtQ >= mtA){
                    std::ifstream in(qFile);
                    int kk; double a, b;
                    while (in >> kk >> a >> b) quantileCache[kk] = {a, b};
                }
            }
            const auto hit = quantileCache.find(k);
            if (hit != quantileCache.end()) return hit->second;
            std::vector<double> vals = Columns(rawT, {Form("raw_%d", k)}, Form("ok_%d==1", k))[0];
            std::pair<double, double> q{0.0, 0.0};
            if (!vals.empty()){ // nth_element, the two values only, in place of a sort of all of them
                const size_t i1 = size_t(0.01 * vals.size()), i99 = std::min(vals.size() - 1, size_t(0.99 * vals.size()));
                std::nth_element(vals.begin(), vals.begin() + i1, vals.end());
                q.first = vals[i1];
                std::nth_element(vals.begin() + i1, vals.begin() + i99, vals.end());
                q.second = vals[i99];
            }
            quantileCache[k] = q;
            const std::string tmp = std::string(qFile) + "." + std::to_string(getpid());
            if (FILE* qf = fopen(tmp.c_str(), "w")){ // written to a temporary name and renamed, two processes never read a half file
                for (const auto& kv : quantileCache) fprintf(qf, "%d\t%.17g\t%.17g\n", kv.first, kv.second.first, kv.second.second);
                fclose(qf);
                gSystem->Rename(tmp.c_str(), qFile);
            }
            return q;
        };
        if (RunPart(3) && rawT && rawT->GetAlias("raw_0") && rawT->GetEntries() > 0){
            std::vector<std::pair<std::string, double>> infoCuts = {{"1_hitZDC", -1.0}, {"0_cuts/1_noCuts", -1.0}, {"0_cuts/2_theta0_1", -1.0}}; // (folder, KE cut in MeV, -1 none)
            for (int pct : {80}) infoCuts.push_back({"2_KECUT" + std::to_string(pct), initialKE * pct / 100.0});
            for (const auto& sc : sampleCuts) infoCuts.push_back({sc.dir, -1.0}); // no KE cut: what each cut does to every variable
            auto lapStart = std::chrono::steady_clock::now();
            const auto lap = [&](const std::string& what){ // seconds since the last lap
                const auto now = std::chrono::steady_clock::now();
                printf(">>> 1_additionalinfo %s took %.1f s\n", what.c_str(), std::chrono::duration<double>(now - lapStart).count()); fflush(stdout);
                lapStart = now;
            };
            for (const auto& ic : infoCuts){
                if (const int ii = int(&ic - &infoCuts[0]); !RunCut(ii < 4 ? 0 : ii - 3) || (ii < 4 && !RunSub(ii))) continue; // P3_CUT: 0 the first four folders (P3_SUB picks one), i the sampleCuts[i - 1] folder
                lapStart = std::chrono::steady_clock::now();
                const SampleCut* sc = nullptr;
                for (const auto& c : sampleCuts) if (c.dir == ic.first) sc = &c;
                const bool nz = sc != nullptr; // the cuts' plots use every fitted proton (AnglesRaw), the others the protons that hit the ZDC
                TTree *tCE = nz ? cutTree(*sc) : events2, *tAll = nz ? tCE : eventsAllEv, *tPart = nz ? tCE : eventsAllPart; // tPart: every recon particle (any pdg)
                const std::string ke = ic.second < 0 ? "" : " && fDEnergy > " + std::to_string(ic.second);
                const bool noTheta = ic.first == "0_cuts/1_noCuts"; // only the ZDC hit and a recon track
                const bool inCuts0 = ic.first.rfind("0_cuts/", 0) == 0; // its plots sit in the folder itself, not in a 1_additionalinfo
                const std::string zb = noTheta ? zdcRecoAll : zdcBase, zbp = noTheta ? zdcRecoPartAll : zdcBasePart, zNote = noTheta ? kNoCutNote : kThetaNote;
                const std::string bCE = nz ? rawCE : zb, bAll = nz ? "PDG == 2212" : zb, bPart = nz ? "1" : zbp, hitTxt = nz ? "" : " that hit ZDC"; // all recon protons are the pdg 2212 ones
                // every histogram of the folder is booked first and filled in one event loop per tree (one for the three samples of a cut, they share the cut file), not one loop per plot
                std::vector<std::pair<std::string, std::string>> reqCE, reqAll, reqPart; // (Draw string, cut) of the CE recon protons, of all recon protons and of all recon particles
                std::vector<int> slot(kInfoVars.size(), -1); // position of the variable in reqCE and reqAll
                std::vector<std::pair<double, double>> range(kInfoVars.size());
                for (size_t iv = 0; iv < kInfoVars.size(); iv++){
                    const auto& v = kInfoVars[iv];
                    const std::string k = std::to_string(v.k);
                    slot[iv] = (int)reqCE.size();
                    double lo = v.lo, hi = v.hi;
                    if (lo >= hi){
                        const auto [q1, q99] = quantiles(v.k); (void)q1;
                        lo = 0.0; hi = std::max(1e-6, 1.2 * q99);
                    }
                    range[iv] = {lo, hi};
                    const std::string sel = ke + " && rawOk(fEvent, TrackID, " + k + ")",
                        exprPre = "rawVal(fEvent, TrackID, " + k + ")>>hInfo", spec = Form("(%d,%g,%g)", v.nBins, lo, hi);
                    reqCE.push_back({P3(exprPre + "CE" + k + spec), P3(bCE + sel)});
                    reqAll.push_back({P3(exprPre + "Recon" + k + spec), P3(bAll + sel)});
                    reqPart.push_back({P3(exprPre + "Part" + k + spec), P3(bPart + sel)});
                }
                lap(ic.first + ": quantiles and booking the requests");
                std::vector<TH1*> hCEs, hAlls, hParts;
                if (tCE == tAll){ // a cut folder: the three samples are rows of one cut file, one event loop
                    std::vector<std::pair<std::string, std::string>> both = reqCE;
                    both.insert(both.end(), reqAll.begin(), reqAll.end());
                    both.insert(both.end(), reqPart.begin(), reqPart.end());
                    const std::vector<TH1*> hs = HistsMT(tCE, both);
                    hCEs.assign(hs.begin(), hs.begin() + reqCE.size()); hAlls.assign(hs.begin() + reqCE.size(), hs.begin() + reqCE.size() + reqAll.size()); hParts.assign(hs.begin() + reqCE.size() + reqAll.size(), hs.end());
                } else {hCEs = HistsMT(tCE, reqCE); hAlls = HistsMT(tAll, reqAll); hParts = HistsMT(tPart, reqPart);}
                lap(ic.first + ": jit and event loop of " + std::to_string(reqCE.size() + reqAll.size() + reqPart.size()) + " histograms");
                for (size_t iv = 0; iv < kInfoVars.size(); iv++){
                    const auto& v = kInfoVars[iv];
                    if (slot[iv] < 0) continue;
                    TH1 *hCE = hCEs[slot[iv]], *hRec = hAlls[slot[iv]], *hPart = hParts[slot[iv]];
                    if (!hCE || !hRec || !hPart) continue;
                    canvas->Clear(); canvas->SetLogy(1);
                    hCE->SetLineColor(kRed); hRec->SetLineColor(kBlue); hPart->SetLineColor(kGreen + 2); hCE->SetLineWidth(2); hRec->SetLineWidth(2); hPart->SetLineWidth(2); hCE->SetStats(0); hRec->SetStats(0); hPart->SetStats(0);
                    hRec->SetTitle((v.label + ", CE recon protons vs all recon protons vs all recon particles" + (nz ? std::string(", ") : std::string(", hit ZDC, " + zNote + ", ")) + (ic.second < 0 ? std::string(sc ? "no KE cut, " + sc->label + ", no ZDC requirement" : "no KE cut") : "ZDC KE > " + std::to_string(int(ic.second)) + " MeV") + ";" + v.label + ";Count").c_str());
                    if (v.name == "nHits") for (int bin = 1; bin <= hRec->GetNbinsX(); bin++) hRec->GetXaxis()->SetBinLabel(bin, bin % 2 == 1 && bin > 1 ? Form("%d", bin - 1) : ""); // labels at 2, 4, 6, 8 hits only
                    hRec->SetMinimum(0.5); hRec->SetMaximum(2.0 * std::max(std::max(hRec->GetMaximum(), hCE->GetMaximum()), hPart->GetMaximum()));
                    hRec->Draw("hist"); hPart->Draw("hist same"); hCE->Draw("hist same");
                    TLegend* leg = new TLegend(0.745, 0.78, 0.99, 0.90); leg->SetTextSize(0.014);
                    leg->AddEntry(hCE, Form("CE recon protons (%.0f)", hCE->GetEntries()), "l");
                    leg->AddEntry(hRec, Form("all recon protons (%.0f)", hRec->GetEntries()), "l");
                    leg->AddEntry(hPart, Form("all recon particles (%.0f)", hPart->GetEntries()), "l");
                    leg->Draw();
                    // flat files named <folder>_<variable> in 1_additionalinfo of the hit ZDC, KE cut and 2 mm cut folders
                    const std::string idir = folder_KEcuts + ic.first + (sc && v.folder == sc->onFolder || inCuts0 ? "/" : "/1_additionalinfo/"); // the variable a cut is on sits beside 1_additionalinfo; 0_cuts has no such folder, its plots sit in it directly
                    gSystem->mkdir(idir.c_str(), true);
                    canvas->SaveAs((idir + v.folder + "_" + v.name + ".pdf").c_str());
                    if (noTheta && v.name == "nParticlesEv") canvas->SaveAs((folder_KEcuts + "4_nParticles/4_nParticles_nParticlesEv.pdf").c_str()); // the general count, no cut but the ZDC hit, beside the 1 and 2 particle cuts
                    canvas->SetLogy(0);
                    delete hCE; delete hRec; delete hPart;
                }
                { // 9, 10: block theta and phi resolution of the CE recon protons of the folder: recon block angle minus the true CE angle (CalcValsReDDIFF of cut_APandT.root), only CE protons have a truth
                    static std::map<std::pair<int, int>, std::pair<double, double>> blockDiff; // (fEvent, CE track id) -> (theta, phi) difference in degrees
                    static bool blockDiffLoaded = false;
                    if (!blockDiffLoaded){
                        blockDiffLoaded = true;
                        ROOT::RDataFrame dd("CalcValsReDDIFF", "cuts/cut_APandT.root");
                        auto evD = dd.Take<int>("fEvent"), tidD = dd.Take<int>("TrackID");
                        auto thD = dd.Take<double>("thetaReD"), phD = dd.Take<double>("phiReD");
                        for (size_t i = 0; i < evD->size(); i++) blockDiff[{(*evD)[i], (*tidD)[i]}] = {(*thD)[i], (*phD)[i]};
                    }
                    const auto keys = Columns(tCE, {P3("fEvent"), P3("ParentID")}, P3(bCE + ke)); // the CE track is the parent of the CE proton
                    const std::string cutNote = (nz ? sc->label + (sc->zdcKE80 ? ", hit ZDC" : ", no ZDC requirement") + ", " + kThetaNote
                                                    : std::string("hit ZDC, ") + zNote + (ic.second < 0 ? ", no KE cut" : ", ZDC KE > 80%"));
                    struct DiffPlot {std::string name, label; double lim; int nBins;};
                    const DiffPlot dps[2] = {{"thetaBlockDiff", "recon block theta - true CE theta (degrees)", 1.0, 80}, {"phiBlockDiff", "recon block phi - true CE phi (degrees)", 180.0, 72}};
                    for (int w = 0; w < 2; w++){
                        TH1D hd(Form("hBlockDiff%d", w), "", dps[w].nBins, -dps[w].lim, dps[w].lim);
                        hd.SetDirectory(nullptr);
                        for (size_t i = 0; i < keys[0].size(); i++){
                            const auto it = blockDiff.find({(int)keys[0][i], (int)keys[1][i]});
                            if (it != blockDiff.end()) hd.Fill(w ? it->second.second : it->second.first);
                        }
                        canvas->Clear(); canvas->SetLogy(0);
                        hd.SetLineColor(kRed); hd.SetLineWidth(2);
                        hd.SetTitle((dps[w].label + ", CE recon protons" + (nz ? "" : " that hit ZDC") + ", " + cutNote + ";" + dps[w].label + ";Count").c_str());
                        hd.SetStats(0);
                        if (w == 1) hd.GetXaxis()->SetNdivisions(508, kFALSE); // phi: -180 to 180 every 45 deg
                        hd.SetMinimum(0); hd.SetMaximum(1.3 * std::max(1.0, hd.GetMaximum()));
                        hd.Draw("hist");
                        TLegend* leg = new TLegend(0.745, 0.80, 0.99, 0.90); leg->SetTextSize(0.014);
                        leg->AddEntry(&hd, Form("CE recon protons (%.0f)", hd.GetEntries()), "l");
                        if (hd.GetEntries() >= 10){
                            hd.Fit("gaus", "Q0");
                            TF1* fg = hd.GetFunction("gaus");
                            TF1* fd = new TF1(Form("fBlockDiffDraw%d", w), "gaus", -dps[w].lim, dps[w].lim);
                            fd->SetParameters(fg->GetParameter(0), fg->GetParameter(1), fg->GetParameter(2)); fd->SetLineColor(kBlack); fd->SetLineWidth(2); fd->Draw("same");
                            leg->AddEntry(fd, Form("Gaussian fit, mean %.3g, sigma %.3g", fg->GetParameter(1), fg->GetParameter(2)), "l");
                        }
                        leg->Draw();
                        canvas->Modified(); canvas->Update();
                        const std::string idir = folder_KEcuts + ic.first + (inCuts0 ? "/" : "/1_additionalinfo/");
                        gSystem->mkdir(idir.c_str(), true);
                        canvas->SaveAs((idir + std::to_string(9 + w) + "_blockDiff_" + dps[w].name + ".pdf").c_str());
                    }
                }
                { // 11, 12: recon minus recon block, theta and phi, CE recon protons (red) and all recon protons (blue) of the folder
                    const std::string cutNote = (nz ? sc->label + (sc->zdcKE80 ? ", hit ZDC" : ", no ZDC requirement") + ", " + kThetaNote
                                                    : std::string("hit ZDC, ") + zNote + (ic.second < 0 ? ", no KE cut" : ", ZDC KE > 80%"));
                    struct ReDiffPlot {std::string name, label, expr, ok; double lim; int nBins;};
                    const ReDiffPlot rps[2] = {
                        {"thetaReconBlockDiff", "recon theta - recon block theta (degrees)", "rawVal(fEvent, TrackID, 9)-rawVal(fEvent, TrackID, 11)", " && rawOk(fEvent, TrackID, 9) && rawOk(fEvent, TrackID, 11)", 1.0, 80},
                        {"phiReconBlockDiff", "recon phi - recon block phi (degrees)", "rawVal(fEvent, TrackID, 14)-rawVal(fEvent, TrackID, 15)", " && rawOk(fEvent, TrackID, 14) && rawOk(fEvent, TrackID, 15)", 180.0, 72}};
                    for (int w = 0; w < 2; w++){
                        const std::string e = w ? "(" + rps[w].expr + ")" : rps[w].expr;
                        const std::string ex = w ? "(" + e + "-360*(" + e + ">180)+360*(" + e + "<-180))" : e; // phi difference into -180..180
                        TH1D hc(Form("hReDiffCE%d", w), "", rps[w].nBins, -rps[w].lim, rps[w].lim), ha(Form("hReDiffAll%d", w), "", rps[w].nBins, -rps[w].lim, rps[w].lim), hq(Form("hReDiffPart%d", w), "", rps[w].nBins, -rps[w].lim, rps[w].lim);
                        hc.SetDirectory(nullptr); ha.SetDirectory(nullptr); hq.SetDirectory(nullptr);
                        const auto colCE = Columns(tCE, {P3(ex)}, P3(bCE + ke + rps[w].ok)), colAll = Columns(tAll, {P3(ex)}, P3(bAll + ke + rps[w].ok)), colPart = Columns(tPart, {P3(ex)}, P3(bPart + ke + rps[w].ok));
                        for (double d : colCE[0]) hc.Fill(d);
                        for (double d : colAll[0]) ha.Fill(d);
                        for (double d : colPart[0]) hq.Fill(d);
                        canvas->Clear(); canvas->SetLogy(1);
                        hc.SetLineColor(kRed); ha.SetLineColor(kBlue); hq.SetLineColor(kGreen + 2); hc.SetLineWidth(2); ha.SetLineWidth(2); hq.SetLineWidth(2); hc.SetStats(0); ha.SetStats(0); hq.SetStats(0);
                        ha.SetTitle((rps[w].label + ", CE recon protons vs all recon protons vs all recon particles, " + cutNote + ";" + rps[w].label + ";Count").c_str());
                        if (w == 1) ha.GetXaxis()->SetNdivisions(508, kFALSE); // phi: -180 to 180 every 45 deg
                        ha.SetMinimum(0.5); ha.SetMaximum(2.0 * std::max(1.0, std::max(std::max(ha.GetMaximum(), hc.GetMaximum()), hq.GetMaximum()))); // log y
                        ha.Draw("hist"); hq.Draw("hist same"); hc.Draw("hist same");
                        TLegend* leg = new TLegend(0.745, 0.76, 0.99, 0.90); leg->SetTextSize(0.011);
                        leg->AddEntry(&hc, Form("CE recon protons (%.0f), mean %.3g, RMS %.3g", hc.GetEntries(), hc.GetMean(), hc.GetRMS()), "l");
                        leg->AddEntry(&ha, Form("all recon protons (%.0f), mean %.3g, RMS %.3g", ha.GetEntries(), ha.GetMean(), ha.GetRMS()), "l");
                        leg->AddEntry(&hq, Form("all recon particles (%.0f), mean %.3g, RMS %.3g", hq.GetEntries(), hq.GetMean(), hq.GetRMS()), "l");
                        leg->Draw();
                        canvas->Modified(); canvas->Update();
                        const std::string idir = folder_KEcuts + ic.first + (inCuts0 ? "/" : "/1_additionalinfo/");
                        gSystem->mkdir(idir.c_str(), true);
                        canvas->SaveAs((idir + std::to_string(11 + w) + "_reconBlockDiff_" + rps[w].name + ".pdf").c_str());
                    }
                }
                lap(ic.first + ": drawing the plots");
            }
        }

    }
    BoxAndCropPlots();
}
