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

// Plot3D has four independent parts that each write their own plots/ folders: 1 position, energy, 2_KECUT (KECUT80 to KECUT95) and 1_hitZDC; 2 the cut ratios and rank csv files of the 6_cuts folders;
// 3 their 1_additionalinfo plots; 4 their phi fits. P3_PART=n runs part n only, so the four can run as four processes side by side (Plot.pbs, plot.sh). Without it all four run.
// Part 2 is one process per cut with P3_CUT=i (0: the 1_hitZDC baseline, 1 to 8: sampleCuts in order), each writing its rank rows to cuts/rank_rows_<i>.tsv; part 5 reads those
// files and writes the rank_*.csv files. Parts 3 and 4 also take P3_CUT=i to do one cut folder only (0: 1_hitZDC and 2_KECUT/KECUT80 to KECUT95). Without P3_PART all of it runs in one process.
static bool RunPart(int n){const char* e = std::getenv("P3_PART"); return !e || !*e || std::atoi(e) == n;}
static bool RunCut(int i){const char* e = std::getenv("P3_CUT"); return !e || !*e || std::atoi(e) == i;}
// P3_SUB (0 to 3): in part 1 the four groups 1_position/1_true, 1_position/2_block, 2_energy and the 80% KE cut plots (1_hitZDC, 2_KECUT, 0_cuts); in parts 3 and 4 of the four folders P3_CUT=0 stands for in parts 3 and 4 (1_hitZDC, 2_KECUT, 0_cuts/1_noCuts, 0_cuts/2_theta0_1) only this one, so they can run as four processes. Without it all four.
// The six counts of AnglesRaw every part 2, 4 and 5 process needs, in one event loop: CE protons / pdg 2212 / any particle with the theta baseline, then CE protons / pdg 2212 / any particle without it.
// Kept in cuts/plot3d_totals.csv (first line the theta baseline, second the counts) and reused while the file is newer than cut_angles_all.root and the baseline is the same; else counted again (written to a
// temporary file and renamed, so processes that start together do not clash).
// the six ratio plots of a list of samples (purity, contamination, background efficiency, signal significance, background rejection efficiency, signal efficiency), one point per sample, in folder
// ce / part: CE recon protons and all recon charged particles of each sample; ceBase / nonBase: the CE recon protons and the non-CE recon particles of the first sample (no cut), the base of the efficiencies
static void DrawCutRatios(TCanvas* canvas, const std::string& folder, const std::vector<std::string>& labels, const std::vector<double>& ce, const std::vector<double>& part, double ceBase, double nonBase, const std::string& what, const std::string& note){
    std::string fileTag = what;
    for (auto& ch : fileTag) if (ch == ' ') ch = '_';
    struct Ratio {const char* file, *label; int color; std::function<double(int)> f;};
    const int n = int(labels.size());
    const Ratio ratios[6] = {
        {"purity", "N_{CE,recon} / N_{recon,particle}", kGreen + 2, [&](int i){return part[i] > 0 ? ce[i] / part[i] : NAN;}},
        {"contamination", "N_{nonCE,recon,particle} / N_{recon,particle}", kGreen + 2, [&](int i){return part[i] > 0 ? (part[i] - ce[i]) / part[i] : NAN;}},
        {"background_efficiency", "N_{nonCE,recon,particle} / N_{nonCE,recon,particle,no cut}", kGreen + 2, [&](int i){return nonBase > 0 ? (part[i] - ce[i]) / nonBase : NAN;}},
        {"signal_significance", "N_{CE,recon} / sqrt(N_{nonCE,recon,particle})", kGreen + 2, [&](int i){return part[i] - ce[i] > 0 ? ce[i] / std::sqrt(part[i] - ce[i]) : NAN;}},
        {"background_rejection_efficiency", "1 - N_{nonCE,recon,particle} / N_{nonCE,recon,particle,no cut}", kGreen + 2, [&](int i){return nonBase > 0 ? 1 - (part[i] - ce[i]) / nonBase : NAN;}},
        {"signal_efficiency", "N_{CE,recon} / N_{CE,recon,no cut}", kRed, [&](int i){return ceBase > 0 ? ce[i] / ceBase : NAN;}}};
    gSystem->mkdir(folder.c_str(), true);
    for (int k = 0; k < 6; k++){
        TGraph g;
        for (int i = 0; i < n; i++){const double v = ratios[k].f(i); if (std::isfinite(v)) g.SetPoint(g.GetN(), i + 1, v);}
        if (g.GetN() < 2) continue;
        double lo = 1e300, hi = -1e300;
        for (int i = 0; i < g.GetN(); i++){lo = std::min(lo, g.GetY()[i]); hi = std::max(hi, g.GetY()[i]);}
        const double pad = 0.1 * std::max(hi - lo, 1e-9);
        canvas->Clear(); canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
        TH1F frame("hCutRatioFrame", (std::string(ratios[k].label) + " vs " + what + ", " + note + ";;" + ratios[k].label).c_str(), n, 0.5, n + 0.5);
        for (int i = 0; i < n; i++) frame.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
        frame.SetStats(0); frame.SetMinimum(lo - pad); frame.SetMaximum(hi + pad);
        frame.SetLineColor(kWhite); frame.Draw("hist"); // empty bins in white: the frame carries the axis labels and the title (the "axis" option draws no title)
        g.SetMarkerStyle(20); g.SetMarkerSize(0.6); g.SetLineWidth(1); g.SetLineColor(ratios[k].color); g.SetMarkerColor(ratios[k].color);
        g.Draw("PL same");
        canvas->SaveAs((folder + std::to_string(k + 1) + "_" + ratios[k].file + "_vs_" + fileTag + ".pdf").c_str());
    }
    canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
}
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
// P3_SEC: splits the slow noCuts folder (0_cuts/1_noCuts, every fitted proton) into processes (Plot.pbs); empty: everything.
//   part 3: 0 to 3 the histogram event loop and variable pdfs of the variables with index % 4 equal to it, 4 the block theta and phi resolution (9, 10), 5 recon minus recon block (11, 12)
//   part 4: the recon block phi fits (2_phi_cosfit), one section only, so 4:0:2 has no P3_SEC split
static bool RunSec(int i){const char* e = std::getenv("P3_SEC"); return !e || !*e || std::atoi(e) == i;}
// the folder of the ZDC KE cut at pct % of the beam KE, under plots/6_cuts: 2_KECUT/KECUT80, KECUT85, KECUT90, KECUT95
static std::string KeCutDir(int pct){return "2_KECUT/KECUT" + std::to_string(pct);}
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
    // eventsAllPart: the same for every recorded particle (charged, any pdg) at the ZDC with a reconstructed track (cut_hitZDC_AllPart.root), the "all recon charged particles" plots
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
                            const std::string& path, const char* opt = "", const std::function<void(TH1*)>& tweak = nullptr, int fitMode = kFitNone){ // fitMode: the unbinned gaussian fit of the plotted values (FitCore.h), before the tweak
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
            if (fitMode != kFitNone){ // the same values and cut as the histogram, read once more for the likelihood (no bins)
                const std::vector<double> xs = Columns(t, SplitDraw(P3(expr)).exprs, P3(cut))[0];
                if (fitMode == kFitCore) FitCore(h, xs); else FitBlocks(h, xs);
            }
            if (allPart1D) h->SetLineColor(kGreen + 2); // all recon charged particles are green, as in the x overlays
            if (tweak) tweak(h);
            if (ceOnly1D || allPart1D) for (TObject* o : *h->GetListOfFunctions()) if (TF1* fit = dynamic_cast<TF1*>(o)) fit->SetLineColor(kBlack); // a red histogram would hide a red fit
            const std::string hName = h->GetName(); // the stats box shows the name: a readable one until the pdf is saved
            h->SetName(StatsTitle(h, title.substr(0, std::min(title.find(','), axes)) + (t == events2 ? ", CE protons" : t == eventsAllEv ? ", all recon protons" : t == eventsAllPart ? ", all recon charged particles" : "")).c_str()); // quantity: the title up to its first comma
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

        // call in drawPlot's tweak: h (all recon protons, blue) with the CE protons (events2, red) and all recon charged particles (eventsAllPart, green) of the same quantity drawn over it, and a legend below the stats box.
        // logY: log y from 0.5, else a linear y 20% above the highest bin.
        auto overlayCEPart = [&](TH1* h, const std::string& exprCE, const std::string& cutCE, const std::string& exprPart, const std::string& cutPart, const char* ceLabel, bool logY){
            TH1* hCE = HistMT(events2, P3(exprCE), P3(cutCE)); hCE->Draw("same");
            TH1* hPart = HistMT(eventsAllPart, P3(exprPart), P3(cutPart)); hPart->Draw("same");
            h->SetLineColor(kBlue); hCE->SetLineColor(kRed); hPart->SetLineColor(kGreen + 2); hCE->SetStats(0); hPart->SetStats(0);
            const double top = std::max({h->GetMaximum(), hCE->GetMaximum(), hPart->GetMaximum()});
            if (logY){h->SetMinimum(0.5); h->SetMaximum(2.0 * top); canvas->SetLogy(1);} else h->SetMaximum(1.2 * top);
            TLegend* leg = new TLegend(0.745, 0.30, 0.99, 0.42); leg->SetTextSize(0.012); // below the stats box
            leg->AddEntry(hPart, "all recon charged particles that hit ZDC", "l");
            leg->AddEntry(h, "all recon protons that hit ZDC", "l");
            leg->AddEntry(hCE, ceLabel, "l");
            leg->Draw();
        };

        const std::string zdcProtons ="boundaryType == 2 && pdg == 2212";
        // position plots: only protons with a reconstructed track
        const std::string zdcReco = zdcProtons + " && isReco(fEvent, TrackID)", zdcRecoAll = zdcReco;
        // the same for every recorded particle (any pdg, eventsAllPart): ZDC hit and a recon track
        const std::string zdcRecoPartAll = "boundaryType == 2 && isReco(fEvent, TrackID)";
        // baseline of every plot and count in 6_cuts: recon theta and recon block theta from 0 to 1 degree (NaN fails), on top of the ZDC hit and the cut
        const std::string thetaBase = "rawVal(fEvent, TrackID, 9) >= 0 && rawVal(fEvent, TrackID, 9) <= 1 && rawVal(fEvent, TrackID, 11) >= 0 && rawVal(fEvent, TrackID, 11) <= 1", zdcBase = zdcRecoAll + " && " + thetaBase, zdcBasePart = zdcRecoPartAll + " && " + thetaBase; // zdcBasePart: every recorded particle (eventsAllPart), not only the protons
        countZdcAll = [&]{
            if (nCERecoAll >= 0) return;
            {const auto c = CountsMT(events2, {P3(zdcRecoAll), P3(zdcProtons)}); nCERecoAll = c[0]; nCEZdcAll = c[1];} // one event loop for both; nRecoAll = CountMT(eventsAllEv, P3(zdcRecoAll)); nRecoAllPart = CountMT(eventsAllPart, P3(zdcRecoPartAll));
        };
        const std::string kNoCutNote = "no theta cut"; // 0_cuts: only the ZDC hit and a recon track
        const std::string kThetaNote = "recon theta and recon block theta 0-1 deg"; // the thetaBase cut, named in the titles of the 6_cuts plots

        // Position plots have no KE cut: 1_position/1_true/ and 1_position/2_block/
        std::vector<std::pair<std::string, std::string>> posKes = {{"", ""}}; // unused folder, extra cut
        std::vector<std::string> posKeLabels = {""};
        if (RunPart(1)) for (size_t pk = 0; pk < posKes.size(); pk++){
            const std::string folder_pos = folder + "/1_position/1_true/", folder_posD = folder + "/1_position/2_block/";
            const std::string zdcReco = zdcRecoAll + posKes[pk].second, zdcRecoD = zdcReco; // the block position plots have the same cut as the true position ones
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
            drawPlot(events2, "fX>>hSliceX200(100,-200,200)", xySlice, "X Position at ZDC, CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_pos + "5_x_fit_200_CE.pdf", "", nullptr, kFitCore);
            drawPlot(events2, "fX>>hSliceX(100,-50,50)", xySlice, "X Position at ZDC, CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_pos + "7_x_fit_50_CE.pdf", "", nullptr, kFitCore);

            // the x slices as one overlay each, no fit: CE recon protons (red), all recon protons (blue), all recon charged particles (green)
            struct XOv {std::string bins, file, tag;};
            for (const auto& o : std::vector<XOv>{{"(100,-200,200)", "9_x_overlay_200.pdf", "200"}, {"(100,-50,50)", "10_x_overlay_50.pdf", "50"}}){
                drawPlot(eventsAllEv, "fX>>hXOv" + o.tag + o.bins, xySlice, "X Position at ZDC, all recon charged particles, all recon protons and CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_pos + o.file, "", [&](TH1* h){
                    overlayCEPart(h, "fX>>hXOvCE" + o.tag + o.bins, xySlice, "fX>>hXOvPart" + o.tag + o.bins, xySliceP, "CE recon protons that hit ZDC", false);
                });
            }


            }
            if (RunSub(1)){ // P3_SUB 1: 2_block
            titleSuffix = ", has recon track, no KE cut" + posKeLabels[pk]; // the cut of the 2_block plots: hit ZDC and a recon track
            // ZDC position, block (fDx, fDy): the ZDC is nCols_detector x nRows_detector blocks over 400 mm. Every bin holds the same whole number of blocks (BlockBins, FitCore.h):
            // 100 bins over 100 mm on a pitch of 400 / 801 mm put an extra block in the bin that starts at x = 0 (a cross at x = 0 and y = 0), and in one bin of ten along the x slices (spikes).
            const int nColsD = (int)ReadSettings().at("nCols_detector"), nRowsD = (int)ReadSettings().at("nRows_detector");
            const std::string binsSlice200 = BlockBins(nColsD, 400.0 / nColsD, 200), binsSlice50 = BlockBins(nColsD, 400.0 / nColsD, 53.5); // x slices: one bin per ZDC cell, the +-50 mm ones with whole cells past 50 (clipped by the frame, clipX)
            // the 4 xy plots: one marker per occupied ZDC cell (nCols_detector x nRows_detector over 400 mm), at the exact cell position (drawXYDots)
            const std::string binsCell = BlockBins(nColsD, 400.0 / nColsD, 200), binsCellY = BlockBins(nRowsD, 400.0 / nRowsD, 200), binsCellZoom = BlockBins(nColsD, 400.0 / nColsD, 53.5), binsCellZoomY = BlockBins(nRowsD, 400.0 / nRowsD, 53.5); // whole cells past +-50: the frame is +-50, the edge cells are clipped
            const std::string cellNote = Form(", 1 square per ZDC cell (%d x %d over 400 mm)", nColsD, nRowsD), zoomNote = cellNote + ", edge cells clipped at 50 mm";
            drawXYDots(eventsAllEv, "fDy:fDx>>h2(" + binsCell + "," + binsCellY + ")", zdcRecoD, "Position at ZDC BLOCK, all recon protons that hit ZDC" + cellNote + ";X Position (mm);Y Position (mm)", folder_posD + "2_xy_recon.pdf");
            drawXYDots(eventsAllEv, "fDy:fDx>>h2(" + binsCellZoom + "," + binsCellZoomY + ")", zdcRecoD, "Position at ZDC BLOCK, all recon protons that hit ZDC" + zoomNote + ";X Position (mm);Y Position (mm)", folder_posD + "4_xy_zoom_recon.pdf", 50);
            drawXYDots(events2, "fDy:fDx>>h2(" + binsCell + "," + binsCellY + ")", zdcRecoD, "Position at ZDC BLOCK, CE recon protons that hit ZDC" + cellNote + ";X Position (mm);Y Position (mm)", folder_posD + "1_xy_CE.pdf");
            drawXYDots(events2, "fDy:fDx>>h2(" + binsCellZoom + "," + binsCellZoomY + ")", zdcRecoD, "Position at ZDC BLOCK, CE recon protons that hit ZDC" + zoomNote + ";X Position (mm);Y Position (mm)", folder_posD + "3_xy_zoom_CE.pdf", 50);

            // the six x fits: ZDC hit and a recon track only
            const std::string xySliceB = zdcReco + " && fabs(fDy) < 2", xySliceBP = zdcRecoP + " && fabs(fDy) < 2";
            drawPlot(events2, "fDx>>hSliceDX200(" + binsSlice200 + ")", xySliceB, "X Position at ZDC BLOCK, CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_posD + "5_x_fit_200_CE.pdf", "", nullptr, kFitBlocks); //whole blocks per bin
            drawPlot(events2, "fDx>>hSliceDX(" + binsSlice50 + ")", xySliceB, "X Position at ZDC BLOCK, CE recon protons that hit ZDC, |Y| < 2 mm, edge cells clipped at 50 mm;X Position (mm);Count", folder_posD + "7_x_fit_50_CE.pdf", "", [&](TH1* h){clipX(h, 50);}, kFitBlocks);
            // the x slices as one overlay each, no fit: CE recon protons (red), all recon protons (blue), all recon charged particles (green)
            for (const auto& o : std::vector<std::array<std::string, 4>>{{"200", binsSlice200, "6_x_overlay_200.pdf", ""}, {"50", binsSlice50, "8_x_overlay_50.pdf", ", edge cells clipped at 50 mm"}}){
                const std::string tag = o[0], bins = "(" + o[1] + ")";
                drawPlot(eventsAllEv, "fDx>>hXOvD" + tag + bins, xySliceB, "X Position at ZDC BLOCK, all recon charged particles, all recon protons and CE recon protons that hit ZDC, |Y| < 2 mm" + o[3] + ";X Position (mm);Count", folder_posD + o[2], "", [&](TH1* h){
                    overlayCEPart(h, "fDx>>hXOvDCE" + tag + bins, xySliceB, "fDx>>hXOvDPart" + tag + bins, xySliceBP, "CE recon protons that hit ZDC", false);
                    if (tag == "50") h->GetXaxis()->SetRangeUser(-50, 50);
                });
            }

            // recon block phi (rawVal 15, 0 to 360 deg) of the same |Y| < 2 mm events as the x fits above, no fit
            const std::string phiOk = " && rawOk(fEvent, TrackID, 15)";
            // one version per x window, like the x fits: |X| < 200 mm (files 9, 10) and |X| < 50 mm (files 11, 12)
            struct PhiWin {int w; std::string nCE, nRecon;};
            for (const auto& pw : std::vector<PhiWin>{{200, "9", "10"}, {50, "11", "12"}}){
                const std::string tag = std::to_string(pw.w), xWin = Form(" && fabs(fDx) < %d", pw.w), winNote = Form(", |X| < %d mm", pw.w);
                drawPlot(events2, "rawVal(fEvent, TrackID, 15)>>hPhiDY2CE" + tag + "(100,0,360)", xySliceB + phiOk + xWin, "Recon block phi, CE recon protons that hit ZDC, |Y| < 2 mm" + winNote + ";Recon block phi (degrees);Count", folder_posD + pw.nCE + "_phiBlock_" + tag + "_CE.pdf");
                // CE recon protons (red), all recon protons (blue) and all recon charged particles (green) on one plot
                drawPlot(eventsAllEv, "rawVal(fEvent, TrackID, 15)>>hPhiDY2AllEv" + tag + "(100,0,360)", xySliceB + phiOk + xWin, "Recon block phi, all recon charged particles, all recon protons and CE recon protons that hit ZDC, |Y| < 2 mm" + winNote + ";Recon block phi (degrees);Count", folder_posD + pw.nRecon + "_phiBlock_overlay_" + tag + ".pdf", "", [&](TH1* h){
                    overlayCEPart(h, "rawVal(fEvent, TrackID, 15)>>hPhiDY2OvCE" + tag + "(100,0,360)", xySliceB + phiOk + xWin, "rawVal(fEvent, TrackID, 15)>>hPhiDY2OvPart" + tag + "(100,0,360)", xySliceBP + phiOk + xWin, "CE recon protons that hit ZDC", false);
                });
            }
            }
            titleSuffix.clear();
        }

        // KE x axes: 0 to KeAxisMax (FitCore.h, 30k MeV at 25 GeV/c) for the 3_KEcuts and 2_energy plots; the tag names the range in k MeV
        struct KERange {double lo, hi; std::string tag;};
        const double keAxisMax = KeAxisMax();
        auto keTag = [](double lo, double hi){return std::string(Form("_%gto%gk", lo / 1000, hi / 1000));};
        const std::vector<KERange> keRanges = {{0, keAxisMax, keTag(0, keAxisMax)}}; // 3_KEcuts: full range only
        const std::vector<KERange> energyRanges = {{0, keAxisMax, keTag(0, keAxisMax)}}; // 2_energy: 0 to 1.2 x the beam KE (KeAxisMax), full range only
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
            // one log y plot per quantity: CE protons (red), all recon protons (blue) and all recon charged particles (green): 1_true_ke, 2_smeared_ke
            for (const auto& q : keQtys){
                for (const auto& r : energyRanges){
                    drawPlot(eventsAllEv, q.branch + ">>hKE" + r.tag.substr(1) + binSpec(r), zdcRecoE, q.title + ", all recon charged particles, all recon protons and CE protons that hit ZDC (reconstructed or not);" + q.axis + ";Count",
                             folder_smear + std::to_string(&q - &keQtys[0] + 1) + "_" + q.file + r.tag + ".pdf", "", [&](TH1* h){
                        overlayCEPart(h, q.branch + ">>hKECE" + r.tag.substr(1) + binSpec(r), zdcProtonsE, q.branch + ">>hKEPart" + r.tag.substr(1) + binSpec(r), zdcRecoPartAll, "CE protons that hit ZDC (reco or not)", true); // every charged particle: contains the protons
                    });
                    canvas->SetLogy(0);
                }
            }

            // ZDC energy smearing (fDEnergy = Gaussian smear of fKEnergy, /zdc/resStochastic and /zdc/resConstant in the sim)
            if (events->GetBranch("fDEnergy")){ // files from before fDEnergy was added don't have it
                gSystem->mkdir(folder_smear.c_str(), true);

                // fKEnergy = true KE, fDEnergy = smeared KE, dE = fDEnergy - fKEnergy
                // the smear sigma is resStochastic * sqrt(KE in GeV) GeV, so the dE axis is +-10000 MeV at the 25 GeV/c beam (24078 MeV KE) and scales with sqrt(KE)
                const double smearLim = 10000.0 * std::sqrt(kBeamKE / 24078.0);

                drawPlot(events2, Form("fDEnergy-fKEnergy>>hSmearDiffCE(100,-%g,%g)", smearLim, smearLim), zdcProtonsE, "Energy Difference at ZDC, CE protons that hit ZDC (reconstructed or not);fDEnergy - fKEnergy (MeV);Count", folder_smear + "5_diff_CE.pdf", "", nullptr, kFitCore);

                // percent error = 100 * (fDEnergy - fKEnergy) / fKEnergy
                const std::string pct = "100*((fDEnergy-fKEnergy)/fKEnergy)"; // bracketed: the C++ of HistMT would divide 100*(a-b) by c, TTree::Draw 100*(a-b/c)
                drawPlot(events2, pct + ">>hSmearPctCE(100,-100,100)", zdcProtonsE, "Percent Error at ZDC, CE protons that hit ZDC (reconstructed or not);Percent Error (%);Count", folder_smear + "7_pct_CE.pdf", "", nullptr, kFitCore);
                

                // the smearing quantities as one log y overlay each: CE protons (red), all recon protons (blue), all recon charged particles (green)
                struct SmearOv {std::string expr, bins, file, axis, title, tag;};
                const std::vector<SmearOv> smearOvs = {{"fDEnergy-fKEnergy", Form("(100,-%g,%g)", smearLim, smearLim), "9_diff_overlay.pdf", "fDEnergy - fKEnergy (MeV)", "Energy Difference at ZDC", "Diff"},
                                                       {pct, "(100,-100,100)", "10_pct_overlay.pdf", "Percent Error (%)", "Percent Error at ZDC", "Pct"}};
                for (const auto& o : smearOvs){
                    drawPlot(eventsAllEv, o.expr + ">>hOv" + o.tag + o.bins, zdcRecoE, o.title + ", all recon charged particles, all recon protons and CE protons that hit ZDC (reconstructed or not);" + o.axis + ";Count",
                             folder_smear + o.file, "", [&](TH1* h){
                        overlayCEPart(h, o.expr + ">>hOvCE" + o.tag + o.bins, zdcProtonsE, o.expr + ">>hOvPart" + o.tag + o.bins, zdcRecoPartAll, "CE protons that hit ZDC (reco or not)", true);
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
        // the counts of the four KE cuts of the loop below (index (pct - 80) / 5), and the totals they are measured against: one event loop per tree for all of them, not one per cut
        std::array<long long, 4> keRealOut{}, keFakeKept{}, kePartKept{};
        if (RunPart(1) && RunSub(3) && !cutEnergies.empty()){
            std::vector<std::string> cReal = {P3(zdcRecoAll)}, cFake = {P3(zdcRecoAll)}, cPart = {P3(zdcRecoPartAll)};
            for (int pct : {80, 85, 90, 95}){
                const std::string cutStr = std::to_string(initialKE * pct / 100.0), br = cutEnergies[0].branch;
                cReal.push_back(P3(zdcRecoAll + " && " + br + " < " + cutStr));
                cFake.push_back(P3(zdcRecoAll + " && " + br + " > " + cutStr));
                cPart.push_back(P3(zdcRecoPartAll + " && " + br + " > " + cutStr));
            }
            const auto nReal = CountsMT(events2, cReal), nFake = CountsMT(eventsAllEv, cFake), nPart = CountsMT(eventsAllPart, cPart);
            for (int i = 0; i < 4; i++){keRealOut[i] = nReal[i + 1]; keFakeKept[i] = nFake[i + 1]; kePartKept[i] = nPart[i + 1];}
        }
        if (RunPart(1) && RunSub(3)) for (int pct : {80, 85, 90, 95}){ // 2_KECUT/KECUT<pct>
            const int pi = (pct - 80) / 5;
            const std::string zdcReco = zdcRecoAll;
            std::string cut = std::to_string(initialKE * pct / 100.0);
            std::string pctStr = std::to_string(pct) + "% of initial KE";
            for (const auto& e : cutEnergies){
                std::string folder_cut = folder_KEcuts + KeCutDir(pct) + "/";
                gSystem->mkdir(folder_cut.c_str(), true);
                const std::string axis = ";" + e.label + " (MeV);Count";
                const double cutVal = initialKE * pct / 100.0;
                // counts over all KE (not only the plotted range): CE protons left out of the cut out of all CE protons, all recon protons kept in out of all recon protons
                const long long nRealOut = keRealOut[pi], nRealTotal = CountMT(events2, P3(zdcReco)); // the totals come from the memo of the CountsMT above
                const long long nFakeKept = keFakeKept[pi], nFakeTotal = CountMT(eventsAllEv, P3(zdcReco));
                const double realOutPct = 100.0 * nRealOut / std::max<long long>(1, nRealTotal);
                const double fakeKeptPct = 100.0 * nFakeKept / std::max<long long>(1, nFakeTotal);
                for (const auto& r : keRanges){
                // real (reconstructed CE protons) and all reconstructed protons of all events on one plot, with the cut and what it does to each
                drawPlot(events2, e.branch + ">>hCutRealOv" + binSpec(r), zdcReco, e.label + " of Protons, CE recon protons vs all recon protons, both hit ZDC, " + pctStr + " cut" + axis, folder_cut + "3_real_vs_recon_" + e.dir + r.tag + ".pdf", "", [&](TH1* h){
                    h->SetLineColor(kRed);
                    TH1* hf = HistMT(eventsAllEv, P3(e.branch + ">>hCutFakeOv" + binSpec(r)), P3(zdcReco)); hf->Draw("same");
                    hf->SetLineColor(kBlue);
                    TH1* hp = HistMT(eventsAllPart, P3(e.branch + ">>hCutPartOv" + binSpec(r)), P3(zdcRecoPartAll)); hp->Draw("same"); // every recon charged particle (any pdg) that hit the ZDC
                    hp->SetLineColor(kGreen + 2);
                    const long long nPartKept = kePartKept[pi], nPartTotal = CountMT(eventsAllPart, P3(zdcRecoPartAll));
                    h->SetMinimum(0.5); h->SetMaximum(2.0 * std::max(h->GetMaximum(), std::max(hf->GetMaximum(), hp->GetMaximum()))); canvas->SetLogy(1); // log y, full x range of the file name, the cut is the dashed line
                    TLine* cl = new TLine(cutVal, 0.5, cutVal, h->GetMaximum()); cl->SetLineStyle(2); cl->Draw();
                    TLegend* leg = new TLegend(0.745, 0.44, 0.99, 0.55); leg->SetTextSize(0.012); // below the stats box; both samples hit the ZDC (title)
                    leg->AddEntry(h, Form("CE recon protons, %lld/%lld = %.1f%% left out", nRealOut, nRealTotal, realOutPct), "l");
                    leg->AddEntry(hf, Form("all recon protons, %lld/%lld = %.1f%% kept in", nFakeKept, nFakeTotal, fakeKeptPct), "l");
                    leg->AddEntry(hp, Form("all recon charged particles, %lld/%lld = %.1f%% kept in", nPartKept, nPartTotal, 100.0 * nPartKept / std::max<long long>(1, nPartTotal)), "l");
                    leg->Draw();
                });
                canvas->SetLogy(0);
                }
            }
        }
        if (RunPart(1) && RunSub(3) && !cutEnergies.empty()){ // 2_KECUT: purity, contamination, background efficiency, signal significance, background rejection and signal efficiency against the ZDC KE cut (80, 85, 90, 95%)
            const std::string folder_ke = folder_KEcuts + "2_KECUT/";
            gSystem->mkdir(folder_ke.c_str(), true);
            // x axis: no cut at all (every fitted track of AnglesRaw, the base of the megaCut ratios), the hit ZDC sample (ZDC hit, recon track, recon block theta 0 to 1 deg), then that sample with the ZDC KE cut
            TTree* baseT = OpenCut("cuts/cut_angles_all.root", "AnglesRaw");
            gROOT->cd();
            const std::array<double, 6> tot = baseT ? RawTotals(baseT, thetaBase) : std::array<double, 6>{};
            const double ceBase = tot[3], nonBase = tot[5] - tot[3];
            const int pcts[4] = {80, 85, 90, 95};
            const char* xLabels[6] = {"no cut", "hit ZDC", "KE>80%", "KE>85%", "KE>90%", "KE>95%"};
            double ceIn[6], partIn[6];
            ceIn[0] = ceBase; partIn[0] = tot[5];
            std::vector<std::string> kes = {""}, cCE, cPart; // the hit ZDC sample, then each KE cut: one event loop per tree for all five
            for (int i = 0; i < 4; i++) kes.push_back(" && " + cutEnergies[0].branch + " > " + std::to_string(initialKE * pcts[i] / 100.0));
            for (const auto& ke : kes){cCE.push_back(P3(zdcBase + ke)); cPart.push_back(P3(zdcBasePart + ke));}
            const auto nCEs = CountsMT(events2, cCE), nParts = CountsMT(eventsAllPart, cPart);
            for (int i = 0; i < 5; i++){ceIn[i + 1] = double(nCEs[i]); partIn[i + 1] = double(nParts[i]);}
            DrawCutRatios(canvas, folder_ke, std::vector<std::string>(xLabels, xLabels + 6), std::vector<double>(ceIn, ceIn + 6), std::vector<double>(partIn, partIn + 6), ceBase, nonBase, "KE cut", "ZDC KE cut on the ZDC energy, recon block theta 0 to 1 deg for hit ZDC and the KE cuts, no cut = every fitted track");
        }

        // Cuts on the reconstructed protons, no KE cut, each in its own folder under 6_cuts: dir, extra Draw cut, label, and the info variable folder the cut is on
        // (those plots sit in the cut folder itself, the other variables in its 1_additionalinfo). The same KE plot as 1_hitZDC with the cut applied goes in 1_additionalinfo.
        const std::string perPlaneCut = " && rawVal(fEvent, TrackID, 19) >= 0 && rawVal(fEvent, TrackID, 19) <= 1"; // at most 1 charged particle in every tracker plane of the event (nMaxChargedPerPlaneEv)
        const std::vector<SampleCut> sampleCuts = [&]{
          std::vector<SampleCut> v = {
            {"3_ResCuts", "rawVal(fEvent, TrackID, 1) < " + std::to_string(kRes4Mm), "cut resDRms (block residual) < 4 mm", "3_ResCuts", "cut_ResCuts"},
            {"4_nParticles/1_particleCut", "rawVal(fEvent, TrackID, 6) <= 1", "cut at most 1 reconstructed particle in the event", "4_nParticles", "cut_nParticles1"},
            {"4_nParticles/2_particleCut", "rawVal(fEvent, TrackID, 6) <= 2", "cut at most 2 reconstructed particles in the event", "4_nParticles", "cut_nParticles2"},
            {"5_nHits", EvenHitsDraw("rawVal(fEvent, TrackID, 10)"), "cut " + EvenHitsText() + " tracker hits", "5_nHits", "cut_nHits"},
          };
          const std::string hitsList = EvenHitsText();
          // 6_megaCut (no ZDC requirement), 7_megaCut2 (ZDC KE > 80%) and 9_megaCut3 (megaCut2 and the even charged tracker hits in the event, all particles added up), each with at most 1 charged particle per tracker plane and
          // 1_1particlecut: at most 1 particle
          auto addMega = [&](int level){ // 0: megaCut, 1: megaCut2, 2: megaCut3
              const bool zdc = level >= 1;
              const std::string top = level == 0 ? "6_megaCut" : level == 1 ? "7_megaCut2" : "9_megaCut3";
              std::string expr = "rawVal(fEvent, TrackID, 6) <= 1 && " + EvenHitsDraw("rawVal(fEvent, TrackID, 10)") + perPlaneCut + " && rawVal(fEvent, TrackID, 9) >= 0 && rawVal(fEvent, TrackID, 9) <= 1 && rawVal(fEvent, TrackID, 11) >= 0 && rawVal(fEvent, TrackID, 11) <= 1";
              if (zdc) expr += " && fDEnergy > " + std::to_string(initialKE * 0.80);
              if (level == 2) expr += " && " + EvenHitsDraw("rawVal(fEvent, TrackID, 17)") + " && rawVal(fEvent, TrackID, 7) == 1";
              const std::string label = std::string(level == 0 ? "MEGA CUT: at most " : level == 1 ? "MEGA CUT 2: ZDC KE > 80%, at most " : "MEGA CUT 3: ZDC KE > 80%, at most ") + "1 reconstructed particle in event AND number of times the reconstructed particle hits tracker [" + hitsList + "] AND at most 1 charged particle per tracker plane, recon theta and recon block theta from 0 degrees to 1 degree" + (level == 2 ? " AND charged tracker plane hits in the event (all particles added up) [" + hitsList + "] AND exactly 1 charged particle at the ZDC in the event" : "");
              const char* file = level == 0 ? "cut_megacut_1" : level == 1 ? "cut_megacut2_1" : "cut_megacut3_1";
              v.push_back({level == 0 ? top : top + "/1_1particlecut", expr, label, top, file, zdc}); // 6_megaCut holds its files directly, no 1_1particlecut subfolder
          };
          addMega(0); addMega(1);
          v.push_back({"8_nZdcParticles/1_1zdcparticlecut", "rawVal(fEvent, TrackID, 7) == 1", "cut exactly 1 charged particle at the ZDC in the event", "8_nZdcParticles", "cut_nZdc1", false}); // P3_CUT 7
          addMega(2); // P3_CUT 8
          v.push_back({"10_nChargedTrackerHits/1_evenHits", EvenHitsDraw("rawVal(fEvent, TrackID, 17)"), "cut " + hitsList + " charged tracker plane hits in the event", "10_nChargedTrackerHits", "cut_nChHits", false}); // P3_CUT 9
          v.push_back({"10_nChargedTrackerHits/2_twoHits", "rawVal(fEvent, TrackID, 17) == 2", "cut exactly 2 charged tracker plane hits in the event", "10_nChargedTrackerHits", "cut_nChHits2", false}); // P3_CUT 10
          return v;
        }();
        struct RankRow {std::string cut, sample; double nCEin, nAllIn, nCEtotal, nCEreconBefore, fom[2], pol[2], polErr[2], nPartIn, nNonProtBefore, nNonPartBefore;}; // nAllIn: recon protons (pdg 2212) in the cut, nPartIn: every recon charged particle in the cut, nNon*Before: the non-CE protons / particles of the sample the cut is made on (the base of the row), so out of cut = before - in;  fom, pol (fitted P), polErr: [0] recon theta / phi, [1] recon block theta / phi
        double nCEreconTotal = 0, nCEreconAll = 0; // nCEreconAll: the same with no cut at all, the theta baseline included (RawTotals tot[3]); every CE proton with a fitted track (AnglesRaw, isFromCE), no ZDC requirement: CE_before of the QualityStudy csv files, the CE_recon_before_cut of the csv files
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
            nCEreconTotal = tot[0]; nCEreconAll = tot[3];
            nNonProtTh = tot[1] - tot[0]; nNonPartTh = tot[2] - tot[0];
            nNonProtAll = tot[4] - tot[3]; nNonPartAll = tot[5] - tot[3];
        }
        P2Log("totals done");
        const std::string rawCE = "isFromCE==1";
        std::vector<RankRow> rankRows; // every cut folder, with and without the 80% KE cut, for the rank_*.csv files below
        // FOM like 5_fomCalc/<theta set>/5_integral_ce_x_A2_dtheta.pdf (PlotCEFraction.cpp): sum over 0.05 deg theta bins of N_CE(bin) / N_neutrons * A_y^2(theta bin center) * 0.05 deg,
        // N_CE = CE recon protons that pass q (zdcRecoAll and the cut), t from theta and A_y from t with the model formula, theta from 0 to 15 deg
        const double nNeutrons = ReadSettings().at("TotalNeutrons"), apScaleFom = ReadSettings().at("APScale");
        auto fomFromCols = [&](const std::vector<std::vector<double>>& c, double* fom){
            // c: raw value and rawOk flag of both thetas (columns 0 to 3, rankStats below reads them with the phi columns in one pass); the rows without a finite value are skipped like the rawOk cut did
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
        auto phiFitCols = [&](TTree* t, const std::vector<std::vector<double>>& cols, int j, const std::string& title, const std::string& path, double* pol, double* polErr){ // cols: theta, phi (degrees) of the protons that pass the fit cut, in entry order (rankStats)
            *pol = 0; *polErr = 0;
            const Long64_t n = cols[0].size();
            AMean aMean;
            const double apScalePhi = ReadSettings().at("APScale"); // settings.csv is read from disk on every call, keep it out of the event loop
            for (Long64_t i = 0; i < n; i++)
                aMean.Add(AEffOfTheta(cols[0][i], apScalePhi));
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
                AddUnbinnedFit(h, "fPhiCosCut", uf, aMean.Mean(), aMean.Err(), !path.empty());
                h->SetMinimum(0); h->SetMaximum(1.1 * h->GetMaximum());
                *pol = uf.p; *polErr = uf.err;
                if (!path.empty()){
                    h->SetName(StatsTitle(h).c_str()); // stats box title; h is deleted below
                    canvas->Update();
                    if (auto* st = (TPaveStats*)h->FindObject("stats")){st->SetX1NDC(0.76); st->SetX2NDC(0.995); st->SetY1NDC(0.55); st->SetY2NDC(0.83);}
                    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
                    tx.DrawLatex(0.13, 0.90, "Unbinned fit  N(#phi) #propto 1 + P #upoint A_{eff} cos#phi (every proton, no bins)");
                    tx.DrawLatex(0.13, 0.865, Form("A_{eff} = %.4f #pm %.4f (fixed, model A at the %s theta),   P = %.2f #pm %.2f", aMean.Mean(), aMean.Err(), j ? "recon block" : "recon", *pol, *polErr));
                    canvas->Modified(); canvas->Update();
                    canvas->SaveAs(path.c_str());
                    WritePlotNumbers(h, path, fitFilesStarted); // stats.csv and fit_values.csv (N0, P, A) next to the pdf
                }
            }
            if (h) delete h;
            canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
        };
        // FOM and the fitted P with its error, for recon theta / phi (j 0) and recon block theta / phi (j 1), of the protons of t that pass q: ONE pass over the sample (columns 0-3 theta, thetaBlock and their rawOk
        // flags, 4-7 phi, phiBlock and theirs), then the rows of each fit are picked here in the same order the cut string "rawOk(phi) && rawOk(theta)" (and for j 1 the block phi window, no spike at 0) gave
        auto rankStats = [&](TTree* t, const std::string& q, double* fom, double* pol, double* polErr){
            std::vector<std::string> exprs;
            for (int k : {9, 11, 14, 15}){exprs.push_back(P3(Form("rawVal(fEvent, TrackID, %d)", k))); exprs.push_back(P3(Form("rawOk(fEvent, TrackID, %d)", k)));} // th, ok, thBlock, okBlock, ph, ok, phBlock, ok
            const auto c = Columns(t, exprs, P3(q));
            fomFromCols(c, fom);
            for (int j = 0; j < 2; j++){
                const size_t iTh = 2 * j, iPh = 4 + 2 * j; // value column of theta and phi of this fit, its ok flag is the next column
                std::vector<std::vector<double>> sel(2);
                for (size_t i = 0; i < c[iTh].size(); i++){
                    if (c[iTh + 1][i] == 0 || c[iPh + 1][i] == 0) continue;
                    if (j && (c[iPh][i] < kBlockPhiLo || c[iPh][i] > kBlockPhiHi)) continue;
                    sel[0].push_back(c[iTh][i]); sel[1].push_back(c[iPh][i]);
                }
                phiFitCols(t, sel, j, "", "", &pol[j], &polErr[j]);
            }
        };
        TTree* ceAllTree = OpenCut("cuts/cut_CEall.root", "ChargeExchange_values"); // every charge exchange the sim made (CutEvents.cpp), counted once
        if (!ceAllTree) return;
        const Long64_t nCEallZdcAll = ceAllTree->GetEntries();
        if (RunPart(2) || RunPart(5)){
        if (RunPart(2) && RunCut(0)) for (const auto& e : cutEnergies) for (const bool ke80 : {false, true}) for (const bool noTheta : {false, true}){ // the uncut samples (1_hitZDC with the theta window, 0_cuts without) are in every ranking as the baseline
            if (noTheta && ke80) continue; // 0_cuts has no KE cut
            // 0_cuts is no cut at all, not even the ZDC hit: every fitted proton of AnglesRaw, like its plots; 1_hitZDC is the ZDC protons with the theta window
            TTree* tB = noTheta ? rawT : events2;
            const std::string q = noTheta ? rawCE : zdcBase + (ke80 ? " && " + e.branch + " > " + std::to_string(initialKE * 0.80) : "");
            const std::string qAll = noTheta ? "PDG == 2212" : q;
            double fom[2], pol[2], polErr[2]; rankStats(tB, q, fom, pol, polErr);
            const std::string qPart = noTheta ? "1" : zdcBasePart + (ke80 ? " && " + e.branch + " > " + std::to_string(initialKE * 0.80) : "");
            rankRows.push_back({noTheta ? "0_cuts (no cut)" : "1_hitZDC (no cut)", ke80 ? "cut and ZDC KE > 80%" : "cut", double(CountMT(tB, P3(q))), double(CountMT(noTheta ? rawT : eventsAllEv, P3(qAll))), double(nCEallZdcAll),
                                nCEreconTotal, {fom[0], fom[1]}, {pol[0], pol[1]}, {polErr[0], polErr[1]}, double(CountMT(noTheta ? rawT : eventsAllPart, P3(qPart))), noTheta ? nNonProtAll : nNonProtTh, noTheta ? nNonPartAll : nNonPartTh});
        }
        for (const auto& sc : sampleCuts) for (const auto& e : cutEnergies) for (const bool ke80 : {false, true}){ // 0_KE_hitZDC (the cut and hitting the ZDC) and 0_KE_KECUT80 (also the 80% KE cut)
            if (!RunPart(2) || !RunCut(int(&sc - &sampleCuts[0]) + 1)) continue;
            P2Log("cut " + sc.dir + (ke80 ? " ke80" : " nokecut") + " start");
            const std::string folder_res = folder_KEcuts + sc.dir + "/1_additionalinfo/";
            const double cutVal = initialKE * 0.80; // the 80% KE cut of 2_KECUT/KECUT80
            const std::string qCut = zdcBase + " && " + sc.expr + (ke80 ? " && " + e.branch + " > " + std::to_string(cutVal) : "");
            const long long nCEq = ke80 ? CountMT(events2, P3(qCut)) : 0, nAllq = ke80 ? CountMT(eventsAllEv, P3(qCut)) : 0; // only the ZDC KE rows and plots use them
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
                const long long nPartc = ke80 ? CountMT(eventsAllPart, P3(zdcBasePart + " && " + sc.expr + " && " + e.branch + " > " + std::to_string(cutVal))) : CountMT(tS, P3("1")); // every recon charged particle of the cut
                const double nonProtOut = nNonProtTh - double(nAllc - nCEc), nonPartOut = nNonPartTh - double(nPartc - nCEc); // the non-CE protons / particles the cut removes
                P2Log("  counts done");
                double fom[2], pol[2], polErr[2]; rankStats(tS, qSCE, fom, pol, polErr);
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
            if ((ke80 && sc.zdcKE80) || (!ke80 && sc.onFolder == "10_nChargedTrackerHits")) for (const auto& r : keRanges){ TH1 *hfKE = nullptr, *hpKE = nullptr; long long nPartKE = 0; /* the all recon protons / particles histograms and the particle count are the same for the log and the linear pdf: made once */ for (const bool lin : {false, true}){ // 0_KE_KECUT80 plots only for the cuts that already hit the ZDC with a KE cut (7_megaCut2, 9_megaCut3), and 0_KE_hitZDC (hit ZDC, no KE cut) for 10_nChargedTrackerHits; the others keep their KE distribution. The ZDC KE plot only with the KE cut: the cut alone needs no ZDC hit. lin: the same plot with a linear y axis, saved as ..._linear.pdf
                drawPlot(events2, e.branch + ">>hResCE" + binSpec(r), qCut, e.label + " of Protons that hit ZDC, " + sc.label + (ke80 ? ", ZDC KE > 80%" : "") + ", " + kThetaNote + ", CE recon protons vs all recon protons;" + e.label + " (MeV);Count", folder_res + (ke80 ? "0_KE_KECUT80_" : "0_KE_hitZDC_") + e.dir + r.tag + (lin ? "_linear" : "") + ".pdf", "", [&](TH1* h){
                    h->SetLineColor(kRed);
                    if (!hfKE){
                        const std::string qCutPart = zdcBasePart + " && " + sc.expr + (ke80 ? " && " + e.branch + " > " + std::to_string(cutVal) : ""); // the same cut on every recon charged particle (any pdg)
                        hfKE = HistMT(eventsAllEv, P3(e.branch + ">>hResAll" + binSpec(r)), P3(qCut));
                        hpKE = HistMT(eventsAllPart, P3(e.branch + ">>hResPart" + binSpec(r)), P3(qCutPart));
                        nPartKE = CountMT(eventsAllPart, P3(qCutPart));
                        hfKE->SetLineColor(kBlue); hpKE->SetLineColor(kGreen + 2);
                    }
                    TH1 *hf = hfKE, *hp = hpKE;
                    hf->Draw("same"); hp->Draw("same");
                    if (lin){ // per bin counts and the CE share of the recon protons and of the recon particles in the cut, beside the linear pdf
                        FILE* bf = fopen((folder_res + (ke80 ? "0_KE_KECUT80_" : "0_KE_hitZDC_") + e.dir + r.tag + "_perbin.csv").c_str(), "w");
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
                    leg->AddEntry(hp, Form("all recon charged particles that hit ZDC%s, %lld", ke80 ? ", KE > 80%" : "", nPartKE), "l");
                    leg->Draw();
                });
                canvas->SetLogy(0);
            } }
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
                {"purity_proton", [](const RankRow& r){return r.nCEin / std::max(1.0, r.nAllIn);}, false},
                {"signal_efficiency", [](const RankRow& r){return r.nCEin / std::max(1.0, r.nCEreconBefore);}, false},
                {"contamination_proton", [](const RankRow& r){return (r.nAllIn - r.nCEin) / std::max(1.0, r.nAllIn);}, true},
                {"FOM_theta", [](const RankRow& r){return r.fom[0];}, false}, {"FOM_thetaBlock", [](const RankRow& r){return r.fom[1];}, false},
                {"P_phiRe", [](const RankRow& r){return r.pol[0];}, false}, {"P_phiReD", [](const RankRow& r){return r.pol[1];}, false},
                {"purity", [](const RankRow& r){return r.nCEin / std::max(1.0, r.nPartIn);}, false},
                {"contamination", [](const RankRow& r){return (r.nPartIn - r.nCEin) / std::max(1.0, r.nPartIn);}, true},
                {"background_efficiency_proton", [](const RankRow& r){return (r.nAllIn - r.nCEin) / std::max(1.0, r.nNonProtBefore);}, true},
                {"background_efficiency", [](const RankRow& r){return (r.nPartIn - r.nCEin) / std::max(1.0, r.nNonPartBefore);}, true},
                {"signal_significance_proton", [](const RankRow& r){return r.nCEin / std::sqrt(std::max(1.0, r.nAllIn - r.nCEin));}, false},
                {"signal_significance", [](const RankRow& r){return r.nCEin / std::sqrt(std::max(1.0, r.nPartIn - r.nCEin));}, false}};
            const char* files[] = {"rank_purity_proton.csv", "rank_signal_efficiency.csv", "rank_contamination_proton.csv",
                                   "rank_FOM_theta.csv", "rank_FOM_thetaBlock.csv", "rank_P_phiRe.csv", "rank_P_phiReD.csv",
                                   "rank_purity.csv", "rank_contamination.csv",
                                   "rank_background_efficiency_proton.csv", "rank_background_efficiency.csv", "rank_signal_significance_proton.csv", "rank_signal_significance.csv"};
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
            if (!rankRows.empty()){ // 4_nParticles: the six ratios of no cut, the 1 and the 2 particle cut, from the rank rows above (the counts of part 2, no new count); no cut = the 0_cuts row
                std::vector<std::string> labels;
                std::vector<double> ce, part;
                double ceBase = 0, nonBase = 0;
                for (const std::pair<const char*, const char*> w : {std::make_pair("0_cuts (no cut)", "no cut"), std::make_pair("4_nParticles/1_particleCut", "1 particle cut"), std::make_pair("4_nParticles/2_particleCut", "2 particle cut")})
                    for (const auto& r : rankRows) if (r.cut == w.first && r.sample == "cut"){
                        if (labels.empty()){ceBase = r.nCEin; nonBase = r.nPartIn - r.nCEin;}
                        labels.push_back(w.second); ce.push_back(r.nCEin); part.push_back(r.nPartIn);
                    }
                if (labels.size() == 3) DrawCutRatios(canvas, folder_KEcuts + "4_nParticles/", labels, ce, part, ceBase, nonBase, "particle cut", "no ZDC requirement, cut samples with recon theta and block theta 0 to 1 deg, no cut = every fitted track");
            }
        }

        } // part 2

        // phi cos fits in every cut folder, like 3_ReconCalcVals/2_phiRe_cosfit.pdf, CE recon protons and all recon protons of the folder on the same graph, at their recon phi and recon block phi:
        // N(phi) = N0 (1 + P <A> cos(phi)), N0 = entries / 50 bins fixed, <A> = mean model A_y at the recon (block) theta of those protons fixed, only P fitted
        if (RunPart(4)) {
            struct PhiSet {std::string dir, expr; TTree* t; std::string base, basePart, note;}; // base: the sample the cut applies to (hit ZDC and reconstructed, or any fitted track for the cuts)
            std::vector<PhiSet> phiSets = {{"1_hitZDC", "", events2, zdcBase, zdcBasePart, kThetaNote}, {KeCutDir(80), " && fDEnergy > " + std::to_string(initialKE * 0.80), events2, zdcBase, zdcBasePart, kThetaNote}, {"0_cuts/1_noCuts", "", rawT, rawCE, "", ""}, {"0_cuts/2_theta0_1", "", events2, zdcBase, zdcBasePart, kThetaNote}};
            for (const auto& sc : sampleCuts) phiSets.push_back({sc.dir, "", cutTree(sc), rawCE, "", ""});
            for (const int family : {2, 3}) for (const int maxP : {1}) for (const int pct : {85, 90, 95}){ // KE scan of 7_megaCut2 and 9_megaCut3, the 1 particle cut (CutPlots.cpp (samples)): P3_CUT 9 to 14, PlotKEOverlay.cpp overlays them on the 80% folder
                TTree* t = OpenCut(Form("cuts/cut_megacut%d_%d_ke%d.root", family, maxP, pct), "AnglesRaw");
                if (t && !t->GetAlias("raw_0")) Plot3DRawAliases(t);
                gROOT->cd();
                if (t) phiSets.push_back({std::string(family == 2 ? "7_megaCut2/" : "9_megaCut3/") + (maxP == 1 ? "1_1particlecut" : "2_2particlecut") + Form("_KE%d", pct), "", t, rawCE, "", ""});
            }
            const int nPlain = int(phiSets.size());
            for (const int pct : {85, 90, 95}) phiSets.push_back({KeCutDir(pct), " && fDEnergy > " + std::to_string(initialKE * pct / 100.0), events2, zdcBase, zdcBasePart, kThetaNote}); // run with KECUT80 (P3_CUT 0, P3_SUB 1)
            for (const auto& ps : phiSets){
                if (const int pi = int(&ps - &phiSets[0]), kx = pi >= nPlain; !RunCut(pi < 4 || kx ? 0 : pi - 3) || ((pi < 4 || kx) && !RunSub(kx ? 1 : pi))) continue; // P3_CUT: 0 the first four sets (P3_SUB picks one) and the KECUT85 to 95 ones (kx, with KECUT80, sub 1), i the sampleCuts[i - 1] set
                if (ps.dir.rfind("6_megaCut", 0) != 0 && ps.dir.rfind("7_megaCut2", 0) != 0 && ps.dir.rfind("9_megaCut3", 0) != 0) continue; // phi cos fits only in the mega cut folders
                gSystem->mkdir((folder_KEcuts + ps.dir + "/2_phi_cosfit").c_str(), true);
                const bool isZdc = ps.t == events2; // the ZDC folders hold the protons that hit the ZDC, the cut folders the fitted protons that pass the cut
                struct Samp {TTree* t; std::string q, name; int color;}; // CE recon protons (red), all recon protons (pdg 2212, blue) and all recon charged particles (any pdg, green) of the folder, drawn together
                const Samp samps[3] = {{ps.t, ps.base + ps.expr, "CE recon protons", kRed},
                                       {isZdc ? eventsAllEv : ps.t, (isZdc ? ps.base : std::string("PDG == 2212")) + ps.expr, "all recon protons", kBlue},
                                       {isZdc ? eventsAllPart : ps.t, (isZdc ? ps.basePart : std::string("1")) + ps.expr, "all recon charged particles", kGreen + 2}};
                const std::string cutTxt = [&]() -> std::string { // every cut the folder applies, for the title
                    if (ps.dir == "0_cuts/1_noCuts") return "no cuts, every fitted proton";
                    if (isZdc) return "hit ZDC, " + ps.note + (ps.expr.empty() ? std::string(", no KE cut") : ", ZDC KE > " + ps.dir.substr(ps.dir.rfind("KECUT") + 5) + "%");
                    for (const auto& c : sampleCuts) if (c.dir == ps.dir) return c.label + (c.zdcKE80 ? ", hit ZDC" : ", no ZDC requirement");
                    if (const size_t k = ps.dir.rfind("_KE"); k != std::string::npos) return "9_megaCut3 1 particle cut, hit ZDC, ZDC KE > " + ps.dir.substr(k + 3) + "%";
                    return ps.dir;
                }();
                const bool noCeil = ps.dir == "0_cuts/1_noCuts"; // no cuts: the start-theta scans count every proton from the start value on, with no upper theta
                const double thHi = noCeil ? INFINITY : 1.0;
                const std::string upTxt = noCeil ? " and up (no upper limit)" : " to 1 deg";
                const std::string where = isZdc ? " that hit ZDC" : " in the cut", dirP = folder_KEcuts + ps.dir + "/2_phi_cosfit/";
                const double apScalePhi = ReadSettings().at("APScale"), pTrue = ReadSettings().at("PPValue");
                for (int j = 1; j < 2; j++){ // j 1: phiReD at thetaReD, the recon block phi fits (the recon phiRe fits, j 0, are not made)
                    const int kT = j ? 11 : 9, kPh = j ? 15 : 14;
                    const std::string rb = j ? "recon block" : "recon", stem = j ? "phiReD" : "phiRe";
                    const std::string dirJ = dirP; // the plots sit in 2_phi_cosfit itself
                    std::vector<double> th[3], ph[3]; // theta and phi of every proton of the sample
                    StartScan scanned[3]; // per sample, the phi histograms and A sums for every start theta, from one pass over the protons (ScanStartTheta, Plot3DCols.h), and the mean A
                    for (int s = 0; s < 3; s++){
                        std::string qs = samps[s].q + Form(" && rawOk(fEvent, TrackID, %d) && rawOk(fEvent, TrackID, %d)", kT, kPh);
                        if (j) qs += std::string(Form(" && rawVal(fEvent, TrackID, %d) >= %g && rawVal(fEvent, TrackID, %d) <= %g", kPh, kBlockPhiLo, kPh, kBlockPhiHi)); // block phi: no spike at 0 (block theta 0 to 1 deg is in the baseline)
                        const std::string colT = Form("rawVal(fEvent, TrackID, %d)", kT), colPh = Form("rawVal(fEvent, TrackID, %d)", kPh);
                        // the rows only feed histograms and sums, so a multi-cluster tree (a cut file of AnglesRaw) is read on all threads in any order
                        auto c = MultiCluster(samps[s].t) ? ColumnsUnordered(samps[s].t, {P3(colT), P3(colPh)}, P3(qs)) : Columns(samps[s].t, {P3(colT), P3(colPh)}, P3(qs));
                        th[s] = std::move(c[0]); ph[s] = std::move(c[1]);
                        scanned[s] = ScanStartTheta(th[s], ph[s], apScalePhi, j ? kBlockPhiLo : 0, j ? kBlockPhiHi : 360, thHi);
                    }
                    // 1: phi distribution and cos fit of the three samples on one plot, N(phi) = N0 (1 + P <A> cos(phi)), N0 = entries / 50 bins fixed, <A> = mean model A of the sample fixed, only P fitted
                    const std::string allNames = "CE recon protons, all recon protons and all recon charged particles";
                    for (const bool ceOnly : {false, true}){ // the CE recon protons get a plot of their own (_CE) with their own y range; the other plot has the all recon protons and all recon charged particles
                        const int sLo = ceOnly ? 0 : 1, sHi = ceOnly ? 1 : 3;
                        TH1D* hP[3] = {nullptr, nullptr, nullptr};
                        double pol[3] = {0, 0, 0}, polErr[3] = {0, 0, 0}, hMax = 0;
                        AMean aMean[3];
                        for (int s = sLo; s < sHi; s++){
                            delete gDirectory->Get(Form("hPhiOv%d", s));
                            hP[s] = new TH1D(Form("hPhiOv%d", s), "", 50, 0, 360);
                            for (double x : ph[s]) hP[s]->Fill(x);
                            aMean[s] = scanned[s].aAll;
                            hP[s]->SetLineColor(samps[s].color); hP[s]->SetLineWidth(2);
                            if (hP[s]->GetEntries() >= 10){
                                const PFit uf = UnbinnedP(ph[s], aMean[s], j ? kBlockPhiLo : 0, j ? kBlockPhiHi : 360); // unbinned fit of every proton's phi (FitCore.h)
                                AddUnbinnedFit(hP[s], Form("fPhiOv%d", s), uf, aMean[s].Mean(), aMean[s].Err(), false); // the fit stays in the histogram for stats.csv and fit_values.csv, drawn below
                                pol[s] = uf.p; polErr[s] = uf.err;
                            }
                            hMax = std::max(hMax, hP[s]->GetMaximum());
                        }
                        canvas->Clear(); canvas->cd();
                        gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.36); gPad->SetBottomMargin(0.13); // the text lines sit above the frame
                        for (int s = sLo; s < sHi; s++){
                            hP[s]->SetStats(0); hP[s]->SetMinimum(0); hP[s]->SetMaximum(1.1 * hMax);
                            hP[s]->SetTitle(Form("%s%s, %s, %s phi distribution;Angle (degrees);Count", (ceOnly ? "CE recon protons" : "all recon protons and all recon charged particles"), where.c_str(), cutTxt.c_str(), rb.c_str()));
                            hP[s]->Draw(s > sLo ? "hist same" : "hist");
                        }
                        for (int s = sLo; s < sHi; s++) if (hP[s]->GetEntries() >= 10){ // the fit curve of each set in its color, a function of its own: a clone of the one in the histogram crashed when the canvas was cleared
                            TF1* fd = new TF1(Form("fPhiDraw%d", s), "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                            fd->SetParameters(hP[s]->GetEntries() / 50, pol[s], aMean[s].Mean()); fd->SetLineColor(samps[s].color); fd->SetLineStyle(2); fd->SetLineWidth(3); fd->Draw("same");
                        }
                        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.025); // small enough for the whole P = ... #pm ... of the longest line to stay on the page
                        tx.DrawLatex(0.18, 0.90, ("Unbinned fit  N(#phi) #propto 1 + P #upoint A_{eff} cos#phi, dashed = fit, A_{eff} fixed at the " + rb + " theta").c_str());
                        for (int s = sLo; s < sHi; s++){
                            tx.SetTextColor(samps[s].color);
                            tx.DrawLatex(0.18, 0.855 - 0.045 * (s - sLo), Form("%s: %.0f entries, A_{eff} = %.4f #pm %.4f,  P = %.3g #pm %.2g", samps[s].name.c_str(), hP[s]->GetEntries(), aMean[s].Mean(), aMean[s].Err(), pol[s], polErr[s]));
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
                    struct Page {TH1D* h; double start, P, eP, A, eA;}; // the phi distribution and fit behind the scan point every 0.1 deg, for 4_..._cosfit_by_start_theta
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
                            if (b % 10 == 0){TH1D* hc = (TH1D*)hS.Clone(Form("hPage%d_%d", s, b)); hc->SetDirectory(nullptr); pages[s].push_back({hc, start, uf.p, uf.err, aSum / hS.GetEntries(), AMean{hS.GetEntries(), aSum, scanned[s].aSum2[b]}.Err()});}
                            scan[s].x.push_back(start); scan[s].y.push_back(uf.p); scan[s].ey.push_back(uf.err);
                            scan[s].kept.push_back(100.0 * hS.GetEntries() / std::max(1.0, s ? double(n) : (noCeil ? nCEreconAll : nCEreconTotal)));
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
                                h->SetTitle(Form("%s%s, %s, %s phi distribution, %s theta %.2f%s;Angle (degrees);Count", allNames.c_str(), where.c_str(), cutTxt.c_str(), rb.c_str(), rb.c_str(), key / 100.0, upTxt.c_str()));
                                h->Draw(first ? "hist" : "hist same"); first = false;
                            }
                            std::vector<TF1*> fds;
                            for (int s = 0; s < 3; s++) if (pgOf[s]){
                                TF1* fd = new TF1(Form("fPageDraw%d_%zu", s, k), "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                                fd->SetParameters(pgOf[s]->h->GetEntries() / 50, pgOf[s]->P, pgOf[s]->A); fd->SetLineColor(samps[s].color); fd->SetLineStyle(2); fd->SetLineWidth(3); fd->Draw("same");
                                fds.push_back(fd);
                            }
                            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.033);
                            tx.DrawLatex(0.18, 0.90, Form("start theta %.2f deg: unbinned fit N(#phi) #propto 1 + P #upoint A_{eff} cos#phi, dashed = fit", key / 100.0));
                            for (int s = 0; s < 3; s++) if (pgOf[s]){
                                tx.SetTextColor(samps[s].color);
                                tx.DrawLatex(0.18, 0.855 - 0.045 * s, Form("%s: %.0f entries, A_{eff} = %.4f #pm %.4f,  P = %.3g #pm %.2g", samps[s].name.c_str(), pgOf[s]->h->GetEntries(), pgOf[s]->A, pgOf[s]->eA, pgOf[s]->P, pgOf[s]->eP));
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
                    for (const int mode : {0, 1, 2, 3, 4, 5}){ // 2, 3: the P scan of one sample per plot, mode / 2 the sample: the CE recon protons (_CE), all recon protons (_proton), all recon charged particles (_particle)
                        const bool pct = mode & 1;
                        const int s = mode / 2;
                        if (s == 1 || scan[s].x.size() <= 1) continue; // no all recon protons plots (s 1), only the CE recon protons and all recon charged particles
                        canvas->Clear(); canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
                        const double yRef = pct ? 0.0 : pTrue;
                        auto conv = [&](double P){return pct ? 100 * (P - pTrue) / pTrue : P;}; // P to the plotted value
                        double ymin = yRef, ymax = yRef; const double xMax = 1.0; // the scan always ends at 1 deg, the upper theta of the fits
                        for (size_t i = 0; i < scan[s].x.size(); i++){ // the error bars are cut at P = -1 and 1
                            ymin = std::min(ymin, conv(std::max(scan[s].y[i] - scan[s].ey[i], -1.0))); ymax = std::max(ymax, conv(std::min(scan[s].y[i] + scan[s].ey[i], 1.0)));
                        }
                        const double pad = 0.05 * (ymax - ymin);
                        TH2D axis("hPstartAxis", Form("Polarization fit P of %s%s, %s, %s phi, %s theta from the start value%s%s;Start theta (deg);%s", samps[s].name.c_str(), where.c_str(), ps.dir.c_str(), rb.c_str(), rb.c_str(), upTxt.c_str(), j ? ", phi 0.5 to 359.5 deg" : "",
                                  pct ? "100 (P_{fit} - P_{true}) / P_{true} (%)" : "P"), 100, 0, xMax, 100, ymin - pad, ymax + pad);
                        axis.SetStats(0);
                        axis.Draw("axis");
                        TLine ref(0, yRef, xMax, yRef); ref.SetLineStyle(3); ref.SetLineWidth(2); ref.SetLineColor(kGray + 2); ref.Draw(); // true P, from the run's PPValue
                        TLegend leg(0.705, 0.76, 0.99, 0.90); leg.SetTextSize(0.014);
                        TGraphAsymmErrors* gs[3] = {nullptr, nullptr, nullptr};
                        const std::string path = dirJ + std::to_string(2 + pct) + "_" + stem + (pct ? "_P_vs_start_theta_pct_error" : "_P_vs_start_theta") + (s == 0 ? "_CE" : s == 1 ? "_proton" : "_particle") + ".pdf";
                        { // the one sample s
                            gs[s] = new TGraphAsymmErrors();
                            for (size_t i = 0; i < scan[s].x.size(); i++){ // the error bar of P stops at -1 and 1, the physical range
                                const double y = scan[s].y[i], lo = std::max(y - scan[s].ey[i], -1.0), hi = std::min(y + scan[s].ey[i], 1.0);
                                gs[s]->SetPoint(i, scan[s].x[i], conv(y));
                                gs[s]->SetPointError(i, 0, 0, conv(y) - conv(lo), conv(hi) - conv(y));
                            }
                            gs[s]->SetLineColorAlpha(samps[s].color, 0.3); gs[s]->SetMarkerColor(samps[s].color); gs[s]->SetMarkerStyle(20); gs[s]->SetMarkerSize(1.0); // see-through error bars, as in PlotKEOverlay.cpp
                            if (gs[s]->GetN()) gs[s]->Draw("P same"); // with the end caps of the error bars
                            leg.AddEntry(gs[s], samps[s].name.c_str(), "lp");
                            for (size_t i = 0; i < scan[s].x.size(); i++) // the points drawn: start theta, value, error, % kept
                                WriteStatsRow(path, "plot,sample,start_theta_deg,value,error,pct_kept", Form("%s,%s,%.4g,%.6g,%.4g,%.4g", gSystem->BaseName(path.c_str()), samps[s].name.c_str(), scan[s].x[i], gs[s]->GetY()[i], pct ? 100 * scan[s].ey[i] / pTrue : scan[s].ey[i], scan[s].kept[i]), fitFilesStarted, "P_vs_start_theta.csv");
                        }
                        leg.Draw();
                        canvas->SaveAs(path.c_str());
                        for (auto* g : gs) delete g; // null for the sample not drawn
                    }
                    if (j == 1){ // 5 to 12: ratios of the counts with block theta from the start value to 1 deg ("in cut") against the start value, CE = the CE recon protons, prot = all recon protons, part = all recon charged particles of the folder
                        const double nonProtBefore = nNonProtAll, nonPartBefore = nNonPartAll, ceTotal = nCEreconAll; // "before cut" here is the run total with no cut at all: every fitted CE proton, non-CE proton and non-CE particle of AnglesRaw (RawTotals), no ZDC requirement and no theta baseline, not the folder's sample // before the start cut: the whole sample
                        struct Ratio {const char* file, *label; int color; std::function<double(double, double, double)> f;}; // f(CE, prot, part) in cut
                        const int fileNum[6] = {5, 6, 8, 10, 11, 12}; // the numbers the particle ratios kept when the all recon protons ones (7, 9) went
                        const Ratio ratios[6] = {
                            {"purity", "N_{CE,recon,in cut} / N_{recon,particle,in cut}", kGreen + 2, [&](double ce, double, double pa){return pa > 0 ? ce / pa : NAN;}},
                            {"contamination", "N_{nonCE,recon,particle,in cut} / N_{recon,particle,in cut}", kGreen + 2, [&](double ce, double, double pa){return pa > 0 ? (pa - ce) / pa : NAN;}},
                            {"background_efficiency", "N_{nonCE,recon,particle,in cut} / N_{nonCE,recon,particle,before cut}", kGreen + 2, [&](double ce, double, double pa){return nonPartBefore > 0 ? (pa - ce) / nonPartBefore : NAN;}},
                            {"signal_significance", "N_{CE,recon,in cut} / sqrt(N_{nonCE,recon,particle,in cut})", kGreen + 2, [&](double ce, double, double pa){return pa - ce > 0 ? ce / std::sqrt(pa - ce) : NAN;}},
                            {"background_rejection_efficiency", "1 - N_{nonCE,recon,particle,in cut} / N_{nonCE,recon,particle,before cut}", kGreen + 2, [&](double ce, double, double pa){return nonPartBefore > 0 ? 1 - (pa - ce) / nonPartBefore : NAN;}},
                            {"signal_efficiency", "N_{CE,recon,in cut} / N_{CE,recon,total}", kRed, [&](double ce, double, double){return ceTotal > 0 ? ce / ceTotal : NAN;}}};
                        std::vector<double> nIn[3]; // counts with theta from the start value to 1 deg, per start value (0.01 deg steps): theta sorted once, two binary searches per start
                        std::vector<double> sortedTh[3];
                        for (int s = 0; s < 3; s++){
                            sortedTh[s] = th[s];
                            std::sort(sortedTh[s].begin(), sortedTh[s].end());
                            const auto last = std::upper_bound(sortedTh[s].begin(), sortedTh[s].end(), thHi);
                            for (int b = 0; b < 100; b++) nIn[s].push_back(double(last - std::lower_bound(sortedTh[s].begin(), last, b * 0.01)));
                        }
                        for (int k = 0; k < 6; k++){
                            TGraph g;
                            for (int b = 0; b < 100; b++){const double v = ratios[k].f(nIn[0][b], nIn[1][b], nIn[2][b]); if (std::isfinite(v)) g.SetPoint(g.GetN(), b * 0.01, v);}
                            if (g.GetN() < 2) continue;
                            const std::string path = dirJ + std::to_string(fileNum[k]) + "_" + stem + "_" + ratios[k].file + "_vs_start_theta.pdf";
                            canvas->Clear(); canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
                            g.SetTitle((std::string(ratios[k].label) + " vs start theta, " + cutTxt + ", " + rb + " theta from the start value" + upTxt + where + ";Start theta (deg);" + ratios[k].label).c_str());
                            DrawPointsLine(g, ratios[k].color, 1.0);
                            TLegend leg(0.705, 0.82, 0.99, 0.90); leg.SetTextSize(0.014); leg.AddEntry(&g, ratios[k].file, "lp"); leg.Draw();
                            canvas->SaveAs(path.c_str());
                            for (int i = 0; i < g.GetN(); i++) WriteStatsRow(path, "plot,start_theta_deg,value", Form("%s,%.4g,%.6g", gSystem->BaseName(path.c_str()), g.GetX()[i], g.GetY()[i]), fitFilesStarted, "ratios_vs_start_theta.csv");
                        }
                        { // 13: the FOM of the CE recon protons with block theta from the start value to 1 deg, same sum as fomStats (A_y^2 at the 0.05 deg bin center * 0.05 deg / N_neutrons)
                            TGraph g;
                            const std::vector<double>& sorted = sortedTh[0]; // the terms summed in theta order, so a start value is a difference of two running sums
                            std::vector<double> running(sorted.size() + 1, 0.0);
                            for (size_t i = 0; i < sorted.size(); i++)
                                running[i + 1] = running[i] + std::pow(APFromT(kBeamMomentum, TFromAngle(kBeamMomentum, (std::floor(sorted[i] / 0.05) + 0.5) * 0.05 * TMath::Pi() / 180), apScaleFom), 2) * 0.05 / nNeutrons;
                            const auto last = std::upper_bound(sorted.begin(), sorted.end(), thHi) - sorted.begin();
                            for (int b = 0; b < 100; b++){
                                const auto first = std::lower_bound(sorted.begin(), sorted.begin() + last, b * 0.01) - sorted.begin();
                                g.SetPoint(g.GetN(), b * 0.01, running[last] - running[first]);
                            }
                            const std::string path = dirJ + "13_" + stem + "_FOM_vs_start_theta.pdf";
                            canvas->Clear(); canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
                            g.SetTitle(("FOM of CE recon protons vs start theta, " + cutTxt + ", " + rb + " theta from the start value" + upTxt + where + ";Start theta (deg);FOM = #Sigma A_{y}^{2} #Delta#theta / N_{neutrons}").c_str());
                            DrawPointsLine(g, kRed, 1.0);
                            canvas->SaveAs(path.c_str());
                            for (int i = 0; i < g.GetN(); i++) WriteStatsRow(path, "plot,start_theta_deg,value", Form("%s,%.4g,%.6g", gSystem->BaseName(path.c_str()), g.GetX()[i], g.GetY()[i]), fitFilesStarted, "ratios_vs_start_theta.csv");
                        }
                    }
                    canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
                }
            }
        }

        // recon protons that hit the ZDC, all KE (no cut): CE protons and every recon proton of all events
        if (RunPart(1) && RunSub(3)) for (const auto& e : cutEnergies) for (const auto& r : keRanges) for (const int var : {0, 1, 2}){ // 1_hitZDC, 0_cuts/1_noCuts (hit ZDC only), 0_cuts/2_theta0_1 (the same sample as 1_hitZDC)
            const bool noTheta = var == 1;
            const std::string hitDir = var == 0 ? "1_hitZDC/" : var == 1 ? "0_cuts/1_noCuts/" : "0_cuts/2_theta0_1/", hb = noTheta ? zdcRecoAll : zdcBase, hbp = noTheta ? zdcRecoPartAll : zdcBasePart, hnote = noTheta ? kNoCutNote : kThetaNote;
            std::string folder_hit = folder_KEcuts + hitDir;
            gSystem->mkdir(folder_hit.c_str(), true);
            drawPlot(events2, e.branch + ">>hHitCE" + binSpec(r), hb, e.label + " of Protons that hit the ZDC, CE recon protons that hit ZDC vs all events recon, " + hnote + ", no KE cut;" + e.label + " (MeV);Count", folder_hit + "1_real_vs_recon_" + e.dir + r.tag + ".pdf", "", [&](TH1* h){
                h->SetLineColor(kRed);
                TH1* hr = HistMT(eventsAllEv, P3(e.branch + ">>hHitRecon" + binSpec(r)), P3(hb)); hr->Draw("same");
                hr->SetLineColor(kBlue);
                TH1* hpr = HistMT(eventsAllPart, P3(e.branch + ">>hHitPart" + binSpec(r)), P3(hbp)); hpr->Draw("same"); // every recon charged particle (any pdg) that hit the ZDC
                hpr->SetLineColor(kGreen + 2);
                h->SetMinimum(0.5); h->SetMaximum(2.0 * std::max(h->GetMaximum(), std::max(hr->GetMaximum(), hpr->GetMaximum()))); canvas->SetLogy(1);
                TLegend* leg = new TLegend(0.745, 0.44, 0.99, 0.55); leg->SetTextSize(0.015);
                leg->AddEntry(h, Form("CE recon protons that hit ZDC, %lld", (long long)h->GetEntries()), "l");
                leg->AddEntry(hr, Form("all recon protons that hit ZDC, %lld", (long long)hr->GetEntries()), "l");
                leg->AddEntry(hpr, Form("all recon charged particles that hit ZDC, %lld", (long long)hpr->GetEntries()), "l");
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
            {"nHits", "5_nHits", "tracker hits", 10, -0.5, MaxTrackerHits() + 2.5, MaxTrackerHits() + 3}, // one bin per count from 0 to 2 * nBlocks + 2
            {"nZdcCharged", "8_nZdcParticles", "charged particles at the ZDC in the event", 7, -0.5, 9.5, 10},
            {"nChargedTrackerHitsEv", "10_nChargedTrackerHits", "charged tracker plane hits in the event", 17, -0.5, 5 * MaxTrackerHits() + 0.5, 5 * MaxTrackerHits() + 1}, // all charged particles of the event added up; 5 per tracker hit slot, 0 to 40 for 4 blocks
            {"thetaBlock", "7_theta", "recon block theta (degrees)", 11, 0, ThetaAxisMax(), 100}, // 0 to 2.3 deg at 25 GeV/c, scaled with 1 / p (FitCore.h)
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
            if (!(q.second > q.first)) q.second = q.first + 1.0; // no values, or all equal (resDRms is 0 for every 2-hit track of a 1-block run): a zero-width axis would not draw
            quantileCache[k] = q;
            const std::string tmp = std::string(qFile) + "." + std::to_string(getpid());
            if (FILE* qf = fopen(tmp.c_str(), "w")){ // written to a temporary name and renamed, two processes never read a half file
                for (const auto& kv : quantileCache) fprintf(qf, "%d\t%.17g\t%.17g\n", kv.first, kv.second.first, kv.second.second);
                fclose(qf);
                gSystem->Rename(tmp.c_str(), qFile);
            }
            return q;
        };
        // P3_PART=0 (Plot.pbs / plot.sh run it first, with all the threads) fills cuts/info_quantiles.tsv once: the part 3 processes start together, and each of them would
        // otherwise find no file and read all of AnglesRaw itself (425 s in each of about 13 processes on 100M events)
        if (const char* partEnv0 = std::getenv("P3_PART"); partEnv0 && *partEnv0 && std::atoi(partEnv0) == 0 && rawT && rawT->GetAlias("raw_0") && rawT->GetEntries() > 0)
            for (const auto& v : kInfoVars) if (v.lo >= v.hi) quantiles(v.k);
        if (RunPart(3) && rawT && rawT->GetAlias("raw_0") && rawT->GetEntries() > 0){
            std::vector<std::pair<std::string, double>> infoCuts = {{"1_hitZDC", -1.0}, {"0_cuts/1_noCuts", -1.0}, {"0_cuts/2_theta0_1", -1.0}}; // (folder, KE cut in MeV, -1 none)
            infoCuts.push_back({KeCutDir(80), initialKE * 0.80});
            for (const auto& sc : sampleCuts) infoCuts.push_back({sc.dir, -1.0}); // no KE cut: what each cut does to every variable
            std::vector<SampleCut> keSamples; // the 85, 90 and 95% KE levels of 7_megaCut2/1_1particlecut and 9_megaCut3/1_1particlecut: the sample of the 80% folder with its KE cut file (cut_megacut<2,3>_1_ke<pct>.root, CutPlots.cpp (samples)), for their 1_additionalinfo
            for (const auto& c : sampleCuts) if (c.dir == "7_megaCut2/1_1particlecut" || c.dir == "9_megaCut3/1_1particlecut") for (const int pct : {85, 90, 95}){
                SampleCut k = c;
                k.dir += Form("_KE%d", pct); k.file += Form("_ke%d", pct);
                if (const size_t p = k.label.find("80%"); p != std::string::npos) k.label.replace(p, 3, std::to_string(pct) + "%");
                keSamples.push_back(k);
            }
            for (const auto& k : keSamples) infoCuts.push_back({k.dir, -1.0}); // run with KECUT80 (P3_CUT 0, P3_SUB 3), like the KECUT85 to 95 ones
            for (int pct : {85, 90, 95}) infoCuts.push_back({KeCutDir(pct), initialKE * pct / 100.0}); // run with the 80% one (P3_CUT 0, P3_SUB 3)
            auto lapStart = std::chrono::steady_clock::now();
            const auto lap = [&](const std::string& what){ // seconds since the last lap
                const auto now = std::chrono::steady_clock::now();
                printf(">>> 1_additionalinfo %s took %.1f s\n", what.c_str(), std::chrono::duration<double>(now - lapStart).count()); fflush(stdout);
                lapStart = now;
            };
            for (const auto& ic : infoCuts){
                if (const int ii = int(&ic - &infoCuts[0]), kx = ii >= 4 + int(sampleCuts.size()); !RunCut(ii < 4 || kx ? 0 : ii - 3) || ((ii < 4 || kx) && !RunSub(kx ? 3 : ii))) continue; // P3_CUT: 0 the first four folders (P3_SUB picks one) and the KECUT85 to 95 ones (kx, with KECUT80, sub 3), i the sampleCuts[i - 1] folder
                lapStart = std::chrono::steady_clock::now();
                const SampleCut* sc = nullptr;
                for (const auto& c : sampleCuts) if (c.dir == ic.first) sc = &c;
                for (const auto& c : keSamples) if (c.dir == ic.first) sc = &c;
                const bool noTheta = ic.first == "0_cuts/1_noCuts"; // no cut at all, not even the ZDC hit: every fitted proton (AnglesRaw). Only 1_real_vs_recon_KEDetector (part 1) there needs the ZDC hit
                const bool nz = sc != nullptr || noTheta; // the cuts' plots and 0_cuts/1_noCuts use every fitted proton (AnglesRaw), the others the protons that hit the ZDC
                TTree *tCE = nz ? (sc ? cutTree(*sc) : rawT) : events2, *tAll = nz ? tCE : eventsAllEv, *tPart = nz ? tCE : eventsAllPart; // tPart: every recon charged particle (any pdg)
                const std::string ke = ic.second < 0 ? "" : " && fDEnergy > " + std::to_string(ic.second);
                const bool inCuts0 = ic.first.rfind("0_cuts/", 0) == 0; // its plots sit in the folder itself, not in a 1_additionalinfo
                const std::string zb = noTheta ? zdcRecoAll : zdcBase, zbp = noTheta ? zdcRecoPartAll : zdcBasePart, zNote = noTheta ? kNoCutNote : kThetaNote;
                const std::string bCE = nz ? rawCE : zb, bAll = nz ? "PDG == 2212" : zb, bPart = nz ? "1" : zbp, hitTxt = nz ? "" : " that hit ZDC"; // all recon protons are the pdg 2212 ones
                // every histogram of the folder is booked first and filled in one event loop per tree (one for the three samples of a cut, they share the cut file), not one loop per plot
                std::vector<std::pair<std::string, std::string>> reqCE, reqAll, reqPart; // (Draw string, cut) of the CE recon protons, of all recon protons and of all recon charged particles
                std::vector<int> slot(kInfoVars.size(), -1); // position of the variable in reqCE and reqAll
                std::vector<std::pair<double, double>> range(kInfoVars.size());
                for (size_t iv = 0; iv < kInfoVars.size(); iv++){
                    if (!RunSec(int(iv % 4))) continue; // P3_SEC 0 to 3: the variables of one group
                    const auto& v = kInfoVars[iv];
                    const std::string k = std::to_string(v.k);
                    slot[iv] = (int)reqCE.size();
                    double lo = v.lo, hi = v.hi;
                    if (lo >= hi){
                        const auto [q1, q99] = quantiles(v.k); (void)q1;
                        lo = 0.0; hi = std::max(1e-6, 1.2 * q99);
                    }
                    range[iv] = {lo, hi};
                    // 0_cuts/1_noCuts nChargedTrackerHitsEv: only the events with at most 1 charged particle in every tracker plane
                    const bool perPlane = noTheta && v.name == "nChargedTrackerHitsEv";
                    const std::string sel = ke + (perPlane ? perPlaneCut : std::string()) + " && rawOk(fEvent, TrackID, " + k + ")",
                        exprPre = "rawVal(fEvent, TrackID, " + k + ")>>hInfo", spec = Form("(%d,%g,%g)", v.nBins, lo, hi);
                    reqCE.push_back({P3(exprPre + "CE" + k + spec), P3(bCE + sel)});
                    reqAll.push_back({P3(exprPre + "Recon" + k + spec), P3(bAll + sel)});
                    reqPart.push_back({P3(exprPre + "Part" + k + spec), P3(bPart + sel)});
                }
                lap(ic.first + ": quantiles and booking the requests");
                std::vector<TH1*> hCEs, hAlls, hParts;
                if (reqCE.empty()){ // P3_SEC 4 or 5: another process does the histograms
                } else if (tCE == tAll){ // a cut folder: the three samples are rows of one cut file, one event loop
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
                    // the residual plot of 3_ResCuts shows the sample before its own cut: every fitted proton of the theta baseline (AnglesRaw), the cut is the dashed line
                    const bool uncut = sc && sc->dir == "3_ResCuts" && v.folder == sc->onFolder;
                    if (uncut){
                        delete hCE; delete hRec; delete hPart;
                        const std::string uk = std::to_string(v.k), usel = " && " + thetaBase + " && rawOk(fEvent, TrackID, " + uk + ")", uex = "rawVal(fEvent, TrackID, " + uk + ")>>hInfoUncut",
                            uspec = Form("(%d,%g,%g)", v.nBins, range[iv].first, range[iv].second);
                        const std::vector<TH1*> us = HistsMT(rawT, {{P3(uex + "CE" + uk + uspec), P3(rawCE + usel)}, {P3(uex + "Recon" + uk + uspec), P3("PDG == 2212" + usel)}, {P3(uex + "Part" + uk + uspec), P3("1" + usel)}});
                        hCE = us[0]; hRec = us[1]; hPart = us[2];
                        if (!hCE || !hRec || !hPart) continue;
                    }
                    canvas->Clear(); canvas->SetLogy(1);
                    hCE->SetLineColor(kRed); hRec->SetLineColor(kBlue); hPart->SetLineColor(kGreen + 2); hCE->SetLineWidth(2); hRec->SetLineWidth(2); hPart->SetLineWidth(2); hCE->SetStats(0); hRec->SetStats(0); hPart->SetStats(0);
                    hRec->SetTitle((v.label + ", CE recon protons vs all recon protons vs all recon charged particles" + (nz ? std::string(", ") : std::string(", hit ZDC, " + zNote + ", ")) + (ic.second < 0 ? std::string(sc ? "no KE cut, " + (uncut ? std::string("no residual cut, dashed line at the 3_ResCuts cut resDRms < 4 mm, ") + kThetaNote : sc->label) + ", no ZDC requirement" : "no KE cut") : "ZDC KE > " + std::to_string(int(ic.second)) + " MeV") + (noTheta && v.name == "nChargedTrackerHitsEv" ? ", at most 1 charged particle per tracker plane" : "") + ";" + v.label + ";Count").c_str());
                    if (v.name == "nHits") for (int bin = 1; bin <= hRec->GetNbinsX(); bin++) hRec->GetXaxis()->SetBinLabel(bin, bin % 2 == 1 && bin > 1 ? Form("%d", bin - 1) : ""); // labels at the even counts only
                    hRec->SetMinimum(0.5); hRec->SetMaximum(2.0 * std::max(std::max(hRec->GetMaximum(), hCE->GetMaximum()), hPart->GetMaximum()));
                    hRec->Draw("hist"); hPart->Draw("hist same"); hCE->Draw("hist same");
                    if (uncut){TLine* rl = new TLine(kRes4Mm, 0.5, kRes4Mm, hRec->GetMaximum()); rl->SetLineStyle(2); rl->Draw();}
                    TLegend* leg = new TLegend(0.745, 0.78, 0.99, 0.90); leg->SetTextSize(0.014);
                    leg->AddEntry(hCE, Form("CE recon protons (%.0f)", hCE->GetEntries()), "l");
                    leg->AddEntry(hRec, Form("all recon protons (%.0f)", hRec->GetEntries()), "l");
                    leg->AddEntry(hPart, Form("all recon charged particles (%.0f)", hPart->GetEntries()), "l");
                    leg->Draw();
                    // flat files named <folder>_<variable> in 1_additionalinfo of the hit ZDC, KE cut and 2 mm cut folders
                    const std::string idir = folder_KEcuts + ic.first + (sc && v.folder == sc->onFolder || inCuts0 ? "/" : "/1_additionalinfo/"); // the variable a cut is on sits beside 1_additionalinfo; 0_cuts has no such folder, its plots sit in it directly
                    gSystem->mkdir(idir.c_str(), true);
                    canvas->SaveAs((idir + v.folder + "_" + v.name + ".pdf").c_str());
                    if (noTheta && v.name == "nParticlesEv") canvas->SaveAs((folder_KEcuts + "4_nParticles/4_nParticles_nParticlesEv.pdf").c_str()); // the general count, no cut but the ZDC hit, beside the 1 and 2 particle cuts
                    if (noTheta && v.name == "nChargedTrackerHitsEv"){gSystem->mkdir((folder_KEcuts + "10_nChargedTrackerHits").c_str(), true); canvas->SaveAs((folder_KEcuts + "10_nChargedTrackerHits/10_nChargedTrackerHits_nChargedTrackerHitsEv.pdf").c_str());} // the general count with no cut, in its own folder
                    if (noTheta && v.name == "nHits"){gSystem->mkdir((folder_KEcuts + "5_nHits").c_str(), true); canvas->SaveAs((folder_KEcuts + "5_nHits/5_nHits_nHits_noCut.pdf").c_str());} // the general count with no cut, beside the cut version
                    if (noTheta && v.name == "nZdcCharged"){gSystem->mkdir((folder_KEcuts + "8_nZdcParticles").c_str(), true); canvas->SaveAs((folder_KEcuts + "8_nZdcParticles/8_nZdcParticles_nZdcCharged.pdf").c_str());} // the general count with no cut, beside the exactly 1 ZDC particle cut
                    canvas->SetLogy(0);
                    delete hCE; delete hRec; delete hPart;
                }
                const std::string cutNote = (sc ? sc->label + (sc->zdcKE80 ? ", hit ZDC" : ", no ZDC requirement") + ", " + kThetaNote : noTheta ? std::string("no cut, no ZDC requirement, ") + zNote
                                                : std::string("hit ZDC, ") + zNote + (ic.second < 0 ? std::string(", no KE cut") : std::string(Form(", ZDC KE > %.0f%%", 100.0 * ic.second / initialKE)))); // every cut of the folder, for the titles of 9 to 12
                if (RunSec(4)){ // 9, 10: block theta and phi resolution of the CE recon protons of the folder: recon block angle minus the true CE angle (CalcValsReDDIFF of cut_APandT.root), only CE protons have a truth
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
                    struct DiffPlot {std::string name, label; double lim; int nBins;};
                    const DiffPlot dps[2] = {{"thetaBlockDiff", "recon block theta - true CE theta (degrees)", 1.0, 80}, {"phiBlockDiff", "recon block phi - true CE phi (degrees)", 180.0, 72}};
                    for (int w = 0; w < 2; w++){
                        TH1D hd(Form("hBlockDiff%d", w), "", dps[w].nBins, -dps[w].lim, dps[w].lim);
                        hd.SetDirectory(nullptr);
                        std::vector<double> diffVals; // the values the unbinned fit takes, the histogram only shows them
                        for (size_t i = 0; i < keys[0].size(); i++){
                            const auto it = blockDiff.find({(int)keys[0][i], (int)keys[1][i]});
                            if (it != blockDiff.end()){const double dv = w ? it->second.second : it->second.first; hd.Fill(dv); diffVals.push_back(dv);}
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
                            // unbinned likelihood fit of a gaussian to the values in -lim to lim (FitCore.h), height in counts per bin so it overlays the histogram
                            std::vector<double> inRange;
                            for (double dv : diffVals) if (dv >= -dps[w].lim && dv <= dps[w].lim) inRange.push_back(dv);
                            double dMean = 0, dVar = 0;
                            for (double dv : inRange) dMean += dv / inRange.size();
                            for (double dv : inRange) dVar += (dv - dMean) * (dv - dMean) / inRange.size();
                            const double dSig = std::max(std::sqrt(dVar), 0.01 * dps[w].lim);
                            UGaussModel gm; gm.lo = -dps[w].lim; gm.hi = dps[w].lim; gm.pLo = {-dps[w].lim, 0.002 * dps[w].lim}; gm.pHi = {dps[w].lim, 2 * dps[w].lim};
                            const UGaussFit gr = FitUGaus(inRange, gm, {dMean, dSig}, {0.1 * dSig, 0.1 * dSig});
                            if (gr.ok){
                                TF1* fd = AddGausCurve(&hd, Form("fBlockDiffDraw%d", w), gr, gm);
                                hd.GetListOfFunctions()->Remove(fd); // drawn by hand below ("hist" does not draw the function list)
                                fd->SetLineColor(kBlack); fd->SetLineWidth(2); fd->Draw("same");
                                leg->AddEntry(fd, Form("Unbinned gaussian fit, mean %.3g, sigma %.3g", fd->GetParameter(1), fd->GetParameter(2)), "l");
                            }
                        }
                        leg->Draw();
                        canvas->Modified(); canvas->Update();
                        const std::string idir = folder_KEcuts + ic.first + (inCuts0 ? "/" : "/1_additionalinfo/");
                        gSystem->mkdir(idir.c_str(), true);
                        canvas->SaveAs((idir + std::to_string(9 + w) + "_blockDiff_" + dps[w].name + ".pdf").c_str());
                    }
                }
                if (RunSec(5)){ // 11, 12: recon minus recon block, theta and phi, CE recon protons (red) and all recon protons (blue) of the folder
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
                        ha.SetTitle((rps[w].label + ", CE recon protons vs all recon protons vs all recon charged particles, " + cutNote + ";" + rps[w].label + ";Count").c_str());
                        if (w == 1) ha.GetXaxis()->SetNdivisions(508, kFALSE); // phi: -180 to 180 every 45 deg
                        ha.SetMinimum(0.5); ha.SetMaximum(2.0 * std::max(1.0, std::max(std::max(ha.GetMaximum(), hc.GetMaximum()), hq.GetMaximum()))); // log y
                        ha.Draw("hist"); hq.Draw("hist same"); hc.Draw("hist same");
                        TLegend* leg = new TLegend(0.745, 0.76, 0.99, 0.90); leg->SetTextSize(0.011);
                        leg->AddEntry(&hc, Form("CE recon protons (%.0f), mean %.3g, RMS %.3g", hc.GetEntries(), hc.GetMean(), hc.GetRMS()), "l");
                        leg->AddEntry(&ha, Form("all recon protons (%.0f), mean %.3g, RMS %.3g", ha.GetEntries(), ha.GetMean(), ha.GetRMS()), "l");
                        leg->AddEntry(&hq, Form("all recon charged particles (%.0f), mean %.3g, RMS %.3g", hq.GetEntries(), hq.GetMean(), hq.GetRMS()), "l");
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
