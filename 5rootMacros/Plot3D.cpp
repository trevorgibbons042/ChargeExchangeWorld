#include <TColor.h>
#include <TH2D.h>
#include <TMarker.h>
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

#include "CutOpen.h"
#include "FitCore.h"
#include "Plot3DCols.h"
#include "QualityCut.h"
#include "TreeColumns.h"

// The maps (gRecoTracks, gRawVals, ...), isReco/rawVal/rawOk and the columns (cut_ZdcCE.root, cut_ZdcAllEv.root) or aliases (AnglesRaw) that replace them in the Draw cuts (P3) are in Plot3DCols.h.

struct SampleCut {std::string dir, expr, label, onFolder, file;}; // see sampleCuts in Plot3D

// Plot3D has four independent parts that each write their own plots/ folders: 1 position, energy, 2_KECUT80 and 1_hitZDC; 2 the cut ratios and rank csv files of the 6_cuts folders;
// 3 their 1_additionalinfo plots; 4 their phi fits. P3_PART=n runs part n only, so the four can run as four processes side by side (Plot.pbs, plot.sh). Without it all four run.
// Part 2 is one process per cut with P3_CUT=i (0: the 1_hitZDC baseline, 1 to 5: sampleCuts in order), each writing its rank rows to cuts/rank_rows_<i>.tsv; part 5 reads those
// files and writes the rank_*.csv files. Parts 3 and 4 also take P3_CUT=i to do one cut folder only (0: 1_hitZDC and 2_KECUT80). Without P3_PART all of it runs in one process.
static bool RunPart(int n){const char* e = std::getenv("P3_PART"); return !e || !*e || std::atoi(e) == n;}
static bool RunCut(int i){const char* e = std::getenv("P3_CUT"); return !e || !*e || std::atoi(e) == i;}

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
    if (!events2 || !eventsAllEv){Error("Plot3D", "missing cuts/cut_ZdcCE.root or cut_hitZDC_AllEv.root (run CutPlots.cpp first)"); return;}


    if (Graph == 1){

        TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
        gPad->SetRightMargin(0.30);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); // room for the rotated Y title and bottom tick labels
        gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X"); // stats/fit box to the right, fewer x labels

        // Every plot goes through this: draw expr from t with cut, set the title (title;xaxis;yaxis), let tweak adjust the
        // histogram (range, fit, ...), save to path. opt is the Draw option, e.g. "colz".
        std::set<std::string> fitFilesStarted; // stats and fit files already truncated in this run
        std::string titleSuffix; // added to the title text (before the axis titles), the KE cut of the position plots
        double nCERecoAll = -1, nCEZdcAll = -1, nRecoAll = -1; // CE protons at the ZDC before any cut (with a reco track, and all), and every reco proton there, counted on the first plot that shows them
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
            if (tweak) tweak(h);
            const std::string hName = h->GetName(); // the stats box shows the name: a readable one until the pdf is saved
            h->SetName(StatsTitle(h, title.substr(0, std::min(title.find(','), axes)) + (t == events2 ? ", CE protons" : t == eventsAllEv ? ", all recon protons" : "")).c_str()); // quantity: the title up to its first comma
            if (t == eventsAllEv || t == events2) countZdcAll(); // the ZDC totals the percent boxes use, counted once
            if (t == eventsAllEv){ // all recon protons: the CE ones (cut_ParentZdc.root) under the same cut
                canvas->cd();
                DrawCEPercent(CountMT(events2, P3(cut)), h->GetEntries(), nCERecoAll, nRecoAll); // lost: vs every reco CE proton at the ZDC
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
        // x:y position maps: log color scale from 1 count up to the fullest bin
        auto drawXY = [&](TTree* t, const std::string& expr, const std::string& cut, const std::string& title, const std::string& path){
            canvas->SetLogz(1);
            drawPlot(t, expr, cut, title, path, "colz", [](TH1* h){h->SetMinimum(1);});
            canvas->SetLogz(0);
        };

        const std::string zdcProtons = "boundaryType == 2 && pdg == 2212";
        // position plots: only protons with a reconstructed track
        const std::string zdcReco = zdcProtons + " && isReco(fEvent, TrackID)", zdcRecoAll = zdcReco;
        // baseline of every plot and count in 6_cuts: recon theta from 0 to 1 degree (NaN fails), on top of the ZDC hit and the cut
        const std::string thetaBase = "rawVal(fEvent, TrackID, 9) >= 0 && rawVal(fEvent, TrackID, 9) <= 1", zdcBase = zdcRecoAll + " && " + thetaBase;
        countZdcAll = [&]{
            if (nCERecoAll >= 0) return;
            nCERecoAll = CountMT(events2, P3(zdcRecoAll)); nCEZdcAll = CountMT(events2, P3(zdcProtons)); nRecoAll = CountMT(eventsAllEv, P3(zdcRecoAll));
        };
        const std::string kThetaNote = "recon theta 0-1 deg"; // the thetaBase cut, named in the titles of the 6_cuts plots
        const std::string zdcRecoDAll = zdcProtons + " && isRecoD(fEvent, TrackID)"; // block position plots: also the block residual cut

        // Position plots have no KE cut: 1_position/1_true/ and 1_position/2_block/
        std::vector<std::pair<std::string, std::string>> posKes = {{"", ""}}; // unused folder, extra cut
        std::vector<std::string> posKeLabels = {""};
        if (RunPart(1)) for (size_t pk = 0; pk < posKes.size(); pk++){
            const std::string folder_pos = folder + "/1_position/1_true/", folder_posD = folder + "/1_position/2_block/";
            const std::string zdcReco = zdcRecoAll + posKes[pk].second, zdcRecoD = zdcRecoDAll + posKes[pk].second;
            const std::string xySlice = zdcReco + " && fabs(fY) < 2", xySliceD = zdcRecoD + " && fabs(fDy) < 2";
            titleSuffix = ", has recon track, no KE cut" + posKeLabels[pk]; // the cut of the 1_true plots: hit ZDC and a recon track
            // ZDC position, true (fX, fY)
            drawXY(eventsAllEv, "fY:fX>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at ZDC, all recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_pos + "2_xy_recon.pdf");
            drawXY(eventsAllEv, "fY:fX>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at ZDC, all recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_pos + "4_xy_zoom_recon.pdf");
            drawXY(events2, "fY:fX>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at ZDC, CE recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_pos + "1_xy_CE.pdf");
            drawXY(events2, "fY:fX>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at ZDC, CE recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_pos + "3_xy_zoom_CE.pdf");

            gStyle->SetOptFit(1111); // probability, chi2/ndf, fit parameters with errors
            drawPlot(events2, "fX>>hSliceX200(100,-200,200)", xySlice, "X Position at ZDC, CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_pos + "5_x_fit_200_CE.pdf", "", [](TH1* h){FitCore(h);});
            drawPlot(events2, "fX>>hSliceX(100,-50,50)", xySlice, "X Position at ZDC, CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_pos + "7_x_fit_50_CE.pdf", "", [](TH1* h){FitCore(h);});
            drawPlot(eventsAllEv, "fX>>hSliceXAllEv200(100,-200,200)", xySlice, "X Position at ZDC, all recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_pos + "6_x_fit_200_recon.pdf", "", [](TH1* h){FitCore(h);});
            drawPlot(eventsAllEv, "fX>>hSliceXAllEv(100,-50,50)", xySlice, "X Position at ZDC, all recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_pos + "8_x_fit_50_recon.pdf", "", [](TH1* h){FitCore(h);});


            titleSuffix = ", has recon track, block residual < 2 mm, no KE cut" + posKeLabels[pk]; // the cut of the 2_block plots: isRecoD
            // ZDC position, block (fDx, fDy): the ZDC is nCols_detector x nRows_detector blocks over 400 mm. Every bin holds the same whole number of blocks (BlockBins, FitCore.h):
            // 100 bins over 100 mm on a pitch of 400 / 801 mm put an extra block in the bin that starts at x = 0 (a cross at x = 0 and y = 0), and in one bin of ten along the x slices (spikes).
            const int nColsD = (int)ReadSettings().at("nCols_detector"), nRowsD = (int)ReadSettings().at("nRows_detector");
            const std::string binsFull = BlockBins(nColsD, 4.0, 200), binsFullY = BlockBins(nRowsD, 4.0, 200), binsZoom = BlockBins(nColsD, 1.0, 50), binsZoomY = BlockBins(nRowsD, 1.0, 50), binsSlice200 = binsFull, binsSlice50 = BlockBins(nColsD, 2.5, 127.27);
            drawXY(eventsAllEv, "fDy:fDx>>h2(" + binsFull + "," + binsFullY + ")", zdcRecoD, "Position at ZDC BLOCK, all recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_posD + "2_xy_recon.pdf");
            drawXY(eventsAllEv, "fDy:fDx>>h2(" + binsZoom + "," + binsZoomY + ")", zdcRecoD, "Position at ZDC BLOCK, all recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_posD + "4_xy_zoom_recon.pdf");
            drawXY(events2, "fDy:fDx>>h2(" + binsFull + "," + binsFullY + ")", zdcRecoD, "Position at ZDC BLOCK, CE recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_posD + "1_xy_CE.pdf");
            drawXY(events2, "fDy:fDx>>h2(" + binsZoom + "," + binsZoomY + ")", zdcRecoD, "Position at ZDC BLOCK, CE recon protons that hit ZDC;X Position (mm);Y Position (mm)", folder_posD + "3_xy_zoom_CE.pdf");

            drawPlot(events2, "fDx>>hSliceDX200(" + binsSlice200 + ")", xySliceD, "X Position at ZDC BLOCK, CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_posD + "5_x_fit_200_CE.pdf", "", [](TH1* h){FitBlocks(h);}); //whole blocks per bin
            drawPlot(events2, "fDx>>hSliceDX(" + binsSlice50 + ")", xySliceD, "X Position at ZDC BLOCK, CE recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_posD + "7_x_fit_50_CE.pdf", "", [](TH1* h){FitBlocks(h);});
            drawPlot(eventsAllEv, "fDx>>hSliceDXAllEv200(" + binsSlice200 + ")", xySliceD, "X Position at ZDC BLOCK, all recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_posD + "6_x_fit_200_recon.pdf", "", [](TH1* h){FitBlocks(h);});
            drawPlot(eventsAllEv, "fDx>>hSliceDXAllEv(" + binsSlice50 + ")", xySliceD, "X Position at ZDC BLOCK, all recon protons that hit ZDC, |Y| < 2 mm;X Position (mm);Count", folder_posD + "8_x_fit_50_recon.pdf", "", [](TH1* h){FitBlocks(h);});
            titleSuffix.clear();
        }

        // KE x axes: 0 to 30k MeV for the 3_KEcuts plots, 10k to 30k MeV for the 2_energy plots
        struct KERange {double lo, hi; std::string tag;};
        const std::vector<KERange> keRanges = {{0, 30000, "_0to30k"}}; // 3_KEcuts: full range only
        const std::vector<KERange> energyRanges = {{10000, 30000, "_10to30k"}}; // 2_energy: 10k to 30k only
        auto binSpec = [](const KERange& r){return std::string(Form("(100,%g,%g)", r.lo, r.hi));};

        // 2_energy/ holds the energy plots, with the KE cut folders ke_cut_<pct>/ next to them
        if (RunPart(1)) {
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
            // CE protons (red) over the recon protons (black) on one log y plot: 1_true_ke, 2_smeared_ke
            for (const auto& q : keQtys){
                for (const auto& r : energyRanges){
                    drawPlot(eventsAllEv, q.branch + ">>hKE" + binSpec(r), zdcRecoE, q.title + ", all recon protons that hit ZDC and CE protons that hit ZDC (reconstructed or not);" + q.axis + ";Count",
                             folder_smear + std::to_string(&q - &keQtys[0] + 1) + "_" + q.file + r.tag + ".pdf", "", [&](TH1* h){
                        TH1* hCE = HistMT(events2, P3(q.branch + ">>hKECE" + binSpec(r)), P3(zdcProtonsE)); hCE->Draw("same");
                        hCE->SetLineColor(kRed); hCE->SetStats(0);
                        h->SetMinimum(0.5); canvas->SetLogy(1);
                        TLegend* leg = new TLegend(0.745, 0.30, 0.99, 0.38); leg->SetTextSize(0.015); // below the stats box
                        leg->AddEntry(h, "all recon protons that hit ZDC", "l");
                        leg->AddEntry(hCE, "CE protons that hit ZDC (reconstructed or not)", "l");
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
                drawPlot(eventsAllEv, "fDEnergy-fKEnergy>>hSmearDiffAllEv(100,-10000,10000)", zdcRecoE, "Energy Difference at ZDC, all recon protons that hit ZDC;fDEnergy - fKEnergy (MeV);Count", folder_smear + "6_diff_recon.pdf", "", [](TH1* h){FitCore(h);});

                // percent error = 100 * (fDEnergy - fKEnergy) / fKEnergy
                const std::string pct = "100*((fDEnergy-fKEnergy)/fKEnergy)"; // bracketed: the C++ of HistMT would divide 100*(a-b) by c, TTree::Draw 100*(a-b/c)
                drawPlot(events2, pct + ">>hSmearPctCE(100,-100,100)", zdcProtonsE, "Percent Error at ZDC, CE protons that hit ZDC (reconstructed or not);Percent Error (%);Count", folder_smear + "7_pct_CE.pdf", "", [](TH1* h){FitCore(h);});
                drawPlot(eventsAllEv, pct + ">>hSmearPctAllEv(100,-100,100)", zdcRecoE, "Percent Error at ZDC, all recon protons that hit ZDC;Percent Error (%);Count", folder_smear + "8_pct_recon.pdf", "", [](TH1* h){FitCore(h);});

            }
            else{Warning("Plot3D", "no fDEnergy branch in cut_CE.root, skipping smearing plots");}
            titleSuffix.clear();
        }

        const double initialKE = kBeamKE; // MeV, beam KE (FitCore.h)
        // KE cuts on the ZDC-smeared energy (fDEnergy) only, in ke_cut_<pct>/. Files from before fDEnergy was added get no KE cuts.
        struct CutEnergy {std::string branch, dir, label;};
        std::vector<CutEnergy> cutEnergies;
        if (events->GetBranch("fDEnergy")) cutEnergies.push_back({"fDEnergy", "KEDetector", "ZDC KE"});
        if (RunPart(1)) for (int pct : {80}){
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
                    h->SetLineColor(kBlue);
                    TH1* hf = HistMT(eventsAllEv, P3(e.branch + ">>hCutFakeOv" + binSpec(r)), P3(zdcReco)); hf->Draw("same");
                    hf->SetLineColor(kRed);
                    h->SetMinimum(0.5); h->SetMaximum(2.0 * std::max(h->GetMaximum(), hf->GetMaximum())); canvas->SetLogy(1); // log y, full x range of the file name, the cut is the dashed line
                    TLine* cl = new TLine(cutVal, 0.5, cutVal, h->GetMaximum()); cl->SetLineStyle(2); cl->Draw();
                    TLegend* leg = new TLegend(0.745, 0.44, 0.99, 0.52); leg->SetTextSize(0.015); // below the stats box
                    leg->AddEntry(h, Form("CE recon protons that hit ZDC, %lld / %lld = %.1f%% left out", nRealOut, nRealTotal, realOutPct), "l");
                    leg->AddEntry(hf, Form("all recon protons that hit ZDC, %lld / %lld = %.1f%% kept in", nFakeKept, nFakeTotal, fakeKeptPct), "l");
                    leg->Draw();
                });
                canvas->SetLogy(0);
                }
            }
        }

        // Cuts on the reconstructed protons, no KE cut, each in its own folder under 6_cuts: dir, extra Draw cut, label, and the info variable folder the cut is on
        // (those plots sit in the cut folder itself, the other variables in its 1_additionalinfo). The same KE plot as 1_hitZDC with the cut applied goes in 1_additionalinfo.
        const std::vector<SampleCut> sampleCuts = {
            {"3_ResCuts", "isRes2(fEvent, TrackID)", "cut resDRms (block residual) < 2 mm", "3_ResCuts", "cut_ResCuts"},
            {"4_nProtons/1_protonCut", "rawVal(fEvent, TrackID, 6) <= 1", "cut at most 1 proton with a track in the event", "4_nProtons", "cut_nProtons1"},
            {"4_nProtons/2_protonCut", "rawVal(fEvent, TrackID, 6) <= 2", "cut at most 2 protons with a track in the event", "4_nProtons", "cut_nProtons2"},
            {"8_nHits", "(rawVal(fEvent, TrackID, 10) == 4 || rawVal(fEvent, TrackID, 10) == 8 || rawVal(fEvent, TrackID, 10) == 12 || rawVal(fEvent, TrackID, 10) == 16)", "cut 4, 8, 12 or 16 tracker hits", "8_nHits", "cut_nHits"},
            {"11_megaCut/1_megaCut", "rawVal(fEvent, TrackID, 6) <= 1 && (rawVal(fEvent, TrackID, 10) == 4 || rawVal(fEvent, TrackID, 10) == 8 || rawVal(fEvent, TrackID, 10) == 12 || rawVal(fEvent, TrackID, 10) == 16) && rawVal(fEvent, TrackID, 9) >= 0 && rawVal(fEvent, TrackID, 9) <= 1",
             "MEGA CUT: 1 Recon Proton in Event AND Number of times Recon Proton hits tracker [4,8,12,16], theta from 0 degrees to 1 degree", "11_megaCut", "cut_megacut"}};
        struct RankRow {std::string cut, sample; double nCEin, nAllIn, nCEtotal, nCEreconBefore, fom[2], pol[2], polErr[2];}; // fom, pol (fitted P), polErr: [0] recon theta / phi, [1] recon block theta / phi
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
        if (rawT && (RunPart(2) || RunPart(4) || RunPart(5))) nCEreconTotal = CountMT(rawT, P3("isFromCE==1 && " + thetaBase));
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
                TF1 fit("fPhiCosCut", "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                fit.SetParNames("N0", "P", "A");
                fit.FixParameter(0, h->GetEntries() / 50); fit.FixParameter(2, aMean); fit.SetParameter(1, 0.5); fit.SetParLimits(1, -1, 1);
                h->Fit(&fit, path.empty() ? "QN" : "Q");
                h->SetMinimum(0); h->SetMaximum(1.1 * h->GetMaximum());
                *pol = fit.GetParameter(1); *polErr = fit.GetParError(1);
                if (!path.empty()){
                    h->SetName(StatsTitle(h).c_str()); // stats box title; h is deleted below
                    canvas->Update();
                    if (auto* st = (TPaveStats*)h->FindObject("stats")){st->SetX1NDC(0.76); st->SetX2NDC(0.995); st->SetY1NDC(0.55); st->SetY2NDC(0.83);}
                    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
                    tx.DrawLatex(0.13, 0.90, "Fit  N(#phi) = N_{0} (1 + P #upoint #LTA#GT cos#phi),   N_{0} = Entries / bins (fixed)");
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
        if (RunPart(2) && RunCut(0)) for (const auto& e : cutEnergies) for (const bool ke80 : {false, true}){ // the uncut sample (1_hitZDC) is in every ranking as the baseline
            const std::string q = zdcBase + (ke80 ? " && " + e.branch + " > " + std::to_string(initialKE * 0.80) : "");
            double fom[2], pol[2], polErr[2]; fomStats(events2, q, fom);
            for (int j = 0; j < 2; j++) phiFit(events2, q, j, "", "", &pol[j], &polErr[j]);
            rankRows.push_back({"1_hitZDC (no cut)", ke80 ? "cut and ZDC KE > 80%" : "cut", double(CountMT(events2, P3(q))), double(CountMT(eventsAllEv, P3(q))), double(nCEallZdcAll),
                                nCEreconTotal, {fom[0], fom[1]}, {pol[0], pol[1]}, {polErr[0], polErr[1]}});
        }
        for (const auto& sc : sampleCuts) for (const auto& e : cutEnergies) for (const bool ke80 : {false, true}){ // 0_KE_hitZDC (the cut and hitting the ZDC) and 0_KE_KECUT80 (also the 80% KE cut)
            if (!RunPart(2) || !RunCut(int(&sc - &sampleCuts[0]) + 1)) continue;
            const std::string folder_res = folder_KEcuts + sc.dir + "/1_additionalinfo/";
            const double cutVal = initialKE * 0.80; // the 80% KE cut of 2_KECUT80
            const std::string qCut = zdcBase + " && " + sc.expr + (ke80 ? " && " + e.branch + " > " + std::to_string(cutVal) : "");
            const long long nCEq = CountMT(events2, P3(qCut)), nAllq = CountMT(eventsAllEv, P3(qCut));
            { // cut_ratios.csv in the cut folder, one row without and one with the 80% KE cut: CE recon protons in the cut (no ZDC requirement without the 80% KE cut) over all recon protons in the cut,
              // over all CE protons the sim made (no reconstruction or ZDC requirement), and over all CE recon protons that hit ZDC before the cut
                const long long nCEallZdc = nCEallZdcAll; // every charge exchange the sim made, one CE proton each
                gSystem->mkdir((folder_KEcuts + sc.dir).c_str(), true);
                FILE* rf = fopen((folder_KEcuts + sc.dir + "/cut_ratios.csv").c_str(), ke80 ? "a" : "w");
                if (!ke80) fprintf(rf, "sample,CE_recon_proton_in_cut,recon_proton_in_cut,CE_recon_proton_in_cut/recon_proton_in_cut,CE_protons_total,CE_recon_before_cut,CE_recon_proton_in_cut/CE_recon_before_cut,nonCE_recon_proton_in_cut,nonCE_recon_proton_in_cut/recon_proton_in_cut,FOM_theta,FOM_thetaBlock,P_phiRe,P_phiRe_err,P_phiReD,P_phiReD_err\n");
                // without the 80% KE cut no ZDC hit is required (AnglesRaw); the ZDC KE cut needs the ZDC, so that row is CE recon protons that hit ZDC
                TTree* tS = ke80 ? events2 : cutTree(sc); // the cut file holds only the rows that pass the cut
                const std::string qS = ke80 ? qCut : "1", qSCE = ke80 ? qCut : rawCE;
                const long long nCEc = CountMT(tS, P3(qSCE)), nAllc = ke80 ? nAllq : CountMT(tS, P3(qS));
                double fom[2], pol[2], polErr[2]; fomStats(tS, qSCE, fom);
                for (int j = 0; j < 2; j++) phiFit(tS, qSCE, j, "", "", &pol[j], &polErr[j]);
                fprintf(rf, "%s,%lld,%lld,%.6f,%lld,%.0f,%.6f,%lld,%.6f,%.6g,%.6g,%.5f,%.5f,%.5f,%.5f\n", ke80 ? "cut and ZDC KE > 80%" : "cut", nCEc, nAllc, double(nCEc) / std::max<long long>(1, nAllc),
                        nCEallZdc, nCEreconTotal, nCEc / std::max(1.0, nCEreconTotal),
                        nAllc - nCEc, double(nAllc - nCEc) / std::max<long long>(1, nAllc), fom[0], fom[1], pol[0], polErr[0], pol[1], polErr[1]);
                fclose(rf);
                rankRows.push_back({sc.dir, ke80 ? "cut and ZDC KE > 80%" : "cut", double(nCEc), double(nAllc), double(nCEallZdc), nCEreconTotal, {fom[0], fom[1]}, {pol[0], pol[1]}, {polErr[0], polErr[1]}});
            }
            if (ke80) for (const auto& r : keRanges){ // the ZDC KE plot only with the KE cut: the cut alone needs no ZDC hit
                drawPlot(events2, e.branch + ">>hResCE" + binSpec(r), qCut, e.label + " of Protons that hit ZDC, " + sc.label + (ke80 ? ", ZDC KE > 80%" : "") + ", " + kThetaNote + ", CE recon protons vs all recon protons;" + e.label + " (MeV);Count", folder_res + "0_KE_KECUT80_" + e.dir + r.tag + ".pdf", "", [&](TH1* h){
                    h->SetLineColor(kBlue);
                    TH1* hf = HistMT(eventsAllEv, P3(e.branch + ">>hResAll" + binSpec(r)), P3(qCut)); hf->Draw("same");
                    hf->SetLineColor(kRed);
                    h->SetMinimum(0.5); h->SetMaximum(2.0 * std::max(h->GetMaximum(), hf->GetMaximum())); canvas->SetLogy(1);
                    TLegend* leg = new TLegend(0.745, 0.44, 0.99, 0.52); leg->SetTextSize(0.015);
                    leg->AddEntry(h, Form("CE recon protons that hit ZDC%s, %lld", ke80 ? " above the 80% KE cut" : "", nCEq), "l");
                    leg->AddEntry(hf, Form("all recon protons that hit ZDC%s, %lld", ke80 ? " above the 80% KE cut" : "", nAllq), "l");
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
                fprintf(rf, "%s\t%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", r.cut.c_str(), r.sample.c_str(), r.nCEin, r.nAllIn, r.nCEtotal, r.nCEreconBefore,
                        r.fom[0], r.fom[1], r.pol[0], r.pol[1], r.polErr[0], r.polErr[1]);
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
                    if (f.size() != 12) continue;
                    RankRow r{f[0], f[1], std::stod(f[2]), std::stod(f[3]), std::stod(f[4]), std::stod(f[5]), {std::stod(f[6]), std::stod(f[7])}, {std::stod(f[8]), std::stod(f[9])}, {std::stod(f[10]), std::stod(f[11])}};
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
                {"P_phiRe", [](const RankRow& r){return r.pol[0];}, false}, {"P_phiReD", [](const RankRow& r){return r.pol[1];}, false}};
            const char* files[] = {"rank_CE_recon_over_recon_proton_in_cut.csv", "rank_CE_recon_over_CE_recon_before_cut.csv", "rank_nonCE_recon_over_recon_proton_in_cut.csv",
                                   "rank_FOM_theta.csv", "rank_FOM_thetaBlock.csv", "rank_P_phiRe.csv", "rank_P_phiReD.csv"};
            for (size_t c = 0; c < 7; c++){
                std::vector<RankRow> sorted = rankRows;
                std::stable_sort(sorted.begin(), sorted.end(), [&](const RankRow& a, const RankRow& b){return cols[c].lowBest ? cols[c].v(a) < cols[c].v(b) : cols[c].v(a) > cols[c].v(b);});
                FILE* rk = fopen((folder_KEcuts + files[c]).c_str(), "w");
                fprintf(rk, "rank,cut,sample");
                for (const auto& col : cols) fprintf(rk, ",%s", col.name);
                fprintf(rk, ",P_phiRe_err,P_phiReD_err,CE_recon_proton_in_cut,recon_proton_in_cut,CE_protons_total,CE_recon_before_cut,nonCE_recon_proton_in_cut\n");
                for (size_t i = 0; i < sorted.size(); i++){
                    const RankRow& r = sorted[i];
                    fprintf(rk, "%zu,%s,%s", i + 1, r.cut.c_str(), r.sample.c_str());
                    for (const auto& col : cols) fprintf(rk, ",%.6g", col.v(r));
                    fprintf(rk, ",%.5f,%.5f,%.0f,%.0f,%.0f,%.0f,%.0f\n", r.polErr[0], r.polErr[1], r.nCEin, r.nAllIn, r.nCEtotal, r.nCEreconBefore, r.nAllIn - r.nCEin);
                }
                fclose(rk);
            }
        }

        } // part 2

        // phi cos fits in every cut folder, like 3_ReconCalcVals/2_phiRe_cosfit.pdf, CE recon protons and all recon protons of the folder on the same graph, at their recon phi and recon block phi:
        // N(phi) = N0 (1 + P <A> cos(phi)), N0 = entries / 50 bins fixed, <A> = mean model A_y at the recon (block) theta of those protons fixed, only P fitted
        if (RunPart(4)) {
            struct PhiSet {std::string dir, expr; TTree* t; std::string base;}; // base: the sample the cut applies to (hit ZDC and reconstructed, or any fitted track for the cuts)
            std::vector<PhiSet> phiSets = {{"1_hitZDC", "", events2, zdcBase}, {"2_KECUT80", " && fDEnergy > " + std::to_string(initialKE * 0.80), events2, zdcBase}};
            for (const auto& sc : sampleCuts) phiSets.push_back({sc.dir, "", cutTree(sc), rawCE});
            for (const auto& ps : phiSets){
                if (const int pi = int(&ps - &phiSets[0]); !RunCut(pi < 2 ? 0 : pi - 1)) continue; // P3_CUT: 0 the first two sets, i the sampleCuts[i - 1] set
                gSystem->mkdir((folder_KEcuts + ps.dir + "/2_phi_cosfit").c_str(), true);
                const bool isZdc = ps.t == events2; // the ZDC folders hold the protons that hit the ZDC, the cut folders the fitted protons that pass the cut
                struct Samp {TTree* t; std::string q, name; int color;}; // CE recon protons (blue) and all recon protons (red) of the folder, drawn on the same graph
                const Samp samps[2] = {{ps.t, ps.base + ps.expr, "CE recon protons", kBlue}, {isZdc ? eventsAllEv : ps.t, (isZdc ? zdcBase : std::string("1")) + ps.expr, "all recon protons", kRed}};
                const std::string cutTxt = [&]() -> std::string { // every cut the folder applies, for the title
                    if (isZdc) return "hit ZDC, " + kThetaNote + (ps.expr.empty() ? ", no KE cut" : ", ZDC KE > 80%");
                    for (const auto& c : sampleCuts) if (c.dir == ps.dir) return c.label + ", no ZDC requirement";
                    return ps.dir;
                }();
                const std::string where = isZdc ? " that hit ZDC" : " in the cut", dirP = folder_KEcuts + ps.dir + "/2_phi_cosfit/";
                const double apScalePhi = ReadSettings().at("APScale"), pTrue = ReadSettings().at("PPValue");
                for (int j = 0; j < 2; j++){ // j 0: phiRe at thetaRe, 1: phiReD at thetaReD
                    const int kT = j ? 11 : 9, kPh = j ? 15 : 14;
                    const std::string rb = j ? "recon block" : "recon", stem = j ? "phiReD" : "phiRe";
                    std::vector<double> th[2], ph[2], aEv[2]; // theta, phi and the model A of every proton of the sample
                    for (int s = 0; s < 2; s++){
                        std::string qs = samps[s].q + Form(" && rawOk(fEvent, TrackID, %d) && rawOk(fEvent, TrackID, %d)", kT, kPh);
                        if (j) qs += Form(" && rawVal(fEvent, TrackID, %d) >= %g && rawVal(fEvent, TrackID, %d) <= %g", kPh, kBlockPhiLo, kPh, kBlockPhiHi); // block phi: no spike at 0
                        const std::string colT = Form("rawVal(fEvent, TrackID, %d)", kT), colPh = Form("rawVal(fEvent, TrackID, %d)", kPh);
                        auto c = Columns(samps[s].t, {P3(colT), P3(colPh)}, P3(qs));
                        th[s] = std::move(c[0]); ph[s] = std::move(c[1]);
                        aEv[s].resize(th[s].size());
                        for (size_t i = 0; i < th[s].size(); i++) aEv[s][i] = APFromT(kBeamMomentum, TFromAngle(kBeamMomentum, (std::floor(th[s][i] / 0.05) + 0.5) * 0.05 * TMath::Pi() / 180), apScalePhi);
                    }
                    // 1: phi distribution and cos fit, N(phi) = N0 (1 + P <A> cos(phi)), N0 = entries / 50 bins fixed, <A> = mean model A of the sample fixed, only P fitted
                    {
                        canvas->Clear(); canvas->Divide(2, 1); // CE recon protons on the left, all recon protons on the right, each on its own linear y axis
                        TH1D* hP[2] = {nullptr, nullptr};
                        double pol[2] = {0, 0}, polErr[2] = {0, 0}, aMean[2] = {0, 0};
                        for (int s = 0; s < 2; s++){
                            delete gDirectory->Get(Form("hPhiOv%d", s));
                            hP[s] = new TH1D(Form("hPhiOv%d", s), "", 50, 0, 360);
                            for (double x : ph[s]) hP[s]->Fill(x);
                            for (double a : aEv[s]) aMean[s] += a / std::max<size_t>(1, aEv[s].size());
                            hP[s]->SetLineColor(samps[s].color); hP[s]->SetLineWidth(2);
                            if (hP[s]->GetEntries() >= 10){
                                TF1 fit(Form("fPhiOv%d", s), "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                                fit.SetParNames("N0", "P", "A");
                                fit.FixParameter(0, hP[s]->GetEntries() / 50); fit.FixParameter(2, aMean[s]); fit.SetParameter(1, 0.5); fit.SetParLimits(1, -1, 1);
                                hP[s]->Fit(&fit, "Q0"); // the fit stays in the histogram for stats.csv and fit_values.csv, drawn below
                                pol[s] = fit.GetParameter(1); polErr[s] = fit.GetParError(1);
                            }
                            canvas->cd(s + 1);
                            gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.17); gPad->SetBottomMargin(0.13);
                            hP[s]->SetTitle(Form("%s%s, %s, %s phi distribution;Angle (degrees);Count", samps[s].name.c_str(), where.c_str(), cutTxt.c_str(), rb.c_str()));
                            hP[s]->SetMinimum(0); hP[s]->SetMaximum(1.1 * hP[s]->GetMaximum());
                            hP[s]->SetStats(0);
                            hP[s]->Draw("hist");
                            if (hP[s]->GetEntries() >= 10){ // the fit curve, a function of its own: a clone of the one in the histogram crashed when the canvas was cleared
                                TF1* fd = new TF1(Form("fPhiDraw%d", s), "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                                fd->SetParameters(hP[s]->GetEntries() / 50, pol[s], aMean[s]); fd->SetLineColor(kBlack); fd->SetLineWidth(3); fd->Draw("same");
                            }
                            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
                            tx.DrawLatex(0.18, 0.88, Form("N(#phi) = N_{0} (1 + P #upoint #LTA#GT cos#phi), %.0f entries", hP[s]->GetEntries()));
                            tx.DrawLatex(0.18, 0.835, Form("#LTA#GT = %.4f (fixed, %s theta),  P = %.3g #pm %.2g", aMean[s], rb.c_str(), pol[s], polErr[s]));
                            if (j) tx.DrawLatex(0.18, 0.79, BlockPhiNote());
                        }
                        const std::string path = dirP + std::to_string(j * 3 + 1) + "_" + stem + "_cosfit.pdf";
                        canvas->cd(); canvas->Modified(); canvas->Update();
                        canvas->SaveAs(path.c_str());
                        for (int s = 0; s < 2; s++){
                            hP[s]->SetName(StatsTitle(hP[s], samps[s].name + where).c_str());
                            WritePlotNumbers(hP[s], path, fitFilesStarted); // stats.csv and fit_values.csv (N0, P, A) next to the pdf
                            delete hP[s];
                        }
                        canvas->Clear();
                    }
                    // 2, 3: P of the same fit for the protons from a start theta up to 5 deg (recon theta, or recon block theta for j 1), start in 0.05 deg steps, and its percent error to the true P (PPValue, dotted line)
                    struct Scan {std::vector<double> x, y, ey, kept;}; // kept: % of all CE recon protons (CE_recon_before_cut), of the sample for all recon protons, that are in the sample with theta from the start value on
                    Scan scan[2];
                    for (int s = 0; s < 2; s++){
                        const size_t n = th[s].size();
                        for (int b = 0; b < 100; b++){
                            const double start = b * 0.05;
                            TH1D hS("hPstart", "", 50, 0, 360);
                            double aSum = 0;
                            for (size_t i = 0; i < n; i++){
                                if (th[s][i] < start || th[s][i] >= 5.0) continue;
                                hS.Fill(ph[s][i]);
                                aSum += aEv[s][i];
                            }
                            if (hS.GetEntries() < 50) break;
                            TF1 fS("fPstart", "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                            fS.FixParameter(0, hS.GetEntries() / 50); fS.FixParameter(2, aSum / hS.GetEntries()); fS.SetParameter(1, 0.5); fS.SetParLimits(1, -1, 1);
                            hS.Fit(&fS, "QN");
                            scan[s].x.push_back(start); scan[s].y.push_back(fS.GetParameter(1)); scan[s].ey.push_back(fS.GetParError(1));
                            scan[s].kept.push_back(100.0 * hS.GetEntries() / std::max(1.0, s ? double(n) : nCEreconTotal));
                        }
                    }
                    for (const bool pct : {false, true}) if (scan[0].x.size() > 1 || scan[1].x.size() > 1){
                        canvas->Clear(); canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
                        const double yRef = pct ? 0.0 : pTrue;
                        double ymin = yRef, ymax = yRef, xMax = 0.05;
                        for (int s = 0; s < 2; s++) for (size_t i = 0; i < scan[s].x.size(); i++){
                            const double v = pct ? 100 * (scan[s].y[i] - pTrue) / pTrue : scan[s].y[i], e = pct ? 100 * scan[s].ey[i] / pTrue : scan[s].ey[i];
                            ymin = std::min(ymin, v - e); ymax = std::max(ymax, v + e); xMax = std::max(xMax, scan[s].x[i] + 0.05);
                        }
                        const double pad = 0.05 * (ymax - ymin);
                        TH2D axis("hPstartAxis", Form("Polarization fit P of %s and %s%s, %s, %s phi, %s theta from the start value to 5 deg%s;Start theta (deg);%s", samps[0].name.c_str(), samps[1].name.c_str(), where.c_str(), ps.dir.c_str(), rb.c_str(), rb.c_str(), j ? ", phi 0.5 to 359.5 deg" : "",
                                  pct ? "100 (P_{fit} - P_{true}) / P_{true} (%)" : "P"), 100, 0, xMax, 100, ymin - pad, ymax + pad);
                        axis.SetStats(0);
                        axis.Draw("axis");
                        TLine ref(0, yRef, xMax, yRef); ref.SetLineStyle(3); ref.SetLineWidth(2); ref.SetLineColor(kGray + 2); ref.Draw(); // true P, from the run's PPValue
                        TLegend leg(0.705, 0.80, 0.99, 0.90); leg.SetTextSize(0.014);
                        TGraphErrors* gs[2] = {nullptr, nullptr};
                        const std::string path = dirP + std::to_string(j * 3 + 2 + pct) + "_" + stem + (pct ? "_P_vs_start_theta_pct_error.pdf" : "_P_vs_start_theta.pdf");
                        for (int s = 1; s >= 0; s--){
                            gs[s] = new TGraphErrors();
                            for (size_t i = 0; i < scan[s].x.size(); i++){
                                gs[s]->SetPoint(i, scan[s].x[i], pct ? 100 * (scan[s].y[i] - pTrue) / pTrue : scan[s].y[i]);
                                gs[s]->SetPointError(i, 0, pct ? 100 * scan[s].ey[i] / pTrue : scan[s].ey[i]);
                            }
                            gs[s]->SetLineColor(samps[s].color); gs[s]->SetMarkerColor(samps[s].color); gs[s]->SetMarkerStyle(20); gs[s]->SetMarkerSize(0.6);
                            if (gs[s]->GetN()) gs[s]->Draw("PZ same");
                            leg.AddEntry(gs[s], samps[s].name.c_str(), "lp");
                            for (size_t i = 0; i < scan[s].x.size(); i++) // the points drawn: start theta, value, error, % kept
                                WriteStatsRow(path, "plot,sample,start_theta_deg,value,error,pct_kept", Form("%s,%s,%.4g,%.6g,%.4g,%.4g", gSystem->BaseName(path.c_str()), samps[s].name.c_str(), scan[s].x[i], gs[s]->GetY()[i], gs[s]->GetEY()[i], scan[s].kept[i]), fitFilesStarted, "P_vs_start_theta.csv");
                        }
                        leg.Draw();
                        canvas->SaveAs(path.c_str());
                        for (auto* g : gs) delete g;
                    }
                    canvas->SetRightMargin(0.30); canvas->SetTopMargin(0.1);
                }
            }
        }

        // recon protons that hit the ZDC, all KE (no cut): CE protons and every recon proton of all events
        if (RunPart(1)) for (const auto& e : cutEnergies) for (const auto& r : keRanges){
            std::string folder_hit = folder_KEcuts + "1_hitZDC/";
            gSystem->mkdir(folder_hit.c_str(), true);
            drawPlot(events2, e.branch + ">>hHitCE" + binSpec(r), zdcBase, e.label + " of Protons that hit the ZDC, CE recon protons that hit ZDC vs all events recon, " + kThetaNote + ", no KE cut;" + e.label + " (MeV);Count", folder_hit + "1_real_vs_recon_" + e.dir + r.tag + ".pdf", "", [&](TH1* h){
                h->SetLineColor(kBlue);
                TH1* hr = HistMT(eventsAllEv, P3(e.branch + ">>hHitRecon" + binSpec(r)), P3(zdcBase)); hr->Draw("same");
                hr->SetLineColor(kRed);
                h->SetMinimum(0.5); h->SetMaximum(2.0 * std::max(h->GetMaximum(), hr->GetMaximum())); canvas->SetLogy(1);
                TLegend* leg = new TLegend(0.745, 0.44, 0.99, 0.52); leg->SetTextSize(0.015);
                leg->AddEntry(h, Form("CE recon protons that hit ZDC, %lld", (long long)h->GetEntries()), "l");
                leg->AddEntry(hr, Form("all recon protons that hit ZDC, %lld", (long long)hr->GetEntries()), "l");
                leg->Draw();
            });
            canvas->SetLogy(0);
        }

        // 1_additionalinfo in every KE cut folder: the cut candidates (6_cuts/3_ResCuts to 8_nHits) of the reconstructed protons that hit the ZDC and pass that KE cut
        // (1_hitZDC: no KE cut), reconstructed CE protons vs all reconstructed protons of all events. Nothing is cut on these variables.
        struct InfoVar {std::string name, folder, label; int k; double lo, hi; int nBins;}; // lo >= hi: range from the 1 to 99% quantiles of all fitted protons
        const std::vector<InfoVar> kInfoVars = {
            {"resDRms", "3_ResCuts", "residual, block position (mm)", 1, 0, 0, 100},
            {"nProtonsEv", "4_nProtons", "protons with a track in the event", 6, -0.5, 19.5, 20},
            {"nHits", "8_nHits", "tracker hits", 10, -0.5, 19.5, 20},
            {"thetaBlock", "7_theta", "recon block theta (degrees)", 11, 0, 2.3, 100}};
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
            std::vector<std::pair<std::string, double>> infoCuts = {{"1_hitZDC", -1.0}}; // (folder, KE cut in MeV, -1 none)
            for (int pct : {80}) infoCuts.push_back({"2_KECUT" + std::to_string(pct), initialKE * pct / 100.0});
            for (const auto& sc : sampleCuts) infoCuts.push_back({sc.dir, -1.0}); // no KE cut: what each cut does to every variable
            auto lapStart = std::chrono::steady_clock::now();
            const auto lap = [&](const std::string& what){ // seconds since the last lap
                const auto now = std::chrono::steady_clock::now();
                printf(">>> 1_additionalinfo %s took %.1f s\n", what.c_str(), std::chrono::duration<double>(now - lapStart).count()); fflush(stdout);
                lapStart = now;
            };
            for (const auto& ic : infoCuts){
                if (const int ii = int(&ic - &infoCuts[0]); !RunCut(ii < 2 ? 0 : ii - 1)) continue; // P3_CUT: 0 the first two folders, i the sampleCuts[i - 1] folder
                lapStart = std::chrono::steady_clock::now();
                const SampleCut* sc = nullptr;
                for (const auto& c : sampleCuts) if (c.dir == ic.first) sc = &c;
                const bool nz = sc != nullptr; // the cuts' plots use every fitted proton (AnglesRaw), the others the protons that hit the ZDC
                TTree *tCE = nz ? cutTree(*sc) : events2, *tAll = nz ? tCE : eventsAllEv;
                const std::string ke = ic.second < 0 ? "" : " && fDEnergy > " + std::to_string(ic.second);
                const std::string bCE = nz ? rawCE : zdcBase, bAll = nz ? "1" : zdcBase, hitTxt = nz ? "" : " that hit ZDC";
                // every histogram of the folder is booked first and filled in one event loop per tree (one for both samples of a cut, they share the cut file), not one loop per plot
                std::vector<std::pair<std::string, std::string>> reqCE, reqAll; // (Draw string, cut) of the CE recon protons and of all recon protons
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
                }
                lap(ic.first + ": quantiles and booking the requests");
                std::vector<TH1*> hCEs, hAlls;
                if (tCE == tAll){
                    std::vector<std::pair<std::string, std::string>> both = reqCE;
                    both.insert(both.end(), reqAll.begin(), reqAll.end());
                    const std::vector<TH1*> hs = HistsMT(tCE, both);
                    hCEs.assign(hs.begin(), hs.begin() + reqCE.size()); hAlls.assign(hs.begin() + reqCE.size(), hs.end());
                } else {hCEs = HistsMT(tCE, reqCE); hAlls = HistsMT(tAll, reqAll);}
                lap(ic.first + ": jit and event loop of " + std::to_string(reqCE.size() + reqAll.size()) + " histograms");
                for (size_t iv = 0; iv < kInfoVars.size(); iv++){
                    const auto& v = kInfoVars[iv];
                    if (slot[iv] < 0) continue;
                    TH1 *hCE = hCEs[slot[iv]], *hRec = hAlls[slot[iv]];
                    if (!hCE || !hRec) continue;
                    canvas->Clear(); canvas->SetLogy(1);
                    hCE->SetLineColor(kBlue); hRec->SetLineColor(kRed); hCE->SetLineWidth(2); hRec->SetLineWidth(2); hCE->SetStats(0); hRec->SetStats(0);
                    hRec->SetTitle((v.label + ", CE recon protons vs all recon protons" + (nz ? std::string(", ") : std::string(", hit ZDC, " + kThetaNote + ", ")) + (ic.second < 0 ? std::string(sc ? "no KE cut, " + sc->label + ", no ZDC requirement" : "no KE cut") : "ZDC KE > " + std::to_string(int(ic.second)) + " MeV") + ";" + v.label + ";Count").c_str());
                    if (v.name == "nHits") for (int bin = 1; bin <= hRec->GetNbinsX(); bin++) hRec->GetXaxis()->SetBinLabel(bin, bin % 4 == 1 && bin > 1 ? Form("%d", bin - 1) : ""); // labels at 4, 8, 12, 16 hits only
                    hRec->SetMinimum(0.5); hRec->SetMaximum(2.0 * std::max(hRec->GetMaximum(), hCE->GetMaximum()));
                    hRec->Draw("hist"); hCE->Draw("hist same");
                    TLegend* leg = new TLegend(0.745, 0.80, 0.99, 0.90); leg->SetTextSize(0.014);
                    leg->AddEntry(hCE, Form("CE recon protons (%.0f)", hCE->GetEntries()), "l");
                    leg->AddEntry(hRec, Form("all recon protons (%.0f)", hRec->GetEntries()), "l");
                    leg->Draw();
                    // flat files named <folder>_<variable> in 1_additionalinfo of the hit ZDC, KE cut and 2 mm cut folders
                    const std::string idir = folder_KEcuts + ic.first + (sc && v.folder == sc->onFolder ? "/" : "/1_additionalinfo/"); // the variable a cut is on sits beside 1_additionalinfo
                    gSystem->mkdir(idir.c_str(), true);
                    canvas->SaveAs((idir + v.folder + "_" + v.name + ".pdf").c_str());
                    canvas->SetLogy(0);
                    delete hCE; delete hRec;
                }
                lap(ic.first + ": drawing the plots");
            }
        }

    }
    BoxAndCropPlots();
}
