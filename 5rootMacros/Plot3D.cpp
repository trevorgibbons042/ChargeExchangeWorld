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
#include <functional>
#include <memory>
#include <set>
#include <TSystem.h>
#include <sstream>
#include <TStyle.h>

#include "FitCore.h"

// (fEvent, TrackID) of every proton with a reconstructed track (cut_angles_all.root, Angles). The position plots only show these protons,
// the ones that can actually be drawn on the ZDC. isReco is called from the TTree::Draw cuts.
std::set<std::pair<int, int>> gRecoTracks;
bool isReco(int fEvent, int trackID){return gRecoTracks.count({fEvent, trackID}) > 0;}

void Plot3D(){
    std::string folder = "plots";

    // step folders on top, particle set folders inside: 1_CEprotons, 2_allprotons (CE events), 3_allprotons_allevents
    std::string folder_pos = folder + "/1_position/1_true/";
    std::string folder_posD = folder + "/1_position/2_block/";

    std::string folder_smear = folder + "/2_energy/";
    std::string folder_KEcuts = folder + "/3_KEcuts/";


    gSystem->mkdir(folder.c_str(), true);
    gSystem->mkdir(folder_pos.c_str(), true);
    gSystem->mkdir(folder_posD.c_str(), true);
    gSystem->mkdir(folder_smear.c_str(), true);
    gSystem->mkdir(folder_KEcuts.c_str(), true);

    double Graph = 1;

    {
        TFile* fReco = TFile::Open("cut_angles_all.root", "READ");
        TTree* reco = fReco ? (TTree*)fReco->Get("Angles") : nullptr;
        if (!reco){Error("Plot3D", "missing Angles tree in cut_angles_all.root (run ReconAngleAll.cpp first)"); return;}
        int ev, tid;
        reco->SetBranchAddress("fEvent", &ev); reco->SetBranchAddress("TrackID", &tid);
        for (Long64_t i = 0; i < reco->GetEntries(); i++){reco->GetEntry(i); gRecoTracks.insert({ev, tid});}
        fReco->Close();
    }

    TFile* File = TFile::Open("cut_CE.root", "READ");
    TFile* File2 = TFile::Open("cut_ParentID.root", "READ");
    TFile* FileFake = TFile::Open("cut_FakeID.root", "READ");

    TTree *events = (TTree*)File->Get("ParticleInfo");
    TTree *ChargeExchange = (TTree*)File->Get("ChargeExchange_values");

    TTree *events2 = (TTree*)File2->Get("ParticleInfo");
    TTree *ChargeExchange2 = (TTree*)File2->Get("ChargeExchange_values");

    TTree *eventsFake = (TTree*)FileFake->Get("ParticleInfo");
    TTree *ChargeExchangeFake = (TTree*)FileFake->Get("ChargeExchange_values");

    // every event, not only the CE events kept in cut_CE.root: third "all protons in all events" set next to CE protons and all protons in CE events
    TFile* FileAllEv = TFile::Open("cut_AllEvProtons.root", "READ"); // ZDC protons of every event, written by CutEvents.cpp (every plot below cuts on them anyway)
    TTree* eventsAllEv = FileAllEv ? (TTree*)FileAllEv->Get("ParticleInfo") : nullptr;
    if (!eventsAllEv){Error("Plot3D", "missing cut_AllEvProtons.root (run CutEvents.cpp first)"); return;}
    const std::string allEvDir = "3_allprotons_allevents/";


    if (Graph == 1){

        TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
        gPad->SetRightMargin(0.30);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); // room for the rotated Y title and bottom tick labels
        gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X"); // stats/fit box to the right, fewer x labels

        // Every plot goes through this: draw expr from t with cut, set the title (title;xaxis;yaxis), let tweak adjust the
        // histogram (range, fit, ...), save to path. opt is the Draw option, e.g. "colz".
        std::set<std::string> fitFilesStarted; // stats and fit files already truncated in this run
        auto drawPlot = [&](TTree* t, const std::string& expr, const std::string& cut, const std::string& title,
                            const std::string& path, const char* opt = "", const std::function<void(TH1*)>& tweak = nullptr){
            gSystem->mkdir(gSystem->GetDirName(path.c_str()), true);
            canvas->Clear();
            t->Draw(expr.c_str(), cut.c_str(), opt);
            TH1* h = t->GetHistogram();
            h->SetTitle(PlotTitle(path, title).c_str());
            if (tweak) tweak(h);
            WritePlotNumbers(h, path, fitFilesStarted); // stats.csv, plus fit_values.csv if h has a fit, next to the plot
            canvas->Update();
            canvas->SaveAs(path.c_str());
        };
        // x:y position maps share a count scale of 0 to 100
        auto drawXY = [&](TTree* t, const std::string& expr, const std::string& cut, const std::string& title, const std::string& path){
            drawPlot(t, expr, cut, title, path, "colz", [](TH1* h){h->SetMinimum(0); h->SetMaximum(100);});
        };

        const std::string zdcProtons = "boundaryType == 2 && pdg == 2212";
        // position plots: only protons with a reconstructed track
        const std::string zdcReco = zdcProtons + " && isReco(fEvent, TrackID)";
        const std::string xySlice = zdcReco + " && fabs(fY) < 2";
        const std::string xySliceD = zdcReco + " && fabs(fDy) < 2";

        // ZDC position, true (fX, fY)
        drawXY(events, "fY:fX>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at ZDC, CE events all protons, reco tracks;X Position (mm);Y Position (mm)", folder_pos + "2_allprotons/1_1_xy_protons.pdf");
        drawXY(events, "fY:fX>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at ZDC, CE events all protons, reco tracks;X Position (mm);Y Position (mm)", folder_pos + "2_allprotons/1_2_xy_protons_zoom.pdf");
        drawXY(eventsAllEv, "fY:fX>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at ZDC, all events all protons, reco tracks;X Position (mm);Y Position (mm)", folder_pos + allEvDir + "1_1_xy_protons.pdf");
        drawXY(eventsAllEv, "fY:fX>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at ZDC, all events all protons, reco tracks;X Position (mm);Y Position (mm)", folder_pos + allEvDir + "1_2_xy_protons_zoom.pdf");
        drawXY(events2, "fY:fX>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at ZDC, CE events CE protons, reco tracks;X Position (mm);Y Position (mm)", folder_pos + "1_CEprotons/1_1_xy_parentid.pdf");
        drawXY(events2, "fY:fX>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at ZDC, CE events CE protons, reco tracks;X Position (mm);Y Position (mm)", folder_pos + "1_CEprotons/1_2_xy_parentid_zoom.pdf");

        gStyle->SetOptFit(1111); // probability, chi2/ndf, fit parameters with errors
        drawPlot(events2, "fX>>hSliceX200(100,-200,200)", xySlice, "X Position at ZDC, CE events CE protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_pos + "1_CEprotons/2_1_x_fit_200.pdf", "", [](TH1* h){FitCore(h);});
        drawPlot(events2, "fX>>hSliceX(100,-50,50)", xySlice, "X Position at ZDC, CE events CE protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_pos + "1_CEprotons/2_2_x_fit_50.pdf", "", [](TH1* h){FitCore(h);});
        drawPlot(events, "fX>>hSliceXAll200(100,-200,200)", xySlice, "X Position at ZDC, CE events all protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_pos + "2_allprotons/2_1_x_fit_200.pdf", "", [](TH1* h){FitCore(h);});
        drawPlot(events, "fX>>hSliceXAll(100,-50,50)", xySlice, "X Position at ZDC, CE events all protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_pos + "2_allprotons/2_2_x_fit_50.pdf", "", [](TH1* h){FitCore(h);});
        drawPlot(eventsAllEv, "fX>>hSliceXAllEv200(100,-200,200)", xySlice, "X Position at ZDC, all events all protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_pos + allEvDir + "2_1_x_fit_200.pdf", "", [](TH1* h){FitCore(h);});
        drawPlot(eventsAllEv, "fX>>hSliceXAllEv(100,-50,50)", xySlice, "X Position at ZDC, all events all protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_pos + allEvDir + "2_2_x_fit_50.pdf", "", [](TH1* h){FitCore(h);});


        // ZDC position, block (fDx, fDy)
        drawXY(events, "fDy:fDx>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at ZDC BLOCK, CE events all protons, reco tracks;X Position (mm);Y Position (mm)", folder_posD + "2_allprotons/1_1_xy_protons.pdf");
        drawXY(events, "fDy:fDx>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at ZDC BLOCK, CE events all protons, reco tracks;X Position (mm);Y Position (mm)", folder_posD + "2_allprotons/1_2_xy_protons_zoom.pdf");
        drawXY(eventsAllEv, "fDy:fDx>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at ZDC BLOCK, all events all protons, reco tracks;X Position (mm);Y Position (mm)", folder_posD + allEvDir + "1_1_xy_protons.pdf");
        drawXY(eventsAllEv, "fDy:fDx>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at ZDC BLOCK, all events all protons, reco tracks;X Position (mm);Y Position (mm)", folder_posD + allEvDir + "1_2_xy_protons_zoom.pdf");
        drawXY(events2, "fDy:fDx>>h2(100,-200,200,100,-200,200)", zdcReco, "Position at ZDC BLOCK, CE events CE protons, reco tracks;X Position (mm);Y Position (mm)", folder_posD + "1_CEprotons/1_1_xy_parentid.pdf");
        drawXY(events2, "fDy:fDx>>h2(100,-50,50,100,-50,50)", zdcReco, "Position at ZDC BLOCK, CE events CE protons, reco tracks;X Position (mm);Y Position (mm)", folder_posD + "1_CEprotons/1_2_xy_parentid_zoom.pdf");

        drawPlot(events2, "fDx>>hSliceDX200(100,-200,200)", xySliceD, "X Position at ZDC BLOCK, CE events CE protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_posD + "1_CEprotons/2_1_x_fit_200.pdf", "", [](TH1* h){FitBlocks(h);}); //one bin per ZDC block
        drawPlot(events2, "fDx>>hSliceDX(100,-127.27,127.27)", xySliceD, "X Position at ZDC BLOCK, CE events CE protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_posD + "1_CEprotons/2_2_x_fit_50.pdf", "", [](TH1* h){FitBlocks(h);}); //7 central ZDC blocks, pitch 400/11 mm
        drawPlot(events, "fDx>>hSliceDXAll200(100,-200,200)", xySliceD, "X Position at ZDC BLOCK, CE events all protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_posD + "2_allprotons/2_1_x_fit_200.pdf", "", [](TH1* h){FitBlocks(h);});
        drawPlot(events, "fDx>>hSliceDXAll(100,-127.27,127.27)", xySliceD, "X Position at ZDC BLOCK, CE events all protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_posD + "2_allprotons/2_2_x_fit_50.pdf", "", [](TH1* h){FitBlocks(h);});
        drawPlot(eventsAllEv, "fDx>>hSliceDXAllEv200(100,-200,200)", xySliceD, "X Position at ZDC BLOCK, all events all protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_posD + allEvDir + "2_1_x_fit_200.pdf", "", [](TH1* h){FitBlocks(h);});
        drawPlot(eventsAllEv, "fDx>>hSliceDXAllEv(100,-127.27,127.27)", xySliceD, "X Position at ZDC BLOCK, all events all protons, reco tracks, |Y| < 2 mm;X Position (mm);Count", folder_posD + allEvDir + "2_2_x_fit_50.pdf", "", [](TH1* h){FitBlocks(h);});

        // Every KE x axis comes in two ranges, each saved as its own file: 0 to 30k MeV and 15k to 30k MeV
        struct KERange {double lo, hi; std::string tag;};
        const std::vector<KERange> keRanges = {{0, 30000, "_0to30k"}, {10000, 30000, "_10to30k"}};
        auto binSpec = [](const KERange& r){return std::string(Form("(100,%g,%g)", r.lo, r.hi));};

        // Energy plots: 2_energy/<set>/ with set = 1_CEprotons, 2_allprotons (CE events), 3_allprotons_allevents
        // True KE (fKEnergy) and smeared KE (fDEnergy) at the ZDC. CE protons and all protons in CE events share the y scale of the taller one.
        struct KESet {TTree* tree; std::string dir, label;};
        const std::vector<KESet> keSets = {{events2, "1_CEprotons/", "CE events CE protons"}, {events, "2_allprotons/", "CE events all protons"}, {eventsAllEv, allEvDir, "all events all protons"}};
        struct KEQty {std::string branch, file, axis, title;};
        std::vector<KEQty> keQtys = {{"fKEnergy", "1_true_ke", "Kinetic Energy (MeV)", "True Kinetic Energy at ZDC"}};
        if (events->GetBranch("fDEnergy")) keQtys.push_back({"fDEnergy", "2_smeared_ke", "Smeared Energy fDEnergy (MeV)", "Smeared Energy at ZDC"}); // files from before fDEnergy was added don't have it
        for (const auto& q : keQtys){
            for (const auto& r : keRanges){
                double ymax = 0;
                for (TTree* t : {events, events2}){
                    t->Draw((q.branch + ">>hKEScale" + binSpec(r)).c_str(), zdcProtons.c_str(), "goff");
                    TH1* hs = (TH1*)gDirectory->Get("hKEScale");
                    ymax = std::max(ymax, 1.1 * hs->GetBinContent(hs->GetMaximumBin()));
                    delete hs;
                }
                for (size_t i = 0; i < keSets.size(); i++){
                    const auto& s = keSets[i];
                    drawPlot(s.tree, q.branch + ">>hKE" + binSpec(r), zdcProtons, q.title + ", " + s.label + ";" + q.axis + ";Count", folder_smear + s.dir + q.file + r.tag + ".pdf", "", [&](TH1* h){if (i < 2) h->SetMaximum(ymax);});
                }
            }
        }

        // ZDC energy smearing (fDEnergy = Gaussian smear of fKEnergy, /zdc/resStochastic and /zdc/resConstant in the sim)
        if (events->GetBranch("fDEnergy")){ // files from before fDEnergy was added don't have it
            gSystem->mkdir(folder_smear.c_str(), true);

            // fKEnergy = true KE, fDEnergy = smeared KE, dE = fDEnergy - fKEnergy

            drawPlot(events, "fDEnergy-fKEnergy>>hSmearDiffP(100,-10000,10000)", zdcProtons, "Energy Difference at ZDC, CE events all protons;fDEnergy - fKEnergy (MeV);Count", folder_smear + "2_allprotons/4_1_diff_protons.pdf", "", [](TH1* h){FitCore(h);});
            drawPlot(events2, "fDEnergy-fKEnergy>>hSmearDiffCE(100,-10000,10000)", zdcProtons, "Energy Difference at ZDC, CE events CE protons;fDEnergy - fKEnergy (MeV);Count", folder_smear + "1_CEprotons/4_1_diff_parentid.pdf", "", [](TH1* h){FitCore(h);});
            drawPlot(eventsAllEv, "fDEnergy-fKEnergy>>hSmearDiffAllEv(100,-10000,10000)", zdcProtons, "Energy Difference at ZDC, all events all protons;fDEnergy - fKEnergy (MeV);Count", folder_smear + allEvDir + "4_1_diff_protons.pdf", "", [](TH1* h){FitCore(h);});

            // percent error = 100 * (fDEnergy - fKEnergy) / fKEnergy
            const std::string pct = "100*(fDEnergy-fKEnergy)/fKEnergy";
            drawPlot(events, pct + ">>hSmearPctP(100,-100,100)", zdcProtons, "Percent Error at ZDC, CE events all protons;Percent Error (%);Count", folder_smear + "2_allprotons/4_2_pct_protons.pdf", "", [](TH1* h){FitCore(h);});
            drawPlot(events2, pct + ">>hSmearPctCE(100,-100,100)", zdcProtons, "Percent Error at ZDC, CE events CE protons;Percent Error (%);Count", folder_smear + "1_CEprotons/4_2_pct_parentid.pdf", "", [](TH1* h){FitCore(h);});
            drawPlot(eventsAllEv, pct + ">>hSmearPctAllEv(100,-100,100)", zdcProtons, "Percent Error at ZDC, all events all protons;Percent Error (%);Count", folder_smear + allEvDir + "4_2_pct_protons.pdf", "", [](TH1* h){FitCore(h);});

        }
        else{Warning("Plot3D", "no fDEnergy branch in cut_CE.root, skipping smearing plots");}

        const double initialKE = 25000; // MeV, beam energy (25 GeV)
        // KE cuts on the true KE (fKEnergy) and on the ZDC-smeared energy (fDEnergy), each in its own folder
        // ke_cut_<pct>/ as 1_* (true KE) and 2_* (ZDC KE). Files from before fDEnergy was added only get the true KE cuts.
        struct CutEnergy {std::string branch, dir, label;};
        std::vector<CutEnergy> cutEnergies = {{"fKEnergy", "KEtrue", "True KE"}};
        if (events->GetBranch("fDEnergy")) cutEnergies.push_back({"fDEnergy", "KEDetector", "ZDC KE"});
        for (int pct : {50, 60, 70, 80, 90}){
            std::string cut = std::to_string(initialKE * pct / 100.0);
            std::string pctStr = std::to_string(pct) + "% of initial KE";
            for (size_t ei = 0; ei < cutEnergies.size(); ei++){
                const auto& e = cutEnergies[ei];
                const std::string n = std::to_string(ei + 1) + "_"; // 1_x = true KE, 2_x = ZDC KE, in one folder per cut
                std::string folder_cut = folder_KEcuts + "ke_cut_" + std::to_string(pct) + "/" + n + e.branch + "/"; // 1_fKEnergy, 2_fDEnergy
                gSystem->mkdir(folder_cut.c_str(), true);
                const std::string axis = ";" + e.label + " (MeV);Count";
                const double cutVal = initialKE * pct / 100.0;
                const std::string realOutCut = zdcReco + " && " + e.branch + " < " + cut;
                const std::string fakeKeptCut = zdcReco + " && " + e.branch + " > " + cut;
                // counts over all KE (not only the plotted range): CE protons left out of the cut out of all CE protons, fake protons kept in out of all fake protons
                const long long nRealOut = events2->GetEntries(realOutCut.c_str()), nRealTotal = events2->GetEntries(zdcReco.c_str());
                const long long nFakeKept = eventsFake->GetEntries(fakeKeptCut.c_str()), nFakeTotal = eventsFake->GetEntries(zdcReco.c_str());
                const double realOutPct = 100.0 * nRealOut / std::max<long long>(1, nRealTotal);
                const double fakeKeptPct = 100.0 * nFakeKept / std::max<long long>(1, nFakeTotal);
                for (const auto& r : keRanges){
                drawPlot(events2, e.branch + ">>hCutReal" + binSpec(r), realOutCut, e.label + " of Protons Left Out with " + pctStr + " cut, CE events CE protons" + axis, folder_cut + "1_CEprotons/1_real_cut_" + e.dir + r.tag + ".pdf", "", [&](TH1*){
                    TPaveText* pt = new TPaveText(0.13, 0.82, 0.55, 0.88, "NDC"); pt->SetFillColor(0); pt->SetTextAlign(12); pt->SetTextSize(0.022);
                    pt->AddText(Form("left out: %lld / %lld CE protons (reco+ZDC) = %.1f%%", nRealOut, nRealTotal, realOutPct));
                    pt->Draw();
                });
                drawPlot(eventsFake, e.branch + ">>hCutFake" + binSpec(r), zdcReco + " && " + e.branch + " > " + cut, e.label + " of Protons Kept In with " + pctStr + " cut, CE events fake protons" + axis, folder_cut + "2_allprotons/2_fake_cut_" + e.dir + r.tag + ".pdf");
                // real (CE protons) and fake (other protons in CE events) on one plot, with the cut and what it does to each
                drawPlot(events2, e.branch + ">>hCutRealOv" + binSpec(r), zdcReco, e.label + " of Protons, real vs fake, " + pctStr + " cut, CE events" + axis, folder_cut + "1_CEprotons/3_real_vs_fake_" + e.dir + r.tag + ".pdf", "", [&](TH1* h){
                    h->SetLineColor(kBlue);
                    eventsFake->Draw((e.branch + ">>hCutFakeOv" + binSpec(r)).c_str(), zdcReco.c_str(), "same");
                    TH1* hf = static_cast<TH1*>(gDirectory->Get("hCutFakeOv"));
                    hf->SetLineColor(kRed);
                    for (TH1* hh : {h, hf}) hh->GetXaxis()->SetRangeUser(std::max(cutVal, r.lo), r.hi); // only the part above the cut is drawn
                    h->SetMaximum(1.15 * std::max(h->GetMaximum(), hf->GetMaximum()));
                    TLegend* leg = new TLegend(0.13, 0.82, 0.42, 0.88); leg->SetTextSize(0.022);
                    leg->AddEntry(h, Form("real CE protons, %lld / %lld = %.1f%% left out", nRealOut, nRealTotal, realOutPct), "l");
                    leg->AddEntry(hf, Form("fake protons, %lld / %lld = %.1f%% kept in", nFakeKept, nFakeTotal, fakeKeptPct), "l");
                    leg->Draw();
                });
                }
            }
        }

        // Real CE protons / all protons per KE bin over the full KE range, in every ke_cut_<pct> folder next to the cut plots, cut marked.
        // Two samples of "all protons": every event (raw output0_t*.root) and CE events only (cut_CE.root). 1_* true KE, 2_* ZDC KE:
        // 4 = real CE / all protons in all events, 5 = real CE / all protons in CE events. Legend gives the overall and above-the-cut percentages.
        {
            const int nFine = 3000, rebin = 30; // 10 MeV bins hold every cut value on an edge, rebinned to 100 plot bins
            auto fineHist = [&](TTree* t, const std::string& branch, const std::string& name, double lo, double hi){
                t->Draw((branch + ">>" + name + "(" + std::to_string(nFine) + "," + std::to_string(lo) + "," + std::to_string(hi) + ")").c_str(), zdcProtons.c_str(), "goff");
                return static_cast<TH1*>(gDirectory->Get(name.c_str()));
            };
            struct Sample {TTree* tree; std::string dir, label;};
            const std::vector<Sample> samples = {{eventsAllEv, "4_real_over_allevents_", "all events"}, {events, "5_real_over_ceevents_", "CE events"}};
            for (size_t ei = 0; ei < cutEnergies.size(); ei++){
                const auto& e = cutEnergies[ei];
                std::vector<TH1*> hRealR; // per x range, 3000 fine bins rebinned to the 100 plot bins; the 0 to 30k one also gives the counts
                for (const auto& r : keRanges) hRealR.push_back(fineHist(events2, e.branch, "hFineReal" + e.dir + r.tag, r.lo, r.hi));
                TH1* hReal = hRealR[0];
                for (size_t si = 0; si < samples.size(); si++){
                    const auto& s = samples[si];
                    std::vector<TH1*> hSampleR;
                    for (const auto& r : keRanges) hSampleR.push_back(fineHist(s.tree, e.branch, "hFineAll" + e.dir + std::to_string(si) + r.tag, r.lo, r.hi));
                    TH1* hSample = hSampleR[0];
                    const double nReal = hReal->Integral(0, nFine + 1), nSample = hSample->Integral(0, nFine + 1);
                    for (int pct : {50, 60, 70, 80, 90}){
                        const double cutVal = initialKE * pct / 100.0;
                        const int firstAbove = static_cast<int>(cutVal / (30000.0 / nFine)) + 1;
                        const double nRealAbove = hReal->Integral(firstAbove, nFine + 1), nSampleAbove = hSample->Integral(firstAbove, nFine + 1);
                        const std::string folder_all = folder_KEcuts + "ke_cut_" + std::to_string(pct) + "/" + std::to_string(ei + 1) + "_" + e.branch + "/";
                        gSystem->mkdir(folder_all.c_str(), true);
                        for (size_t ri = 0; ri < keRanges.size(); ri++){
                        const auto& r = keRanges[ri];
                        drawPlot(events2, e.branch + ">>hRatioKE" + binSpec(r), zdcProtons,e.label + " of Protons, real CE / all protons in " + s.label + ", " + std::to_string(pct) + "% of initial KE cut;" + e.label + " (MeV);real CE / all protons", folder_all + "2_allprotons/" + s.dir + e.dir + r.tag + ".pdf", "", [&](TH1* h){
                            h->Reset();
                            TH1* hr = static_cast<TH1*>(hRealR[ri]->Clone()); hr->Rebin(rebin);
                            TH1* hs = static_cast<TH1*>(hSampleR[ri]->Clone()); hs->Rebin(rebin);
                            h->Add(hr); h->Divide(hs);
                            delete hr; delete hs;
                            h->GetXaxis()->SetRangeUser(std::max(cutVal, r.lo), r.hi); // only the part above the cut is drawn
                            h->SetMinimum(0); h->SetMaximum(1.15 * h->GetMaximum());
                            TPaveText* pt = new TPaveText(0.13, 0.82, 0.42, 0.88, "NDC"); pt->SetFillColor(0); pt->SetTextAlign(12); pt->SetTextSize(0.022);
                            pt->AddText(Form("all KE: %.0f/%.0f = %.1f%%", nReal, nSample, 100.0 * nReal / std::max(1.0, nSample)));
                            pt->AddText(Form("above cut: %.0f/%.0f = %.1f%%", nRealAbove, nSampleAbove, 100.0 * nRealAbove / std::max(1.0, nSampleAbove)));
                            pt->Draw();
                        });
                        }
                    }
                }
            }
        }

        // Angles at the ZDC in degrees, protons only, computed from the ParticleInfo position at the ZDC boundary.
        // True position = (fX, fY, fZ), block position = (fDx, fDy, fDz).
        // events = all protons in CE events (cut_CE.root), events2 = CE protons, same ParentID (cut_ParentID.root).
        // Zoom = theta < zoomThetaMax; the phi zoom only keeps those small-theta protons.
        const double zoomThetaMax = 0.5; // degrees
        const std::string zdcProtonCut = "boundaryType == 2 && pdg == 2212";
        const std::string zoomLabel = ", theta < " + std::to_string(zoomThetaMax).substr(0, 3) + " deg";
        auto drawAngle = [&](TTree* t, const std::string& expr, const std::string& cut, const std::string& title, const std::string& file){
            drawPlot(t, expr, cut, title + ";Angle (degrees);Count", folder_pos + file);
        };
        struct PosSet {std::string x, y, z, where, tag; int thetaNum;};
        const std::vector<PosSet> posSets = {{"fX", "fY", "fZ", "ZDC", "", 1}, {"fDx", "fDy", "fDz", "ZDC block", "_block", 3}};
        const bool doTheta = false; // theta/phi plots skipped for now
        if (doTheta) for (const auto& p : posSets) {
            const std::string thetaExpr = "TMath::ACos(" + p.z + "/sqrt(" + p.x + "*" + p.x + "+" + p.y + "*" + p.y + "+" + p.z + "*" + p.z + "))*TMath::RadToDeg()";
            const std::string phiExpr = "fmod(TMath::ATan2(" + p.y + "," + p.x + ")*TMath::RadToDeg()+360,360)";
            const std::string zoomCut = zdcProtonCut + " && " + thetaExpr + " < " + std::to_string(zoomThetaMax);
            const std::string th = "3_", ph = "4_", n = std::to_string(p.thetaNum), n2 = std::to_string(p.thetaNum + 1);
            const std::string hs = p.tag.empty() ? "" : "Block";

            drawAngle(events,  thetaExpr + ">>hTheta" + hs + "All(100,0,1.5)", zdcProtonCut, "Theta at " + p.where + ", CE events all protons", "2_allprotons/" + th + n + "_theta" + p.tag + ".pdf");
            drawAngle(events2, thetaExpr + ">>hTheta" + hs + "CE(100,0,1.5)", zdcProtonCut, "Theta at " + p.where + ", CE events CE protons", "1_CEprotons/" + th + n + "_theta" + p.tag + "_parentid.pdf");
            drawAngle(events,  thetaExpr + ">>hTheta" + hs + "ZoomAll(100,0,0.5)", zdcProtonCut, "Theta at " + p.where + ", CE events all protons" + zoomLabel, "2_allprotons/" + th + n2 + "_theta" + p.tag + "_zoom.pdf");
            drawAngle(events2, thetaExpr + ">>hTheta" + hs + "ZoomCE(100,0,0.5)", zdcProtonCut, "Theta at " + p.where + ", CE events CE protons" + zoomLabel, "1_CEprotons/" + th + n2 + "_theta" + p.tag + "_parentid_zoom.pdf");

            drawAngle(events,  phiExpr + ">>hPhi" + hs + "ZDCAll(100,0,360)", zdcProtonCut, "Phi at " + p.where + ", CE events all protons", "2_allprotons/" + ph + n + "_phi" + p.tag + ".pdf");
            drawAngle(events2, phiExpr + ">>hPhi" + hs + "ZDCCE(100,0,360)", zdcProtonCut, "Phi at " + p.where + ", CE events CE protons", "1_CEprotons/" + ph + n + "_phi" + p.tag + "_parentid.pdf");
            drawAngle(events,  phiExpr + ">>hPhi" + hs + "ZDCZoomAll(100,0,360)", zoomCut, "Phi at " + p.where + ", CE events all protons" + zoomLabel, "2_allprotons/" + ph + n2 + "_phi" + p.tag + "_zoom.pdf");
            drawAngle(events2, phiExpr + ">>hPhi" + hs + "ZDCZoomCE(100,0,360)", zoomCut, "Phi at " + p.where + ", CE events CE protons" + zoomLabel, "1_CEprotons/" + ph + n2 + "_phi" + p.tag + "_parentid_zoom.pdf");

            drawAngle(eventsAllEv, thetaExpr + ">>hTheta" + hs + "AllEv(100,0,1.5)", zdcProtonCut, "Theta at " + p.where + ", all events all protons", allEvDir + th + n + "_theta" + p.tag + ".pdf");
            drawAngle(eventsAllEv, thetaExpr + ">>hTheta" + hs + "ZoomAllEv(100,0,0.5)", zdcProtonCut, "Theta at " + p.where + ", all events all protons" + zoomLabel, allEvDir + th + n2 + "_theta" + p.tag + "_zoom.pdf");
            drawAngle(eventsAllEv, phiExpr + ">>hPhi" + hs + "ZDCAllEv(100,0,360)", zdcProtonCut, "Phi at " + p.where + ", all events all protons", allEvDir + ph + n + "_phi" + p.tag + ".pdf");
            drawAngle(eventsAllEv, phiExpr + ">>hPhi" + hs + "ZDCZoomAllEv(100,0,360)", zoomCut, "Phi at " + p.where + ", all events all protons" + zoomLabel, allEvDir + ph + n2 + "_phi" + p.tag + "_zoom.pdf");
        }

    }
}
