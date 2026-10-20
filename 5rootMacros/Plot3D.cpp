#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraph2D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TROOT.h>

#include <algorithm>
#include <functional>
#include <memory>
#include <TSystem.h>
#include <sstream>
#include <TStyle.h>

#include "FitCore.h"

void Plot3D(){
    std::string folder = "plots";

    std::string folder_pos = folder + "/1_posZDC/1_true/";
    std::string folder_posD = folder + "/1_posZDC/2_block/";

    std::string folder_energy = folder + "/4_energy/";


    gSystem->mkdir(folder.c_str(), true);
    gSystem->mkdir(folder_pos.c_str(), true);
    gSystem->mkdir(folder_posD.c_str(), true);
    gSystem->mkdir(folder_energy.c_str(), true);

    double Graph = 1;

    TFile* File = TFile::Open("cut_CE.root", "READ");
    TFile* File2 = TFile::Open("cut_ParentID.root", "READ");
    TFile* FileFake = TFile::Open("cut_FakeID.root", "READ");

    TTree *events = (TTree*)File->Get("ParticleInfo");
    TTree *ChargeExchange = (TTree*)File->Get("ChargeExchange_values");

    TTree *events2 = (TTree*)File2->Get("ParticleInfo");
    TTree *ChargeExchange2 = (TTree*)File2->Get("ChargeExchange_values");

    TTree *eventsFake = (TTree*)FileFake->Get("ParticleInfo");
    TTree *ChargeExchangeFake = (TTree*)FileFake->Get("ChargeExchange_values");


    if (Graph == 1){

        TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
        gPad->SetRightMargin(0.30);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); // room for the rotated Y title and bottom tick labels
        gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetNdivisions(505, "X"); // stats/fit box to the right, fewer x labels

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
        const std::string xySlice = zdcProtons + " && fabs(fY) < 2";
        const std::string xySliceD = zdcProtons + " && fabs(fDy) < 2";

        // ZDC position, true (fX, fY)
        drawXY(events, "fY:fX>>h2(100,-200,200,100,-200,200)", zdcProtons, "Position at ZDC, CE events all protons;X Position (mm);Y Position (mm)", folder_pos + "2_allprotons/1_1_xy_protons.pdf");
        drawXY(events, "fY:fX>>h2(100,-50,50,100,-50,50)", zdcProtons, "Position at ZDC, CE events all protons;X Position (mm);Y Position (mm)", folder_pos + "2_allprotons/1_2_xy_protons_zoom.pdf");
        drawXY(events2, "fY:fX>>h2(100,-200,200,100,-200,200)", zdcProtons, "Position at ZDC, CE events CE protons;X Position (mm);Y Position (mm)", folder_pos + "1_CEprotons/1_1_xy_parentid.pdf");
        drawXY(events2, "fY:fX>>h2(100,-50,50,100,-50,50)", zdcProtons, "Position at ZDC, CE events CE protons;X Position (mm);Y Position (mm)", folder_pos + "1_CEprotons/1_2_xy_parentid_zoom.pdf");

        gStyle->SetOptFit(1111); // probability, chi2/ndf, fit parameters with errors
        drawPlot(events2, "fX>>hSliceX200(100,-200,200)", xySlice, "X Position at ZDC, CE events CE protons, |Y| < 2 mm;X Position (mm);Count", folder_pos + "1_CEprotons/2_1_x_fit_200.pdf", "", [](TH1* h){FitCore(h);});
        drawPlot(events2, "fX>>hSliceX(100,-50,50)", xySlice, "X Position at ZDC, CE events CE protons, |Y| < 2 mm;X Position (mm);Count", folder_pos + "1_CEprotons/2_2_x_fit_50.pdf", "", [](TH1* h){FitCore(h);});
        drawPlot(events, "fX>>hSliceXAll200(100,-200,200)", xySlice, "X Position at ZDC, CE events all protons, |Y| < 2 mm;X Position (mm);Count", folder_pos + "2_allprotons/2_1_x_fit_200.pdf", "", [](TH1* h){FitCore(h);});
        drawPlot(events, "fX>>hSliceXAll(100,-50,50)", xySlice, "X Position at ZDC, CE events all protons, |Y| < 2 mm;X Position (mm);Count", folder_pos + "2_allprotons/2_2_x_fit_50.pdf", "", [](TH1* h){FitCore(h);});


        // ZDC position, block (fDx, fDy)
        drawXY(events, "fDy:fDx>>h2(100,-200,200,100,-200,200)", zdcProtons, "Position at ZDC BLOCK, CE events all protons;X Position (mm);Y Position (mm)", folder_posD + "2_allprotons/1_1_xy_protons.pdf");
        drawXY(events, "fDy:fDx>>h2(100,-50,50,100,-50,50)", zdcProtons, "Position at ZDC BLOCK, CE events all protons;X Position (mm);Y Position (mm)", folder_posD + "2_allprotons/1_2_xy_protons_zoom.pdf");
        drawXY(events2, "fDy:fDx>>h2(100,-200,200,100,-200,200)", zdcProtons, "Position at ZDC BLOCK, CE events CE protons;X Position (mm);Y Position (mm)", folder_posD + "1_CEprotons/1_1_xy_parentid.pdf");
        drawXY(events2, "fDy:fDx>>h2(100,-50,50,100,-50,50)", zdcProtons, "Position at ZDC BLOCK, CE events CE protons;X Position (mm);Y Position (mm)", folder_posD + "1_CEprotons/1_2_xy_parentid_zoom.pdf");

        drawPlot(events2, "fDx>>hSliceDX200(100,-200,200)", xySliceD, "X Position at ZDC BLOCK, CE events CE protons, |Y| < 2 mm;X Position (mm);Count", folder_posD + "1_CEprotons/2_1_x_fit_200.pdf", "", [](TH1* h){FitBlocks(h);}); //one bin per ZDC block
        drawPlot(events2, "fDx>>hSliceDX(100,-127.27,127.27)", xySliceD, "X Position at ZDC BLOCK, CE events CE protons, |Y| < 2 mm;X Position (mm);Count", folder_posD + "1_CEprotons/2_2_x_fit_50.pdf", "", [](TH1* h){FitBlocks(h);}); //7 central ZDC blocks, pitch 400/11 mm
        drawPlot(events, "fDx>>hSliceDXAll200(100,-200,200)", xySliceD, "X Position at ZDC BLOCK, CE events all protons, |Y| < 2 mm;X Position (mm);Count", folder_posD + "2_allprotons/2_1_x_fit_200.pdf", "", [](TH1* h){FitBlocks(h);});
        drawPlot(events, "fDx>>hSliceDXAll(100,-127.27,127.27)", xySliceD, "X Position at ZDC BLOCK, CE events all protons, |Y| < 2 mm;X Position (mm);Count", folder_posD + "2_allprotons/2_2_x_fit_50.pdf", "", [](TH1* h){FitBlocks(h);});

        // Energy plots: smearing/<set>/<quantity>/ with set = 2_allprotons, 1_CEprotons and quantity = ke, 2d_vs_true, diff_ratio_pct_error
        const std::string folder_smear = folder_energy + "1_smearning/";
        // True kinetic energy at the ZDC
        double zoomMax = 0; // y max of the proton zoom, reused for the same-ParentID plot so the two share a scale
        drawPlot(events, "fKEnergy", zdcProtons, "Kinetic Energy at ZDC, CE events all protons;Kinetic Energy (MeV);Count", folder_smear + "2_allprotons/1_true_ke_full.pdf");
        // The CE protons plot below shares the y scale of the all protons plot zoomed on the ~24000 MeV peak
        events->Draw("fKEnergy>>hKEScale", zdcProtons.c_str(), "goff");
        TH1* hScale = (TH1*)gDirectory->Get("hKEScale");
        hScale->GetXaxis()->SetRangeUser(20000, 28000);
        zoomMax = 1.1 * hScale->GetBinContent(hScale->GetMaximumBin());
        drawPlot(events2, "fKEnergy", zdcProtons, "Kinetic Energy at ZDC, CE events CE protons;Kinetic Energy (MeV);Count", folder_smear + "1_CEprotons/1_true_ke_full.pdf", "",
                 [&](TH1* h){h->SetMaximum(zoomMax);});

        // ZDC energy smearing (fDEnergy = Gaussian smear of fKEnergy, /zdc/resStochastic and /zdc/resConstant in the sim)
        if (events->GetBranch("fDEnergy")){ // files from before fDEnergy was added don't have it
            gSystem->mkdir(folder_smear.c_str(), true);

            // fKEnergy = true KE, fDEnergy = smeared KE, dE = fDEnergy - fKEnergy

            drawPlot(events, "fDEnergy-fKEnergy>>hSmearDiffP(100,-10000,10000)", zdcProtons, "Energy Difference at ZDC, CE events all protons;fDEnergy - fKEnergy (MeV);Count", folder_smear + "2_allprotons/4_1_diff_protons.pdf", "", [](TH1* h){FitCore(h);});
            drawPlot(events2, "fDEnergy-fKEnergy>>hSmearDiffCE(100,-10000,10000)", zdcProtons, "Energy Difference at ZDC, CE events CE protons;fDEnergy - fKEnergy (MeV);Count", folder_smear + "1_CEprotons/4_1_diff_parentid.pdf", "", [](TH1* h){FitCore(h);});

            // percent error = 100 * (fDEnergy - fKEnergy) / fKEnergy
            const std::string pct = "100*(fDEnergy-fKEnergy)/fKEnergy";
            drawPlot(events, pct + ">>hSmearPctP(100,-100,100)", zdcProtons, "Percent Error at ZDC, CE events all protons;Percent Error (%);Count", folder_smear + "2_allprotons/4_2_pct_protons.pdf", "", [](TH1* h){FitCore(h);});
            drawPlot(events2, pct + ">>hSmearPctCE(100,-100,100)", zdcProtons, "Percent Error at ZDC, CE events CE protons;Percent Error (%);Count", folder_smear + "1_CEprotons/4_2_pct_parentid.pdf", "", [](TH1* h){FitCore(h);});

            // Same KE plots for CE protons and all protons: beam energy zoom, and the range of the CE protons' true KE at the ZDC
            events2->SetEstimate(events2->GetEntries() + 1);
            const Long64_t nCE = events2->Draw("fKEnergy", zdcProtons.c_str(), "goff");
            const std::string range = nCE > 0 ? Form("(100,%g,%g)", *std::min_element(events2->GetV1(), events2->GetV1() + nCE),
                                                                      *std::max_element(events2->GetV1(), events2->GetV1() + nCE)) : "";
            struct KESet {TTree* tree; std::string dir, label;};
            for (const auto& set : {KESet{events2, "1_CEprotons/", "CE protons"}, KESet{events, "2_allprotons/", "all protons"}}){
                const std::string where = ", CE events " + set.label;
                const std::string axisTrue = ";Kinetic Energy (MeV);Count", axisSmear = ";Smeared Energy fDEnergy (MeV);Count";
                drawPlot(set.tree, "fKEnergy>>hTrueKEZoom(100,20000,30000)", zdcProtons, "True Kinetic Energy at ZDC" + where + axisTrue, folder_smear + set.dir + "2_1_true_ke_20to30GeV.pdf");
                drawPlot(set.tree, "fDEnergy>>hSmearedKEZoom(100,20000,30000)", zdcProtons, "Smeared Energy at ZDC" + where + axisSmear, folder_smear + set.dir + "2_2_smeared_ke_20to30GeV.pdf");
                if (range.empty()) continue;
                drawPlot(set.tree, "fKEnergy>>hTrueKERange" + range, zdcProtons, "True Kinetic Energy at ZDC" + where + ", CE proton KE range" + axisTrue, folder_smear + set.dir + "3_1_true_ke_ce_range.pdf");
                drawPlot(set.tree, "fDEnergy>>hSmearedKERange" + range, zdcProtons, "Smeared Energy at ZDC" + where + ", CE proton KE range" + axisSmear, folder_smear + set.dir + "3_2_smeared_ke_ce_range.pdf");
            }
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
                std::string folder_cut = folder_energy + "ke_cut_" + std::to_string(pct) + "/";
                gSystem->mkdir(folder_cut.c_str(), true);
                const std::string axis = ";" + e.label + " (MeV);Count";
                drawPlot(events2, e.branch + ">>hCutReal(100,10000,30000)", zdcProtons + " && " + e.branch + " < " + cut, e.label + " of Protons Left Out with " + pctStr + " cut, CE events CE protons" + axis, folder_cut + n + "1_real_cut_" + e.dir + ".pdf");
                drawPlot(eventsFake, e.branch + ">>hCutFake(100,10000,30000)", zdcProtons + " && " + e.branch + " > " + cut, e.label + " of Protons Kept In with " + pctStr + " cut, CE events fake protons" + axis, folder_cut + n + "2_fake_cut_" + e.dir + ".pdf");
                // real (CE protons) and fake (other protons in CE events) on one plot, with the cut and what it does to each
                const double cutVal = initialKE * pct / 100.0;
                const std::string realOutCut = zdcProtons + " && " + e.branch + " < " + cut;
                const std::string fakeKeptCut = zdcProtons + " && " + e.branch + " > " + cut;
                const double realOutPct = 100.0 * events2->GetEntries(realOutCut.c_str()) / std::max<Long64_t>(1, events2->GetEntries(zdcProtons.c_str()));
                const double fakeKeptPct = 100.0 * eventsFake->GetEntries(fakeKeptCut.c_str()) / std::max<Long64_t>(1, eventsFake->GetEntries(zdcProtons.c_str()));
                drawPlot(events2, e.branch + ">>hCutRealOv(100,10000,30000)", zdcProtons, e.label + " of Protons, real vs fake, " + pctStr + " cut, CE events" + axis, folder_cut + n + "3_real_vs_fake_" + e.dir + ".pdf", "", [&](TH1* h){
                    h->SetLineColor(kBlue);
                    eventsFake->Draw((e.branch + ">>hCutFakeOv(100,10000,30000)").c_str(), zdcProtons.c_str(), "same");
                    TH1* hf = static_cast<TH1*>(gDirectory->Get("hCutFakeOv"));
                    hf->SetLineColor(kRed);
                    h->SetMaximum(1.15 * std::max(h->GetMaximum(), hf->GetMaximum()));
                    TLine* line = new TLine(cutVal, 0, cutVal, h->GetMaximum());
                    line->SetLineStyle(2); line->Draw();
                    TLegend* leg = new TLegend(0.13, 0.76, 0.50, 0.88); leg->SetTextSize(0.03);
                    leg->AddEntry(h, Form("real CE protons, %.1f%% left out", realOutPct), "l");
                    leg->AddEntry(hf, Form("fake protons, %.1f%% kept in", fakeKeptPct), "l");
                    leg->Draw();
                });
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
        for (const auto& p : posSets) {
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
        }

    }
}
