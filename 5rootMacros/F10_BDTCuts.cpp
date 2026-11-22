#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "FitCore.h"

// 10_BDTcuts: the signal vs background BDT of TmvaTrain.cpp / TmvaApply.cpp, cut on its score. Run in the run folder (build/root/<run>/) after F9_PlotReconNew.cpp: root -l -b -q 'F10_BDTCuts.cpp'
// 1. The polarization fits and the signal and background efficiencies of the 36 BDT cuts -0.50 to 0.20 in steps of 0.02 (plots/10_BDTcuts/scan_csv/BDT<cut>/scan_stats.csv of F9_PlotReconNew.cpp, each with the BDT cut as its only cut):
//    1_P_vs_BDTcut.pdf, 2_signalEff_vs_bkgEff.pdf, 3_cut_metrics/ (one pdf each of signal efficiency, background efficiency, background rejection, purity, contamination, significance), 4_FOM_vs_BDTcut.pdf and bdt_cut_table.csv.
// 2. The TMVA training plots, drawn from resultTMVAtraining_SigBkg.root (TmvaTrain.cpp) into plots/10_BDTcuts/TMVA/: BDT response with the overtraining check, ROC curve, cut efficiencies, the input variables (input_variables/) and the correlation matrices.
// Signal = CE only tracks (evClass 0, every fitted hit a CE child), background = no CE child hit (evClass 2); CE mix (evClass 1) is in neither. (The TMVA training itself used isFromCE == 1 rows as signal.) Both read only the csv and root files, no sample is read again.
namespace {
// figure of merit as F5_FomCalc.cpp: N_S / N_neutrons * A_eff^2 (N_S = signal tracks passing the cut, A_eff the mean at the recon block theta of those tracks)
double Fom(const CutRow& r){return r.totN > 0 ? r.nS / r.totN * r.meanA * r.meanA : 0.0;}

int SignalColor(){return CbCE();}
int BackgroundColor(){return CbAll();}

void DrawBdtCuts(TCanvas& c, const std::string& dir){
    std::vector<CutRow> rows;
    for (int hund = -50; hund <= 20; hund += 2){ // the cuts in hundredths, as F9_PlotReconNew.cpp
        CutRow r;
        r.cut = hund / 100.0;
        if (ReadCutRow(Form("%s/scan_csv/BDT%.2f/scan_stats.csv", dir.c_str(), r.cut), r)) rows.push_back(r);
    }
    if (rows.empty()){Warning("F10_BDTCuts", "no %s/scan_csv/BDT<cut>/scan_stats.csv (run F9_PlotReconNew.cpp with cuts/tmva_scores.root)", dir.c_str()); return;}

    {
        std::ofstream t(dir + "/bdt_cut_table.csv");
        t.precision(10);
        t << "cut,sAll,bAll,nS,nB,sigEff,bkgEff,bkgRejection,purity,contamination,significance,FOM,nFitCEMix,pCEMix,errCEMix,nFitAll,pAll,errAll,nFitCE,pCE,errCE,PPValue\n";
        for (const CutRow& r : rows)
            t << r.cut << "," << r.sAll << "," << r.bAll << "," << r.nS << "," << r.nB << "," << (r.sAll > 0 ? r.nS / r.sAll : 0.0) << "," << (r.bAll > 0 ? r.nB / r.bAll : 0.0) << "," << (r.bAll > 0 ? 1.0 - r.nB / r.bAll : 0.0) << "," << (r.nS + r.nB > 0 ? r.nS / (r.nS + r.nB) : 0.0) << "," << (r.nS + r.nB > 0 ? r.nB / (r.nS + r.nB) : 0.0) << "," << (r.nS + r.nB > 0 ? r.nS / std::sqrt(r.nS + r.nB) : 0.0) << "," << Fom(r) << ","
              << r.nFitMix << "," << r.pMix << "," << r.eMix << "," << r.nFitAll << "," << r.pAll << "," << r.eAll << "," << r.nFitCE << "," << r.pCE << "," << r.eCE << "," << r.ppv << "\n";
    }

    // 1_P_vs_BDTcut: the P of the CE, CE + CE mix and total recon fits (the unbinned phi fit of F9_PlotReconNew.cpp, error with the beam A_P error) against the BDT cut; the dotted line is the true P (PPValue)
    {
        const std::string path = dir + "/1_P_vs_BDTcut.pdf";
        const char* name[3] = {"CE recon protons (fitted at the true CE angle)", "CE + CE mix recon protons", "total recon particles"};
        const int color[3] = {CbCE(), CbProton(), CbAll()}, marker[3] = {33, 20, 21};
        const double shift[3] = {-0.005, 0.0, 0.005};
        TGraphErrors* g[3];
        double lo = rows[0].ppv, hi = rows[0].ppv;
        for (int s = 0; s < 3; s++){
            g[s] = new TGraphErrors();
            for (const CutRow& r : rows){
                const double n = s == 0 ? r.nFitCE : s == 1 ? r.nFitMix : r.nFitAll, p = s == 0 ? r.pCE : s == 1 ? r.pMix : r.pAll, e = s == 0 ? r.eCE : s == 1 ? r.eMix : r.eAll;
                if (n < 10) continue; // no fit
                const int i = g[s]->GetN();
                g[s]->SetPoint(i, r.cut + shift[s], p); g[s]->SetPointError(i, 0, e);
                lo = std::min(lo, p - e); hi = std::max(hi, p + e);
            }
            g[s]->SetLineColor(color[s]); g[s]->SetMarkerColor(color[s]); g[s]->SetMarkerStyle(marker[s]); g[s]->SetMarkerSize(1.3); g[s]->SetLineWidth(2);
        }
        const double pad = 0.1 * (hi - lo > 0 ? hi - lo : 0.1);
        const double xLo = rows.front().cut - 0.05, xHi = rows.back().cut + 0.05;
        TH1F frame("hFrameP", "", 1, xLo, xHi);
        frame.SetStats(0); frame.SetMinimum(-1.0); frame.SetMaximum(1.0); // P is between -1 and 1; points with a larger error bar run off the frame
        frame.SetTitle(PlotTitle(path, "polarization fit P vs BDT score cut, BDT score >= cut is the only cut (no other cut);BDT score cut;polarization fit P").c_str());
        c.Clear(); c.SetLogx(false); c.SetLogy(false);
        frame.Draw("axis");
        TLine lt(xLo, rows[0].ppv, xHi, rows[0].ppv); lt.SetLineStyle(3); lt.SetLineWidth(3); lt.Draw();
        TLegend leg(0.12, 0.13, 0.70, 0.27); // the low part of the frame is empty
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        for (int s = 0; s < 3; s++){g[s]->Draw("PZ same"); leg.AddEntry(g[s], name[s], "pl");}
        TLine* dummy = new TLine(); dummy->SetLineStyle(3); dummy->SetLineWidth(3);
        leg.AddEntry(dummy, Form("true P (PPValue) = %.3g", rows[0].ppv), "l");
        leg.Draw();
        c.SaveAs(path.c_str());
        for (int s = 0; s < 3; s++) delete g[s];
    }

    // 2_signalEff_vs_bkgEff: the share of the signal and of the background recon tracks that pass each BDT cut (of every track, no other cut); the points are labelled with the cut
    {
        const std::string path = dir + "/2_signalEff_vs_bkgEff.pdf";
        TGraph g;
        double xMin = 1.0;
        for (const CutRow& r : rows){
            if (r.sAll <= 0 || r.bAll <= 0) continue;
            const double bkg = r.nB / r.bAll;
            if (bkg > 0) xMin = std::min(xMin, bkg);
            g.SetPoint(g.GetN(), bkg > 0 ? bkg : 1e-9, r.nS / r.sAll);
        }
        TH1F frame("hFrameE", "", 1, 0.5 * xMin, 1.5);
        frame.SetStats(0); frame.SetMinimum(0.0); frame.SetMaximum(1.08);
        frame.SetTitle(PlotTitle(path, "signal efficiency vs background efficiency of the BDT score cuts -0.50 to 0.20, no other cut;background efficiency (background tracks passing / all background tracks);signal efficiency (signal tracks passing / all signal tracks)").c_str());
        c.Clear(); c.SetLogx(true); c.SetLogy(false);
        frame.Draw("axis");
        g.SetLineColor(SignalColor()); g.SetMarkerColor(SignalColor()); g.SetMarkerStyle(20); g.SetMarkerSize(1.3); g.SetLineWidth(2);
        g.Draw("PL same");
        TLatex tl; tl.SetTextSize(0.028); tl.SetTextAlign(13);
        size_t j = 0;
        for (const CutRow& r : rows){
            if (r.sAll <= 0 || r.bAll <= 0) continue;
            double x, y; g.GetPoint((int)j++, x, y);
            if (std::lround(r.cut * 100) % 10 != 0) continue; // a label every 0.1 only, the 0.02 steps would overlap
            tl.DrawLatex(x * 1.12, y - 0.015, Form("%.2f", r.cut));
        }
        c.SaveAs(path.c_str());
        c.SetLogx(false);
    }
    // 2a_signalEff_vs_bkgRejection: the ROC curve of the cuts, background rejection (1 - background efficiency) against signal efficiency, the points labelled with the cut every 0.1
    {
        const std::string path = dir + "/2a_signalEff_vs_bkgRejection.pdf";
        TGraph g;
        std::vector<double> cutOfPoint;
        for (const CutRow& r : rows){
            if (r.sAll <= 0 || r.bAll <= 0) continue;
            g.SetPoint(g.GetN(), r.nS / r.sAll, 1.0 - r.nB / r.bAll);
            cutOfPoint.push_back(r.cut);
        }
        TH1F frame("hFrameRoc", "", 1, 0.0, 1.05);
        frame.SetStats(0); frame.SetMinimum(0.0); frame.SetMaximum(1.08);
        frame.SetTitle(PlotTitle(path, "background rejection vs signal efficiency of the BDT score cuts -0.50 to 0.20, signal = CE only tracks (evClass 0), background = non-CE tracks (evClass 2), no other cut;signal efficiency;background rejection (1 - background efficiency)").c_str());
        c.Clear(); c.SetLogx(false); c.SetLogy(false); c.SetGrid();
        frame.Draw("axis");
        g.SetLineColor(SignalColor()); g.SetMarkerColor(SignalColor()); g.SetMarkerStyle(20); g.SetMarkerSize(1.0); g.SetLineWidth(2);
        g.Draw("PL same");
        TLatex tl; tl.SetTextSize(0.028); tl.SetTextAlign(13);
        for (int j = 0; j < g.GetN(); j++){
            if (std::lround(cutOfPoint[j] * 100) % 10 != 0) continue; // a label every 0.1 only
            double x, y; g.GetPoint(j, x, y);
            tl.DrawLatex(x + 0.01, y - 0.015, Form("%.2f", cutOfPoint[j]));
        }
        c.SaveAs(path.c_str());
        c.SetGrid(0, 0);
    }

    // 3_cut_metrics/: one pdf per quantity against the BDT cut, from the tracks passing it (no other cut): signal efficiency, background efficiency (log y), background rejection (1 - background efficiency),
    // purity S / (S + B), contamination B / (S + B) and the significance S / sqrt(S + B); S = CE only tracks (evClass 0) passing, B = non-CE tracks (evClass 2) passing
    {
        gSystem->mkdir((dir + "/3_cut_metrics").c_str(), true);
        struct Metric {const char* file; const char* title; const char* yTitle; bool logY; int marker;};
        const Metric mt[6] = {
            {"1_signalEfficiency", "signal efficiency, CE only tracks passing / all CE only tracks", "signal efficiency", false, 20},
            {"2_backgroundEfficiency", "background efficiency, non-CE tracks passing / all non-CE tracks", "background efficiency", true, 21},
            {"3_backgroundRejection", "background rejection, 1 - background efficiency", "background rejection", false, 33},
            {"4_purity", "purity S / (S + B)", "purity", false, 22},
            {"5_contamination", "contamination B / (S + B)", "contamination", false, 23},
            {"6_significance", "significance S / #sqrt{S + B}", "significance", false, 34}};
        const int color[6] = {CbCE(), CbAll(), CbProton(), kBlack, kMagenta + 1, kOrange + 7};
        for (int q = 0; q < 6; q++){
            const std::string path = dir + "/3_cut_metrics/" + mt[q].file + "_vs_BDTcut.pdf";
            TGraph g;
            double top = 0;
            for (const CutRow& r : rows){
                const double sb = r.nS + r.nB;
                const double y = q == 0 ? (r.sAll > 0 ? r.nS / r.sAll : 0) : q == 1 ? (r.bAll > 0 ? r.nB / r.bAll : 0) : q == 2 ? (r.bAll > 0 ? 1.0 - r.nB / r.bAll : 0) : q == 3 ? (sb > 0 ? r.nS / sb : 0) : q == 4 ? (sb > 0 ? r.nB / sb : 0) : (sb > 0 ? r.nS / std::sqrt(sb) : 0);
                if (mt[q].logY && y <= 0) continue; // log axis
                g.SetPoint(g.GetN(), r.cut, y); top = std::max(top, y);
            }
            if (g.GetN() == 0) continue;
            double ymin = 0.0, ymax = q < 5 ? 1.08 : 1.15 * top;
            if (mt[q].logY){ymin = 1e-8; ymax = 3.0;}
            TH1F frame(Form("hFrameM%d", q), "", 1, rows.front().cut - 0.05, rows.back().cut + 0.05);
            frame.SetStats(0); frame.SetMinimum(ymin); frame.SetMaximum(ymax);
            frame.SetTitle(PlotTitle(path, std::string(mt[q].title) + " vs BDT score cut, BDT score >= cut is the only cut (no other cut), signal = CE only tracks (evClass 0), background = non-CE tracks (evClass 2);BDT score cut;" + mt[q].yTitle).c_str());
            c.Clear(); c.SetLogx(false); c.SetLogy(mt[q].logY);
            frame.Draw("axis");
            g.SetLineColor(color[q]); g.SetMarkerColor(color[q]); g.SetMarkerStyle(mt[q].marker); g.SetMarkerSize(1.3); g.SetLineWidth(2);
            g.Draw("PL same");
            c.SaveAs(path.c_str());
            c.SetLogy(false);
        }
    }

    // 4_FOM_vs_BDTcut: the figure of merit of F5_FomCalc.cpp, N_S / N_neutrons x A_eff^2, of the signal tracks passing each cut
    {
        const std::string path = dir + "/4_FOM_vs_BDTcut.pdf";
        TGraph g;
        double top = 0;
        for (const CutRow& r : rows){g.SetPoint(g.GetN(), r.cut, Fom(r)); top = std::max(top, Fom(r));}
        TH1F frame("hFrameFom", "", 1, rows.front().cut - 0.05, rows.back().cut + 0.05);
        frame.SetStats(0); frame.SetMinimum(0.0); frame.SetMaximum(1.15 * (top > 0 ? top : 1.0));
        frame.SetTitle(PlotTitle(path, "FOM vs BDT score cut, FOM = N_{S} / N_{neutrons} #times A_{eff}^{2}, N_{S} = signal tracks passing, BDT score >= cut is the only cut (no other cut);BDT score cut;FOM").c_str());
        c.Clear(); c.SetLogx(false); c.SetLogy(false);
        frame.Draw("axis");
        g.SetLineColor(CbCE()); g.SetMarkerColor(CbCE()); g.SetMarkerStyle(20); g.SetMarkerSize(1.3); g.SetLineWidth(2);
        g.Draw("PL same");
        c.SaveAs(path.c_str());
    }
}

// a pdf of one histogram pair, signal and background, normalised to unit area
void DrawPair(TCanvas& c, TH1* s, TH1* b, const std::string& path, const std::string& title){
    TH1* hs = (TH1*)s->Clone(Form("%s_cs", s->GetName())); TH1* hb = (TH1*)b->Clone(Form("%s_cb", b->GetName()));
    hs->SetDirectory(nullptr); hb->SetDirectory(nullptr);
    if (hs->Integral() > 0) hs->Scale(1.0 / hs->Integral());
    if (hb->Integral() > 0) hb->Scale(1.0 / hb->Integral());
    hs->SetStats(0); hb->SetStats(0);
    hs->SetLineColor(SignalColor()); hs->SetLineWidth(2); hb->SetLineColor(BackgroundColor()); hb->SetLineWidth(2);
    hs->SetTitle(PlotTitle(path, title).c_str());
    hs->SetMaximum(1.35 * std::max(hs->GetMaximum(), hb->GetMaximum())); hs->SetMinimum(0.0);
    c.Clear(); c.SetLogy(false);
    hs->Draw("hist"); hb->Draw("hist same");
    TLegend leg(0.62, 0.78, 0.96, 0.90);
    leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(hs, "signal (fit used a CE child hit)", "l"); leg.AddEntry(hb, "background (no CE child hit)", "l");
    leg.Draw();
    c.SaveAs(path.c_str());
    delete hs; delete hb;
}

void DrawTmva(TCanvas& c, const char* file, const std::string& dir){
    TFile f(file);
    if (f.IsZombie()){Warning("F10_BDTCuts", "no %s (run TmvaTrain.cpp), no TMVA plots", file); return;}
    TDirectory* bd = (TDirectory*)f.Get("dataset/Method_BDT/BDT");
    if (!bd){Warning("F10_BDTCuts", "no dataset/Method_BDT/BDT in %s, no TMVA plots", file); return;}
    const std::string tdir = dir + "/TMVA";
    gSystem->mkdir((tdir + "/input_variables").c_str(), true);
    const std::string sample = "TMVA test sample (half of the TMVA input events), no cuts";
    gStyle->SetOptStat(0);

    // 1_BDT_response: the BDT score of the test sample (histograms) and the training sample (markers); the Kolmogorov-Smirnov probabilities say whether the two agree (overtraining check)
    TH1D *s = (TH1D*)bd->Get("MVA_BDT_S"), *b = (TH1D*)bd->Get("MVA_BDT_B"), *ts = (TH1D*)bd->Get("MVA_BDT_Train_S"), *tb = (TH1D*)bd->Get("MVA_BDT_Train_B");
    if (s && b && ts && tb){
        const std::string path = tdir + "/1_BDT_response.pdf";
        s->SetStats(0);
        s->SetLineColor(SignalColor()); s->SetFillColorAlpha(SignalColor(), 0.3); s->SetLineWidth(2);
        b->SetLineColor(BackgroundColor()); b->SetFillColorAlpha(BackgroundColor(), 0.3); b->SetLineWidth(2);
        ts->SetMarkerStyle(20); ts->SetMarkerColor(SignalColor()); ts->SetLineColor(SignalColor()); tb->SetMarkerStyle(21); tb->SetMarkerColor(BackgroundColor()); tb->SetLineColor(BackgroundColor());
        s->SetTitle(PlotTitle(path, "BDT response, test sample (histograms) and training sample (markers), no cuts;BDT score;(1/N) dN/dx").c_str());
        s->SetMaximum(1.5 * std::max(std::max(s->GetMaximum(), b->GetMaximum()), std::max(ts->GetMaximum(), tb->GetMaximum()))); s->SetMinimum(0.0);
        c.Clear(); c.SetLogy(false);
        s->Draw("hist"); b->Draw("hist same"); ts->Draw("E same"); tb->Draw("E same");
        TLegend leg(0.12, 0.74, 0.55, 0.90);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(s, "signal, test", "f"); leg.AddEntry(b, "background, test", "f"); leg.AddEntry(ts, "signal, training", "p"); leg.AddEntry(tb, "background, training", "p");
        leg.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026);
        tx.DrawLatex(0.58, 0.86, Form("KS test vs training  signal %.3f, background %.3f", s->KolmogorovTest(ts), b->KolmogorovTest(tb)));
        c.SaveAs(path.c_str());
    }

    // 2_ROC: background rejection against signal efficiency (the test sample, the training sample dashed); the ROC integral is the area under the curve
    TH1D *rej = (TH1D*)bd->Get("MVA_BDT_rejBvsS"), *trej = (TH1D*)bd->Get("MVA_BDT_trainingRejBvsS");
    if (rej){
        const std::string path = tdir + "/2_ROC.pdf";
        rej->SetStats(0); rej->SetLineColor(SignalColor()); rej->SetLineWidth(3);
        rej->SetTitle(PlotTitle(path, "ROC curve of the BDT, test sample (solid) and training sample (dashed), no cuts;signal efficiency;background rejection").c_str());
        rej->SetMinimum(0.0); rej->SetMaximum(1.1);
        c.Clear(); c.SetLogy(false); c.SetGrid();
        rej->Draw("hist");
        if (trej){trej->SetLineColor(kBlack); trej->SetLineStyle(2); trej->SetLineWidth(2); trej->Draw("hist same");}
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.03);
        std::string roc = Form("ROC integral, test %.3f", rej->Integral("width")); // not nested in one Form: it reuses one buffer
        if (trej) roc += Form(", training %.3f", trej->Integral("width"));
        tx.DrawLatex(0.15, 0.30, roc.c_str());
        c.SaveAs(path.c_str());
        c.SetGrid(0, 0);
    }

    // 3_cut_efficiencies: the share of the test signal and background above each BDT cut
    TH1D *es = (TH1D*)bd->Get("MVA_BDT_effS"), *eb = (TH1D*)bd->Get("MVA_BDT_effB");
    if (es && eb){
        const std::string path = tdir + "/3_cut_efficiencies.pdf";
        es->SetStats(0); es->SetLineColor(SignalColor()); es->SetLineWidth(3); eb->SetLineColor(BackgroundColor()); eb->SetLineWidth(3);
        es->SetTitle(PlotTitle(path, "BDT cut efficiencies, " + sample + ";BDT score cut;efficiency (events with a score above the cut / all events)").c_str());
        es->SetMinimum(0.0); es->SetMaximum(1.15);
        c.Clear(); c.SetLogy(false);
        es->Draw("hist"); eb->Draw("hist same");
        TLegend leg(0.62, 0.78, 0.96, 0.90);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(es, "signal efficiency", "l"); leg.AddEntry(eb, "background efficiency", "l");
        leg.Draw();
        c.SaveAs(path.c_str());
    }

    // 5a_classifier_cut_efficiencies (TMVAGui 5a, mvaeffs): the signal and background efficiency, the signal purity S / (S + B), efficiency x purity and the significance S / sqrt(S + B) against the BDT cut,
    // with the events of the test sample as S and B (the efficiencies times the number of test signal and background events); the significance is drawn scaled to its maximum of 1, the maximum is marked
    if (es && eb && s && b){
        const std::string path = tdir + "/5a_classifier_cut_efficiencies.pdf";
        const double nSig = s->GetEntries(), nBkg = b->GetEntries();
        TH1D* pu = (TH1D*)es->Clone("hMvaPurity"); TH1D* ep = (TH1D*)es->Clone("hMvaEffPur"); TH1D* sg = (TH1D*)es->Clone("hMvaSignif");
        pu->SetDirectory(nullptr); ep->SetDirectory(nullptr); sg->SetDirectory(nullptr);
        double sgMax = 0, cutAtMax = 0;
        for (int i = 1; i <= es->GetNbinsX(); i++){
            const double S = es->GetBinContent(i) * nSig, B = eb->GetBinContent(i) * nBkg;
            const double pur = S + B > 0 ? S / (S + B) : 0.0, sig = S + B > 0 ? S / std::sqrt(S + B) : 0.0;
            pu->SetBinContent(i, pur); ep->SetBinContent(i, es->GetBinContent(i) * pur); sg->SetBinContent(i, sig);
            if (sig > sgMax){sgMax = sig; cutAtMax = es->GetBinCenter(i);}
        }
        if (sgMax > 0) sg->Scale(1.0 / sgMax);
        TH1D* se = (TH1D*)es->Clone("hMvaEffS"); TH1D* be = (TH1D*)eb->Clone("hMvaEffB");
        se->SetDirectory(nullptr); be->SetDirectory(nullptr);
        const int col[5] = {SignalColor(), BackgroundColor(), kBlack, kMagenta + 1, CbProton()};
        TH1D* h[5] = {se, be, pu, ep, sg};
        const char* nm[5] = {"signal efficiency", "background efficiency", "signal purity S/(S+B)", "signal efficiency #times purity", "significance S/#sqrt{S+B} (scaled to its maximum)"};
        for (int q = 0; q < 5; q++){h[q]->SetStats(0); h[q]->SetLineColor(col[q]); h[q]->SetLineWidth(3);}
        se->SetTitle(PlotTitle(path, "classifier cut efficiencies vs BDT cut, " + sample + ";BDT score cut;efficiency, purity, scaled significance").c_str());
        se->SetMinimum(0.0); se->SetMaximum(1.35);
        c.Clear(); c.SetLogy(false);
        se->Draw("hist");
        for (int q = 1; q < 5; q++) h[q]->Draw("hist same");
        TLegend leg(0.12, 0.76, 0.70, 0.90);
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetNColumns(1);
        for (int q = 0; q < 5; q++) leg.AddEntry(h[q], nm[q], "l");
        leg.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026);
        tx.DrawLatex(0.14, 0.70, Form("maximum significance %.1f at BDT cut %.3f (S = %.0f test signal, B = %.0f test background events)", sgMax, cutAtMax, nSig, nBkg));
        c.SaveAs(path.c_str());
        for (int q = 0; q < 5; q++) delete h[q];
    }

    // 4_ and 5_correlation: the linear correlations of the input variables
    const char* corr[2] = {"CorrelationMatrixS", "CorrelationMatrixB"};
    const char* corrName[2] = {"signal", "background"};
    for (int k = 0; k < 2; k++){
        TH2F* h = (TH2F*)f.Get(Form("dataset/%s", corr[k]));
        if (!h) continue;
        const std::string path = Form("%s/%d_correlation_%s.pdf", tdir.c_str(), k + 4, corrName[k]);
        h->SetStats(0); h->SetMarkerSize(1.0);
        h->SetTitle(PlotTitle(path, std::string("linear correlation of the BDT input variables, ") + corrName[k] + ", TMVA training sample, no cuts").c_str());
        gStyle->SetPaintTextFormat("4.0f");
        c.Clear(); c.SetLogy(false); c.SetLeftMargin(0.22); c.SetBottomMargin(0.22); c.SetRightMargin(0.16);
        h->Draw("COLZ TEXT");
        c.SaveAs(path.c_str());
        c.SetLeftMargin(0.13); c.SetBottomMargin(0.13); c.SetRightMargin(0.04);
    }

    // input_variables: every input of the BDT, signal and background normalised to unit area
    TDirectory* iv = (TDirectory*)f.Get("dataset/InputVariables_Id");
    if (iv){
        int n = 0;
        TIter next(iv->GetListOfKeys());
        while (TKey* key = (TKey*)next()){
            const std::string nm = key->GetName(), tag = "__Signal_Id"; // TMVA names them <variable>__Signal_Id and <variable>__Background_Id (Transformations=I)
            if (nm.size() <= tag.size() || nm.compare(nm.size() - tag.size(), tag.size(), tag) != 0) continue;
            const std::string var = nm.substr(0, nm.size() - tag.size());
            TH1* hs = (TH1*)iv->Get(nm.c_str());
            TH1* hb = (TH1*)iv->Get((var + "__Background_Id").c_str());
            if (!hs || !hb) continue;
            DrawPair(c, hs, hb, Form("%s/input_variables/%d_%s.pdf", tdir.c_str(), ++n, var.c_str()), var + ", signal vs background, TMVA training sample (unit area), no cuts;" + var + ";fraction of events");
        }
    }
}
}

void F10_BDTCuts(const char* tmvaFile = "resultTMVAtraining_SigBkg.root"){
    const std::string dir = "plots/10_BDTcuts";
    gSystem->mkdir(dir.c_str(), true);
    gStyle->SetOptStat(0);
    TCanvas c("cBdtCuts", "", 900, 650);
    c.SetRightMargin(0.04);
    DrawBdtCuts(c, dir);
    DrawTmva(c, tmvaFile, dir);
    printf(">>> F10_BDTCuts: plots/10_BDTcuts (1_P_vs_BDTcut, 2_signalEff_vs_bkgEff, 3_cut_metrics, 4_FOM_vs_BDTcut, bdt_cut_table.csv) and plots/10_BDTcuts/TMVA\n");
}
