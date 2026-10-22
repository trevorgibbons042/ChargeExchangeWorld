#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "FitCore.h"

// 2D maps over theta ranges [lo, hi): x = lower bound, y = upper bound (hi > lo always, so only the upper triangle is filled). True KE or detector KE
// (energy argument, folders 1_KEtrue / 2_KEDetector, trees with the _DE suffix), one set of
// maps per KE cut (ke_all, ke_gt50 ... ke_gt90, the same cuts as PlotPolarization.cpp and PlotAngle.cpp).
//   FOM maps, one set per theta type, plots/6_fomCalc/1_KEtrue/FOM_calc/<1_reconTheta, 2_reconBlockTheta>/3_2drangemaps/<1_allprotons, 2_ceprotons>/<ke>/6_fom_2d.pdf (2_ceprotons
//   has only that one, 1_allprotons also 7_pct_CE_2d.pdf, 8_fom_x_CE_2d.pdf), values in 1_csv/fom_2d_<protons>_<ke>.csv:
//     FOM = sum of eff_x_AP2 over the fom.csv bins whose center is in [lo, hi) (1_csv/<ke>.csv, PlotAngle.cpp), % CE = N_CE / N of those bins,
//     FOM x CE = FOM times N_CE / N, the last point of the eff_x_AP2_x_CE curve in PlotFomCumulative.cpp for that range.
//   P maps, one per theta type and AP source, plots/6_fomCalc/1_KEtrue/3_polarization/<1_reconTheta, 2_reconBlockTheta>/total_fits/3_fitranges/<1_apCE, 2_apRecon,
//   3_apReconBlock>/<1_ce_protons, 2_all_protons>/<ke>/ (all protons: their phi is fit with the same AP, so the unpolarized non-CE protons pull P down): P_fit_2d.pdf and .csv, and P_pct_error_2d.pdf = |100 * (P_fit - PPValue) / PPValue| (PPValue from settings.csv),
//   log color scale from 0.01 to 100:
//     N0 * (1 + AP * P * cos(phi)) fit to the CE protons with theta in [lo, hi), AP fixed to their mean AP from the chosen source (same fit as FitSlice
//     in PlotPolarization.cpp). From cut_APandT.root: the theta set's tree (thetaRe, phiRe or thetaReD, phiReD) with the AP source tree as friend.
static const double kStep = 0.05;    // deg, bin width of both axes (0 to 5 deg: 100 x 100 bins, like the other 2D histograms)
static const double kThetaMaxFP = 5.0; // deg, same range as PlotPolarization.cpp
static const int kPhiBinsFP = 100;
static const size_t kMinNFP = 10;     // CE protons needed in a range for its fit

// Fits N0 * (1 + AP * P * cos(phi)) to a phi histogram (0 to 360 deg), AP fixed. false when the fit fails.
static bool FitRange(const std::array<double, kPhiBinsFP>& counts, double ap, double& p, double& pErr){
    TH1D h("hFomPolPhi", "", kPhiBinsFP, 0, 360);
    h.SetDirectory(nullptr);
    double n = 0;
    for (int b = 0; b < kPhiBinsFP; b++){h.SetBinContent(b + 1, counts[b]); n += counts[b];}
    h.SetEntries(n);
    TF1 f("fFomPol", "[0]*(1+[1]*[2]*cos(x*TMath::DegToRad()))", 0, 360);
    f.SetParameters(n / kPhiBinsFP, 0.5, ap);
    f.FixParameter(2, ap);
    f.SetParLimits(1, -1, 1);
    if (h.Fit(&f, "LQ0") != 0) return false;
    p = f.GetParameter(1); pErr = f.GetParError(1);
    return true;
}

struct Map {std::string file, title; TH2D h;};
static Map MakeMap(const std::string& file, const std::string& title, int nSlices){
    Map m{file, title, TH2D(("hMap" + file).c_str(), "", nSlices, 0, kThetaMaxFP, nSlices, 0, kThetaMaxFP)};
    m.h.SetDirectory(nullptr);
    return m;
}

// Draws the map with its values outside the graph (box in the right margin) and the marked cells on it, and writes the stats.csv row. Normally the min and max
// cells are marked; with a reference value (ref, not NaN) only the cell closest to it is, for the P_fit map against the true polarization.
static void SaveMap(TCanvas& c, Map& m, const std::string& path, const std::string& what, bool logz = false, double ref = std::nan("")){
    m.h.SetTitle((m.title + " over recon #theta range, " + what + ";lower #theta bound (deg);upper #theta bound (deg)").c_str());
    c.SetLogz(logz);
    const bool hasRef = !std::isnan(ref);
    int nCells = 0, minX = 0, minY = 0, maxX = 0, maxY = 0, clX = 0, clY = 0;
    double mn = 0, mx = 0, clV = 0, clD = 0;
    for (int bx = 1; bx <= m.h.GetNbinsX(); bx++) for (int by = 1; by <= m.h.GetNbinsY(); by++){
        const double v = m.h.GetBinContent(bx, by);
        if (v == 0) continue;
        if (nCells == 0 || v < mn){mn = v; minX = bx; minY = by;}
        if (nCells == 0 || v > mx){mx = v; maxX = bx; maxY = by;}
        if (hasRef && (nCells == 0 || std::fabs(v - ref) < clD)){clD = std::fabs(v - ref); clV = v; clX = bx; clY = by;}
        nCells++;
    }
    if (mn < 0){ // empty cells would draw in the color of 0 once the scale reaches below it: put them under the scale so they stay blank
        for (int bx = 1; bx <= m.h.GetNbinsX(); bx++) for (int by = 1; by <= m.h.GetNbinsY(); by++)
            if (m.h.GetBinContent(bx, by) == 0) m.h.SetBinContent(bx, by, mn - 1.0);
        m.h.SetMinimum(mn); m.h.SetMaximum(mx);
    }
    m.h.Draw("COLZ");
    auto lo = [&](int bx){return m.h.GetXaxis()->GetBinLowEdge(bx);};
    auto hi = [&](int by){return m.h.GetYaxis()->GetBinUpEdge(by);};
    TPaveText box(0.815, 0.68, 0.995, 0.86, "NDC"); // below the title, right of the palette labels
    box.SetBorderSize(1); box.SetFillColor(0); box.SetTextAlign(12); box.SetTextSize(0.022);
    box.AddText(m.file.c_str());
    box.AddText(Form("filled cells = %d", nCells));
    if (hasRef){
        box.AddText(Form("true value = %.4g", ref));
        box.AddText(Form("closest = %.4g", clV));
        if (nCells > 0) box.AddText(Form("at %.2f to %.2f deg", lo(clX), hi(clY)));
    } else {
        box.AddText(Form("min = %.4g", mn));
        box.AddText(Form("max = %.4g", mx));
        if (nCells > 0){
            box.AddText(Form("min at %.2f to %.2f deg", lo(minX), hi(minY)));
            box.AddText(Form("max at %.2f to %.2f deg", lo(maxX), hi(maxY)));
        }
    }
    box.DrawClone();
    // the marked cells on the graph itself: an outline around the cell (so its own color stays visible) and its value beside it (flipped to the left near the right edge)
    auto mark = [&](int bx, int by, const char* tag, double v){
        const double x = m.h.GetXaxis()->GetBinCenter(bx), y = m.h.GetYaxis()->GetBinCenter(by);
        const double x1 = m.h.GetXaxis()->GetBinLowEdge(bx), x2 = m.h.GetXaxis()->GetBinUpEdge(bx), y1 = m.h.GetYaxis()->GetBinLowEdge(by), y2 = m.h.GetYaxis()->GetBinUpEdge(by);
        const double pad = kStep; // one cell of room around it, so the marked cell itself is not covered by the line
        TBox outer(x1 - pad, y1 - pad, x2 + pad, y2 + pad); outer.SetFillStyle(0); outer.SetLineColor(kBlack); outer.SetLineWidth(3); outer.DrawClone();
        TBox inner(x1 - pad, y1 - pad, x2 + pad, y2 + pad); inner.SetFillStyle(0); inner.SetLineColor(kWhite); inner.SetLineWidth(1); inner.DrawClone();
        const bool left = x > 0.7 * kThetaMaxFP;
        const double dx = 0.12, w = 1.4, h = 0.3;
        TPaveText label(left ? x - dx - w : x + dx, y - h / 2, left ? x - dx : x + dx + w, y + h / 2);
        label.SetBorderSize(1); label.SetFillColor(kWhite); label.SetTextSize(0.022); label.SetTextColor(kBlack);
        label.AddText(Form("%s %.4g", tag, v));
        label.DrawClone();
    };
    if (nCells > 0){
        if (hasRef) mark(clX, clY, "closest", clV);
        else {mark(minX, minY, "min", mn); mark(maxX, maxY, "max", mx);}
    }
    // stats.csv next to the plot: the theta range (lower, upper bound) of the min and max, or of the cell closest to the reference; the columns that do not apply stay empty
    static std::set<std::string> statsStarted;
    if (nCells > 0){
        const std::string name = gSystem->BaseName(path.c_str());
        WriteStatsRow(path, "plot,filled_cells,min,min_lo_deg,min_hi_deg,max,max_lo_deg,max_hi_deg,true_value,closest,closest_lo_deg,closest_hi_deg",
                      hasRef ? Form("%s,%d,,,,,,,%.8g,%.8g,%.2f,%.2f", name.c_str(), nCells, ref, clV, lo(clX), hi(clY))
                             : Form("%s,%d,%.8g,%.2f,%.2f,%.8g,%.2f,%.2f,,,,", name.c_str(), nCells, mn, lo(minX), hi(minY), mx, lo(maxX), hi(maxY)), statsStarted);
    }
    c.SaveAs(path.c_str());
    c.SetLogz(0);
}

struct KECut {std::string dir, cut, label;};

void PlotFomPolarization(std::string energy = "fKEnergy"){ // energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the input trees and plot folder
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); // colorblind-friendly, many steps so 10% and 100% differ clearly
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    const std::string energyDir = (energy == "fDEnergy") ? "2_KEDetector" : "1_KEtrue";
    const std::string keLabel = (energy == "fDEnergy") ? "ZDC KE" : "true KE";
    const int nSlices = (int)std::round(kThetaMaxFP / kStep);
    const double initialKE = 25000; // MeV
    std::vector<KECut> kes = {{"ke_all", "", "no KE cut"}};
    for (int pct : {50, 60, 70, 80, 90}) kes.push_back({"ke_gt" + std::to_string(pct), "finalKE > " + std::to_string((int)(initialKE * pct / 100.0)), "KE cut " + std::to_string(pct) + "%"});
    TCanvas c("cFomMaps", "", 1500, 1200);
    c.SetRightMargin(0.30); c.SetLeftMargin(0.12); c.SetBottomMargin(0.12);

    // FOM maps: per fom.csv bin, its center, eff_x_AP2 of all protons (the csv column), the same for CE protons only (N_CE / total neutrons times
    // mean_AP^2) and N, N_CE
    const double totalNeutrons = ReadSettings().at("TotalNeutrons");
    struct FomBin {double center, effAP2All, effAP2CE, n, nCE;};
    for (const std::string& setDir : {"1_reconTheta", "2_reconBlockTheta"}){ // the FOM tables of both theta types, PlotAngle.cpp
        const std::string fomDir = "plots/6_fomCalc/" + energyDir + "/FOM_calc/" + setDir + "/";
        const std::string posLabel = setDir == "1_reconTheta" ? "recon theta" : "recon block theta";
        gSystem->mkdir((fomDir + "1_csv").c_str(), true);
        for (const auto& ke : kes){
            std::vector<FomBin> fomBins;
            std::ifstream in(fomDir + "1_csv/" + ke.dir + ".csv");
            if (!in){printf("PlotFomPolarization: %s1_csv/%s.csv missing (run PlotAngle.cpp first)\n", fomDir.c_str(), ke.dir.c_str()); continue;}
            std::string line; std::getline(in, line);
            while (std::getline(in, line)){
                std::vector<double> v; std::stringstream ss(line); std::string f;
                while (std::getline(ss, f, ',')) v.push_back(std::stod(f));
                if (v.size() < 8 || v[1] - v[0] > 0.5) continue; // skip the total row
                fomBins.push_back({0.5 * (v[0] + v[1]), v[7], v[3] / totalNeutrons * v[5] * v[5], v[2], v[3]});
            }
            // 1_allprotons: the csv FOM (N / total neutrons), with the % CE and CE-weighted maps; 2_ceprotons: the FOM of the CE protons only (N_CE / total neutrons),
            // where % CE and the CE weight would be 100% / 1, so only the FOM map
            for (const bool ceOnly : {false, true}){
                const std::string typeDir = ceOnly ? "2_ceprotons" : "1_allprotons";
                std::vector<Map> fomMaps = {MakeMap("6_fom_2d", ceOnly ? "FOM of CE protons = #varepsilon_{CE} A_{y}^{2}" : "FOM = #varepsilon A_{y}^{2}", nSlices)};
                if (!ceOnly){
                    fomMaps.push_back(MakeMap("7_pct_CE_2d", "CE protons (%)", nSlices));
                    fomMaps.push_back(MakeMap("8_fom_x_CE_2d", "FOM #times N_{CE} / N_{all protons}", nSlices));
                }
                const std::string mapDir = fomDir + "3_2drangemaps/" + typeDir + "/" + ke.dir + "/";
                gSystem->mkdir(mapDir.c_str(), true);
                FILE* fomOut = fopen((fomDir + "1_csv/fom_2d_" + typeDir + "_" + ke.dir + ".csv").c_str(), "w");
                if (fomOut) fprintf(fomOut, "theta_lo,theta_hi,N,N_CE,pct_CE,FOM,FOM_x_CE\n");
                for (int i = 0; i < nSlices; i++) for (int j = i + 1; j <= nSlices; j++){ // j > i: the upper bound is always above the lower
                    const double lo = i * kStep, hi = j * kStep;
                    double run = 0, nAll = 0, nCEbins = 0;
                    for (const auto& b : fomBins){
                        if (b.center < lo || b.center >= hi) continue;
                        run += ceOnly ? b.effAP2CE : b.effAP2All; nAll += b.n; nCEbins += b.nCE;
                    }
                    if (nAll == 0) continue;
                    const double pctCE = 100.0 * nCEbins / nAll, fomCE = run * nCEbins / nAll;
                    // x bin i + 1 holds lo, y bin j holds hi (bin j spans (j - 1, j) * kStep)
                    fomMaps[0].h.SetBinContent(i + 1, j, run);
                    if (!ceOnly){
                        fomMaps[1].h.SetBinContent(i + 1, j, pctCE);
                        fomMaps[2].h.SetBinContent(i + 1, j, fomCE);
                    }
                    if (fomOut) fprintf(fomOut, "%g,%g,%.0f,%.0f,%.6g,%.8g,%.8g\n", lo, hi, nAll, nCEbins, pctCE, run, ceOnly ? 0.0 : fomCE);
                }
                if (fomOut) fclose(fomOut);
                for (auto& m : fomMaps) SaveMap(c, m, mapDir + m.file + ".pdf", keLabel + ", " + posLabel + ", " + ke.label);
            }
        }
    }

    // P maps, one per theta type, AP source and KE cut
    const double pRef = ReadSettings().at("PPValue");
    TFile* ceFile = TFile::Open("cut_APandT.root", "READ");
    if (!ceFile){printf("PlotFomPolarization: cut_APandT.root missing\n"); return;}
    TFile* allFile = TFile::Open("cut_angles_all.root", "READ");
    TTree* all = allFile ? (TTree*)allFile->Get(("Angles" + suffix).c_str()) : nullptr;
    if (!all){printf("PlotFomPolarization: Angles tree in cut_angles_all.root missing\n"); return;}
    all->SetEstimate(all->GetEntries() + 1);
    struct ThetaSet {std::string dir, name, tree, theta, phi, allTheta, allPhi;}; // all* : theta and phi (rad) of every proton in the Angles tree
    const std::vector<std::string> typeDirs = {"1_ce_protons", "2_all_protons"}, typeNames = {"CE protons", "all protons"};
    struct ApSrc {std::string dir, name, tree, branch;}; // the AP of each proton; theta and phi come from the theta set, rows line up between the trees
    const std::vector<ThetaSet> sets = {{"1_reconTheta", "recon theta", "CalcValsRe" + suffix, "thetaRe", "phiRe", "outTheta", "outPhi"},
                                        {"2_reconBlockTheta", "recon block theta", "CalcValsReD" + suffix, "thetaReD", "phiReD", "outDTheta", "outDPhi"}};
    const std::vector<ApSrc> apSrcs = {{"1_apCE", "AP from CE", "CE" + suffix, "apCE"},
                                       {"2_apRecon", "AP from recon", "CalcValsRe" + suffix, "apRe"},
                                       {"3_apReconBlock", "AP from recon block", "CalcValsReD" + suffix, "apReD"}};
    for (const auto& set : sets) for (const auto& apSrc : apSrcs){
        TTree* ce = (TTree*)ceFile->Get(set.tree.c_str());
        if (!ce || !ceFile->Get(apSrc.tree.c_str())){printf("PlotFomPolarization: tree %s or %s missing\n", set.tree.c_str(), apSrc.tree.c_str()); continue;}
        const bool apInOwnTree = apSrc.tree == set.tree;
        if (!apInOwnTree && !ce->GetFriend(apSrc.tree.c_str())) ce->AddFriend(apSrc.tree.c_str());
        ce->SetEstimate(ce->GetEntries() + 1);
        for (const auto& ke : kes){
            const Long64_t nCE = ce->Draw((set.theta + ":" + set.phi + ":" + (apInOwnTree ? "" : apSrc.tree + ".") + apSrc.branch).c_str(), ke.cut.c_str(), "goff");
            // per theta slice, the phi counts of the CE protons and of all protons, and the AP sum of the CE protons (the AP of both fits, as in PlotPolarization.cpp)
            std::vector<std::array<double, kPhiBinsFP>> phiCE(nSlices), phiAll(nSlices);
            std::vector<double> apSum(nSlices, 0.0), apN(nSlices, 0.0);
            auto binPhi = [&](std::vector<std::array<double, kPhiBinsFP>>& counts, double theta, double phi){
                if (theta < 0 || theta >= kThetaMaxFP) return -1;
                const int s = (int)(theta / kStep), pb = (int)(phi / 360.0 * kPhiBinsFP);
                if (pb >= 0 && pb < kPhiBinsFP) counts[s][pb]++;
                return s;
            };
            for (Long64_t k = 0; k < nCE; k++){
                const int s = binPhi(phiCE, ce->GetVal(0)[k], ce->GetVal(1)[k]);
                if (s >= 0){apSum[s] += ce->GetVal(2)[k]; apN[s]++;}
            }
            const std::string deg = "*180/TMath::Pi()";
            const Long64_t nAll = all->Draw((set.allTheta + deg + ":" + set.allPhi + deg).c_str(), ke.cut.c_str(), "goff");
            for (Long64_t k = 0; k < nAll; k++) binPhi(phiAll, all->GetVal(0)[k], all->GetVal(1)[k]);

            for (int type = 0; type < 2; type++){
                const auto& phiCounts = type == 0 ? phiCE : phiAll;
                const std::string dir = "plots/6_fomCalc/" + energyDir + "/3_polarization/" + set.dir + "/total_fits/3_fitranges/" + apSrc.dir + "/" + typeDirs[type] + "/" + ke.dir + "/";
                gSystem->mkdir(dir.c_str(), true);
                Map pMap = MakeMap("P_fit_2d", "P_{fit}", nSlices);
                Map pctMap = MakeMap("P_pct_error_2d", "|P_{fit} percent error vs PPValue| (%)", nSlices);
                pctMap.h.SetMinimum(0.01); pctMap.h.SetMaximum(100); // log color scale, so the floor is above 0
                FILE* out = fopen((dir + "P_fit_2d.csv").c_str(), "w");
                if (out) fprintf(out, "theta_lo,theta_hi,N,AP_mean,P,P_err,abs_pct_err\n");
                for (int i = 0; i < nSlices; i++){
                    std::array<double, kPhiBinsFP> counts{}; // slices i to j-1 summed, extended by one slice per step of j
                    double apS = 0, apCount = 0, nPhi = 0;
                    for (int j = i + 1; j <= nSlices; j++){
                        for (int b = 0; b < kPhiBinsFP; b++){counts[b] += phiCounts[j - 1][b]; nPhi += phiCounts[j - 1][b];}
                        apS += apSum[j - 1]; apCount += apN[j - 1];
                        double p = 0, pErr = 0;
                        if (apCount == 0 || nPhi < kMinNFP || std::fabs(apS / apCount) < 1e-12 || !FitRange(counts, apS / apCount, p, pErr)) continue;
                        const double absPct = std::fabs(100.0 * (p - pRef) / pRef);
                        pMap.h.SetBinContent(i + 1, j, p);
                        pctMap.h.SetBinContent(i + 1, j, absPct);
                        if (out) fprintf(out, "%g,%g,%.0f,%.8g,%.8g,%.8g,%.8g\n", i * kStep, j * kStep, nPhi, apS / apCount, p, pErr, absPct);
                    }
                }
                if (out) fclose(out);
                for (Map* m : {&pMap, &pctMap}) SaveMap(c, *m, dir + m->file + ".pdf", typeNames[type] + ", " + apSrc.name + ", " + set.name + ", " + keLabel + ", " + ke.label, m == &pctMap, m == &pMap ? pRef : std::nan(""));
            }
        }
    }
}
