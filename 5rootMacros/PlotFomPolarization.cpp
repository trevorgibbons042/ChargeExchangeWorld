#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMath.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <TROOT.h>
#include "FitCore.h"
#include "CutFiles.h"
#include "FomCommon.h"

// 2D maps over theta ranges [lo, hi): x = lower bound, y = upper bound (hi > lo always, so only the upper triangle is filled). True KE or detector KE
// (energy argument, trees with the _DE suffix), one set of
// maps for the one cut, ke_gt80 (no KE cut, KeCuts in CutFiles.h). The FOM is only in PlotCEPerBin.cpp (cumulative FOM plot).
//   P maps, one per theta type (AP: the model equation at the CE t, same for both), plots/5_fomCalc/<1_reconTheta, 2_reconBlockTheta>/<ke>/ (recon block theta: block residual cut), files tagged apCE_<CE or recon> (recon protons: their
//   phi is fit with A_eff = the mean A(recon theta) over all recon protons in the range, so the unpolarized non-CE protons pull P down): 5_P_fit_2d_<tag>.pdf and .csv, and 6_P_pct_error_2d_<tag>.pdf = |100 * (P_fit - PPValue) / PPValue| (PPValue from settings.csv),
//   log color scale from 0.01 to 100:
//     N0 * (1 + AP * P * cos(phi)) fit to the CE protons with theta in [lo, hi), AP fixed to their mean AP from the chosen source (same fit as FitSlice
//     in PlotPolarization.cpp). From cut_APandT.root: the theta set's tree (thetaRe, phiRe or thetaReD, phiReD) with the AP source tree as friend.
static const double kStep = 0.05;    // deg, bin width of both axes (0 to 5 deg: 100 x 100 bins, like the other 2D histograms)
static const double kThetaMaxFP = 5.0; // deg, same range as PlotPolarization.cpp
static const int kPhiBinsFP = 100;
static const double kLowEdge = 0.3;  // deg, lower edge of the second wide fit plot
static const size_t kMinNFP = 10;     // CE protons needed in a range for its fit

// Fits N0 * (1 + AP * P * cos(phi)) to a phi histogram (0 to 360 deg), AP fixed. false when the fit fails.
// Same binned Poisson likelihood as the TF1 "L" fit in PlotPolarization.cpp, solved directly because the maps need about 5000 fits each (Minuit took ~80 s).
// The bin-center cosines sum to 0 over 0 to 360 deg, so N0 = total / bins for any c = AP * P and the likelihood in c separates:
// maximize sum n_b ln(1 + c cos_b) with Newton steps; error of c = 1 / sqrt(sum n_b cos_b^2 / (1 + c cos_b)^2). P stays in [-1, 1] like the old fit's limits.
static bool FitRange(const std::array<double, kPhiBinsFP>& counts, double ap, double& p, double& pErr){
    static const std::array<double, kPhiBinsFP> cosB = []{
        std::array<double, kPhiBinsFP> c{};
        for (int b = 0; b < kPhiBinsFP; b++) c[b] = std::cos((b + 0.5) * 2 * TMath::Pi() / kPhiBinsFP);
        return c;
    }();
    auto curvature = [&](double c, double* grad){ // -d2/dc2 of the log likelihood, and its gradient in grad
        double g = 0, h = 0;
        for (int b = 0; b < kPhiBinsFP; b++){
            if (counts[b] == 0) continue;
            const double d = cosB[b] / (1 + c * cosB[b]);
            g += counts[b] * d; h += counts[b] * d * d;
        }
        if (grad) *grad = g;
        return h;
    };
    const double cMax = std::fabs(ap) * (1 - 1e-9);
    double c = 0;
    for (int it = 0; it < 100; it++){
        double g = 0;
        const double h = curvature(c, &g);
        if (h <= 0) return false;
        const double cNew = std::clamp(c + g / h, -cMax, cMax);
        const bool done = std::fabs(cNew - c) < 1e-12;
        c = cNew;
        if (done) break;
    }
    const double h = curvature(c, nullptr);
    if (h <= 0) return false;
    p = c / ap; pErr = 1 / std::sqrt(h) / std::fabs(ap);
    return true;
}

struct Map {std::string file, title; TH2D h;};
static Map MakeMap(const std::string& file, const std::string& title, int nSlices){
    Map m{file, title, TH2D(("hMap" + file).c_str(), "", nSlices, 0, kThetaMaxFP, nSlices, 0, kThetaMaxFP)};
    m.h.SetDirectory(nullptr);
    return m;
}



using PhiCounts = std::vector<std::array<double, kPhiBinsFP>>; // per theta slice, the phi histogram (0 to 360 deg)
// Adds one proton (theta, phi in deg) to counts; returns its theta slice, or -1 when theta is outside 0 to kThetaMaxFP.
static int BinPhi(PhiCounts& counts, double theta, double phi){
    if (theta < 0 || theta >= kThetaMaxFP) return -1;
    const int s = (int)(theta / kStep), pb = (int)(phi / 360.0 * kPhiBinsFP);
    if (pb >= 0 && pb < kPhiBinsFP) counts[s][pb]++;
    return s;
}

// Per theta slice, the phi histogram of the protons and the sum of their A. The model summed over the protons of a range is
// N0_tot (1 + A_eff P cos(phi)) with A_eff = sum A / N, the plain mean of A over the protons that land in the range, so a range needs only these sums.
struct PhiAP {
    PhiCounts phi; std::vector<double> apSum, apN;
    PhiAP(int nSlices) : phi(nSlices), apSum(nSlices, 0.0), apN(nSlices, 0.0) {}
};
// A of every proton is the one at the center of its theta slice (t from theta, APFromT), so apSum of a slice is N_slice * A_slice
static void AddProton(PhiAP& d, double theta, double phi, double apScale){
    const int s = BinPhi(d.phi, theta, phi);
    if (s >= 0){d.apSum[s] += APFromT(kBeamMomentum, TFromAngle(kBeamMomentum, (s + 0.5) * kStep * TMath::Pi() / 180), apScale); d.apN[s]++;}
}

void PlotFomPolarization(std::string energy = "fKEnergy"){ // energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the input trees and plot folder
    ROOT::EnableImplicitMT(NThreads()); // implicit MT for the file reads (NThreads: ANA_NTHREADS, FitCore.h)
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); // colorblind-friendly, many steps so 10% and 100% differ clearly
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    const std::string keLabel = (energy == "fDEnergy") ? "ZDC KE" : "true KE";
    const int nSlices = (int)std::round(kThetaMaxFP / kStep);
    const std::vector<KeCut> kes = {KeCuts(keLabel).front()}; // CutFiles.h, the one cut ke_gt80 (no KE cut, no residual cut folders)
    TCanvas c("cFomMaps", "", 1500, 1200);
    c.SetRightMargin(0.30); c.SetLeftMargin(0.12); c.SetBottomMargin(0.12);

    // P maps, one per theta type, AP source and KE cut
    const double pRef = ReadSettings().at("PPValue"), apScale = ReadSettings().at("APScale");
    TFile* ceFile = TFile::Open("cuts/cut_APandT.root", "READ");
    if (!ceFile){printf("PlotFomPolarization: cut_APandT.root missing\n"); return;}
    TFile* allFile = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree* all = allFile ? (TTree*)allFile->Get(kAngleTree.c_str()) : nullptr; // every fitted proton; the cut lists "Angles" + suffix hold those that pass the quality cut
    if (!all || !all->GetBranch("passQuality")){printf("PlotFomPolarization: AnglesRaw with passQuality in cut_angles_all.root missing\n"); return;}
    // block: the block position set, which only takes the block AP and the block residual cut (kBlockCut, QualityCut.h); the recon set takes AP from CE and from recon
    struct ThetaSet {std::string dir, name, tree, theta, phi, allTheta, allPhi; bool block;}; // all* : theta and phi (rad) of every proton with a reconstructed path in the Angles tree
    const std::vector<std::string> typeTags = {"CE", "recon"}, typeNames = {"CE recon protons", "all recon protons"};
    // the one AP source: the model equation (APFromT, FitCore.h) at the CE t (tCE) of each CE proton; theta and phi come from the theta set, rows line up between the trees
    struct ApSrc {std::string tag, name, tree, branch;};
    std::vector<ThetaSet> sets = {{"1_reconTheta", "recon theta", "CalcValsRe" + suffix, "thetaRe", "phiRe", "outTheta", "outPhi", false},
                                        {"2_reconBlockTheta", "recon block theta", "CalcValsReD" + suffix, "thetaReD", "phiReD", "outDTheta", "outDPhi", true}};
    if (FomBlockOnly()) sets.erase(sets.begin()); // the cut runs besides mega only have the block position set
    const std::vector<ApSrc> apSrcs = {{"apCE", "AP from CE t", "CE" + suffix, "tCE"}};
    // phi counts of protons with reconstructed paths per theta set and KE cut: they do not depend on the AP source, so one pass over the (large) all-protons tree fills them all
    std::vector<std::vector<PhiAP>> phiAll(sets.size(), std::vector<PhiAP>(kes.size(), PhiAP(nSlices))); // A of each proton from its recon theta (APFromT, FitCore.h)
    {
        std::vector<double> th(sets.size(), 0.0), ph(sets.size(), 0.0);
        // the entries of every KE cut (block set: with the block residual cut), CutTrees.cpp
        std::vector<std::vector<TEntryList*>> lists(sets.size());
        for (size_t s = 0; s < sets.size(); s++) for (const auto& ke : kes) lists[s].push_back(CutList(ke.dir, all, sets[s].block, "Angles" + suffix));
        all->SetBranchStatus("*", 0);
        for (size_t s = 0; s < sets.size(); s++){
            all->SetBranchStatus(sets[s].allTheta.c_str(), 1); all->SetBranchStatus(sets[s].allPhi.c_str(), 1);
            all->SetBranchAddress(sets[s].allTheta.c_str(), &th[s]); all->SetBranchAddress(sets[s].allPhi.c_str(), &ph[s]);
        }
        const double toDeg = 180 / TMath::Pi();
        for (Long64_t e = 0, n = all->GetEntries(); e < n; e++){
            all->GetEntry(e);
            for (size_t s = 0; s < sets.size(); s++){
                for (size_t k = 0; k < kes.size(); k++)
                    if (lists[s][k] && lists[s][k]->Contains(e)) AddProton(phiAll[s][k], th[s] * toDeg, ph[s] * toDeg, apScale);
            }
        }
    }
    for (size_t si = 0; si < sets.size(); si++) for (const auto& apSrc : apSrcs){
        const auto& set = sets[si];
                TTree* ce = (TTree*)ceFile->Get(set.tree.c_str());
        if (!ce || !ceFile->Get(apSrc.tree.c_str())){printf("PlotFomPolarization: tree %s or %s missing\n", set.tree.c_str(), apSrc.tree.c_str()); continue;}
        const bool apInOwnTree = apSrc.tree == set.tree;
        if (!apInOwnTree && !ce->GetFriend(apSrc.tree.c_str())) ce->AddFriend(apSrc.tree.c_str());
        ce->SetEstimate(ce->GetEntries() + 1);
        for (size_t ki = 0; ki < kes.size(); ki++){
            const auto& ke = kes[ki];
            ApplyCut(ce, ke.dir, set.block);
            // the CE protons take their own CE angles (thetaCE, phiCE, deg) from the CE tree, not the recon or recon block ones of the theta set
            const std::string ceCols = (apInOwnTree ? "" : apSrc.tree + ".");
            const Long64_t nCE = ce->Draw((ceCols + "thetaCE:" + ceCols + "phiCE:" + ceCols + apSrc.branch + ":" + ceCols + "apCE:" + set.theta + ":" + set.phi).c_str(), "", "goff");
            ce->SetEntryList(nullptr);
            // per theta slice, the phi counts of the CE protons and the AP sum of the CE protons (the AP of both fits, as in PlotPolarization.cpp)
            PhiAP phiCE(nSlices);
            for (Long64_t k = 0; k < nCE; k++) AddProton(phiCE, ce->GetVal(0)[k], ce->GetVal(1)[k], apScale);
            PhiAP phiCERecon(nSlices); // the same CE protons at their recon angles (set.theta, set.phi), for the all-theta phi fit
            for (Long64_t k = 0; k < nCE; k++) AddProton(phiCERecon, ce->GetVal(4)[k], ce->GetVal(5)[k], apScale);
            double apSimMean = 0; // mean of the sim A (LastAP, apCE) of the CE protons, the A of the CE phi fits
            for (Long64_t k = 0; k < nCE; k++) apSimMean += ce->GetVal(3)[k] / std::max<Long64_t>(nCE, 1);

            for (int type = 0; type < 2; type++){
                // file numbers continue PlotPolarization.cpp's (CE 2k-1, recon 2k): kinds 5, 6 in the recon theta AP folders, 3, 4 in the block KE folder (no total fits there)
                auto fileNum = [type](int kind){return std::to_string(2 * (kind - 1) + 1 + type) + "_";};
                const PhiAP& protons = type == 0 ? phiCE : phiAll[si][ki];
                const std::string dir = FomRoot() + "/" + set.dir + "/" + (set.block ? ke.dir + "/" : "") + kFitDir; // the position set has no KE cut folder // with the total and slice phi fits of PlotPolarization.cpp
                const std::string tag = typeTags[type];
                gSystem->mkdir(dir.c_str(), true);
                Map pMap = MakeMap(fileNum(set.block ? 3 : 5) + "P_fit_2d", "P_{fit}, fit N_{0} (1 + #LTAP#GT P cos#phi), #LTAP#GT fixed", nSlices);
                Map pctMap = MakeMap(fileNum(set.block ? 4 : 6) + "P_pct_error_2d", "|100 (P_{fit} - PPValue) / PPValue| (%)", nSlices);
                pctMap.h.SetMinimum(0); pctMap.h.SetMaximum(100); // linear color scale 0 to 100 %, larger errors draw in the top color
                FILE* out = fopen((dir + pMap.file + "_" + tag + ".csv").c_str(), "w");
                if (out) fprintf(out, "theta_lo,theta_hi,N,AP_mean,P,P_err,abs_pct_err,b,b_max\n");
                for (int i = 0; i < nSlices; i++){
                    std::array<double, kPhiBinsFP> counts{}; // slices i to j-1 summed, extended by one slice per step of j
                    double apS = 0, apCount = 0, nPhi = 0;
                    for (int j = i + 1; j <= nSlices; j++){
                        for (int b = 0; b < kPhiBinsFP; b++){counts[b] += protons.phi[j - 1][b]; nPhi += protons.phi[j - 1][b];}
                        apS += protons.apSum[j - 1]; apCount += protons.apN[j - 1];
                        double p = 0, pErr = 0;
                        if (apCount == 0 || nPhi < kMinNFP || std::fabs(apS / apCount) < 1e-12 || !FitRange(counts, apS / apCount, p, pErr)) continue;
                        const double absPct = std::fabs(100.0 * (p - pRef) / pRef);
                        pMap.h.SetBinContent(i + 1, j, p);
                        pctMap.h.SetBinContent(i + 1, j, absPct);
                        if (out) fprintf(out, "%g,%g,%.0f,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n", i * kStep, j * kStep, nPhi, apS / apCount, p, pErr, absPct, p * apS / kPhiBinsFP, apS / kPhiBinsFP); // b = P sum(N A) per phi bin, |b| <= b_max
                    }
                }
                if (out) fclose(out);
                for (const int s0 : {0}){ // phi histogram of the range s0 * kStep to kThetaMaxFP (all thetas, then from kLowEdge) (20 bins) with the fit of FitRange (same fit as the from0p3 point at theta_max = kThetaMaxFP) and the input P
                    std::array<double, kPhiBinsFP> counts{};
                    double apS = 0, apCount = 0, nPhi = 0;
                    const PhiAP& fitProtons = type == 0 ? phiCERecon : protons; // CE protons with a recon track at their recon angles, recon protons at theirs
                    for (int s = s0; s < nSlices; s++){
                        for (int b = 0; b < kPhiBinsFP; b++){counts[b] += fitProtons.phi[s][b]; nPhi += fitProtons.phi[s][b];}
                        apS += fitProtons.apSum[s]; apCount += fitProtons.apN[s];
                    }
                    double p = 0, pErr = 0; (void)p; (void)pErr;
                    if (apCount > 0 && nPhi >= kMinNFP && std::fabs(apS / apCount) > 1e-12){
                        // N(phi) = N0 (1 + P <A> cos(phi)) per bin, N0 = entries / bins and <A> fixed, only P is fitted. CE protons take the sim A (mean apCE), recon protons the A from t(theta) at the slice centers
                        const double aFix = type == 0 ? apSimMean : apS / apCount;
                        const int nFitBins = kPhiBinsFP / 2; // 50 bins
                        TH1D h("hPhiRange", "", nFitBins, 0, 360);
                        for (int b = 0; b < kPhiBinsFP; b++) h.AddBinContent(b / 2 + 1, counts[b]);
                        h.SetEntries(nPhi);
                        h.SetTitle(((type == 0 ? std::string("CE recon protons only (not all CE protons)") : typeNames[type]) + ", " + set.name + " angles" + ", phi fit over all theta" + ";Angle (degrees);Count").c_str());
                        h.GetXaxis()->SetNdivisions(505);
                        TF1 fFit("fPhiRange", "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                        fFit.SetParNames("N0", "P", "A");
                        fFit.FixParameter(0, nPhi / nFitBins); fFit.FixParameter(2, aFix); fFit.SetParameter(1, 0.5); fFit.SetParLimits(1, -1, 1);
                        gStyle->SetOptStat(1111); gStyle->SetOptFit(1111);
                        c.Clear(); c.SetLogz(0); c.SetLogy(0);
                        c.SetRightMargin(0.25); c.SetTopMargin(0.17);
                        TH1D* hd = (TH1D*)h.DrawCopy();
                        hd->SetStats(true); // the histogram was built while the stat box was off for the maps
                        hd->Fit(&fFit, "Q");
                        hd->SetMinimum(0); hd->SetMaximum(1.1 * hd->GetMaximum());
                        hd->SetName("Fit results"); // the first line of the stats box is the histogram name
                        c.Update();
                        if (auto* st = (TPaveStats*)hd->FindObject("stats")){st->SetX1NDC(0.76); st->SetX2NDC(0.995); st->SetY1NDC(0.55); st->SetY2NDC(0.83);}
                        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.028);
                        tx.DrawLatex(0.13, 0.90, "Fit  N(#phi) = N_{0} (1 + P #upoint #LTA#GT cos#phi),   N_{0} = Entries / bins (fixed)");
                        tx.DrawLatex(0.13, 0.865, Form("#LTA#GT = %.4f (fixed, %s),   P = %.2f #pm %.2f", aFix, type == 0 ? "mean sim A (LastAP)" : "mean A from t(#theta)", fFit.GetParameter(1), fFit.GetParError(1)));
                        c.Modified(); c.Update();
                        c.SaveAs((dir + fileNum(set.block ? 7 : 9) + "phi_fit_all_theta_" + tag + ".pdf").c_str());
                        gStyle->SetOptStat(0); c.SetRightMargin(0.30); c.SetTopMargin(0.1);
                    }
                }
            }
        }
    }
}
