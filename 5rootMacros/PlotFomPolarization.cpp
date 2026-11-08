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
#include "TreeColumns.h"

// 2D maps over theta ranges [lo, hi): x = lower bound, y = upper bound (hi > lo always, so only the upper triangle is filled). True KE or detector KE
// (energy argument, trees with the _DE suffix), one set of
// maps for the one cut, ke_gt80 (no KE cut, KeCuts in CutFiles.h). The FOM is only in PlotCEPerBin.cpp (cumulative FOM plot).
//   P maps, one per theta type (AP: the model equation at the CE t, same for both), plots/5_fomCalc/<1_reconTheta_estneu, 2_reconBlockTheta_estneu>/<ke>/, files tagged apCE_<CE or recon> (recon protons: their
//   phi is fit with A_eff = the mean A(RECONTHETA_ESTNEU) over all recon protons in the range, so the unpolarized non-CE protons pull P down): 5_P_fit_2d_<tag>.pdf and .csv, and 6_P_pct_error_2d_<tag>.pdf = |100 * (P_fit - PPValue) / PPValue| (PPValue from settings.csv),
//   log color scale from 0.01 to 100:
//     N0 * (1 + AP * P * cos(phi)) fit to the CE protons with theta in [lo, hi), AP fixed to their mean AP from the chosen source (same fit as FitSlice
//     in PlotPolarization.cpp). From cut_APandT.root: the theta set's tree (thetaRe, phiRe or thetaReD, phiReD) with the AP source tree as friend.
static const double kStep = 0.05;    // deg, bin width of both axes (0 to 5 deg: 100 x 100 bins, like the other 2D histograms)
static const double kThetaMaxFP = 5.0; // deg, same range as PlotPolarization.cpp
static const int kPhiBinsFP = 100;
static const double kLowEdge = 0.3;  // deg, lower edge of the second wide fit plot
static const size_t kMinNFP = 10;     // CE protons needed in a range for its fit

// Unbinned fit of N(phi) ~ 1 + AP * P * cos(phi) to the protons of a range, AP fixed: from the sums of cos^k(phi) of the protons (UnbinnedPFromMoments, FitCore.h), so a range of slices needs only the sums
// of its slices and the maps can fit their ~5000 ranges. P stays in [-1, 1]. block: the phi range of the block fits (kBlockPhiLo to kBlockPhiHi). false when the fit fails.
static bool FitRange(const double* mom, double n, const AMean& ap, bool block, double& p, double& pErr){
    const PFit r = UnbinnedPFromMoments(mom, n, ap, block ? kBlockPhiLo : 0, block ? kBlockPhiHi : 360);
    p = r.p; pErr = r.err;
    return r.ok && r.err > 0;
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
    PhiCounts phi; std::vector<double> apSum, apN, apSum2; // apSum2: the sums of A^2, for the spread of A (AMean)
    std::vector<double> apLo, apHi, dSum, dN; // per theta slice, the sums of A at thetaMin and thetaMax of the est neutron theta range, and of sin(h) / h with its count (PhiRange): the errors by source
    double totalApLo = 0, totalApHi = 0, totalDSum = 0, totalDN = 0; // the same for every proton at any theta
    std::vector<double> apLo2, apHi2, dSum2, dN2; double totalApLo2 = 0, totalApHi2 = 0, totalDSum2 = 0, totalDN2 = 0; // the same for the block grid range alone (the second range of AMean and PhiRange, block set only)
    std::vector<std::array<double, kPhiMoments>> mom; // per theta slice, the sums of cos^k(phi) of its protons (the unbinned fit, AddCosMoments)
    std::array<double, kPhiBinsFP> total{}; double totalApSum = 0, totalApSum2 = 0, totalN = 0; // every proton at any theta (no theta cut), for the all-theta phi fits
    std::array<double, kPhiMoments> totalMom{};
    PhiAP(int nSlices = 0) : phi(nSlices), apSum(nSlices, 0.0), apN(nSlices, 0.0), apSum2(nSlices, 0.0), apLo(nSlices, 0.0), apHi(nSlices, 0.0), dSum(nSlices, 0.0), dN(nSlices, 0.0), apLo2(nSlices, 0.0), apHi2(nSlices, 0.0), dSum2(nSlices, 0.0), dN2(nSlices, 0.0), mom(nSlices) {for (auto& m : mom) m.fill(0.0);}
};
// A of every proton is the one at the center of its theta slice (t from theta, APFromT), so apSum of a slice is N_slice * A_slice
// thMin, thMax (deg) and phHalf (rad): the est neutron range of the proton, NaN when it has none (the true CE angles): then A at both ends is the central A and the phi range is left out.
static void AddProton(PhiAP& d, double theta, double phi, double apScale, bool block, double thMin = NAN, double thMax = NAN, double phHalf = NAN, double thMin2 = NAN, double thMax2 = NAN, double phHalf2 = NAN){
    if (block && !InBlockPhi(phi)) return; // block phi: the spike at 0 is left out of every fit
    const int pt = (int)(phi / 360.0 * kPhiBinsFP);
    const double aAny = AEffOfTheta(theta, apScale); // the A of the proton at its own theta, for the total and for its theta slice
    const bool hasR = std::isfinite(thMin) && std::isfinite(thMax);
    const double aLo = hasR ? AEffOfTheta(thMin, apScale) : aAny, aHi = hasR ? AEffOfTheta(thMax, apScale) : aAny;
    PhiRange dr; if (std::isfinite(phHalf)) dr.Add(phHalf);
    const bool hasR2 = std::isfinite(thMin2) && std::isfinite(thMax2); // the block grid range alone
    const double aLo2 = hasR2 ? AEffOfTheta(thMin2, apScale) : aAny, aHi2 = hasR2 ? AEffOfTheta(thMax2, apScale) : aAny;
    PhiRange dr2; if (std::isfinite(phHalf2)) dr2.Add(phHalf2);
    if (pt >= 0 && pt < kPhiBinsFP){
        d.total[pt]++; d.totalN++; AddCosMoments(d.totalMom.data(), phi);
        d.totalApSum += aAny; d.totalApSum2 += aAny * aAny; d.totalApLo += aLo; d.totalApHi += aHi; d.totalDSum += dr.sum; d.totalDN += dr.n;
        d.totalApLo2 += aLo2; d.totalApHi2 += aHi2; d.totalDSum2 += dr2.sum; d.totalDN2 += dr2.n;
    }
    const int s = BinPhi(d.phi, theta, phi);
    if (s >= 0) AddCosMoments(d.mom[s].data(), phi);
    if (s >= 0){d.apSum[s] += aAny; d.apSum2[s] += aAny * aAny; d.apN[s]++; d.apLo[s] += aLo; d.apHi[s] += aHi; d.dSum[s] += dr.sum; d.dN[s] += dr.n; d.apLo2[s] += aLo2; d.apHi2[s] += aHi2; d.dSum2[s] += dr2.sum; d.dN2[s] += dr2.n;}
}
// Adds the counts and sums of s to d (same number of slices).
static void AddPhiAP(PhiAP& d, const PhiAP& s){
    for (size_t i = 0; i < d.phi.size(); i++){
        for (int b = 0; b < kPhiBinsFP; b++) d.phi[i][b] += s.phi[i][b];
        d.apSum[i] += s.apSum[i]; d.apN[i] += s.apN[i]; d.apSum2[i] += s.apSum2[i];
        d.apLo[i] += s.apLo[i]; d.apHi[i] += s.apHi[i]; d.dSum[i] += s.dSum[i]; d.dN[i] += s.dN[i];
        d.apLo2[i] += s.apLo2[i]; d.apHi2[i] += s.apHi2[i]; d.dSum2[i] += s.dSum2[i]; d.dN2[i] += s.dN2[i];
        for (int k = 0; k < kPhiMoments; k++) d.mom[i][k] += s.mom[i][k];
    }
    for (int k = 0; k < kPhiMoments; k++) d.totalMom[k] += s.totalMom[k];
    for (int b = 0; b < kPhiBinsFP; b++) d.total[b] += s.total[b];
    d.totalApSum += s.totalApSum; d.totalApSum2 += s.totalApSum2; d.totalN += s.totalN;
    d.totalApLo += s.totalApLo; d.totalApHi += s.totalApHi; d.totalDSum += s.totalDSum; d.totalDN += s.totalDN;
    d.totalApLo2 += s.totalApLo2; d.totalApHi2 += s.totalApHi2; d.totalDSum2 += s.totalDSum2; d.totalDN2 += s.totalDN2;
}
// AddProton for every row of cols (theta, phi in rad, columns 0 and 1; toDeg converts; with 5 columns also thetaMin, thetaMax (rad) and the phi half width (rad), columns 2 to 4, with 8 the block grid range alone, columns 5 to 7) on all threads (MapRows, FitCore.h): each row range fills its own PhiAP, added into d in range order after.
// The counts are exact; the AP sums can differ from a serial loop in the last digits.
static void AddProtons(PhiAP& d, const std::vector<std::vector<double>>& cols, double toDeg, double apScale, bool block){
    const size_t nSlices = d.apSum.size();
    const auto parts = MapRows((Long64_t)cols[0].size(), [&](Long64_t lo, Long64_t hi){
        PhiAP part((int)nSlices);
        for (Long64_t e = lo; e < hi; e++){
            if (cols.size() >= 8) AddProton(part, cols[0][e] * toDeg, cols[1][e] * toDeg, apScale, block, cols[2][e] * toDeg, cols[3][e] * toDeg, cols[4][e], cols[5][e] * toDeg, cols[6][e] * toDeg, cols[7][e]);
            else if (cols.size() >= 5) AddProton(part, cols[0][e] * toDeg, cols[1][e] * toDeg, apScale, block, cols[2][e] * toDeg, cols[3][e] * toDeg, cols[4][e]);
            else AddProton(part, cols[0][e] * toDeg, cols[1][e] * toDeg, apScale, block);
        }
        return part;
    });
    for (const PhiAP& part : parts) AddPhiAP(d, part);
}

void PlotFomPolarization(std::string energy = "fKEnergy"){ // energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the input trees and plot folder
    FastJit();
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
    TTree* all = allFile ? (TTree*)allFile->Get(kAngleTree) : nullptr; // every fitted proton
    if (!all){printf("PlotFomPolarization: AnglesRaw in cut_angles_all.root missing\n"); return;}
    // block: the block position set, which only takes the block AP; the recon set takes AP from CE and from recon
    struct ThetaSet {std::string dir, name, tree, theta, phi, allTheta, allPhi, thMin, thMax, phHalf, allThMin, allThMax, allPhHalf, thMin2, thMax2, phHalf2, allThMin2, allThMax2, allPhHalf2; bool block;}; // ...2: the block grid range columns (block set only) // thMin, thMax (deg), phHalf (rad): the est neutron range columns of the CE tree; allThMin, allThMax, allPhHalf (rad) those of AnglesRaw // all* : theta and phi (rad) of every proton with a reconstructed path in the Angles tree
    // types: CE recon protons, all recon protons (pdg 2212) and all recon charged particles (every charged particle of AnglesRaw, any pdg)
    const std::vector<std::string> typeTags = {"CE", "recon", "reconParticles"}, typeNames = {"CE recon protons", "all recon protons", "all recon charged particles"};
    // the one AP source: the model equation (APFromT, FitCore.h) at the CE t (tCE) of each CE proton; theta and phi come from the theta set, rows line up between the trees
    struct ApSrc {std::string tag, name, tree, branch;};
    std::vector<ThetaSet> sets = {{"1_reconTheta_estneu", "RECONTHETA_ESTNEU", "CalcValsRe" + suffix, "thetaRe", "phiRe", "thetaRelN", "phiRelN", "thetaReMin", "thetaReMax", "phiReHalf", "thetaRelNMin", "thetaRelNMax", "phiRelNHalf", "", "", "", "", "", "", false},
                                        {"2_reconBlockTheta_estneu", "RECONTHETABL_ESTNEU", "CalcValsReD" + suffix, "thetaReD", "phiReD", "thetaRelND", "phiRelND", "thetaReDMin", "thetaReDMax", "phiReDHalf", "thetaRelNDMin", "thetaRelNDMax", "phiRelNDHalf", "thetaReDBMin", "thetaReDBMax", "phiReDBHalf", "thetaRelNDBMin", "thetaRelNDBMax", "phiRelNDBHalf", true}};
    if (FomBlockOnly()) sets.erase(sets.begin()); // the cut runs besides mega only have the block position set
    const std::vector<ApSrc> apSrcs = {{"apCE", "AP from CE t", "CE" + suffix, "tCE"}};
    // phi counts of protons with reconstructed paths per theta set and KE cut: they do not depend on the AP source, so one read of the (large) all-protons tree per theta set and KE cut fills them all
    std::vector<std::vector<PhiAP>> phiAll(sets.size(), std::vector<PhiAP>(kes.size(), PhiAP(nSlices))); // A of each proton from its RECONTHETA_ESTNEU (APFromT, FitCore.h)
    auto phiAllPart = phiAll; // the same for every recon charged particle (any pdg), not only the protons
    {
        const double toDeg = 180 / TMath::Pi();
        for (size_t s = 0; s < sets.size(); s++) for (size_t k = 0; k < kes.size(); k++){
            // no KE or residual cut here: every recon proton (pdg 2212) and every recon charged particle (the Angles view); one multithreaded read (ColumnsUnordered, TreeColumns.h) of their angles each, rows in any order (AddProton only adds counts and AP sums)
            const auto anglesPart = ColumnsUnordered(all, sets[s].block ? std::vector<std::string>{sets[s].allTheta, sets[s].allPhi, sets[s].allThMin, sets[s].allThMax, sets[s].allPhHalf, sets[s].allThMin2, sets[s].allThMax2, sets[s].allPhHalf2} : std::vector<std::string>{sets[s].allTheta, sets[s].allPhi, sets[s].allThMin, sets[s].allThMax, sets[s].allPhHalf}, AngleSelection("Angles"));
            AddProtons(phiAllPart[s][k], anglesPart, toDeg, apScale, sets[s].block);
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
            ce->SetEntryList(nullptr); // no cut: every CE proton with a recon track
            // the CE protons take their own CE angles (thetaCE, phiCE, deg) from the CE tree, not the recon or recon block ones of the theta set
            const std::string ceCols = (apInOwnTree ? "" : apSrc.tree + ".");
            std::vector<std::string> ceExprs = {ceCols + "thetaCE", ceCols + "phiCE", ceCols + apSrc.branch, ceCols + "apCE", set.theta, set.phi, set.thMin, set.thMax, set.phHalf};
            if (set.block) for (const std::string& e : {set.thMin2, set.thMax2, set.phHalf2}) ceExprs.push_back(e); // the block grid range alone
            const std::vector<std::vector<double>> ceV = Columns(ce, ceExprs); // CE angles, AP, then the recon angles and their est neutron range
            const Long64_t nCE = (Long64_t)ceV[0].size();
            ce->SetEntryList(nullptr);
            // per theta slice, the phi counts of the CE protons and the AP sum of the CE protons (the AP of both fits, as in PlotPolarization.cpp)
            PhiAP phiCE(nSlices);
            PhiAP phiCERecon(nSlices); // the same CE protons at their recon angles (set.theta, set.phi), for the all-theta phi fit
            for (Long64_t k = 0; k < nCE; k++){
                AddProton(phiCE, ceV[0][k], ceV[1][k], apScale, set.block);
                if (set.block) AddProton(phiCERecon, ceV[4][k], ceV[5][k], apScale, set.block, ceV[6][k], ceV[7][k], ceV[8][k], ceV[9][k], ceV[10][k], ceV[11][k]);
                else AddProton(phiCERecon, ceV[4][k], ceV[5][k], apScale, set.block, ceV[6][k], ceV[7][k], ceV[8][k]);
            }

            struct Ovl {TH1D* h = nullptr; double nPhi = 0; AMean aFix; bool ok = false; PFit uf; PhiRange pr, pr2;}; // uf: the unbinned fit of the range
            Ovl ovl[2][3]; // per range and type: the phi counts of the fit plot, collected here and drawn together after the type loop
            for (int type = 0; type < 3; type++){
                if (type == 1) continue; // no all recon protons (pdg 2212): only the CE recon protons and all recon charged particles
                // file numbers continue PlotPolarization.cpp's (CE 2k-1, recon 2k): kinds 5, 6 in the RECONTHETA_ESTNEU AP folders, 3, 4 in the block KE folder (no total fits there)
                // the all recon charged particles files come after those: 20 + kind
                auto fileNum = [type](int kind){return std::to_string(type == 2 ? 20 + kind : 2 * (kind - 1) + 1 + type) + "_";};
                const PhiAP& protons = type == 0 ? phiCE : type == 1 ? phiAll[si][ki] : phiAllPart[si][ki];
                const std::string dir = FomRoot() + "/" + set.dir + "/" + kFitDir; // // with the total and slice phi fits of PlotPolarization.cpp
                const std::string tag = typeTags[type];
                gSystem->mkdir(dir.c_str(), true);
                Map pMap = MakeMap(fileNum(set.block ? 3 : 5) + "P_fit_2d", "P_{fit}, unbinned fit 1 + A_{eff} P cos#phi, A_{eff} = #Sigma N_{i} AP(#theta_{i}) / N_{0} fixed", nSlices);
                Map pctMap = MakeMap(fileNum(set.block ? 4 : 6) + "P_pct_error_2d", "|100 (P_{fit} - PPValue) / PPValue| (%)", nSlices);
                pctMap.h.SetMinimum(0); pctMap.h.SetMaximum(100); // linear color scale 0 to 100 %, larger errors draw in the top color
                FILE* out = fopen((dir + pMap.file + "_" + tag + ".csv").c_str(), "w");
                if (out) fprintf(out, "theta_lo,theta_hi,N,AP_mean,P,P_err,abs_pct_err,b,b_max\n");
                for (int i = 0; i < nSlices; i++){
                    std::array<double, kPhiMoments> mom{}; // sums of cos^k(phi) of slices i to j-1, extended by one slice per step of j
                    AMean aR; // the A of slices i to j-1
                    double nPhi = 0;
                    for (int j = i + 1; j <= nSlices; j++){
                        for (int b = 0; b < kPhiBinsFP; b++) nPhi += protons.phi[j - 1][b];
                        for (int k = 0; k < kPhiMoments; k++) mom[k] += protons.mom[j - 1][k];
                        aR += AMean{protons.apN[j - 1], protons.apSum[j - 1], protons.apSum2[j - 1]};
                        double p = 0, pErr = 0;
                        if (aR.n == 0 || nPhi < kMinNFP || std::fabs(aR.Mean()) < 1e-12 || !FitRange(mom.data(), nPhi, aR, set.block, p, pErr)) continue;
                        const double absPct = std::fabs(100.0 * (p - pRef) / pRef);
                        pMap.h.SetBinContent(i + 1, j, p);
                        pctMap.h.SetBinContent(i + 1, j, absPct);
                        if (out) fprintf(out, "%g,%g,%.0f,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n", i * kStep, j * kStep, nPhi, aR.Mean(), p, pErr, absPct, p * aR.sum / kPhiBinsFP, aR.sum / kPhiBinsFP); // b = P sum(N A) per phi bin, |b| <= b_max
                    }
                }
                if (out) fclose(out);
                for (const int range : {0, 1}){ // range 0: all thetas, range 1: theta 0 to 1 deg; phi histogram (all thetas, then from kLowEdge) (20 bins) with the fit of FitRange (same fit as the from0p3 point at theta_max = kThetaMaxFP) and the input P
                    std::array<double, kPhiBinsFP> counts{};
                    std::array<double, kPhiMoments> mom{};
                    AMean aR; PhiRange prR, prR2; // the A of the protons of the range, and their est neutron phi range and block grid phi range
                    double nPhi = 0;
                    const PhiAP& fitProtons = type == 0 ? phiCERecon : protons; // CE protons with a recon track at their recon angles, recon protons at theirs
                    if (range == 0){counts = fitProtons.total; mom = fitProtons.totalMom; nPhi = fitProtons.totalN; aR = AMean{fitProtons.totalN, fitProtons.totalApSum, fitProtons.totalApSum2, fitProtons.totalApLo, fitProtons.totalApHi, fitProtons.totalApLo2, fitProtons.totalApHi2}; prR.n = fitProtons.totalDN; prR.sum = fitProtons.totalDSum; prR2.n = fitProtons.totalDN2; prR2.sum = fitProtons.totalDSum2;} // no theta cut
                    else for (int sl = 0; sl < std::min(nSlices, (int)std::lround(1.0 / kStep)); sl++){ // the slices of theta 0 to 1 deg
                        for (int b = 0; b < kPhiBinsFP; b++){counts[b] += fitProtons.phi[sl][b]; nPhi += fitProtons.phi[sl][b];}
                        for (int k = 0; k < kPhiMoments; k++) mom[k] += fitProtons.mom[sl][k];
                        aR += AMean{fitProtons.apN[sl], fitProtons.apSum[sl], fitProtons.apSum2[sl], fitProtons.apLo[sl], fitProtons.apHi[sl], fitProtons.apLo2[sl], fitProtons.apHi2[sl]}; prR.n += fitProtons.dN[sl]; prR.sum += fitProtons.dSum[sl]; prR2.n += fitProtons.dN2[sl]; prR2.sum += fitProtons.dSum2[sl];
                    }
                    double p = 0, pErr = 0; (void)p; (void)pErr;
                    if (aR.n > 0 && nPhi >= kMinNFP && std::fabs(aR.Mean()) > 1e-12){
                        // N(phi) = N0 (1 + P <A> cos(phi)) per bin, N0 = entries / bins and <A> fixed, only P is fitted. A_eff = sum AP(theta_i) / N of the protons of the fit (aR), theta_i the own theta of each proton in the frame of the fitted phi
                        const AMean aFix = aR;
                        const int nFitBins = kPhiBinsFP / 2; // 50 bins
                        TH1D* h = new TH1D(Form("hPhiRange%d_%d", range, type), "", nFitBins, 0, 360);
                        h->SetDirectory(nullptr);
                        for (int b = 0; b < kPhiBinsFP; b++) h->AddBinContent(b / 2 + 1, counts[b]);
                        h->SetEntries(nPhi);
                        ovl[range][type] = {h, nPhi, aFix, true, UnbinnedPFromMoments(mom.data(), nPhi, aFix, set.block ? kBlockPhiLo : 0, set.block ? kBlockPhiHi : 360, &prR, &prR2), prR, prR2}; // drawn with the other two types after this loop
                    }
                }
            }
            // the three types of each range on one plot: CE recon protons (red), all recon protons (blue), all recon charged particles (green), each with its own fit of the same color
            for (const int range : {0, 1}){
                if (!ovl[range][0].ok && !ovl[range][1].ok && !ovl[range][2].ok) continue;
                const std::string dir = FomRoot() + "/" + set.dir + "/" + kFitDir;
                const int colors[3] = {kRed, kBlue, kGreen + 2};
                double pFit[3] = {0, 0, 0}, pErr[3] = {0, 0, 0}, hMax = 0;
                c.Clear(); c.SetLogz(0); c.SetLogy(0);
                c.SetRightMargin(0.25); c.SetTopMargin(0.36); c.SetLeftMargin(0.14);
                gStyle->SetOptStat(0); gStyle->SetOptFit(0);
                for (int type = 0; type < 3; type++) if (ovl[range][type].ok){
                    TH1D* h = ovl[range][type].h;
                    AddUnbinnedFit(h, Form("fPhiRange%d_%d", range, type), ovl[range][type].uf, ovl[range][type].aFix.Mean(), ovl[range][type].aFix.Err(), false); // the unbinned fit of the range (sums of cos^k phi), the histogram only shows it
                    pFit[type] = ovl[range][type].uf.p; pErr[type] = ovl[range][type].uf.err;
                    h->SetLineColor(colors[type]); h->SetLineWidth(2); h->SetStats(0);
                    hMax = std::max(hMax, h->GetMaximum());
                }
                bool first = true;
                for (int type = 0; type < 3; type++) if (ovl[range][type].ok){
                    TH1D* h = ovl[range][type].h;
                    h->SetMinimum(0); h->SetMaximum(1.1 * hMax);
                    h->SetTitle((std::string("CE recon protons (only the ones with a recon track) and all recon charged particles, ") + set.name + " angles" + (range == 0 ? ", phi fit over all theta, no cuts (no KE, residual or theta cut)" : ", phi fit over theta 0 to 1 deg (no KE or residual cut)") + ";Angle (degrees);Count").c_str());
                    h->GetXaxis()->SetNdivisions(505);
                    h->Draw(first ? "hist" : "hist same"); first = false;
                }
                for (int type = 0; type < 3; type++) if (ovl[range][type].ok){
                    TF1* fd = new TF1(Form("fPhiRangeDraw%d_%d", range, type), "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                    fd->SetParameters(ovl[range][type].nPhi / ovl[range][type].h->GetNbinsX(), pFit[type], ovl[range][type].aFix.Mean());
                    fd->SetLineColor(colors[type]); fd->SetLineStyle(2); fd->SetLineWidth(3); fd->Draw("same");
                }
                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.025); // small enough for the whole P = ... #pm ... of the longest line to stay on the page
                tx.DrawLatex(0.15, 0.90, "Unbinned fit  N(#phi) #propto 1 + P #upoint A_{eff} cos#phi (every proton), dashed = fit");

tx.DrawLatex(0.15, 0.865, AEffEq());
                for (int type = 0; type < 3; type++) if (ovl[range][type].ok){
                    tx.SetTextColor(colors[type]);
                    tx.DrawLatex(0.15, 0.828 - 0.045 * (type == 2 ? 1 : 0), Form("%s: %.0f entries, A_{eff} = %.4f #pm %.4f,  P = %.2f #pm %.2f", typeNames[type].c_str(), ovl[range][type].nPhi, ovl[range][type].aFix.Mean(), ovl[range][type].aFix.Err(), pFit[type], pErr[type]));
                }
                tx.SetTextColor(kBlack);
                if (set.block) tx.DrawLatex(0.15, 0.72, BlockPhiNote());
                { // the errors of A_eff and P by source of each fit, in the right margin (CE recon protons and all recon charged particles: two small boxes)
                    int nBox = 0;
                    for (int type = 0; type < 3; type++) if (ovl[range][type].ok){
                        DrawErrorBySource(ovl[range][type].aFix, ovl[range][type].uf, true, 0.755, 0.06 + 0.22 * nBox, 0.995, 0.27 + 0.22 * nBox, typeNames[type].c_str(), 0.0155);
                        nBox++;
                    }
                }
                c.Modified(); c.Update();
                c.SaveAs((dir + std::to_string(2 * ((range == 0 ? (set.block ? 7 : 9) : (set.block ? 8 : 10)) - 1) + 1) + "_" + (range == 0 ? "phi_fit_all_theta" : "phi_fit_theta_0_to_1") + ".pdf").c_str());
                if (ovl[range][0].ok){ // the CE recon protons alone, their few entries vanish next to the recon protons and particles above: own pdf, the next number after the overlay
                    TH1D* h = ovl[range][0].h;
                    c.Clear(); c.SetTopMargin(0.27);
                    h->SetLineColor(kRed); h->SetLineWidth(2); h->SetStats(0); h->SetMinimum(0); h->SetMaximum(-1111); h->SetMaximum(1.1 * h->GetMaximum()); // -1111: drop the overlay maximum, the axis follows the CE bins
                    h->SetTitle((std::string("CE recon protons only (not all CE protons), ") + set.name + " angles" + (range == 0 ? ", phi fit over all theta, no cuts (no KE, residual or theta cut)" : ", phi fit over theta 0 to 1 deg (no KE or residual cut)") + ";Angle (degrees);Count").c_str());
                    h->Draw("hist");
                    TF1* fd = new TF1(Form("fPhiRangeCE%d", range), "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
                    fd->SetParameters(ovl[range][0].nPhi / h->GetNbinsX(), pFit[0], ovl[range][0].aFix.Mean()); fd->SetLineColor(kBlack); fd->SetLineWidth(3); fd->Draw("same");
                    DrawErrorBySource(ovl[range][0].aFix, ovl[range][0].uf, true); // the errors of A_eff and P by source
                    TLatex txc; txc.SetNDC(); txc.SetTextSize(0.03);
                    txc.DrawLatex(0.15, 0.90, "Unbinned fit  N(#phi) #propto 1 + P #upoint A_{eff} cos#phi (every proton)");
                    txc.DrawLatex(0.15, 0.855, Form("%.0f entries, %s = %.4f #pm %.4f (per-proton #theta),  P =%.2f #pm %.2f", ovl[range][0].nPhi, AEffEq(), ovl[range][0].aFix.Mean(), ovl[range][0].aFix.Err(), pFit[0], pErr[0]));
                    if (set.block) txc.DrawLatex(0.15, 0.81, BlockPhiNote());
                    c.Modified(); c.Update();
                    c.SaveAs((dir + std::to_string(2 * ((range == 0 ? (set.block ? 7 : 9) : (set.block ? 8 : 10)) - 1) + 2) + "_" + (range == 0 ? "phi_fit_all_theta" : "phi_fit_theta_0_to_1") + "_CE.pdf").c_str());
                    delete fd;
                }
                for (int type = 0; type < 3; type++) delete ovl[range][type].h;
                c.SetRightMargin(0.30); c.SetTopMargin(0.1); c.SetLeftMargin(0.12);
            }
        }
    }
}
