#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "FitCore.h"
#include "CutFiles.h"

// P error budget of one run: every statistical, systematic, bias and dilution item of /Users/trevorg04/Desktop/P_fit_error_sources.md that one sim run can give, one txt per value in
// plots/11_errorBudget/, and P_summary.txt with P = P_fit +- sigma_stat +- delta_syst, the biases and the dilution.
// The sample is the total recon fit of F9_PlotReconNew.cpp's 6_megaCut (every recon track after the KE 80%, even planes, theta 0 to 1 degree and 1 hit per plane cuts, the block phi fit in the est neutron frame),
// so P_fit and errFit here equal pAll and errAll of plots/6_ReconNew_1/6_megaCut/scan_stats.csv. The variations (cuts, fit range) are refits of the same rows.
// Items that need another sim or another ReconAngleAll run (efficiency, P = 0, seed road, plane weights, kTubeSigmas, H/C target mix) are written as "not computed".
namespace {
constexpr double kRad = 180.0 / M_PI;
struct Row {double thMin, thMax, phHalf, th2Min, th2Max, ph2Half, thetaRelND, phiRelND, zdcKEAllD; int eventID, evClass, nPlanesUsed, nPlanesHit, nTrackerHits;};
struct Cuts {double ke = 0.80, thLo = 0.0, thHi = 1.0, phiTrim = 0.0; bool even = true, avgOne = true;}; // ke: fraction of the beam KE, phiTrim: degrees cut off both ends of the phi fit range
struct Sel {double ph, thDeg, a, a1lo, a1hi, a2lo, a2hi, h1, h2; bool r1, r2; int cls; long ceIdx;}; // ph: est neutron phi (deg), cls: evClass, ceIdx: row of the CE truth of the event (-1 none)
struct Result {PFit uf; AMean a; PhiRange q1, q2; size_t n = 0;};
struct MeanErr {double mean = 0, sd = 0, sem = 0, n = 0;};
MeanErr Stats(const std::vector<double>& v){
    MeanErr m; m.n = (double)v.size();
    if (v.empty()) return m;
    for (double x : v) m.mean += x;
    m.mean /= m.n;
    if (v.size() > 1){double s = 0; for (double x : v) s += (x - m.mean) * (x - m.mean); m.sd = std::sqrt(s / (m.n - 1)); m.sem = m.sd / std::sqrt(m.n);}
    return m;
}
double Dof(bool has, double h){return !has || h <= 1e-9 ? 1.0 : (h < M_PI ? std::sin(h) / h : 0.0);} // sin(h)/h of PhiRange::Add
}

void PErrorBudget(){
    TFile* f = TFile::Open("cuts/cut_angles_all.root");
    TTree* t = f ? (TTree*)f->Get("AnglesRaw") : nullptr;
    if (!t){Error("PErrorBudget", "no cuts/cut_angles_all.root (run ReconAngleAll.cpp)"); return;}
    ROOT::EnableImplicitMT(NThreads());
    const auto cols = ColumnsUnordered(t, {"thetaRelNDMin", "thetaRelNDMax", "phiRelNDHalf", "thetaRelNDBMin", "thetaRelNDBMax", "phiRelNDBHalf", "thetaRelND", "phiRelND", "zdcKEAllD", "fEvent", "evClass", "nPlanesUsed", "nPlanesHitEv", "nTrackerHitsEv"});
    std::vector<Row> rows(cols[0].size());
    for (size_t i = 0; i < rows.size(); i++){
        size_t k = 0;
        Row& x = rows[i];
        x.thMin = cols[k++][i] * kRad; x.thMax = cols[k++][i] * kRad; x.phHalf = cols[k++][i]; x.th2Min = cols[k++][i] * kRad; x.th2Max = cols[k++][i] * kRad; x.ph2Half = cols[k++][i];
        x.thetaRelND = cols[k++][i]; x.phiRelND = cols[k++][i]; x.zdcKEAllD = cols[k++][i];
        x.eventID = (int)cols[k++][i]; x.evClass = (int)cols[k++][i]; x.nPlanesUsed = (int)cols[k++][i]; x.nPlanesHit = (int)cols[k++][i]; x.nTrackerHits = (int)cols[k++][i];
    }
    const auto settings = ReadSettings();
    const double apScale = settings.at("APScale"), beamKE = settings.at("BeamKE"), pTrue = settings.at("PPValue");
    const int maxPlanes = MaxTrackerHits();

    // the true CE of each event (ChargeExchange_values, the carbon CE rows): true phi and theta in the true neutron frame (deg), LastAP, LastT (model t), LastTheta (lab theta)
    TTree* ceT = OpenCut("cuts/cut_CE.root", "ChargeExchange_values");
    gROOT->cd();
    if (!ceT){Error("PErrorBudget", "no cuts/cut_CE.root"); return;}
    const auto ce = Columns(ceT, {"fEvent", "CEPhiNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "CEThetaNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "LastAP", "LastT", "LastTheta"});
    std::map<long long, size_t> byEvent;
    for (size_t i = 0; i < ce[0].size(); i++) byEvent.emplace((long long)ce[0][i], i);

    // the rows of a cut set (F9_PlotReconNew.cpp passes + inFit): every recon track, evClass 0 to 2
    const auto select = [&](const Cuts& c){
        std::vector<Sel> s;
        const double keCut = c.ke * beamKE;
        for (const Row& x : rows){
            if (x.evClass < 0 || x.evClass > 2) continue;
            if (!(std::isfinite(x.zdcKEAllD) && x.zdcKEAllD > keCut)) continue;
            if (c.even && !EvenHits(x.nPlanesUsed, maxPlanes)) continue;
            const double th = x.thetaRelND * kRad, ph = x.phiRelND * kRad;
            if (!(th >= c.thLo && th <= c.thHi)) continue; // NAN fails too
            if (c.avgOne && !(x.nPlanesHit > 0 && x.nTrackerHits == x.nPlanesHit)) continue;
            if (!(ph > kBlockPhiLo + c.phiTrim && ph < kBlockPhiHi - c.phiTrim)) continue;
            Sel q;
            q.ph = ph; q.thDeg = th; q.a = AEffOfTheta(th, apScale);
            q.r1 = std::isfinite(x.thMin) && std::isfinite(x.thMax); q.r2 = std::isfinite(x.th2Min) && std::isfinite(x.th2Max);
            q.a1lo = q.r1 ? AEffOfTheta(x.thMin, apScale) : q.a; q.a1hi = q.r1 ? AEffOfTheta(x.thMax, apScale) : q.a;
            q.a2lo = q.r2 ? AEffOfTheta(x.th2Min, apScale) : q.a; q.a2hi = q.r2 ? AEffOfTheta(x.th2Max, apScale) : q.a;
            q.h1 = x.phHalf; q.h2 = x.ph2Half; q.cls = x.evClass;
            const auto it = byEvent.find((long long)x.eventID);
            q.ceIdx = it == byEvent.end() ? -1 : (long)it->second;
            s.push_back(q);
        }
        return s;
    };
    // the signed A_eff bias of the recon theta: CE recon rows, A at the recon theta minus A at the true CE theta (F9 smearBias), with the error of that mean
    const auto smearOf = [&](const std::vector<Sel>& s){
        std::vector<double> d;
        for (const Sel& q : s) if (q.cls == 0 && q.ceIdx >= 0) d.push_back(q.a - AEffOfTheta(ce[2][q.ceIdx], apScale));
        return Stats(d);
    };
    const auto fit = [&](const std::vector<Sel>& s, int maxCls, double smear, double trim){
        Result r;
        std::vector<double> ph;
        for (const Sel& q : s){
            if (q.cls > maxCls) continue;
            ph.push_back(q.ph); r.a.Add(q.a, q.a1lo, q.a1hi, q.a2lo, q.a2hi);
            if (q.r1) r.q1.Add(q.h1);
            if (q.r2) r.q2.Add(q.h2);
        }
        r.a.smearBias = smear; r.n = ph.size();
        if (ph.size() >= 10) r.uf = UnbinnedP(ph, r.a, kBlockPhiLo + trim, kBlockPhiHi - trim, nullptr, true, &r.q1, &r.q2);
        return r;
    };

    // ---------- baseline: the megaCut total recon fit ----------
    const Cuts base{};
    const std::vector<Sel> s0 = select(base);
    const MeanErr sm0 = smearOf(s0);
    const Result B = fit(s0, 2, sm0.mean, 0.0);
    const PFit& u = B.uf;
    if (!u.ok){Error("PErrorBudget", "baseline fit failed, %zu rows", B.n); return;}
    const double P = u.p, A = B.a.Mean(), N = (double)B.n;

    const std::string dir = "plots/11_errorBudget";
    gSystem->mkdir(dir.c_str(), true);
    const auto put = [&](const char* name, const std::string& body){std::ofstream o(dir + "/" + name + ".txt"); o << body;};
    const auto fmt = [](const char* f, auto... a){return std::string(Form(f, a...));}; // every argument is a double

    // ---------- statistical ----------
    const double sigFit = u.errFit, sigAStat = u.errAStat;
    put("sigma_P_fit", fmt("sigma_P_fit (likelihood error errFit, statistics of the phi distribution)\nvalue = %.5f\nsample = total recon, megaCut, N = %.0f fit tracks\ncheck sqrt(2/N)/|A_eff| = %.5f (A_eff = %.5f)\n", sigFit, N, std::sqrt(2.0 / N) / std::fabs(A), A));
    put("sigma_AP_stat", fmt("sigma_AP_stat (spread of AP(theta_i) over the fit tracks, AMean::StatErr = std(A)/sqrt(N))\nvalue on A_eff = %.6f (A_eff = %.5f, relative %.5f)\nvalue on P = %.5f\nnote it may double count with sigma_P_fit, since A_eff comes from the same tracks that are fit\n", B.a.StatErr(), A, B.a.StatErr() / std::fabs(A), sigAStat));
    // bootstrap of the fit rows
    const int nBoot = 200;
    std::vector<double> bp, bpr;
    {
        TRandom3 rng(20260101);
        std::vector<const Sel*> pick;
        for (const Sel& q : s0) pick.push_back(&q);
        for (int b = 0; b < nBoot; b++){
            std::vector<double> ph(pick.size());
            AMean a; PhiRange q1, q2;
            for (size_t i = 0; i < pick.size(); i++){
                const Sel& q = *pick[(size_t)(rng.Rndm() * pick.size())];
                ph[i] = q.ph; a.Add(q.a, q.a1lo, q.a1hi, q.a2lo, q.a2hi);
                if (q.r1) q1.Add(q.h1);
                if (q.r2) q2.Add(q.h2);
            }
            a.smearBias = sm0.mean;
            const PFit v = UnbinnedP(ph, a, kBlockPhiLo, kBlockPhiHi, nullptr, false, &q1, &q2);
            if (v.ok){bp.push_back(v.p); bpr.push_back(v.pRaw);}
        }
    }
    const MeanErr bs = Stats(bp), bsr = Stats(bpr);
    put("sigma_P_bootstrap", fmt("sigma_P_bootstrap (%.0f resamples with replacement of the %.0f fit rows, refit each, corrected P)\nspread of P = %.5f (mean %.5f)\nspread of P before correction = %.5f\nratio to errFit = %.3f\n", bs.n, N, bs.sd, bs.mean, bsr.sd, sigFit != 0 ? bs.sd / sigFit : 0.0));
    // MC statistics of the dilution and smear corrections
    std::vector<double> d1v, d2v, ad;
    for (const Sel& q : s0){d1v.push_back(Dof(q.r1, q.h1)); d2v.push_back(Dof(q.r2, q.h2)); ad.push_back(q.a * Dof(q.r1, q.h1));}
    const MeanErr d1 = Stats(d1v), d2 = Stats(d2v);
    put("sigma_D_phi_MC", fmt("sigma_D_phi_MC (statistics of the mean D = sin(h)/h over the fit tracks)\nD_phi = %.5f +- %.5f, relative %.5f -> on P %.5f\nD_block = %.5f +- %.5f, relative %.5f -> on P %.5f\n", d1.mean, d1.sem, d1.sem / d1.mean, std::fabs(P) * d1.sem / d1.mean, d2.mean, d2.sem, d2.sem / d2.mean, std::fabs(P) * d2.sem / d2.mean));
    put("sigma_smearBias_MC", fmt("sigma_smearBias_MC (statistics of the smearBias mean over the CE recon rows, N = %.0f)\nsmearBias = %+.6f +- %.6f (A_eff units)\non P = %.5f\n", sm0.n, sm0.mean, sm0.sem, std::fabs(P) * sm0.sem / std::fabs(A)));
    // purity
    double nCls[3] = {0, 0, 0};
    for (const Sel& q : s0) for (int k = q.cls; k < 3; k++) nCls[k]++;
    const double f0 = nCls[0] / nCls[2], f1 = nCls[1] / nCls[2];
    put("sigma_purity", fmt("sigma_purity (binomial error of the CE fraction of the fit sample, N = %.0f)\nCE only fraction f0 = %.5f +- %.5f (%.0f rows)\nCE + CE mix fraction f1 = %.5f +- %.5f (%.0f rows)\n", nCls[2], f0, std::sqrt(f0 * (1 - f0) / nCls[2]), nCls[0], f1, std::sqrt(f1 * (1 - f1) / nCls[2]), nCls[1]));

    // ---------- systematic ----------
    const double rel2 = B.n > 0 && A != 0 ? std::fabs(B.a.sHi2 - B.a.sLo2) / (2 * N) / std::fabs(A) : 0.0;
    const double dModel = u.errAModel, dBeam = u.errABeam, dRange = u.errARange, dRange2 = std::fabs(P) * rel2;
    put("delta_AP_model", fmt("delta_AP_model (errors of the fitted AP model coefficients A = 0.7079 +- 0.0303, k = 0.3860 +- 0.0523, no covariance)\nrelative on A_eff = %.5f\nvalue on P = %.5f\n", kModelAFrac, dModel));
    put("delta_AP_beam", fmt("delta_AP_beam (beam momentum spread %.1f MeV/c of %.1f MeV/c)\nrelative on A_eff = %.5f\nvalue on P = %.5f\n", kBeamMomentumSD, kBeamMomentum, kBeamAFrac, dBeam));
    put("delta_AP_range", fmt("delta_AP_range (est neutron theta range of each track: half the A_eff difference between the two ends, AMean::RangeErr)\nvalue on A_eff = %.6f\nvalue on P = %.5f\n", B.a.RangeErr(), dRange));
    put("delta_AP_range_block", fmt("delta_AP_range_block (block grid theta range, the second range sLo2 and sHi2; computed here, AddAErrToP does not carry it)\nrelative on A_eff = %.5f\nvalue on P = %.5f\n", rel2, dRange2));
    // t mismatch: CE rows of the fit sample, A_eff of the analysis (transverse t at the true CE theta) against the model AP of the CE (LastAP)
    std::vector<double> aCode, aModel, tModel, tCode;
    for (const Sel& q : s0) if (q.cls == 0 && q.ceIdx >= 0){
        const double thN = ce[2][q.ceIdx];
        aCode.push_back(AEffOfTheta(thN, apScale)); aModel.push_back(ce[3][q.ceIdx]);
        tModel.push_back(ce[4][q.ceIdx]); tCode.push_back(TFromAngle(kBeamMomentum, thN / kRad));
    }
    const MeanErr mAC = Stats(aCode), mAM = Stats(aModel), mTM = Stats(tModel), mTC = Stats(tCode);
    const double rT = mAC.mean != 0 ? mAM.mean / mAC.mean : 1.0; // P_code = P_true * A_model / A_code
    const double dTmis = std::fabs(P) * std::fabs(rT - 1.0);
    // the error on that ratio from the per-track difference A_model,i - A_analysis,i (the two means use the same tracks, so their separate errors overstate it): ratio - 1 = mean(d) / mean(A_analysis)
    std::vector<double> aDiff;
    for (size_t i = 0; i < aCode.size(); i++) aDiff.push_back(aModel[i] - aCode[i]);
    const MeanErr mD = Stats(aDiff);
    const double sigRatio = mAC.mean != 0 ? mD.sem / std::fabs(mAC.mean) : 0.0;
    // the same with t from the block theta in the est neutron frame (what the fit really uses) against the model AP of the CE: this holds the t definition and the theta smear together
    std::vector<double> aBlk;
    for (const Sel& q : s0) if (q.cls == 0 && q.ceIdx >= 0) aBlk.push_back(q.a);
    const MeanErr mAB = Stats(aBlk);
    const double rTB = mAB.mean != 0 ? mAM.mean / mAB.mean : 1.0;
    put("delta_AP_tmismatch_block", fmt("delta_AP_tmismatch_block (analysis t from the est n frame block theta against the model AP of the CE, %.0f CE recon rows; t definition and theta smear together)\nA_eff block theta = %.5f +- %.5f, A_eff model (LastAP) = %.5f +- %.5f\nratio A_model / A_block = %.4f (the t-only ratio is %.4f)\nvalue on P = %.5f\n", mAB.n, mAB.mean, mAB.sem, mAM.mean, mAM.sem, rTB, rT, std::fabs(P) * std::fabs(rTB - 1.0)));
    put("delta_AP_tmismatch", fmt("delta_AP_tmismatch (analysis t = (p sin theta)^2 at the true CE theta against the model t and AP of the CE, %.0f CE recon rows of the sample)\nA_eff analysis = %.5f +- %.5f, A_eff model (LastAP) = %.5f +- %.5f\nmean t model = %.6g, mean t analysis = %.6g (GeV^2 if LastT is in GeV^2)\nratio A_model / A_analysis = %.4f -> P_code = P_true x that\nvalue on P = %.5f\nerror on the ratio from the per-track A_model - A_analysis = %.5f -> on P %.5f\n", mAC.n, mAC.mean, mAC.sem, mAM.mean, mAM.sem, mTM.mean, mTC.mean, rT, dTmis, sigRatio, std::fabs(P) * sigRatio));
    // cut and fit range variations: refit the same rows with one setting changed. sd = sqrt(|err_var^2 - err_base^2|) is the statistical size of the shift (nested samples); delta = the largest sqrt(max(0, shift^2 - sd^2)), the shift beyond statistics
    struct Var {const char* label; Cuts c;};
    const auto variation = [&](const char* name, const char* what, std::vector<Var> vars){
        std::string txt = std::string(name) + " (" + what + ")\nbase P = " + std::string(Form("%.5f", P)) + " +- " + std::string(Form("%.5f", u.err)) + ", N = " + std::string(Form("%.0f", N)) + "\n";
        double dMax = 0;
        for (const Var& v : vars){
            const std::vector<Sel> sv = select(v.c);
            const Result r = fit(sv, 2, smearOf(sv).mean, v.c.phiTrim);
            const double dP = r.uf.p - P, sd = std::sqrt(std::fabs(r.uf.err * r.uf.err - u.err * u.err));
            const double sig = std::sqrt(std::max(0.0, dP * dP - sd * sd)); // the part of the shift beyond its statistical size
            if (r.uf.ok) dMax = std::max(dMax, sig);
            txt += std::string(Form("%s: N = %zu, P = %.5f +- %.5f, P - P_base = %+.5f, statistical size of the shift %.5f, shift beyond statistics %.5f\n", v.label, r.n, r.uf.p, r.uf.err, dP, sd, sig));
        }
        txt += std::string(Form("value on P (largest shift beyond statistics) = %.5f\n", dMax));
        put(name, txt);
        return dMax;
    };
    Cuts c;
    c = base; c.ke = 0.70; Cuts ke70 = c; c = base; c.ke = 0.90; Cuts ke90 = c;
    const double dKE = variation("delta_KEcut", "KE cut 70% and 90% of the beam KE instead of 80%", {{"KE 70%", ke70}, {"KE 90%", ke90}});
    c = base; c.thHi = 0.8; Cuts th08 = c; c = base; c.thHi = 1.2; Cuts th12 = c;
    const double dTh = variation("delta_theta_cut", "theta cut 0 to 0.8 and 0 to 1.2 degree instead of 0 to 1", {{"theta 0 to 0.8", th08}, {"theta 0 to 1.2", th12}});
    c = base; c.even = false; Cuts noEven = c;
    const double dPl = variation("delta_planes_cut", "even planes cut dropped", {{"no even planes cut", noEven}});
    c = base; c.avgOne = false; Cuts noAvg = c;
    const double dHit = variation("delta_hits_cut", "one hit per plane cut dropped", {{"no 1 hit per plane cut", noAvg}});
    c = base; c.phiTrim = 1.0; Cuts trim1 = c;
    const double dFr = variation("delta_fit_range", "phi fit range 1 to 359 degree instead of 0 to 360", {{"phi 1 to 359", trim1}});
    const char* notComputed[][2] = {
        {"delta_efficiency", "tracker efficiency 0.95 needs another sim (/detector/trackerEfficiency 0.95) and the same analysis"},
        {"delta_false_asym", "the false asymmetry needs a P = 0 sim and +0.7 and -0.7 sims (half sum = false asymmetry, half difference = P)"},
        {"delta_seed_road", "the 5 mrad scatter road (kSeedScatterRad, SeedExtend.h) needs ReconAngleAll.cpp rerun with the road changed"},
        {"delta_plane_weights", "the plane weights 1, 1, 0.5, 0.25 need ReconAngleAll.cpp rerun with other weights"},
        {"delta_kTubeSigmas", "kTubeSigmas = 2 (tube holds 92% of the CE vertices) needs ReconAngleAll.cpp rerun with 1 and 3"},
        {"delta_target_HC", "the CH2 target (CE on H with t slope 14.5 and on C with 76) needs a sim with the target material split or a CE target column"}};
    for (const auto& nc : notComputed) put(nc[0], std::string(nc[0]) + "\nnot computed: " + nc[1] + "\n");

    // ---------- bias and dilution ----------
    const double Dphi = B.q1.D(), Dblk = B.q2.D();
    put("D_phi", fmt("D_phi (est neutron phi range dilution, mean sin(h)/h over the fit tracks, corrected: P = P_raw / D)\nD = %.5f +- %.5f\nshift of P = %+.5f\n", Dphi, d1.sem, u.biasPhi));
    put("delta_phi_range", fmt("delta_phi_range (est neutron phi range: the uncertainty of the phi dilution correction, taken as its full size |P| x |1/D - 1|, D = mean sin(h)/h over the fit tracks)\nD = %.5f\nvalue on P = %.5f\n", Dphi, std::fabs(u.biasPhi)));
    put("D_block", fmt("D_block (block grid phi range dilution, mean sin(h)/h of the block range, corrected: P = P / D)\nD = %.5f +- %.5f\nshift of P = %+.5f\n", Dblk, d2.sem, u.biasPhi2));
    // actual dilution: CE recon rows, est neutron phi against the true neutron frame phi of the CE
    std::vector<double> cosD, dModelRows;
    for (const Sel& q : s0) if (q.cls == 0 && q.ceIdx >= 0 && q.r1){
        cosD.push_back(std::cos((q.ph - ce[1][q.ceIdx]) / kRad)); dModelRows.push_back(Dof(q.r1, q.h1));
    }
    const MeanErr mCos = Stats(cosD), mDm = Stats(dModelRows);
    put("D_actual", fmt("D_actual (CE recon rows: mean cos(phi_est - phi_true in the true neutron frame), the dilution really present; N = %.0f)\nD_actual = %.5f +- %.5f\nD of the code (mean sin(h)/h of the same rows) = %.5f +- %.5f\nD_actual / D_code = %.4f (1 = the code correction is right)\n", mCos.n, mCos.mean, mCos.sem, mDm.mean, mDm.sem, mDm.mean != 0 ? mCos.mean / mDm.mean : 0.0));
    // fits of the four samples for the purity table
    const Result rAll = B, rMix = fit(s0, 1, sm0.mean, 0.0), rCE = fit(s0, 0, sm0.mean, 0.0);
    PFit ufTrue; double nTrue = 0;
    {
        std::vector<double> ph; AMean a;
        for (const Sel& q : s0) if (q.cls == 0 && q.ceIdx >= 0){ph.push_back(ce[1][q.ceIdx]); const double a0 = ce[3][q.ceIdx]; a.Add(a0, a0, a0, a0, a0);}
        nTrue = (double)ph.size();
        if (ph.size() >= 10) ufTrue = UnbinnedP(ph, a, kBlockPhiLo, kBlockPhiHi, nullptr, true, nullptr, nullptr);
    }
    put("D_purity", fmt("D_purity (dilution by the non-CE tracks in the sample: P of each sample against the true P = %.2f)\n", pTrue) +
        fmt("CE only, true CE angle and LastAP: N = %.0f, P = %.5f +- %.5f\n", nTrue, ufTrue.p, ufTrue.err) +
        fmt("CE only, recon est angle: N = %.0f, P = %.5f +- %.5f\n", (double)rCE.n, rCE.uf.p, rCE.uf.err) +
        fmt("CE + CE mix, recon: N = %.0f, P = %.5f +- %.5f\n", (double)rMix.n, rMix.uf.p, rMix.uf.err) +
        fmt("total recon: N = %.0f, P = %.5f +- %.5f\n", (double)rAll.n, rAll.uf.p, rAll.uf.err) +
        fmt("CE only fraction of the total sample f0 = %.4f, CE + CE mix f1 = %.4f; if the non-CE tracks are flat in phi the total P is the CE P times f\n", f0, f1));
    // R: P_code = P_true R, R = <A d> / (<A> <d>) over the fit tracks (A and d of one track are correlated through theta)
    const auto Rof = [&](int maxCls){
        double n = 0, sa = 0, sd = 0, sad = 0;
        for (const Sel& q : s0){if (q.cls > maxCls) continue; const double d = Dof(q.r1, q.h1); n++; sa += q.a; sd += d; sad += q.a * d;}
        return n > 0 && sa != 0 && sd != 0 ? (sad / n) / ((sa / n) * (sd / n)) : 1.0;
    };
    const double R2 = Rof(2), R0 = Rof(0);
    put("R_corr", fmt("R_corr (product of means: R = <A_i d_i> / (<A><d>), P_code = P_true x R to first order)\nR total recon = %.5f -> P / R = %.5f (shift %+.5f)\nR CE only = %.5f\n", R2, P / R2, P / R2 - P, R0));
    put("bias_smear", fmt("bias_smear (recon theta smear on A_eff, smearBias = mean AP(theta_rec) - mean AP(theta_true) over the CE recon rows, N = %.0f; corrected in the fit)\nsmearBias = %+.6f +- %.6f (A_eff units), %+.3f%% of A_eff\nshift of P = %+.5f\n", sm0.n, sm0.mean, sm0.sem, 100.0 * sm0.mean / A, u.biasSmear));
    put("bias_P_true_angles", fmt("bias_P_true_angles (the fit model itself: CE rows at the true CE angle and LastAP, no recon, no dilution)\nP = %.5f +- %.5f, true P = %.5f, P - P_true = %+.5f (%.2f sigma), N = %.0f\n", ufTrue.p, ufTrue.err, pTrue, ufTrue.p - pTrue, ufTrue.err > 0 ? (ufTrue.p - pTrue) / ufTrue.err : 0.0, nTrue));

    // ---------- totals ----------
    const double sigTheta = std::fabs(P) * sm0.sem / std::fabs(A), sigPhi = std::fabs(P) * d1.sem / d1.mean; // statistics of the smearBias and D means
    const double sigStat = std::sqrt(sigFit * sigFit + sigAStat * sigAStat + sigTheta * sigTheta + sigPhi * sigPhi);
    const std::vector<std::pair<const char*, double>> deltas = {{"AP model coefficients", dModel}, {"beam momentum spread", dBeam}, {"est neutron theta range", dRange}, {"block grid theta range", dRange2}, {"est neutron phi range", std::fabs(u.biasPhi)}, {"t definition mismatch", dTmis}, {"KE cut", dKE}, {"theta cut", dTh}, {"even planes cut", dPl}, {"1 hit per plane cut", dHit}, {"phi fit range", dFr}};
    double dSq = 0;
    for (const auto& d : deltas) dSq += d.second * d.second;
    const double dSyst = std::sqrt(dSq), sigTot = std::sqrt(sigStat * sigStat + dSq);
    const double pull = sigTot > 0 ? (P - pTrue) / sigTot : 0.0;
    put("bias_closure", fmt("bias_closure (corrected P against the true P of the run)\nP = %.5f, P_true = %.5f, P - P_true = %+.5f\nsigma_stat = %.5f, delta_syst = %.5f, total = %.5f\npull = (P - P_true) / total = %+.3f (one run, a pull needs many runs: error on the mean 1/sqrt(N runs), on the RMS 1/sqrt(2 N runs))\nP - P_true in sigma_stat only = %+.3f\n", P, pTrue, P - pTrue, sigStat, dSyst, sigTot, pull, sigStat > 0 ? (P - pTrue) / sigStat : 0.0));
    std::string sum = fmt("P_summary, run %.0f blocks, %.0f neutrons, true P = %.3f, beam KE %.0f MeV\nsample: total recon, 6_megaCut (KE %.0f%%, even planes, theta 0 to 1 degree, 1 hit per plane), N = %.0f fit tracks\n\n", settings.at("nBlocks"), settings.at("TotalNeutrons"), pTrue, beamKE, 100 * base.ke, N);
    sum += fmt("P = %.4f  +-  %.4f (stat)  +-  %.4f (syst)\n", P, sigStat, dSyst);
    if (std::fabs(u.pRaw) > 0.9999) sum += fmt("WARNING the fit is at the +-1 bound (P before correction = %.4f): the likelihood maximum is outside [-1, 1], so errFit and the bootstrap (%.4f) are not Gaussian errors, and the corrected P can exceed 1\n", u.pRaw, bs.sd);
    sum += "\n";
    sum += fmt("STATISTICAL sigma_stat = %.4f = errFit %.4f (+) AP spread %.4f (+) theta %.4f (+) phi %.4f; bootstrap %.4f\n", sigStat, sigFit, sigAStat, sigTheta, sigPhi, bs.sd);
    sum += "SYSTEMATIC delta_syst, each on P (quadrature):\n";
    for (const auto& d : deltas) sum += fmt("  %.5f  ", d.second) + d.first + "\n";
    sum += "  not computed (need other runs) delta_efficiency, delta_false_asym, delta_seed_road, delta_plane_weights, delta_kTubeSigmas, delta_target_HC\n\n";
    sum += "BIASES AND DILUTION (applied in P, not in the errors):\n";
    sum += fmt("  P before any correction = %.4f\n", u.pRaw);
    sum += fmt("  smear of recon theta on A_eff: smearBias %+.5f, shift of P %+.4f\n", sm0.mean, u.biasSmear);
    sum += fmt("  est neutron phi dilution D_phi = %.4f, shift of P %+.4f (actual D = %.4f +- %.4f)\n", Dphi, u.biasPhi, mCos.mean, mCos.sem);
    sum += fmt("  block grid phi dilution D_block = %.4f, shift of P %+.4f\n", Dblk, u.biasPhi2);
    sum += fmt("  total dilution D = %.4f\n", u.dil);
    sum += fmt("  product of means R = %.4f (P / R = %.4f, not applied)\n", R2, P / R2);
    sum += fmt("  non-CE purity: CE only fraction %.3f, CE + mix %.3f (not corrected)\n", f0, f1);
    sum += fmt("  t mismatch A_model / A_analysis = %.4f (not applied)\n", rT);
    sum += fmt("  fit model bias with true angles: P - P_true = %+.4f +- %.4f\n\n", ufTrue.p - pTrue, ufTrue.err);
    sum += fmt("CLOSURE: P - P_true = %+.4f, pull = %+.3f (total error %.4f)\n", P - pTrue, pull, sigTot);
    put("P_summary", sum);
    printf("%s", sum.c_str());
    printf(">>> wrote %s/*.txt\n", dir.c_str());
}
