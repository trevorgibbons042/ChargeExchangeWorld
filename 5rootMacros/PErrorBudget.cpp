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

// P error budget of one run: the statistical, systematic, bias and dilution items of /Users/trevorg04/Desktop/P_fit_error_sources.md that one sim run can give, a txt per item that feeds the tree in
// plots/11_errorBudget/ (sigma_P_fit, sigma_theta, delta_*, D_block, D_track, bias_closure), P_summary.txt with P = P_fit +- sigma_stat +- delta_syst, and P_table.md, the same as a tree.
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
// np CE cross section parameters A, alpha, g0, g1, c0, c1: the values of G4ChargeExchangeNP.cc (Lyubovitsky, Barton et al. 1976) and the default errors of Ddelta_CEXS
constexpr double kCEXSNom[6] = {1828.5, 0.20096, 0.03129, -0.0084418, 7.6817, -0.5507};
constexpr double kCEXSErr[6] = {27.73083269, 0.006177454, 0.354852944, 0.129440399, 2.640354537, 0.804150501};
constexpr bool kCEXSVary[6] = {true, true, false, false, false, true}; // g0, g1, c0 errors blow the cross section up to x0.46 to x2.68, so they stay at their nominal value
// the cross section of GetCrossSection without the Z factor (it cancels in a ratio), s in GeV^2, s0 = 10 GeV^2
double CEXSShape(const double* p, double s){
    const double sS0 = s / 10.0, lg = std::log(sS0);
    const double G = std::max(0.0, 1 + p[2] + p[3] * lg), C = std::max(1e-16, 1 + p[4] + p[5] * lg);
    return p[0] * std::pow(sS0, 2 * p[1] - 2) * (1 + G) / C;
}
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
    // P_result.csv: item,quantity,value,error,note, one row per number of the txt files (commas in a note become semicolons, error empty when none)
    std::string csvRows;
    const auto csvRow = [&](const std::string& item, const std::string& q, double v, double e, std::string note){
        std::replace(note.begin(), note.end(), ',', ';');
        csvRows += item + "," + q + "," + std::string(Form("%.6g", v)) + "," + (std::isnan(e) ? std::string() : std::string(Form("%.6g", e))) + "," + note + "\n";
    };
    const auto row = [&](const std::string& item, const std::string& q, double v, const std::string& note = ""){csvRow(item, q, v, NAN, note);};
    const auto rowE = [&](const std::string& item, const std::string& q, double v, double e, const std::string& note = ""){csvRow(item, q, v, e, note);};
    const auto rowNote = [&](const std::string& item, std::string note){std::replace(note.begin(), note.end(), ',', ';'); csvRows += item + ",status,,," + note + "\n";};

    // ---------- statistical ----------
    const double sigFit = u.errFit;
    put("sigma_P_fit", fmt("sigma_P_fit (likelihood error errFit, statistics of the phi distribution)\nvalue = %.5f\nsample = total recon, megaCut, N = %.0f fit tracks\ncheck sqrt(2/N)/|A_eff| = %.5f (A_eff = %.5f)\n", sigFit, N, std::sqrt(2.0 / N) / std::fabs(A), A));
    // bootstrap of the fit rows
    const int nBoot = 200;
    std::vector<double> bp;
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
            if (v.ok) bp.push_back(v.p);
        }
    }
    const MeanErr bs = Stats(bp);
    // purity
    double nCls[3] = {0, 0, 0};
    for (const Sel& q : s0) for (int k = q.cls; k < 3; k++) nCls[k]++;
    const double f0 = nCls[0] / nCls[2], f1 = nCls[1] / nCls[2];

    // ---------- systematic ----------
    const double dModel = u.errAModel, dBeam = u.errABeam, dRange = u.errARange;
    put("delta_AP_model", fmt("delta_AP_model (errors of the fitted AP model coefficients A = 0.7079 +- 0.0303, k = 0.3860 +- 0.0523, no covariance)\nrelative on A_eff = %.5f\nvalue on P = %.5f\n", kModelAFrac, dModel));
    put("delta_AP_beam", fmt("delta_AP_beam (beam momentum spread %.1f MeV/c of %.1f MeV/c)\nrelative on A_eff = %.5f\nvalue on P = %.5f\n", kBeamMomentumSD, kBeamMomentum, kBeamAFrac, dBeam));
    put("delta_AP_range", fmt("delta_AP_range (est neutron theta range of each track: half the A_eff difference between the two ends, AMean::RangeErr)\nvalue on A_eff = %.6f\nvalue on P = %.5f\n", B.a.RangeErr(), dRange));
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
        rowE(name, "base P", P, u.err); row(name, "base N", N);
        for (const Var& v : vars){
            const std::vector<Sel> sv = select(v.c);
            const Result r = fit(sv, 2, smearOf(sv).mean, v.c.phiTrim);
            const double dP = r.uf.p - P, sd = std::sqrt(std::fabs(r.uf.err * r.uf.err - u.err * u.err));
            const double sig = std::sqrt(std::max(0.0, dP * dP - sd * sd)); // the part of the shift beyond its statistical size
            if (r.uf.ok) dMax = std::max(dMax, sig);
            row(name, std::string(v.label) + " N", (double)r.n); rowE(name, std::string(v.label) + " P", r.uf.p, r.uf.err); row(name, std::string(v.label) + " P - P_base", dP);
            row(name, std::string(v.label) + " statistical size of the shift", sd); row(name, std::string(v.label) + " shift beyond statistics", sig);
            txt += std::string(Form("%s: N = %zu, P = %.5f +- %.5f, P - P_base = %+.5f, statistical size of the shift %.5f, shift beyond statistics %.5f\n", v.label, r.n, r.uf.p, r.uf.err, dP, sd, sig));
        }
        txt += std::string(Form("value on P (largest shift beyond statistics) = %.5f\n", dMax));
        row(name, "value on P", dMax, "largest shift beyond statistics");
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
    for (const auto& nc : notComputed){put(nc[0], std::string(nc[0]) + "\nnot computed: " + nc[1] + "\n"); rowNote(nc[0], std::string("not computed: ") + nc[1]);}

    // ---------- bias and dilution ----------
    const double Dblk = B.q2.D();
    std::vector<double> d2v;
    for (const Sel& q : s0) d2v.push_back(Dof(q.r2, q.h2));
    const MeanErr d2 = Stats(d2v);
    put("D_block", fmt("D_block (block grid phi range dilution, mean sin(h)/h of the block range, shown, not applied)\nD = %.5f +- %.5f\nshift of P = %+.5f\n", Dblk, d2.sem, u.biasPhi2));
    // actual dilution: CE recon rows, est neutron phi against the true neutron frame phi of the CE
    std::vector<double> cosD;
    for (const Sel& q : s0) if (q.cls == 0 && q.ceIdx >= 0 && q.r1) cosD.push_back(std::cos((q.ph - ce[1][q.ceIdx]) / kRad));
    const MeanErr mCos = Stats(cosD);

    // ---------- totals, in the structure of the error list: P fit uncertainty (sigma_P, delta_P), Dilution (D), Bias ----------
    // Dtrack: the true tracks (2+ tracker hits) that sit in a pure reconstructed track over all true tracks with 2+ tracker hits, no cuts (F8_TrackerFits.cpp writes the two counts from the sim output)
    // Dtrack is the CE child protons only (nPureTrueCE of nTrueTracksCE); the all particle ratio (nPureTrue of nTrueTr) is only in D_track.txt
    double nPureTrue = 0, nTrueTr = 0, nPureTrueCE = 0, nTrueTrCE = 0;
    { std::ifstream in("plots/8_trackerfits/track_purity.csv"); std::string hdr; char c1 = 0, c2 = 0, c3 = 0; if (!(std::getline(in, hdr) && (in >> nPureTrue >> c1 >> nTrueTr >> c2 >> nPureTrueCE >> c3 >> nTrueTrCE))) nTrueTrCE = 0; }
    const bool haveTrk = nTrueTrCE > 0;
    const double Dtrk = haveTrk ? nPureTrueCE / nTrueTrCE : 0.0, sigTrk = haveTrk ? std::sqrt(Dtrk * (1 - Dtrk) / nTrueTrCE) : 0.0; // binomial error
    put("D_track", haveTrk ? fmt("D_track (CE child protons with 2+ tracker hits that sit in a pure reconstructed track, all hits one track ID, no cuts)\nD_track = %.5f +- %.5f (binomial), %.0f of %.0f\nall particles with 2+ tracker hits: %.5f, %.0f of %.0f\n", Dtrk, sigTrk, nPureTrueCE, nTrueTrCE, nTrueTr > 0 ? nPureTrue / nTrueTr : 0.0, nPureTrue, nTrueTr) : std::string("D_track\nnot computed: no plots/8_trackerfits/track_purity.csv with the CE columns (run F8_TrackerFits.cpp)\n"));
    // sigma_phi: statistics of the actual phi dilution. delta_phi: the est n frame phi dilution really present, counted once here (D_actual) and not again as Ddelta_nphi
    const double Dact = mCos.mean;
    const double sigPhi = d2.mean != 0 ? std::fabs(P) * d2.sem / std::fabs(d2.mean) : 0.0; // std(D_block)/sqrt(N) over the fit tracks, like sigma_theta
    // delta_phi: half the P span of the est n frame phi range, every fit track's phi moved by -h and by +h (h = phiRelNDHalf), refit, like delta_theta is half the A_eff span of the theta range
    const auto phiShift = [&](double sgn){
        std::vector<double> ph;
        for (const Sel& q : s0) ph.push_back(q.ph + (q.r1 ? sgn * q.h1 * kRad : 0.0));
        return UnbinnedP(ph, B.a, kBlockPhiLo, kBlockPhiHi, nullptr, false, &B.q1, &B.q2);
    };
    const PFit pPhiLo = phiShift(-1.0), pPhiHi = phiShift(+1.0);
    const double dPhi = pPhiLo.ok && pPhiHi.ok ? std::fabs(pPhiHi.p - pPhiLo.p) / 2 : 0.0;
    const double sigTheta = u.errAStat;
    put("delta_phi", fmt("delta_phi (est n frame phi range, every fit track's phi moved by -h and +h, h = phiRelNDHalf, refit; half the P span)\nP at -h = %.5f, P at +h = %.5f\nvalue on P = %.5f\nsigma_phi (std(D_block)/sqrt(N) over the fit tracks, D_block = %.5f +- %.5f, N = %.0f) on P = %.5f\n", pPhiLo.p, pPhiHi.p, dPhi, d2.mean, d2.sem, d2.n, sigPhi));
    put("sigma_theta", fmt("sigma_theta (statistical error of A_eff over the fit tracks, std(A)/sqrt(N), AMean::StatErr, N = %.0f)\nvalue on A_eff = %.6f (A_eff = %.5f)\nvalue on P = %.5f\nmay overlap errFit, both come from the same tracks\n", N, B.a.StatErr(), A, sigTheta));
    // sigma_stat = errFit (+) sigma_phi (+) sigma_theta. delta_syst = delta_phi, delta_theta, the AP model and beam terms. The AP spread, t mismatch, cut variations and R stay in their own txt files, not in the sums.
    const double sigStat = std::sqrt(sigFit * sigFit + sigPhi * sigPhi + sigTheta * sigTheta);
    const std::vector<std::pair<std::string, double>> deltas = {{"delta_phi est n frame phi range, half the P span", dPhi}, {"delta_theta est n frame theta range", dRange}, {"delta_model A, k", dModel}, {"delta_beam beam momentum spread", dBeam}};
    double dSq = 0;
    for (const auto& d : deltas) dSq += d.second * d.second;
    const double dSyst = std::sqrt(dSq), sigTot = std::sqrt(sigStat * sigStat + dSq);
    const double pull = sigTot > 0 ? (P - pTrue) / sigTot : 0.0;
    put("bias_closure", fmt("bias_closure (corrected P against the true P of the run)\nP = %.5f, P_true = %.5f, P - P_true = %+.5f\nsigma_stat = %.5f, delta_syst = %.5f, total = %.5f\npull = (P - P_true) / total = %+.3f (one run, a pull needs many runs: error on the mean 1/sqrt(N runs), on the RMS 1/sqrt(2 N runs))\nP - P_true in sigma_stat only = %+.3f\n", P, pTrue, P - pTrue, sigStat, dSyst, sigTot, pull, sigStat > 0 ? (P - pTrue) / sigStat : 0.0));
    const double sigCon = std::sqrt(f0 * (1 - f0) / nCls[2]);
    // Ddelta_CEXS: sigma at the beam s over all 2^6 corners of the parameter errors, smallest and biggest sigma as ratios r to the nominal, the CE + CE mix count scaled by r, half the f1 span
    const double mN = 0.939565, mP = 0.938272, sBeam = mP * mP + 2 * (beamKE / 1000.0 + mN) * mP + mN * mN;
    const double sigNom = CEXSShape(kCEXSNom, sBeam);
    double rMin = 1, rMax = 1, pLo[6], pHi[6];
    std::copy(kCEXSNom, kCEXSNom + 6, pLo); std::copy(kCEXSNom, kCEXSNom + 6, pHi);
    for (int m = 0; m < 64; m++){
        double p[6];
        bool skip = false;
        for (int i = 0; i < 6; i++){
            if (!kCEXSVary[i] && ((m >> i) & 1)) skip = true;
            p[i] = kCEXSVary[i] ? kCEXSNom[i] + ((m >> i) & 1 ? kCEXSErr[i] : -kCEXSErr[i]) : kCEXSNom[i];
        }
        if (skip) continue;
        const double r = CEXSShape(p, sBeam) / sigNom;
        if (r < rMin){rMin = r; std::copy(p, p + 6, pLo);}
        if (r > rMax){rMax = r; std::copy(p, p + 6, pHi);}
    }
    const double nOther = nCls[2] - nCls[0];
    const auto f0At = [&](double r){return r * nCls[0] / (r * nCls[0] + nOther);};
    const double dCEXS = std::fabs(f0At(rMax) - f0At(rMin)) / 2, dNCE = nCls[0] * (rMax - rMin) / 2;
    put("delta_CEXS", fmt("delta_CEXS (np CE cross section A, alpha, c1 over their errors, all 8 corners, g0, g1, c0 held at nominal, beam KE %.0f MeV, s = %.4f GeV^2)\nsigma_min / sigma = %.4f, sigma_max / sigma = %.4f\nCE recon (CE only) N = %.0f, from %.0f to %.0f (+- %.0f)\nf0 = %.5f, at sigma_min %.5f, at sigma_max %.5f\nDdelta_CEXS = half the f0 span = %.5f\nparameters at sigma_min: %.5g %.5g %.5g %.5g %.5g %.5g\nparameters at sigma_max: %.5g %.5g %.5g %.5g %.5g %.5g\n", beamKE, sBeam, rMin, rMax, nCls[0], nCls[0] * rMin, nCls[0] * rMax, dNCE, f0, f0At(rMin), f0At(rMax), dCEXS, pLo[0], pLo[1], pLo[2], pLo[3], pLo[4], pLo[5], pHi[0], pHi[1], pHi[2], pHi[3], pHi[4], pHi[5]));
    const bool atBound = std::fabs(u.pRaw) > 0.9999;

    std::string sum = fmt("P_summary, run %.0f blocks, %.0f neutrons, true P = %.3f, beam KE %.0f MeV\nsample: total recon, 6_megaCut (KE %.0f%%, even planes, theta 0 to 1 degree, 1 hit per plane), N = %.0f fit tracks\n\n", settings.at("nBlocks"), settings.at("TotalNeutrons"), pTrue, beamKE, 100 * base.ke, N);
    sum += fmt("P = %.4f  +-  %.4f (stat)  +-  %.4f (syst)\n", P, sigStat, dSyst);
    if (atBound) sum += fmt("WARNING the fit is at the +-1 bound (P before correction = %.4f): the likelihood maximum is outside [-1, 1], so errFit and the bootstrap (%.4f) are not Gaussian errors, and the corrected P can exceed 1\n", u.pRaw, bs.sd);
    sum += "\nP FIT UNCERTAINTY\n";
    sum += fmt("1. sigma_P (stat) = %.4f = errFit %.4f (+) sigma_phi %.4f (+) sigma_AP %.4f; bootstrap %.4f\n", sigStat, sigFit, sigPhi, sigTheta, bs.sd);
    sum += fmt("   a. sigma_phi = %.5f (std(D_block)/sqrt(N) over the fit tracks, N = %.0f)\n", sigPhi, d2.n);
    sum += fmt("   b. sigma_AP = %.5f (sigma_t = sigma_theta, std(A)/sqrt(N) over the fit tracks, N = %.0f)\n", sigTheta, N);
    sum += fmt("2. delta_P (syst) = %.4f, each on P (quadrature):\n", dSyst);
    for (const auto& d : deltas) sum += fmt("   %.5f  ", d.second) + d.first + "\n";
    sum += "   not computed (need other runs) efficiency, seed road, plane weights, kTubeSigmas, H/C target mix\n";
    sum += "\nDILUTION (shown, not in the +- sums)\n";
    sum += fmt("1. Dcon: CE only fraction f0 = %.4f, CE + CE mix f1 = %.4f (not corrected)\n", f0, f1);
    sum += fmt("   a. Dsigma_con = %.4f (binomial error of f0)\n", sigCon);
    sum += fmt("   b. Ddelta_con = Ddelta_CEXS (A, alpha, c1) = %.4f (sigma x %.3f to x %.3f, CE only count %.0f +- %.0f)\n", dCEXS, rMin, rMax, nCls[0], dNCE);
    if (haveTrk){
        sum += fmt("2. Dtrack = %.4f (CE protons with a pure reconstructed track %.0f of %.0f with 2+ tracker hits, not corrected)\n", Dtrk, nPureTrueCE, nTrueTrCE);
        sum += fmt("   a. Dsigma_track = %.4f (binomial)\n", sigTrk);
    } else sum += "2. Dtrack: not computed (no plots/8_trackerfits/track_purity.csv, run F8_TrackerFits.cpp)\n";
    sum += "3. Dphi\n";
    sum += fmt("   Ddelta_phi block grid phi range D_block = %.4f, shift of P %+.4f (not applied)\n", Dblk, u.biasPhi2);
    sum += fmt("   Ddelta_nphi est n frame phi - true n frame phi, D_actual = %.4f +- %.4f (not applied)\n", Dact, mCos.sem);
    sum += fmt("   total applied dilution D = %.4f\n", u.dil);
    sum += fmt("4. DAP, Dtheta: P before any correction = %.4f\n", u.pRaw);
    sum += fmt("   Ddelta_ntheta est n frame theta - true n frame theta (smearBias) %+.5f, would shift P by %+.4f (not applied)\n", sm0.mean, P * (A / (A - sm0.mean) - 1.0));
    sum += "   Ddelta_theta block theta - true theta smear: not computed\n";
    sum += fmt("\nCLOSURE: P - P_true = %+.4f, pull = %+.3f (total error %.4f)\n", P - pTrue, pull, sigTot);
    put("P_summary", sum);

    // P_table.md, the same numbers as a markdown table
    // errors on A_eff first (the AP terms, A_eff units and relative), then the errors on P; "-" where a term acts on P only
    const double aSigTh = B.a.StatErr(), aRng = B.a.RangeErr(), aMod = std::fabs(A) * kModelAFrac, aBeam = B.a.BeamErr();
    const double aT = std::sqrt(aRng * aRng + aBeam * aBeam), pT = std::sqrt(dRange * dRange + dBeam * dBeam); // delta_t (syst): delta_theta and delta_beam act through t; the stat part sigma_t is sigma_theta
    const double aTot = std::sqrt(aT * aT + aMod * aMod);
    const double pAP = std::sqrt(pT * pT + dModel * dModel);
    const auto aCell = [&](double a){return std::string(Form("%.5f (%.2f%%)", a, 100.0 * a / std::fabs(A)));};
    // three sections (P fit uncertainty on A_eff and on P, Dilution, Bias), each a tree in a code block: a label with its branch characters and one value
    using Rows = std::vector<std::pair<std::string, std::string>>;
    const auto pv = [&](double v){return std::string(Form("%.4f", v));};
    const auto width = [](const std::string& s){size_t n = 0; for (unsigned char ch : s) if ((ch & 0xC0) != 0x80) n++; return n;}; // characters of a UTF-8 string
    const auto block = [&](const Rows& r){
        size_t w = 0;
        for (const auto& x : r) w = std::max(w, width(x.first));
        std::string s = "```\n";
        for (const auto& x : r) s += x.first + (x.second.empty() ? std::string() : std::string(w + 3 - width(x.first), ' ') + x.second) + "\n";
        return s + "```\n\n";
    };
    std::string md = std::string("# P error tree, ") + gSystem->BaseName(gSystem->pwd()) + fmt(" (%.0f neutrons, P_true %.2f, %.0f blocks, %.0f MeV)\n\n", settings.at("TotalNeutrons"), pTrue, settings.at("nBlocks"), beamKE);
    md += fmt("**P = %.3f ± %.3f (stat) ± %.3f (syst)**, pull %+.2f.", P, sigStat, dSyst, pull);
    if (atBound) md += " The raw fit sits at the +1 or -1 bound, so errFit and the total are not Gaussian at this N.";
    md += "\n\nA parent is the quadrature of its children.\n\n";

    md += fmt("## 1. P fit uncertainty\n\n### On P\n\nP = %.4f. Absolute error, the relative error in brackets.\n\n", P);
    const auto pCell = [&](double v){return std::string(Form("%.4f (%.2f%%)", v, 100.0 * v / std::fabs(P)));};
    md += block({{"σ_P (stat)  [on P]", pCell(sigStat)},
                 {"├─ errFit (likelihood)", pCell(sigFit)},
                 {"├─ σ_φ  stat of φ over the tracks", pCell(sigPhi)},
                 {"└─ σ_AP  stat of θ over the tracks (σ_t)", pCell(sigTheta)},
                 {"δ_P (syst)  [on P]", pCell(dSyst)},
                 {"├─ δ_φ  est n frame φ range", pCell(dPhi)},
                 {"└─ δ_AP  [on P, |P| × the relative A_eff error]", pCell(pAP)}});
    md += fmt("### On A_eff\n\nA_eff = %.5f ± %.5f (stat) ± %.5f (syst). Absolute error, the relative error in brackets.\n\n", A, aSigTh, aTot);
    md += block({{"AP  [on A_eff]", fmt("± %.5f (stat) ± %.5f (syst)", aSigTh, aTot)},
                 {"├─ σ_AP (stat)", aCell(aSigTh)},
                 {"│  └─ σ_t", aCell(aSigTh)},
                 {"│     └─ σ_θ  stat of θ over the tracks", aCell(aSigTh)},
                 {"└─ δ_AP (syst)", aCell(aTot)},
                 {"   ├─ δ_t", aCell(aT)},
                 {"   │  ├─ δ_θ  est n frame θ range", aCell(aRng)},
                 {"   │  └─ δ_beam  beam spread", aCell(aBeam)},
                 {"   └─ δ_model  A, k", aCell(aMod)}});
    md += "- δ_AP and its terms are the same in both trees: P carries the relative A_eff error. σ_AP is σ_t, the stat of θ over the tracks, std(A)/√N, and σ_φ is std(D_block)/√N; σ_AP may overlap errFit, both come from the same tracks. δ_AP is the syst side, δ_t and δ_model in quadrature.\n";
    md += "- δ_φ is half the P span when every track's est n frame φ moves by −h and +h (h the φ range half-width), refit. δ_t is δ_θ and δ_beam in quadrature (both act through t). δ_θ is half the A_eff span of the θ range, δ_model the A, k envelope (no covariance), δ_beam the ±1σ beam momentum spread propagated.\n";
    md += "- Not computed: efficiency, seed road, plane weights, kTubeSigmas, H/C target mix.\n\n";

    md += "## 2. Dilution (not in the ± sums)\n\nThe relative error in brackets: Dσ_con and Dδ_con against f0, Dσ_track against |Dtrack|, Dδ_θ against |A_eff|.\n\n";
    const auto rCell = [&](double v, double base, const char* f){return std::string(Form(f, v)) + Form(" (%.2f%%)", 100.0 * std::fabs(v) / std::fabs(base));};
    md += block({{fmt("Dcon  f0 = %.3f (CE only), f1 = %.3f (CE + CE mix), not applied", f0, f1), ""},
                 {"├─ Dσ_con  binomial error of f0", rCell(sigCon, f0, "%.4f")},
                 {"└─ Dδ_con", rCell(dCEXS, f0, "%.4f")},
                 {"   └─ Dδ_CEXS  A, α, c1", rCell(dCEXS, f0, "%.4f")},
                 {"Dtrack  not applied", haveTrk ? pv(Dtrk) : std::string("not computed")},
                 {"└─ Dσ_track", haveTrk ? rCell(sigTrk, Dtrk, "%.4f") : std::string("not computed")},
                 {"Dφ", ""},
                 {"└─ Dδ_φ  block grid φ range D_block, not applied", pv(Dblk)},
                 {"DAP  not applied", ""},
                 {"└─ Dθ", ""},
                 {"   └─ Dδ_θ  θ_block − θ_true smear (smearBias)", rCell(sm0.mean, A, "%+.5f")}});
    md += haveTrk ? fmt("- Dtrack = CE child protons with a pure reconstructed track (all hits one track ID) over the CE child protons with 2 or more tracker hits, no cuts (F8_TrackerFits.cpp), %.0f of %.0f; Dσ_track is binomial.\n", nPureTrueCE, nTrueTrCE) : std::string("- Dtrack not computed: run F8_TrackerFits.cpp for plots/8_trackerfits/track_purity.csv.\n");
    md += fmt("- Dδ_CEXS is half the f0 span when the CE only count is scaled by σ_min/σ = %.3f and σ_max/σ = %.3f, the np CE cross section over all 8 corners of A, α, c1 at their errors (g0, g1, c0 held at nominal), at the beam s. CE only recon N = %.0f ± %.0f.\n", rMin, rMax, nCls[0], dNCE);
    md += fmt("- Dδ_φ would shift P by %+.4f; it is measured and not applied (kApplyBlockPhiDilution in FitCore.h).\n", u.biasPhi2);

    md += fmt("- Dδ_θ would shift P by %+.4f; it is measured and not applied (kApplySmearBias in FitCore.h).\n", P * (A / (A - sm0.mean) - 1.0));
    { std::ofstream o(dir + "/P_table.md"); o << md; }
    // P_result.csv: every number of the txt files; the rows of the totals first, then the variations and the not computed items collected above. F9_PlotReconNew.cpp reads the item P rows (P +- stat +- syst on the total recon fit plot of 6_megaCut)
    const std::string earlyRows = csvRows;
    csvRows.clear();
    row("P", "P", P); row("P", "sigma_stat", sigStat); row("P", "delta_syst", dSyst); row("P", "sigma_total", sigTot); row("P", "P_true", pTrue); row("P", "P - P_true", P - pTrue); row("P", "pull", pull); row("P", "P_before_correction", u.pRaw);
    row("P_stat", "errFit", sigFit); row("P_stat", "sigma_phi", sigPhi); row("P_stat", "sigma_AP", sigTheta); row("P_stat", "bootstrap sd", bs.sd, "200 resamples");
    row("P_syst", "delta_phi", dPhi); row("P_syst", "delta_AP", pAP); row("P_syst", "delta_t", pT); row("P_syst", "delta_theta", dRange); row("P_syst", "delta_beam", dBeam); row("P_syst", "delta_model", dModel);
    row("A_eff", "A_eff", A); row("A_eff", "sigma_AP", aSigTh); row("A_eff", "sigma_t", aSigTh); row("A_eff", "sigma_theta", aSigTh); row("A_eff", "delta_AP", aTot); row("A_eff", "delta_t", aT); row("A_eff", "delta_theta", aRng); row("A_eff", "delta_beam", aBeam); row("A_eff", "delta_model", aMod);
    row("sigma_P_fit", "errFit", sigFit); row("sigma_P_fit", "N", N); row("sigma_P_fit", "check sqrt(2/N)/|A_eff|", std::sqrt(2.0 / N) / std::fabs(A));
    row("sigma_theta", "value on A_eff", B.a.StatErr()); row("sigma_theta", "value on P", sigTheta);
    row("delta_AP_model", "relative on A_eff", kModelAFrac); row("delta_AP_model", "value on P", dModel);
    row("delta_AP_beam", "beam momentum spread MeV/c", kBeamMomentumSD); row("delta_AP_beam", "beam momentum MeV/c", kBeamMomentum); row("delta_AP_beam", "relative on A_eff", kBeamAFrac); row("delta_AP_beam", "value on P", dBeam);
    row("delta_AP_range", "value on A_eff", B.a.RangeErr()); row("delta_AP_range", "value on P", dRange);
    row("delta_AP_tmismatch_block", "CE recon rows", mAB.n); rowE("delta_AP_tmismatch_block", "A_eff block theta", mAB.mean, mAB.sem); rowE("delta_AP_tmismatch_block", "A_eff model (LastAP)", mAM.mean, mAM.sem);
    row("delta_AP_tmismatch_block", "ratio A_model / A_block", rTB); row("delta_AP_tmismatch_block", "t-only ratio", rT); row("delta_AP_tmismatch_block", "value on P", std::fabs(P) * std::fabs(rTB - 1.0));
    row("delta_AP_tmismatch", "CE recon rows", mAC.n); rowE("delta_AP_tmismatch", "A_eff analysis", mAC.mean, mAC.sem); rowE("delta_AP_tmismatch", "A_eff model (LastAP)", mAM.mean, mAM.sem);
    row("delta_AP_tmismatch", "mean t model", mTM.mean); row("delta_AP_tmismatch", "mean t analysis", mTC.mean); row("delta_AP_tmismatch", "ratio A_model / A_analysis", rT); row("delta_AP_tmismatch", "value on P", dTmis);
    row("delta_AP_tmismatch", "error on the ratio", sigRatio); row("delta_AP_tmismatch", "error on the ratio on P", std::fabs(P) * sigRatio);
    row("delta_phi", "P at -h", pPhiLo.p); row("delta_phi", "P at +h", pPhiHi.p); row("delta_phi", "value on P", dPhi); row("delta_phi", "sigma_phi on P", sigPhi); rowE("delta_phi", "D_block", d2.mean, d2.sem); row("delta_phi", "N", d2.n);
    rowE("D_block", "D", Dblk, d2.sem, "block grid phi range dilution; not applied"); row("D_block", "shift of P", u.biasPhi2);
    if (haveTrk){
        rowE("D_track", "D CE protons", Dtrk, sigTrk, "binomial; pure track over CE protons with 2+ tracker hits"); row("D_track", "CE protons with a pure track", nPureTrueCE); row("D_track", "CE protons with 2+ tracker hits", nTrueTrCE);
        row("D_track", "D all particles", nTrueTr > 0 ? nPureTrue / nTrueTr : 0.0); row("D_track", "particles with a pure track", nPureTrue); row("D_track", "particles with 2+ tracker hits", nTrueTr);
    } else rowNote("D_track", "not computed: no plots/8_trackerfits/track_purity.csv");
    row("delta_CEXS", "beam KE MeV", beamKE); row("delta_CEXS", "s GeV^2", sBeam); row("delta_CEXS", "sigma_min / sigma", rMin); row("delta_CEXS", "sigma_max / sigma", rMax);
    row("delta_CEXS", "CE only N", nCls[0]); row("delta_CEXS", "N at sigma_min", nCls[0] * rMin); row("delta_CEXS", "N at sigma_max", nCls[0] * rMax); row("delta_CEXS", "N error", dNCE);
    row("delta_CEXS", "f0", f0); row("delta_CEXS", "f0 at sigma_min", f0At(rMin)); row("delta_CEXS", "f0 at sigma_max", f0At(rMax)); row("delta_CEXS", "value (half the f0 span)", dCEXS);
    const char* cexsName[6] = {"A", "alpha", "g0", "g1", "c0", "c1"};
    for (int i = 0; i < 6; i++){row("delta_CEXS", std::string("param at sigma_min ") + cexsName[i], pLo[i]); row("delta_CEXS", std::string("param at sigma_max ") + cexsName[i], pHi[i]);}
    row("Dcon", "f0", f0, "CE only fraction"); row("Dcon", "f1", f1, "CE + CE mix fraction"); row("Dcon", "Dsigma_con", sigCon, "binomial error of f0"); row("Dcon", "Ddelta_con", dCEXS);
    row("Dphi", "Ddelta_phi D_block", Dblk, "not applied"); row("Dphi", "Ddelta_phi shift of P", u.biasPhi2); rowE("Dphi", "Ddelta_nphi D_actual", Dact, mCos.sem, "not applied"); row("Dphi", "total applied D", u.dil);
    row("Dtheta", "smearBias", sm0.mean, "A_eff shift of the recon theta smear; not applied"); row("Dtheta", "would shift P by", P * (A / (A - sm0.mean) - 1.0));
    { std::ofstream o(dir + "/P_result.csv"); o << "item,quantity,value,error,note\n" << csvRows << earlyRows; }
    printf("%s", sum.c_str());
    printf(">>> wrote %s/*.txt and P_table.md\n", dir.c_str());
}
