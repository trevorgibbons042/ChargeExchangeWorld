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

// P error budget of one run: the statistical, systematic, bias and dilution items of /Users/trevorg04/Desktop/P_fit_error_sources.md that one sim run can give, written to
// plots/11_errorBudget/: P_result.csv (item,quantity,value,error,note, one row per number: P = P_fit +- sigma_stat +- delta_syst, sigma_*, delta_*, D_*, closure) and P_table.md, the same as a tree. The summary text is printed to the terminal only.
// The sample is the total recon fit of F9_PlotReconNew.cpp's 6_megaCut (every recon track after the KE 90%, theta 0 to 1 degree, 1 hit per plane, at least 3 planes and exactly 1 recon track in the event cuts, the block phi fit in the est neutron frame),
// so P_fit and errFit here equal pAll and errAll of plots/6_ReconNew_1/6_megaCut/scan_stats.csv. The variations (cuts, fit range) are refits of the same rows.
// Items that need another sim or another ReconAngleAll run (efficiency, P = 0, seed road, plane weights, kTubeSigmas, H/C target mix) are written as "not computed".
namespace {
constexpr double kRad = 180.0 / M_PI;
struct Row {double thMin, thMax, phHalf, th2Min, th2Max, ph2Half, thetaRelND, phiRelND, zdcKEAllD; int eventID, evClass, nPlanesUsed, nPlanesHit, nTrackerHits, nRecoEv;};
struct Cuts {double ke = 0.90, thLo = 0.0, thHi = 1.0, phiTrim = 0.0; bool even = false, avgOne = true;}; // even: the even planes cut, not in the megaCut (F9)
 // ke: fraction of the beam KE, phiTrim: degrees cut off both ends of the phi fit range
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
    // AnglesRaw and cut_CE are opened without the Events friend: every slot of an RDataFrame opens its own copy of the tree, and the friend's index on fEvent costs about 0.9 GB per copy at 100M events,
    // so the RAM grew with the threads (killed at 64 and 128). The Events columns used here are read once, serially, into a table sorted by fEvent that every slot shares (CalcAPandT.cpp does the same).
    TFile* anglesFile = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree* t = anglesFile && !anglesFile->IsZombie() ? (TTree*)anglesFile->Get("AnglesRaw") : nullptr;
    if (!t){Error("PErrorBudget", "no cuts/cut_angles_all.root (run ReconAngleAll.cpp)"); return;}
    struct EvRow {int fEvent; double ke, planesHit, trackerHits, nReco;};
    std::vector<EvRow> evTable; // sorted by fEvent
    double tally[4] = {0, 0, 0, 0}; // sums over the events of nPureTrueEv, nTrueEv, nPureCEEv, nTrueCEEv (Dtrack, below), no cuts
    { // serial, before ROOT::EnableImplicitMT: with it on, every GetEntry of the 38M events farms its decompression out to the pool and the loop does not finish
        TFile evFile("cuts/cut_events.root", "READ");
        TTree* ev = evFile.IsZombie() ? nullptr : (TTree*)evFile.Get("Events");
        if (!ev){Error("PErrorBudget", "no cuts/cut_events.root (rerun ReconAngleAll.cpp)"); return;}
        const char* names[] = {"fEvent", "zdcKEAllD", "nPlanesHitEv", "nTrackerHitsEv", "nRecoEv", "nPureTrueEv", "nTrueEv", "nPureCEEv", "nTrueCEEv"};
        TLeaf* leaf[9];
        ev->SetBranchStatus("*", 0);
        for (int i = 0; i < 9; i++){
            ev->SetBranchStatus(names[i], 1); leaf[i] = ev->GetLeaf(names[i]);
            if (!leaf[i]){Error("PErrorBudget", "no %s in the Events tree (rerun ReconAngleAll.cpp)", names[i]); return;}
        }
        const Long64_t nEv = ev->GetEntries();
        evTable.reserve(nEv);
        for (Long64_t i = 0; i < nEv; i++){
            ev->GetEntry(i);
            evTable.push_back({(int)leaf[0]->GetValue(), leaf[1]->GetValue(), leaf[2]->GetValue(), leaf[3]->GetValue(), leaf[4]->GetValue()});
            for (int k = 0; k < 4; k++) tally[k] += leaf[5 + k]->GetValue();
        }
        std::sort(evTable.begin(), evTable.end(), [](const EvRow& a, const EvRow& b){return a.fEvent < b.fEvent;});
    }
    ROOT::EnableImplicitMT(NThreads());
    // one pass that fills each Row by its entry number (no column vectors next to rows, and the rows in the order of AnglesRaw)
    std::vector<Row> rows((size_t)t->GetEntries());
    ROOT::RDataFrame dfRaw(*t);
    auto dfRows = WithEntryNumber(dfRaw, dfRaw.GetNSlots()).Define("__ev", [&evTable](int e){ // the Events columns of the row's event: zdcKEAllD, nPlanesHitEv, nTrackerHitsEv, nRecoEv
        const auto it = std::lower_bound(evTable.begin(), evTable.end(), e, [](const EvRow& r, int ev){return r.fEvent < ev;});
        return it != evTable.end() && it->fEvent == e ? ROOT::RVecD{it->ke, it->planesHit, it->trackerHits, it->nReco} : ROOT::RVecD{0., 0., 0., 0.};
    }, {"fEvent"}).Define("__row", "ROOT::RVecD{(double)thetaRelNDMin, (double)thetaRelNDMax, (double)phiRelNDHalf, (double)thetaRelNDBMin, (double)thetaRelNDBMax, (double)phiRelNDBHalf, (double)thetaRelND, (double)phiRelND, __ev[0], (double)fEvent, (double)evClass, (double)nPlanesUsed, __ev[1], __ev[2], __ev[3]}");
    dfRows.Foreach([&](ULong64_t entry, const ROOT::RVecD& v){
        size_t k = 0;
        Row& x = rows[(size_t)entry];
        x.thMin = v[k++] * kRad; x.thMax = v[k++] * kRad; x.phHalf = v[k++]; x.th2Min = v[k++] * kRad; x.th2Max = v[k++] * kRad; x.ph2Half = v[k++];
        x.thetaRelND = v[k++]; x.phiRelND = v[k++]; x.zdcKEAllD = v[k++];
        x.eventID = (int)v[k++]; x.evClass = (int)v[k++]; x.nPlanesUsed = (int)v[k++]; x.nPlanesHit = (int)v[k++]; x.nTrackerHits = (int)v[k++]; x.nRecoEv = (int)v[k++]; // nRecoEv: the fits (rows) of the event, from the Events tree
    }, {"entryNo", "__row"});
    const auto settings = ReadSettings();
    const double apScale = settings.at("APScale"), beamKE = settings.at("BeamKE"), pTrue = settings.at("PPValue");
    const int maxPlanes = MaxTrackerHits();

    // the true CE of each event (ChargeExchange_values, the carbon CE rows): true phi and theta in the true neutron frame (deg), LastAP, LastT (model t), LastTheta (lab theta)
    TFile* ceFile = TFile::Open("cuts/cut_CE.root", "READ"); // no Events friend, see above
    TTree* ceT = ceFile && !ceFile->IsZombie() ? (TTree*)ceFile->Get("ChargeExchange_values") : nullptr;
    gROOT->cd();
    if (!ceT){Error("PErrorBudget", "no cuts/cut_CE.root"); return;}
    const auto ce = Columns(ceT, {"fEvent", "CEPhiNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "CEThetaNeu(LastTheta, LastPhiLab, CENDirX, CENDirY, CENDirZ) * 180 / TMath::Pi()", "LastAP", "LastT", "LastTheta"});
    std::map<long long, size_t> byEvent;
    for (size_t i = 0; i < ce[0].size(); i++) byEvent.emplace((long long)ce[0][i], i);

    // select() runs 8 times (the baseline and 7 variations). The rows that pass the loosest cut set of those (KE 85%, theta to 1.2 degree, no avgOne, phi untrimmed) are listed once in the order of rows,
    // so each call scans that list and not every row. Keep it at least as loose as every Cuts used below.
    EventCut looseCut; looseCut.ke = 0.85 * beamKE; looseCut.theta = true; looseCut.thHi = 1.2;
    std::vector<const Row*> loose;
    for (const Row& x : rows){
        if (x.evClass < 0 || x.evClass > 2 || x.nPlanesUsed < 3 || x.nRecoEv != 1) continue;
        if (!PassEventCut(looseCut, x.zdcKEAllD, x.nPlanesUsed, x.nRecoEv, x.thetaRelND * kRad, x.nPlanesHit, x.nTrackerHits)) continue;
        loose.push_back(&x);
    }
    // the rows of a cut set (F9_PlotReconNew.cpp passes + inFit): every recon track, evClass 0 to 2
    const auto select = [&](const Cuts& c){
        std::vector<Sel> s;
        EventCut ec; ec.ke = c.ke * beamKE; ec.even = c.even; ec.maxPlanes = maxPlanes; ec.theta = true; ec.thLo = c.thLo; ec.thHi = c.thHi; ec.avgOne = c.avgOne; ec.firstTrack = true;
        for (const Row* xp : loose){
            const Row& x = *xp;
            if (x.evClass < 0 || x.evClass > 2) continue;
            const double th = x.thetaRelND * kRad;
            if (!PassEventCut(ec, x.zdcKEAllD, x.nPlanesUsed, x.nRecoEv, th, x.nPlanesHit, x.nTrackerHits)) continue; // F9 megaCut: no 2 plane tracks, only events with 1 recon track
            const double ph = x.phiRelND * kRad;
            if (!(ph > kBlockPhiLo + c.phiTrim && ph < kBlockPhiHi - c.phiTrim)) continue;
            Sel q;
            const AEffIn e = AEffInOf(th, x.thMin, x.thMax, x.th2Min, x.th2Max, apScale);
            q.ph = ph; q.thDeg = th; q.a = e.a; q.r1 = e.r1; q.r2 = e.r2; q.a1lo = e.a1lo; q.a1hi = e.a1hi; q.a2lo = e.a2lo; q.a2hi = e.a2hi;
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
    const auto fmt = [](const char* f, auto... a){return std::string(Form(f, a...));}; // every argument is a double
    // P_result.csv: item,quantity,value,error,note, one row per number of the budget (commas in a note become semicolons, error empty when none)
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
    // bootstrap of the fit rows
    const int nBoot = 200;
    std::vector<double> bp;
    {
        // the 200 resamples are independent (each has its own seed, 20260101 + b) and UnbinnedP(..., parallel = false) is serial, so they run on all the threads; the results come back in b order
        std::vector<unsigned> ids(nBoot);
        std::iota(ids.begin(), ids.end(), 0u);
        ROOT::TThreadExecutor ex(NThreads());
        const auto boots = ex.Map([&](unsigned b){
            TRandom3 rng(20260101 + b);
            const size_t np = s0.size();
            std::vector<double> ph(np);
            AMean a; PhiRange q1, q2;
            for (size_t i = 0; i < np; i++){
                const Sel& q = s0[(size_t)(rng.Rndm() * np)];
                ph[i] = q.ph; a.Add(q.a, q.a1lo, q.a1hi, q.a2lo, q.a2hi);
                if (q.r1) q1.Add(q.h1);
                if (q.r2) q2.Add(q.h2);
            }
            a.smearBias = sm0.mean;
            return UnbinnedP(ph, a, kBlockPhiLo, kBlockPhiHi, nullptr, false, &q1, &q2);
        }, ids);
        for (const PFit& v : boots) if (v.ok) bp.push_back(v.p);
    }
    const MeanErr bs = Stats(bp);
    // purity
    double nCls[3] = {0, 0, 0};
    for (const Sel& q : s0) for (int k = q.cls; k < 3; k++) nCls[k]++;
    const double f0 = nCls[0] / nCls[2], f1 = nCls[1] / nCls[2];

    // ---------- systematic ----------
    const double dModel = u.errAModel, dBeam = u.errABeam, dRange = u.errARange;
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
    // cut and fit range variations: refit the same rows with one setting changed. sd = sqrt(|err_var^2 - err_base^2|) is the statistical size of the shift (nested samples); delta = the largest sqrt(max(0, shift^2 - sd^2)), the shift beyond statistics
    struct Var {const char* label; Cuts c;};
    const auto variation = [&](const char* name, const char* what, std::vector<Var> vars){
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
        }
        row(name, "value on P", dMax, "largest shift beyond statistics");
        return dMax;
    };
    Cuts c;
    c = base; c.ke = 0.85; Cuts ke85 = c; c = base; c.ke = 0.95; Cuts ke95 = c;
    const double dKE = variation("delta_KEcut", "KE cut 85% and 95% of the beam KE instead of 90%", {{"KE 85%", ke85}, {"KE 95%", ke95}});
    c = base; c.thHi = 0.8; Cuts th08 = c; c = base; c.thHi = 1.2; Cuts th12 = c;
    const double dTh = variation("delta_theta_cut", "theta cut 0 to 0.8 and 0 to 1.2 degree instead of 0 to 1", {{"theta 0 to 0.8", th08}, {"theta 0 to 1.2", th12}});
    c = base; c.even = true; Cuts noEven = c;
    const double dPl = variation("delta_planes_cut", "even planes cut added", {{"even planes cut", noEven}});
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
    for (const auto& nc : notComputed) rowNote(nc[0], std::string("not computed: ") + nc[1]);

    // ---------- bias and dilution ----------
    const double Dblk = B.q2.D();
    std::vector<double> d2v;
    for (const Sel& q : s0) d2v.push_back(Dof(q.r2, q.h2));
    const MeanErr d2 = Stats(d2v);
    // actual dilution: CE recon rows, est neutron phi against the true neutron frame phi of the CE
    std::vector<double> cosD;
    for (const Sel& q : s0) if (q.cls == 0 && q.ceIdx >= 0 && q.r1) cosD.push_back(std::cos((q.ph - ce[1][q.ceIdx]) / kRad));
    const MeanErr mCos = Stats(cosD);

    // ---------- totals, in the structure of the error list: P fit uncertainty (sigma_P, delta_P), Dilution (D), Bias ----------
    // Dtrack: the true tracks (2+ tracker hits) that sit in a pure reconstructed track over all true tracks with 2+ tracker hits, no cuts (the tally columns of the Events tree, ReconAngleAll.cpp, made from the sim output)
    // Dtrack is the CE child protons only (nPureTrueCE of nTrueTracksCE); the all particle ratio (nPureTrue of nTrueTr) is only in the D_track rows of P_result.csv
    const double nPureTrue = tally[0], nTrueTr = tally[1], nPureTrueCE = tally[2], nTrueTrCE = tally[3]; // the sums of the tally columns of the Events tree (ReconAngleAll.cpp), no cuts, read with the table above
    const bool haveTrk = nTrueTrCE > 0;
    const double Dtrk = haveTrk ? nPureTrueCE / nTrueTrCE : 0.0, sigTrk = haveTrk ? std::sqrt(Dtrk * (1 - Dtrk) / nTrueTrCE) : 0.0; // binomial error
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
    // sigma_stat = errFit (+) sigma_phi (+) sigma_theta. delta_syst = delta_phi, delta_theta, the AP model and beam terms. The AP spread, t mismatch, cut variations and R stay in their own rows of P_result.csv, not in the sums.
    const double sigStat = std::sqrt(sigFit * sigFit + sigPhi * sigPhi + sigTheta * sigTheta);
    const std::vector<std::pair<std::string, double>> deltas = {{"delta_phi est n frame phi range, half the P span", dPhi}, {"delta_theta est n frame theta range", dRange}, {"delta_model A, k", dModel}, {"delta_beam beam momentum spread", dBeam}};
    double dSq = 0;
    for (const auto& d : deltas) dSq += d.second * d.second;
    const double dSyst = std::sqrt(dSq), sigTot = std::sqrt(sigStat * sigStat + dSq);
    const double pull = sigTot > 0 ? (P - pTrue) / sigTot : 0.0;
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
    const bool atBound = std::fabs(u.pRaw) > 0.9999;

    std::string sum = fmt("P_summary, run %.0f blocks, %.0f neutrons, true P = %.3f, beam KE %.0f MeV\nsample: total recon, 6_megaCut (KE %.0f%%, theta 0 to 1 degree, 1 hit per plane, 3+ planes, 1 recon track per event), N = %.0f fit tracks\n\n", settings.at("nBlocks"), settings.at("TotalNeutrons"), pTrue, beamKE, 100 * base.ke, N);
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
    } else sum += "2. Dtrack: not computed (no true tracks with 2+ tracker hits in cuts/cut_events.root, rerun ReconAngleAll.cpp)\n";
    sum += "3. Dphi\n";
    sum += fmt("   Ddelta_phi block grid phi range D_block = %.4f, shift of P %+.4f (not applied)\n", Dblk, u.biasPhi2);
    sum += fmt("   Ddelta_nphi est n frame phi - true n frame phi, D_actual = %.4f +- %.4f (not applied)\n", Dact, mCos.sem);
    sum += fmt("   total applied dilution D = %.4f\n", u.dil);
    sum += fmt("4. DAP, Dtheta: P before any correction = %.4f\n", u.pRaw);
    sum += fmt("   Ddelta_ntheta est n frame theta - true n frame theta (smearBias) %+.5f, would shift P by %+.4f (not applied)\n", sm0.mean, P * (A / (A - sm0.mean) - 1.0));
    sum += "   Ddelta_theta block theta - true theta smear: not computed\n";
    sum += fmt("\nCLOSURE: P - P_true = %+.4f, pull = %+.3f (total error %.4f)\n", P - pTrue, pull, sigTot);

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
    // P with the dilutions fixed: / f0 (Dcon); then also / Dtrack, / D_block (Dphi) and the A_eff smear shift A / (A - smearBias) (DAP). stat and syst scale with P; the D errors (Dsigma_con, Ddelta_con, Dsigma_track, relative) are added to syst in quadrature
    const double fCon = f0 > 0 ? 1.0 / f0 : 1.0;
    const double rCon = f0 > 0 ? std::sqrt(sigCon * sigCon + dCEXS * dCEXS) / f0 : 0.0;
    const bool useTrk = haveTrk && Dtrk > 0;
    const double fAll = fCon * (useTrk ? 1.0 / Dtrk : 1.0) / (Dblk > 0 ? Dblk : 1.0) * (A != sm0.mean ? A / (A - sm0.mean) : 1.0);
    const double rAll = std::sqrt(rCon * rCon + (useTrk ? (sigTrk / Dtrk) * (sigTrk / Dtrk) : 0.0));
    const double pCon = P * fCon, sCon = sigStat * fCon, yCon = std::sqrt(dSyst * dSyst * fCon * fCon + pCon * pCon * rCon * rCon);
    const double pFix = P * fAll, sFix = sigStat * fAll, yFix = std::sqrt(dSyst * dSyst * fAll * fAll + pFix * pFix * rAll * rAll);
    const auto pullOf = [&](double p, double s, double y){return s > 0 || y > 0 ? (p - pTrue) / std::sqrt(s * s + y * y) : 0.0;};
    md += "\n\n";
    md += fmt("With Dcon fixed (P / f0, f0 = %.3f): **P = %.3f ± %.3f (stat) ± %.3f (syst)**, pull %+.2f.\n\n", f0, pCon, sCon, yCon, pullOf(pCon, sCon, yCon));
    md += fmt("With Dcon, Dtrack, Dφ and DAP fixed (P / (f0 · Dtrack · D_block) and the A_eff smear shift, factor %.3f on P): **P = %.3f ± %.3f (stat) ± %.3f (syst)**, pull %+.2f.", fAll, pFix, sFix, yFix, pullOf(pFix, sFix, yFix));
    if (!useTrk) md += " Dtrack is not computed and left out.";
    md += " The Dσ_con, Dδ_con and Dσ_track errors are in the syst.";
    md += "\n\nA parent is the quadrature of its children.\n\n";

    md += fmt("## 1. P fit uncertainty\n\n### On P\n\nP = %.4f. Absolute error, the relative error in brackets.\n\n", P);
    const auto pCell = [&](double v){return std::string(Form("%.4f (%.2f%%)", v, 100.0 * v / std::fabs(P)));};
    md += block({{"σ_P (stat)  [on P]", pCell(sigStat)},
                 {"├─ errFit (likelihood)", pCell(sigFit)},
                 {"├─ σ_D  stat of the block φ dilution D over the tracks", pCell(sigPhi)},
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
    md += haveTrk ? fmt("- Dtrack = CE child protons with a pure reconstructed track (all hits one track ID) over the CE child protons with 2 or more tracker hits, no cuts (F8_TrackerFits.cpp), %.0f of %.0f; Dσ_track is binomial.\n", nPureTrueCE, nTrueTrCE) : std::string("- Dtrack not computed: no true tracks with 2+ tracker hits in cuts/cut_events.root, rerun ReconAngleAll.cpp.\n");
    md += fmt("- Dδ_CEXS is half the f0 span when the CE only count is scaled by σ_min/σ = %.3f and σ_max/σ = %.3f, the np CE cross section over all 8 corners of A, α, c1 at their errors (g0, g1, c0 held at nominal), at the beam s. CE only recon N = %.0f ± %.0f.\n", rMin, rMax, nCls[0], dNCE);
    md += fmt("- Dδ_φ would shift P by %+.4f; it is measured and not applied (kApplyBlockPhiDilution in FitCore.h).\n", u.biasPhi2);

    md += fmt("- Dδ_θ would shift P by %+.4f; it is measured and not applied (kApplySmearBias in FitCore.h).\n", P * (A / (A - sm0.mean) - 1.0));
    { std::ofstream o(dir + "/P_table.md"); o << md; }
    // P_result.csv: every number of the budget; the rows of the totals first, then the variations and the not computed items collected above. F9_PlotReconNew.cpp reads the item P rows (P +- stat +- syst on the total recon fit plot of 6_megaCut)
    const std::string earlyRows = csvRows;
    csvRows.clear();
    row("P", "P", P); row("P", "sigma_stat", sigStat); row("P", "delta_syst", dSyst); row("P", "sigma_total", sigTot); row("P", "P_true", pTrue); row("P", "P - P_true", P - pTrue); row("P", "pull", pull); row("P", "P_before_correction", u.pRaw);
    row("P_Dcon_fixed", "P", pCon, "P / f0"); row("P_Dcon_fixed", "sigma_stat", sCon); row("P_Dcon_fixed", "delta_syst", yCon, "includes Dsigma_con and Ddelta_con"); row("P_Dcon_fixed", "pull", pullOf(pCon, sCon, yCon));
    row("P_all_fixed", "P", pFix, "P / (f0 * Dtrack * D_block) and the A_eff smear shift"); row("P_all_fixed", "factor on P", fAll); row("P_all_fixed", "sigma_stat", sFix); row("P_all_fixed", "delta_syst", yFix, "includes Dsigma_con; Ddelta_con and Dsigma_track"); row("P_all_fixed", "pull", pullOf(pFix, sFix, yFix));
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
    } else rowNote("D_track", "not computed: no true tracks with 2+ tracker hits in cuts/cut_events.root");
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
    printf(">>> wrote %s/P_table.md and P_result.csv\n", dir.c_str());
}
