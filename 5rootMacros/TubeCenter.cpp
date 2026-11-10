#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include "FitCore.h"
#include "CutOpen.h"
#include "TreeColumns.h"

// The tube centre of every reconstructed particle of each 6_cuts sample that has an AnglesRaw cut file (all rows of cuts/<file>.root, which carry the tube end points: CutPlots.cpp (samples), rawCols):
// the mean of the two tube end points (tubeX0/tubeX1, tubeY0/tubeY1, tubeZ0/tubeZ1), the point SmearMatrix.cpp takes as the estimated neutron direction from the beam origin. Four pages in the
// 1_additionalinfo of each folder:
//   15 tube centre y vs x, 16 r = sqrt(x^2 + y^2) vs the theta of the estimated neutron direction, 17 r vs its phi, 18 r vs the tube centre z.
// Est neutron direction: theta = acos(z / |c|), phi = atan2(y, x) (0 to 360 deg) of the tube centre c. Rows without a tube (NaN) are left out.
// Only the rows of the cut are read (one multithreaded pass of each cut file), never the whole of AnglesRaw. TUBE_CUT = one folder of the list (default: all of them), smear_cuts.sh runs it.
struct TubeSpec {const char* dir; const char* file; int kePct; const char* labelDir;}; // kePct replaces the 80% of the cuts.txt of labelDir; file "" = no cut file of its own: the reconstructed ZDC protons (cut_ZdcAllEv.root) with ZDC KE above kePct % of the beam KE (0: all of them), their rows of cut_angles_all.root
static const std::vector<TubeSpec>& TubeSpecs(){
    static const std::vector<TubeSpec> specs = {
        {"1_hitZDC", "", 0, ""},
        {"2_KECUT/KECUT80", "", 80, ""},
        {"2_KECUT/KECUT85", "", 85, ""},
        {"2_KECUT/KECUT90", "", 90, ""},
        {"2_KECUT/KECUT95", "", 95, ""},
        {"3_ResCuts", "cut_ResCuts", 0, "3_ResCuts"},
        {"4_nZdcParticles/1_1zdcparticlecut", "cut_nZdc1", 0, "4_nZdcParticles/1_1zdcparticlecut"},
        {"6_megaCut3/1_1particlecut", "cut_megacut3_1", 0, "6_megaCut3/1_1particlecut"},
        {"6_megaCut3/1_1particlecut_KE85", "cut_megacut3_1_ke85", 85, "6_megaCut3/1_1particlecut"},
        {"6_megaCut3/1_1particlecut_KE90", "cut_megacut3_1_ke90", 90, "6_megaCut3/1_1particlecut"},
        {"6_megaCut3/1_1particlecut_KE95", "cut_megacut3_1_ke95", 95, "6_megaCut3/1_1particlecut"},
        {"5_avgTrackerHitsPerPlane", "cut_avgHitsPlane", 0, "5_avgTrackerHitsPerPlane"}};
    return specs;
}

// the cut as the titles name it: the first line of cuts.txt of the folder after "<dir>: " (Plot3D.cpp writes it), the folder name when there is none
static std::string TubeCutLabel(const TubeSpec& s){
    std::string label = s.dir;
    if (!*s.file) return s.kePct ? Form("reconstructed ZDC protons with ZDC KE > %d%% of the beam KE", s.kePct) : "all reconstructed ZDC protons, no cuts";
    if (FILE* f = fopen(Form("plots/6_cuts/%s/cuts.txt", s.labelDir), "r")){
        char buf[4096];
        if (fgets(buf, sizeof buf, f)){std::string l = buf; while (!l.empty() && (l.back() == '\n' || l.back() == '\r')) l.pop_back(); const size_t p = l.find(": "); label = p == std::string::npos ? l : l.substr(p + 2);}
        fclose(f);
    }
    if (s.kePct) if (const size_t p = label.find("80%"); p != std::string::npos) label.replace(p, 3, std::to_string(s.kePct) + "%");
    return label;
}

static void SavePage(TH2D& h, const std::string& path, const std::string& title){
    TCanvas c("cTube", "Plots", 1500, 1200);
    c.SetRightMargin(0.15); c.SetLeftMargin(0.14); c.SetBottomMargin(0.13);
    h.SetTitle(PlotTitle(path, title).c_str());
    h.SetStats(0);
    h.GetZaxis()->SetTitle("Particles");
    h.GetZaxis()->SetTitleOffset(1.4);
    h.Draw("COLZ");
    c.SaveAs(path.c_str());
}

static std::pair<double, double> Range(const std::vector<double>& v){
    if (v.empty()) return {0, 1};
    const auto mm = std::minmax_element(v.begin(), v.end());
    const double pad = 0.02 * (*mm.second - *mm.first) + 1e-6;
    return {*mm.first - pad, *mm.second + pad};
}

// the six tube end point columns of the rows of the sample (empty: the cut file is missing or has no tube columns); nRows = the rows read. The order of the rows does not matter for the histograms: ColumnsUnordered reads on all threads in any order
static std::vector<std::vector<double>> TubeColumns(const TubeSpec& s, long long& nRows){
    const std::vector<std::string> tube = {"tubeX0", "tubeY0", "tubeZ0", "tubeX1", "tubeY1", "tubeZ1"};
    nRows = 0;
    if (*s.file){
        TTree* t = OpenCut(std::string("cuts/") + s.file + ".root", "AnglesRaw");
        if (!t){printf("TubeCenter %s: cuts/%s.root (AnglesRaw) missing, skipped\n", s.dir, s.file); return {};}
        if (!t->GetBranch("tubeX0")){printf("TubeCenter %s: cuts/%s.root has no tube columns (rerun CutPlots.cpp samples), skipped\n", s.dir, s.file); return {};}
        nRows = t->GetEntries();
        return ColumnsUnordered(t, tube);
    }
    // (fEvent, TrackID) of the reconstructed ZDC protons above the KE cut (CutPlots.cpp keKeys: the highest ZDC KE if a key repeats), then their rows of AnglesRaw
    TTree* z = OpenCut("cuts/cut_ZdcAllEv.root", "ParticleInfo");
    TTree* ang = OpenCut("cuts/cut_angles_all.root", "AnglesRaw");
    if (!z || !ang || !ang->GetBranch("tubeX0")){printf("TubeCenter %s: cut_ZdcAllEv.root or cut_angles_all.root (with tube columns) missing, skipped\n", s.dir); return {};}
    const auto k = ColumnsUnordered(z, {"fEvent", "TrackID", "fDEnergy"}, "isReco == 1");
    std::map<std::pair<int, int>, double> keOf;
    for (size_t i = 0; i < k[0].size(); i++){double& v = keOf.emplace(std::make_pair((int)k[0][i], (int)k[1][i]), k[2][i]).first->second; if (k[2][i] > v) v = k[2][i];}
    const double cutVal = kBeamKE * s.kePct / 100.0;
    std::vector<std::string> cols = {"fEvent", "TrackID"}; cols.insert(cols.end(), tube.begin(), tube.end());
    const auto all = ColumnsUnordered(ang, cols);
    std::vector<std::vector<double>> a(6);
    for (size_t i = 0; i < all[0].size(); i++){
        const auto it = keOf.find({(int)all[0][i], (int)all[1][i]});
        if (it == keOf.end() || !(s.kePct == 0 || it->second > cutVal)) continue;
        for (int c = 0; c < 6; c++) a[c].push_back(all[2 + c][i]);
    }
    nRows = (long long)a[0].size();
    return a;
}

static void TubePages(const TubeSpec& s){
    long long nRows = 0;
    const auto a = TubeColumns(s, nRows);
    if (a.empty()) return;
    const std::string out = std::string("plots/6_cuts/") + s.dir + "/1_additionalinfo/";
    gSystem->mkdir(out.c_str(), true);
    std::vector<double> x, y, z, r, th, ph;
    x.reserve(a[0].size()); y.reserve(a[0].size()); z.reserve(a[0].size()); r.reserve(a[0].size()); th.reserve(a[0].size()); ph.reserve(a[0].size());
    for (size_t i = 0; i < a[0].size(); i++){
        const double cx = 0.5 * (a[0][i] + a[3][i]), cy = 0.5 * (a[1][i] + a[4][i]), cz = 0.5 * (a[2][i] + a[5][i]);
        if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(cz)) continue;
        const double m = std::sqrt(cx * cx + cy * cy + cz * cz);
        x.push_back(cx); y.push_back(cy); z.push_back(cz); r.push_back(std::hypot(cx, cy));
        th.push_back(std::acos(std::max(-1.0, std::min(1.0, cz / m))) * TMath::RadToDeg());
        double p = std::atan2(cy, cx) * TMath::RadToDeg(); if (p < 0) p += 360;
        ph.push_back(p);
    }
    const size_t n = x.size();
    printf("TubeCenter %s: %zu of %lld rows have a tube\n", s.dir, n, nRows);
    if (!n) return;
    const std::string what = "all recon particles of " + std::string(s.dir) + " (" + TubeCutLabel(s) + "), tube centre, " + Form("%zu rows", n);
    const auto rx = Range(x), ry = Range(y), rr = Range(r), rz = Range(z), rt = Range(th);

    { TH2D h("hTubeXY", "", 100, rx.first, rx.second, 100, ry.first, ry.second);
      for (size_t i = 0; i < n; i++) h.Fill(x[i], y[i]);
      SavePage(h, out + "15_tubeCenter_xy.pdf", what + ";Tube centre X (mm);Tube centre Y (mm)"); }
    { TH2D h("hTubeRTh", "", 100, rt.first, rt.second, 100, rr.first, rr.second);
      for (size_t i = 0; i < n; i++) h.Fill(th[i], r[i]);
      SavePage(h, out + "16_tubeCenter_r_vs_estNeuTheta.pdf", what + ";Estimated neutron theta (degrees);r = #sqrt{X^{2} + Y^{2}} of the tube centre (mm)"); }
    { TH2D h("hTubeRPh", "", 72, 0, 360, 100, rr.first, rr.second);
      for (size_t i = 0; i < n; i++) h.Fill(ph[i], r[i]);
      SavePage(h, out + "17_tubeCenter_r_vs_estNeuPhi.pdf", what + ";Estimated neutron phi (degrees);r = #sqrt{X^{2} + Y^{2}} of the tube centre (mm)"); }
    { TH2D h("hTubeRZ", "", 100, rz.first, rz.second, 100, rr.first, rr.second);
      for (size_t i = 0; i < n; i++) h.Fill(z[i], r[i]);
      SavePage(h, out + "18_tubeCenter_r_vs_z.pdf", what + ";Tube centre Z (mm);r = #sqrt{X^{2} + Y^{2}} of the tube centre (mm)"); }
}

void TubeCenter(){
    ROOT::EnableImplicitMT(NThreads());
    FastJit();
    gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255);
    const char* only = std::getenv("TUBE_CUT");
    for (const TubeSpec& s : TubeSpecs()) if (!only || !*only || std::string(only) == s.dir) TubePages(s);
}
