#pragma once
#include <memory>
#include <TAxis.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>

#include <TLegend.h>
#include <TList.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <set>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>
#include <thread>
#include <functional>
#include <numeric>
#include <cstdlib>
#include <sys/stat.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TTree.h>

// Every plotting macro includes this header: plain (not bold) Helvetica for every title, label, stats box, legend and text box.
// A function, not a global with a lambda initializer: ACLiC cannot load a macro with a runtime-initialized global on LCRC (Failed to materialize symbols).
// NThreads() calls it, and every macro calls NThreads() at its start.
inline void PlainFonts(){
    gStyle->SetTextFont(42); gStyle->SetStatFont(42); gStyle->SetLegendFont(42);
    gStyle->SetTitleFont(42, "XYZ"); gStyle->SetLabelFont(42, "XYZ"); gStyle->SetTitleFont(42, ""); // "": the pad title
}

// Beam of the run: the mean momentum and kinetic energy the generator (5generator.cc) wrote to settings.csv (BeamMomentum, BeamKE, MeV), set by /randGun/Energy
// and /randGun/EnergySD. Without settings.csv or the keys (an older run) the old beam of 25 GeV/c, which is a KE of sqrt(p^2 + m^2) - m = 24078 MeV, not 25 GeV.
inline std::map<std::string, double> ReadSettings(const std::string& path); // defined below, its default path is added there
static const double kNeutronMass = 939.565;   // MeV
inline double BeamSetting(const char* name, double fallback){
    try {
        const auto s = ReadSettings("settings.csv");
        const auto it = s.find(name);
        if (it != s.end() && it->second > 0) return it->second;
    } catch (const std::exception&) {}
    if (std::string(name) == "BeamMomentum") std::fprintf(stderr, "Warning: no BeamMomentum in settings.csv, using the 25 GeV/c beam of the first runs\n");
    return fallback;
}
static const double kBeamMomentum = BeamSetting("BeamMomentum", 25000.0); // MeV/c
static const double kBeamKE = BeamSetting("BeamKE", std::sqrt(25000.0 * 25000.0 + kNeutronMass * kNeutronMass) - kNeutronMass); // MeV
static const double kBeamMomentumSD = BeamSetting("BeamMomentumSD", 0.0); // MeV/c, the spread of the beam momentum (/randGun/EnergySD), 0 without one

// Axis limits that follow the beam; the fixed values of the first runs were made for the 25 GeV/c beam. KeAxisMax: 1.2 x the beam KE rounded up to 5 GeV (30000 MeV at 25 GeV/c).
// Theta at a given t goes as 1 / p, so the theta axes and theta windows scale by ThetaScale (2.3 deg at 25 GeV/c). t itself does not scale with p: the t axes stay fixed.
inline double KeAxisMax(){return std::ceil(1.2 * kBeamKE / 5000.0 - 1e-9) * 5000.0;} // MeV
inline double ThetaScale(){return 25000.0 / kBeamMomentum;}
inline double ThetaAxisMax(){return 2.3 * ThetaScale();} // deg

// The beam spread as relative errors, the same for every event. t = (p sin(theta))^2 gives sigma_t / t = 2 sigma_p / p. A_P is -kAPA (1 - exp(-kAPk p)) p sin(theta) (APFromT), so
// sigma_A / |A| = sigma_p (1 / p + kAPk e^(-kAPk p) / (1 - e^(-kAPk p))), the theta cancels. P is the fit amplitude over A, so P carries the same fraction as A (PFit, AddAErrToP).
// pBeam, sdP in MeV/c; kAPk is per GeV in the model.
static const double kAPA = 0.7079, kAPk = 0.3860;
inline double BeamTFrac(double pBeam, double sdP){return pBeam > 0 ? 2 * sdP / pBeam : 0;}
inline double BeamAFrac(double pBeam, double sdP){
    if (pBeam <= 0) return 0;
    const double k = kAPk / 1000.0, e = std::exp(-k * pBeam); // per MeV
    return sdP * (1 / pBeam + k * e / (1 - e));
}
static const double kBeamTFrac = BeamTFrac(kBeamMomentum, kBeamMomentumSD);
static const double kBeamAFrac = BeamAFrac(kBeamMomentum, kBeamMomentumSD);

// Transverse-only t in GeV^2: (p_n sin(theta))^2, theta = proton lab polar angle (rad), p_n = lab beam momentum
// (MeV). The full invariant also has an energy term, which the measured final KE (after energy loss) overwhelms.
inline double TFromAngle(double pn, double theta){
    const double pT = pn * std::sin(theta) / 1000.0;
    return pT * pT;
}

// Same formula as the model (G4ChargeExchange.cc, per GeV there), including its minus sign. pbeam in MeV, t in GeV^2.
// AP times apScale (the sim's /CE/APScale), clamped to -1.0 at the low end like the model.
inline double APFromT(double pbeam, double t, double apScale){
    return std::max(-1.0, apScale * -kAPA * (1 - std::exp(-kAPk * pbeam / 1000.0)) * std::sqrt(std::max(t, 0.0)));
}

// Run settings the sim wrote to settings.csv (a row of names, a row of values), as name -> value (.at("APScale"), ...).
// Throws if the file is missing or the two rows do not match.
inline std::map<std::string, double> ReadSettings(const std::string& path = "settings.csv"){
    // memoised by path and mtime+size (called per plot): returns a copy, so callers see no difference
    static std::map<std::string, std::pair<std::pair<long long, long long>, std::map<std::string, double>>> memo;
    struct stat st;
    const bool haveStat = (::stat(path.c_str(), &st) == 0);
    const std::pair<long long, long long> stamp = haveStat ? std::make_pair((long long)st.st_mtime, (long long)st.st_size) : std::make_pair(-1LL, -1LL);
    if (haveStat){
        auto it = memo.find(path);
        if (it != memo.end() && it->second.first == stamp) return it->second.second;
    }
    std::ifstream in(path);
    std::string names, values;
    if (!std::getline(in, names) || !std::getline(in, values))
        throw std::runtime_error("cannot read " + path + " (copy it from the sim's root/<folder>/ next to the ROOT files)");
    std::map<std::string, double> settings;
    std::istringstream nameStream(names), valueStream(values);
    std::string name, value;
    while (std::getline(nameStream, name, ',') && std::getline(valueStream, value, ',')){
        if (!name.empty() && name[0] == '#') name.erase(0, 1); // # marks a default value, the name is the same
        settings[name] = std::stod(value);
    }
    if (settings.empty()) throw std::runtime_error(path + " has no settings");
    if (haveStat) memo[path] = {stamp, settings};
    return settings;
}

// The tracker hit cuts (5_nHits, 6_megaCut, 7_megaCut2, 9_megaCut3) keep the even counts 2, 4, ... 2 * nBlocks (1 block: 2; 3 blocks: 2, 4, 6; 7 blocks: 2 to 14), one pair of hits per block.
// MaxTrackerHits reads settings.csv, so call it once outside an event loop and pass the value to EvenHits.
inline int MaxTrackerHits(){return 2 * (int)std::lround(ReadSettings().at("nBlocks"));}
inline bool EvenHits(int n, int maxHits){return n >= 2 && n <= maxHits && n % 2 == 0;}
// the Draw string of the same cut on expression e, e.g. (e == 2 || e == 4 || e == 6)
inline std::string EvenHitsDraw(const std::string& e){
    std::string s = "(";
    for (int n = 2, maxHits = MaxTrackerHits(); n <= maxHits; n += 2) s += (n > 2 ? " || " : "") + e + " == " + std::to_string(n);
    return s + ")";
}
// the ratio and FOM graphs: round markers (size 0.6) joined by a width 1 line in one color, drawn "APL" with x from 0 to xMax
inline void DrawPointsLine(TGraph& g, int color, double xMax){
    g.SetMarkerStyle(20); g.SetMarkerSize(0.6); g.SetLineWidth(1); g.SetLineColor(color); g.SetMarkerColor(color);
    g.Draw("APL");
    g.GetXaxis()->SetLimits(0, xMax);
}
// the counts for a title, "2, 4, 6"
inline std::string EvenHitsText(){
    std::string s;
    for (int n = 2, maxHits = MaxTrackerHits(); n <= maxHits; n += 2) s += (n > 2 ? ", " : "") + std::to_string(n);
    return s;
}

// Threads for ROOT::EnableImplicitMT: env ANA_NTHREADS, else the hardware thread count (plot.sh splits the cores between its lanes).
inline unsigned NThreads(){
    static const bool fonts = (PlainFonts(), true); (void)fonts; // first call sets the fonts of every plot
    if (const char* e = std::getenv("ANA_NTHREADS")){
        const long v = std::atol(e);
        if (v > 0) return (unsigned)v;
    }
    const unsigned h = std::thread::hardware_concurrency();
    return h ? h : 1;
}

// Runs f(lo, hi) on contiguous row ranges of [0, n) on NThreads() threads (8 ranges per thread) and returns the results in range order. The row loops that follow a Columns or
// DrawPlan read use it, so they run on all the cores instead of one. f only reads shared data; each range builds and returns its own partial result.
#include <ROOT/TThreadExecutor.hxx>
template <class F>
inline auto MapRows(Long64_t n, F f) -> std::vector<decltype(f(Long64_t(0), Long64_t(0)))>
{
    const unsigned nChunk = 8 * NThreads();
    std::vector<unsigned> ids(nChunk);
    std::iota(ids.begin(), ids.end(), 0u);
    ROOT::TThreadExecutor ex(NThreads());
    return ex.Map([&](unsigned c){return f(n * c / nChunk, n * (c + 1) / nChunk);}, ids);
}

// Rows collected one entry at a time, per RDataFrame slot, as one vector per column (RowsAppend), and the slots appended column by column (RowsMerge, each slot column freed at once).
// The values of one row go into the same position of every column, so the columns stay aligned whatever the order of the rows. Used by ColumnsUnordered (TreeColumns.h) and DrawPlan.
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
using RowCols = std::vector<std::vector<double>>;
inline void RowsAppend(RowCols& acc, const ROOT::RVecD& v){
    if (acc.size() < v.size()) acc.resize(v.size());
    for (size_t c = 0; c < v.size(); c++) acc[c].push_back(v[c]);
}
inline void RowsMerge(std::vector<RowCols>& slots){
    RowCols& dst = slots[0];
    size_t ncol = 0; for (auto& s : slots) ncol = std::max(ncol, s.size());
    dst.resize(ncol);
    for (size_t c = 0; c < ncol; c++){
        size_t tot = 0; for (auto& s : slots) if (c < s.size()) tot += s[c].size();
        dst[c].reserve(tot);
        for (size_t i = 1; i < slots.size(); i++) if (c < slots[i].size()){
            dst[c].insert(dst[c].end(), slots[i][c].begin(), slots[i][c].end());
            std::vector<double>().swap(slots[i][c]); // clear() would keep the capacity
        }
    }
}

// A tree with more than one cluster (a multithreaded Snapshot writes one cluster per task). rdfentry_ is wrong under implicit MT on such a tree (ROOT 6.38: 48 to 76% of the rows of a 5-cluster
// AnglesRaw at 4 and 8 threads), so the macros take the entry number from WithEntryNumber below. MultiCluster only tells the order-free readers (ColumnsUnordered, DrawPlan::unordered) which trees they are for.
inline bool MultiCluster(TTree* t){
    auto it = t->GetClusterIterator(0);
    const Long64_t n = t->GetEntries();
    return n > 0 && it.Next() < n && it.Next() < n; // a second cluster starts before the last entry
}
// The entry number of each row in its tree, exact under implicit MT, as a ULong64_t column col (default "entryNo"); use it instead of rdfentry_ (tested on a 5-cluster AnglesRaw: rdfentry_ labelled 48 to 76% of
// the rows wrong at 4 and 8 threads, this was exact in 18 of 18 runs). Every task starts at s.EntryRange().first and counts the rows it reads. The counter has to see every row of the tree, so the column is
// forced by an always true Filter right here, before any Filter of the caller drops rows. Call it on a node of an RDataFrame built after ROOT::EnableImplicitMT; nSlots is df.GetNSlots().
// For a TChain the number is the one inside the file of the row, which orders the rows of one file (and so of one event) correctly.
inline ROOT::RDF::RNode WithEntryNumber(ROOT::RDF::RNode df, unsigned nSlots, const std::string& col = "entryNo"){
    auto count = std::make_shared<std::vector<Long64_t>>(std::max(1u, nSlots) + 1, 0);
    return df.DefinePerSample(col + "_b", [count](unsigned slot, const ROOT::RDF::RSampleInfo& s){(*count)[slot] = 0; return (Long64_t)s.EntryRange().first;})
             .DefineSlot(col, [count](unsigned slot, Long64_t first){return (ULong64_t)(first + (*count)[slot]++);}, {col + "_b"})
             .Filter([](ULong64_t){return true;}, {col});
}

// z (mm) of the carbon blocks and tracker planes of the run (2construction.cc). stackSpread = 1 (settings.csv): the blocks sit in evenly spaced slots between z = 0 and the ZDC
// front face (7700 mm), the face being the last slot (block i at i * 7700 / (nBlocks + 1)); a tracker plane dz behind every block and dz in front of every block but the first,
// and one plane dz in front of the ZDC (last in planeZ). Here dz is the gap between the edge of the block and the edge of the 0.5 mm tracker (the ZDC face for the last plane).
// A run without the key is the old layout, stack centred at 4 m, elements dz apart centre to centre: each block is followed by nBlocksTracking
// tracker planes; frontTrackers (settings.csv, 1 when the run has none, an older sim) puts one more group of planes before the first block.
struct StackLayout {std::vector<double> blockZ, planeZ;};
inline StackLayout StackZ(const std::map<std::string, double>& s){
    const int nBlocks = (int)std::lround(s.at("nBlocks")), nTrk = (int)std::lround(s.at("nBlocksTracking"));
    const double dzMm = s.at("dz") * 1000.0;
    if (s.count("stackSpread") && std::lround(s.at("stackSpread")) == 1){
        const double zdcFrontMm = 7700.0, pitch = zdcFrontMm / (nBlocks + 1);
        const double offset = s.at("CarbonThickness") * 1000.0 / 2.0 + dzMm + 0.25; // block centre to tracker centre: half block + gap + half tracker (0.5 mm)
        StackLayout sp;
        for (int i = 0; i < nBlocks; i++){
            const double zc = (i + 1) * pitch;
            sp.blockZ.push_back(zc);
            if (i > 0) sp.planeZ.push_back(zc - offset);
            sp.planeZ.push_back(zc + offset);
        }
        sp.planeZ.push_back(zdcFrontMm - dzMm - 0.25);
        return sp;
    }
    const int front = s.count("frontTrackers") ? (int)std::lround(s.at("frontTrackers")) : 1;
    const int nElements = (nBlocks + front) * nTrk + nBlocks;
    const double zStart = 4000.0 - (nElements - 1) * dzMm / 2.0;
    StackLayout st; int element = 0;
    auto planes = [&](){for (int j = 0; j < nTrk; j++, element++) st.planeZ.push_back(zStart + element * dzMm);};
    if (front) planes();
    for (int i = 0; i < nBlocks; i++){st.blockZ.push_back(zStart + element * dzMm); element++; planes();}
    return st;
}

// The charge exchange selection of every stage: a CE with its vertex inside one of the carbon blocks (CE in the tracker silicon is left out).
// The block z positions mirror MyDetectorConstruction (2construction.cc): the blocks and planes as StackZ gives them. fZ < 7700 is the front face of the ZDC.
inline std::string CarbonCeCut(){
    static std::string memo; static std::map<std::string, double> memoSettings;
    const auto s = ReadSettings();
    if (!memo.empty() && s == memoSettings) return memo;
    memoSettings = s;
    const double halfMm = s.at("CarbonThickness") * 1000.0 / 2.0;
    std::string windows;
    for (const double zc : StackZ(s).blockZ){
        char w[96]; std::snprintf(w, sizeof w, "%s(fZ>%.3f && fZ<%.3f)", windows.empty() ? "" : " || ", zc - halfMm, zc + halfMm);
        windows += w;
    }
    memo = "didChargeExchange>0 && fZ<7700 && (" + windows + ")";
    return memo;
}

// Plot title in the shared "N. Quantity, Label" format. title is "title;xaxis;yaxis". N is the leading digit groups of the
// file name joined by dots (1_1_xy_all.pdf -> 1.1), omitted if the file has none. Spaces around every ';' segment are trimmed.
inline std::string PlotTitle(const std::string& path, const std::string& title){
    std::string file = path.substr(path.find_last_of('/') + 1);
    std::string num;
    for (size_t i = 0; i < file.size() && std::isdigit((unsigned char)file[i]); ){
        size_t j = i;
        while (j < file.size() && std::isdigit((unsigned char)file[j])) j++;
        num += (num.empty() ? "" : ".") + file.substr(i, j - i);
        i = (j + 1 < file.size() && file[j] == '_') ? j + 1 : file.size();
    }
    std::string out = num.empty() ? "" : num + ". ";
    size_t start = 0;
    for (bool first = true; start <= title.size(); first = false){
        size_t end = title.find(';', start);
        if (end == std::string::npos) end = title.size();
        size_t a = title.find_first_not_of(' ', start), b = title.find_last_not_of(' ', end ? end - 1 : 0);
        std::string seg = (a == std::string::npos || a >= end) ? "" : title.substr(a, b - a + 1);
        out += (first ? "" : ";") + seg;
        start = end + 1;
    }
    return out;
}

// Seed from median and IQR (robust to sparse bins and far-away outliers),
// then refit on a window around the fitted mean until it settles.
// A fit drawn over a sub-range stops at the range edge. Adds dotted continuations of f from each edge out to the
// histogram axis limits, so it is visible where the fit curve would keep going. Appended after the fit function, so the
// stats box still reports the real fit.
inline void DrawFitTails(TH1* h, TF1* f){
    double lo = f->GetXmin(), hi = f->GetXmax();
    auto tail = [&](double a, double b, const char* tag){
        if (b <= a) return;
        TF1* t = (TF1*)f->Clone(Form("%s_%s", f->GetName(), tag));
        t->SetRange(a, b);
        t->SetLineStyle(kDotted);
        t->SetLineColor(f->GetLineColor());
        h->GetListOfFunctions()->Add(t);
    };
    tail(h->GetXaxis()->GetXmin(), lo, "tailLo");
    tail(hi, h->GetXaxis()->GetXmax(), "tailHi");
}

// The fit equation in a box in the right margin under the stats box (the fit plots use a 0.30 right margin), added to h's
// function list so it is drawn with h. Parameter names match the stats box.
static constexpr const char* kGausEq ="Constant #times exp(-(x - Mean)^{2} / 2 Sigma^{2})";
// On every plot of all recon protons: the CE protons among them under the same cuts, bottom right below the stats and fit boxes.
// nCEAll: the CE recon protons before any cut, for the share of them the plot's cuts lose (left out when not given), in its own box (DrawCELost).
// nRecon < 0: a CE protons only plot, only the CE lost box.
// nReconAll: every recon proton before any cut, for the share of them in the plot's cut region (left out when not given or nRecon < 0).
// Call it after the legend (and any stats-like info box) is drawn on gPad: the CE lost box is placed under it.
inline void DrawCELost(double nCE, double nCEAll);
inline void DrawCEPercent(double nCE, double nRecon, double nCEAll = -1, double nReconAll = -1){
    const bool keptLine = nRecon >= 0 && nReconAll > 0;
    const int lines = (nRecon >= 0) + keptLine;
    if (nCEAll > 0) DrawCELost(nCE, nCEAll);
    if (lines == 0) return;
    TPaveText* pt = new TPaveText(0.745, 0.02, 0.99, 0.02 + 0.04 * lines, "NDC"); // right of the graph and its tick labels on canvases with a 0.28 or 0.30 right margin
    pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetTextAlign(12); pt->SetTextSize(0.016);
    if (nRecon >= 0) pt->AddText(Form("CE / recon = %.0f / %.0f = %.2f%%", nCE, nRecon, nRecon > 0 ? 100.0 * nCE / nRecon : 0.0));
    if (keptLine) pt->AddText(Form("recon in cut / all = %.0f / %.0f = %.2f%%", nRecon, nReconAll, 100.0 * nRecon / nReconAll));
    pt->Draw();
}

// "CE lost by cuts" in its own one line box right under the legend on gPad, same width and text size; without a legend where the legend would sit:
// under the stats box (the histogram's, or a stats-like text box in the right margin), else under where gStyle puts the stats box.
inline void DrawCELost(double nCE, double nCEAll){
    TLegend* leg = nullptr; TH1* h = nullptr; TPave* info = nullptr;
    for (TObject* o : *gPad->GetListOfPrimitives()){
        if (auto* l = dynamic_cast<TLegend*>(o)) leg = l; // the last legend drawn
        else if (auto* hh = dynamic_cast<TH1*>(o)){if (!h || h->TestBit(TH1::kNoStats)) h = hh;} // the first histogram with a stats box (a DrawFrame frame has none)
        else if (auto* p = dynamic_cast<TPave*>(o)) if (p->GetX1NDC() >= 0.7 && p->GetY1NDC() >= 0.4 && strcmp(p->GetName(), "fitEquation")) info = p; // a fit info box in the right margin (PlotPolarization)
    }
    double x1 = 0.745, x2 = 0.99, top = gStyle->GetStatY() - gStyle->GetStatH(), size = 0.016;
    if (leg){x1 = leg->GetX1NDC(); x2 = leg->GetX2NDC(); top = leg->GetY1NDC(); if (leg->GetTextSize() > 0) size = leg->GetTextSize();}
    else {
        TPave* st = h && !h->TestBit(TH1::kNoStats) ? (TPave*)h->FindObject("stats") : nullptr;
        if (h && !h->TestBit(TH1::kNoStats) && !st){gPad->Modified(); gPad->Update(); st = (TPave*)h->FindObject("stats");} // ROOT makes the stats box when it paints
        if (st) top = st->GetY1NDC();
        else if (info) top = info->GetY1NDC();
    }
    TPaveText* pt = new TPaveText(x1, top - 0.05, x2, top - 0.01, "NDC"); // 0.01 gap: cut_side_box.py tells the frames apart
    pt->SetName("ceLost");
    pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetTextAlign(12); pt->SetTextSize(size);
    pt->AddText(Form("CE lost by cuts = %.0f / %.0f = %.1f%%", nCEAll - nCE, nCEAll, 100.0 * (nCEAll - nCE) / nCEAll));
    pt->Draw();
}

// The fit's stats box goes in the right margin too: ROOT widens a stats box with fit values into the graph, but keeps the place of one that
// already exists, so it is made here at the margin's width. The equation box sits above the CE box of DrawCEPercent (at most 0.02 to 0.10).
inline void AddFitEquation(TH1* h, const std::vector<std::string>& lines){
    if (!h->GetListOfFunctions()->FindObject("stats")){
        TPaveStats* st = new TPaveStats(0.745, 0.58, 0.99, 0.935, "brNDC");
        st->SetName("stats"); st->SetParent(h);
        st->SetOptStat(gStyle->GetOptStat()); st->SetOptFit(gStyle->GetOptFit());
        st->SetBorderSize(1); st->SetFillColor(0); st->SetTextFont(42);
        h->GetListOfFunctions()->Add(st);
    }
    TPaveText* eq = new TPaveText(0.745, 0.17, 0.99, 0.17 + 0.04 * (lines.size() + 1), "NDC");
    eq->SetName("fitEquation");
    eq->SetBorderSize(1); eq->SetFillColor(0); eq->SetTextAlign(12); eq->SetTextSize(0.017);
    eq->AddText("Fit function:");
    for (const auto& l : lines) eq->AddText(l.c_str());
    h->GetListOfFunctions()->Add(eq);
}

// Unbinned fit of the polarization P to the phi of the protons (degrees, from loDeg to hiDeg): every proton enters the likelihood with its own phi, no bins. The density of a proton is
// (1 + P A cos phi) / N(P), N(P) = (hi - lo) + P A (sin hi - sin lo) (angles in rad), A the model A of the sample (fixed, aFix) or of the proton itself (aEach). P is found with Newton's method,
// limited to -1 to 1 (a polarization), and its error is 1 / sqrt(-d2lnL/dP2) at that P. The fit that replaced the chi2 fit of a 50 bin histogram, which the 45 degree grid of the recon block phi
// (spikes) pushed to the limit of P.
// err is the error on P that every caller quotes: the likelihood error (errFit) and the error from A (errA = |P| * the relative error of A, see AMean) in quadrature.
// A fixed to the mean of a sample has the error of that mean (the spread of A over the theta distribution of the sample, std(A) / sqrt(N)) and the beam spread (kBeamAFrac, see BeamAFrac); the A of
// each proton (aEach) has the beam spread only.
struct PFit {double p = 0, err = 0, errFit = 0, errA = 0; bool ok = false;};
inline void AddAErrToP(PFit& r, double aRelErr){
    r.errFit = r.err;
    r.errA = std::fabs(r.p) * aRelErr;
    r.err = std::sqrt(r.errFit * r.errFit + r.errA * r.errA);
}

// The mean A of a sample with its error: A = sum / n over the protons, StatErr = std(A) / sqrt(n) (the error of the mean over the theta distribution of the sample), and Err adds the beam
// spread (kBeamAFrac) in quadrature, relative to the mean: this is the "+-" of the fixed A of a fit, and RelErr is what the P of that fit carries (AddAErrToP).
// Add per proton, or add the sums of several parts (+=): n, sum and sum^2 add over slices.
struct AMean {
    double n = 0, sum = 0, sum2 = 0;
    void Add(double a){n++; sum += a; sum2 += a * a;}
    AMean& operator+=(const AMean& o){n += o.n; sum += o.sum; sum2 += o.sum2; return *this;}
    double Mean() const {return n > 0 ? sum / n : 0.0;}
    double StatErr() const {
        if (n < 2) return 0.0;
        const double m = Mean(), var = std::max(0.0, (sum2 / n - m * m) * n / (n - 1));
        return std::sqrt(var / n);
    }
    double RelErr() const {
        const double m = std::fabs(Mean()), s = m > 0 ? StatErr() / m : 0.0;
        return std::sqrt(s * s + kBeamAFrac * kBeamAFrac);
    }
    double Err() const {return std::fabs(Mean()) * RelErr();}
};
inline AMean AMeanOf(const std::vector<double>& v){AMean a; for (double x : v) a.Add(x); return a;}
// P of an unbinned fit by Newton's method from 0, limited to -1 to 1: derivs(p, d1, d2) gives dlnL/dP and d2lnL/dP2. Stops when P moves less than tol, or when d2 is not negative.
template <class D>
inline double NewtonP(D derivs, double tol){
    double p = 0, d1, d2;
    for (int it = 0; it < 100; it++){
        derivs(p, d1, d2);
        if (!(d2 < 0)) break;
        const double next = std::max(-1.0, std::min(1.0, p - d1 / d2));
        const bool done = std::fabs(next - p) < tol;
        p = next;
        if (done) break;
    }
    return p;
}
// cosPhi: cos(phi) of every proton (a: its A, or null for aFix). parallel: the sums over the protons run on all threads (not when the caller is already one of the threads).
inline PFit UnbinnedPFromCos(const double* cosPhi, Long64_t n, const double* a, const AMean& aMean, double loDeg, double hiDeg, bool parallel = true){
    PFit r;
    const double aFix = aMean.Mean();
    if (n < 10 || (!a && std::fabs(aFix) < 1e-12)) return r;
    const double lo = loDeg * TMath::Pi() / 180, hi = hiDeg * TMath::Pi() / 180, span = hi - lo, s = std::sin(hi) - std::sin(lo);
    struct Sums {double d1 = 0, d2 = 0;};
    const auto sums = [&](double p){
        const auto range = [&](Long64_t from, Long64_t to){
            Sums t;
            for (Long64_t i = from; i < to; i++){
                const double ai = a ? a[i] : aFix, ci = ai * cosPhi[i], q = 1 + p * ci, as = ai * s, nr = span + p * as;
                t.d1 += ci / q - as / nr;
                t.d2 += -ci * ci / (q * q) + as * as / (nr * nr);
            }
            return t;
        };
        Sums tot;
        if (parallel && n > 200000) for (const Sums& t : MapRows(n, range)){tot.d1 += t.d1; tot.d2 += t.d2;}
        else tot = range(0, n);
        return tot;
    };
    const double p = NewtonP([&](double pp, double& d1, double& d2){const Sums t = sums(pp); d1 = t.d1; d2 = t.d2;}, 1e-10);
    const Sums t = sums(p);
    r.p = p; r.err = t.d2 < 0 ? 1 / std::sqrt(-t.d2) : 0; r.ok = true;
    AddAErrToP(r, a ? kBeamAFrac : aMean.RelErr()); // the A of each proton: the beam spread only
    return r;
}
// aMean: the A fixed for every proton (its error goes into P), or, with aEach, the A of each proton (aMean is then not used: pass AMean())
inline PFit UnbinnedP(const std::vector<double>& phiDeg, const AMean& aMean, double loDeg = 0, double hiDeg = 360, const std::vector<double>* aEach = nullptr, bool parallel = true){
    std::vector<double> cs(phiDeg.size());
    for (size_t i = 0; i < cs.size(); i++) cs[i] = std::cos(phiDeg[i] * TMath::Pi() / 180);
    return UnbinnedPFromCos(cs.data(), (Long64_t)cs.size(), aEach ? aEach->data() : nullptr, aMean, loDeg, hiDeg, parallel);
}

// The curve N0 (1 + P A cos phi) of the unbinned P in the histogram h (50 bins over 0 to 360, N0 = entries / bins), where the stats box and WritePlotNumbers (fit_values.csv) look for a fit.
// Without chi2: chi2 and ndf stay 0. draw: the curve is drawn with the histogram; otherwise it only sits in the function list, like a fit with option "0".
// A is the fixed A with its error (the stats box shows A +- aErr instead of +- 0).
inline TF1* AddUnbinnedFit(TH1* h, const char* name, const PFit& r, double aFix, double aErr, bool draw){
    TF1* f = new TF1(name, "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
    f->SetParNames("N0", "P", "A");
    f->SetParameters(h->GetEntries() / h->GetNbinsX(), r.p, aFix);
    f->SetParError(1, r.err);
    f->SetParError(2, aErr);
    f->SetNDF(0); f->SetChisquare(0);
    if (!draw) f->SetBit(TF1::kNotDraw);
    h->GetListOfFunctions()->Add(f);
    return f;
}

// The same unbinned fit from sums, for the maps that fit thousands of theta ranges (PlotFomPolarization.cpp): ln(1 + x) = x - x^2/2 + ... with x = P A cos(phi) never exceeds |0.12| (|A| < 0.2),
// so the log likelihood of any set of protons is the polynomial sum_k (-1)^(k+1) (P A)^k m_k / k - n ln N(P), m_k = sum over the protons of cos^k(phi), and the sums of cos^k add over slices.
// Eight terms: the left out term is below 1e-8 of n. Same density, range and limit of P as UnbinnedP.
inline constexpr int kPhiMoments = 8;
inline void AddCosMoments(double* m, double phiDeg){
    const double c = std::cos(phiDeg * TMath::Pi() / 180);
    double pw = 1;
    for (int k = 0; k < kPhiMoments; k++){pw *= c; m[k] += pw;}
}
inline PFit UnbinnedPFromMoments(const double* m, double n, const AMean& aMean, double loDeg = 0, double hiDeg = 360){
    PFit r;
    const double aFix = aMean.Mean();
    if (n < 10 || std::fabs(aFix) < 1e-12) return r;
    const double lo = loDeg * TMath::Pi() / 180, hi = hiDeg * TMath::Pi() / 180, span = hi - lo, as = aFix * (std::sin(hi) - std::sin(lo));
    const auto derivs = [&](double p, double& d1, double& d2){
        d1 = 0; d2 = 0;
        double ak = 1; // A^k
        for (int k = 1; k <= kPhiMoments; k++){
            ak *= aFix;
            const double sign = k % 2 ? 1 : -1;
            d1 += sign * ak * std::pow(p, k - 1) * m[k - 1];
            if (k > 1) d2 += sign * (k - 1) * ak * std::pow(p, k - 2) * m[k - 1];
        }
        const double nr = span + p * as;
        d1 -= n * as / nr; d2 += n * as * as / (nr * nr);
    };
    const double p = NewtonP(derivs, 1e-12);
    double d1, d2;
    derivs(p, d1, d2);
    r.p = p; r.err = d2 < 0 ? 1 / std::sqrt(-d2) : 0; r.ok = true;
    AddAErrToP(r, aMean.RelErr());
    return r;
}

// The RDataFrame string expressions (Filter, Define, Histo of a Draw string) are compiled by cling when the macro runs, with the optimizer on, and that was most of the time of the plot macros:
// they read small trees, so the code of an expression runs for a few hundred thousand rows only. ".O 0" halved the time of Plot3D (ANA_JIT_O1=1 keeps the default). Not for the macros that read the big trees.
inline void FastJit(){if (!std::getenv("ANA_JIT_O1")) gInterpreter->ProcessLine(".O 0");}

// Bin string "n,lo,hi" for a plot of the ZDC block positions (fDx, fDy): the ZDC is a grid of nCols cells over 400 mm, and every bin holds the same whole number of cells.
// A bin width that is not a multiple of the pitch (a plain 100 bins over 100 mm on a pitch of 400 / 801 mm) holds one cell more in some bins than in the others:
// a bin that starts on a cell center, like the one at x = 0, holds the extra cell, which is a line at x = 0 and y = 0 in the 2D plots and spikes in the slices.
// The bin is k cells wide, the closest width to targetMm, and the bins are centered on the grid: an odd nCols has a cell at 0, so n and k are odd, an even one has an edge at 0, so n is even.
inline std::string BlockBins(int nCols, double targetMm, double halfRangeMm){
    const double pitch = 400.0 / nCols;
    int k = std::max(1, (int)std::lround(targetMm / pitch));
    if (nCols % 2 == 1 && k % 2 == 0) k = (k > 1 && std::fabs((k - 1) * pitch - targetMm) <= std::fabs((k + 1) * pitch - targetMm)) ? k - 1 : k + 1;
    int n = std::max(1, (int)std::lround(2.0 * halfRangeMm / (k * pitch)));
    if (nCols % 2 == 1 && n % 2 == 0) n++;
    if (nCols % 2 == 0 && n % 2 == 1) n++;
    const double lo = 0.5 * n * k * pitch;
    return Form("%d,%.6f,%.6f", n, -lo, lo);
}

// ---- Unbinned fits of a gaussian (the position, energy and residual plots) ----
// The raw values x of the plot enter the likelihood one by one, no bins. In window mode (cell == 0) the density of a value in [lo, hi] is
//   (1 - fb) [ (1 - f2) G(x; mu, s) + f2 G(x; mu2, s2) ] + fb / (hi - lo),
// each G a gaussian normalised over the window (erf), the second gaussian (nGaus 2) and the flat part (flat) optional. In cell mode (cell > 0, the values sit on a grid of that pitch,
// the ZDC blocks) the probability of a value is the gaussian probability of its cell, normalised over [lo, hi].
// Parameters, in this order: mu, s, [mu2, s2, f2 if nGaus is 2], [fb if flat]. The minimiser is a plain Nelder-Mead (no Minuit2: it is not in every ROOT build) and the errors come from the
// inverse of the finite-difference Hessian of -ln L at the minimum, like the P fit's 1 / sqrt(-d2lnL/dP2).
constexpr int kFitNone = 0, kFitCore = 1, kFitBlocks = 2; // the gaussian fit a plot gets (Plot3D.cpp drawPlot): none, FitCore (window, background, second gaussian), FitBlocks (block cells)
inline double GausCdf(double z){return 0.5 * std::erfc(-z / 1.4142135623730951);}

struct UGaussModel {
    double lo = 0, hi = 0, cell = 0; int nGaus = 1; bool flat = false;
    std::vector<double> pLo, pHi; // bounds of the parameters
};
inline double UGaussNll(const std::vector<double>& x, const UGaussModel& m, const std::vector<double>& p){
    for (size_t k = 0; k < p.size(); k++) if (!(p[k] >= m.pLo[k] && p[k] <= m.pHi[k])) return 1e300;
    const double mu = p[0], s = p[1];
    double mu2 = 0, s2 = 1, f2 = 0, fb = 0;
    size_t ix = 2;
    if (m.nGaus == 2){mu2 = p[2]; s2 = p[3]; f2 = p[4]; ix = 5;}
    if (m.flat) fb = p[ix];
    const double z1 = GausCdf((m.hi - mu) / s) - GausCdf((m.lo - mu) / s);
    const double z2 = m.nGaus == 2 ? GausCdf((m.hi - mu2) / s2) - GausCdf((m.lo - mu2) / s2) : 1.0;
    if (!(z1 > 1e-12) || !(z2 > 1e-12)) return 1e300;
    const double inv = 1.0 / std::sqrt(2.0 * TMath::Pi()), flatDens = fb / (m.hi - m.lo);
    double nll = 0;
    for (double v : x){
        double pdf;
        if (m.cell > 0) pdf = (GausCdf((v + 0.5 * m.cell - mu) / s) - GausCdf((v - 0.5 * m.cell - mu) / s)) / z1;
        else {
            const double a = (v - mu) / s;
            pdf = inv * std::exp(-0.5 * a * a) / (s * z1);
            if (m.nGaus == 2){const double b = (v - mu2) / s2; pdf = (1 - f2) * pdf + f2 * inv * std::exp(-0.5 * b * b) / (s2 * z2);}
            pdf = (1 - fb) * pdf + flatDens;
        }
        nll -= std::log(std::max(pdf, 1e-300));
    }
    return nll;
}

inline std::vector<double> NelderMead(const std::function<double(const std::vector<double>&)>& f, const std::vector<double>& x0, const std::vector<double>& step){
    const size_t n = x0.size();
    std::vector<double> best = x0;
    for (int restart = 0; restart < 3; restart++){ // a restart from the result gets out of a stalled simplex
        std::vector<std::vector<double>> s(n + 1, best);
        for (size_t i = 0; i < n; i++) s[i + 1][i] += step[i] * (restart ? 0.3 : 1.0);
        std::vector<double> fv(n + 1);
        for (size_t i = 0; i <= n; i++) fv[i] = f(s[i]);
        for (int it = 0; it < 3000; it++){
            std::vector<size_t> idx(n + 1);
            std::iota(idx.begin(), idx.end(), 0);
            std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){return fv[a] < fv[b];});
            std::vector<std::vector<double>> s2(n + 1); std::vector<double> f2(n + 1);
            for (size_t i = 0; i <= n; i++){s2[i] = s[idx[i]]; f2[i] = fv[idx[i]];}
            s = s2; fv = f2;
            double spread = 0;
            for (size_t i = 1; i <= n; i++) for (size_t k = 0; k < n; k++) spread = std::max(spread, std::fabs(s[i][k] - s[0][k]) / (std::fabs(s[0][k]) + step[k]));
            if (spread < 1e-8 && std::fabs(fv[n] - fv[0]) < 1e-9 * (std::fabs(fv[0]) + 1.0)) break;
            std::vector<double> c(n, 0.0);
            for (size_t i = 0; i < n; i++) for (size_t k = 0; k < n; k++) c[k] += s[i][k] / n;
            auto along = [&](double t){std::vector<double> q(n); for (size_t k = 0; k < n; k++) q[k] = c[k] + t * (s[n][k] - c[k]); return q;}; // t = -1 reflect, -2 expand, 0.5 contract inside, -0.5 contract outside
            const std::vector<double> xr = along(-1.0); const double fr = f(xr);
            if (fr < fv[0]){
                const std::vector<double> xe = along(-2.0); const double fe = f(xe);
                if (fe < fr){s[n] = xe; fv[n] = fe;} else {s[n] = xr; fv[n] = fr;}
            } else if (fr < fv[n - 1]){
                s[n] = xr; fv[n] = fr;
            } else {
                const std::vector<double> xc = along(fr < fv[n] ? -0.5 : 0.5); const double fc = f(xc);
                if (fc < std::min(fr, fv[n])){s[n] = xc; fv[n] = fc;}
                else for (size_t i = 1; i <= n; i++){for (size_t k = 0; k < n; k++) s[i][k] = s[0][k] + 0.5 * (s[i][k] - s[0][k]); fv[i] = f(s[i]);}
            }
        }
        size_t b = 0;
        for (size_t i = 1; i <= n; i++) if (fv[i] < fv[b]) b = i;
        best = s[b];
    }
    return best;
}

// Errors of the parameters p at the minimum of f (-ln L): sqrt of the diagonal of the inverse Hessian, from central differences with steps h. 0 where the Hessian cannot be inverted
// (a parameter at its bound).
inline std::vector<double> NllErrors(const std::function<double(const std::vector<double>&)>& f, const std::vector<double>& p, const std::vector<double>& h){
    const size_t n = p.size();
    std::vector<double> err(n, 0.0);
    std::vector<std::vector<double>> H(n, std::vector<double>(2 * n, 0.0)); // [Hessian | identity], inverted in place
    const double f0 = f(p);
    auto at = [&](size_t i, double di, size_t j, double dj){std::vector<double> q = p; q[i] += di; q[j] += dj; return f(q);};
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j <= i; j++){
            double v;
            if (i == j) v = (at(i, h[i], i, 0) - 2 * f0 + at(i, -h[i], i, 0)) / (h[i] * h[i]);
            else v = (at(i, h[i], j, h[j]) - at(i, h[i], j, -h[j]) - at(i, -h[i], j, h[j]) + at(i, -h[i], j, -h[j])) / (4 * h[i] * h[j]);
            if (!(std::fabs(v) < 1e290)) return err;
            H[i][j] = H[j][i] = v;
        }
        H[i][n + i] = 1.0;
    }
    for (size_t c = 0; c < n; c++){ // Gauss-Jordan with row pivoting
        size_t piv = c;
        for (size_t r = c + 1; r < n; r++) if (std::fabs(H[r][c]) > std::fabs(H[piv][c])) piv = r;
        if (std::fabs(H[piv][c]) < 1e-300) return err;
        std::swap(H[c], H[piv]);
        const double d = H[c][c];
        for (size_t k = 0; k < 2 * n; k++) H[c][k] /= d;
        for (size_t r = 0; r < n; r++) if (r != c){const double m = H[r][c]; for (size_t k = 0; k < 2 * n; k++) H[r][k] -= m * H[c][k];}
    }
    for (size_t i = 0; i < n; i++) err[i] = H[i][n + i] > 0 ? std::sqrt(H[i][n + i]) : 0.0;
    return err;
}

struct UGaussFit {bool ok = false; std::vector<double> p, err; double n = 0;};
inline UGaussFit FitUGaus(const std::vector<double>& x, const UGaussModel& m, const std::vector<double>& p0, const std::vector<double>& step){
    UGaussFit r;
    r.n = (double)x.size();
    const auto nll = [&](const std::vector<double>& p){return UGaussNll(x, m, p);};
    r.p = NelderMead(nll, p0, step);
    if (!(nll(r.p) < 1e299)) return r;
    std::vector<double> h(step.size());
    for (size_t k = 0; k < h.size(); k++) h[k] = step[k] / 20;
    r.err = NllErrors(nll, r.p, h);
    r.ok = true;
    for (double v : r.p) if (!std::isfinite(v)) r.ok = false;
    return r;
}

// The fitted curve as a TF1 of the names of the old binned fits (Constant, Mean, Sigma, [Constant2, Mean2, Sigma2], [Background]) in the function list of h, so the stats box, the fit equation box
// and fit_values.csv read it like before. Height in counts per bin of h (it overlays the histogram), chi2 and NDF 0 (unbinned). Constant = mass * binWidth / (sigma sqrt(2 pi) Z).
inline TF1* AddGausCurve(TH1* h, const char* name, const UGaussFit& r, const UGaussModel& m){
    const double bw = h->GetBinWidth(1), nEff = r.n, root2pi = std::sqrt(2.0 * TMath::Pi());
    const double fb = m.flat ? r.p[m.nGaus == 2 ? 5 : 2] : 0.0, fbErr = m.flat ? r.err[m.nGaus == 2 ? 5 : 2] : 0.0;
    const double f2 = m.nGaus == 2 ? r.p[4] : 0.0;
    const char* form = m.nGaus == 2 ? "gaus(0) + gaus(3) + pol0(6)" : m.flat ? "gaus(0) + pol0(3)" : "gaus";
    TF1* f = new TF1(name, form, m.lo, m.hi);
    auto gauss = [&](int slot, double mu, double s, double muErr, double sErr, double mass){
        const double z = GausCdf((m.hi - mu) / s) - GausCdf((m.lo - mu) / s);
        const double c = mass * nEff * bw / (s * root2pi * z);
        f->SetParameter(slot, c); f->SetParError(slot, c / std::sqrt(std::max(nEff, 1.0)));
        f->SetParameter(slot + 1, mu); f->SetParError(slot + 1, muErr);
        f->SetParameter(slot + 2, s); f->SetParError(slot + 2, sErr);
    };
    gauss(0, r.p[0], r.p[1], r.err[0], r.err[1], (1 - fb) * (1 - f2));
    if (m.nGaus == 2){
        gauss(3, r.p[2], r.p[3], r.err[2], r.err[3], (1 - fb) * f2);
        f->SetParNames("Constant", "Mean", "Sigma", "Constant2", "Mean2", "Sigma2", "Background");
    } else if (m.flat) f->SetParNames("Constant", "Mean", "Sigma", "Background");
    else f->SetParNames("Constant", "Mean", "Sigma");
    if (m.flat){
        const int slot = m.nGaus == 2 ? 6 : 3;
        f->SetParameter(slot, nEff * bw * fb / (m.hi - m.lo)); f->SetParError(slot, nEff * bw * fbErr / (m.hi - m.lo));
    }
    f->SetNDF(0); f->SetChisquare(0);
    h->GetListOfFunctions()->Add(f);
    return f;
}

// Pearson chi2 of h against the curve f per degree of freedom, bins with at least 1 expected count: only used to decide whether the core fit needs a second gaussian, never to fit.
inline double CurveChi2PerNdf(TH1* h, TF1* f, int nPar){
    double chi2 = 0; int nb = 0;
    for (int b = 1; b <= h->GetNbinsX(); b++){
        const double x = h->GetBinCenter(b);
        if (x < f->GetXmin() || x > f->GetXmax()) continue;
        const double e = f->Eval(x);
        if (e < 1) continue;
        chi2 += (h->GetBinContent(b) - e) * (h->GetBinContent(b) - e) / e; nb++;
    }
    return nb > nPar ? chi2 / (nb - nPar) : 0.0;
}

// Block positions (fDx) are discrete cell centers (one bin of h per block, bin width = block pitch). Unbinned fit of the gaussian to the values x, each value with the probability of
// its cell (gaussian integrated over the block width), over the whole axis of h, seeded from the mean and RMS. No windowed refit, since a block cannot be split.
inline void FitBlocks(TH1* h, const std::vector<double>& x){
    const double bw = h->GetBinWidth(1), lo = h->GetXaxis()->GetXmin(), hi = h->GetXaxis()->GetXmax();
    std::vector<double> v;
    for (double a : x) if (a >= lo && a <= hi) v.push_back(a);
    if (v.size() < 3) return;
    double mean = 0, ms = 0;
    for (double a : v) mean += a / v.size();
    for (double a : v) ms += (a - mean) * (a - mean) / v.size();
    const double sigma = std::max(std::sqrt(ms), 0.5 * bw);
    UGaussModel m; m.lo = lo; m.hi = hi; m.cell = bw; m.pLo = {lo, 0.05 * bw}; m.pHi = {hi, hi - lo};
    const UGaussFit r = FitUGaus(v, m, {mean, sigma}, {0.1 * sigma, 0.1 * sigma});
    if (!r.ok) return;
    AddGausCurve(h, "fBlocks", r, m);
    AddFitEquation(h, {kGausEq, "unbinned fit", "integrated over each block"});
}

// Unbinned fit of the core of the values x of h: a gaussian on a window around the peak (mean +- 2.5 sigma, refitted until it settles), plus a flat background when the plot has one (phi
// difference wrapped at +-180), plus a second wider gaussian when that does not describe the shoulders. h is only used for the seed and the axis.
inline void FitCore(TH1* h, const std::vector<double>& x){
    if (h->GetEntries() < 3) return;

    // Seed from the peak bin and its FWHM, not the quantiles: a flat background (e.g. phi difference, wrapped at
    // +-180) inflates the IQR until the window covers the whole histogram and the fit goes broad.
    double bw = h->GetBinWidth(1);
    // The seed (peak bin, background, FWHM) is read from a lightly smoothed copy: with few entries one noisy bin
    // looked like a one-bin-wide peak and the fit locked onto it.
    std::unique_ptr<TH1> hs(static_cast<TH1*>(h->Clone("fitCoreSeed")));
    hs->SetDirectory(nullptr);
    hs->Smooth(2);
    int pk = hs->GetMaximumBin();
    double mean = hs->GetBinCenter(pk);

    // Flat background level = median bin content (robust to the peak). Only when it is a real fraction of the peak
    // (e.g. phi difference wrapped at +-180) is the background fitted; otherwise the plain gaussian is kept.
    std::vector<double> contents;
    for (int b = 1; b <= h->GetNbinsX(); b++) contents.push_back(hs->GetBinContent(b));
    std::nth_element(contents.begin(), contents.begin() + contents.size() / 2, contents.end());
    double bkg = contents[contents.size() / 2];
    bool useBkg = bkg > 0.05 * hs->GetBinContent(pk);
    if (!useBkg) bkg = 0;

    double half = bkg + 0.5 * (hs->GetBinContent(pk) - bkg);
    int lo = pk, hi = pk;
    while (lo > 1 && hs->GetBinContent(lo - 1) > half) lo--;
    while (hi < h->GetNbinsX() && hs->GetBinContent(hi + 1) > half) hi++;
    double sigma = std::max((h->GetBinCenter(hi) - h->GetBinCenter(lo) + bw) / 2.355, bw);

    const double axLo = h->GetXaxis()->GetXmin(), axHi = h->GetXaxis()->GetXmax();
    const auto inRange = [&](double a, double b){std::vector<double> v; for (double t : x) if (t >= a && t <= b) v.push_back(t); return v;};
    // the window fits: mean and sigma are refitted on mean +- halfWin until they stop moving
    UGaussFit r; UGaussModel m; m.flat = useBkg;
    for (int i = 0; i < 4; i++){
        const double halfWin = useBkg ? std::max(3.0 * sigma, 1.5 * bw) : std::max(2.5 * sigma, 1.5 * bw);
        m.lo = std::max(axLo, mean - halfWin); m.hi = std::min(axHi, mean + halfWin);
        const std::vector<double> v = inRange(m.lo, m.hi);
        if (v.size() < 5) return;
        m.pLo = {m.lo, 0.05 * bw}; m.pHi = {m.hi, m.hi - m.lo};
        std::vector<double> p0 = {mean, sigma}, step = {0.1 * sigma, 0.1 * sigma};
        if (useBkg){
            const double gMass = (hs->GetBinContent(pk) - bkg) * sigma * std::sqrt(2.0 * TMath::Pi()) / bw, fMass = bkg * (m.hi - m.lo) / bw;
            p0.push_back(std::min(0.9, std::max(0.01, fMass / (gMass + fMass)))); step.push_back(0.05);
            m.pLo.push_back(0.0); m.pHi.push_back(0.95);
        }
        r = FitUGaus(v, m, p0, step);
        if (!r.ok) return;
        const double mm = r.p[0], s = std::fabs(r.p[1]);
        if (s < bw || !std::isfinite(mm) || !std::isfinite(s)) return;
        const bool settled = std::fabs(mm - mean) < 0.05 * sigma && std::fabs(s - sigma) < 0.05 * sigma;
        mean = mm; sigma = s;
        if (settled) break;
    }
    TF1* g = AddGausCurve(h, useBkg ? "fCoreBkg" : "fGaus", r, m);
    if (!useBkg){
        AddFitEquation(h, {kGausEq, "unbinned fit"});
        DrawFitTails(h, g);
        return;
    }
    // A flat background cannot take up shoulders around the core, so the window fit pushes it up to the shoulder level.
    // When the fit is poor, refit the whole range with a second, wider gaussian for the shoulders: the flat part is then
    // set by the far wings, and the narrow gaussian stays the core (Constant, Mean, Sigma).
    if (CurveChi2PerNdf(h, g, 4) > 3){
        const double mm = r.p[0], s = std::fabs(r.p[1]), fb0 = r.p[2], range = axHi - axLo;
        UGaussModel m2; m2.lo = axLo; m2.hi = axHi; m2.nGaus = 2; m2.flat = true;
        m2.pLo = {axLo, bw, mm - 4 * s, 1.5 * s, 0.0, 0.0}; m2.pHi = {axHi, 2 * s, mm + 4 * s, 0.5 * range, 1.0, 0.95};
        const UGaussFit r2 = FitUGaus(inRange(axLo, axHi), m2, {mm, s, mm, 4 * s, 0.3, 0.3 * fb0}, {0.1 * s, 0.1 * s, s, 0.3 * s, 0.05, 0.02});
        if (r2.ok){
            h->GetListOfFunctions()->Remove(g); delete g;
            g = AddGausCurve(h, "fCore2Bkg", r2, m2);
        }
    }
    if (std::string(g->GetName()) == "fCore2Bkg")
        AddFitEquation(h, {kGausEq, "+ Constant2 #times exp(-(x - Mean2)^{2} / 2 Sigma2^{2})", "+ Background", "unbinned fit"});
    else
        AddFitEquation(h, {kGausEq, "+ Background", "unbinned fit"});
    DrawFitTails(h, g);
}

// The first line of a stats box is the histogram's name: "Fit results" when h has a fit function, else what (or "Statistics").
// Set it with h->SetName right before the canvas is painted and saved, and put the old name back after, so Draw ">>name" and gDirectory->Get still find h.
// Recon block phi sits on the block grid: tracks along the x axis pile up at exactly 0 (a spike in the first phi bin). Every fit of the recon block phi
// (the thetaReconBlock sets) only uses the protons with phi from kBlockPhiLo to kBlockPhiHi degrees, and says so on its plot.
constexpr double kBlockPhiLo = 0.5, kBlockPhiHi = 359.5;
inline bool InBlockPhi(double phi){return phi >= kBlockPhiLo && phi <= kBlockPhiHi;}
inline const char* BlockPhiNote(){return "block phi fit: 0.5 to 359.5 deg only (the spike at 0 is left out)";}

inline std::string StatsTitle(TH1* h, const std::string& what = ""){
    for (TObject* o : *h->GetListOfFunctions()) if (dynamic_cast<TF1*>(o)) return "Fit results";
    return what.empty() ? "Statistics" : what;
}

// Appends the fit of h (first fit function on it) to fit_values.csv in the folder of plotPath. The first write to a
// file in a run truncates it, so reruns do not pile up.
inline void WriteFitValues(TH1* h, const std::string& plotPath, std::set<std::string>& started){
    TF1* f = nullptr;
    for (TObject* o : *h->GetListOfFunctions()) {
        if ((f = dynamic_cast<TF1*>(o))) break;
    }
    if (!f) return;
    const std::string csv = std::string(gSystem->GetDirName(plotPath.c_str())) + "/fit_values.csv";
    FILE* outCsv = fopen(csv.c_str(), started.insert(csv).second ? "w" : "a");
    if (!outCsv) return;
    if (started.insert(csv + "#header").second) fprintf(outCsv, "plot,function,entries,chi2,ndf,parameter,value,error\n");
    for (int i = 0; i < f->GetNpar(); i++) {
        fprintf(outCsv, "%s,%s,%.0f,%.6g,%d,%s,%.8g,%.4g\n", gSystem->BaseName(plotPath.c_str()), f->GetName(), h->GetEntries(), f->GetChisquare(), f->GetNDF(), f->GetParName(i), f->GetParameter(i), f->GetParError(i));
    }
    fclose(outCsv);
}

// Appends entries, mean and standard deviation of h (x, and y too for 2D histograms) to stats.csv in the folder of
// plotPath. Same truncate-on-first-write rule as WriteFitValues.
inline void WriteHistStats(TH1* h, const std::string& plotPath, std::set<std::string>& started){
    const std::string file = std::string(gSystem->GetDirName(plotPath.c_str())) + "/stats.csv";
    FILE* out = fopen(file.c_str(), started.insert(file).second ? "w" : "a");
    if (!out) return;
    if (started.insert(file + "#header").second) fprintf(out, "plot,entries,mean,sd,mean_y,sd_y\n");
    fprintf(out, "%s,%.0f,%.8g,%.8g", gSystem->BaseName(plotPath.c_str()), h->GetEntries(), h->GetMean(1), h->GetStdDev(1));
    if (h->GetDimension() == 2) fprintf(out, ",%.8g,%.8g\n", h->GetMean(2), h->GetStdDev(2));
    else fprintf(out, ",,\n");
    fclose(out);
}

// Appends one row (the numbers drawn on a graph, not a histogram) to stats.csv in the folder of plotPath, after a header
// row written once per file. Same truncate-on-first-write rule as WriteHistStats; header and row are comma-joined strings.
// name: the csv in the plot's folder, for a second table with its own header next to stats.csv
inline void WriteStatsRow(const std::string& plotPath, const std::string& header, const std::string& row, std::set<std::string>& started, const std::string& name = "stats.csv"){
    const std::string file = std::string(gSystem->GetDirName(plotPath.c_str())) + "/" + name;
    FILE* out = fopen(file.c_str(), started.insert(file).second ? "w" : "a");
    if (!out) return;
    if (started.insert(file + "#header").second) fprintf(out, "%s\n", header.c_str());
    fprintf(out, "%s\n", row.c_str());
    fclose(out);
}

// Every number worth keeping from one saved plot: stats.csv always, fit_values.csv when it has a fit.
inline void WritePlotNumbers(TH1* h, const std::string& plotPath, std::set<std::string>& started){
    WriteHistStats(h, plotPath, started);
    WriteFitValues(h, plotPath, started);
}

// The side box with the cut numbers on every 6_cuts pdf and the crop of every pdf to its content (cut_side_box.py, crop_pdfs.py in the run folder, plot.sh runs them last
// too): called at the end of the macros that write 6_cuts plots, so a single macro run leaves finished plots.
inline void BoxAndCropPlots(){
    if (std::getenv("PLOT_PIPELINE")) return; // plot.sh / Plot.pbs run both scripts after all lanes finished; running them here races with other lanes
    for (const char* script : {"cut_side_box.py", "crop_pdfs.py"})
        if (!gSystem->AccessPathName(script)) gSystem->Exec((std::string("python3 -I ") + script + " plots").c_str());
}
