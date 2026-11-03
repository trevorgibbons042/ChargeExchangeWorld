#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <TTree.h>

#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>

#include "BlockSquares.h"
#include "CutOpen.h"
#include "DrawPlan.h"
#include "FitCore.h"
#include "CutFiles.h"

// Protons whose recon block phi is off from the CE phi by absLo to absHi degrees either way, -absHi to -absLo or absLo to absHi (CalcValsReDDIFF.phiReD in cut_APandT.root).
// The block phi diff window is the only cut: no KE cut, no block residual cut. The xy plots need a ZDC position, so they hold the ZDC hits only.
// Their ZDC positions come from cut_ParentZdc.root, matched on (fEvent, ParentID).
// energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the input trees like PlotAngle.cpp. Next to each window plot the same for every recon
// proton, and the block fit residual, R^2 and impact of every recon proton and the CE ones.
void PlotDiffPhiCuts(double absLo = 50, double absHi = 180, std::string energy = "fKEnergy"){
    FastJit();
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    const std::string folder = std::string("plots/4_CEblock/3_diffphiBlock/");
    gSystem->mkdir(folder.c_str(), true);
    gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X");

    TFile* fAP = TFile::Open("cuts/cut_APandT.root", "READ");
    TTree* diff = (TTree*)fAP->Get(("CalcValsReDDIFF" + suffix).c_str());
    // ZDC protons of the CE children, cuts/cut_ParentZdc.root (CutOpen.h): every row passes the ZDC proton cut
    TTree* pos = OpenCut("cuts/cut_ParentZdc.root", "ParticleInfo");
    if (!diff || !pos){Error("PlotDiffPhiCuts", "missing tree in cut_APandT.root or cut_ParentZdc.root"); return;}

    // Every column pull of this macro is booked here and collected in one multithreaded pass per tree (DrawPlan.h), the render code below reads the columns in entry order.
    ROOT::EnableImplicitMT(NThreads());
    DrawPlan plan;
    plan.unordered = true; // nothing below depends on the row order, so AnglesRaw (many clusters) is read on all threads
    // CE events in the phi diff window, either sign: the only cut
    auto cDiff = plan.C(diff, {"fEvent", "TrackID"}, Form("fabs(phiReD) >= %g && fabs(phiReD) <= %g", absLo, absHi));
    // Every recon proton (the rows of AnglesRaw in cut_angles_all.root, all events): (fEvent, TrackID) for the positions,
    // the KE of the energy (finalKE / finalKED), isFromCE and the block line fit quality (ReconAngleAll.cpp)
    TFile* fAng = TFile::Open("cuts/cut_angles_all.root", "READ");
    TTree* ang = fAng ? (TTree*)fAng->Get(kAngleTree) : nullptr;
    if (!ang || !ang->GetBranch("r2D")){Error("PlotDiffPhiCuts", "AnglesRaw with r2D missing in cut_angles_all.root (rerun ReconAngleAll.cpp)"); return;}
    const std::string angSel = AngleSelection("Angles"); // the same rows for both column sets, so they line up
    auto cAng = plan.C(ang, {"fEvent", "TrackID", KeVar(suffix), "isFromCE", "resDRms", "outPhi", "outDPhi", "ParentID", "PDG"}, angSel);
    DrawPlan::Cols cEx;
    if (ang->GetBranch("nProtonsEv")) cEx = plan.C(ang, {"resRms", "nProtonsEv", "resHalfRatioD", "outTheta*180/TMath::Pi()", "outDTheta*180/TMath::Pi()", "nHits"}, angSel);
    // ZDC protons of the CE children (pos) and of every event (allEv), the cut files: every row passes the ZDC proton cut
    TTree* allEv = OpenCut("cuts/cut_AllEvProtons.root", "ParticleInfo");
    if (!allEv){Error("PlotDiffPhiCuts", "cuts/cut_AllEvProtons.root missing (run CutEvents.cpp)"); return;}
    auto cPos = plan.C(pos, {"fEvent", "ParentID", "fDx", "fDy", "fX", "fY", "TrackID"});
    auto cAllEv = plan.C(allEv, {"fEvent", "TrackID", "fDx", "fDy", "fX", "fY"});
    plan.Run();
    gROOT->cd();
    const Long64_t nSel = cDiff.N();
    const time_t tKeys = time(nullptr);
    // Key sets: sorted vectors of packed (fEvent, TrackID) keys, not std::set, which is slow and large for every recon proton of the run. count() takes the pair the old sets took.
    struct KeySet {
        std::vector<uint64_t> k;
        static uint64_t Pack(long long ev, long long tr){return (uint64_t(uint32_t(ev)) << 32) | uint32_t(tr);}
        void Finish(){std::sort(k.begin(), k.end()); k.erase(std::unique(k.begin(), k.end()), k.end());}
        size_t count(const std::pair<long long, long long>& p) const {return std::binary_search(k.begin(), k.end(), Pack(p.first, p.second));}
        size_t size() const {return k.size();}
    };
    KeySet keys;
    for (Long64_t i = 0; i < nSel; i++){
        keys.k.push_back(KeySet::Pack(llround(cDiff.V(0)[i]), llround(cDiff.V(1)[i])));
    }
    keys.Finish();

    const Long64_t nAng = cAng.N();
    std::vector<const double*> angCols(9); // views of the DrawPlan columns, no copy
    for (int c = 0; c < 9; c++) angCols[c] = cAng.V(c);
    // every recon proton (all events, CE ones included) in the same block phi diff window: the CE protons keep their block phi - CE phi
    // window above, the non-CE ones have no CE phi so their diff is recon block phi - recon tracker phi (outPhi), wrapped to -180 to 180
    // like CalcAPandT.cpp. ceWinKeys: the window CE protons as (fEvent, proton TrackID), the Angles tree key: AnglesRaw rows from a CE neutron in keys (ParentID)
    KeySet reconKeys, ceWinKeys;
    struct AngKeys {std::vector<uint64_t> ce, recon;};
    const auto angParts = MapRows(nAng, [&](Long64_t lo, Long64_t hi){ // row ranges on all threads (MapRows, FitCore.h), each collects its own keys
        AngKeys part;
        for (Long64_t i = lo; i < hi; i++){
            if (angCols[3][i] > 0.5){
                if (keys.count({llround(angCols[0][i]), llround(angCols[7][i])})) part.ce.push_back(KeySet::Pack(llround(angCols[0][i]), llround(angCols[1][i])));
                continue;
            }
            const double dPhi = std::remainder((angCols[6][i] - angCols[5][i]) * 180 / TMath::Pi(), 360.0);
            if (fabs(dPhi) >= absLo && fabs(dPhi) <= absHi) part.recon.push_back(KeySet::Pack(llround(angCols[0][i]), llround(angCols[1][i])));
        }
        return part;
    });
    for (const auto& part : angParts){
        ceWinKeys.k.insert(ceWinKeys.k.end(), part.ce.begin(), part.ce.end());
        reconKeys.k.insert(reconKeys.k.end(), part.recon.begin(), part.recon.end());
    }
    ceWinKeys.Finish(); reconKeys.Finish();

    // positions at the ZDC, block (fDx, fDy) and true (fX, fY): the window protons from cut_ParentZdc.root (cut_APandT.root TrackID is the
    // CE neutron's track, so it matches the proton's ParentID there), every recon proton from cut_AllEvProtons.root (matched on TrackID)
    // the block maps have one bin per ZDC cell (nCols_detector x nRows_detector over 400 mm), drawn as one square per cell (BlockSquares.h); the true ones keep 100x100
    const int nColsD = (int)ReadSettings().at("nCols_detector"), nRowsD = (int)ReadSettings().at("nRows_detector");
    TH2D hBlock("hBlock", "", nColsD, -200, 200, nRowsD, -200, 200), hTrue("hTrue", "", 100, -200, 200, 100, -200, 200);
    TH2D hBlockAll("hBlockAll", "", nColsD, -200, 200, nRowsD, -200, 200), hTrueAll("hTrueAll", "", 100, -200, 200, 100, -200, 200);
    // the rows of a column set whose (fEvent, id) is in a key set, found on all threads; the histograms are then filled here in row order, so the stats sums are those of a serial fill
    auto matching = [&](const DrawPlan::Cols& c, const KeySet& set){
        const auto parts = MapRows(c.N(), [&](Long64_t lo, Long64_t hi){
            std::vector<Long64_t> rows;
            for (Long64_t i = lo; i < hi; i++) if (set.count({llround(c.V(0)[i]), llround(c.V(1)[i])})) rows.push_back(i);
            return rows;
        });
        std::vector<Long64_t> rows;
        for (const auto& p : parts) rows.insert(rows.end(), p.begin(), p.end());
        return rows;
    };
    for (Long64_t i : matching(cPos, keys)){
        hBlock.Fill(cPos.V(2)[i], cPos.V(3)[i]);
        hTrue.Fill(cPos.V(4)[i], cPos.V(5)[i]);
    }
    for (Long64_t i : matching(cAllEv, reconKeys)){
        hBlockAll.Fill(cAllEv.V(2)[i], cAllEv.V(3)[i]);
        hTrueAll.Fill(cAllEv.V(4)[i], cAllEv.V(5)[i]);
    }
    printf(">>> PlotDiffPhiCuts key sets and ZDC fills took %.0f s\n", difftime(time(nullptr), tKeys)); fflush(stdout);
    hBlockAll.Add(&hBlock);
    hTrueAll.Add(&hTrue);

    TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
    gPad->SetRightMargin(0.30);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);

    const std::string cutLabel = Form("CE recon protons, block phi diff -%g to -%g and %g to %g deg, no other cut, N = %lld", absHi, absLo, absLo, absHi, nSel);
    std::set<std::string> numbersStarted;
    auto save = [&](TH1* h, const std::string& title, const std::string& file, const char* opt, bool squares = false){
        h->SetTitle(PlotTitle(file, title).c_str());
        canvas->Clear();
        h->Draw(opt);
        if (squares) DrawBlockSquares(h, 0, false);
        WritePlotNumbers(h, folder + file, numbersStarted);
        canvas->Update();
        canvas->SaveAs((folder + file).c_str());
    };
    const std::string xy = ";X Position (mm);Y Position (mm)";
    save(&hBlock, "Position at ZDC BLOCK, " + cutLabel + ", hit ZDC (position), 1 square per ZDC cell" + xy, "1_xy_block.pdf", "colz", true);
    const std::string reconLabel = Form("all recon protons, block phi diff -%g to -%g and %g to %g deg (CE: block - CE, non-CE: block - tracker), no other cut, N = %lld", absHi, absLo, absLo, absHi,
                                            nSel + (Long64_t)reconKeys.size());
    save(&hBlockAll, "Position at ZDC BLOCK, " + reconLabel + ", 1 square per ZDC cell" + xy, "2_xy_block_recon.pdf", "colz", true);

    // Recon histograms on log y: all recon particles (green, every row of AnglesRaw), all recon protons (blue, pdg 2212) and the CE recon protons in the phi diff window (red)
    auto drawOverlay = [&](TH1D& hRecon, TH1D& hProt, TH1D& hCE, TH1D* hWin, const std::string& label, const std::string& what, const std::string& file,
                           const std::vector<std::string>& eqLines){
        hRecon.SetLineColor(kGreen + 2); hProt.SetLineColor(kBlue); hCE.SetLineColor(kRed);
        for (TH1D* h : {&hRecon, &hProt, &hCE}) h->SetLineWidth(2);
        if (hWin){hWin->SetLineColor(kRed); hWin->SetLineWidth(2);}
        hRecon.SetStats(0);
        hRecon.SetMinimum(0.5); // log y: the recon protons outnumber the CE ones by far
        hRecon.SetMaximum(2.0 * std::max({hRecon.GetMaximum(), hProt.GetMaximum(), hCE.GetMaximum(), hWin ? hWin->GetMaximum() : 0.0}));
        hRecon.SetTitle(PlotTitle(file, label + ", " + what + ";" + label + ";Count").c_str());
        canvas->Clear();
        canvas->SetLogy(1);
        hRecon.Draw(); hProt.Draw("same"); hCE.Draw("same"); if (hWin) hWin->Draw("same");
        TLegend leg(0.745, 0.72, 0.99, 0.88); leg.SetTextSize(0.014);
        leg.AddEntry(&hRecon, Form("all recon particles, no cuts (%.0f)", hRecon.GetEntries()), "l");
        leg.AddEntry(&hProt, Form("all recon protons, no cuts (%.0f)", hProt.GetEntries()), "l");
        leg.AddEntry(&hCE, Form("CE recon protons in window (%.0f, %.2f%%)", hCE.GetEntries(), hRecon.GetEntries() > 0 ? 100.0 * hCE.GetEntries() / hRecon.GetEntries() : 0.0), "l");
        if (hWin) leg.AddEntry(hWin, Form("phi diff window (%.0f)", hWin->GetEntries()), "l");
        leg.Draw();
        TPaveText eq(0.745, 0.70 - 0.05 * (eqLines.size() + 1), 0.99, 0.70, "NDC");
        eq.SetBorderSize(1); eq.SetFillColor(0); eq.SetTextAlign(12); eq.SetTextSize(0.016);
        for (const auto& l : eqLines) eq.AddText(l.c_str());
        if (!eqLines.empty()) eq.Draw();
        WritePlotNumbers(&hRecon, folder + file, numbersStarted);
        canvas->Update();
        canvas->SaveAs((folder + file).c_str());
        canvas->SetLogy(0);
    };

    // Per row, once for all the overlays below: bit 1 = a proton (pdg 2212), bit 2 = a CE proton in the window (ceWinKeys). Row ranges on all threads, each writes its own part.
    std::vector<unsigned char> rowFlag(nAng);
    MapRows(nAng, [&](Long64_t lo, Long64_t hi){
        for (Long64_t i = lo; i < hi; i++)
            rowFlag[i] = (llround(angCols[8][i]) == 2212 ? 1 : 0) | (ceWinKeys.count({llround(angCols[0][i]), llround(angCols[1][i])}) ? 2 : 0);
        return 0;
    });
    // Fills hRecon (every row with a finite x), hProt (the protons) and hCE (the CE ones in the window) from the column x. Each row range fills its own copies on all threads (MapRows), added
    // into the three histograms in range order after: the bin contents are exact, the mean and RMS sums can differ from a serial fill in the last digits.
    struct Hist3 {std::shared_ptr<TH1D> recon, prot, ce;};
    auto fill3 = [&](const double* x, TH1D& hRecon, TH1D& hProt, TH1D& hCE){
        const bool addDir = TH1::AddDirectoryStatus(); TH1::AddDirectory(false); // the copies of the ranges are not registered in the current directory
        const auto parts = MapRows(nAng, [&](Long64_t lo, Long64_t hi){
            Hist3 p{std::make_shared<TH1D>(), std::make_shared<TH1D>(), std::make_shared<TH1D>()};
            hRecon.Copy(*p.recon); hProt.Copy(*p.prot); hCE.Copy(*p.ce); // same binning, still empty
            for (Long64_t i = lo; i < hi; i++){
                if (!std::isfinite(x[i])) continue;
                p.recon->Fill(x[i]);
                if (rowFlag[i] & 1) p.prot->Fill(x[i]);
                if (rowFlag[i] & 2) p.ce->Fill(x[i]);
            }
            return p;
        });
        TH1::AddDirectory(addDir);
        for (const Hist3& p : parts){hRecon.Add(p.recon.get()); hProt.Add(p.prot.get()); hCE.Add(p.ce.get());}
    };

    // KE: every recon proton and the CE ones in the phi diff window (finalKE of the Angles tree)
    {
        TH1D hKERecon("hKERecon", "", 100, 0, 30000), hKEProt("hKEProt", "", 100, 0, 30000), hKECE("hKECE", "", 100, 0, 30000);
        fill3(angCols[2], hKERecon, hKEProt, hKECE); // no ZDC hit: no KE, not filled
        drawOverlay(hKERecon, hKEProt, hKECE, nullptr, "Kinetic Energy at ZDC (MeV)", Form("all recon particles, all recon protons and CE recon protons, block phi diff -%g to -%g and %g to %g deg", absHi, absLo, absLo, absHi), "5_ke.pdf", {});
    }

    // Block line fit quality of every recon proton and the CE ones in the phi diff window: 6_resDRms.pdf
    struct QualVar {std::string label, file; double lo, hi; std::vector<std::string> eq;}; // eq: short lines, they sit in the right margin
    const std::string line = "fit to block hits: #hat{x} = x_{0} + b_{x}z, #hat{y} = y_{0} + b_{y}z";
    const std::string ssRes = "SS_{res} = #Sigma_{i} (x_{i} - #hat{x}_{i})^{2} + (y_{i} - #hat{y}_{i})^{2}";
    const std::vector<QualVar> vars = { // same order as column 4 of angCols
        {"block residual (mm)", "6_resDRms.pdf", 0, 5, {line, "residual = #sqrt{SS_{res} / n}", ssRes}}};
    for (size_t v = 0; v < vars.size(); v++){
        TH1D hRecon(Form("hRecon%zu", v), "", 100, vars[v].lo, vars[v].hi), hProt(Form("hProt%zu", v), "", 100, vars[v].lo, vars[v].hi), hCE(Form("hCE%zu", v), "", 100, vars[v].lo, vars[v].hi);
        fill3(angCols[4 + v], hRecon, hProt, hCE);
        drawOverlay(hRecon, hProt, hCE, nullptr, vars[v].label, Form("all recon particles, all recon protons and CE recon protons, block phi diff -%g to -%g and %g to %g deg", absHi, absLo, absLo, absHi), vars[v].file, vars[v].eq);
    }

    // The other cut candidates (6_cuts/3_ResCuts to 5_nHits), shown only, nothing is cut on them: the true position residual, protons per event,
    // residual of the 2nd half over the 1st half of the tracker hits, recon theta, recon block theta and tracker hits. 9_ to 17_
    if (ang->GetBranch("nProtonsEv")){
        const Long64_t nEx = cEx.N();
        if (nEx != nAng){Error("PlotDiffPhiCuts", "extra columns have %lld rows, not %lld", nEx, nAng); return;}
        struct ExVar {std::string label, file; double lo, hi; int nBins; std::vector<std::string> eq;};
        const std::string lineT = "fit to hits: #hat{x} = x_{0} + b_{x}z, #hat{y} = y_{0} + b_{y}z";
        const std::vector<ExVar> exVars = {
            {"residual, true position (mm)", "9_resRms.pdf", 0, 5, 100, {lineT, "residual = #sqrt{SS_{res} / n}"}},
            {"reconstructed particles in the event", "12_nParticles.pdf", -0.5, 19.5, 20, {}},
            {"block residual, hits after first 2 / first 2 hits (>2 hits)", "14_resHalf.pdf", 0, 10, 100, {}},
            {"recon theta (degrees)", "15_theta.pdf", 0, 2.3, 100, {}},
            {"recon block theta (degrees)", "16_thetaBlock.pdf", 0, 2.3, 100, {}},
            {"tracker hits", "17_nHits.pdf", -0.5, 10.5, 11, {}}};
        std::vector<const double*> exCols(exVars.size());
        for (size_t c = 0; c < exVars.size(); c++) exCols[c] = cEx.V(c);
        for (size_t v = 0; v < exVars.size(); v++){
            TH1D hRecon(Form("hExRecon%zu", v), "", exVars[v].nBins, exVars[v].lo, exVars[v].hi), hProt(Form("hExProt%zu", v), "", exVars[v].nBins, exVars[v].lo, exVars[v].hi), hCE(Form("hExCE%zu", v), "", exVars[v].nBins, exVars[v].lo, exVars[v].hi);
            fill3(exCols[v], hRecon, hProt, hCE);
            drawOverlay(hRecon, hProt, hCE, nullptr, exVars[v].label, Form("all recon particles, all recon protons and CE recon protons, block phi diff -%g to -%g and %g to %g deg", absHi, absLo, absLo, absHi), exVars[v].file, exVars[v].eq);
        }
    }
    BoxAndCropPlots();
}
