#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <set>
#include <string>
#include <utility>

#include "FitCore.h"

// Protons whose recon block phi is off from the CE phi by phiLo to phiHi degrees (CalcValsReDDIFF.phiReD in cut_APandT.root).
// Their ZDC positions come from cut_ParentID.root, matched on (fEvent, ParentID); their KE is finalKE from cut_APandT.root.
// energy: "fKEnergy" (true KE) or "fDEnergy" (detector energy), picks the input trees like PlotAngle.cpp.
void PlotDiffPhiCuts(double phiLo = -180, double phiHi = -50, std::string energy = "fKEnergy"){
    const std::string suffix = (energy == "fDEnergy") ? "_DE" : "";
    const std::string folder = std::string("plots/4_ReconCalcVals/3_diffphi_cuts/") + (energy == "fDEnergy" ? "2_KEDetector" : "1_KEtrue") + "/";
    gSystem->mkdir(folder.c_str(), true);
    gStyle->SetStatX(0.99); gStyle->SetStatW(0.19); gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255); gStyle->SetNdivisions(505, "X");

    TFile* fAP = TFile::Open("cut_APandT.root", "READ");
    TTree* diff = (TTree*)fAP->Get(("CalcValsReDDIFF" + suffix).c_str());
    TFile* fPos = TFile::Open("cut_ParentID.root", "READ");
    TTree* pos = (TTree*)fPos->Get("ParticleInfo");
    if (!diff || !pos){Error("PlotDiffPhiCuts", "missing tree in cut_APandT.root or cut_ParentID.root"); return;}

    // events in the phi diff window
    diff->SetEstimate(diff->GetEntries() + 1);
    Long64_t nSel = diff->Draw("fEvent:TrackID:finalKE", Form("phiReD >= %g && phiReD < %g", phiLo, phiHi), "goff");
    std::set<std::pair<long long, long long>> keys;
    TH1D hKE("hKE", "", 100, 0, 30000);
    for (Long64_t i = 0; i < nSel; i++){
        keys.insert({llround(diff->GetV1()[i]), llround(diff->GetV2()[i])});
        hKE.Fill(diff->GetV3()[i]);
    }

    // their positions at the ZDC: block (fDx, fDy) and true (fX, fY)
    pos->SetEstimate(pos->GetEntries() + 1);
    // cut_APandT.root TrackID is the CE neutron's track, so it matches the proton's ParentID here
    Long64_t nPos = pos->Draw("fEvent:ParentID:fDx:fDy:fX:fY", "boundaryType == 2 && pdg == 2212", "goff");
    TH2D hBlock("hBlock", "", 100, -200, 200, 100, -200, 200);
    TH2D hTrue("hTrue", "", 100, -200, 200, 100, -200, 200);
    for (Long64_t i = 0; i < nPos; i++){
        if (!keys.count({llround(pos->GetVal(0)[i]), llround(pos->GetVal(1)[i])})) continue;
        hBlock.Fill(pos->GetVal(2)[i], pos->GetVal(3)[i]);
        hTrue.Fill(pos->GetVal(4)[i], pos->GetVal(5)[i]);
    }

    TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
    gPad->SetRightMargin(0.30);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);

    const std::string cutLabel = Form("CE events CE protons, block phi diff %g to %g deg, N = %lld", phiLo, phiHi, nSel);
    std::set<std::string> numbersStarted;
    auto save = [&](TH1* h, const std::string& title, const std::string& file, const char* opt){
        h->SetTitle(PlotTitle(file, title).c_str());
        canvas->Clear();
        h->Draw(opt);
        WritePlotNumbers(h, folder + file, numbersStarted);
        canvas->Update();
        canvas->SaveAs((folder + file).c_str());
    };
    save(&hBlock, "Position at ZDC BLOCK, " + cutLabel + ";X Position (mm);Y Position (mm)", "1_xy_block.pdf", "colz");
    save(&hTrue, "Position at ZDC, " + cutLabel + ";X Position (mm);Y Position (mm)", "2_xy.pdf", "colz");
    save(&hKE, "Kinetic Energy at ZDC, " + cutLabel + ";Kinetic Energy (MeV);Count", "3_ke.pdf", "");
}
