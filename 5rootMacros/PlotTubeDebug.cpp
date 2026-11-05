#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>
#include <TBox.h>
#include <TGraph.h>
#include <TLine.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "FitCore.h"

// Troubleshooting plot of the vertex tube of ReconAngleAll.cpp (the tube* and *N branches of AnglesRaw). Not part of run.sh or plot.sh: run it by hand from a run folder, after ReconAngleAll.
//   root -l -b -q 'PlotTubeDebug.cpp'                    two example CE protons (different tube blocks, 5+ tracker hits)
//   root -l -b -q 'PlotTubeDebug.cpp(26079, 15)'         one track: event 26079, TrackID 15 (any proton of AnglesRaw, not only CE ones)
//   root -l -b -q 'PlotTubeDebug.cpp(-1, -1, 3, 6)'      3 examples with 6+ tracker hits
// Per track, x(z) and y(z) side by side: the full stack (carbon blocks gray, tracker planes green dashed, tracker block hits, the block-position fit dashed blue, the tube pink, the true vertex a red star)
// and below it a zoom on the tube block, with the tube at its true width. The text on the right gives the true vertex, the tube, the neutron theta and phi range and the true neutron direction.
// The pdf goes to tube_debug/ in the run folder. python3 crop_pdfs.py tube_debug crops it to its content (ROOT puts the figure in a corner of an A4 page).

// The tracker block hits of a track from the raw ParticleInfo: u = block position (fDx for a = 0, fDy for a = 1) at z = fDz.
// ut, zt: the same hits at their true position (fX or fY at fZ).
static void TubeHits(TChain& pi, int ev, int tr, int a, std::vector<double>& u, std::vector<double>& z, std::vector<double>& ut, std::vector<double>& zt){
    u.clear(); z.clear(); ut.clear(); zt.clear();
    const std::string cut = "fEvent==" + std::to_string(ev) + "&&TrackID==" + std::to_string(tr) + "&&boundaryType==0";
    pi.Draw(a ? "fDy:fDz" : "fDx:fDz", cut.c_str(), "goff");
    for (Long64_t i = 0; i < pi.GetSelectedRows(); i++){u.push_back(pi.GetV1()[i]); z.push_back(pi.GetV2()[i]);}
    pi.Draw(a ? "fY:fZ" : "fX:fZ", cut.c_str(), "goff");
    for (Long64_t i = 0; i < pi.GetSelectedRows(); i++){ut.push_back(pi.GetV1()[i]); zt.push_back(pi.GetV2()[i]);}
}

void PlotTubeDebug(int event = -1, int trackID = -1, int nTracks = 2, int minHits = 5){
    const auto s = ReadSettings();
    const StackLayout layout = StackZ(s);
    const double half = s.at("CarbonThickness") * 500.0, zdc = 7700.0;
    const double rad2d = 180.0 / M_PI;
    constexpr double kTubeSigmas = 2.0; // the tube radius is this many sigma of the fit, the same number as kTubeSigmas in ReconAngleAll.cpp

    TFile fa("cuts/cut_angles_all.root");
    TTree* t = (TTree*)fa.Get("AnglesRaw");
    if (!t){Error("PlotTubeDebug", "no AnglesRaw in cuts/cut_angles_all.root (run ReconAngleAll.cpp)"); return;}
    if (!t->GetBranch("tubeRadius")){Error("PlotTubeDebug", "AnglesRaw has no tube branches (rerun ReconAngleAll.cpp)"); return;}
    TChain pi("ParticleInfo");
    pi.Add("output0_t*.root");

    int fromCE, evn, trk, par, blk, nh;
    double z0, z1, x0, y0, x1, y1, rad, xv, yv, zv, x0D, bxD, y0D, byD, thLo, thHi, phMid, phHalf;
    t->SetBranchAddress("isFromCE", &fromCE); t->SetBranchAddress("fEvent", &evn); t->SetBranchAddress("TrackID", &trk); t->SetBranchAddress("ParentID", &par);
    t->SetBranchAddress("tubeBlock", &blk); t->SetBranchAddress("nHits", &nh);
    t->SetBranchAddress("tubeZ0", &z0); t->SetBranchAddress("tubeZ1", &z1); t->SetBranchAddress("tubeX0", &x0); t->SetBranchAddress("tubeY0", &y0);
    t->SetBranchAddress("tubeX1", &x1); t->SetBranchAddress("tubeY1", &y1); t->SetBranchAddress("tubeRadius", &rad);
    t->SetBranchAddress("xVertexOut", &xv); t->SetBranchAddress("yVertexOut", &yv); t->SetBranchAddress("zVertexOut", &zv);
    t->SetBranchAddress("x0DOut", &x0D); t->SetBranchAddress("bxDOut", &bxD); t->SetBranchAddress("y0DOut", &y0D); t->SetBranchAddress("byDOut", &byD);
    t->SetBranchAddress("thetaNMin", &thLo); t->SetBranchAddress("thetaNMax", &thHi); t->SetBranchAddress("phiNMid", &phMid); t->SetBranchAddress("phiNHalf", &phHalf);

    // the true neutron direction at the CE: (event, CE neutron TrackID) -> direction, the neutron TrackID is the ParentID of its CE proton
    std::map<std::pair<int, int>, std::array<double, 3>> nDir;
    if (TFile* fc = TFile::Open("cuts/cut_CE.root")){
        if (TTree* ce = (TTree*)fc->Get("ChargeExchange_values")){
            int e, id; double dx, dy, dz;
            ce->SetBranchAddress("fEvent", &e); ce->SetBranchAddress("TrackID", &id);
            ce->SetBranchAddress("CENDirX", &dx); ce->SetBranchAddress("CENDirY", &dy); ce->SetBranchAddress("CENDirZ", &dz);
            for (Long64_t i = 0; i < ce->GetEntries(); i++){ce->GetEntry(i); nDir[{e, id}] = {dx, dy, dz};}
        }
        fc->Close();
    }

    // the tracks to draw: the one asked for, or examples, each with its own tube block
    std::vector<Long64_t> picks;
    {
        int lastBlk = -1;
        for (Long64_t i = 0; i < t->GetEntries() && (int)picks.size() < (event >= 0 ? 6 : nTracks); i++){
            t->GetEntry(i);
            if (!(rad > 0)) continue; // no tube: no tracker hit after a block
            if (event >= 0){if (evn == event && (trackID < 0 || trk == trackID)) picks.push_back(i);}
            else if (fromCE && nh >= minHits && blk != lastBlk){picks.push_back(i); lastBlk = blk;}
        }
    }
    if (picks.empty()){Error("PlotTubeDebug", "no track with a tube found for these arguments"); return;}
    const int nEx = (int)picks.size();

    // the full-view panels of every track share one axis range per x and per y, the zoom panels one span (x and y alike, each starts where its data is)
    double fullMax[2] = {0, 0}, zoomSpan[2] = {0, 0};
    for (const Long64_t i : picks){
        t->GetEntry(i);
        for (int a = 0; a < 2; a++){
            std::vector<double> u, z, ut, zt; TubeHits(pi, evn, trk, a, u, z, ut, zt);
            const double u0 = a ? y0 : x0, u1 = a ? y1 : x1, uv = a ? yv : xv;
            double m = std::fabs(uv); for (double uu : u) m = std::max(m, std::fabs(uu));
            fullMax[a] = std::max(fullMax[a], 1.12 * m + 0.2);
            double lo = std::min(u0, u1) - rad, hi = std::max(u0, u1) + rad;
            lo = std::min(lo, uv); hi = std::max(hi, uv);
            for (size_t k = 0; k < u.size(); k++) if (z[k] >= z0 - 20 && z[k] <= z1 + 40){lo = std::min(lo, u[k]); hi = std::max(hi, u[k]);}
            zoomSpan[a] = std::max(zoomSpan[a], 1.25 * (hi - lo));
        }
    }
    zoomSpan[0] = zoomSpan[1] = std::max(zoomSpan[0], zoomSpan[1]);

    gStyle->SetOptStat(0);
    const double canvasH = 1150.0 * nEx, footFrac = 100.0 / canvasH;
    TCanvas c("c", "", 2000, (int)canvasH);
    auto main = new TPad("m", "", 0, footFrac, 0.62, 1); main->Draw(); main->SetFillStyle(0); main->Divide(2, 2 * nEx);
    auto tp = new TPad("tp", "", 0.62, footFrac, 1, 1); tp->Draw(); tp->SetFillStyle(0);

    for (int k = 0; k < nEx; k++){
        t->GetEntry(picks[k]);
        for (int zoom = 0; zoom < 2; zoom++) for (int a = 0; a < 2; a++){
            main->cd(4 * k + 2 * zoom + a + 1);
            gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.13);
            std::vector<double> u, z, ut, zt; TubeHits(pi, evn, trk, a, u, z, ut, zt);
            const char* nm = a ? "y" : "x";
            const double u0 = a ? y0 : x0, u1 = a ? y1 : x1, uv = a ? yv : xv, ui = a ? y0D : x0D, bi = a ? byD : bxD;
            const double zlo = zoom ? z0 - 20 : 0, zhi = zoom ? z1 + 40 : zdc;
            double lo, hi;
            if (zoom){
                double l = std::min(u0, u1) - rad, h = std::max(u0, u1) + rad;
                l = std::min(l, uv); h = std::max(h, uv);
                for (size_t j = 0; j < u.size(); j++) if (z[j] >= zlo && z[j] <= zhi){l = std::min(l, u[j]); h = std::max(h, u[j]);}
                const double mid = 0.5 * (l + h); lo = mid - zoomSpan[a] / 2; hi = mid + zoomSpan[a] / 2;
            } else {lo = -fullMax[a]; hi = fullMax[a];}
            auto fr = new TH1F("", "", 10, zlo, zhi); fr->SetMinimum(lo); fr->SetMaximum(hi);
            fr->SetTitle(Form(";z (mm);%s (mm)", nm));
            for (auto ax : {fr->GetXaxis(), fr->GetYaxis()}){ax->SetTitleSize(0.055); ax->SetLabelSize(0.05);}
            fr->GetYaxis()->SetTitleOffset(1.1); fr->Draw("axis");
            for (double zc : layout.blockZ){auto b = new TBox(zc - half, lo, zc + half, hi); b->SetFillColor(kGray); b->SetLineColor(kGray); b->Draw();}
            for (double zp : layout.planeZ){
                if (zp < zlo || zp > zhi) continue;
                auto pl = new TLine(zp, lo, zp, hi); pl->SetLineColor(kGreen + 2); pl->SetLineStyle(7); pl->SetLineWidth(2); pl->Draw();
            }
            fr->Draw("axis same");
            // the tube: its true width in the zoom, a thick line in the full view (0.3 mm is under a pixel there)
            if (zoom){
                double bz[4] = {z0, z1, z1, z0}, bu[4] = {u0 + rad, u1 + rad, u1 - rad, u0 - rad};
                auto band = new TGraph(4, bz, bu); band->SetFillColor(kRed - 9); band->SetLineColor(kRed - 9); band->Draw("f same");
            }
            if (!zoom){auto tl = new TGraph(2); tl->SetPoint(0, z0, u0); tl->SetPoint(1, z1, u1); tl->SetLineColor(kRed - 9); tl->SetLineWidth(8); tl->Draw("l same");}
            auto fl = new TGraph(2); fl->SetPoint(0, zlo, ui + bi * zlo); fl->SetPoint(1, zhi, ui + bi * zhi); fl->SetLineColor(kBlue + 1); fl->SetLineStyle(2); fl->SetLineWidth(2); fl->Draw("l same");
            // true hit positions: white dots with a black ring, under the black block hits
            if (!ut.empty()){
                auto ring = new TGraph((int)ut.size(), zt.data(), ut.data()); ring->SetMarkerStyle(20); ring->SetMarkerColor(kBlack); ring->SetMarkerSize(zoom ? 2.0 : 1.3); ring->Draw("p same");
                auto dot = new TGraph((int)ut.size(), zt.data(), ut.data()); dot->SetMarkerStyle(20); dot->SetMarkerColor(kWhite); dot->SetMarkerSize(zoom ? 1.5 : 0.9); dot->Draw("p same");
            }
            if (!u.empty()){auto g = new TGraph((int)u.size(), z.data(), u.data()); g->SetMarkerStyle(20); g->SetMarkerSize(zoom ? 0.9 : 0.6); g->Draw("p same");}
            auto v = new TMarker(zv, uv, 29); v->SetMarkerColor(kRed + 1); v->SetMarkerSize(2.5); v->Draw();
            const std::string zoomTxt = zoom ? std::string(Form(", zoom on block %d", blk)) : "";
            TLatex tt; tt.SetNDC(); tt.SetTextAlign(22); tt.SetTextSize(0.048);
            tt.DrawLatex(0.55, 0.935, Form("Event %d TrackID %d : %s(z)%s", evn, trk, nm, zoomTxt.c_str()));
        }

        // the numbers of the track, right of its panels
        tp->cd();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.03);
        double y = 1.0 - 0.015 - double(k) / nEx; const double dy = 0.0165 * 2.0 / nEx;
        auto line = [&](const std::string& str, bool bold = false){tx.SetTextFont(bold ? 62 : 42); tx.DrawLatex(0.02, y, str.c_str()); y -= dy;};
        const double xlo = std::min(x0, x1) - rad, xhi = std::max(x0, x1) + rad, ylo = std::min(y0, y1) - rad, yhi = std::max(y0, y1) + rad;
        const double fz = (zv - z0) / (z1 - z0);
        line(Form("Event %d TrackID %d (%s, %d tracker hits)", evn, trk, fromCE ? "CE proton" : "not from a CE", nh), true);
        line(Form("True vertex (sim truth): x = %.3f, y = %.3f, z = %.1f mm", xv, yv, zv));
        line(Form("Tube = CE block %d, z from %.0f to %.0f mm", blk, z0, z1), true);
        line(Form("  centre at front face: x = %.3f, y = %.3f mm", x0, y0));
        line(Form("  centre at back face:  x = %.3f, y = %.3f mm", x1, y1));
        line(Form("  radius %.3f mm (#sigma = %.3f mm, x %.0f)", rad, rad / kTubeSigmas, kTubeSigmas));
        line(Form("  x range %.3f to %.3f mm,  y range %.3f to %.3f mm", xlo, xhi, ylo, yhi));
        line(Form("  true vertex distance from the tube line: %.3f mm", std::hypot(xv - (x0 + (x1 - x0) * fz), yv - (y0 + (y1 - y0) * fz))));
        line("Incoming neutron (from the beam spot, 2#sigma, to the tube)", true);
        line(Form("  #theta range %.4f#circ to %.4f#circ", thLo * rad2d, thHi * rad2d));
        if (phHalf >= M_PI - 1e-9) line("  #phi range: all (tube covers the beam axis)");
        else line(Form("  #phi range %.1f#circ #pm %.1f#circ (mid #pm half width)", phMid * rad2d, phHalf * rad2d));
        const auto nd = nDir.find({evn, par});
        if (fromCE && nd != nDir.end()){
            double th = std::atan2(std::hypot(nd->second[0], nd->second[1]), nd->second[2]), ph = std::atan2(nd->second[1], nd->second[0]); if (ph < 0) ph += 2 * M_PI;
            const bool inTh = th >= thLo && th <= thHi, inPh = phHalf >= M_PI - 1e-9 || std::fabs(std::remainder(ph - phMid, 2 * M_PI)) <= phHalf;
            line("Actual neutron direction at the CE (sim truth)", true);
            line(Form("  #theta = %.4f#circ (%s the range), #phi = %.1f#circ (%s the range)", th * rad2d, inTh ? "inside" : "OUTSIDE", ph * rad2d, inPh ? "inside" : "OUTSIDE"));
        }

        // what every line, point and shade of the panels is
        auto key = [](TObject* o, const char* label, const char* opt, TLegend* lg){lg->AddEntry(o, label, opt);};
        auto lg = new TLegend(0.02, y - 7 * dy - 0.005, 0.98, y - 0.005);
        lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetTextSize(0.028); lg->SetMargin(0.1);
        auto kBlock = new TGraph(); kBlock->SetFillColor(kGray); kBlock->SetLineColor(kGray);
        auto kPlane = new TGraph(); kPlane->SetLineColor(kGreen + 2); kPlane->SetLineStyle(7); kPlane->SetLineWidth(2);
        auto kHit = new TGraph(); kHit->SetMarkerStyle(20); kHit->SetMarkerSize(1.2);
        auto kFit = new TGraph(); kFit->SetLineColor(kBlue + 1); kFit->SetLineStyle(2); kFit->SetLineWidth(2);
        auto kBand = new TGraph(); kBand->SetFillColor(kRed - 9); kBand->SetLineColor(kRed - 9);
        auto kVertex = new TGraph(); kVertex->SetMarkerStyle(29); kVertex->SetMarkerColor(kRed + 1); kVertex->SetMarkerSize(1.8);
        key(kBlock, "carbon block", "f", lg);
        key(kPlane, "tracker plane (z of the silicon)", "l", lg);
        key(kHit, "tracker block hit (block position)", "p", lg);
        auto kTrue = new TGraph(); kTrue->SetMarkerStyle(24); kTrue->SetMarkerSize(1.2);
        key(kTrue, "tracker hit, true position (white dot)", "p", lg);
        key(kFit, "block-position fit (the tube centre)", "l", lg);
        key(kBand, Form("tube, #pm %.0f#sigma (true width in the zoom)", kTubeSigmas), "f", lg);
        key(kVertex, "true vertex (sim truth)", "p", lg);
        lg->Draw();
    }

    c.cd(0);
    TLatex lt; lt.SetNDC(); lt.SetTextSize(0.0115);
    lt.DrawLatex(0.02, 46.0 / canvasH, "Zoom rows: every x and y panel covers the same distance in mm. Full-view rows share one x range and one y range. The tube is a thick line in the full view.");
    lt.DrawLatex(0.02, 14.0 / canvasH, "Tube from the block-position fit, see ReconAngleAll.cpp. No other cuts.");
    gSystem->mkdir("tube_debug", kTRUE);
    const std::string out = event >= 0 ? "tube_debug/tube_debug_" + std::to_string(event) + (trackID >= 0 ? "_" + std::to_string(trackID) : "") + ".pdf" : "tube_debug/tube_debug.pdf";
    c.SaveAs(out.c_str());
    printf(">>> PlotTubeDebug wrote %s (python3 crop_pdfs.py tube_debug crops it)\n", out.c_str()); fflush(stdout);
}
