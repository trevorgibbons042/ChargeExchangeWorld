#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>
#include <TBox.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "../FitCore.h"
#include "../TreeColumns.h"
#include "../SeedExtend.h"
#include "../CeTracks.h"

// One event of the event fit (ReconAngleAll.cpp, file cuts/cut_angles_all.root): x(z) and y(z) side by side with the carbon blocks (gray) and the tracker planes (green dashed),
// the legend in its own column on the right, listing only what that pdf shows.
// Every charged tracker hit of the event is drawn as its track ID (TrackerHitsAll, any particle). The hits the event fit takes (SeedExtend.h, the same call as ReconAngleAll.cpp) get a black dot
// with the ID in black above it; the other hits are the ID in red. Look at the IDs of the dots to see whether the line follows one track or jumps between tracks.
// The blue dashed line is the weighted line fit of ReconAngleAll.cpp through the dots (x0DOut + bxDOut z), read from AnglesRaw, not refitted. The pink line is the vertex tube
// (tube* branches of AnglesRaw, the fit across the carbon block upstream of the first tracker hit). The black dot is the centre of the tube.
// Two pdfs, rows top to bottom:
//   tube_debug/event_fit_<event>_top.pdf   row 0: the full view of the hits and the fit alone, no tube, no tube centre, no est neutron line
//                                          (in the pdf below the tube is a thick line in row 1, its width is under a pixel there)
//   tube_debug/event_fit_<event>.pdf       row 1: the same full view with the est neutron line (black: from the beam spot x = y = z = 0 through the centre of the tube, on to the end of the row)
//                                          row 2: the view of the tube block only (20 mm before the tube to 20 mm past it), the tube at its true width; its range follows the tube, the fit line,
//                                          the neutron line and the used dots, not the 2+ hit particles
// An event without a tube has the first pdf only, named event_fit_<event>.pdf. No cuts.
// Run by hand from a run folder after ReconAngleAll.cpp; the copy steps leave debug/ out, so copy 5rootMacros/debug/ into the run folder as debug/ (the ../ includes expect the headers one folder up) or run it by its path (needs settings.csv, output0_t*.root and cuts/cut_angles_all.root).
//   root -l -b -q 'PlotTrackerDots.cpp'          the first event with 2+ hits on a plane and 3+ planes used
//   root -l -b -q 'PlotTrackerDots.cpp(387)'     event 387
void PlotTrackerDots(int event = -1){
    ROOT::EnableImplicitMT(NThreads());
    const auto s = ReadSettings();
    const StackLayout layout = StackZ(s);
    const double half = s.at("CarbonThickness") * 500.0, zdc = StackZ(s).zdcFront; // z of the ZDC front face
    const double zStart = s.count("zdcFrontZ") && s.count("blockRegion") ? zdc - s.at("blockRegion") * 1000.0 : 0.0; // the full view starts where the block region does (the gun is at z = 0, the drift before the region is empty)

    // the fit of the event: x = x0 + bx z, y = y0 + by z, and its tube
    double x0, bx, y0, by, tz0 = 0, tz1 = 0, tx0 = 0, tx1 = 0, ty0 = 0, ty1 = 0, trad = 0; int nMax, nUsed, ev, tblk = -1;
    {
        TFile fa("cuts/cut_angles_all.root");
        TTree* t = (TTree*)fa.Get("AnglesRaw");
        if (!t){Error("PlotTrackerDots", "no AnglesRaw in cuts/cut_angles_all.root (run ReconAngleAll.cpp)"); return;}
        t->SetBranchAddress("fEvent", &ev); t->SetBranchAddress("nMaxChargedPerPlaneEv", &nMax); t->SetBranchAddress("nPlanesUsed", &nUsed);
        t->SetBranchAddress("x0DOut", &x0); t->SetBranchAddress("bxDOut", &bx); t->SetBranchAddress("y0DOut", &y0); t->SetBranchAddress("byDOut", &by);
        t->SetBranchAddress("tubeBlock", &tblk); t->SetBranchAddress("tubeZ0", &tz0); t->SetBranchAddress("tubeZ1", &tz1);
        t->SetBranchAddress("tubeX0", &tx0); t->SetBranchAddress("tubeX1", &tx1); t->SetBranchAddress("tubeY0", &ty0); t->SetBranchAddress("tubeY1", &ty1); t->SetBranchAddress("tubeRadius", &trad);
        bool found = false;
        for (Long64_t i = 0; i < t->GetEntries() && !found; i++){
            t->GetEntry(i);
            found = event >= 0 ? ev == event : (nMax >= 2 && nUsed >= 3);
        }
        if (!found){Error("PlotTrackerDots", "no event found in AnglesRaw for these arguments"); return;}
        event = ev;
    }
    const bool hasTube = trad > 0;
    const int nParts = hasTube ? 2 : 1;

    // every charged tracker hit of the event with its track ID, and which of them the event fit uses
    struct Hit {double z, x, y; long id; int trk; bool ce;}; // trk: the track SeedExtend gives the hit, 0 the one the event fit uses, -1 no track
    std::vector<Hit> hits;
    {
        TChain th("TrackerHitsAll");
        th.Add("output0_t*.root");
        const auto c = Columns(&th, {"fDx", "fDy", "fDz", "TrackID", "ParentID"}, "fEvent==" + std::to_string(event) + " && fZ < fDz"); // the hits that entered through the front face, as the finder sees them
        std::vector<double> hx, hy, hz; std::vector<long> hid, hpar;
        std::map<int, std::set<int>> ceTracks; // a CE hit is a hit of a CE child proton: its ParentID is a CE neutron of the event (CeTracks.h)
        CarbonCeTracks(layout.blockZ, half, ceTracks, "PlotTrackerDots");
        const auto ceOfEvent = ceTracks.find(event);
        for (size_t i = 0; i < c[0].size(); i++) if (c[2][i] != 0.0){hx.push_back(c[0][i]); hy.push_back(c[1][i]); hz.push_back(c[2][i]); hid.push_back(std::lround(c[3][i])); hpar.push_back(std::lround(c[4][i]));}
        const double hitSigmaMm = 400.0 / s.at("nCols_blocks") / std::sqrt(12.0);
        const int trackerSkip = (s.count("trackerEfficiency") && s.at("trackerEfficiency") < 1.0) ? 1 : 0; // trackers that miss hits (sim /detector/trackerEfficiency below 1): a track may pass one plane with no hit
        // more tracks: SeedExtend again on the hits no track has taken yet, until it finds none (at most 8)
        std::vector<int> trk(hz.size(), -1);
        std::vector<long> planesAll; // every plane the event has a hit on, the hits of the tracks taken off too
        for (double zz : hz) planesAll.push_back(std::lround(zz * 1000.0));
        std::sort(planesAll.begin(), planesAll.end()); planesAll.erase(std::unique(planesAll.begin(), planesAll.end()), planesAll.end());
        for (int k = 0; k < 8; k++){
            std::vector<double> rx, ry, rz; std::vector<size_t> at;
            for (size_t i = 0; i < hz.size(); i++) if (trk[i] < 0){rx.push_back(hx[i]); ry.push_back(hy[i]); rz.push_back(hz[i]); at.push_back(i);}
            const auto found = SeedExtend(rx, ry, rz, planesAll, hitSigmaMm, trackerSkip);
            if (found.empty()) break;
            for (size_t j : found) trk[at[j]] = k;
        }
        for (size_t i = 0; i < hz.size(); i++) hits.push_back({hz[i], hx[i], hy[i], hid[i], trk[i], ceOfEvent != ceTracks.end() && ceOfEvent->second.count((int)hpar[i])});
    }
    if (hits.empty()){Error("PlotTrackerDots", "no tracker hits for event %d in output0_t*.root", event); return;}
    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b){return a.z < b.z;});

    // where the line of each track runs: back to the downstream face of the nearest carbon block upstream of its first tracker hit (the block it could have come out of), on to its last tracker hit
    double zBack[8], zLast[8], zFirst[8];
    for (int k = 0; k < 8; k++){zFirst[k] = INFINITY; zLast[k] = -INFINITY;}
    for (const auto& h : hits) if (h.trk >= 0){zFirst[h.trk] = std::min(zFirst[h.trk], h.z); zLast[h.trk] = std::max(zLast[h.trk], h.z);}
    for (int k = 0; k < 8; k++){
        zBack[k] = zFirst[k];
        double zb = -INFINITY;
        for (double zc : layout.blockZ) if (zc < zFirst[k]) zb = std::max(zb, zc);
        if (std::isfinite(zb)) zBack[k] = zb + half;
    }

    // the z range of each row, and for the zoom row (2) the u range: one span for x and y, each centred on its own data
    // row 2 is the tube block only: from just before the tube to just past it, so the tube at its true width is visible
    double zlo[3] = {zStart, zStart, tz0 - 20}, zhi[3] = {zdc, zdc, tz1 + 20}, mid[3][2] = {{0, 0}, {0, 0}, {0, 0}}, span[3] = {0, 0, 0};
    constexpr int row = 2; // only the zoom row has a span of its own
    if (hasTube) for (int a = 0; a < 2; a++){
        const double u0 = a ? ty0 : tx0, u1 = a ? ty1 : tx1, ui = a ? y0 : x0, bi = a ? by : bx;
        double lo = std::min(u0, u1) - trad, hi = std::max(u0, u1) + trad;
        const double nSlope = (u0 + u1) / (tz0 + tz1); // the est neutron line, 0 at z = 0 through the tube centre
        for (double zz : {zlo[row], zhi[row]}){lo = std::min({lo, ui + bi * zz, nSlope * zz}); hi = std::max({hi, ui + bi * zz, nSlope * zz});}
        for (const auto& h : hits) if (h.trk == 0 && h.z >= zlo[row] && h.z <= zhi[row]){const double u = a ? h.y : h.x; lo = std::min(lo, u); hi = std::max(hi, u);} // the used dots only
        mid[row][a] = 0.5 * (lo + hi); span[row] = std::max(span[row], 1.25 * (hi - lo));
    }

    gStyle->SetOptStat(0);
    // the colour and marker shape of each track: track 0 (blue, like its fit line) is the one the event fit uses
    const int kTrkMark[8] = {20, 21, 22, 23, 33, 34, 29, 47};
    const int kTrkCol[8] = {kBlue + 1, kMagenta + 1, kCyan + 2, kOrange + 7, kViolet + 1, kGreen + 3, kAzure + 1, kPink + 1};
    gStyle->SetLineWidth(1); gStyle->SetFrameLineWidth(1);
    gSystem->mkdir("tube_debug", kTRUE);
    for (int part = 0; part < nParts; part++){
        // part 0 draws row 0 alone, part 1 rows 1 and 2
        const int nR = part == 0 ? 1 : 2;
        const int rowOf[2] = {part == 0 ? 0 : 1, 2};
        TCanvas cv(Form("cv%d", part), "", 2600, 650 * nR);
        auto main = new TPad("m", "", 0, 0, 0.62, 1); main->Draw(); main->SetFillStyle(0); main->Divide(2, nR);
        auto lp = new TPad("lp", "", 0.62, 0, 1, 1); lp->Draw(); lp->SetFillStyle(0);
        bool anyOther = false, anyCE = false;
        for (int k = 0; k < nR; k++) for (int a = 0; a < 2; a++){
            const int row = rowOf[k];
            main->cd(2 * k + a + 1);
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.12);
            const char* nm = a ? "y" : "x";
            const bool zoom = row == 2; // the zoom row; rows 0 and 1 are the full view
            double lo, hi;
            if (zoom){lo = mid[row][a] - span[row] / 2; hi = mid[row][a] + span[row] / 2;}
            else {
                double m = 0;
                for (const auto& h : hits) m = std::max(m, std::fabs(a ? h.y : h.x));
                hi = 1.2 * m + 0.2; lo = -hi;
            }
            auto fr = new TH1F("", "", 10, zlo[row], zhi[row]); fr->SetMinimum(lo); fr->SetMaximum(hi);
            fr->SetTitle(Form(";z (mm);%s (mm)", nm));
            for (auto ax : {fr->GetXaxis(), fr->GetYaxis()}){ax->SetTitleSize(0.055); ax->SetLabelSize(0.05); ax->SetAxisColor(kBlack); ax->SetTickLength(0.02);}
            if (zoom) fr->GetXaxis()->SetNdivisions(505);
            fr->GetYaxis()->SetTitleOffset(1.2); fr->Draw("axis");
            for (double zc : layout.blockZ){
                if (zc + half < zlo[row] || zc - half > zhi[row]) continue; // a box fully outside the z range is not clipped right in the pdf: it shows as a gray bar inside the panel
                auto b = new TBox(zc - half, lo, zc + half, hi); b->SetFillColor(kGray); b->SetLineColor(kGray); b->SetLineWidth(0); b->Draw();}
            for (double zp : layout.planeZ){
                if (zp < zlo[row] || zp > zhi[row]) continue;
                auto pl = new TLine(zp, lo, zp, hi); pl->SetLineColor(kGreen + 2); pl->SetLineStyle(7); pl->SetLineWidth(1); pl->Draw();
            }
            fr->Draw("axis same");
            // the vertex tube, under the fit line and the dots: its true width in the zoom row, a thick line in the full view
            if (hasTube && part == 1){
                const double u0 = a ? ty0 : tx0, u1 = a ? ty1 : tx1;
                if (zoom){
                    double bz[4] = {tz0, tz1, tz1, tz0}, bu[4] = {u0 + trad, u1 + trad, u1 - trad, u0 - trad};
                    auto band = new TGraph(4, bz, bu); band->SetFillColor(kRed - 9); band->SetLineColor(kRed - 9); band->Draw("f same");
                } else {auto tl = new TGraph(2); tl->SetPoint(0, tz0, u0); tl->SetPoint(1, tz1, u1); tl->SetLineColor(kRed - 9); tl->SetLineWidth(6); tl->Draw("l same");}
            }
            // the fit line (track 0) runs from its carbon block to its last tracker hit, see zBack above
            const double zA = std::max(zlo[row], zBack[0]), zB = std::min(zhi[row], zLast[0]);
            // the fit line, thin and dashed
            const double ui = a ? y0 : x0, bi = a ? by : bx;
            auto fl = new TGraph(2); fl->SetPoint(0, zA, ui + bi * zA); fl->SetPoint(1, zB, ui + bi * zB);
            fl->SetLineColorAlpha(kBlue + 1, 0.6); fl->SetLineStyle(1); fl->SetLineWidth(1); fl->Draw("l same");
            // the est neutron line (not in row 0): from the beam spot (0, 0) through the tube centre, on to the end of the row
            if (hasTube && row >= 1){
                const double nSlope = ((a ? ty0 + ty1 : tx0 + tx1)) / (tz0 + tz1), zs = std::max(zlo[row], 0.0);
                auto nl = new TGraph(2); nl->SetPoint(0, zs, nSlope * zs); nl->SetPoint(1, zhi[row], nSlope * zhi[row]);
                nl->SetLineColor(kBlack); nl->SetLineWidth(1); nl->Draw("l same");
            }
            // the black dot: the centre of the tube
            if (hasTube && part == 1){
                auto cd = new TGraph(1); cd->SetPoint(0, 0.5 * (tz0 + tz1), 0.5 * (a ? ty0 + ty1 : tx0 + tx1));
                cd->SetMarkerStyle(20); cd->SetMarkerColor(kBlack); cd->SetMarkerSize(zoom ? 1.4 : 0.9); cd->Draw("p same");
            }
            // the hits: a hit of a track is a dot with its track ID in the colour of that track, any other hit is its ID in red. Each track other than 0 is joined by a thin line
            // through its hits (track 0 has the fit line).
            std::vector<std::vector<double>> tz(8), tu(8);
            std::vector<double> oz, ou, cez, ceu; // the hits of no track, and the CE hits
            TLatex idTxt; idTxt.SetTextSize(zoom ? 0.032 : 0.026);
            for (const auto& h : hits){
                if (h.z < zlo[row] || h.z > zhi[row]) continue;
                const double u = a ? h.y : h.x;
                if (u < lo || u > hi) continue;
                if (h.ce){cez.push_back(h.z); ceu.push_back(u);}
                if (h.trk >= 0){tz[h.trk].push_back(h.z); tu[h.trk].push_back(u); idTxt.SetTextColor(kTrkCol[h.trk]);}
                else {anyOther = true; oz.push_back(h.z); ou.push_back(u); idTxt.SetTextColor(kRed + 1);}
                std::string label = std::to_string(h.id), tight; // the digits of an ID pulled together, a "1" in a wide digit cell reads as a gap ("11" looked like "1 1")
                for (size_t d = 0; d < label.size(); d++) tight += d == 0 ? std::string(1, label[d]) : "#kern[-0.12]{" + std::string(1, label[d]) + "}";
                idTxt.SetTextAlign(11); idTxt.DrawLatex(h.z + 0.01 * (zhi[row] - zlo[row]), u + 0.015 * (hi - lo), tight.c_str()); // up and to the right of the dot, clear of it
            }
            if (!cez.empty()){anyCE = true; auto cg = new TGraph((int)cez.size(), cez.data(), ceu.data()); cg->SetMarkerStyle(24); cg->SetMarkerColor(CbCE()); cg->SetMarkerSize(zoom ? 2.2 : 1.6); cg->Draw("p same");}
            if (!oz.empty()){auto og = new TGraph((int)oz.size(), oz.data(), ou.data()); og->SetMarkerStyle(20); og->SetMarkerColor(kRed + 1); og->SetMarkerSize(zoom ? 1.0 : 0.7); og->Draw("p same");}
            for (int k = 0; k < 8; k++) if (!tz[k].empty()){
                auto g = new TGraph((int)tz[k].size(), tz[k].data(), tu[k].data());
                g->SetMarkerStyle(kTrkMark[k]); g->SetMarkerColor(kTrkCol[k]); g->SetMarkerSize(zoom ? 1.2 : 0.9);
                g->Draw("p same");
                // the line of a track other than 0: the least squares line through all its hits, dashed in its colour (track 0 has the fit line of AnglesRaw)
                if (k == 0) continue;
                double n = 0, sz = 0, su = 0, szz = 0, szu = 0;
                for (const auto& h : hits) if (h.trk == k){const double u = a ? h.y : h.x; n++; sz += h.z; su += u; szz += h.z * h.z; szu += h.z * u;}
                const double den = n * szz - sz * sz;
                if (n < 2 || den == 0) continue;
                const double bk = (n * szu - sz * su) / den, ak = (su - bk * sz) / n;
                auto ln = new TGraph(2); const double kA = std::max(zlo[row], zBack[k]), kB = std::min(zhi[row], zLast[k]); // this track's own ends
                ln->SetPoint(0, kA, ak + bk * kA); ln->SetPoint(1, kB, ak + bk * kB);
                ln->SetLineColorAlpha(kTrkCol[k], 0.6); ln->SetLineStyle(1); ln->SetLineWidth(1); ln->Draw("l same");
            }
            TLatex tt; tt.SetNDC(); tt.SetTextAlign(22); tt.SetTextSize(0.05);
            if (zoom) tt.DrawLatex(0.55, 0.94, Form("Event %d : %s(z), zoom on the tube block, no cuts", event, nm));
            else tt.DrawLatex(0.55, 0.94, Form("Event %d : %s(z), track IDs of the hits, no cuts", event, nm));
        }

        // the legend column: a key for each thing drawn in this pdf, empty graphs styled like the drawn objects (the white dot of a used hit is shown as an open ring, white on white would not show)
        lp->cd();
        auto kBlock = new TGraph(); kBlock->SetFillColor(kGray); kBlock->SetLineColor(kGray);
        auto kPlane = new TGraph(); kPlane->SetLineColor(kGreen + 2); kPlane->SetLineStyle(7); kPlane->SetLineWidth(1);
        auto kUsed = new TGraph(); kUsed->SetMarkerStyle(20); kUsed->SetMarkerSize(0.9);
        auto kIgn = new TGraph(); kIgn->SetMarkerStyle(20); kIgn->SetMarkerColor(kRed + 1); kIgn->SetMarkerSize(0.6);
        auto kFit = new TGraph(); kFit->SetLineColor(kBlue + 1); kFit->SetLineStyle(1); kFit->SetLineWidth(1);
        auto kTube = new TGraph(); kTube->SetFillColor(kRed - 9); kTube->SetLineColor(kRed - 9); kTube->SetLineWidth(6);
        auto kNeu = new TGraph(); kNeu->SetLineColor(kBlack); kNeu->SetLineWidth(1);
        auto kCen = new TGraph(); kCen->SetMarkerStyle(20); kCen->SetMarkerColor(kBlack); kCen->SetMarkerSize(0.9);
        // TLegend shrinks all its text to fit the longest label in the width and spreads the entries over the box height: short labels, and a box sized from the entries,
        // so the text is the same size in pixels whatever the number of rows of the canvas
        std::vector<std::pair<TGraph*, std::pair<const char*, const char*>>> keys;
        keys.push_back({kBlock, {"carbon block", "f"}});
        keys.push_back({kPlane, {"tracker position", "l"}});
        kUsed->SetMarkerStyle(kTrkMark[0]); kUsed->SetMarkerColor(kTrkCol[0]); keys.push_back({kUsed, {"track 0, in the fit (ID)", "p"}});
        int nTrk = 0; for (const auto& h : hits) nTrk = std::max(nTrk, h.trk + 1);
        for (int k = 1; k < nTrk; k++){auto g = new TGraph(); g->SetMarkerStyle(kTrkMark[k]); g->SetMarkerSize(0.9); g->SetMarkerColor(kTrkCol[k]); g->SetLineColor(kTrkCol[k]); g->SetLineStyle(1); keys.push_back({g, {Form("track %d, line (ID)", k), "lp"}});}
        if (anyCE){auto kc = new TGraph(); kc->SetMarkerStyle(24); kc->SetMarkerColor(CbCE()); kc->SetMarkerSize(1.0); keys.push_back({kc, {"CE hit (ring)", "p"}});}
        if (anyOther) keys.push_back({kIgn, {"other hit (track ID, red)", "p"}});
        keys.push_back({kFit, {"weighted line fit", "l"}});
        const char* tubeLabel = Form("vertex tube (block %d, %.1f mm)", tblk, trad);
        if (hasTube && part == 1) keys.push_back({kTube, {tubeLabel, "f"}});
        if (hasTube && part == 1) keys.push_back({kNeu, {"est neutron direction", "l"}}); // from (0, 0) through the tube centre
        if (hasTube && part == 1) keys.push_back({kCen, {"tube centre", "p"}});
        const double boxH = keys.size() * 52.0 / (650.0 * nR);
        auto lg = new TLegend(0.02, 0.5 - boxH / 2, 0.98, 0.5 + boxH / 2);
        lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetTextSize(0.0723 / nR); lg->SetMargin(0.12);
        for (const auto& k : keys) lg->AddEntry(k.first, k.second.first, k.second.second);
        lg->Draw();

        cv.SaveAs(Form("tube_debug/event_fit_%d%s.pdf", event, hasTube && part == 0 ? "_top" : ""));
    }
}
