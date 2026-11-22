#include <TCanvas.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include "FitCore.h"
#include "TreeColumns.h"
#include "CutFiles.h"

// The numbers are the Events tree of cuts/cut_events.root (ReconAngleAll.cpp runs the finder below on TrackerHitsAll once and writes the per event tally); this macro does not read the sim output any more.
// plots/8_trackerfits: how well the seed and extend track finding (SeedExtend.h) rebuilds the tracks of an event, from TrackerHitsAll (every charged track that entered a tracker), no cuts.
// Every event with a tracker hit: SeedExtend runs on the tracker hits (block positions, no track ID) and its hits are taken off, then it runs again on the rest, until it finds no track.
// A reconstructed track is pure when all its hits carry one track ID, mixed when they carry 2 or more. A true track (a track ID with 2 or more tracker hits) is missing when none of its hits is in any
// reconstructed track of the event, pure or mixed. A CE hit is a hit of a CE child proton (its ParentID is a CE neutron of the event, CeTracks.h), a CE track a track with such hits only.
//   4_reconstructable_tracks.pdf  every true track with 2 or more tracker hits, what the finding could rebuild: all particles, and the CE ones (CE child protons)
// One pdf per kind, the x axis the number of tracks of that kind in an event, the y axis the number of events, both axes the same in the four pdfs:
//   1_pure_tracks.pdf     pure reconstructed tracks: all, and the CE ones (the CE proton rebuilt on its own)
//   2_mixed_tracks.pdf    mixed reconstructed tracks: the ones with no CE hit in them (green), and the CE particles that have no pure track and sit in mixed tracks (blue)
//   3_missing_tracks.pdf  true tracks with 2 or more hits that no reconstructed track (pure or mixed) holds: all, and the CE ones
//   5_max_hits_per_track.pdf  per event the most tracker hits of any one charged particle (the longest track): all events, CE (a CE particle has a pure reconstructed track) and CE mix (its hits only in mixed reconstructed tracks), the same CE and CE mix as 1_pure_tracks and 2_mixed_tracks
//   6_last_tracker_of_track.pdf  for every true track with 2 or more hits, the last (highest z) tracker it hits, 1 = the first tracker: all particles, and the CE child protons
//   11_max_hits_clean_chain.pdf  the events of 5_max_hits_per_track by chain length: the longest run of trackers in a row with exactly 1 charged hit each, any start (2,1,1,1 gives 3)
//   9_maxhits_vs_last_tracker_CE.pdf, 10_maxhits_vs_last_tracker_all.pdf  2D, per event the tracker hits of the longest track (5_max_hits_per_track) against the last tracker of that track (6_last_tracker_of_track): CE events (a CE particle has a pure track), and all events with a track of 2 or more hits
//   12_chain_vs_last_tracker_CE.pdf, 13_chain_vs_last_tracker_all.pdf  2D, the clean chain length of 11_max_hits_clean_chain against the last tracker of the longest track, CE and all events
//   1_indivtrackers/  one pdf per tracker plane (named by its order and z): the charged hits on that plane per event, all, CE (an event whose CE particle has a pure reconstructed track) and CE + CE mix (CE, and the events whose CE particle is only in mixed reconstructed tracks), only the events with at least 1 hit on the plane (no 0 bin)
void F8_TrackerFits(){
    FastJit();
    ROOT::EnableImplicitMT(NThreads());
    gStyle->SetOptStat(0);

    // The tally of every event is a row of the Events tree (cuts/cut_events.root, ReconAngleAll.cpp runs the finder on the raw hits once and writes it): this macro only draws it.
    TTree* ev = EventsTree();
    if (!ev){Error("F8_TrackerFits", "no Events tree in cuts/cut_events.root (run ReconAngleAll.cpp)"); return;}
    const auto settings = ReadSettings();
    // the tracker planes of the run (z in um), the plane before the ZDC last: planeAll all of them; below `planes` keeps the planes some event hits, as F8 always did, and the plane numbers of ReconAngleAll are renumbered to it
    std::vector<long> planeAll;
    for (double z : StackZ(settings).planeZ) planeAll.push_back(std::lround(z * 1000.0));
    const size_t nAll = planeAll.size();
    if (ev->GetEntries() == 0 || nAll == 0){Error("F8_TrackerFits", "no events or no tracker planes"); return;}

    long long nTrueTracks = 0; // the true tracks with 2 or more tracker hits, over all events
    std::vector<int> nPure, nPureCE, nMixedCE, nMixedNoCE, nMissing, nMissingCE, nTrue, nTrueCE; // per event, as the Events columns nPureEv, nPureCEEv, nMixedCEEv, nMixedNoCEEv, nMissingEv, nMissingCEEv, nTrueEv, nTrueCEEv
    std::vector<std::vector<int>> nHitAll(nAll); // per plane of the run, the front face tracker hits of each event (nHitPlane<i>)
    std::vector<int> lastPlane, lastPlaneCE, lastPlaneMix; // per true track with 2+ hits (all events), the number of its last (highest z) tracker plane, 1 = the first tracker: all tracks, and the tracks of the events of class CE (a CE particle with a pure track) and CE mix (a CE child proton, pure track or only inside mixed tracks)
    long long nEvLast = 0, nEvLastCE = 0, nEvLastMix = 0; // the events that gave at least 1 track to lastPlane, lastPlaneCE, lastPlaneMix (the events of the 6_last_tracker_of_track legend)
    std::vector<int> chainLen, maxHitLast, maxHit, evClass; // per event: the longest run of trackers in a row with exactly 1 charged hit; the last tracker of its longest track; the most tracker hits of any one track ID; its class (trkClassEv): 0 no CE particle with 2+ hits, 1 CE, 2 CE mix
    {
        // the scalar columns of every event in entry order, one event loop on all the threads (Columns, TreeColumns.h), the hits of every plane of the run with them
        std::vector<std::string> names = {"trkClassEv", "nPureEv", "nPureCEEv", "nMixedCEEv", "nMixedNoCEEv", "nMissingEv", "nMissingCEEv", "nTrueEv", "nTrueCEEv", "maxHitEv", "maxHitLastEv", "chainLenEv"};
        for (size_t p = 0; p < nAll; p++) names.push_back(Form("nHitPlane%zu", p));
        const auto cols = Columns(ev, names);
        const auto asInts = [](const std::vector<double>& v){std::vector<int> r(v.size()); for (size_t i = 0; i < v.size(); i++) r[i] = (int)std::lround(v[i]); return r;};
        evClass = asInts(cols[0]); nPure = asInts(cols[1]); nPureCE = asInts(cols[2]); nMixedCE = asInts(cols[3]); nMixedNoCE = asInts(cols[4]); nMissing = asInts(cols[5]); nMissingCE = asInts(cols[6]);
        nTrue = asInts(cols[7]); nTrueCE = asInts(cols[8]); maxHit = asInts(cols[9]); maxHitLast = asInts(cols[10]); chainLen = asInts(cols[11]);
        for (size_t p = 0; p < nAll; p++) nHitAll[p] = asInts(cols[12 + p]);
        for (int t : nTrue) nTrueTracks += t;

        // the true tracks of every event (trueLastPlane, trueKind: one entry per true track with 2+ hits) only go into histograms, so the order of the events does not matter: each thread fills its own tally, the tallies are added up
        // kind 0 not CE, 1 CE with no pure track, 2 CE with a pure track
        struct LastTally {std::vector<int> all, ce, mix; long long nEv = 0, nEvCE = 0, nEvMix = 0;};
        ROOT::RDataFrame dfEv(*ev);
        std::vector<LastTally> tally(dfEv.GetNSlots());
        dfEv.ForeachSlot([&](unsigned slot, const ROOT::RVec<int>& trueLast, const ROOT::RVec<int>& trueKind){
            LastTally& t = tally[slot];
            bool anyCE = false, anyMix = false;
            for (size_t k = 0; k < trueLast.size(); k++){
                t.all.push_back(trueLast[k]);
                if (trueKind[k] >= 1){t.mix.push_back(trueLast[k]); anyMix = true;}
                if (trueKind[k] == 2){t.ce.push_back(trueLast[k]); anyCE = true;}
            }
            if (!trueLast.empty()){t.nEv++; if (anyCE) t.nEvCE++; if (anyMix) t.nEvMix++;} // mix = CE + CE mix, so nEvLastMix >= nEvLastCE
        }, {"trueLastPlane", "trueKind"});
        for (const LastTally& t : tally){
            lastPlane.insert(lastPlane.end(), t.all.begin(), t.all.end()); lastPlaneCE.insert(lastPlaneCE.end(), t.ce.begin(), t.ce.end()); lastPlaneMix.insert(lastPlaneMix.end(), t.mix.begin(), t.mix.end());
            nEvLast += t.nEv; nEvLastCE += t.nEvCE; nEvLastMix += t.nEvMix;
        }
    }
    // the planes some event hits, in z order: the numbers of the last tracker (1 = first) and the chain are over these, as before the tally moved into ReconAngleAll.cpp
    std::vector<long> planes; std::vector<int> planeNo(nAll + 1, 0); // planeNo[old number] = the new number, 0 for a plane no event hits
    std::vector<std::vector<int>> nHit;
    for (size_t p = 0; p < nAll; p++) if (std::any_of(nHitAll[p].begin(), nHitAll[p].end(), [](int n){return n > 0;})){
        planes.push_back(planeAll[p]); planeNo[p + 1] = (int)planes.size(); nHit.push_back(std::move(nHitAll[p]));
    }
    const size_t nPl = planes.size();
    if (nPl == 0){Error("F8_TrackerFits", "no tracker hits"); return;}
    for (auto* v : {&lastPlane, &lastPlaneCE, &lastPlaneMix, &maxHitLast}) for (int& n : *v) if (n > 0) n = planeNo[n];
    for (size_t e = 0; e < chainLen.size(); e++){ // the longest run of the planes some event hits in a row with exactly 1 charged hit
        int run = 0, best = 0;
        for (size_t p = 0; p < nPl; p++){run = nHit[p][e] == 1 ? run + 1 : 0; best = std::max(best, run);}
        chainLen[e] = best;
    }

    gSystem->mkdir("plots/8_trackerfits", true);
    std::set<std::string> numberFiles;
    const long long nEvents = (long long)nPure.size();
    // a plot: two counts per event drawn over each other, in the colours of F9_PlotReconNew.cpp (ClassColor): all and no CE green, CE vermillion, CE mix (a CE hit among others) blue
    struct Kind {const char* file; const char* what; const char* label1; const std::vector<int>* v1; int col1; const char* label2; const std::vector<int>* v2; int col2; const char* xlab;};
    const Kind kinds[4] = {
        {"1_pure_tracks",    "Pure reconstructed tracks (all hits one track ID)", "all", &nPure, CbAll(), "CE particles", &nPureCE, CbCE(), "Reconstructed tracks per event"},
        {"2_mixed_tracks",   "Mixed reconstructed tracks (hits of 2+ track IDs)", "no CE hit", &nMixedNoCE, CbAll(), "CE particles, in mixed tracks only", &nMixedCE, CbProton(), "Reconstructed tracks per event"},
        {"3_missing_tracks", "Missing tracks (2+ hits, no hit in any reconstructed track)", "all", &nMissing, CbAll(), "CE", &nMissingCE, CbCE(), "Missing tracks per event (true, not reconstructed)"},
        {"4_reconstructable_tracks", "Reconstructable tracks (true tracks with 2+ hits)", "all particles", &nTrue, CbAll(), "CE particles", &nTrueCE, CbCE(), "Reconstructable tracks per event"},
    };
    // the three plots share both axes: the x range reaches the highest count of any of them, the y range the tallest bin of any of them
    int top = 1;
    for (const Kind& k : kinds) top = std::max({top, *std::max_element(k.v1->begin(), k.v1->end()), *std::max_element(k.v2->begin(), k.v2->end())});
    auto fill = [&](const char* name, const std::vector<int>& v){
        TH1D* h = new TH1D(name, "", top, 0.5, top + 0.5); // events with no track of the kind are not drawn: their bin would be the whole sample
        for (int n : v) if (n >= 1) h->Fill(n);
        return h;
    };
    std::vector<TH1D*> h1, h2;
    double tallest = 1.0;
    for (const Kind& k : kinds){
        h1.push_back(fill(Form("h1_%s", k.file), *k.v1)); h2.push_back(fill(Form("h2_%s", k.file), *k.v2));
        tallest = std::max({tallest, h1.back()->GetMaximum(), h2.back()->GetMaximum()});
    }
    TCanvas canvas("c", "", 1200, 900);
    canvas.SetLogy(1);
    for (size_t i = 0; i < 4; i++){
        const Kind& k = kinds[i];
        TH1D *a = h1[i], *b = h2[i];
        const long long n1 = std::accumulate(k.v1->begin(), k.v1->end(), 0LL), n2 = std::accumulate(k.v2->begin(), k.v2->end(), 0LL); // tracks over all events
        a->SetLineColor(k.col1); a->SetLineWidth(2);
        b->SetLineColor(k.col2); b->SetLineWidth(2);
        const std::string path = std::string("plots/8_trackerfits/") + k.file + ".pdf";
        a->SetTitle(PlotTitle(path, Form("%s per event, %lld true tracks with 2+ hits in %lld events, no cuts;%s;Events", k.what, nTrueTracks, nEvLast, k.xlab)).c_str());
        a->SetMinimum(0.5); a->SetMaximum(3.0 * tallest);
        canvas.cd();
        a->Draw("hist"); b->Draw("hist same");
        TLegend lg(0.55, 0.74, 0.88, 0.88);
        lg.SetBorderSize(0); lg.SetFillStyle(0);
        lg.AddEntry(a, Form("%s: %lld tracks, %.0f events", k.label1, n1, a->GetEntries()), "l");
        lg.AddEntry(b, Form("%s: %lld tracks, %.0f events", k.label2, n2, b->GetEntries()), "l");
        lg.Draw();
        canvas.Update();
        WritePlotNumbers(a, path, numberFiles);
        canvas.SaveAs(path.c_str());
    }
    // plots/8_trackerfits/1_indivtrackers: one pdf per tracker plane, the charged hits on the plane per event (all, and the CE child protons), only the events with at least 1 hit on that plane (n >= 1, no 0 bin).
    // Both axes are the same in all the planes.
    gSystem->mkdir("plots/8_trackerfits/1_indivtrackers", true);
    int topHits = 1;
    for (size_t p = 0; p < nPl; p++) topHits = std::max({topHits, *std::max_element(nHit[p].begin(), nHit[p].end())});
    std::vector<TH1D*> hAll, hCE, hMix; // hits on the plane per event: all events, events whose hits on the plane are all CE child protons, events with a CE child proton and other hits there
    double tallestHits = 1.0;
    for (size_t p = 0; p < nPl; p++){
        hAll.push_back(new TH1D(Form("hitAll_%zu", p), "", topHits, 0.5, topHits + 0.5));
        hCE.push_back(new TH1D(Form("hitCE_%zu", p), "", topHits, 0.5, topHits + 0.5));
        hMix.push_back(new TH1D(Form("hitMix_%zu", p), "", topHits, 0.5, topHits + 0.5));
        for (size_t e = 0; e < nHit[p].size(); e++){
            const int n = nHit[p][e];
            if (n >= 1) hAll.back()->Fill(n);
            if (n >= 1 && evClass[e] == 1) hCE.back()->Fill(n);
            if (n >= 1 && evClass[e] >= 1) hMix.back()->Fill(n); // CE + CE mix, every event with a CE particle
        }
        tallestHits = std::max(tallestHits, hAll.back()->GetMaximum());
    }
    for (size_t p = 0; p < nPl; p++){
        const double zMm = planes[p] / 1000.0;
        TH1D *a = hAll[p], *b = hCE[p], *m = hMix[p];
        const long long n1 = std::accumulate(nHit[p].begin(), nHit[p].end(), 0LL); // hits over all events
        a->SetLineColor(CbAll()); a->SetLineWidth(2);
        b->SetLineColor(CbCE()); b->SetLineWidth(2);
        m->SetLineColor(CbProton()); m->SetLineWidth(2);
        const std::string path = Form("plots/8_trackerfits/1_indivtrackers/%zu_tracker_z%ldmm.pdf", p + 1, std::lround(zMm));
        const long long nOnPlane = std::count_if(nHit[p].begin(), nHit[p].end(), [](int n){return n >= 1;}); // events with a charged hit on this plane
        a->SetTitle(PlotTitle(path, Form("Charged hits on tracker %zu at z = %.1f m per event, %lld events with a hit on it, no cuts;Charged hits on the tracker per event;Events", p + 1, zMm / 1000.0, nOnPlane)).c_str());
        a->SetMinimum(0.5); a->SetMaximum(3.0 * tallestHits);
        canvas.cd();
        a->Draw("hist"); m->Draw("hist same"); b->Draw("hist same"); // CE (vermillion) drawn last, on top of CE + CE mix (blue)
        TLegend lg(0.5, 0.70, 0.88, 0.88);
        lg.SetBorderSize(0); lg.SetFillStyle(0);
        lg.AddEntry(a, Form("all charged: %lld hits, mean %.2f per event", n1, a->GetMean()), "l");
        lg.AddEntry(b, Form("CE: %.0f events", b->GetEntries()), "l");
        lg.AddEntry(m, Form("CE + CE mix: %.0f events", m->GetEntries()), "l");
        lg.Draw();
        canvas.Update();
        WritePlotNumbers(a, path, numberFiles);
        canvas.SaveAs(path.c_str());
    }
    // 5_max_hits_per_track.pdf: per event the most tracker hits of any one charged particle (the longest track), all events, the CE events (a CE particle has a pure reconstructed track) and the CE + CE mix events together (a CE particle with a pure reconstructed track, or only inside mixed tracks)
    {
        const int topMax = std::max(1, *std::max_element(maxHit.begin(), maxHit.end()));
        TH1D *a = new TH1D("maxHitAll", "", topMax, 0.5, topMax + 0.5), *b = new TH1D("maxHitCE", "", topMax, 0.5, topMax + 0.5), *m = new TH1D("maxHitCEMix", "", topMax, 0.5, topMax + 0.5); // b: CE, m: CE + CE mix
        long long tracksAll = 0, tracksCE = 0, tracksCEMix = 0; // the true tracks with 2+ hits of the events in each line
        for (size_t e = 0; e < maxHit.size(); e++){
            if (maxHit[e] < 2) continue; // only the events with a track of 2+ hits, the events of 6_last_tracker_of_track
            a->Fill(maxHit[e]); tracksAll += nTrue[e];
            if (evClass[e] == 1){b->Fill(maxHit[e]); tracksCE += nTrueCE[e];} // CE: a CE particle has a pure track, CE + CE mix (blue): every event with a CE particle, as in 1_indivtrackers
            if (evClass[e] >= 1){m->Fill(maxHit[e]); tracksCEMix += nTrueCE[e];}
        }
        a->SetLineColor(CbAll()); a->SetLineWidth(2);
        b->SetLineColor(CbCE()); b->SetLineWidth(2);
        m->SetLineColor(CbProton()); m->SetLineWidth(2);
        const std::string path = "plots/8_trackerfits/5_max_hits_per_track.pdf";
        a->SetTitle(PlotTitle(path, Form("Most trackers entered by one charged particle per event, %lld events with a track of 2+ hits, no cuts;Trackers the longest charged track entered, 1 hit per tracker;Events", nEvLast)).c_str());
        a->SetMinimum(0.5); a->SetMaximum(100.0 * std::max({a->GetMaximum(), b->GetMaximum(), m->GetMaximum()})); // 2 decades of headroom on the log y axis: the legend sits above every curve
        const double keepBottom = canvas.GetBottomMargin();
        canvas.cd();
        canvas.SetBottomMargin(0.14); a->GetXaxis()->SetTitleOffset(1.2); // the x title below the tick labels, not over them
        a->Draw("hist"); m->Draw("hist same"); b->Draw("hist same"); // CE (vermillion) drawn last, on top of CE + CE mix (blue)
        TLegend lg(0.12, 0.74, 0.55, 0.88); // inside the frame: the labels end before the right edge of the axes
        lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.026);
        lg.AddEntry(a, Form("all: %lld tracks, %.0f events", tracksAll, a->GetEntries()), "l");
        lg.AddEntry(b, Form("CE: %lld tracks, %.0f events", tracksCE, b->GetEntries()), "l");
        lg.AddEntry(m, Form("CE + CE mix: %lld tracks, %.0f events", tracksCEMix, m->GetEntries()), "l");
        lg.Draw();
        canvas.Update();
        WritePlotNumbers(a, path, numberFiles);
        canvas.SaveAs(path.c_str());
        canvas.SetBottomMargin(keepBottom);
    }
    // 6_last_tracker_of_track.pdf: for every true track with 2 or more tracker hits, the number of the last (highest z) tracker it hits, 1 = the first tracker (the numbers of 1_indivtrackers); all particles and the CE child protons
    {
        TH1D *a = new TH1D("lastPlaneAll", "", (int)nPl, 0.5, nPl + 0.5), *b = new TH1D("lastPlaneCE", "", (int)nPl, 0.5, nPl + 0.5), *m = new TH1D("lastPlaneMix", "", (int)nPl, 0.5, nPl + 0.5);
        for (int n : lastPlane) a->Fill(n);
        for (int n : lastPlaneCE) b->Fill(n);
        for (int n : lastPlaneMix) m->Fill(n);
        m->SetLineColor(CbProton()); m->SetLineWidth(2); // the colours of the recon plots (ClassColor, F9_PlotReconNew.cpp): CE vermillion, CE mix blue, all green
        a->SetLineColor(CbAll()); a->SetLineWidth(2);
        b->SetLineColor(CbCE()); b->SetLineWidth(2);
        const std::string path = "plots/8_trackerfits/6_last_tracker_of_track.pdf";
        a->SetTitle(PlotTitle(path, Form("Last tracker of a track per track, %lld true tracks with 2+ hits in %lld events, no cuts;Last tracker the track hit;Tracks", nTrueTracks, nEvLast)).c_str());
        a->SetMinimum(0.5); a->SetMaximum(40.0 * std::max({a->GetMaximum(), b->GetMaximum(), m->GetMaximum()})); // headroom for the legend over the bars
        canvas.cd();
        a->Draw("hist"); m->Draw("hist same"); b->Draw("hist same"); // CE (vermillion) drawn last, on top of CE mix (blue)
        TLegend lg(0.12, 0.74, 0.55, 0.90);
        lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.018);
        lg.AddEntry(a, Form("all recon particles: %zu tracks, %lld events", lastPlane.size(), nEvLast), "l");
        lg.AddEntry(b, Form("CE: %zu tracks, %lld events", lastPlaneCE.size(), nEvLastCE), "l");
        lg.AddEntry(m, Form("CE + CE mix: %zu tracks, %lld events", lastPlaneMix.size(), nEvLastMix), "l");
        lg.Draw();
        canvas.Update();
        WritePlotNumbers(a, path, numberFiles);
        canvas.SaveAs(path.c_str());
    }
    // 9 and 10: per event, the tracker hits of the longest track (5_max_hits_per_track) against the last tracker of that track (6_last_tracker_of_track); CE = a CE particle has a pure track, all = every event with a track of 2+ hits
    {
        const int topMax = std::max(1, *std::max_element(maxHit.begin(), maxHit.end()));
        const double keepRight = canvas.GetRightMargin();
        gStyle->SetPalette(kViridis);
        canvas.cd(); canvas.SetLogy(0); canvas.SetRightMargin(0.15);
        for (int kind = 0; kind < 2; kind++){ // 0 CE, 1 all
            TH2D* h = new TH2D(Form("maxVsLast%d", kind), "", topMax, 0.5, topMax + 0.5, (int)nPl, 0.5, nPl + 0.5);
            for (size_t e = 0; e < maxHit.size(); e++){
                if (maxHit[e] < 2 || (kind == 0 && evClass[e] != 1)) continue;
                h->Fill(maxHit[e], maxHitLast[e]);
            }
            const std::string path = kind == 0 ? "plots/8_trackerfits/9_maxhits_vs_last_tracker_CE.pdf" : "plots/8_trackerfits/10_maxhits_vs_last_tracker_all.pdf";
            h->SetTitle(PlotTitle(path, Form("Longest track of the event, %s, %.0f events, no cuts;Trackers the longest charged track entered, 1 hit per tracker;Last tracker of that track", kind == 0 ? "CE" : "all recon", h->GetEntries())).c_str());
            canvas.SetLogz(kind);
            h->SetMinimum(kind == 1 ? 0.5 : 0.0);
            h->Draw("colz text");
            canvas.Update();
            WritePlotNumbers(h, path, numberFiles);
            canvas.SaveAs(path.c_str());
        }
        canvas.SetLogz(0); canvas.SetLogy(1); canvas.SetRightMargin(keepRight);
    }
    // 11_max_hits_clean_chain.pdf: the events of 5_max_hits_per_track by chain length, the trackers in a row from tracker 1 with exactly 1 charged hit each (1,1,1,2 gives 3, 2,1,1 gives 0)
    {
        const int topChain = std::max(1, *std::max_element(chainLen.begin(), chainLen.end()));
        TH1D *a = new TH1D("chainAll", "", topChain + 1, -0.5, topChain + 0.5), *b = new TH1D("chainCE", "", topChain + 1, -0.5, topChain + 0.5), *m = new TH1D("chainCEMix", "", topChain + 1, -0.5, topChain + 0.5);
        long long evAll = 0, evCE = 0, evMix = 0, trAll = 0, trCE = 0, trMix = 0; // events and true tracks with 2+ hits, as in 5_max_hits_per_track
        for (size_t e = 0; e < chainLen.size(); e++){
            if (maxHit[e] < 2) continue;
            a->Fill(chainLen[e]); evAll++; trAll += nTrue[e];
            if (evClass[e] == 1){b->Fill(chainLen[e]); evCE++; trCE += nTrueCE[e];}
            if (evClass[e] >= 1){m->Fill(chainLen[e]); evMix++; trMix += nTrueCE[e];}
        }
        a->SetLineColor(CbAll()); a->SetLineWidth(2);
        b->SetLineColor(CbCE()); b->SetLineWidth(2);
        m->SetLineColor(CbProton()); m->SetLineWidth(2);
        const std::string path = "plots/8_trackerfits/11_max_hits_clean_chain.pdf";
        a->SetTitle(PlotTitle(path, Form("Longest run of trackers in a row with exactly 1 charged hit, %lld events with a track of 2+ hits, no cuts;Longest run of trackers with 1 hit each;Events", evAll)).c_str());
        a->SetMinimum(0.5); a->SetMaximum(100.0 * std::max({a->GetMaximum(), b->GetMaximum(), m->GetMaximum()}));
        const double keepBottom = canvas.GetBottomMargin();
        canvas.cd();
        canvas.SetBottomMargin(0.14); a->GetXaxis()->SetTitleOffset(1.2);
        a->Draw("hist"); m->Draw("hist same"); b->Draw("hist same"); // CE (vermillion) drawn last, on top of CE + CE mix (blue)
        TLegend lg(0.12, 0.74, 0.55, 0.88);
        lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.026);
        lg.AddEntry(a, Form("all: %lld tracks, %lld events", trAll, evAll), "l");
        lg.AddEntry(b, Form("CE: %lld tracks, %lld events", trCE, evCE), "l");
        lg.AddEntry(m, Form("CE + CE mix: %lld tracks, %lld events", trMix, evMix), "l");
        lg.Draw();
        canvas.Update();
        WritePlotNumbers(a, path, numberFiles);
        canvas.SaveAs(path.c_str());
        canvas.SetBottomMargin(keepBottom);
    }
    // 12 and 13: the clean chains of 11_max_hits_clean_chain (1 charged hit on each tracker from tracker 1, one particle) as in 9 and 10: the chain length against the last tracker of the longest track of the event (6_last_tracker_of_track); CE = a CE particle has a pure track, all = every clean chain event
    {
        const int topChain = std::max(2, *std::max_element(chainLen.begin(), chainLen.end()));
        const double keepRight = canvas.GetRightMargin();
        gStyle->SetPalette(kViridis);
        canvas.cd(); canvas.SetLogy(0); canvas.SetRightMargin(0.15);
        for (int kind = 0; kind < 2; kind++){ // 0 CE, 1 all
            TH2D* h = new TH2D(Form("chainVsLast%d", kind), "", topChain, 0.5, topChain + 0.5, (int)nPl, 0.5, nPl + 0.5);
            for (size_t e = 0; e < chainLen.size(); e++){
                if (chainLen[e] < 2 || (kind == 0 && evClass[e] != 1)) continue;
                h->Fill(chainLen[e], maxHitLast[e]);
            }
            const std::string path = kind == 0 ? "plots/8_trackerfits/12_chain_vs_last_tracker_CE.pdf" : "plots/8_trackerfits/13_chain_vs_last_tracker_all.pdf";
            h->SetTitle(PlotTitle(path, Form("Clean chains, %s, %.0f events, no cuts;Trackers of the clean chain, 1 hit per tracker (the first tracker alone is 1 hit);Last tracker of the longest track", kind == 0 ? "CE" : "all recon", h->GetEntries())).c_str());
            canvas.SetLogz(kind);
            h->SetMinimum(kind == 1 ? 0.5 : 0.0);
            h->Draw("colz text");
            canvas.Update();
            WritePlotNumbers(h, path, numberFiles);
            canvas.SaveAs(path.c_str());
        }
        canvas.SetLogz(0); canvas.SetLogy(1); canvas.SetRightMargin(keepRight);
    }
    BoxAndCropPlots();
}
