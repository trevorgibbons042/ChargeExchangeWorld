#include <TError.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/TSeq.hxx>
#include <ROOT/TThreadExecutor.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_set>
#include <vector>

#include "FitCore.h"
#include "QualityCut.h"

// Turns a cuts/cut_angles_all.root of the old 5-tree format (AnglesRaw, Angles, AnglesCE, Angles_DE, AnglesCE_DE) into the one-tree format of ReconAngleAll.cpp (AnglesRaw with
// inCE and finalKED) without refitting a single track. Run after CutEvents.cpp (needs cuts/cut_CE.root).
//   inCE         the event has a carbon CE: its fEvent is in cuts/cut_CE.root (ChargeExchange_values)
//   finalKED     detector KE of the track end: AnglesRaw has only the true KE, so it comes from Angles_DE, whose rows must be the rows of AnglesRaw (same order, same count);
//                a file whose Angles_DE was cut by the removed quality test has fewer rows and is refused
// On any mismatch the old file stays. On success it is kept as cuts/cut_angles_all_old5.root (delete it when the new file is checked).
// All the threads work on it, and no file is merged. AnglesRaw is cut into one block of consecutive entries per thread; a first pass gives every row its finalKED (a table with one value per row);
// then one multithreaded RDataFrame Snapshot copies the tree with the two new branches into the output file. The rows reach the file in the order the Snapshot buffers arrive, not in the order of the old file; nothing after this
// reads them by position (the entry lists, CalcAPandT.cpp and the plots go by fEvent, TrackID and the entry number of the new file).
struct ConvBlock {
    Long64_t first = 0, last = 0;   // the entries [first, last) of AnglesRaw
    Long64_t nCE = 0;               // rows of the block in a CE event
    bool ok = true;
};

void ConvertAngles(){
    const char* inPath = "cuts/cut_angles_all.root";
    const char* outPath = "cuts/cut_angles_all.converted.root";
    const char* oldPath = "cuts/cut_angles_all_old5.root";
    ROOT::EnableThreadSafety();
    const unsigned nThreads = NThreads();

    std::unordered_set<int> ceEvents; // events with a carbon CE
    {
        TFile fCE("cuts/cut_CE.root", "READ");
        TTree* ce = fCE.IsZombie() ? nullptr : (TTree*)fCE.Get("ChargeExchange_values");
        if (!ce){Error("ConvertAngles", "cuts/cut_CE.root or its ChargeExchange_values tree is missing (run CutEvents.cpp)"); return;}
        int ev = 0;
        ce->SetBranchAddress("fEvent", &ev);
        for (Long64_t i = 0, n = ce->GetEntries(); i < n; i++){ce->GetEntry(i); ceEvents.insert(ev);}
    }

    Long64_t nRaw = 0, nAng = 0, nAngCE = 0, nDE = 0;
    std::vector<std::string> columns; // the branches of AnglesRaw
    {
        TFile in(inPath, "READ");
        TTree* raw = in.IsZombie() ? nullptr : (TTree*)in.Get("AnglesRaw");
        TTree* ang = in.IsZombie() ? nullptr : (TTree*)in.Get("Angles");
        TTree* angCE = in.IsZombie() ? nullptr : (TTree*)in.Get("AnglesCE");
        TTree* angDE = in.IsZombie() ? nullptr : (TTree*)in.Get("Angles_DE");
        if (!raw){Error("ConvertAngles", "no AnglesRaw tree in %s", inPath); return;}
        if (raw->GetBranch("inCE")){printf("[ConvertAngles] %s is already in the one-tree format, nothing to do\n", inPath); return;}
        if (!ang || !angCE || !angDE){Error("ConvertAngles", "%s lacks Angles, AnglesCE or Angles_DE: not the old format", inPath); return;}
        nRaw = raw->GetEntries(); nAng = ang->GetEntries(); nAngCE = angCE->GetEntries(); nDE = angDE->GetEntries();
        for (auto* br : *raw->GetListOfBranches()) columns.push_back(br->GetName());
    }

    const Long64_t nBlocks = std::max<Long64_t>(1, std::min<Long64_t>(nThreads, nRaw));
    std::vector<ConvBlock> blocks(nBlocks);
    for (Long64_t k = 0; k < nBlocks; k++){blocks[k].first = k * nRaw / nBlocks; blocks[k].last = (k + 1) * nRaw / nBlocks;}
    ROOT::TThreadExecutor pool(nThreads);

    if (nDE != nRaw){Error("ConvertAngles", "Angles_DE has %lld rows, AnglesRaw %lld: the old file was cut by a quality test that no longer exists. Old file kept", nDE, nRaw); return;}
    // finalKED of every row; each block reads its own stretch of Angles_DE
    std::vector<double> finalKED(nRaw, NAN);
    pool.Foreach([&](unsigned k){
        ConvBlock& b = blocks[k];
        if (!b.ok) return;
        TFile f(inPath, "READ");
        TTree* raw = f.IsZombie() ? nullptr : (TTree*)f.Get("AnglesRaw");
        TTree* angDE = f.IsZombie() ? nullptr : (TTree*)f.Get("Angles_DE");
        if (!raw || !angDE){b.ok = false; return;}
        double finalKEde = 0;
        int fEvent = 0;
        raw->SetBranchStatus("*", 0);
        raw->SetBranchStatus("fEvent", 1);
        raw->SetBranchAddress("fEvent", &fEvent);
        angDE->SetBranchStatus("*", 0); // only finalKE of Angles_DE is read, not its other 37 branches
        angDE->SetBranchStatus("finalKE", 1);
        angDE->SetBranchAddress("finalKE", &finalKEde); // in Angles_DE, finalKE holds the detector KE
        for (Long64_t i = b.first; i < b.last; i++){
            raw->GetEntry(i);
            angDE->GetEntry(i);
            finalKED[i] = finalKEde;
            if (ceEvents.count(fEvent)) b.nCE++;
        }
    }, ROOT::TSeqU((unsigned)nBlocks));

    Long64_t nCE = 0;
    for (const auto& b : blocks){
        if (!b.ok){Error("ConvertAngles", "block %lld to %lld failed. Old file kept", b.first, b.last); return;}
        nCE += b.nCE;
    }
    if (nCE != nAngCE){
        Error("ConvertAngles", "rows in CE events %lld do not match AnglesCE %lld. Old file kept", nCE, nAngCE);
        return;
    }

    // one multithreaded Snapshot (it writes the file through a TBufferMerger of its own, no merge pass afterwards)
    ROOT::EnableImplicitMT(nThreads);
    {
        ROOT::RDataFrame df("AnglesRaw", inPath);
        auto withNew = df.Define("finalKED", [&finalKED](ULong64_t e){return finalKED[e];}, {"rdfentry_"})
                         .Define("inCE", [&ceEvents](int ev){return (int)ceEvents.count(ev);}, {"fEvent"});
        std::vector<std::string> outCols = columns;
        for (const char* c : {"finalKED", "inCE"}) outCols.push_back(c);
        withNew.Snapshot("AnglesRaw", outPath, outCols);
    }
    {
        TFile conv(outPath, "READ");
        TTree* t = conv.IsZombie() ? nullptr : (TTree*)conv.Get("AnglesRaw");
        if (!t || t->GetEntries() != nRaw){Error("ConvertAngles", "%s has %lld rows, not %lld. Old file kept", outPath, t ? (long long)t->GetEntries() : -1LL, nRaw); gSystem->Unlink(outPath); return;}
    }
    if (gSystem->Rename(inPath, oldPath) != 0 || gSystem->Rename(outPath, inPath) != 0){Error("ConvertAngles", "could not move the files into place (%s, %s)", oldPath, outPath); return;}
    printf("[ConvertAngles] %lld rows, %lld of them in CE events; old file kept as %s\n", nRaw, nCE, oldPath);
}
