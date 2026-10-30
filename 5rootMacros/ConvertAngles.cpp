#include <TError.h>
#include <TFile.h>
#include <TFileMerger.h>
#include <TSystem.h>
#include <TTree.h>
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
// passQuality, inCE and finalKED) without refitting a single track. Run after CutEvents.cpp (needs cuts/cut_CE.root).
//   passQuality  zOrderOK and PassesQuality, the test ReconAngleAll.cpp fills Angles with
//   inCE         the event has a carbon CE: its fEvent is in cuts/cut_CE.root (ChargeExchange_values)
//   finalKED     detector KE of the track end: AnglesRaw has only the true KE, so it comes from Angles_DE, which holds the passQuality rows in the same order; NaN for the other rows
// The result must reproduce the old trees: the passQuality rows number the entries of Angles, the passQuality and inCE rows those of AnglesCE, and Angles_DE is used up exactly.
// On any mismatch the old file stays. On success it is kept as cuts/cut_angles_all_old5.root (delete it when the new file is checked).
// All the threads work on it: AnglesRaw is cut into one block of consecutive entries per thread. A first pass counts the passQuality rows of every block, which tells each block
// where its rows start in Angles_DE; a second pass reads, converts and compresses every block into its own file; the files are then merged in block order, so the rows keep
// the order of the old file.
struct ConvBlock {
    Long64_t first = 0, last = 0;   // the entries [first, last) of AnglesRaw
    Long64_t deFirst = 0;           // the Angles_DE row of the block's first passQuality row
    Long64_t nPass = 0, nPass2 = 0, nPassCE = 0; // passQuality rows (first pass, second pass) and those in a CE event
    int compression = 0;
    std::string tmp;
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
    {
        TFile in(inPath, "READ");
        TTree* raw = in.IsZombie() ? nullptr : (TTree*)in.Get("AnglesRaw");
        TTree* ang = in.IsZombie() ? nullptr : (TTree*)in.Get("Angles");
        TTree* angCE = in.IsZombie() ? nullptr : (TTree*)in.Get("AnglesCE");
        TTree* angDE = in.IsZombie() ? nullptr : (TTree*)in.Get("Angles_DE");
        if (!raw){Error("ConvertAngles", "no AnglesRaw tree in %s", inPath); return;}
        if (raw->GetBranch("passQuality")){printf("[ConvertAngles] %s is already in the one-tree format, nothing to do\n", inPath); return;}
        if (!ang || !angCE || !angDE){Error("ConvertAngles", "%s lacks Angles, AnglesCE or Angles_DE: not the old format", inPath); return;}
        nRaw = raw->GetEntries(); nAng = ang->GetEntries(); nAngCE = angCE->GetEntries(); nDE = angDE->GetEntries();
    }

    const Long64_t nBlocks = std::max<Long64_t>(1, std::min<Long64_t>(nThreads, nRaw));
    std::vector<ConvBlock> blocks(nBlocks);
    for (Long64_t k = 0; k < nBlocks; k++){blocks[k].first = k * nRaw / nBlocks; blocks[k].last = (k + 1) * nRaw / nBlocks;}
    ROOT::TThreadExecutor pool(nThreads);

    // pass 1: the passQuality rows of every block (only the four branches the test reads)
    pool.Foreach([&](unsigned k){
        ConvBlock& b = blocks[k];
        TFile f(inPath, "READ");
        TTree* t = f.IsZombie() ? nullptr : (TTree*)f.Get("AnglesRaw");
        if (!t){b.ok = false; return;}
        double resRms = 0, impact = 0, impactD = 0;
        int zOrderOK = 0;
        t->SetBranchStatus("*", 0);
        for (const char* br : {"resRms", "impact", "impactD", "zOrderOK"}) t->SetBranchStatus(br, 1);
        t->SetBranchAddress("resRms", &resRms); t->SetBranchAddress("impact", &impact); t->SetBranchAddress("impactD", &impactD); t->SetBranchAddress("zOrderOK", &zOrderOK);
        for (Long64_t i = b.first; i < b.last; i++){
            t->GetEntry(i);
            if (zOrderOK && PassesQuality(resRms, impact, impactD)) b.nPass++;
        }
    }, ROOT::TSeqU((unsigned)nBlocks));
    Long64_t rowsBefore = 0;
    for (auto& b : blocks){b.deFirst = rowsBefore; rowsBefore += b.nPass;}

    // pass 2: every block into its own file
    pool.Foreach([&](unsigned k){
        ConvBlock& b = blocks[k];
        if (!b.ok) return;
        TFile f(inPath, "READ");
        TTree* raw = f.IsZombie() ? nullptr : (TTree*)f.Get("AnglesRaw");
        TTree* angDE = f.IsZombie() ? nullptr : (TTree*)f.Get("Angles_DE");
        if (!raw || !angDE){b.ok = false; return;}
        double resRms = 0, impact = 0, impactD = 0, finalKEde = 0;
        int zOrderOK = 0, fEvent = 0;
        raw->SetBranchAddress("resRms", &resRms); raw->SetBranchAddress("impact", &impact); raw->SetBranchAddress("impactD", &impactD);
        raw->SetBranchAddress("zOrderOK", &zOrderOK); raw->SetBranchAddress("fEvent", &fEvent);
        angDE->SetBranchStatus("*", 0); // only finalKE of Angles_DE is read, not its other 37 branches
        angDE->SetBranchStatus("finalKE", 1);
        angDE->SetBranchAddress("finalKE", &finalKEde); // in Angles_DE, finalKE holds the detector KE

        b.tmp = Form("cuts/tmp_conv_%u.root", k);
        TFile out(b.tmp.c_str(), "RECREATE");
        b.compression = out.GetCompressionSettings();
        TTree* conv = raw->CloneTree(0); // same branches, filled from raw's addresses below
        conv->SetName("AnglesRaw");
        int passQuality = 0, inCE = 0;
        double finalKED = NAN;
        conv->Branch("finalKED", &finalKED);
        conv->Branch("passQuality", &passQuality);
        conv->Branch("inCE", &inCE);

        Long64_t deRow = b.deFirst;
        for (Long64_t i = b.first; i < b.last; i++){
            raw->GetEntry(i);
            passQuality = zOrderOK && PassesQuality(resRms, impact, impactD) ? 1 : 0;
            inCE = ceEvents.count(fEvent) ? 1 : 0;
            finalKED = NAN;
            if (passQuality){
                if (deRow >= nDE){b.ok = false; return;}
                angDE->GetEntry(deRow++);
                finalKED = finalKEde;
                b.nPass2++;
                if (inCE) b.nPassCE++;
            }
            conv->Fill();
        }
        conv->Write();
        out.Close();
    }, ROOT::TSeqU((unsigned)nBlocks));

    auto removeTmp = [&](){for (const auto& b : blocks) if (!b.tmp.empty()) gSystem->Unlink(b.tmp.c_str());};
    Long64_t nPass = 0, nPassCE = 0;
    for (const auto& b : blocks){
        if (!b.ok || b.nPass != b.nPass2){
            Error("ConvertAngles", "block %lld to %lld failed (more quality rows than Angles_DE entries, or the two passes disagree: QualityCut.h changed since ReconAngleAll.cpp ran?). Old file kept",
                  b.first, b.last);
            removeTmp(); return;
        }
        nPass += b.nPass; nPassCE += b.nPassCE;
    }
    if (nPass != nAng || nPass != nDE || nPassCE != nAngCE){
        Error("ConvertAngles", "rows do not match the old trees: quality %lld vs Angles %lld, Angles_DE %lld; quality and CE %lld vs AnglesCE %lld. Old file kept",
              nPass, nAng, nDE, nPassCE, nAngCE);
        removeTmp(); return;
    }

    // the blocks in order into one file: the baskets are copied, not recompressed, when the compression is the same
    {
        TFileMerger merger;
        merger.OutputFile(outPath, "RECREATE", blocks[0].compression);
        for (const auto& b : blocks) merger.AddFile(b.tmp.c_str());
        if (!merger.Merge()){Error("ConvertAngles", "merging the blocks into %s failed. Old file kept", outPath); removeTmp(); return;}
    }
    removeTmp();
    {
        TFile conv(outPath, "READ");
        TTree* t = conv.IsZombie() ? nullptr : (TTree*)conv.Get("AnglesRaw");
        if (!t || t->GetEntries() != nRaw){Error("ConvertAngles", "%s has %lld rows, not %lld. Old file kept", outPath, t ? (long long)t->GetEntries() : -1LL, nRaw); return;}
    }
    if (gSystem->Rename(inPath, oldPath) != 0 || gSystem->Rename(outPath, inPath) != 0){Error("ConvertAngles", "could not move the files into place (%s, %s)", oldPath, outPath); return;}
    printf("[ConvertAngles] %lld rows, %lld pass the quality cut, %lld of them in CE events; old file kept as %s\n", nRaw, nPass, nPassCE, oldPath);
}
