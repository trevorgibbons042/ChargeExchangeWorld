#include <TFile.h>
#include <TError.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <unordered_set>

#include "FitCore.h"
#include "QualityCut.h"

// Turns a cuts/cut_angles_all.root of the old 5-tree format (AnglesRaw, Angles, AnglesCE, Angles_DE, AnglesCE_DE) into the one-tree format of ReconAngleAll.cpp (AnglesRaw with
// passQuality, inCE and finalKED) without refitting a single track. Run after CutEvents.cpp (needs cuts/cut_CE.root).
//   passQuality  zOrderOK and PassesQuality, the test ReconAngleAll.cpp fills Angles with
//   inCE         the event has a carbon CE: its fEvent is in cuts/cut_CE.root (ChargeExchange_values)
//   finalKED     detector KE of the track end: AnglesRaw has only the true KE, so it comes from Angles_DE, which holds the passQuality rows in the same order; NaN for the other rows
// The result must reproduce the old trees: the passQuality rows number the entries of Angles, the passQuality and inCE rows those of AnglesCE, and Angles_DE is used up exactly.
// On any mismatch the old file stays. On success it is kept as cuts/cut_angles_all_old5.root (delete it when the new file is checked).
void ConvertAngles(){
    const char* inPath = "cuts/cut_angles_all.root";
    const char* outPath = "cuts/cut_angles_all.converted.root";
    const char* oldPath = "cuts/cut_angles_all_old5.root";
    ROOT::EnableImplicitMT(NThreads()); // parallel compression of the output baskets

    std::unordered_set<int> ceEvents; // events with a carbon CE
    {
        TFile fCE("cuts/cut_CE.root", "READ");
        TTree* ce = fCE.IsZombie() ? nullptr : (TTree*)fCE.Get("ChargeExchange_values");
        if (!ce){Error("ConvertAngles", "cuts/cut_CE.root or its ChargeExchange_values tree is missing (run CutEvents.cpp)"); return;}
        int ev = 0;
        ce->SetBranchAddress("fEvent", &ev);
        for (Long64_t i = 0, n = ce->GetEntries(); i < n; i++){ce->GetEntry(i); ceEvents.insert(ev);}
    }

    TFile in(inPath, "READ");
    TTree* raw = in.IsZombie() ? nullptr : (TTree*)in.Get("AnglesRaw");
    TTree* ang = in.IsZombie() ? nullptr : (TTree*)in.Get("Angles");
    TTree* angCE = in.IsZombie() ? nullptr : (TTree*)in.Get("AnglesCE");
    TTree* angDE = in.IsZombie() ? nullptr : (TTree*)in.Get("Angles_DE");
    if (!raw){Error("ConvertAngles", "no AnglesRaw tree in %s", inPath); return;}
    if (raw->GetBranch("passQuality")){printf("[ConvertAngles] %s is already in the one-tree format, nothing to do\n", inPath); return;}
    if (!ang || !angCE || !angDE){Error("ConvertAngles", "%s lacks Angles, AnglesCE or Angles_DE: not the old format", inPath); return;}

    double resRms = 0, impact = 0, impactD = 0, finalKEde = 0;
    int zOrderOK = 0, fEvent = 0;
    raw->SetBranchAddress("resRms", &resRms); raw->SetBranchAddress("impact", &impact); raw->SetBranchAddress("impactD", &impactD);
    raw->SetBranchAddress("zOrderOK", &zOrderOK); raw->SetBranchAddress("fEvent", &fEvent);
    angDE->SetBranchStatus("*", 0); // only finalKE of Angles_DE is read, not its other 37 branches
    angDE->SetBranchStatus("finalKE", 1);
    angDE->SetBranchAddress("finalKE", &finalKEde); // in Angles_DE, finalKE holds the detector KE

    TFile out(outPath, "RECREATE");
    TTree* conv = raw->CloneTree(0); // same branches, filled from raw's addresses below
    conv->SetName("AnglesRaw");
    int passQuality = 0, inCE = 0;
    double finalKED = NAN;
    conv->Branch("finalKED", &finalKED);
    conv->Branch("passQuality", &passQuality);
    conv->Branch("inCE", &inCE);

    Long64_t nPass = 0, nPassCE = 0;
    const Long64_t nRaw = raw->GetEntries(), nDE = angDE->GetEntries();
    for (Long64_t i = 0; i < nRaw; i++){
        raw->GetEntry(i);
        passQuality = zOrderOK && PassesQuality(resRms, impact, impactD) ? 1 : 0;
        inCE = ceEvents.count(fEvent) ? 1 : 0;
        finalKED = NAN;
        if (passQuality){
            if (nPass >= nDE){Error("ConvertAngles", "more quality rows than Angles_DE entries (%lld): QualityCut.h changed since ReconAngleAll.cpp ran? Old file kept", nDE); return;}
            angDE->GetEntry(nPass);
            finalKED = finalKEde;
            nPass++;
            if (inCE) nPassCE++;
        }
        conv->Fill();
    }
    conv->Write();
    out.Close();

    if (nPass != ang->GetEntries() || nPass != nDE || nPassCE != angCE->GetEntries()){
        Error("ConvertAngles", "rows do not match the old trees: quality %lld vs Angles %lld, Angles_DE %lld; quality and CE %lld vs AnglesCE %lld. Old file kept, %s is not used",
              nPass, (long long)ang->GetEntries(), (long long)nDE, nPassCE, (long long)angCE->GetEntries(), outPath);
        return;
    }
    in.Close();
    if (gSystem->Rename(inPath, oldPath) != 0 || gSystem->Rename(outPath, inPath) != 0){Error("ConvertAngles", "could not move the files into place (%s, %s)", oldPath, outPath); return;}
    printf("[ConvertAngles] %lld rows, %lld pass the quality cut, %lld of them in CE events; old file kept as %s\n", nRaw, nPass, nPassCE, oldPath);
}
