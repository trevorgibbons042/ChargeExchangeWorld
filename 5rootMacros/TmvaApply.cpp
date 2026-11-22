// Applies the BDT of TmvaTrain.cpp (dataset/weights/SigBkg_BDT.weights.xml) to every row of AnglesRaw and writes the score to a small tree TMVAScores:
// fEvent and trackNo (the key of a row), isFromCE and BDT. Every row is scored; a NaN input goes in as kTmvaMissing (-999, TmvaVars.h) like in TmvaSplit.cpp. The AnglesRaw tree itself is not copied. The rows of TMVAScores are in the order of AnglesRaw.
// Run in the run folder (build/root/<run>/): root -l -b -q 'TmvaApply.cpp("cuts/cut_angles_all.root", "cuts/tmva_scores.root")'
// Multithreaded: NThreads() workers (FitCore.h), each scores its own slice of the rows with its own TMVA::Reader and its own open file (a TTree and a Reader are not shared between threads).
// Adapted from TMVA_apply_for_Grace.C: no std::filesystem (cling on the cluster), the variables come from TmvaVars.h, no CopyTree of the input.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "TmvaVars.h"
#include "CutFiles.h"

void TmvaApply(const char* input = "cuts/cut_angles_all.root", const char* output = "cuts/tmva_scores.root"){
    const auto t0 = std::chrono::steady_clock::now();
    const char* weights = "dataset/weights/SigBkg_BDT.weights.xml";
    if (gSystem->AccessPathName(input)){Error("TmvaApply", "no %s", input); return;}
    if (gSystem->AccessPathName(weights)){Error("TmvaApply", "no %s (run TmvaTrain.cpp first)", weights); return;}

    const std::vector<TmvaVar> vars = TmvaVariables();
    Long64_t n = 0;
    {
        TFile in(input);
        TTree* raw = (TTree*)in.Get("AnglesRaw");
        if (!raw){Error("TmvaApply", "no AnglesRaw in %s", input); return;}
        AttachEvents(raw); // zdcKEAllD and nChargedParticle are columns of the Events tree (cuts/cut_events.root), a friend of AnglesRaw
        for (const auto& x : vars) if (!raw->GetBranch(x.name.c_str())){Error("TmvaApply", "AnglesRaw has no %s", x.name.c_str()); return;}
        n = raw->GetEntries();
    }

    ROOT::EnableThreadSafety();
    TMVA::Tools::Instance();
    // one slice of the rows per worker: [lo, hi), its own float inputs and Reader (the Reader keeps the addresses of the floats, so the slices are made before and never moved) and its own results
    struct Slice {
        Long64_t lo = 0, hi = 0;
        std::vector<Float_t> in;
        TMVA::Reader* reader = nullptr;
        std::vector<Int_t> fEvent, trackNo, isFromCE;
        std::vector<Double_t> bdt;
    };
    const unsigned nWork = (unsigned)std::max<Long64_t>(1, std::min<Long64_t>((Long64_t)NThreads(), n / 1000 + 1));
    std::vector<Slice> slices(nWork);
    // the Readers are made one after the other, BookMVA is not safe to run on several threads at once
    for (unsigned k = 0; k < nWork; k++){
        Slice& s = slices[k];
        s.lo = n * k / nWork; s.hi = n * (k + 1) / nWork;
        s.in.assign(vars.size(), 0.f);
        s.reader = new TMVA::Reader("!Color:!Silent");
        for (std::size_t i = 0; i < vars.size(); i++) s.reader->AddVariable(vars[i].name.c_str(), &s.in[i]);
        s.reader->BookMVA("BDT", weights);
    }
    std::cout << "==> " << nWork << " threads scoring " << n << " rows" << std::endl;

    std::vector<std::thread> pool;
    for (unsigned k = 0; k < nWork; k++){
        pool.emplace_back([&, k]{
            Slice& s = slices[k];
            TFile in(input);
            TTree* raw = (TTree*)in.Get("AnglesRaw");
            TFile evFile("cuts/cut_events.root"); // a file handle of this thread: the shared one of EventsTree() is not used across threads
            if (TTree* ev = (TTree*)evFile.Get("Events")) raw->AddFriend(ev); // the Events columns of the variables (zdcKEAllD, nChargedParticle), indexed on fEvent
            // the branches keep their own type, the reader takes floats
            std::vector<Double_t> asDouble(vars.size(), 0.0);
            std::vector<Int_t> asInt(vars.size(), 0);
            Int_t fEvent = 0, trackNo = 0, isFromCE = 0;
            raw->SetBranchStatus("*", 0);
            for (std::size_t i = 0; i < vars.size(); i++){
                raw->SetBranchStatus(vars[i].name.c_str(), 1);
                if (vars[i].type == 'D') raw->SetBranchAddress(vars[i].name.c_str(), &asDouble[i]);
                else raw->SetBranchAddress(vars[i].name.c_str(), &asInt[i]);
            }
            for (const char* b : {"fEvent", "trackNo", "isFromCE"}) raw->SetBranchStatus(b, 1);
            raw->SetBranchAddress("fEvent", &fEvent);
            raw->SetBranchAddress("trackNo", &trackNo);
            raw->SetBranchAddress("isFromCE", &isFromCE);
            const std::size_t rows = (std::size_t)(s.hi - s.lo);
            s.fEvent.reserve(rows); s.trackNo.reserve(rows); s.isFromCE.reserve(rows); s.bdt.reserve(rows);
            for (Long64_t e = s.lo; e < s.hi; e++){
                raw->GetEntry(e);
                for (std::size_t i = 0; i < vars.size(); i++) s.in[i] = (Float_t)TmvaValue(vars[i].type == 'D' ? asDouble[i] : (double)asInt[i]); // NaN becomes kTmvaMissing like in TmvaSplit.cpp
                s.fEvent.push_back(fEvent); s.trackNo.push_back(trackNo); s.isFromCE.push_back(isFromCE);
                s.bdt.push_back(s.reader->EvaluateMVA("BDT"));
            }
        });
    }
    for (auto& t : pool) t.join();
    for (auto& s : slices){delete s.reader; s.reader = nullptr;}

    // the small output tree is filled here, one slice after the other, so its rows are in the order of AnglesRaw
    TFile outFile(output, "RECREATE");
    TTree scores("TMVAScores", "BDT score of every row of AnglesRaw");
    Int_t fEvent = 0, trackNo = 0, isFromCE = 0;
    Double_t bdt = 0;
    scores.Branch("fEvent", &fEvent);
    scores.Branch("trackNo", &trackNo);
    scores.Branch("isFromCE", &isFromCE);
    scores.Branch("BDT", &bdt);
    for (const auto& s : slices){
        for (std::size_t r = 0; r < s.bdt.size(); r++){
            fEvent = s.fEvent[r]; trackNo = s.trackNo[r]; isFromCE = s.isFromCE[r]; bdt = s.bdt[r];
            scores.Fill();
        }
    }
    scores.Write();
    outFile.Close();
    std::cout << "==> Wrote " << n << " BDT scores to " << output << " in " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s" << std::endl;
}
