// Splits AnglesRaw of a run folder into a signal file (CE only, evClass == 0: every fitted hit a CE child proton) and a background file (isFromCE == 0) that TmvaTrain.cpp reads.
// Run in the run folder (build/root/<run>/): root -l -b -q 'TmvaSplit.cpp("cuts/cut_angles_all.root", 2)'
// Writes cuts/cut_tmva_signal.root (tree TmvaSignal) and cuts/cut_tmva_background.root (tree TmvaBackground), each with only the row key (fEvent, trackNo) and the detector-level inputs of TmvaVars.h:
// no truth column, not even the label isFromCE (the file says which class a row is in).
// Every row is kept; a NaN input is written as kTmvaMissing (-999, TmvaVars.h) because TMVA cannot take a NaN. prescale keeps only the events with fEvent % prescale == 0, all the tracks of such an event, so the files stay a size TMVA can load
// (the full 80M rows of a run are about 10 GB); 1 keeps every event. The default 2 trains the BDT on the even events; the odd events are never trained on, and F9_PlotReconNew.cpp runs its BDT cut fits on them only (it reads the prescale from cuts/tmva_prescale.txt, written here). The row order is not kept, the key of a row is fEvent + trackNo.
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"

#include "TmvaVars.h"

void TmvaSplit(const char* input = "cuts/cut_angles_all.root", int prescale = 2, int bkgThin = 3){
    if (gSystem->AccessPathName(input)){Error("TmvaSplit", "no %s", input); return;}
    if (prescale < 1) prescale = 1;
    if (bkgThin < 1) bkgThin = 1; // keeps 1 background event in bkgThin of the events prescale kept; every signal event stays (the signal is about 1 in 5000 of the background)

    const std::vector<TmvaVar> vars = TmvaVariables();
    std::vector<std::string> columns = {"fEvent", "trackNo"}; // evClass and isFromCE only pick the rows below, it is not written
    {
        TFile in(input);
        TTree* raw = (TTree*)in.Get("AnglesRaw");
        if (!raw){Error("TmvaSplit", "no AnglesRaw in %s", input); return;}
        for (const auto& x : vars){
            if (!raw->GetBranch(x.name.c_str())){Error("TmvaSplit", "AnglesRaw has no %s (run ReconAngleAll.cpp of this version)", x.name.c_str()); return;}
            columns.push_back(x.name);
        }
    }

    ROOT::EnableImplicitMT(NThreads());
    ROOT::RDataFrame df("AnglesRaw", input);
    ROOT::RDF::RNode kept = df.Filter(Form("fEvent %% %d == 0", prescale));
    for (const auto& x : vars) if (x.type == 'D') kept = kept.Redefine(x.name, Form("TMath::Finite(%s) ? %s : %g", x.name.c_str(), x.name.c_str(), kTmvaMissing)); // NaN becomes kTmvaMissing, no row is dropped
    auto signal = kept.Filter("evClass == 0"); // CE only: every fitted hit is a CE child proton. CE mix events (evClass 1) are in neither class, a mixed event is not a clean signal or background
    auto background = kept.Filter(Form("isFromCE == 0 && fEvent %% %d == 0", prescale * bkgThin));
    // one pass over AnglesRaw per file (a Snapshot runs when it is booked); the row count is read back from the file it wrote
    auto signalFile = signal.Snapshot("TmvaSignal", "cuts/cut_tmva_signal.root", columns);
    auto backgroundFile = background.Snapshot("TmvaBackground", "cuts/cut_tmva_background.root", columns);
    const ULong64_t nSig = *signalFile->Count(), nBkg = *backgroundFile->Count();
    std::cout << "==> cuts/cut_tmva_signal.root " << nSig << " rows, cuts/cut_tmva_background.root " << nBkg << " rows (prescale " << prescale << ", background 1 in " << bkgThin << " of those)" << std::endl;
    {std::ofstream out("cuts/tmva_prescale.txt"); out << prescale << "\n";} // the events with fEvent % prescale == 0 trained the BDT, F9_PlotReconNew.cpp leaves them out of the BDT cut fits
}
