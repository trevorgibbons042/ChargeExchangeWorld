// Splits AnglesRaw of a run folder into a signal file (the fit used a CE child hit, isFromCE == 1) and a background file (isFromCE == 0) that TmvaTrain.cpp reads.
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
#include "ROOT/RDFHelpers.hxx" // RunGraphs
#include "ROOT/RSnapshotOptions.hxx"

#include "TmvaVars.h"
#include "CutFiles.h"

void TmvaSplit(const char* input = "cuts/cut_angles_all.root", int prescale = 2){
    if (gSystem->AccessPathName(input)){Error("TmvaSplit", "no %s", input); return;}
    if (prescale < 1) prescale = 1;

    const std::vector<TmvaVar> vars = TmvaVariables();
    std::vector<std::string> columns = {"fEvent", "trackNo"}; // isFromCE only picks the rows below, it is not written
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
    TFile* fin = TFile::Open(input, "READ");
    TTree* rawTree = fin ? (TTree*)fin->Get("AnglesRaw") : nullptr;
    if (!rawTree){Error("TmvaSplit", "no AnglesRaw in %s", input); return;}
    AttachEvents(rawTree); // zdcKEAllD and nChargedParticle are columns of the Events tree (cuts/cut_events.root), a friend of AnglesRaw
    ROOT::RDataFrame df(*rawTree);
    ROOT::RDF::RNode kept = df.Filter(Form("fEvent %% %d == 0", prescale));
    for (const auto& x : vars) if (x.type == 'D') kept = kept.Redefine(x.name, Form("TMath::Finite(%s) ? %s : %g", x.name.c_str(), x.name.c_str(), kTmvaMissing)); // NaN becomes kTmvaMissing, no row is dropped
    auto signal = kept.Filter("isFromCE == 1");
    auto background = kept.Filter("isFromCE == 0");
    // both Snapshots are lazy and RunGraphs runs them in one pass over AnglesRaw (an eager Snapshot runs when it is booked, one pass each); the row counts are booked in the same pass
    ROOT::RDF::RSnapshotOptions lazy; lazy.fLazy = true;
    auto signalFile = signal.Snapshot("TmvaSignal", "cuts/cut_tmva_signal.root", columns, lazy);
    auto backgroundFile = background.Snapshot("TmvaBackground", "cuts/cut_tmva_background.root", columns, lazy);
    auto nSigRes = signal.Count(), nBkgRes = background.Count(); // booked with the Snapshots: the same pass counts the rows, no read of the files just written
    ROOT::RDF::RunGraphs({signalFile, backgroundFile, nSigRes, nBkgRes});
    const ULong64_t nSig = *nSigRes, nBkg = *nBkgRes;
    std::cout << "==> cuts/cut_tmva_signal.root " << nSig << " rows, cuts/cut_tmva_background.root " << nBkg << " rows (prescale " << prescale << ")" << std::endl;
    {std::ofstream out("cuts/tmva_prescale.txt"); out << prescale << "\n";} // the events with fEvent % prescale == 0 trained the BDT, F9_PlotReconNew.cpp leaves them out of the BDT cut fits
}
