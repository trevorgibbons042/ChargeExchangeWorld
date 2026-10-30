#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>

#include <cstdio>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include "FitCore.h"
#include <ROOT/RSnapshotOptions.hxx>
#include <ROOT/RDFHelpers.hxx> // RunGraphs, RResultHandle

void CutEvents(){
    ROOT::EnableImplicitMT(NThreads());
    
    TChain *events = new TChain("ParticleInfo");
    TChain *ChargeExchange = new TChain("ChargeExchange_values");
    
    TString filename = "output0_t*.root";
    events->Add(filename);
    ChargeExchange->Add(filename);
    
    Long64_t Sync = ChargeExchange->Draw("fEvent:TrackID", CarbonCeCut().c_str(), "goff");

    // Lookup tables instead of one huge "a || b || ..." string per cut: those strings were JIT-compiled by cling
    // (serial, slow) and tested term by term for every entry.
    std::unordered_map<int, std::vector<int>> ceTracksInEvent;
    std::unordered_set<long long> ceEventParent;
    for (Long64_t i = 0; i < Sync; i++) {
        int ev = (int)ChargeExchange->GetV1()[i];
        int id = (int)ChargeExchange->GetV2()[i];
        ceTracksInEvent[ev].push_back(id);
        ceEventParent.insert(((long long)ev << 32) | (unsigned int)id);
    }

    // event has a charge exchange
    auto inEvent = [&](int ev){return ceTracksInEvent.count(ev) > 0;};
    // track is a child of a charge exchange track in this event
    auto isParent = [&](int ev, int parent){return ceEventParent.count(((long long)ev << 32) | (unsigned int)parent) > 0;};

    ROOT::RDataFrame df(*events);
    ROOT::RDataFrame dfCEAll(*ChargeExchange);

    // The cuts chain off each other (see CutOpen.h): the big sim output is read once, in stage 1, for the three base files; the other two cuts are filters of the small cut_CE.
    //   stage 1 (sim output):  cut_CE (events with a carbon CE), cut_AllEvProtons (ZDC protons of every event), cut_CEall (every CE)
    //   stage 2 (cut_CE):      cut_ParentZdc (ZDC protons that are children of a carbon CE track), cut_ParentTracker (the same at the first tracker)
    //   stage 3 (cut_CEall):   the carbon CE rows, added as the second tree (UPDATE) of cut_CE
    // Each stage is lazy: all the cuts of a stage are written in one pass over their input. The result handles own the open output files, so they are cleared (files closed)
    // before the next stage reads those files; left open, cut_CE.root came out without ChargeExchange_values.
    ROOT::RDF::RSnapshotOptions lazy;
    lazy.fLazy = true;
    const std::string pInfo = events->GetName(), ceTree = ChargeExchange->GetName();

    {
        std::vector<ROOT::RDF::RResultHandle> handles;
        handles.emplace_back(df.Filter(inEvent, {"fEvent"}).Snapshot(pInfo, "cuts/cut_CE.root", "", lazy));
        handles.emplace_back(df.Filter("boundaryType == 2 && pdg == 2212").Snapshot(pInfo, "cuts/cut_AllEvProtons.root", "", lazy));
        handles.emplace_back(dfCEAll.Snapshot(ceTree, "cuts/cut_CEall.root", "", lazy)); // every CE of the sim (the whole ChargeExchange_values tree, one row per CE)
        ROOT::RDF::RunGraphs(handles);
    }
    {
        std::vector<ROOT::RDF::RResultHandle> handles;
        ROOT::RDataFrame dfCEfile(pInfo, "cuts/cut_CE.root"); // stays alive until RunGraphs
        auto children = dfCEfile.Filter(isParent, {"fEvent", "ParentID"});
        handles.emplace_back(children.Filter("boundaryType == 2 && pdg == 2212").Snapshot(pInfo, "cuts/cut_ParentZdc.root", "", lazy));
        handles.emplace_back(children.Filter("boundaryType == 0 && pdg == 2212").Snapshot(pInfo, "cuts/cut_ParentTracker.root", "", lazy));
        ROOT::RDF::RunGraphs(handles);
    }
    ROOT::RDF::RSnapshotOptions update;
    update.fMode = "UPDATE";
    ROOT::RDataFrame dfCEallRows(ceTree, "cuts/cut_CEall.root");
    dfCEallRows.Filter(CarbonCeCut()).Snapshot(ceTree, "cuts/cut_CE.root", "", update); // the carbon-block CE rows, beside ParticleInfo
}
