#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TFileMerger.h>

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
    auto dfCE = dfCEAll.Filter(CarbonCeCut());

    // Lazy snapshots: all cuts of one dataframe are written in a single pass over its input, and the three
    // dataframes run at the same time. Each tree goes to its own temp file (two snapshots must not write into
    // one file at once); the temp files are merged into the final cut files afterwards.
    ROOT::RDF::RSnapshotOptions lazy;
    lazy.fLazy = true;

    const std::vector<std::string> cutFiles = {"cuts/cut_CE.root", "cuts/cut_ParentID.root"};
    const std::vector<std::string> treeNames = {events->GetName(), ChargeExchange->GetName()};
    auto tmpName = [](const std::string& tree, const std::string& cutFile){return "cuts/tmp_" + tree + "_" + cutFile.substr(5);}; // cutFile is cuts/<name>

    auto selectedEvent = df.Filter(inEvent, {"fEvent"});
    auto selectedParent = df.Filter(isParent, {"fEvent", "ParentID"});


    std::vector<ROOT::RDF::RResultHandle> handles;
    auto queue = [&](auto& node, size_t tree, size_t cut){
        handles.emplace_back(node.Snapshot(treeNames[tree], tmpName(treeNames[tree], cutFiles[cut]), "", lazy));
    };
    queue(selectedEvent, 0, 0);   queue(selectedParent, 0, 1);
    queue(dfCE, 1, 0);            queue(dfCE, 1, 1);
    // Single-tree cut files (no merge step), see CutOpen.h: ZDC protons of every event, the ZDC and first tracker protons of the CE children, and the CE outside the carbon blocks.
    handles.emplace_back(df.Filter("boundaryType == 2 && pdg == 2212").Snapshot(treeNames[0], "cuts/cut_AllEvProtons.root", "", lazy));
    handles.emplace_back(selectedParent.Filter("boundaryType == 2 && pdg == 2212").Snapshot(treeNames[0], "cuts/cut_ParentZdc.root", "", lazy));
    handles.emplace_back(selectedParent.Filter("boundaryType == 0 && pdg == 2212").Snapshot(treeNames[0], "cuts/cut_ParentTracker.root", "", lazy));
    handles.emplace_back(dfCEAll.Filter("!(" + CarbonCeCut() + ")").Snapshot(treeNames[1], "cuts/cut_CEtracker.root", "", lazy));
    handles.emplace_back(dfCEAll.Snapshot(treeNames[1], "cuts/cut_CEall.root", "", lazy)); // every CE of the sim (the whole ChargeExchange_values tree, one row per CE)
    ROOT::RDF::RunGraphs(handles);
    ROOT::DisableImplicitMT(); // merging with implicit MT on crashed in TTree::FlushBasketsImpl

    for (size_t cut = 0; cut < cutFiles.size(); cut++) {
        {
            TFileMerger merger;
            merger.OutputFile(cutFiles[cut].c_str(), "RECREATE");
            for (const auto& tree : treeNames) merger.AddFile(tmpName(tree, cutFiles[cut]).c_str());
            if (!merger.Merge()) std::cerr << "CutEvents: merge into " << cutFiles[cut] << " failed\n";
        }
        for (const auto& tree : treeNames) std::remove(tmpName(tree, cutFiles[cut]).c_str());
    }
}
