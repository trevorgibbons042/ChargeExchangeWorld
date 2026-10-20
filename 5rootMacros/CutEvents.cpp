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
#include <ROOT/RSnapshotOptions.hxx>

void CutEvents(){
    ROOT::EnableImplicitMT(128);
    
    TChain *events = new TChain("ParticleInfo");
    TChain *ChargeExchange = new TChain("ChargeExchange_values");
    
    TString filename = "output0_t*.root";
    events->Add(filename);
    ChargeExchange->Add(filename);
    
    Long64_t Sync = ChargeExchange->Draw("fEvent:TrackID", "didChargeExchange>0 && fZ<7000", "goff");

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
    // event has a charge exchange track that is not this track's parent
    auto isFake = [&](int ev, int parent){
        auto it = ceTracksInEvent.find(ev);
        if (it == ceTracksInEvent.end()) return false;
        for (int id : it->second) {if (id != parent) return true;}
        return false;
    };

    ROOT::RDataFrame df(*events);
    ROOT::RDataFrame dfCEAll(*ChargeExchange);
    auto dfCE = dfCEAll.Filter("fZ<7000");

    // Lazy snapshots: all cuts of one dataframe are written in a single pass over its input, and the three
    // dataframes run at the same time. Each tree goes to its own temp file (two snapshots must not write into
    // one file at once); the temp files are merged into the final cut files afterwards.
    ROOT::RDF::RSnapshotOptions lazy;
    lazy.fLazy = true;

    const std::vector<std::string> cutFiles = {"cut_CE.root", "cut_ParentID.root", "cut_FakeID.root"};
    const std::vector<std::string> treeNames = {events->GetName(), ChargeExchange->GetName()};
    auto tmpName = [](const std::string& tree, const std::string& cutFile){return "tmp_" + tree + "_" + cutFile;};

    auto selectedEvent = df.Filter(inEvent, {"fEvent"});
    auto selectedParent = df.Filter(isParent, {"fEvent", "ParentID"});
    auto selectedFake = df.Filter(isFake, {"fEvent", "ParentID"});


    std::vector<ROOT::RDF::RResultHandle> handles;
    auto queue = [&](auto& node, size_t tree, size_t cut){
        handles.emplace_back(node.Snapshot(treeNames[tree], tmpName(treeNames[tree], cutFiles[cut]), "", lazy));
    };
    queue(selectedEvent, 0, 0);   queue(selectedParent, 0, 1);   queue(selectedFake, 0, 2);
    queue(dfCE, 1, 0);            queue(dfCE, 1, 1);             queue(dfCE, 1, 2);
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
