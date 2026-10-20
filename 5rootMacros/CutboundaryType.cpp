
#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TTree.h>
#include <TFileMerger.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <unordered_set>
#include <TStyle.h>

void CutboundaryType(){
    ROOT::EnableImplicitMT(128);

    TFile* File = TFile::Open("cut_CE.root", "READ");
    //TFile* File = TFile::Open("cut_ParentID.root", "READ");
    
    TTree *events = (TTree*)File->Get("ParticleInfo");
    TTree *ChargeExchange = (TTree*)File->Get("ChargeExchange_values");
    
    Long64_t Sync = events->Draw("fEvent:TrackID", "boundaryType == 2", "goff");
    
    // Set lookup instead of one huge "(a && b) || ..." string: cling JIT-compiled that string serially and it was
    // tested term by term for every entry.
    std::unordered_set<long long> hitTracks;
    for (Long64_t i = 0; i < Sync; i++)
        hitTracks.insert(((long long)(int)events->GetV1()[i] << 32) | (unsigned int)(int)events->GetV2()[i]);

    auto hasHit = [&](int ev, int id){return hitTracks.count(((long long)ev << 32) | (unsigned int)id) > 0;};
    auto hasHitProton = [&](int ev, int id, int boundary, int pdg){return boundary == 2 && pdg == 2212 && hasHit(ev, id);};

    ROOT::RDataFrame df(*events);
    ROOT::RDF::RSnapshotOptions lazy;
    lazy.fLazy = true;

    auto selectedEvent = df.Filter(hasHit, {"fEvent", "TrackID"});
    auto selectedEvent2 = df.Filter(hasHitProton, {"fEvent", "TrackID", "boundaryType", "pdg"});


    // Lazy snapshots: both cuts of a dataframe are written in one pass over its input and the two dataframes run
    // at the same time. Each tree goes to its own temp file (snapshots must not write into one file at once);
    // the temp files are merged into the final cut files afterwards.
    const std::vector<std::string> cutFiles = {"cut_EOL.root", "cut_EOL2.root"};
    const std::vector<std::string> treeNames = {events->GetName()};
    auto tmpName = [](const std::string& tree, const std::string& cutFile){return "tmp_" + tree + "_" + cutFile;};

    std::vector<ROOT::RDF::RResultHandle> handles;
    handles.emplace_back(selectedEvent.Snapshot(treeNames[0], tmpName(treeNames[0], cutFiles[0]), "", lazy));
    handles.emplace_back(selectedEvent2.Snapshot(treeNames[0], tmpName(treeNames[0], cutFiles[1]), "", lazy));
    ROOT::RDF::RunGraphs(handles);
    ROOT::DisableImplicitMT(); // merging with implicit MT on crashed in TTree::FlushBasketsImpl

    for (const auto& cutFile : cutFiles) {
        {
            TFileMerger merger;
            merger.OutputFile(cutFile.c_str(), "RECREATE");
            for (const auto& tree : treeNames) merger.AddFile(tmpName(tree, cutFile).c_str());
            if (!merger.Merge()) std::cerr << "CutboundaryType: merge into " << cutFile << " failed\n";
        }
        for (const auto& tree : treeNames) std::remove(tmpName(tree, cutFile).c_str());
    }
}
