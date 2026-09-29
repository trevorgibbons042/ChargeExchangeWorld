
#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <memory>
#include <filesystem>
#include <sstream>
#include <TStyle.h>

void CutboundaryType(){
    ROOT::EnableImplicitMT();

    TFile* File = TFile::Open("cut_CE.root", "READ");
    
    TTree *events = (TTree*)File->Get("ParticleInfo");
    TTree *ChargeExchange = (TTree*)File->Get("ChargeExchange_values");
    
    Long64_t Sync = events->Draw("fEvent:TrackID", "boundaryType == 2", "goff");
    
    TString cut_ID = "0";
    for (Long64_t i = 0; i < Sync; i++)
        cut_ID += Form(
            " || (TrackID==%.0f && fEvent==%.0f)",
            events->GetV2()[i],
            events->GetV1()[i]);
    
    ROOT::RDataFrame df(*events);
    ROOT::RDF::RSnapshotOptions update;
    update.fMode = "UPDATE";
    
    TString cut_ID2 = "(" + cut_ID + ") && boundaryType == 2 && pdg == 2212";
    auto selectedEvent = df.Filter(cut_ID.Data());
    auto selectedEvent2 = df.Filter(cut_ID2.Data());
    
    selectedEvent.Snapshot(events->GetName(),"cut_EOL.root");
    selectedEvent2.Snapshot(events->GetName(),"cut_EOL2.root");
}
