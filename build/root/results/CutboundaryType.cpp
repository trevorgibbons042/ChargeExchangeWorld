
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
    //TFile* File = TFile::Open("cut_ParentID.root", "READ");
    
    TTree *events = (TTree*)File->Get("ParticleInfo");
    TTree *ChargeExchange = (TTree*)File->Get("ChargeExchange_values");
    TTree *CalcVals = (TTree*)File->Get("CalcVals");
    
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
    
    ROOT::RDataFrame dfCV(*CalcVals);
    auto selectedEventCV = dfCV.Filter(cut_ID.Data());
    auto selectedEventCV2 = df.Filter(cut_ID2.Data());
    
    selectedEvent.Snapshot(events->GetName(),"cut_EOL.root");
    selectedEventCV.Snapshot(CalcVals->GetName(),"cut_EOL.root", "", update);
    selectedEvent2.Snapshot(events->GetName(),"cut_EOL2.root");
    selectedEventCV2.Snapshot(CalcVals->GetName(),"cut_EOL2.root","", update);
}
