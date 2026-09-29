#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>

void CutEvents(){
    ROOT::EnableImplicitMT();
    
    TChain *events = new TChain("ParticleInfo");
    TChain *ChargeExchange = new TChain("ChargeExchange_values");
    
    for (int i = 0; i < 128; i++){
        TString filename = Form("output0_t%d.root",i);
        events->Add(filename);
        ChargeExchange->Add(filename);
    }
    
    Long64_t Sync = ChargeExchange->Draw("fEvent:TrackID", "didChargeExchange>0", "goff");

    TString cut_EventSync = "0";

    for (Long64_t i = 0; i < Sync; i++)
        cut_EventSync += Form(" || fEvent==%.0f", ChargeExchange->GetV1()[i]);

    TString cut_ParentID = "0";

    for (Long64_t i = 0; i < Sync; i++)
        cut_ParentID += Form(
            " || (ParentID==%.0f && fEvent==%.0f)",
            ChargeExchange->GetV2()[i],
            ChargeExchange->GetV1()[i]
        );
    
    ROOT::RDataFrame df(*events);
    ROOT::RDataFrame dfCE(*ChargeExchange);
    ROOT::RDF::RSnapshotOptions update;
    update.fMode = "UPDATE";

    auto selectedEvent = df.Filter(cut_EventSync.Data());
    auto selectedParent = df.Filter(cut_ParentID.Data());
    
    selectedEvent.Snapshot(events->GetName(),"cut_CE.root");
    dfCE.Snapshot("ChargeExchange_values","cut_CE.root","",update);
    dfCE.Snapshot("CalcVals","cut_CE.root","",update);
    
    selectedParent.Snapshot(events->GetName(),"cut_ParentID.root");
    dfCE.Snapshot("ChargeExchange_values","cut_ParentID.root","",update);
    dfCE.Snapshot("CalcVals","cut_ParentID.root","",update);
}

