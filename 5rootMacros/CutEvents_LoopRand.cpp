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
    TChain *CalcVals = new TChain("CalcVals");
    
    for (double loopvar = 0; loopvar < 2.1; loopvar = loopvar + .1){
        for (int i = 0; i < 101; i++){
            for (int j = 0; j < 128; j++){
                TString filename = Form("LoopRand_%.1f/output%d_t%d.root",loopvar, i, j);
                events->Add(filename);
                ChargeExchange->Add(filename);
                CalcVals->Add(filename);
            }
        }
        
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
    ROOT::RDataFrame dfCV(*CalcVals);

    ROOT::RDF::RSnapshotOptions update;
    update.fMode = "UPDATE";

    auto selectedEvent = df.Filter(cut_EventSync.Data());
    auto selectedParent = df.Filter(cut_ParentID.Data());

    auto selectedEventCV = dfCV.Filter(cut_EventSync.Data());
    auto selectedParentCV = dfCV.Filter(cut_ParentID.Data());
    
    selectedEvent.Snapshot(events->GetName(),"cut_CE.root");
    dfCE.Snapshot("ChargeExchange_values","cut_CE.root","",update);
    selectedEventCV.Snapshot(CalcVals->GetName(),"cut_CE.root", update);
    
    selectedParent.Snapshot(events->GetName(),"cut_ParentID.root");
    dfCE.Snapshot("ChargeExchange_values","cut_ParentID.root","",update);
    selectedParentCV.Snapshot(CalcVals->GetName(),"cut_ParentID.root", update);
}
