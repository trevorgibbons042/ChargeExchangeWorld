#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TROOT.h>

#include <algorithm>
#include <memory>

void Plot3D(const char* filename = "/Users/trevorg04/G4ChargeExchange/build/output0.root"){
    TFile *File = TFile::Open(filename, "READ");
    
    TTree *events = (TTree*)File->Get("Hits");
    TTree *ChargeExchange = (TTree*)File->Get("ChargeExchange_values");
    ChargeExchange->AddFriend(events, "h");
    events->AddFriend(ChargeExchange, "ce");
    
    //gStyle->SetPalette(kRainBow);
    //events->SetMarkerStyle(20);
    //events->SetMarkerSize(1.5);
    //ChargeExchange->SetMarkerStyle(20);
    //ChargeExchange->SetMarkerSize(1.5);
    
    double Graph = 4;
    
    if (Graph == 4){
        TCanvas* canvas = new TCanvas("canvas", "Plots", 1200, 1200);
        canvas->Clear();
        canvas->Divide(3,3);
        
        Long64_t EventSync = ChargeExchange->Draw("fEvent", "didChargeExchange>0", "goff");
        TString cut0 = "0";
        for (Long64_t i = 0; i < EventSync; i++){
            cut0 += Form(" || fEvent==%.0f", ChargeExchange->GetV1()[i]);
        }
        
        TString cut1 = "(" + cut0 + ") && pdg==2212";
        
        Long64_t ParentSync = ChargeExchange->Draw("TrackID", "didChargeExchange>0", "goff");
        TString cut2 = "0";
        for (Long64_t i = 0; i < ParentSync; i++){
            ChargeExchange->Draw("TrackID", "didChargeExchange>0", "goff");
            cut2 += Form(" || (ParentID==%.0f ", ChargeExchange->GetV1()[i]);
            ChargeExchange->Draw("fEvent", "didChargeExchange>0", "goff");
            cut2 += Form("&& fEvent == %.0f) ", ChargeExchange->GetV1()[i]);
        }
        
        TString cut3 = "(" + cut2 + ") && pdg==2212";
        
    
        canvas->cd(1);
        events->Draw("fZ:fY:fX", cut0.Data(), "");
        events->GetHistogram()->SetTitle("Events with Charge Exchange, all Particles");
    
        canvas->cd(2);
        ChargeExchange->Draw("fZ", "", "");
        ChargeExchange->GetHistogram()->SetTitle("Where ChargeExchange happens (neutron)");
        
        canvas->cd(3);
        events->Draw("fKEnergy", cut0.Data(), "");
        events->GetHistogram()->SetTitle("Energy of all Particles in Charge Exchange");
        
        canvas->cd(4);
        events->Draw("fZ:fY:fX", cut3.Data(), "");
        events->GetHistogram()->SetTitle("Events with Charge Exchange, same ParentID and pdg==2212");
        
        canvas->cd(5);
        events->Draw("fY:fX", cut3.Data(), "colz");
        events->GetHistogram()->SetTitle("Events with Charge Exchange, same ParentID and pdg==2212");
        
        canvas->cd(6);
        events->Draw("fKEnergy", cut3.Data(), "");
        events->GetHistogram()->SetTitle("Events with Charge Exchange, same ParentID and pdg==2212");
        
        canvas->cd(7);
        events->Draw("fKEnergy", "pdg == 2212 && fKEnergy > 18800", "");
        events->GetHistogram()->SetTitle("Events with Charge Exchange, pdg==2212");

        
        canvas->Update();
    }
    
    if (Graph == 3){
        TCanvas* canvas = new TCanvas("canvas", "Plots", 1200, 1200);
        canvas->Clear();
        canvas->Divide(1,1);
        
        canvas->cd(1);
        ChargeExchange->Draw("didChargeExchange", "didChargeExchange>0", "");
        
        canvas->Update();
    }
    
    if (Graph == 2){
        TCanvas* canvas = new TCanvas("canvas", "Plots", 1200, 1200);
        canvas->Clear();
        canvas->Divide(3, 3);
        
        canvas->cd(1);
        ChargeExchange->Draw("phi_new", "", "");
        canvas->cd(2);
        events->Draw("fY:fX", "fX > -10 && fX < 10 && fY > -10 && fY < 10 && fZ > 490", "colz");
        canvas->cd(3);
        ChargeExchange->Draw("AnalyzingPower", "", "");
        canvas->cd(4);
        ChargeExchange->Draw("momentumCMS", "", "");
        canvas->cd(5);
        ChargeExchange->Draw("theta_init", "", "");
        canvas->cd(6);
        events->Draw("fZ:fY:fX", "fZ > -510 && fZ < 510", "");
        
        
        canvas->Update();
    }
    
    if (Graph == 1){
        TCanvas* canvas = new TCanvas("canvas", "Plots", 1200, 900);
        canvas->Clear();
        canvas->Divide(5, 4);
        
        double minimum_lv2x = 1000;
        double minimum_fX = 500;
        TString cut = Form("fX < %.17g", minimum_fX);
        TString cut2 = Form("fParentID > %.17g", 1.0);
        TString cut3 = Form("lv2_output_x < %.17g", minimum_lv2x);
        
        canvas->cd(1);
        ChargeExchange->Draw("lastcostT_active", cut.Data(), "");
        canvas->cd(2);
        ChargeExchange->Draw("lastcostT_active", cut2.Data(), "");
        canvas->cd(3);
        ChargeExchange->Draw("lastcostT_active", "", "");
        canvas->cd(4);
        ChargeExchange->Draw("sint_output", cut.Data(), "");
        canvas->cd(5);
        ChargeExchange->Draw("sint_output", cut2.Data(), "");
        canvas->cd(6);
        ChargeExchange->Draw("sint_output", "", "");
        canvas->cd(7);
        ChargeExchange->Draw("lv2_output_x", cut.Data(), "");
        canvas->cd(8);
        ChargeExchange->Draw("lv2_output_x", cut2.Data(), "");
        canvas->cd(9);
        ChargeExchange->Draw("lv2_output_x", "", "");
        canvas->cd(10);
        ChargeExchange->Draw("fParentID", cut.Data(), "");
        canvas->cd(11);
        ChargeExchange->Draw("fParentID", cut2.Data(), "");
        canvas->cd(12);
        ChargeExchange->Draw("fParentID", "", "");
        canvas->cd(13);
        ChargeExchange->Draw("fEvent", cut.Data(), "");
        canvas->cd(14);
        ChargeExchange->Draw("fEvent", cut2.Data(), "");
        canvas->cd(15);
        events->Draw("fX", cut.Data(), "");
        canvas->cd(16);
        events->Draw("fX", cut2.Data(), "");
        canvas->cd(17);
        events->Draw("fX", "", "");
        canvas->cd(18);
        events->Draw("fZ:fY:fX", cut.Data(), "");
        canvas->cd(19);
        events->Draw("fZ:fY:fX", cut2.Data(), "");
        canvas->cd(20);
        events->Draw("fZ:fY:fX", "", "");
        
        canvas->Update();
    }
}
