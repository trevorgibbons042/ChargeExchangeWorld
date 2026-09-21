#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TROOT.h>
#include <TString.h>

#include <algorithm>
#include <memory>

void Plot3D(){
    int Graph = 6;
    int maxFiles = 21;
    std::string filename_firstpart = "/Users/trevorg04/G4ChargeExchange/build/Polarization_Tests_otherphysics/output";
    
    TCanvas* canvas = new TCanvas("canvas", "Plots", 1200, 1200);
    canvas->Clear();
    canvas->Divide(5,5);
    
    
    if (Graph == 6){
        for (int i = 0; i <= maxFiles; ++i) {
            std::string filename = filename_firstpart
            + std::to_string(i)
            + ".root";
        TFile* File = TFile::Open(filename.c_str(), "READ");
            
        TTree *events = (TTree*)File->Get("ChargeExchange_values");
        canvas->cd(i+1);
        events->Draw("fEvent", "", "");
        delete File;
        }
    canvas->Update();
    }
    
    if (Graph == 5){
        for (int i = 0; i <= maxFiles; ++i) {
            std::string filename = filename_firstpart
            + std::to_string(i)
            + ".root";
        TFile* File = TFile::Open(filename.c_str(), "READ");
            
        TTree *events = (TTree*)File->Get("Hits");
        events->AddFriend("ChargeExchange_values = ChargeExchange_values");
        canvas->cd(i+1);
        events->Draw("fX", "fX > -100 && fX < 100 && fY > -100 && fY < 100 && fZ > 499 && pdg == 2212 && ChargeExchangeYES > 0", "colz");
        delete File;
        }
    canvas->Update();
    }
    
    if (Graph == 4){
        for (int i = 0; i <= maxFiles; ++i) {
            std::string filename = filename_firstpart
            + std::to_string(i)
            + ".root";
        TFile* File = TFile::Open(filename.c_str(), "READ");
            
        TTree *events = (TTree*)File->Get("Hits");
        canvas->cd(i+1);
        events->AddFriend("ChargeExchange_values = ChargeExchange_values");
        events->Draw("fY:fX", "fX > -100 && fX < 100 && fY > -100 && fY < 100 && fZ > 499 && pdg == 2212 && ChargeExchangeYES == 1", "colz");
        delete File;
        }
    canvas->Update();
    }
    
    if (Graph == 3){
        for (int i = 0; i <= maxFiles; ++i) {
            std::string filename = filename_firstpart
            + std::to_string(i)
            + ".root";
        TFile* File = TFile::Open(filename.c_str(), "READ");
            
        TTree *events = (TTree*)File->Get("ChargeExchange_values");
        canvas->cd(i+1);
        events->Draw("theta_init", "theta_init >0 && theta_init < 0.001", "colz");
        delete File;
        }
    canvas->Update();
    }
    
    if (Graph == 2){
        for (int i = 0; i <= maxFiles; ++i) {
            std::string filename = filename_firstpart
            + std::to_string(i)
            + ".root";
        TFile* File = TFile::Open(filename.c_str(), "READ");
            
        TTree *events = (TTree*)File->Get("ChargeExchange_values");
        canvas->cd(i+1);
        events->Draw("phi_new", "", "");
        delete File;
        }
    canvas->Update();
    }

    if (Graph == 1){
        for (int i = 0; i <= maxFiles; ++i) {
            std::string filename = filename_firstpart
            + std::to_string(i)
            + ".root";
        TFile* File = TFile::Open(filename.c_str(), "READ");
            
        TTree *events = (TTree*)File->Get("Hits");
        canvas->cd(i+1);
        events->Draw("fZ:fY:fX", "fZ > 499 && pdg == 2212", "");
        //fX < 100 && fX > -100 && fY < 100 && fY > -100
        //pdg ==
        delete File;
        }
    canvas->Update();
    }
}
