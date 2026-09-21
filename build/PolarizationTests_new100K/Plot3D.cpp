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
    int Graph = 1;
    std::string filename_firstpart = "/Users/trevorg04/G4ChargeExchange/build/PolarizationTests_new100K/output";
    
    TCanvas* canvas = new TCanvas("canvas", "Plots", 1200, 1200);
    canvas->Clear();
    
    if (Graph == 3){
        canvas->Divide(3,2);
        for (int i = 0; i <= 5; ++i) {
            std::string filename =
            filename_firstpart
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
        canvas->Divide(3,2);
        for (int i = 0; i <= 5; ++i) {
            std::string filename =
            filename_firstpart
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
        canvas->Divide(3,2);
        for (int i = 0; i <= 5; ++i) {
            std::string filename =
            filename_firstpart
            + std::to_string(i)
            + ".root";
        TFile* File = TFile::Open(filename.c_str(), "READ");
            
        TTree *events = (TTree*)File->Get("Hits");
        canvas->cd(i+1);
        events->Draw("fY:fX", "fX < 10 && fX > -10 && fY > -10 && fY < 10", "colz");
        delete File;
        }
    canvas->Update();
    }
}
