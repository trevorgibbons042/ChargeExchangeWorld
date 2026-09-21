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

void Plot3D(){
    int Graph = 1;
    int maxFiles = 21;
    TCanvas* canvas = new TCanvas("canvas", "Plots", 1200, 1200);
    canvas->Clear();
    canvas->Divide(4,2);
    
    if (Graph == 5){
        for (int i = 0; i <= maxFiles; ++i) {
            std::string filename =
            "/Users/trevorg04/G4ChargeExchange/build/PolarizationTests_percent/output"
            + std::to_string(i)
            + ".root";
        TFile* File = TFile::Open(filename.c_str(), "READ");
            
        TTree *events = (TTree*)File->Get("Hits");
        canvas->cd(i+1);
        events->Draw("fEvent");
        delete File;
        }
    canvas->Update();
    }
    
    if (Graph == 4){
        for (int i = 0; i <= maxFiles; ++i) {
            std::string filename =
            "/Users/trevorg04/G4ChargeExchange/build/PolarizationTests_percent/output"
            + std::to_string(i)
            + ".root";
        TFile* File = TFile::Open(filename.c_str(), "READ");
            
        TTree *events = (TTree*)File->Get("Hits");
        canvas->cd(i+1);
        events->Draw("fX", "fX < 10 && fX > -10 && fY > -1 && fY < 1", "");
        delete File;
        }
    canvas->Update();
    }
    
    if (Graph == 3){
        for (int i = 0; i <= maxFiles; ++i) {
            std::string filename =
            "/Users/trevorg04/G4ChargeExchange/build/PolarizationTests_percent/output"
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
            std::string filename =
            "/Users/trevorg04/G4ChargeExchange/build/PolarizationTests_percent/output"
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
        for (int i = 0; i <= maxFiles/3; ++i) {
            std::string filename =
            "/Users/trevorg04/G4ChargeExchange/build/PolarizationTests_percent/output"
            + std::to_string(3*i)
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
