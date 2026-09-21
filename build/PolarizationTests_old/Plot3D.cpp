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
    double Graph = 3;
    TCanvas* canvas = new TCanvas("canvas", "Plots", 1200, 1200);
    canvas->Clear();
    
    if (Graph == 3){
        canvas->Divide(5, 3);
        for (int i = 0; i <= 14; ++i) {
            std::string filename =
            "/Users/trevorg04/G4ChargeExchange/build/PolarizationTests_old/output"
            + std::to_string(i)
            + ".root";
        TFile* File = TFile::Open(filename.c_str(), "READ");
            
        TTree *events = (TTree*)File->Get("Hits");
        //TTree *events = (TTree*)File->Get("ChargeExchange_values");
        canvas->cd(i+1);
        events->Draw("fY:fX", "fX < 10 && fX > -10 && fY < 10 && fY > -10", "colz");
        delete File;
        }
    canvas->Update();
    }
}
