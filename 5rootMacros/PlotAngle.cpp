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
#include <filesystem>
#include <sstream>
#include <TStyle.h>

void PlotAngle(){
    std::string folder = "pngs_angle";
    std::filesystem::create_directories(folder);
    
    double Graph = 1;
    
    TFile* File = TFile::Open("cut_angles.root", "READ");

    TTree *events = (TTree*)File->Get("Angles");

    if (Graph == 1){
        TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
        gPad->SetRightMargin(0.30);
        
        canvas->Clear();
        events->Draw("Angle", "", "colz");
        events->GetHistogram()->SetTitle("Reconsturced Angle (theta) Distribution");
        events->GetHistogram()->SetMinimum(0);
        events->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder + "/1.png").c_str());
        
        canvas->Clear();
        events->Draw("Angle", "PDG==2212", "colz");
        events->GetHistogram()->SetTitle("Reconsturced Angle (theta) Distribution");
        events->GetHistogram()->SetMinimum(0);
        events->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder + "/2.png").c_str());
    }
}
