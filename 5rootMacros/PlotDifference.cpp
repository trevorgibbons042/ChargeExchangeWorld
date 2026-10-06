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

void PlotDifference(){
    std::string folder = "pngs_difference";
    std::filesystem::create_directories(folder);
    
    double Graph = 1;
    
    TFile* File = TFile::Open("angle_differences.root", "READ");

    TTree *events = (TTree*)File->Get("AngleDifference");

    if (Graph == 1){
        TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
        gPad->SetRightMargin(0.30);
        
        canvas->Clear();
        events->Draw("deltaTheta");
        events->GetHistogram()->SetTitle("Difference in theta");
        canvas->Update();
        canvas->SaveAs((folder + "/1.png").c_str());
        
        canvas->Clear();
        events->Draw("deltaThetaD");
        events->GetHistogram()->SetTitle("Difference in theta (DETECTOR BLOCK)");
        canvas->Update();
        canvas->SaveAs((folder + "/1_2.png").c_str());
        
        canvas->Clear();
        events->Draw("deltaTheta", "PDG==2212");
        events->GetHistogram()->SetTitle("Difference in theta, PDG==2212");
        canvas->Update();
        canvas->SaveAs((folder + "/2.png").c_str());
        
        canvas->Clear();
        events->Draw("deltaThetaD", "PDG==2212");
        events->GetHistogram()->SetTitle("Difference in theta (DETECTOR BLOCK), PDG==2212");
        canvas->Update();
        canvas->SaveAs((folder + "/2_2.png").c_str());
        
        canvas->Clear();
        events->Draw("deltaPhi");
        events->GetHistogram()->SetTitle("Difference in phi");
        canvas->Update();
        canvas->SaveAs((folder + "/3.png").c_str());
        
        canvas->Clear();
        events->Draw("deltaPhiD");
        events->GetHistogram()->SetTitle("Difference in phi (DETECTOR BLOCK)");
        canvas->Update();
        canvas->SaveAs((folder + "/3_2.png").c_str());
        
        canvas->Clear();
        events->Draw("deltaPhi", "PDG==2212");
        events->GetHistogram()->SetTitle("Difference in phi, PDG==2212");
        canvas->Update();
        canvas->SaveAs((folder + "/4.png").c_str());
        
        canvas->Clear();
        events->Draw("deltaPhiD", "PDG==2212");
        events->GetHistogram()->SetTitle("Difference in phi (DETECTOR BLOCK), PDG==2212");
        canvas->Update();
        canvas->SaveAs((folder + "/4_2.png").c_str());
    }
}
