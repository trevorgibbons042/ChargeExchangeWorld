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
        events->Draw("outTheta", "", "colz");
        events->GetHistogram()->SetTitle("Reconsturced Angle (theta) Distribution");
        canvas->Update();
        canvas->SaveAs((folder + "/1.png").c_str());
        
        canvas->Clear();
        events->Draw("outDTheta", "", "colz");
        events->GetHistogram()->SetTitle("Reconsturced Angle (theta) Distribution (DETECTOR BLOCK)");
        canvas->Update();
        canvas->SaveAs((folder + "/1_2.png").c_str());
        
        canvas->Clear();
        events->Draw("outTheta", "PDG==2212", "colz");
        events->GetHistogram()->SetTitle("Reconsturced Angle (theta) Distribution, PDG == 2212");
        canvas->Update();
        canvas->SaveAs((folder + "/2.png").c_str());
        
        canvas->Clear();
        events->Draw("outDTheta", "PDG==2212", "colz");
        events->GetHistogram()->SetTitle("Reconsturced Angle (theta) Distribution (DETECTOR BLOCK), PDG == 2212");
        canvas->Update();
        canvas->SaveAs((folder + "/2_2.png").c_str());
        
        canvas->Clear();
        events->Draw("outPhi", "", "colz");
        events->GetHistogram()->SetTitle("Reconsturced Angle (phi) Distribution");
        canvas->Update();
        canvas->SaveAs((folder + "/3.png").c_str());
        
        canvas->Clear();
        events->Draw("outDPhi", "", "colz");
        events->GetHistogram()->SetTitle("Reconsturced Angle (phi) Distribution (DETECTOR BLOCK)");
        canvas->Update();
        canvas->SaveAs((folder + "/3_2.png").c_str());
        
        canvas->Clear();
        events->Draw("outPhi", "PDG==2212", "colz");
        events->GetHistogram()->SetTitle("Reconsturced Angle (phi) Distribution, PDG == 2212");
        canvas->Update();
        canvas->SaveAs((folder + "/4.png").c_str());
        
        canvas->Clear();
        events->Draw("outDPhi", "PDG==2212", "colz");
        events->GetHistogram()->SetTitle("Reconsturced Angle (phi) (DETECTOR BLOCK) Distribution, PDG == 2212");
        canvas->Update();
        canvas->SaveAs((folder + "/4_2.png").c_str());
    }
}
