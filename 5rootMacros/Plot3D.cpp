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

void Plot3D(){
    std::string folder = "pngs";
    std::filesystem::create_directories(folder);
    
    double Graph = 1;
    
    TFile* File = TFile::Open("cut_CE.root", "READ");
    TFile* File2 = TFile::Open("cut_ParentID.root", "READ");

    TTree *events = (TTree*)File->Get("ParticleInfo");
    TTree *ChargeExchange = (TTree*)File->Get("ChargeExchange_values");
    TTree *CalcVals = (TTree*)File->Get("CalcVals");

    TTree *events2 = (TTree*)File2->Get("ParticleInfo");
    TTree *ChargeExchange2 = (TTree*)File2->Get("ChargeExchange_values");
    TTree *CalcVals2 = (TTree*)File2->Get("CalcVals");

    
    if (Graph == 1){
        
        TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
        gPad->SetRightMargin(0.30);
        
        canvas->Clear();
        events->Draw("fY:fX>>h2(100,-500,500,100,-500,500)", "boundaryType == 2", "colz");
        events->GetHistogram()->SetTitle("Position at end of World, Charge Exchange Events, All Particles;X Position (cm);Y Position (cm)");
        events->GetHistogram()->SetMinimum(0);
        events->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder + "/1.png").c_str());
        
        canvas->Clear();
        events->Draw("fY:fX>>h2(100,-500,500,100,-500,500)", "boundaryType == 2 && pdg == 2212", "colz");
        events->GetHistogram()->SetTitle("Position at end of World, Charge Exchange Events, Only Protons;X Position (cm);Y Position (cm)");
        events->GetHistogram()->SetMinimum(0);
        events->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder + "/2.png").c_str());
        
        canvas->Clear();
        events->Draw("fY:fX>>h2(100,-10,10,100,-10,10)", "boundaryType == 2 && pdg == 2212", "colz");
        events->GetHistogram()->SetTitle("Position at end of World, Charge Exchange Events, Only Protons;X Position (cm);Y Position (cm)");
        events->GetHistogram()->SetMinimum(0);
        events->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder + "/2_2.png").c_str());
        
        canvas->Clear();
        ChargeExchange->Draw("fZ>>hX(100,0,2500)", "", "");
        ChargeExchange->GetHistogram()->SetTitle("Z position where Charge Exchange Happens; Z Position (cm); Count");
        canvas->Update();
        canvas->SaveAs((folder + "/3.png").c_str());
        
        canvas->Clear();
        events2->Draw("fEvent>>hX(50,0,1000000)", "boundaryType == 2 && pdg == 2212", "");
        events2->GetHistogram()->SetTitle("Events where Charge Exchange happens, Hits the end of the world, Only Protons, same ParentID; Event Number ; Count");
        canvas->Update();
        canvas->SaveAs((folder + "/4.png").c_str());
        
        canvas->Clear();
        events2->Draw("fY:fX>>h2(100,-500,500,100,-500,500)", "boundaryType == 2 && pdg == 2212", "colz");
        events2->GetHistogram()->SetTitle("Position at end of World, Charge Exchange Events, Only Protons, same ParentID; X Position (cm);Y Position (cm)");
        events2->GetHistogram()->SetMinimum(0);
        events2->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder + "/5.png").c_str());
        
        canvas->Clear();
        events2->Draw("fY:fX>>h2(100,-20,20,100,-20,20)", "boundaryType == 2 && pdg == 2212", "colz");
        events2->GetHistogram()->SetTitle("Position at end of World, Charge Exchange Events, Only Protons, same ParentID; X Position (cm);Y Position (cm)");
        events2->GetHistogram()->SetMinimum(0);
        events2->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder + "/6.png").c_str());
        
        canvas->Clear();
        events->Draw("fKEnergy>>hX(100,0,25000)", "boundaryType == 2", "");
        events->GetHistogram()->SetTitle("Kinetic Energy at end of World, Charge Exchange Events, All Particles;Kinetic Energy (MeV); Count");
        events->GetHistogram()->GetYaxis()->SetRangeUser(0, 830);
        canvas->Update();
        canvas->SaveAs((folder + "/7.png").c_str());
        
        canvas->Clear();
        events->Draw("fKEnergy>>hX(100,0,25000)", "boundaryType == 2 && pdg == 2212", "");
        events->GetHistogram()->SetTitle("Kinetic Energy at end of World, Charge Exchange Events, Only Protons;Kinetic Energy (MeV); Count");
        canvas->Update();
        canvas->SaveAs((folder + "/8.png").c_str());
        
        canvas->Clear();
        events->Draw("fKEnergy>>hX(100,10000,25000)", "boundaryType == 2 && pdg == 2212", "");
        events->GetHistogram()->SetTitle("Kinetic Energy at end of World, Charge Exchange Events, Only Protons;Kinetic Energy (MeV); Count");
        events->GetHistogram()->GetYaxis()->SetRangeUser(0, 830);
        canvas->Update();
        canvas->SaveAs((folder + "/8_2.png").c_str());
        
        canvas->Clear();
        events2->Draw("fKEnergy>>hX(100,10000,25000)", "boundaryType == 2 && pdg == 2212", "");
        events2->GetHistogram()->SetTitle("Kinetic Energy at end of World, Charge Exchange Events, Only Protons, Same ParentID; Kinetic Energy (MeV); Count");
        events2->GetHistogram()->GetYaxis()->SetRangeUser(0, 20000);
        canvas->Update();
        canvas->SaveAs((folder + "/9.png").c_str());
        
        canvas->Clear();
        events2->Draw("fKEnergy>>hX(100,10000,25000)", "boundaryType == 2 && pdg == 2212", "");
        events2->GetHistogram()->SetTitle("Kinetic Energy at end of World, Charge Exchange Events, Only Protons, Same ParentID; Kinetic Energy (MeV); Count");
        canvas->Update();
        canvas->SaveAs((folder + "/9_2.png").c_str());
        
        canvas->Clear();
        ChargeExchange->Draw("LastPhi", "", "");
        ChargeExchange->GetHistogram()->SetTitle("LastPhi of all Charge Exchange");
        canvas->Update();
        canvas->SaveAs((folder + "/10.png").c_str());
        
        canvas->Clear();
        ChargeExchange->Draw("LastTheta", "", "");
        ChargeExchange->GetHistogram()->SetTitle("LastTheta of all Charge Exchange");
        canvas->Update();
        canvas->SaveAs((folder + "/10_2.png").c_str());
        
        canvas->Clear();
        events2->Draw("phiCalcXY0detector0Pos", "pdg == 2212", "");
        events2->GetHistogram()->SetTitle("Phi calculated from x,y=0 at end of world, Charge Exchange Events, Only Protons, Same ParentID;");
        canvas->Update();
        canvas->SaveAs((folder + "/11.png").c_str());

        canvas->Clear();
        events2->Draw("phiCalcXY0detectorRealPos", "pdg == 2212", "");
        events2->GetHistogram()->SetTitle("Phi calculated from detector position at end of world, Charge Exchange Events, Only Protons, Same ParentID;");
        canvas->Update();
        canvas->SaveAs((folder + "/11_2.png").c_str());

        canvas->Clear();
        events2->Draw("thetaCalcXY0detector0Pos", "pdg == 2212", "");
        events2->GetHistogram()->SetTitle("Theta calculated from x,y=0 at end of world, Charge Exchange Events, Only Protons, Same ParentID;");
        canvas->Update();
        canvas->SaveAs((folder + "/12.png").c_str());

        canvas->Clear();
        events2->Draw("thetaCalcXY0detectorRealPos", "pdg == 2212", "");
        events2->GetHistogram()->SetTitle("Theta calculated from detector position at end of world, Charge Exchange Events, Only Protons, Same ParentID;");
        canvas->Update();
        canvas->SaveAs((folder + "/12_2.png").c_str());
    }
}
