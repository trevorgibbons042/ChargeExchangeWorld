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
    
    std::string folder_pos = folder + "/pos/";
    std::string folder_posD = folder + "/pos/posD/";
    
    std::string folder_energy = folder + "/energy/";
    std::string folder_energydiff = folder + "/energy/difference";
    
    std::string folder_angles = folder + "/angle/";
    
    std::filesystem::create_directories(folder);
    std::filesystem::create_directories(folder_pos);
    std::filesystem::create_directories(folder_posD);
    std::filesystem::create_directories(folder_energy);
    std::filesystem::create_directories(folder_energydiff);
    std::filesystem::create_directories(folder_angles);
    
    double Graph = 1;
    
    TFile* File = TFile::Open("cut_CE.root", "READ");
    TFile* File2 = TFile::Open("cut_ParentID.root", "READ");
    TFile* FileFake = TFile::Open("cut_FakeID.root", "READ");

    TTree *events = (TTree*)File->Get("ParticleInfo");
    TTree *ChargeExchange = (TTree*)File->Get("ChargeExchange_values");
    TTree *CalcVals = (TTree*)File->Get("CalcVals");

    TTree *events2 = (TTree*)File2->Get("ParticleInfo");
    TTree *ChargeExchange2 = (TTree*)File2->Get("ChargeExchange_values");
    TTree *CalcVals2 = (TTree*)File2->Get("CalcVals");
    
    TTree *eventsFake = (TTree*)FileFake->Get("ParticleInfo");
    TTree *ChargeExchangeFake = (TTree*)File2->Get("ChargeExchange_values");
    TTree *CalcValsFake = (TTree*)FileFake->Get("CalcVals");

    
    if (Graph == 1){
        
        TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
        gPad->SetRightMargin(0.30);
        
        canvas->Clear();
        events->Draw("fY:fX>>h2(100,-500,500,100,-500,500)", "boundaryType == 2", "colz");
        events->GetHistogram()->SetTitle("Position at end of World, Charge Exchange Events, All Particles;X Position (cm);Y Position (cm)");
        events->GetHistogram()->SetMinimum(0);
        events->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder_pos + "1.png").c_str());
        
        canvas->Clear();
        events->Draw("fY:fX>>h2(100,-500,500,100,-500,500)", "boundaryType == 2 && pdg == 2212", "colz");
        events->GetHistogram()->SetTitle("Position at end of World, Charge Exchange Events, Only Protons;X Position (cm);Y Position (cm)");
        events->GetHistogram()->SetMinimum(0);
        events->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder_pos + "2.png").c_str());
        
        canvas->Clear();
        events->Draw("fY:fX>>h2(100,-10,10,100,-10,10)", "boundaryType == 2 && pdg == 2212", "colz");
        events->GetHistogram()->SetTitle("Position at end of World, Charge Exchange Events, Only Protons;X Position (cm);Y Position (cm)");
        events->GetHistogram()->SetMinimum(0);
        events->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder_pos + "2_2.png").c_str());
        
        canvas->Clear();
        events2->Draw("fY:fX>>h2(100,-500,500,100,-500,500)", "boundaryType == 2 && pdg == 2212", "colz");
        events2->GetHistogram()->SetTitle("Position at end of World, Charge Exchange Events, Only Protons, same ParentID; X Position (cm);Y Position (cm)");
        events2->GetHistogram()->SetMinimum(0);
        events2->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder_pos + "3.png").c_str());
        
        canvas->Clear();
        events2->Draw("fY:fX>>h2(100,-20,20,100,-20,20)", "boundaryType == 2 && pdg == 2212", "colz");
        events2->GetHistogram()->SetTitle("Position at end of World, Charge Exchange Events, Only Protons, same ParentID; X Position (cm);Y Position (cm)");
        events2->GetHistogram()->SetMinimum(0);
        events2->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder_pos + "3_2.png").c_str());
        
        
        
        
        
        canvas->Clear();
        events->Draw("fDy:fDx>>h2(100,-500,500,100,-500,500)", "boundaryType == 2", "colz");
        events->GetHistogram()->SetTitle("Position at end of World (DETECTOR BLOCK), Charge Exchange Events, All Particles;X Position (cm);Y Position (cm)");
        events->GetHistogram()->SetMinimum(0);
        events->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder_posD + "1.png").c_str());
        
        canvas->Clear();
        events->Draw("fDy:fDx>>h2(100,-500,500,100,-500,500)", "boundaryType == 2 && pdg == 2212", "colz");
        events->GetHistogram()->SetTitle("Position at end of World (DETECTOR BLOCK), Charge Exchange Events, Only Protons;X Position (cm);Y Position (cm)");
        events->GetHistogram()->SetMinimum(0);
        events->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder_posD + "2.png").c_str());
        
        canvas->Clear();
        events->Draw("fDy:fDx>>h2(100,-10,10,100,-10,10)", "boundaryType == 2 && pdg == 2212", "colz");
        events->GetHistogram()->SetTitle("Position at end of World (DETECTOR BLOCK), Charge Exchange Events, Only Protons;X Position (cm);Y Position (cm)");
        events->GetHistogram()->SetMinimum(0);
        events->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder_posD + "2_2.png").c_str());
        
        canvas->Clear();
        events2->Draw("fDy:fDx>>h2(100,-500,500,100,-500,500)", "boundaryType == 2 && pdg == 2212", "colz");
        events2->GetHistogram()->SetTitle("Position at end of World (DETECTOR BLOCK), Charge Exchange Events, Only Protons, same ParentID; X Position (cm);Y Position (cm)");
        events2->GetHistogram()->SetMinimum(0);
        events2->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder_posD + "3.png").c_str());
        
        canvas->Clear();
        events2->Draw("fDy:fDx>>h2(100,-20,20,100,-20,20)", "boundaryType == 2 && pdg == 2212", "colz");
        events2->GetHistogram()->SetTitle("Position at end of World (DETECTOR BLOCK), Charge Exchange Events, Only Protons, same ParentID; X Position (cm);Y Position (cm)");
        events2->GetHistogram()->SetMinimum(0);
        events2->GetHistogram()->SetMaximum(100);
        canvas->Update();
        canvas->SaveAs((folder_posD + "3_2.png").c_str());
        
        
        
        
        
        
        
        
        canvas->Clear();
        ChargeExchange->Draw("fZ", "", "");
        ChargeExchange->GetHistogram()->SetTitle("Z position where Charge Exchange Happens; Z Position (cm); Count");
        canvas->Update();
        canvas->SaveAs((folder_pos + "4.png").c_str());
        
        canvas->Clear();
        events2->Draw("fEvent", "boundaryType == 2 && pdg == 2212", "");
        events2->GetHistogram()->SetTitle("Events where Charge Exchange happens, Hits the end of the world, Only Protons, same ParentID; Event Number ; Count");
        canvas->Update();
        canvas->SaveAs((folder_pos + "4_2.png").c_str());
        
        
        
        
        canvas->Clear();
        events->Draw("fKEnergy", "boundaryType == 2", "");
        events->GetHistogram()->SetTitle("Kinetic Energy at end of World, Charge Exchange Events, All Particles;Kinetic Energy (MeV); Count");
        events->GetHistogram()->GetYaxis()->SetRangeUser(0, 830);
        canvas->Update();
        canvas->SaveAs((folder_energy + "1.png").c_str());
        
        canvas->Clear();
        events->Draw("fKEnergy", "boundaryType == 2 && pdg == 2212", "");
        events->GetHistogram()->SetTitle("Kinetic Energy at end of World, Charge Exchange Events, Only Protons;Kinetic Energy (MeV); Count");
        canvas->Update();
        canvas->SaveAs((folder_energy + "2.png").c_str());
        
        canvas->Clear();
        events->Draw("fKEnergy", "boundaryType == 2 && pdg == 2212", "");
        events->GetHistogram()->SetTitle("Kinetic Energy at end of World, Charge Exchange Events, Only Protons;Kinetic Energy (MeV); Count");
        canvas->Update();
        canvas->SaveAs((folder_energy + "2_2.png").c_str());
        
        canvas->Clear();
        events2->Draw("fKEnergy", "boundaryType == 2 && pdg == 2212", "");
        events2->GetHistogram()->SetTitle("Kinetic Energy at end of World, Charge Exchange Events, Only Protons, Same ParentID; Kinetic Energy (MeV); Count");
        canvas->Update();
        canvas->SaveAs((folder_energy + "3.png").c_str());
        
        canvas->Clear();
        events2->Draw("fKEnergy", "boundaryType == 2 && pdg == 2212 && fKEnergy < 10000", "");
        events2->GetHistogram()->SetTitle("Real CE protons left out with 10000 MeV cut");
        canvas->Update();
        canvas->SaveAs((folder_energy + "4.png").c_str());
        
        canvas->Clear();
        eventsFake->Draw("fKEnergy", "boundaryType == 2 && pdg == 2212 && fKEnergy > 10000", "");
        eventsFake->GetHistogram()->SetTitle("Fake CE protons kept in with 10000 MeV cut");
        canvas->Update();
        canvas->SaveAs((folder_energy + "4_2.png").c_str());
        
        

        
        
        
        
        
        
        canvas->Clear();
        ChargeExchange->Draw("LastTheta", "", "");
        ChargeExchange->GetHistogram()->SetTitle("LastTheta of all Charge Exchange");
        canvas->Update();
        canvas->SaveAs((folder_angles + "3.png").c_str());
        
        canvas->Clear();
        ChargeExchange->Draw("LastPhi", "", "");
        ChargeExchange->GetHistogram()->SetTitle("LastPhi of all Charge Exchange");
        canvas->Update();
        canvas->SaveAs((folder_angles + "4.png").c_str());
        
        canvas->Clear();
        CalcVals->Draw("thetaCalcXY0detector0Pos", "pdg == 2212", "");
        CalcVals->GetHistogram()->SetTitle("Theta calculated from x,y=0 at end of world, Charge Exchange Events");
        canvas->Update();
        canvas->SaveAs((folder_angles + "1.png").c_str());

        canvas->Clear();
        CalcVals2->Draw("thetaCalcXY0detector0Pos", "pdg == 2212", "");
        CalcVals2->GetHistogram()->SetTitle("Theta calculated from detector position at end of world, Charge Exchange Events, Only Protons, Same ParentID;");
        canvas->Update();
        canvas->SaveAs((folder_angles + "1_2.png").c_str());
        canvas->Clear();
        CalcVals->Draw("phiCalcXY0detector0Pos", "pdg == 2212", "");
        CalcVals->GetHistogram()->SetTitle("Phi calculated from x,y=0 at end of world, Charge Exchange Events");
        canvas->Update();
        canvas->SaveAs((folder_angles + "2.png").c_str());

        canvas->Clear();
        CalcVals2->Draw("phiCalcXY0detector0Pos", "pdg == 2212", "");
        CalcVals2->GetHistogram()->SetTitle("Phi calculated from detector position at end of world, Charge Exchange Events, Only Protons, Same ParentID;");
        canvas->Update();
        canvas->SaveAs((folder_angles + "2_2.png").c_str());
         
    }
}
