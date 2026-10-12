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
    std::string foldermain = "pngs/pngs_reconangle/";
    std::filesystem::create_directories(foldermain);
    
    std::string folder = foldermain + "CEPROTONS/NOCUT";
    std::filesystem::create_directories(folder);
    
    std::string folder2 = foldermain + "CEPROTONS/KECUTBIG";
    std::filesystem::create_directories(folder2);
    
    std::string folder3 = foldermain + "CEPROTONS/KECUTSMOL";
    std::filesystem::create_directories(folder3);
    
    double Graph = 1;
    
    TFile* File = TFile::Open("cut_APandT.root", "READ");

    TTree *events = (TTree*)File->Get("CE");
    TTree *events2 = (TTree*)File->Get("CalcValsRe");
    TTree *events3 = (TTree*)File->Get("CalcValsReD");
    TTree *events4 = (TTree*)File->Get("CalcValsReDIFF");
    TTree *events5 = (TTree*)File->Get("CalcValsReDDIFF");

    if (Graph == 1){
        TCanvas* canvas = new TCanvas("canvas", "Plots", 1500, 1200);
        gPad->SetRightMargin(0.30);
        
        canvas->Clear();
        events->Draw("thetaCE", "", "");
        events->GetHistogram()->SetTitle("CE theta Distribution");
        canvas->Update();
        canvas->SaveAs((folder + "/1.png").c_str());
        
        canvas->Clear();
        events->Draw("phiCE", "", "");
        events->GetHistogram()->SetTitle("CE phi Distribution");
        canvas->Update();
        canvas->SaveAs((folder + "/2.png").c_str());
        
        canvas->Clear();
        events2->Draw("thetaRe", "", "");
        events2->GetHistogram()->SetTitle("Recon theta Distribution");
        canvas->Update();
        canvas->SaveAs((folder + "/3.png").c_str());
        
        canvas->Clear();
        events2->Draw("phiRe", "", "");
        events2->GetHistogram()->SetTitle("Recon phi Distribution");
        canvas->Update();
        canvas->SaveAs((folder + "/4.png").c_str());
        
        canvas->Clear();
        events3->Draw("thetaReD", "", "");
        events3->GetHistogram()->SetTitle("Recon block theta Distribution");
        canvas->Update();
        canvas->SaveAs((folder + "/5.png").c_str());
        
        canvas->Clear();
        events3->Draw("phiReD", "", "");
        events3->GetHistogram()->SetTitle("Recon block phi Distribution");
        canvas->Update();
        canvas->SaveAs((folder + "/6.png").c_str());
        
        canvas->Clear();
        events4->Draw("thetaReD", "", "");
        events4->GetHistogram()->SetTitle("Difference between CE and Recon Theta");
        events4->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder + "/7.png").c_str());
        
        canvas->Clear();
        events4->Draw("phiReD", "", "");
        events4->GetHistogram()->SetTitle("Difference between CE and Recon Phi");
        events4->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder + "/8.png").c_str());
        
        canvas->Clear();
        events5->Draw("thetaReD", "", "");
        events5->GetHistogram()->SetTitle("Difference between CE and Recon Block Theta");
        events5->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder + "/9.png").c_str());
        
        canvas->Clear();
        events5->Draw("phiReD", "", "");
        events5->GetHistogram()->SetTitle("Difference between CE and Recon Block Phi");
        events5->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder + "/10.png").c_str());
        
        canvas->Clear();
        events->Draw("thetaCE", "finalKE > 10000", "");
        events->GetHistogram()->SetTitle("CE theta Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folder2 + "/1.png").c_str());
        
        canvas->Clear();
        events->Draw("phiCE", "finalKE > 10000", "");
        events->GetHistogram()->SetTitle("CE phi Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folder2 + "/2.png").c_str());
        
        canvas->Clear();
        events2->Draw("thetaRe", "finalKE > 10000", "");
        events2->GetHistogram()->SetTitle("Recon theta Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folder2 + "/3.png").c_str());
        
        canvas->Clear();
        events2->Draw("phiRe", "finalKE > 10000", "");
        events2->GetHistogram()->SetTitle("Recon phi Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folder2 + "/4.png").c_str());
        
        canvas->Clear();
        events3->Draw("thetaReD", "finalKE > 10000", "");
        events3->GetHistogram()->SetTitle("Recon block theta Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folder2 + "/5.png").c_str());
        
        canvas->Clear();
        events3->Draw("phiReD", "finalKE > 10000", "");
        events3->GetHistogram()->SetTitle("Recon block phi Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folder2 + "/6.png").c_str());
        
        canvas->Clear();
        events4->Draw("thetaReD", "finalKE > 10000", "");
        events4->GetHistogram()->SetTitle("Difference between CE and Recon Theta, finalKE > 10000");
        events4->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder2 + "/7.png").c_str());
        
        canvas->Clear();
        events4->Draw("phiReD", "finalKE > 10000", "");
        events4->GetHistogram()->SetTitle("Difference between CE and Recon Phi, finalKE > 10000");
        events4->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder2 + "/8.png").c_str());
        
        canvas->Clear();
        events5->Draw("thetaReD", "finalKE > 10000", "");
        events5->GetHistogram()->SetTitle("Difference between CE and Recon Block Theta, finalKE > 10000");
        events5->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder2 + "/9.png").c_str());
        
        canvas->Clear();
        events5->Draw("phiReD", "finalKE > 10000", "");
        events5->GetHistogram()->SetTitle("Difference between CE and Recon Block Phi, finalKE > 10000");
        events5->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder2 + "/10.png").c_str());
        
        
        canvas->Clear();
        events->Draw("thetaCE", "finalKE < 10000", "");
        events->GetHistogram()->SetTitle("CE theta Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folder3 + "/1.png").c_str());
        
        canvas->Clear();
        events->Draw("phiCE", "finalKE < 10000", "");
        events->GetHistogram()->SetTitle("CE phi Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folder3 + "/2.png").c_str());
        
        canvas->Clear();
        events2->Draw("thetaRe", "finalKE < 10000", "");
        events2->GetHistogram()->SetTitle("Recon theta Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folder3 + "/3.png").c_str());
        
        canvas->Clear();
        events2->Draw("phiRe", "finalKE < 10000", "");
        events2->GetHistogram()->SetTitle("Recon phi Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folder3 + "/4.png").c_str());
        
        canvas->Clear();
        events3->Draw("thetaReD", "finalKE < 10000", "");
        events3->GetHistogram()->SetTitle("Recon block theta Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folder3 + "/5.png").c_str());
        
        canvas->Clear();
        events3->Draw("phiReD", "finalKE < 10000", "");
        events3->GetHistogram()->SetTitle("Recon block phi Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folder3 + "/6.png").c_str());
        
        canvas->Clear();
        events4->Draw("thetaReD", "finalKE < 10000", "");
        events4->GetHistogram()->SetTitle("Difference between CE and Recon Theta, finalKE < 10000");
        events4->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder3 + "/7.png").c_str());
        
        canvas->Clear();
        events4->Draw("phiReD", "finalKE < 10000", "");
        events4->GetHistogram()->SetTitle("Difference between CE and Recon Phi, finalKE < 10000");
        events4->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder3 + "/8.png").c_str());
        
        canvas->Clear();
        events5->Draw("thetaReD", "finalKE < 10000", "");
        events5->GetHistogram()->SetTitle("Difference between CE and Recon Block Theta, finalKE < 10000");
        events5->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder3 + "/9.png").c_str());
        
        canvas->Clear();
        events5->Draw("phiReD", "finalKE < 10000", "");
        events5->GetHistogram()->SetTitle("Difference between CE and Recon Block Phi, finalKE < 10000");
        events5->GetHistogram()->Fit("gaus");
        canvas->Update();
        canvas->SaveAs((folder3 + "/10.png").c_str());

        //All protons and CE protons, theta/phi and block theta/phi
        TFile* FileAll = TFile::Open("cut_angles_all.root", "READ");
        TTree *all = (TTree*)FileAll->Get("Angles");
        TTree *allCE = (TTree*)FileAll->Get("AnglesCE");

        std::string folderALLPROTONSNOCUT = foldermain + "ANGLES/ALLPROTONS/NOCUT";
        std::filesystem::create_directories(folderALLPROTONSNOCUT);
        std::string folderALLPROTONSKECUTBIG = foldermain + "ANGLES/ALLPROTONS/KECUTBIG";
        std::filesystem::create_directories(folderALLPROTONSKECUTBIG);
        std::string folderALLPROTONSKECUTSMOL = foldermain + "ANGLES/ALLPROTONS/KECUTSMOL";
        std::filesystem::create_directories(folderALLPROTONSKECUTSMOL);
        std::string folderCEPROTONS_FROMCENOCUT = foldermain + "ANGLES/CEPROTONS_FROMCE/NOCUT";
        std::filesystem::create_directories(folderCEPROTONS_FROMCENOCUT);
        std::string folderCEPROTONS_FROMCEKECUTBIG = foldermain + "ANGLES/CEPROTONS_FROMCE/KECUTBIG";
        std::filesystem::create_directories(folderCEPROTONS_FROMCEKECUTBIG);
        std::string folderCEPROTONS_FROMCEKECUTSMOL = foldermain + "ANGLES/CEPROTONS_FROMCE/KECUTSMOL";
        std::filesystem::create_directories(folderCEPROTONS_FROMCEKECUTSMOL);
        std::string folderCEEVENTPROTONSNOCUT = foldermain + "ANGLES/CEEVENTPROTONS/NOCUT";
        std::filesystem::create_directories(folderCEEVENTPROTONSNOCUT);
        std::string folderCEEVENTPROTONSKECUTBIG = foldermain + "ANGLES/CEEVENTPROTONS/KECUTBIG";
        std::filesystem::create_directories(folderCEEVENTPROTONSKECUTBIG);
        std::string folderCEEVENTPROTONSKECUTSMOL = foldermain + "ANGLES/CEEVENTPROTONS/KECUTSMOL";
        std::filesystem::create_directories(folderCEEVENTPROTONSKECUTSMOL);

        canvas->Clear();
        all->Draw("outTheta*180/TMath::Pi()", "", "");
        all->GetHistogram()->SetTitle("All Protons Recon theta Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSNOCUT + "/1.png").c_str());
        
        canvas->Clear();
        all->Draw("outPhi*180/TMath::Pi()", "", "");
        all->GetHistogram()->SetTitle("All Protons Recon phi Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSNOCUT + "/2.png").c_str());
        
        canvas->Clear();
        all->Draw("outDTheta*180/TMath::Pi()", "", "");
        all->GetHistogram()->SetTitle("All Protons Recon block theta Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSNOCUT + "/3.png").c_str());
        
        canvas->Clear();
        all->Draw("outDPhi*180/TMath::Pi()", "", "");
        all->GetHistogram()->SetTitle("All Protons Recon block phi Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSNOCUT + "/4.png").c_str());
        
        canvas->Clear();
        all->Draw("outTheta*180/TMath::Pi()", "finalKE > 10000", "");
        all->GetHistogram()->SetTitle("All Protons Recon theta Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSKECUTBIG + "/1.png").c_str());
        
        canvas->Clear();
        all->Draw("outPhi*180/TMath::Pi()", "finalKE > 10000", "");
        all->GetHistogram()->SetTitle("All Protons Recon phi Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSKECUTBIG + "/2.png").c_str());
        
        canvas->Clear();
        all->Draw("outDTheta*180/TMath::Pi()", "finalKE > 10000", "");
        all->GetHistogram()->SetTitle("All Protons Recon block theta Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSKECUTBIG + "/3.png").c_str());
        
        canvas->Clear();
        all->Draw("outDPhi*180/TMath::Pi()", "finalKE > 10000", "");
        all->GetHistogram()->SetTitle("All Protons Recon block phi Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSKECUTBIG + "/4.png").c_str());
        
        canvas->Clear();
        all->Draw("outTheta*180/TMath::Pi()", "finalKE < 10000", "");
        all->GetHistogram()->SetTitle("All Protons Recon theta Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSKECUTSMOL + "/1.png").c_str());
        
        canvas->Clear();
        all->Draw("outPhi*180/TMath::Pi()", "finalKE < 10000", "");
        all->GetHistogram()->SetTitle("All Protons Recon phi Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSKECUTSMOL + "/2.png").c_str());
        
        canvas->Clear();
        all->Draw("outDTheta*180/TMath::Pi()", "finalKE < 10000", "");
        all->GetHistogram()->SetTitle("All Protons Recon block theta Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSKECUTSMOL + "/3.png").c_str());
        
        canvas->Clear();
        all->Draw("outDPhi*180/TMath::Pi()", "finalKE < 10000", "");
        all->GetHistogram()->SetTitle("All Protons Recon block phi Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderALLPROTONSKECUTSMOL + "/4.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outTheta*180/TMath::Pi()", "isFromCE==1", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon theta Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCENOCUT + "/1.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outPhi*180/TMath::Pi()", "isFromCE==1", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon phi Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCENOCUT + "/2.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDTheta*180/TMath::Pi()", "isFromCE==1", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon block theta Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCENOCUT + "/3.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDPhi*180/TMath::Pi()", "isFromCE==1", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon block phi Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCENOCUT + "/4.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outTheta*180/TMath::Pi()", "isFromCE==1 && finalKE > 10000", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon theta Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCEKECUTBIG + "/1.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outPhi*180/TMath::Pi()", "isFromCE==1 && finalKE > 10000", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon phi Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCEKECUTBIG + "/2.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDTheta*180/TMath::Pi()", "isFromCE==1 && finalKE > 10000", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon block theta Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCEKECUTBIG + "/3.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDPhi*180/TMath::Pi()", "isFromCE==1 && finalKE > 10000", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon block phi Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCEKECUTBIG + "/4.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outTheta*180/TMath::Pi()", "isFromCE==1 && finalKE < 10000", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon theta Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCEKECUTSMOL + "/1.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outPhi*180/TMath::Pi()", "isFromCE==1 && finalKE < 10000", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon phi Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCEKECUTSMOL + "/2.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDTheta*180/TMath::Pi()", "isFromCE==1 && finalKE < 10000", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon block theta Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCEKECUTSMOL + "/3.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDPhi*180/TMath::Pi()", "isFromCE==1 && finalKE < 10000", "");
        allCE->GetHistogram()->SetTitle("CE Proton Recon block phi Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderCEPROTONS_FROMCEKECUTSMOL + "/4.png").c_str());

        //Protons in CE events, may or may not come from the CE process
        canvas->Clear();
        allCE->Draw("outTheta*180/TMath::Pi()", "", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon theta Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSNOCUT + "/1.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outPhi*180/TMath::Pi()", "", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon phi Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSNOCUT + "/2.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDTheta*180/TMath::Pi()", "", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon block theta Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSNOCUT + "/3.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDPhi*180/TMath::Pi()", "", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon block phi Distribution, no KE cut");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSNOCUT + "/4.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outTheta*180/TMath::Pi()", "finalKE > 10000", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon theta Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSKECUTBIG + "/1.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outPhi*180/TMath::Pi()", "finalKE > 10000", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon phi Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSKECUTBIG + "/2.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDTheta*180/TMath::Pi()", "finalKE > 10000", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon block theta Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSKECUTBIG + "/3.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDPhi*180/TMath::Pi()", "finalKE > 10000", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon block phi Distribution, finalKE > 10000");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSKECUTBIG + "/4.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outTheta*180/TMath::Pi()", "finalKE < 10000", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon theta Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSKECUTSMOL + "/1.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outPhi*180/TMath::Pi()", "finalKE < 10000", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon phi Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSKECUTSMOL + "/2.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDTheta*180/TMath::Pi()", "finalKE < 10000", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon block theta Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSKECUTSMOL + "/3.png").c_str());
        
        canvas->Clear();
        allCE->Draw("outDPhi*180/TMath::Pi()", "finalKE < 10000", "");
        allCE->GetHistogram()->SetTitle("CE Event Proton Recon block phi Distribution, finalKE < 10000");
        canvas->Update();
        canvas->SaveAs((folderCEEVENTPROTONSKECUTSMOL + "/4.png").c_str());
        

    }
}
