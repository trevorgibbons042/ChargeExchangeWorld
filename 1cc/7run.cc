#include "7run.hh"
#include <filesystem>
#include <sstream>

MyRunAction::MyRunAction()
{
    fMessenger = new G4GenericMessenger(this, "/output/", "Output file settings");
    fMessenger->DeclareProperty("folder", outputFolder, "Folder for ROOT output files");
    
    //CREATING EMPTY TABLES
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    
    //6detector.cc
    man->CreateNtuple("ParticleInfo", "Position when Entering/Leaving Volume");
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleIColumn("TrackID");
    man->CreateNtupleIColumn("ParentID");
    man->CreateNtupleIColumn("copyNo");
    man->CreateNtupleIColumn("pdg");
    man->CreateNtupleIColumn("boundaryType");
    man->CreateNtupleDColumn("fX");
    man->CreateNtupleDColumn("fY");
    man->CreateNtupleDColumn("fZ");
    man->CreateNtupleDColumn("fKEnergy");
    man->CreateNtupleDColumn("fDx");
    man->CreateNtupleDColumn("fDy");
    man->CreateNtupleDColumn("fDz");
    man->FinishNtuple(0);

    //9stepping.cc
    man->CreateNtuple("ChargeExchange_values", "ChargeExchange values");
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleDColumn("TrackID");
    man->CreateNtupleDColumn("didChargeExchange");
    man->CreateNtupleDColumn("fX");
    man->CreateNtupleDColumn("fY");
    man->CreateNtupleDColumn("fZ");
    man->CreateNtupleDColumn("LastPhi");
    man->CreateNtupleDColumn("LastTheta");
    man->CreateNtupleDColumn("LastMomentum");
    man->FinishNtuple(1);

    //9stepping.cc
    man->CreateNtuple("CalcVals", "Calculated Values");
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleIColumn("TrackID");
    man->CreateNtupleIColumn("ParentID");
    man->CreateNtupleIColumn("pdg");
    man->CreateNtupleDColumn("phiCalcXY0detector0Pos");
    //man->CreateNtupleDColumn("phiCalcXY0detectorRealPos");
    man->CreateNtupleDColumn("thetaCalcXY0detector0Pos");
    //man->CreateNtupleDColumn("thetaCalcXY0detectorRealPos");
    //man->CreateNtupleDColumn("phiCalcWithMomDirec");
    //man->CreateNtupleDColumn("phiCalcWithMomDirecdetectorRealPos");
    man->FinishNtuple(2);
};

MyRunAction::~MyRunAction()
{delete fMessenger;};

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->Reset();

    G4int runID = run->GetRunID();
    std::stringstream strRunID;
    strRunID << runID;

    std::string folder = "root/" + outputFolder;

    std::filesystem::create_directories(folder);
    std::string fileName = "root/" + outputFolder + "/output" + strRunID.str() + ".root";
    
    G4cout << "Opening: " << fileName << G4endl;
    man->OpenFile(fileName);
    G4cout << "Past open: " << G4endl;
};

void MyRunAction::EndOfRunAction(const G4Run*)
{
    
    ////Writes into root file
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    G4cout << "writing: " << G4endl;
    man->Write();

    G4cout << "closing: " << G4endl;
    man->CloseFile();

    G4cout << "End of Event" << G4endl<< G4endl<< G4endl<< G4endl<< G4endl<< G4endl;
};

