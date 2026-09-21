#include "7run.hh"

MyRunAction::MyRunAction()
{
    //CREATING EMPTY TABLES
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    ////Position of Photons (in detector, with efficiency, for 6detector.cc)
    /*
    man->CreateNtuple("Hits", "Hits"); //name (internal ntuple name), title (descriptive title)
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleDColumn("fX");
    man->CreateNtupleDColumn("fY");
    man->CreateNtupleDColumn("fZ");
    //man->CreateNtupleDColumn("fWlenPhoton");
    man->FinishNtuple(0); //all this only creates the sections, no info yet
    */
    ////Position of Photon (true position, for 6detector.cc)
    //man->CreateNtuple("Photons", "Photons");
    //man->CreateNtupleIColumn("fEvent");
    //man->CreateNtupleDColumn("fX");
    //man->CreateNtupleDColumn("fY");
    //man->CreateNtupleDColumn("fZ");
    //man->CreateNtupleDColumn("fWlenPhoton");
    //man->FinishNtuple(1);

    ////Energy Deposition (for 8event.cc)
    //man->CreateNtuple("Scoring", "Event");
    //man->CreateNtupleDColumn("fEdep");
    //man->FinishNtuple(2);

    man->CreateNtuple("Hits", "Hits");
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleDColumn("fX");
    man->CreateNtupleDColumn("fY");
    man->CreateNtupleDColumn("fZ");
    man->CreateNtupleIColumn("pdg");
    man->CreateNtupleIColumn("fKEnergy");
    man->CreateNtupleIColumn("TrackID");
    man->CreateNtupleIColumn("ParentID");
    man->FinishNtuple(0);

    man->CreateNtuple("ChargeExchange_values", "ChargeExchange values");
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleDColumn("LastPhi");
    man->CreateNtupleDColumn("LastTheta");
    man->CreateNtupleDColumn("LastMomentum");
    man->CreateNtupleDColumn("fX");
    man->CreateNtupleDColumn("fY");
    man->CreateNtupleDColumn("fZ");
    man->CreateNtupleDColumn("didChargeExchange");
    man->CreateNtupleDColumn("TrackID");
    man->FinishNtuple(1);
};

MyRunAction::~MyRunAction()
{};

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->Reset();

    G4int runID = run->GetRunID();
    std::stringstream strRunID;
    strRunID << runID;

    std::string fileName = "output"+ strRunID.str() + ".root"; //puts root into custom folder
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
