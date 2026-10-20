#include "7run.hh"
#include "3physics.hh"
#include "G4ProcessTable.hh"
#include "G4HadronicProcess.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4ChargeExchangeNP.hh"
#include "G4ChargeExchange.hh"
#include <filesystem>
#include <sstream>

MyRunAction::MyRunAction()
{
    fMessenger = new G4GenericMessenger(this, "/output/", "Output file settings");
    fMessenger->DeclareProperty("folder", outputFolder, "Folder for ROOT output files");
    fMessenger->DeclareProperty("print", printStatus, "Print run status lines (opening/writing/closing/end of run)");

    fCEMessenger = new G4GenericMessenger(this, "/CE/", "Charge exchange settings");
    fCEMessenger->DeclareProperty("NPCrossSectionFactor", NPCrossSectionFactor, "Factor on the neutronChargeExNP cross section (1 = regular value)");

    //CREATING EMPTY TABLES
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    //man->SetNtupleMerging(true);
    
    //6detector.cc
    man->CreateNtuple("ParticleInfo", "Position when Entering/Leaving Volume");
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleIColumn("TrackID");
    man->CreateNtupleIColumn("ParentID");
    man->CreateNtupleIColumn("pdg");
    man->CreateNtupleIColumn("boundaryType");
    man->CreateNtupleDColumn("fX");
    man->CreateNtupleDColumn("fY");
    man->CreateNtupleDColumn("fZ");
    man->CreateNtupleDColumn("fKEnergy");
    man->CreateNtupleDColumn("fDx");
    man->CreateNtupleDColumn("fDy");
    man->CreateNtupleDColumn("fDz");
    man->CreateNtupleDColumn("fDEnergy");
    man->FinishNtuple(0);

    //9stepping.cc
    man->CreateNtuple("ChargeExchange_values", "ChargeExchange values");
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleIColumn("TrackID");
    man->CreateNtupleDColumn("didChargeExchange");
    man->CreateNtupleDColumn("fX");
    man->CreateNtupleDColumn("fY");
    man->CreateNtupleDColumn("fZ");
    man->CreateNtupleDColumn("LastPhi");
    man->CreateNtupleDColumn("LastTheta");
    man->CreateNtupleDColumn("LastMomentum");
    man->CreateNtupleDColumn("LastT");
    man->CreateNtupleDColumn("LastAP");
    man->FinishNtuple(1);

};

MyRunAction::~MyRunAction()
{delete fMessenger; delete fCEMessenger;};

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->Reset();

    //the CE NP cross section dataset is built per thread inside G4ChargeExchangePhysics, so set its factor here on every thread
    auto* ceProcess = dynamic_cast<G4HadronicProcess*>(G4ProcessTable::GetProcessTable()->FindProcess("neutronChargeExNP", "neutron"));
    if (ceProcess != nullptr){
        for (auto* dataSet : ceProcess->GetCrossSectionDataStore()->GetDataSetList()){
            if (auto* npXS = dynamic_cast<G4ChargeExchangeNP*>(dataSet)){npXS->SetCrossSectionFactor(NPCrossSectionFactor);}
        }
    }
    MyPhysicsList::ApplyAPScale();

    G4int runID = run->GetRunID();
    std::stringstream strRunID;
    strRunID << runID;

    std::string folder = "root/" + outputFolder;

    std::filesystem::create_directories(folder);
    std::string fileName = "root/" + outputFolder + "/output" + strRunID.str() + ".root";
    
    if (printStatus){G4cout << "Opening: " << fileName << G4endl;}
    man->OpenFile(fileName);
    if (printStatus){G4cout << "Past open: " << G4endl;}
};

void MyRunAction::EndOfRunAction(const G4Run*)
{
    
    ////Writes into root file
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    if (printStatus){G4cout << "writing: " << G4endl;}
    man->Write();

    if (printStatus){G4cout << "closing: " << G4endl;}
    man->CloseFile();

    if (printStatus){G4cout << "End of Run" << G4endl<< G4endl<< G4endl<< G4endl<< G4endl<< G4endl;}
};

