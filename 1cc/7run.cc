#include "7run.hh"
#include "3physics.hh"
#include "G4ProcessTable.hh"
#include "G4HadronicProcess.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4ChargeExchangeNP.hh"
#include "G4ChargeExchange.hh"
#include "2construction.hh"
#include "5generator.hh"
#include "9stepping.hh"
#include "G4RunManager.hh"
#include "G4Neutron.hh"
#include "G4ProcessManager.hh"
#include <filesystem>
#include <fstream>
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
    man->CreateNtupleDColumn("LastPhiLab");
    man->CreateNtupleDColumn("LastThetaModel");
    man->CreateNtupleDColumn("nProtonSec");
    man->CreateNtupleDColumn("PinTheta");
    man->CreateNtupleDColumn("ModelTheta");
    man->CreateNtupleDColumn("CMTheta");
    man->CreateNtupleSColumn("CEVolume"); // physical volume of the post-step point of the CE step (physMat carbon, physTrackSlab silicon, physWorld vacuum)
    man->CreateNtupleSColumn("CEPreVolume");
    man->CreateNtupleSColumn("CECreator"); // process that created the CE track ("primary" for the beam neutron)
    man->CreateNtupleDColumn("CEPolX"); // polarization of the CE track in the world frame
    man->CreateNtupleDColumn("CEPolY");
    man->CreateNtupleDColumn("CEPolZ");
    man->CreateNtupleDColumn("CEVertexX"); // where the CE neutron was made (track->GetVertexPosition(), the gun position for the beam neutron)
    man->CreateNtupleDColumn("CEVertexY");
    man->CreateNtupleDColumn("CEVertexZ");
    man->CreateNtupleDColumn("CENDirX"); // direction of the CE neutron going into the CE (pre-step momentum direction, world frame)
    man->CreateNtupleDColumn("CENDirY");
    man->CreateNtupleDColumn("CENDirZ");
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

void MyRunAction::EndOfRunAction(const G4Run* run)
{
    //settings the analysis macros read (5rootMacros/FitCore.h ReadSettings), written once from the master
    if (IsMaster()){
        std::vector<std::pair<std::string, double>> values = {
            {"APScale", MyPhysicsList::GetAPScale()}, {"NPCrossSectionFactor", NPCrossSectionFactor},
            {"TotalNeutrons", run->GetNumberOfEventToBeProcessed()},
            {"PPValue", MyPrimaryGenerator::GetPPValue()}, {"PPSigma", MyPrimaryGenerator::GetPPSigma()},
            {"resStochastic", MySteppingAction::GetResStochastic()}, {"resConstant", MySteppingAction::GetResConstant()}};
        auto* detector = static_cast<const MyDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        for (const auto& value : detector->GetSettings()){values.push_back(value);}
        //incident neutron momentum (MeV/c) and kinetic energy (MeV), mean and standard deviation over the primaries, the analysis macros take the beam from here (FitCore.h)
        G4double meanMomentum, sdMomentum, meanKE, sdKE;
        MyPrimaryGenerator::GetBeamStats(meanMomentum, sdMomentum, meanKE, sdKE);
        values.push_back({"BeamMomentum", meanMomentum}); values.push_back({"BeamMomentumSD", sdMomentum});
        values.push_back({"BeamKE", meanKE}); values.push_back({"BeamKESD", sdKE});
        //the nophysics macros inactivate neutronInelastic with the rest of the hadronic physics
        auto* neutronInelastic = G4ProcessTable::GetProcessTable()->FindProcess("neutronInelastic", "neutron");
        values.push_back({"NoPhysics", neutronInelastic != nullptr && !G4Neutron::Neutron()->GetProcessManager()->GetProcessActivation(neutronInelastic)});

        //a name starting with # still has its built-in default value (keep in sync with the defaults in the classes above)
        const std::map<std::string, double> defaults = {
            {"APScale", 1}, {"NPCrossSectionFactor", 1}, {"PPValue", 0}, {"PPSigma", 0}, {"resStochastic", 0.30}, {"resConstant", 0},
            {"nBlocks", 4}, {"nBlocksTracking", 2}, {"nCols_blocks", 11}, {"nRows_blocks", 11}, {"nCols_detector", 11}, {"nRows_detector", 11},
            {"dz", 0.25}, {"CarbonThickness", 0.2}, {"randPosDetector", -1}, {"randPosScaling", 1}, {"NoPhysics", 0}};

        std::string names, row;
        for (const auto& [name, value] : values){
            auto defaultValue = defaults.find(name);
            names += (names.empty() ? "" : ",") + std::string(defaultValue != defaults.end() && defaultValue->second == value ? "#" : "") + name;
            row += (row.empty() ? "" : ",") + std::to_string(value);
        }
        std::filesystem::create_directories("root/" + outputFolder);
        std::ofstream("root/" + outputFolder + "/settings.csv") << names << "\n" << row << "\n";
    }

    ////Writes into root file
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    if (printStatus){G4cout << "writing: " << G4endl;}
    man->Write();

    if (printStatus){G4cout << "closing: " << G4endl;}
    man->CloseFile();

    if (printStatus){G4cout << "End of Run" << G4endl<< G4endl<< G4endl<< G4endl<< G4endl<< G4endl;}
};

