#include "9stepping.hh"
#include "G4VVisManager.hh"
#include "G4Text.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4StepStatus.hh"

#include <cmath>
#include "G4ChargeExchange.hh"
#include "G4HadronicProcess.hh"
#include "G4ChargeExchangeNP.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{
    fEventAction = eventAction;

    fMessenger = new G4GenericMessenger(this, "/zdc/", "ZDC energy resolution");
    fMessenger->DeclareProperty("resStochastic", fResStochastic, "Stochastic term a in a/sqrt(E[GeV])");
    fMessenger->DeclareProperty("resConstant", fResConstant, "Constant term b");
};

MySteppingAction::~MySteppingAction(){
    delete fMessenger;
};

void MySteppingAction::UserSteppingAction(const G4Step *step){    
    ////Go the volume where the step is happening
    //G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    ////Get the detector geometry
    //const MyDetectorConstruction *detectorConstruction = 
    //static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    
    ////Set scoring volume (the volume that records energy data) equal to detector geometry
    //G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();
    //if(volume!=fScoringVolume){return;}; //ignore every step unless it happens inside this chosen volume

    ////Get Energy Deposited in the material during the step
    //G4double edep = step->GetTotalEnergyDeposit();
    //fEventAction->AddEdep(edep);

    G4Track* track = step->GetTrack();
    G4int id = track->GetTrackID();
    G4int ParentID = track->GetParentID();
    G4double pdg = track->GetParticleDefinition()->GetPDGEncoding();

    auto* post = step->GetPostStepPoint();

    //tracker cell entry: the step ends on the cell boundary but starts in the world, so the SD never sees it.
    //checking pre-in-world also skips cell-to-cell crossings inside the tracker
    if (post->GetStepStatus() == fGeomBoundary && (pdg == 2212 || pdg == 2112)){
        const G4VPhysicalVolume* postVolume = post->GetPhysicalVolume();
        const G4VPhysicalVolume* preVolume = step->GetPreStepPoint()->GetPhysicalVolume();
        const G4LogicalVolume* trackLV = static_cast<const MyDetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetTrackVolume();

        if (postVolume && postVolume->GetLogicalVolume() == trackLV &&
            preVolume && preVolume->GetLogicalVolume() != trackLV){
            G4AnalysisManager *man = G4AnalysisManager::Instance();
            G4ThreeVector Pos = post->GetPosition();
            G4ThreeVector detectorCenter = post->GetTouchable()->GetTranslation(); //global center, cells sit inside a slab

            man->FillNtupleIColumn(0, 0, G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
            man->FillNtupleIColumn(0, 1, id);
            man->FillNtupleIColumn(0, 2, ParentID);
            man->FillNtupleIColumn(0, 3, static_cast<G4int>(pdg));
            man->FillNtupleIColumn(0, 4, 0);
            man->FillNtupleDColumn(0, 5, Pos.x());
            man->FillNtupleDColumn(0, 6, Pos.y());
            man->FillNtupleDColumn(0, 7, Pos.z());
            man->FillNtupleDColumn(0, 8, post->GetKineticEnergy());
            man->FillNtupleDColumn(0, 9, detectorCenter.x());
            man->FillNtupleDColumn(0, 10, detectorCenter.y());
            man->FillNtupleDColumn(0, 11, detectorCenter.z());
            man->FillNtupleDColumn(0, 12, post->GetKineticEnergy());
            man->AddNtupleRow(0);
        }
    }

    //cheap process-type check first, dynamic_cast only for hadronic steps
    const G4VProcess* definedProcess = post->GetProcessDefinedStep();
    const G4HadronicProcess* process = nullptr;
    if (definedProcess != nullptr && definedProcess->GetProcessType() == fHadronic){
        process = dynamic_cast<const G4HadronicProcess*>(definedProcess);
    }

    if (process != nullptr && process->GetProcessName() == "neutronChargeExNP"){
        G4AnalysisManager *man = G4AnalysisManager::Instance();
        G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

        const G4ChargeExchange* modelNP = 
        dynamic_cast<const G4ChargeExchange*>
        (process->GetHadronicInteraction());

        
        if (modelNP != nullptr){
            G4ThreeVector reactionPosition = step->GetPostStepPoint()->GetPosition();

            const auto* secondaries = step->GetSecondaryInCurrentStep();
            G4double thetaCE = NAN;

            for (const auto* secondary : *secondaries) {
                if (secondary->GetDefinition()->GetPDGEncoding() != 2212) {continue;}
                thetaCE =(secondary->GetMomentum()).theta();
            }
            
            man->FillNtupleIColumn(1, 0, evt);
            man->FillNtupleIColumn(1, 1, id);
            man->FillNtupleDColumn(1, 2, 1);
            man->FillNtupleDColumn(1, 3, reactionPosition.getX());
            man->FillNtupleDColumn(1, 4, reactionPosition.getY());
            man->FillNtupleDColumn(1, 5, reactionPosition.getZ());
            man->FillNtupleDColumn(1, 6, modelNP->GetLastphinew());
            man->FillNtupleDColumn(1, 7, thetaCE);
            man->FillNtupleDColumn(1, 8, modelNP->GetLastmomentum());
            man->FillNtupleDColumn(1, 9, modelNP->GetLastT());
            man->FillNtupleDColumn(1, 10, modelNP->GetLastAP());
            man->AddNtupleRow(1);
        }
    }

    auto* postVol = post->GetPhysicalVolume();
    if (post->GetStepStatus() == fGeomBoundary && (pdg == 2212 || pdg == 2112)
        && postVol && postVol->GetName() == "physDetector"){
        G4AnalysisManager *man = G4AnalysisManager::Instance();
        G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
        G4ThreeVector Pos = post->GetPosition();
        G4double KE       = post->GetKineticEnergy();
        G4ThreeVector detectorCenter = post->GetTouchable()->GetTranslation(); //global center, cells sit inside a slab

        man->FillNtupleIColumn(0, 0, evt);
        man->FillNtupleIColumn(0, 1, id);
        man->FillNtupleIColumn(0, 2, ParentID);
        man->FillNtupleIColumn(0, 3, pdg);
        man->FillNtupleIColumn(0, 4, 2);
        man->FillNtupleDColumn(0, 5, Pos.x());
        man->FillNtupleDColumn(0, 6, Pos.y());
        man->FillNtupleDColumn(0, 7, Pos.z());
        man->FillNtupleDColumn(0, 8, KE);
        man->FillNtupleDColumn(0, 9, detectorCenter.x());
        man->FillNtupleDColumn(0, 10, detectorCenter.y());
        man->FillNtupleDColumn(0, 11, detectorCenter.z());

        G4double relRes = std::sqrt(fResStochastic*fResStochastic/(KE/GeV) + fResConstant*fResConstant);
        man->FillNtupleDColumn(0, 12, std::max(0., G4RandGauss::shoot(KE, KE*relRes)));
        man->AddNtupleRow(0);

        track->SetTrackStatus(fStopAndKill);
    }

    //labeling
    /*
        G4String label = "T" + std::to_string(id) + " " + track->GetParticleDefinition()->GetParticleName();
        G4ThreeVector vertex = track->GetVertexPosition();
        G4ThreeVector dir = track->GetMomentumDirection();
        G4ThreeVector labelPosition;
        
        if(track->GetParentID() == 0){labelPosition = vertex + dir * -0.2*m;}
        else if(id == 2){labelPosition = vertex + dir * 0.5*m;}
        else{labelPosition = vertex + dir * 0.5*m;}

        G4Text text(label, labelPosition);
        text.SetScreenSize(14);

        G4double charge = track->GetParticleDefinition()->GetPDGCharge();
        
        if (charge > 0.0) {
            G4VisAttributes att(G4Colour(0.0, 0.0, 1.0));
            text.SetVisAttributes(att);}
        else if (charge < 0.0) {
            G4VisAttributes att(G4Colour(1.0, 0.0, 0.0));
            text.SetVisAttributes(att);}
        else{
            G4VisAttributes att(G4Colour(0.0, 1.0, 0.0));
            text.SetVisAttributes(att);}
        
        if (vis && charge > 0.0) {vis->Draw(text);}
        */
};