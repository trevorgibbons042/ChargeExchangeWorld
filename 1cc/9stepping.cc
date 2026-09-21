#include "9stepping.hh"
#include "G4VVisManager.hh"
#include "G4Text.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4StepStatus.hh"

#include "G4ChargeExchange.hh"
#include "G4HadronicProcess.hh"
#include "G4ChargeExchangeNP.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{
    fEventAction = eventAction;
};

MySteppingAction::~MySteppingAction(){
};

void MySteppingAction::UserSteppingAction(const G4Step *step){
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    ////Go the volume where the step is happening
    //G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    ////Get the detector geometry
    //const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    
    ////Set scoring volume (the volume that records energy data) equal to detector geometry
    //G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();
    //if(volume!=fScoringVolume){return;}; //ignore every step unless it happens inside this chosen volume

    ////Get Energy Deposited in the material during the step
    //G4double edep = step->GetTotalEnergyDeposit();
    //fEventAction->AddEdep(edep);

    G4Track* track = step->GetTrack();
    G4int id = track->GetTrackID();
    G4double charge = track->GetParticleDefinition()->GetPDGCharge();
    G4double pdg = track->GetParticleDefinition()->GetPDGEncoding();
    auto* vis = G4VVisManager::GetConcreteInstance();
    G4double Kenergy = track->GetKineticEnergy();
    G4double ParentID = track->GetParentID();
    G4VisAttributes textVis;

    // Only the original primary particle
    auto* pre  = step->GetPreStepPoint();
    auto* post = step->GetPostStepPoint();
    
    const G4HadronicProcess* process =
    dynamic_cast<const G4HadronicProcess*>
    (step->GetPostStepPoint()->GetProcessDefinedStep());

    //World Boundary Values
    if (post->GetStepStatus() == fWorldBoundary){
        G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
        G4ThreeVector posTrue = post->GetPosition();

        man->FillNtupleIColumn(0, 0, evt); //for "Hits"
        man->FillNtupleDColumn(0, 1, posTrue[0]); //for fX
        man->FillNtupleDColumn(0, 2, posTrue[1]); //for fY
        man->FillNtupleDColumn(0, 3, posTrue[2]); //for fZ
        man->FillNtupleIColumn(0, 4, pdg);
        man->FillNtupleIColumn(0, 5, Kenergy);
        man->FillNtupleIColumn(0, 6, id);
        man->FillNtupleIColumn(0, 7, ParentID);
        man->AddNtupleRow(0);
    }

    //Charge Exchange Process Values
    if(process != nullptr){
        const G4ChargeExchange* modelNP = 
        dynamic_cast<const G4ChargeExchange*>
        (process->GetHadronicInteraction());
    }
        
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    
    if (process != nullptr && process->GetProcessName() == "neutronChargeExNP"){

        const G4ChargeExchange* modelNP = 
        dynamic_cast<const G4ChargeExchange*>
        (process->GetHadronicInteraction());

        if (modelNP != nullptr){
            G4ThreeVector reactionPosition = step->GetPostStepPoint()->GetPosition();
            man->FillNtupleIColumn(1, 0, evt);
            man->FillNtupleDColumn(1, 1, modelNP->GetLastphinew());
            man->FillNtupleDColumn(1, 2, modelNP->GetLasttheta());
            man->FillNtupleDColumn(1, 3, modelNP->GetLastmomentum());
            man->FillNtupleDColumn(1, 4, reactionPosition.getX());
            man->FillNtupleDColumn(1, 5, reactionPosition.getY());
            man->FillNtupleDColumn(1, 6, reactionPosition.getZ());
            man->FillNtupleDColumn(1, 7, 1);
            man->FillNtupleDColumn(1, 8, id);
            man->AddNtupleRow(1);
        }
    }


    /*
    if (process->GetProcessName() == "neutronChargeExNP" && track->GetDefinition()->GetPDGEncoding() == 2112) {
        G4cout << "Killing incoming primary after momentum change"
               << " TrackID=" << track->GetTrackID()
               << G4endl;
        track->SetTrackStatus(fStopAndKill);
    }
    else if (process == nullptr) {return;}
    */

    /*
    if (process->GetProcessName() == "hadElastic" && deltaP > 0) {
        G4cout << "Killing incoming primary after elastic"
            << " TrackID=" << track->GetTrackID()
            << G4endl;
        track->SetTrackStatus(fStopAndKill);
    }
    */

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