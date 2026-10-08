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
    auto* vis = G4VVisManager::GetConcreteInstance();
    G4double Kenergy = track->GetKineticEnergy();
    G4VisAttributes textVis;

    // Only the original primary particle
    auto* pre  = step->GetPreStepPoint();
    auto* post = step->GetPostStepPoint();
    
    const G4HadronicProcess* process =
    dynamic_cast<const G4HadronicProcess*>
    (step->GetPostStepPoint()->GetProcessDefinedStep());

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