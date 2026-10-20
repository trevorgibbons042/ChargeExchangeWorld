#include "10tracking.hh"

#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

MyTrackingAction::MyTrackingAction()
{};

MyTrackingAction::~MyTrackingAction()
{};


void MyTrackingAction::PreUserTrackingAction(const G4Track* track)
{
    G4int pdg =track->GetParticleDefinition()->GetPDGEncoding();

    if (pdg == 2212 || pdg == 2112){
        G4AnalysisManager* man = G4AnalysisManager::Instance();

        G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
        G4int trackID  = track->GetTrackID();
        G4int parentID = track->GetParentID();

        G4ThreeVector Pos = track->GetPosition();
        G4ThreeVector Mom = track->GetMomentum();
        G4double KE = track->GetKineticEnergy();

        man->FillNtupleIColumn(0, evt);
        man->FillNtupleIColumn(1, trackID);
        man->FillNtupleIColumn(2, parentID);
        man->FillNtupleIColumn(3, pdg);
        man->FillNtupleIColumn(4, 3);
        man->FillNtupleDColumn(5, Pos.x());
        man->FillNtupleDColumn(6, Pos.y());
        man->FillNtupleDColumn(7, Pos.z());
        //man->FillNtupleDColumn(9, Mom.x());
        //man->FillNtupleDColumn(10, Mom.y());
        //man->FillNtupleDColumn(11, Mom.z());
        man->FillNtupleDColumn(8, KE);
        man->FillNtupleDColumn(12, KE);
        man->AddNtupleRow(0);
    }
};
