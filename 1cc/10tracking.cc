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
    G4AnalysisManager* man = G4AnalysisManager::Instance();

    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4int trackID  = track->GetTrackID();
    G4int parentID = track->GetParentID();

    G4int pdg =track->GetParticleDefinition()->GetPDGEncoding();
    G4ThreeVector Pos = track->GetPosition();
    G4ThreeVector Mom = track->GetMomentum();
    G4double KE = track->GetKineticEnergy();

    man->FillNtupleIColumn(0, evt);
    man->FillNtupleIColumn(1, trackID);
    man->FillNtupleIColumn(2, parentID);
    man->FillNtupleIColumn(4, pdg);
    man->FillNtupleIColumn(5, 3);
    man->FillNtupleDColumn(6, Pos.x());
    man->FillNtupleDColumn(7, Pos.y());
    man->FillNtupleDColumn(8, Pos.z());
    //man->FillNtupleDColumn(9, Mom.x());
    //man->FillNtupleDColumn(10, Mom.y());
    //man->FillNtupleDColumn(11, Mom.z());
    man->FillNtupleDColumn(9, KE);
    man->AddNtupleRow(0);
};