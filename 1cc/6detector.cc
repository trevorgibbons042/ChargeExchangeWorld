#include "6detector.hh"
#include "2construction.hh"
#include <cmath>

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{};

MySensitiveDetector::~MySensitiveDetector()
{};

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist){
    //only boundary crossings of protons/neutrons are recorded, bail out cheaply otherwise
    G4StepPoint *pre = aStep->GetPreStepPoint();
    G4StepPoint *post = aStep->GetPostStepPoint();
    if (post->GetStepStatus() != fGeomBoundary){return true;}

    G4Track *track = aStep->GetTrack();
    G4int pdg = track->GetParticleDefinition()->GetPDGEncoding();
    if (pdg != 2212 && pdg != 2112){return true;}

    G4AnalysisManager *man = G4AnalysisManager::Instance();

    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4int trackID = track->GetTrackID();
    G4int parentID = track->GetParentID();

    const G4VTouchable *touchable = pre->GetTouchable();
    G4ThreeVector detectorCenter = touchable->GetTranslation(); //global center, cells sit inside a slab

    auto preVolume  = pre->GetPhysicalVolume();
    auto postVolume = post->GetPhysicalVolume();
    const bool preInTrack  = preVolume  && preVolume->GetLogicalVolume()  == trackLV;
    const bool postInTrack = postVolume && postVolume->GetLogicalVolume() == trackLV;

    if (postInTrack && !preInTrack){
        G4ThreeVector Pos = pre->GetPosition();
        G4ThreeVector Mom = pre->GetMomentum();
        G4double KE       = pre->GetKineticEnergy();

        man->FillNtupleIColumn(0, evt);
        man->FillNtupleIColumn(1, trackID);
        man->FillNtupleIColumn(2, parentID);
        man->FillNtupleIColumn(3, pdg);
        man->FillNtupleIColumn(4, 0);
        man->FillNtupleDColumn(5, Pos.x());
        man->FillNtupleDColumn(6, Pos.y());
        man->FillNtupleDColumn(7, Pos.z());
        man->FillNtupleDColumn(8, KE);
        man->FillNtupleDColumn(9, detectorCenter.x());
        man->FillNtupleDColumn(10, detectorCenter.y());
        man->FillNtupleDColumn(11, detectorCenter.z());
        man->FillNtupleDColumn(12, KE);
        man->AddNtupleRow(0);
    }

    if (preInTrack && !postInTrack){
        G4ThreeVector Pos = post->GetPosition();
        G4ThreeVector Mom = post->GetMomentum();
        G4double KE       = post->GetKineticEnergy();

        man->FillNtupleIColumn(0, evt);
        man->FillNtupleIColumn(1, trackID);
        man->FillNtupleIColumn(2, parentID);
        man->FillNtupleIColumn(3, pdg);
        man->FillNtupleIColumn(4, 1);
        man->FillNtupleDColumn(5, Pos.x());
        man->FillNtupleDColumn(6, Pos.y());
        man->FillNtupleDColumn(7, Pos.z());
        man->FillNtupleDColumn(8, KE);
        man->FillNtupleDColumn(9, detectorCenter.x());
        man->FillNtupleDColumn(10, detectorCenter.y());
        man->FillNtupleDColumn(11, detectorCenter.z());
        man->FillNtupleDColumn(12, KE);
        man->AddNtupleRow(0);
    }

return true;
}
