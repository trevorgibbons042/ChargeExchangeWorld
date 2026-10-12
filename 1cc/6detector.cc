#include "6detector.hh"
#include "2construction.hh"
#include <cmath>

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{};

MySensitiveDetector::~MySensitiveDetector()
{};

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist){  
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4Track *track = aStep->GetTrack();
    G4int trackID = track->GetTrackID();
    G4int parentID = track->GetParentID();
    G4int pdg = track->GetParticleDefinition()->GetPDGEncoding();

    G4StepPoint *pre = aStep->GetPreStepPoint();
    G4StepPoint *post = aStep->GetPostStepPoint();

    const G4VTouchable *touchable = pre->GetTouchable();
    G4VPhysicalVolume *physVol = touchable->GetVolume();
    G4int copyNo = touchable->GetCopyNumber();
    G4ThreeVector detectorCenter = physVol->GetTranslation();

    G4String volumeName = physVol->GetName();
    auto preVolume  = pre->GetPhysicalVolume();
    auto postVolume = post->GetPhysicalVolume();
    G4String preVolumeName  = "";
    G4String postVolumeName = "";
    if (preVolume) {preVolumeName = preVolume->GetName();}
    if (postVolume) {postVolumeName = postVolume->GetName();}

    if (post->GetStepStatus() == fGeomBoundary && postVolumeName == "physTrack" && preVolumeName != "physTrack" 
    && (pdg == 2212 || pdg == 2112)){
        G4ThreeVector Pos = pre->GetPosition();
        G4ThreeVector Mom = pre->GetMomentum();
        G4double KE       = pre->GetKineticEnergy();

        man->FillNtupleIColumn(0, evt);
        man->FillNtupleIColumn(1, trackID);
        man->FillNtupleIColumn(2, parentID);
        man->FillNtupleIColumn(3, copyNo);
        man->FillNtupleIColumn(4, pdg);
        man->FillNtupleIColumn(5, 0);
        man->FillNtupleDColumn(6, Pos.x());
        man->FillNtupleDColumn(7, Pos.y());
        man->FillNtupleDColumn(8, Pos.z());
        man->FillNtupleDColumn(9, KE);
        man->FillNtupleDColumn(10, detectorCenter.x());
        man->FillNtupleDColumn(11, detectorCenter.y());
        man->FillNtupleDColumn(12, detectorCenter.z());
        man->AddNtupleRow(0);
    }

    if (post->GetStepStatus() == fGeomBoundary && preVolumeName == "physTrack" && postVolumeName != "physTrack"
    && (pdg == 2212 || pdg == 2112)){
        G4ThreeVector Pos = post->GetPosition();
        G4ThreeVector Mom = post->GetMomentum();
        G4double KE       = post->GetKineticEnergy();
    
        man->FillNtupleIColumn(0, evt);
        man->FillNtupleIColumn(1, trackID);
        man->FillNtupleIColumn(2, parentID);
        man->FillNtupleIColumn(3, copyNo);
        man->FillNtupleIColumn(4, pdg);
        man->FillNtupleIColumn(5, 1);
        man->FillNtupleDColumn(6, Pos.x());
        man->FillNtupleDColumn(7, Pos.y());
        man->FillNtupleDColumn(8, Pos.z());
        man->FillNtupleDColumn(9, KE);
        man->FillNtupleDColumn(10, detectorCenter.x());
        man->FillNtupleDColumn(11, detectorCenter.y());
        man->FillNtupleDColumn(12, detectorCenter.z());
        man->AddNtupleRow(0);
    }

return true;
}