#include "6detector.hh"
#include "2construction.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{};

MySensitiveDetector::~MySensitiveDetector()
{};

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist){  
    G4cout << "Testing Sensitive Detector" << G4endl;
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4Track *track = aStep->GetTrack();
    G4int trackID = track->GetTrackID();
    G4int parentID = track->GetParentID();
    G4int pdg = track->GetParticleDefinition()->GetPDGEncoding();

    G4StepPoint *pre = aStep->GetPreStepPoint();
    G4StepPoint *post = aStep->GetPostStepPoint();

    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int copyNo = touchable->GetCopyNumber();
    G4VPhysicalVolume *physVol = touchable->GetVolume();
    G4String volumeName = physVol->GetName();

    if (pre->GetStepStatus() == fGeomBoundary && 
    volumeName == "physMat"){
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
        man->FillNtupleDColumn(9, Mom.x());
        man->FillNtupleDColumn(10, Mom.y());
        man->FillNtupleDColumn(11, Mom.z());
        man->FillNtupleDColumn(12, KE);
        man->AddNtupleRow(0);
    }

    if (post->GetStepStatus() == fGeomBoundary && 
    volumeName == "physMat"){
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
        man->FillNtupleDColumn(9, Mom.x());
        man->FillNtupleDColumn(10, Mom.y());
        man->FillNtupleDColumn(11, Mom.z());
        man->FillNtupleDColumn(12, KE);
        man->AddNtupleRow(0);
    }

    if (pre->GetStepStatus() == fGeomBoundary && 
    volumeName == "physDetector"){
        G4ThreeVector Pos = pre->GetPosition();
        G4ThreeVector Mom = pre->GetMomentum();
        G4double KE       = pre->GetKineticEnergy();

        man->FillNtupleIColumn(0, evt);
        man->FillNtupleIColumn(1, trackID);
        man->FillNtupleIColumn(2, parentID);
        man->FillNtupleIColumn(3, copyNo);
        man->FillNtupleIColumn(4, pdg);
        man->FillNtupleIColumn(5, 2);
        man->FillNtupleDColumn(6, Pos.x());
        man->FillNtupleDColumn(7, Pos.y());
        man->FillNtupleDColumn(8, Pos.z());
        man->FillNtupleDColumn(9, Mom.x());
        man->FillNtupleDColumn(10, Mom.y());
        man->FillNtupleDColumn(11, Mom.z());
        man->FillNtupleDColumn(12, KE);
        man->AddNtupleRow(0);

        auto detConstruction = static_cast<const MyDetectorConstruction*>(
        G4RunManager::GetRunManager() ->GetUserDetectorConstruction());
        G4ThreeVector detectorPos = detConstruction->GetDetectorPosition();

        man->FillNtupleIColumn(2, 0, evt);
        man->FillNtupleIColumn(2, 1, trackID);
        man->FillNtupleIColumn(2, 2, parentID);
        man->FillNtupleIColumn(2, 3, pdg);
        man->FillNtupleDColumn(2, 4, std::atan2(Pos.y(), Pos.x()));
        man->FillNtupleDColumn(2, 5, std::atan2(Pos.y()-detectorPos.y(), Pos.x()-detectorPos.x()));
        man->FillNtupleDColumn(
            2, 6, std::atan2(std::sqrt(Pos.x()*Pos.x() + Pos.y()*Pos.y()), Pos.z())
        );
        man->FillNtupleDColumn(
            2, 7, std::atan2(
                std::sqrt(
            (Pos.x()-detectorPos.x())*(Pos.x()-detectorPos.x())
            + (Pos.y()-detectorPos.y())*(Pos.y()-detectorPos.y())), 
            (Pos.z()-detectorPos.z())
            )
        );
        man->AddNtupleRow(2);

        track->SetTrackStatus(fStopAndKill);
    }

return true;
}