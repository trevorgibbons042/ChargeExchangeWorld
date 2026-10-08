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

    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int copyNo = touchable->GetCopyNumber();
    G4ThreeVector detectorCenter = touchable->GetTranslation();

    G4VPhysicalVolume *physVol = touchable->GetVolume();
    G4String volumeName = physVol->GetName();

    auto preVolume  = pre->GetPhysicalVolume();
    auto postVolume = post->GetPhysicalVolume();
    G4String preVolumeName  = "";
    G4String postVolumeName = "";
    if (preVolume) {preVolumeName = preVolume->GetName();}
    if (postVolume) {postVolumeName = postVolume->GetName();}

    if (post->GetStepStatus() == fGeomBoundary && postVolumeName == "physTrack" && preVolumeName != "physTrack"){
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

    if (post->GetStepStatus() == fGeomBoundary && preVolumeName == "physTrack" && postVolumeName != "physTrack"){
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
        man->FillNtupleDColumn(9, KE);
        man->FillNtupleDColumn(10, detectorCenter.x());
        man->FillNtupleDColumn(11, detectorCenter.y());
        man->FillNtupleDColumn(12, detectorCenter.z());
        man->AddNtupleRow(0);

        //auto detConstruction = static_cast<const MyDetectorConstruction*>(
        //G4RunManager::GetRunManager() ->GetUserDetectorConstruction());
        //G4ThreeVector detectorPos = detConstruction->GetDetectorPosition();
        
        G4double phi = std::atan2(Pos.y(), Pos.x());
        if (phi < 0) {phi += 2.0 * CLHEP::pi;}

        G4double theta = std::acos(Pos.z()/std::sqrt(Pos.x()*Pos.x() + Pos.y()*Pos.y() + Pos.z()*Pos.z()));

        man->FillNtupleIColumn(2, 0, evt);
        man->FillNtupleIColumn(2, 1, trackID);
        man->FillNtupleIColumn(2, 2, parentID);
        man->FillNtupleIColumn(2, 3, pdg);
        man->FillNtupleDColumn(2, 4, phi);
        man->FillNtupleDColumn(2, 5, theta);
        man->AddNtupleRow(2);

        track->SetTrackStatus(fStopAndKill);
    }

return true;
}