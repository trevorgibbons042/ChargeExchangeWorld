#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4TrackStatus.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsOrderedFreeVector.hh" //detector efficency stuff

//root stuff
#include "G4AnalysisManager.hh"
#include "G4VNtupleManager.hh"
#include "G4Event.hh"
#include "G4Track.hh"

//one row of the ParticleInfo ntuple (booked in 7run.cc)
//boundaryType: 0 tracker entry, 1 tracker exit, 2 ZDC, 3 track start
inline void FillParticleInfo(const G4Track *track, G4int boundaryType, const G4ThreeVector &pos, G4double KE,
                             const G4ThreeVector &detectorCenter, G4double measuredKE)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->FillNtupleIColumn(0, 0, G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
    man->FillNtupleIColumn(0, 1, track->GetTrackID());
    man->FillNtupleIColumn(0, 2, track->GetParentID());
    man->FillNtupleIColumn(0, 3, track->GetParticleDefinition()->GetPDGEncoding());
    man->FillNtupleIColumn(0, 4, boundaryType);
    man->FillNtupleDColumn(0, 5, pos.x());
    man->FillNtupleDColumn(0, 6, pos.y());
    man->FillNtupleDColumn(0, 7, pos.z());
    man->FillNtupleDColumn(0, 8, KE);
    man->FillNtupleDColumn(0, 9, detectorCenter.x());
    man->FillNtupleDColumn(0, 10, detectorCenter.y());
    man->FillNtupleDColumn(0, 11, detectorCenter.z());
    man->FillNtupleDColumn(0, 12, measuredKE);
    man->AddNtupleRow(0);
}

class MySensitiveDetector : public G4VSensitiveDetector
{
public:
    MySensitiveDetector(G4String);
    ~MySensitiveDetector();
    void SetTrackVolume(const G4LogicalVolume *lv) {trackLV = lv;} //logicTrack, compared by pointer in ProcessHits
private:
    const G4LogicalVolume *trackLV = nullptr;
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *) override;

    //adding in detector efficency
    G4PhysicsOrderedFreeVector *quEff;
};

#endif