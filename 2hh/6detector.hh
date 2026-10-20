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