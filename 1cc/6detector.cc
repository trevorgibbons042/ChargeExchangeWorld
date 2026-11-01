#include "6detector.hh"
#include "2construction.hh"
#include <cmath>

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{};

MySensitiveDetector::~MySensitiveDetector()
{};

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist){
    //only boundary crossings of charged particles are recorded, bail out cheaply otherwise
    G4StepPoint *pre = aStep->GetPreStepPoint();
    G4StepPoint *post = aStep->GetPostStepPoint();
    if (post->GetStepStatus() != fGeomBoundary){return true;}

    G4Track *track = aStep->GetTrack();
    if (track->GetParticleDefinition()->GetPDGCharge() == 0.){return true;}

    G4ThreeVector detectorCenter = pre->GetTouchable()->GetTranslation(); //global center, cells sit inside a slab

    auto preVolume  = pre->GetPhysicalVolume();
    auto postVolume = post->GetPhysicalVolume();
    const bool preInTrack  = preVolume  && preVolume->GetLogicalVolume()  == trackLV;
    const bool postInTrack = postVolume && postVolume->GetLogicalVolume() == trackLV;

    if (postInTrack && !preInTrack){
        FillParticleInfo(track, 0, pre->GetPosition(), pre->GetKineticEnergy(), detectorCenter, pre->GetKineticEnergy());
    }

    //leaving the tracker is not recorded: a tracker hit is the entry only (the entry from the world is also filled in 9stepping.cc)

return true;
}
