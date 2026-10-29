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

    G4ThreeVector detectorCenter = pre->GetTouchable()->GetTranslation(); //global center, cells sit inside a slab

    auto preVolume  = pre->GetPhysicalVolume();
    auto postVolume = post->GetPhysicalVolume();
    const bool preInTrack  = preVolume  && preVolume->GetLogicalVolume()  == trackLV;
    const bool postInTrack = postVolume && postVolume->GetLogicalVolume() == trackLV;

    if (postInTrack && !preInTrack){
        FillParticleInfo(track, 0, pre->GetPosition(), pre->GetKineticEnergy(), detectorCenter, pre->GetKineticEnergy());
    }

    if (preInTrack && !postInTrack){
        FillParticleInfo(track, 1, post->GetPosition(), post->GetKineticEnergy(), detectorCenter, post->GetKineticEnergy());
    }

return true;
}
