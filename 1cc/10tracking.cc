#include "10tracking.hh"
#include "6detector.hh"

MyTrackingAction::MyTrackingAction()
{};

MyTrackingAction::~MyTrackingAction()
{};


void MyTrackingAction::PreUserTrackingAction(const G4Track* track)
{
    if (track->GetParticleDefinition()->GetPDGCharge() != 0.){
        //no detector at the track start, so the detector center columns are 0
        FillParticleInfo(track, 3, track->GetPosition(), track->GetKineticEnergy(), G4ThreeVector(), track->GetKineticEnergy());
    }
};
