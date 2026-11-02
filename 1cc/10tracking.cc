#include "10tracking.hh"
#include "6detector.hh"

MyTrackingAction::MyTrackingAction()
{};

MyTrackingAction::~MyTrackingAction()
{};


void MyTrackingAction::PreUserTrackingAction(const G4Track* track)
{
    gTrackRows.clear(); //the rows of the previous track were flushed in PostUserTrackingAction, this is only a safety net
    if (track->GetParticleDefinition()->GetPDGCharge() != 0.){
        //no detector at the track start, so the detector center columns are 0
        FillParticleInfo(track, 3, track->GetPosition(), track->GetKineticEnergy(), G4ThreeVector(), track->GetKineticEnergy());
    }
};

void MyTrackingAction::PostUserTrackingAction(const G4Track* track)
{
    FlushParticleInfo(track); //writes the track's ParticleInfo rows if it hit enough tracker planes or is a CE proton
};
