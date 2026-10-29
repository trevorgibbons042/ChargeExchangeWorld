#include "10tracking.hh"
#include "6detector.hh"

MyTrackingAction::MyTrackingAction()
{};

MyTrackingAction::~MyTrackingAction()
{};


void MyTrackingAction::PreUserTrackingAction(const G4Track* track)
{
    G4int pdg =track->GetParticleDefinition()->GetPDGEncoding();

    if (pdg == 2212 || pdg == 2112){
        //no detector at the track start, so the detector center columns are 0
        FillParticleInfo(track, 3, track->GetPosition(), track->GetKineticEnergy(), G4ThreeVector(), track->GetKineticEnergy());
    }
};
