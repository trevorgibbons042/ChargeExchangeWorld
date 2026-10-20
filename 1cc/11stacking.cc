#include "11stacking.hh"

#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include <cstdlib>

MyStackingAction::MyStackingAction()
{
    fMessenger = new G4GenericMessenger(this, "/stack/", "Stacking Properties");
    fMessenger->DeclareProperty("killExtras", killExtras, "Kill gamma, e+-, mu+-, neutrino secondaries");
};

MyStackingAction::~MyStackingAction()
{
    delete fMessenger;
};

G4ClassificationOfNewTrack MyStackingAction::ClassifyNewTrack(const G4Track* track)
{
    if (killExtras && track->GetParentID() > 0){
        //kill gamma, e+-, mu+-, neutrinos; hadrons stay so pion-made p/n are kept
        G4int apdg = std::abs(track->GetParticleDefinition()->GetPDGEncoding());
        if (apdg == 22 || apdg == 11 || apdg == 13 || apdg == 12 || apdg == 14 || apdg == 16){return fKill;}
    }
    return fUrgent;
};
