#ifndef STACKING_HH
#define STACKING_HH

#include "G4UserStackingAction.hh"
#include "G4GenericMessenger.hh"

class G4Track;

//kills gamma, e+-, mu+-, neutrinos secondaries (hadrons kept), toggle with /stack/killExtras true|false
class MyStackingAction : public G4UserStackingAction
{
public:
    MyStackingAction();
    ~MyStackingAction();

    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track) override;

private:
    G4GenericMessenger *fMessenger;
    G4bool killExtras = true;
};

#endif
