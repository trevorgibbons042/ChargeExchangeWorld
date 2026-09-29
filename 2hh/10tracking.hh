#ifndef TRACKING_HH
#define TRACKING_HH

#include "G4UserTrackingAction.hh"

class G4Track;

class MyTrackingAction : public G4UserTrackingAction
{
public:
    MyTrackingAction();
    ~MyTrackingAction();

    void PreUserTrackingAction(const G4Track* track) override;
};

#endif