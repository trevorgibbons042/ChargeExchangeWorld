#ifndef STEPPING_HH
#define STEPPING_HH

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4GenericMessenger.hh"

#include "2construction.hh"
#include "8event.hh"

//root stuff
//#include "g4root.hh"
#include "G4AnalysisManager.hh"
#include "G4VNtupleManager.hh"

class MySteppingAction : public G4UserSteppingAction
{
public:
    MySteppingAction(MyEventAction* eventAction);
    ~MySteppingAction();

    virtual void UserSteppingAction(const G4Step*);

private:
    MyEventAction *fEventAction;
    std::set<G4int> labeledTracks;

    //ZDC energy resolution: sigma/E = sqrt((a/sqrt(E[GeV]))^2 + b^2)
    G4GenericMessenger *fMessenger;
    G4double fResStochastic = 0.30;
    G4double fResConstant = 0.0;
};

#endif