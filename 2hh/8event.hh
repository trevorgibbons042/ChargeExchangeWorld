#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Run.hh"

#include "7run.hh"

//root stuff
//#include "g4root.hh"
#include "G4AnalysisManager.hh"
#include "G4VNtupleManager.hh"

class MyEventAction : public G4UserEventAction
{
public:
    MyEventAction(MyRunAction*);
    ~MyEventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    //make fEdep be the total amount of energy deposited over many steps (aka a Event)
    void AddEdep(G4double edep){fEdep += edep;};

private:
    G4double fEdep;
};

#endif