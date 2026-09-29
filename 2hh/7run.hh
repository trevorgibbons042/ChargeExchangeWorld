#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"

//root stuff
#include "G4AnalysisManager.hh"
#include "G4VNtupleManager.hh"
#include "G4GenericMessenger.hh"

class MyRunAction : public G4UserRunAction
{
public:
    MyRunAction();
    ~MyRunAction();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
private:
    G4GenericMessenger *fMessenger;
    G4String outputFolder = "output";
};

#endif