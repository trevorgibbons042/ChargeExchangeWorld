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
    G4GenericMessenger *fCEMessenger;
    G4double NPCrossSectionFactor = 1.0; //multiplies the neutronChargeExNP cross section, /CE/NPCrossSectionFactor
    G4String outputFolder = "output";
    G4bool printStatus = true; //opening/writing/closing/end of run lines, /output/print false turns them off
};

#endif