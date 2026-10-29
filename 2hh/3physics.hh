#ifndef PHYSICS_HH
#define PHYSICS_HH

#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"

#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "FTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_ATL.hh"
#include "G4ChargeExchangePhysics.hh"
#include "G4GenericMessenger.hh"

class MyPhysicsList : public G4VModularPhysicsList
{
public:
    MyPhysicsList();
    ~MyPhysicsList();
    //void ConstructParticle();
    static void ApplyAPScale(); //pushes /CE/APScale into the calling thread's NP model, call from each thread at the start of a run
    static G4double GetAPScale() {return APScale;}
private:
    G4GenericMessenger *fCEMessenger;
    static G4bool noFTFNeutronCE; //true turns off neutron quark exchange inside FTF, /CE/noFTFNeutronCE, static like APScale
    static G4double APScale; //multiplies the analyzing power (clamped at -1.0), /CE/APScale, static so worker threads read the master's value
};

#endif