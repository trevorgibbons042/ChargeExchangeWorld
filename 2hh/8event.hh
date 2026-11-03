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

    //called once per particle at the ZDC entry (9stepping.cc), every particle counts, charged or not, reconstructed or not
    void AddZdcParticle(G4bool charged, G4double KE){
        fNZdcAll++; fZdcKEAll += KE;
        if (charged){fNZdcCharged++; fZdcKECharged += KE;}
    };

private:
    G4double fEdep;
    G4int fNZdcAll = 0, fNZdcCharged = 0;       //particles that reached the ZDC in this event
    G4double fZdcKEAll = 0., fZdcKECharged = 0.; //their summed kinetic energy at the ZDC entry (MeV)
};

#endif