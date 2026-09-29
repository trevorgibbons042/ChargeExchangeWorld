#ifndef ACTION_HH
#define ACTION_HH

#include "G4VUserActionInitialization.hh"

#include "5generator.hh"
#include "7run.hh"
#include "8event.hh"
#include "9stepping.hh"

class MyActionInitialization : public G4VUserActionInitialization
{
public:
    MyActionInitialization();
    ~MyActionInitialization();
    virtual void Build() const;
    virtual void BuildForMaster() const;
};

#endif