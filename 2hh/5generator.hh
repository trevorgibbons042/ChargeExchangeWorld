#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleTable.hh"
#include "G4GenericMessenger.hh"

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
    MyPrimaryGenerator();
    ~MyPrimaryGenerator();
    virtual void GeneratePrimaries(G4Event*);
    void PrintPolarization();
    void PrintPPIndex();

private:
    G4ParticleGun *fParticleGun;
    G4GenericMessenger* fMessenger;
    std::vector<G4ThreeVector>  polarizationTable = {
    G4ThreeVector(0,1,0),
    G4ThreeVector(0,-1,0),
    G4ThreeVector(1,0,0),
    G4ThreeVector(-1,0,0),
    G4ThreeVector(0,0,1),
    G4ThreeVector(0,0,-1),

    G4ThreeVector(0.540302305868,0.841470984808,0),
    G4ThreeVector(-0.999135150273,0.0415806624333,0),
    G4ThreeVector(-0.818277111064,-0.574823946533,0),
    G4ThreeVector(-0.772764487556,0.634692875943,0),

    G4ThreeVector(0.2675,0.152,0.95149),

    G4ThreeVector(-0.49672,0.2675,0.82566),
    G4ThreeVector(0.60643,-0.58337,0.5403),
    G4ThreeVector(0.68033,0.5403,-0.49521),

    G4ThreeVector(0.53512,-0.41615,-0.73516),
    G4ThreeVector(-0.72017,0.62201,-0.30733),
    G4ThreeVector(-0.13274,-0.83094,0.5403),

    G4ThreeVector(-0.08132,-0.50905,-0.85689)};

    G4double PPTable[21] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
                        -0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1.0};

    G4int PolarizationIndex = -1;
    G4int PPIndex = -1;
    G4int randPos = -1;
    G4int randPAngle = -1;
};

#endif