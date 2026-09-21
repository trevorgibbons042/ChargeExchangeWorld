#ifndef CONSTRUCTION_HH //if cons_hh is not defined
#define CONSTRUCTION_HH //then define it lol

//p much any class (including ones in .cc) u add needs a .hh file LOL, lucky that the classes = .hh file name 
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4SystemOfUnits.hh" //adds meters and such
#include "G4GenericMessenger.hh" //adds messengers

#include "6detector.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();

    virtual G4VPhysicalVolume *Construct();

    ////Needs to be accessed by 9stepping.cc
    //G4LogicalVolume *GetScoringVolume() const {return fScoringVolume;}

private:
    ////Put variables in here for when you want it to change
    virtual void ConstructSDandField(); //sensiive detector mangetic field
    
    ////messengar stuff
    G4GenericMessenger *fMessenger;
    G4int nCols, nRows;

    ////Defining World
    G4NistManager *nist;
    G4Material *worldMat, *CMat;
    G4MaterialPropertiesTable *mptWorld;
    G4Box *solidWorld;
    G4LogicalVolume *logicWorld;
    G4VPhysicalVolume *physWorld;

    ////Detector Defining
    //G4LogicalVolume *fScoringVolume, *logicDetector;

    void DefineMaterial();
};

#endif
