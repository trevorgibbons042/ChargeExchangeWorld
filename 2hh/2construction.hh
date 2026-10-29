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
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4SystemOfUnits.hh" //adds meters and such
#include "G4GenericMessenger.hh" //adds messengers

#include "6detector.hh"

//places cell copyNo (row = copyNo / nCols, col = copyNo % nCols) on an nCols x nRows grid centered in its slab
class GridParameterisation : public G4VPVParameterisation
{
public:
    GridParameterisation(G4int cols, G4double cellX, G4double cellY, G4double totalX, G4double totalY)
        : fCols(cols), fCellX(cellX), fCellY(cellY), fTotalX(totalX), fTotalY(totalY) {}

    virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
    {
        G4int row = copyNo / fCols;
        G4int col = copyNo % fCols;
        physVol->SetTranslation(G4ThreeVector(-fTotalX / 2.0 + fCellX / 2.0 + col * fCellX,
                                              -fTotalY / 2.0 + fCellY / 2.0 + row * fCellY, 0.));
    }

    //every cell has the same size, so the solid is never resized
    virtual void ComputeDimensions(G4Box &, const G4int, const G4VPhysicalVolume *) const {}

private:
    G4int fCols;
    G4double fCellX, fCellY, fTotalX, fTotalY;
};

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();

    virtual G4VPhysicalVolume *Construct();
    G4ThreeVector GetDetectorPosition() const{return detectorPosition;}
    const G4LogicalVolume *GetTrackVolume() const {return logicTrack;} //read per call, geometry rebuilds change the pointer

    //geometry values for settings.csv, names match the /detector/ commands
    std::vector<std::pair<std::string, double>> GetSettings() const
    {
        return {{"nBlocks", nBlocks}, {"nBlocksTracking", nBlocksTracking},
                {"nCols_blocks", nCols}, {"nRows_blocks", nRows},
                {"nCols_detector", nCols_detector}, {"nRows_detector", nRows_detector},
                {"dz", dz}, {"CarbonThickness", CarbonThickness},
                {"randPosDetector", randPosDetector}, {"randPosScaling", randPosScaling}};
    }

    ////Needs to be accessed by 9stepping.cc
    //G4LogicalVolume *GetScoringVolume() const {return fScoringVolume;}

private:
    ////Put variables in here for when you want it to change
    virtual void ConstructSDandField(); //sensiive detector mangetic field
    
    ////messengar stuff
    G4GenericMessenger *fMessenger;
    G4int nCols, nRows, nBlocks, nBlocksTracking;
    G4int nCols_detector, nRows_detector;

    ////Defining World
    G4NistManager *nist;
    G4Material *worldMat, *CMat, *DMat;
    G4MaterialPropertiesTable *mptWorld;
    
    G4Box *solidWorld, *solidMat, *solidDetector, *solidTrack;
    G4LogicalVolume *logicWorld, *logicMat, *logicDetector, *logicTrack;
    G4VPhysicalVolume *physWorld, *physMat, *physDetector, *physTrack;
    G4ThreeVector detectorPosition;

    G4int randPosDetector = -1;
    G4double randPosScaling = 1;
    G4double CarbonThickness, dz;

    ////Detector Defining
    G4LogicalVolume *fScoringVolume;

    void DefineMaterial();
};

#endif
