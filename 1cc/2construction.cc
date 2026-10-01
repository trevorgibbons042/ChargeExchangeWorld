#include "2construction.hh"
#include "Randomize.hh"
#include "G4SDManager.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
    ////Messengar. allows us to make customized gui commands (like change number of rows on detector)    
    fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");
    fMessenger->DeclareProperty("nCols", nCols, "Number of columns");
    fMessenger->DeclareProperty("nRows", nRows, "Number of rows");
    fMessenger->DeclareProperty("nBlocks", nBlocks, "Number of Blocks");
    fMessenger->DeclareProperty("randPosDetector", randPosDetector, "Is detector in random position?");
    fMessenger->DeclareProperty("randPosScaling", randPosScaling, "randPos Scaling Factor");

    //default values
    nBlocks = 20;
    nCols = 9;
    nRows = 9;

    DefineMaterial();
}

MyDetectorConstruction::~MyDetectorConstruction()
{
    delete fMessenger;
}

void MyDetectorConstruction::DefineMaterial(){
    //sets up material manager
    nist = G4NistManager::Instance();

    //world mat and properties (refractive index)
    worldMat = nist->FindOrBuildMaterial("G4_Galactic"); //get air from manager, need a world first before detector
    G4double energyWorld[2] = {1.239841939*eV/0.9, 1.239841939*eV/0.2}; //E/wavelength
    G4double rindexWorld[2] = {1.1, 1.1}; //this is not true cuz of dispersion
    G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
    mptWorld->AddProperty("RINDEX", energyWorld, rindexWorld, 2);
    worldMat->SetMaterialPropertiesTable(mptWorld);

    //Charge Exchange material
    CMat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
}

G4VPhysicalVolume *MyDetectorConstruction::Construct() //defines volume and material
{
    //resets geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    G4GeometryManager::GetInstance()->CloseGeometry();

    //world volume calculations, all values are half of width
    G4double xWorld = 0.5*m;
    G4double yWorld = 0.5*m;
    G4double zWorld = 26*m;
    solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    physWorld = new G4PVPlacement(0, G4ThreeVector(0. ,0. ,0. ), logicWorld, "physWorld", 0, false, 0, false); 

    //material
    G4double totalX = 1.0 * m;
    G4double totalY = 1.0 * m;
    G4double cellX = totalX / nCols;
    G4double cellY = totalY / nRows;

    solidMat = new G4Box("solidMat", cellX / 2.0,cellY / 2.0,0.05 * m);
    logicMat = new G4LogicalVolume(solidMat, CMat, "logicMat");

    for (G4int i = 0; i < nBlocks; i++) {
        G4double z = 25.0 * m * (i + 1) / nBlocks;
        for (G4int row = 0; row < nRows; row++) {
            for (G4int col = 0; col < nCols; col++) {
                G4double x = -totalX / 2.0+ cellX / 2.0 + col * cellX;
                G4double y = -totalY / 2.0 + cellY / 2.0 + row * cellY;
                G4int copyNo =i*nRows*nCols + row*nCols + col;
                physMat = new G4PVPlacement(0, G4ThreeVector(x, y, z),logicMat,"physMat",
                logicWorld,false,copyNo, false);
            }
        }
    }

    G4double randX = 0;
    G4double randY = 0;
    G4double randZ = 0;

    if (randPosDetector == 1){
        randX = randPosScaling*G4RandGauss::shoot(0, 0.5*mm);
        randY = randPosScaling*G4RandGauss::shoot(0, 0.5*mm);
        randZ = randPosScaling*G4RandGauss::shoot(0, 0.5*cm);
    }

    G4ThreeVector detectorPos = G4ThreeVector(randX,randY, randZ+25.2*m);

    solidDetector = new G4Box("solidDetector", 0.49*m, 0.49*m, 0.1*m);
    logicDetector = new G4LogicalVolume(solidDetector, CMat, "logicDetector");
    physDetector = new G4PVPlacement(0, detectorPos, logicDetector,"physDetector",
                    logicWorld,false,0, false);

    return physWorld;
}

//this for making the detector (need here b/c we gonna use it to find boundaries)
void MyDetectorConstruction::ConstructSDandField(){
    MySensitiveDetector* sensDet = new MySensitiveDetector("SensitiveDetector");

    G4SDManager::GetSDMpointer()->AddNewDetector(sensDet);

    logicMat->SetSensitiveDetector(sensDet);
    logicDetector->SetSensitiveDetector(sensDet);
}