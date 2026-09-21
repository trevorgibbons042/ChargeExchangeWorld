#include "2construction.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
    ////Messengar. allows us to make customized gui commands (like change number of rows on detector)    
    //fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");
    //fMessenger->DeclareProperty("nCols", nCols, "Number of columns");
    //fMessenger->DeclareProperty("nRows", nRows, "Number of rows");
    //nCols = 10; //default values
    //nRows = 10;

    DefineMaterial();
}

MyDetectorConstruction::~MyDetectorConstruction()
{}

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
    CMat = nist->FindOrBuildMaterial("G4_C");
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
    G4double zWorld = 0.5*m;
    solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    physWorld = new G4PVPlacement(0, G4ThreeVector(0. ,0. ,0. ), logicWorld, "physWorld", 0, false, 0, true); 

    //material
    G4Box *solidMat = new G4Box("solidMat", 0.3*m,0.2*m,0.1*m);
    G4LogicalVolume *logicMat = new G4LogicalVolume(solidMat, CMat, "logicMat");
    G4VPhysicalVolume *physMat = new G4PVPlacement(0, G4ThreeVector(0*m,0*m,0.25*m), logicMat, "physMat", logicWorld, false, 0, true);

    ////for 9stepping.cc
    //fScoringVolume = logicRadiator;

    return physWorld;
}

//this for making the detector (need here b/c we gonna use it to find boundaries)
void MyDetectorConstruction::ConstructSDandField(){
    ////For 6detector.cc
    //MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");
    //logicDetector->SetSensitiveDetector(sensDet);
}
