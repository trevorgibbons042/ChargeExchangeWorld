#include "2construction.hh"
#include "Randomize.hh"
#include "G4SDManager.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
    ////Messengar. allows us to make customized gui commands (like change number of rows on detector)    
    fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");
    fMessenger->DeclareProperty("nBlocks", nBlocks, "Number of Blocks");
    fMessenger->DeclareProperty("nBlocksTracking", nBlocksTracking, "Number of Tracking Blocks");
    fMessenger->DeclareProperty("nCols_blocks", nCols, "Number of columns for tracking");
    fMessenger->DeclareProperty("nRows_blocks", nRows, "Number of rows for tracking");
    fMessenger->DeclareProperty("nCols_detector", nCols_detector, "Number of columns for detector");
    fMessenger->DeclareProperty("nRows_detector", nRows_detector, "Number of rows for detector");
    fMessenger->DeclareProperty("dz", dz, "Distance between trackers and blocks (in meters)");
    fMessenger->DeclareProperty("CarbonThickness", CarbonThickness, "Thickness of carbons (Meters)");
    fMessenger->DeclareProperty("randPosDetector", randPosDetector, "Is detector in random position?");
    fMessenger->DeclareProperty("randPosScaling", randPosScaling, "randPos Scaling Factor");

    //default values
    nBlocks = 4;
    nBlocksTracking = 2;
    nCols = 11;
    nRows = 11;
    nCols_detector= 11;
    nRows_detector = 11;
    dz = 0.25;
    CarbonThickness = 0.2;


    DefineMaterial();
}

MyDetectorConstruction::~MyDetectorConstruction()
{
    delete fMessenger;
}

void MyDetectorConstruction::DefineMaterial(){
    //sets up material manager
    nist = G4NistManager::Instance();

    worldMat = nist->FindOrBuildMaterial("G4_Galactic");
    CMat = nist->FindOrBuildMaterial("G4_C");
    DMat = nist->FindOrBuildMaterial("G4_Si");
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
    G4double zWorld = 9*m;
    solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    physWorld = new G4PVPlacement(0, G4ThreeVector(0. ,0. ,0. ), logicWorld, "physWorld", 0, false, 0, false); 

    //material
    G4double zCenter = 4 * m;

    solidMat = new G4Box("solidMat", 0.5 * m,0.5 * m, (CarbonThickness * m)/2);
    logicMat = new G4LogicalVolume(solidMat, CMat, "logicMat");

    G4double totalX = 0.4 * m;
    G4double totalY = 0.4 * m;
    G4double cellX = totalX / nCols;
    G4double cellY = totalY / nRows;
    G4double detectorThickness = 0.002 * m;

    solidTrack = new G4Box("solidTrack",cellX / 2.0,cellY / 2.0,detectorThickness / 2.0);
    logicTrack = new G4LogicalVolume(solidTrack, DMat, "logicTrack");

    G4int nTrackerPlanes = (nBlocks + 1) * nBlocksTracking;
    G4int nElements = nTrackerPlanes + nBlocks;

    G4double zStart = zCenter - ((nElements - 1) * (dz * m)) / 2.0;

    G4int elementIndex = 0;
    G4int trackerIndex = 0;

    for (G4int i = 0; i < nBlocks; i++) {
        for (G4int j = 0; j < nBlocksTracking; j++) {
            G4double zTracker = zStart + elementIndex * (dz * m);
            for (G4int row = 0; row < nRows; row++) {
                for (G4int col = 0; col < nCols; col++) {

                    G4double x = -totalX / 2.0 + cellX / 2.0 + col * cellX;
                    G4double y =-totalY / 2.0 + cellY / 2.0 + row * cellY;
                    G4int copyNo = trackerIndex * nRows * nCols + row * nCols + col;

                    physTrack = new G4PVPlacement( nullptr,
                        G4ThreeVector(x, y, zTracker), logicTrack, "physTrack",
                        logicWorld,false,copyNo,false);
                }
            }
            trackerIndex++;
            elementIndex++;
        }
        G4double zCarbon = zStart + elementIndex * (dz * m);
        physMat = new G4PVPlacement(nullptr,
            G4ThreeVector(0, 0, zCarbon), logicMat, "physMat",
            logicWorld,false,i,false);

        elementIndex++;
    }


    // ============================================================
    // Final tracker group after the last carbon
    // ============================================================

    for (G4int j = 0; j < nBlocksTracking; j++) {

        G4double zTracker =
            zStart + elementIndex * (dz * m);

        for (G4int row = 0; row < nRows; row++) {
            for (G4int col = 0; col < nCols; col++) {

                G4double x =
                    -totalX / 2.0
                    + cellX / 2.0
                    + col * cellX;

                G4double y =
                    -totalY / 2.0
                    + cellY / 2.0
                    + row * cellY;

                G4int copyNo =
                    trackerIndex * nRows * nCols
                    + row * nCols
                    + col;

                physTrack = new G4PVPlacement(
                    nullptr,
                    G4ThreeVector(x, y, zTracker),
                    logicTrack,
                    "physTrack",
                    logicWorld,
                    false,
                    copyNo,
                    false
                );
            }
        }

        trackerIndex++;
        elementIndex++;
    }
    

    G4ThreeVector detectorPos = G4ThreeVector(0,0, 8*m);

    G4double detectorX = 0.4*m;
    G4double detectorY = 0.4*m;
    G4double detectorZ = 0.6*m;
    G4double cellXd = detectorX / nCols_detector;
    G4double cellYd = detectorY / nRows_detector;

    solidDetector = new G4Box("solidDetector", cellXd / 2.0, cellYd / 2.0, detectorZ / 2.0);
    logicDetector = new G4LogicalVolume(solidDetector, DMat, "logicDetector");
    
    for (G4int row = 0; row < nRows_detector; row++){
        for (G4int col = 0; col < nCols_detector; col++){
            G4double x = detectorPos.x()- detectorX / 2.0+ cellXd / 2.0+ col * cellXd;
            G4double y = detectorPos.y()- detectorY / 2.0+ cellYd / 2.0+ row * cellYd;
            G4double z = detectorPos.z();
            G4ThreeVector cellPos = G4ThreeVector(x, y, z);
            G4int copyNumber = row * nCols_detector + col;

            physDetector = new G4PVPlacement(0, cellPos, logicDetector,"physDetector",
                                            logicWorld,false,copyNumber, false);
        }
    }

    return physWorld;
}

//this for making the detector (need here b/c we gonna use it to find boundaries)
void MyDetectorConstruction::ConstructSDandField(){
    MySensitiveDetector* sensDet = new MySensitiveDetector("SensitiveDetector");

    G4SDManager::GetSDMpointer()->AddNewDetector(sensDet);

    logicMat->SetSensitiveDetector(sensDet);
    logicDetector->SetSensitiveDetector(sensDet);
    logicTrack->SetSensitiveDetector(sensDet);
}