#include "2construction.hh"
#include "Randomize.hh"
#include "G4SDManager.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
    ////Messengar. allows us to make customized gui commands (like change number of rows on detector)    
    fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");
    fMessenger->DeclareProperty("nBlocks", nBlocks, "Number of Blocks");
    fMessenger->DeclareProperty("nBlocksTracking", nBlocksTracking, "Number of Tracking Blocks");
    fMessenger->DeclareProperty("patchCells", patchCells, "Display only, build just the N x N detector cells at the lower left corner (0 = all cells)");

    fMessenger->DeclareProperty("nCols_blocks", nCols, "Number of columns for tracking");
    fMessenger->DeclareProperty("nRows_blocks", nRows, "Number of rows for tracking");
    fMessenger->DeclareProperty("nCols_detector", nCols_detector, "Number of columns for detector");
    fMessenger->DeclareProperty("nRows_detector", nRows_detector, "Number of rows for detector");

    fMessenger->DeclareProperty("dz", dz, "Distance between trackers and blocks (in meters)");
    fMessenger->DeclareProperty("CarbonThickness", CarbonThickness, "Thickness of carbons (Meters)");
    
    fMessenger->DeclareProperty("randPosDetector", randPosDetector, "Is detector in random position?");
    fMessenger->DeclareProperty("randPosScaling", randPosScaling, "randPos Scaling Factor");

    //default values
    nBlocks = 3; //default 3 blocks with their trackers, /detector/nBlocks 0 leaves only the ZDC
    nBlocksTracking = 2;
    nCols = 11;
    nRows = 11;
    nCols_detector= 5;
    nRows_detector = 5;
    dz = 0.025;
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
    CMat = nist->FindOrBuildMaterial("G4_POLYETHYLENE"); //CH2 blocks, 0.94 g/cm3
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
    //ZDC: 0.6 m deep, centred at z = 8 m, its front face at z = 7.7 m is the last of the evenly spaced slots (stackSpread = 1 in settings.csv)
    G4ThreeVector detectorPos = G4ThreeVector(0,0, 8*m);
    G4double detectorX = 0.4*m;
    G4double detectorY = 0.4*m;
    G4double detectorZ = 0.6*m;

    solidMat = new G4Box("solidMat", 0.5 * m,0.5 * m, (CarbonThickness * m)/2);
    logicMat = new G4LogicalVolume(solidMat, CMat, "logicMat");

    G4double totalX = 0.4 * m;
    G4double totalY = 0.4 * m;
    G4double cellX = totalX / nCols;
    G4double cellY = totalY / nRows;
    G4double detectorThickness = 500 * um;

    //each grid of cells is one parameterised volume inside a slab, so a plane costs one volume, not nCols*nRows
    solidTrack = new G4Box("solidTrack",cellX / 2.0,cellY / 2.0,detectorThickness / 2.0);
    logicTrack = new G4LogicalVolume(solidTrack, DMat, "logicTrack");

    G4Box *solidTrackSlab = new G4Box("solidTrackSlab", totalX / 2.0, totalY / 2.0, detectorThickness / 2.0);
    G4LogicalVolume *logicTrackSlab = new G4LogicalVolume(solidTrackSlab, DMat, "logicTrackSlab");
    physTrack = new G4PVParameterised("physTrack", logicTrack, logicTrackSlab, kUndefined, nRows * nCols,
                    new GridParameterisation(nCols, cellX, cellY, totalX, totalY));

    //the blocks sit in evenly spaced slots between z = 0 and the ZDC front face, the ZDC face being the last slot. Every block has a tracker dz behind it
    //and, except the first, one dz in front of it; one more tracker sits dz in front of the ZDC. dz is the gap between the edge of the block and the edge of the
    //tracker (the ZDC face and the tracker for the last one). nBlocksTracking is not used by this layout.
    G4double zdcFront = detectorPos.z() - detectorZ / 2.0;
    G4double pitch = zdcFront / (nBlocks + 1);
    G4double trackerOffset = (CarbonThickness * m) / 2.0 + dz * m + detectorThickness / 2.0; //block centre to tracker centre
    G4int trackerIndex = 0;
    auto placeTracker = [&](G4double zTracker) {
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zTracker), logicTrackSlab, "physTrackSlab",
                    logicWorld, false, trackerIndex, false);
        trackerIndex++;
    };

    for (G4int i = 0; i < nBlocks; i++) {
        G4double zCarbon = (i + 1) * pitch;
        physMat = new G4PVPlacement(nullptr,
            G4ThreeVector(0, 0, zCarbon), logicMat, "physMat",
            logicWorld,false,i,false);

        if (i > 0) placeTracker(zCarbon - trackerOffset);
        placeTracker(zCarbon + trackerOffset);
    }
    if (nBlocks > 0) placeTracker(zdcFront - dz * m - detectorThickness / 2.0); //no blocks, no trackers

    G4double cellXd = detectorX / nCols_detector;
    G4double cellYd = detectorY / nRows_detector;

    solidDetector = new G4Box("solidDetector", cellXd / 2.0, cellYd / 2.0, detectorZ / 2.0);
    logicDetector = new G4LogicalVolume(solidDetector, DMat, "logicDetector");

    //patchCells > 0 builds only the first patchCells x patchCells cells at the lower left corner of the grid (same pitch), so the GUI can draw it
    G4bool patch = patchCells > 0 && patchCells < nCols_detector && patchCells < nRows_detector;
    G4int detCols = patch ? patchCells : nCols_detector;
    G4int detRows = patch ? patchCells : nRows_detector;
    G4double slabX = detCols * cellXd;
    G4double slabY = detRows * cellYd;
    G4ThreeVector slabPos = detectorPos;
    if (patch) slabPos += G4ThreeVector(-detectorX / 2.0 + slabX / 2.0, -detectorY / 2.0 + slabY / 2.0, 0);

    G4Box *solidDetectorSlab = new G4Box("solidDetectorSlab", slabX / 2.0, slabY / 2.0, detectorZ / 2.0);
    G4LogicalVolume *logicDetectorSlab = new G4LogicalVolume(solidDetectorSlab, DMat, "logicDetectorSlab");
    physDetector = new G4PVParameterised("physDetector", logicDetector, logicDetectorSlab, kUndefined,
                    detRows * detCols,
                    new GridParameterisation(detCols, cellXd, cellYd, slabX, slabY));
    new G4PVPlacement(nullptr, slabPos, logicDetectorSlab, "physDetectorSlab", logicWorld, false, 0, false);

    return physWorld;
}

//this for making the detector (need here b/c we gonna use it to find boundaries)
void MyDetectorConstruction::ConstructSDandField(){
    //reuse the detector after /run/reinitializeGeometry, registering it again triggers DET1010
    auto* sensDet = static_cast<MySensitiveDetector*>(G4SDManager::GetSDMpointer()->FindSensitiveDetector("SensitiveDetector", false));
    if (!sensDet) {
        sensDet = new MySensitiveDetector("SensitiveDetector");
        G4SDManager::GetSDMpointer()->AddNewDetector(sensDet);
    }
    sensDet->SetTrackVolume(logicTrack); //refreshed each call, geometry rebuilds change the pointer

    //logicMat (CH2) is not sensitive: its hits are never recorded and cost a call per step
    logicDetector->SetSensitiveDetector(sensDet);
    logicTrack->SetSensitiveDetector(sensDet);
}