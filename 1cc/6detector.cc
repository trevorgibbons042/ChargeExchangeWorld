#include "6detector.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
    ////Detector Efficency Input
    ////1st colm: wavelength, 2nd colm: efficiency percentage (based on tests)
    //quEff = new G4PhysicsOrderedFreeVector();
    //std::ifstream datafile;
    //datafile.open("4dat/1eff.dat");
    //while(1){
    //    G4double wlenPhoton, quef;
    //    datafile >> wlenPhoton >> quef;
    //    if(datafile.eof()){break;};
    //    G4cout << wlenPhoton << " " << quef << G4endl;
    //    quEff->InsertValues(wlenPhoton, quef/100);
    //};
    //datafile.close();

};

MySensitiveDetector::~MySensitiveDetector()
{};

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
    ////Gives info on track when entering volume and kills photons    
    G4Track *track = aStep->GetTrack();

    //track->SetTrackStatus(fStopAndKill);

    ////Defining when photons enter and exit
    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

    ////1. finds wavelength of the photon 
    //G4ThreeVector momPhoton = preStepPoint->GetMomentum();
    //G4double magmomPhoton = momPhoton.mag();
    //G4double wlenPhoton = (1.239841939*eV)/(magmomPhoton)*(1E+03);

    ////2. Photon true position

    ////3. Copy number (detector integer ID)
    //const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable(); //gets which are touching
    //G4int copyNo = touchable->GetCopyNumber(); //puts that number into copyNo
    //G4cout << "Copy number: " << copyNo << G4endl;

    ////4. Position of Detector hit (from copy number)
    //G4VPhysicalVolume *physVol = touchable->GetVolume(); //gets which are touching
    //G4ThreeVector posDetector = physVol->GetTranslation(); //get which number detector it touched
    //G4cout << "Detector position:" << posDetector << G4endl;

    ////Root table pt2, NOTE - You need to create empty tables in 7run.cc
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    
    ////Position of Photons (in detector, with efficiency, check 7run.cc)
    //if(G4UniformRand() < quEff->Value(wlenPhoton)){
    //man->FillNtupleIColumn(0, 0, evt); //FillNtupleIColumn()
    //man->FillNtupleDColumn(0, 1, posDetector[0]); //for fX
    //man->FillNtupleDColumn(0, 2, posDetector[1]); //for fY
    //man->FillNtupleDColumn(0, 3, posDetector[2]); //for fZ
    //man->AddNtupleRow(0);
    //};

    /*
    ////Position of Photon (true position, check 7run.cc)
    man->FillNtupleIColumn(0, 0, evt); //for "Hits"
    man->FillNtupleDColumn(0, 1, posTrue[0]); //for fX
    man->FillNtupleDColumn(0, 2, posTrue[1]); //for fY
    man->FillNtupleDColumn(0, 3, posTrue[2]); //for fZ
    //man->FillNtupleDColumn(0, 4, wlenPhoton); //for wavelength
    man->AddNtupleRow(0);
    */

return true;
}