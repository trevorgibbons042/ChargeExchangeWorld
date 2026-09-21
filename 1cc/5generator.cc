#include "5generator.hh"

#include "Randomize.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
    fMessenger = new G4GenericMessenger(this, "/Polarization/", "Gun Properties");
    fMessenger->DeclareProperty("PolarizationIndex", PolarizationIndex, "Polarization Vector Index Number");
    fMessenger->DeclareProperty("PPIndex", PPIndex, "Polarization Percentage Index Number");
    fMessenger->DeclareMethod("PrintPolarization", &MyPrimaryGenerator::PrintPolarization, "Printing Polarization from Polarization Index");
    fMessenger->DeclareMethod("PrintPPIndex", &MyPrimaryGenerator::PrintPPIndex, "Print PPIndex");
    
    //number of particles per event (can do 1 run with bunch of events tho)
    fParticleGun = new G4ParticleGun(1);

    //type of particle
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName = "neutron";
    G4ParticleDefinition *particle = particleTable->FindParticle("neutron");

    //inital position and direction
    G4ThreeVector pos(0. ,0. , 0. );
    G4ThreeVector mom(0., 0. , 1. );

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);

    //sets momentum
    fParticleGun->SetParticleMomentum(25*GeV);
    fParticleGun->SetParticleDefinition(particle);
};

MyPrimaryGenerator::~MyPrimaryGenerator()
{
    delete fParticleGun;
    delete fMessenger;
};

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
    //put here if u want all runs to be the same, only ran with run starts
    G4ThreeVector polarization;
    G4double PPNumber;

    if (PolarizationIndex > -1){
        polarization = polarizationTable.at(PolarizationIndex);
        fParticleGun->SetParticlePolarization(polarization);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    } else if (PPIndex > -1){
        
        G4double PPNumber = PPTable[PPIndex];

        if (G4UniformRand() < std::abs(PPNumber)){
            if (PPNumber > 0){
                polarization = G4ThreeVector(0,1,0);
                fParticleGun->SetParticlePolarization(polarization);
                fParticleGun->GeneratePrimaryVertex(anEvent);
            }
            else{
                polarization = G4ThreeVector(0,-1,0); 
                fParticleGun->SetParticlePolarization(polarization);
                fParticleGun->GeneratePrimaryVertex(anEvent);}
        }
        else{
            G4double RandomTheta = (G4UniformRand()*CLHEP::twopi - CLHEP::twopi/2);
            G4double RandomPhi = G4UniformRand()*CLHEP::twopi;

            polarization = G4ThreeVector(
            std::sin(RandomTheta)*std::cos(RandomPhi), 
            std::sin(RandomTheta)*std::sin(RandomPhi), 
            std::cos(RandomTheta));

            fParticleGun->SetParticlePolarization(polarization);
            fParticleGun->GeneratePrimaryVertex(anEvent);
        }
    } 
    else {polarization = G4ThreeVector(0,0,0);
        fParticleGun->SetParticlePolarization(polarization);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
};

void MyPrimaryGenerator::PrintPolarization(){
    const G4ThreeVector& polarization = polarizationTable.at(PolarizationIndex);
    G4cout << "PolarizationIndex = " << PolarizationIndex
           << ", polarization = " << polarization
           << G4endl;
}

void MyPrimaryGenerator::PrintPPIndex(){
    G4cout << "PPIndex = " << PPIndex << G4endl;
    G4cout << "Percentage = " << PPTable[PPIndex] << G4endl;
}