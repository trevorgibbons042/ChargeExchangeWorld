#include "5generator.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
    fMessenger = new G4GenericMessenger(this, "/Polarization/", "Gun Properties");
    fMessenger->DeclareProperty("PolarizationIndex", PolarizationIndex, "Polarization Vector Index Number");
    fMessenger->DeclareProperty("PPIndex", PPIndex, "Polarization Percentage Index Number");
    fMessenger->DeclareMethod("PrintPolarization", &MyPrimaryGenerator::PrintPolarization, "Printing Polarization from Polarization Index");
    fMessenger->DeclareMethod("PrintPPIndex", &MyPrimaryGenerator::PrintPPIndex, "Print PPIndex");

    fMessenger = new G4GenericMessenger(this, "/randGun/", "random Gun Properties");
    fMessenger->DeclareProperty("randPAngle", randPAngle, "Is random angle on?");
    fMessenger->DeclareProperty("randPos", randPos, "Is random position on?");
    fMessenger->DeclareProperty("EnergyBool", EnergyBool, "Is random energy on?");
    fMessenger->DeclareProperty("EnergySD", EnergySD, "Standard Deviation on Energy");

    //number of particles per event (can do 1 run with bunch of events tho)
    fParticleGun = new G4ParticleGun(1);

    //type of particle
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName = "neutron";
    G4ParticleDefinition *particle = particleTable->FindParticle("neutron");
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
    G4ThreeVector polarization(0,0,0);

    G4double theta = 0;
    G4double phi = 0;
    G4double randX = 0;
    G4double randY = 0;
    G4double randZ = 0;

    if (randPAngle == 1){
    theta = .005*G4UniformRand();
    phi = twopi*G4UniformRand();
    }


    G4ThreeVector mom(std::sin(theta)*std::cos(phi),
                std::sin(theta)*std::sin(phi),
                cos(theta));

    if (randPos == 1){
    randX = G4RandGauss::shoot(0, 0.011/2)*mm;
    randY = G4RandGauss::shoot(0, 0.006/2)*mm;
    randZ = G4RandGauss::shoot(0, 3/2)*cm;
    }

    G4ThreeVector pos(randX, randY, randZ);

    G4double Energy = 25*GeV;

    if (EnergyBool == true){
        Energy = 25*GeV*G4RandGauss::shoot(0,EnergySD);
    }
    fParticleGun->SetParticleMomentum(25*GeV);

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);

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