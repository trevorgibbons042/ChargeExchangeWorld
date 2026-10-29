#include "5generator.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include <algorithm>
#include <cmath>
#include <mutex>

G4double MyPrimaryGenerator::PPValue = 0;
G4double MyPrimaryGenerator::PPSigma = 0;

namespace {
    //sums over every primary of every worker thread: n, momentum, momentum^2, KE, KE^2 (MeV/c, MeV)
    std::mutex beamMutex;
    G4double beamN = 0, beamP = 0, beamP2 = 0, beamKE = 0, beamKE2 = 0;
}

void MyPrimaryGenerator::GetBeamStats(G4double& meanMomentum, G4double& sdMomentum, G4double& meanKE, G4double& sdKE)
{
    std::lock_guard<std::mutex> lock(beamMutex);
    meanMomentum = beamN > 0 ? beamP/beamN : 0;
    meanKE = beamN > 0 ? beamKE/beamN : 0;
    sdMomentum = beamN > 0 ? std::sqrt(std::max(0.0, beamP2/beamN - meanMomentum*meanMomentum)) : 0;
    sdKE = beamN > 0 ? std::sqrt(std::max(0.0, beamKE2/beamN - meanKE*meanKE)) : 0;
}

MyPrimaryGenerator::MyPrimaryGenerator()
{
    fMessenger = new G4GenericMessenger(this, "/Polarization/", "Gun Properties");
    fMessenger->DeclareProperty("PolarizationIndex", PolarizationIndex, "Polarization Vector Index Number");
    fMessenger->DeclareMethod("PPValue", &MyPrimaryGenerator::SetPPValue, "Polarization fraction, -1.0 to 1.0");
    fMessenger->DeclareMethod("PPSigma", &MyPrimaryGenerator::SetPPSigma, "Gaussian sigma on polarization fraction, 0 or more; each event is clamped to -1.0 to 1.0");
    fMessenger->DeclareMethod("PrintPolarization", &MyPrimaryGenerator::PrintPolarization, "Printing Polarization from Polarization Index");
    fMessenger->DeclareMethod("PrintPPValue", &MyPrimaryGenerator::PrintPPValue, "Print PPValue and PPSigma");

    fRandMessenger = new G4GenericMessenger(this, "/randGun/", "random Gun Properties");
    fRandMessenger->DeclareProperty("randPAngle", randPAngle, "Is random angle on?");
    fRandMessenger->DeclareProperty("randPos", randPos, "Is random position on?");
    fRandMessenger->DeclareProperty("EnergySD", EnergySD, "Standard Deviation on Energy");

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
    delete fRandMessenger;
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

    G4double Momentum = 25*GeV;

    if (EnergySD > 0){
        Momentum = 25*GeV*G4RandGauss::shoot(1,EnergySD);
    }
    //SetParticleMomentum prints a message on every call after the first, so convert the momentum to kinetic energy here.
    G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
    const G4double kineticEnergy = std::sqrt(Momentum*Momentum + mass*mass) - mass;
    fParticleGun->SetParticleEnergy(kineticEnergy);
    {
        std::lock_guard<std::mutex> lock(beamMutex);
        const G4double p = Momentum/MeV, ke = kineticEnergy/MeV;
        beamN++; beamP += p; beamP2 += p*p; beamKE += ke; beamKE2 += ke*ke;
    }

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);

    if (PolarizationIndex > -1){
        polarization = polarizationTable.at(PolarizationIndex);
        fParticleGun->SetParticlePolarization(polarization);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    } else if (PPOn){

        G4double PPNumber = PPValue;
        if (PPSigma > 0){
            PPNumber = std::clamp(G4RandGauss::shoot(PPValue, PPSigma), -1.0, 1.0);
        }

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

void MyPrimaryGenerator::SetPPValue(G4double value){
    if (value < -1.0 || value > 1.0){
        G4cerr << "PPValue " << value << " outside [-1, 1], ignored" << G4endl;
        return;
    }
    PPValue = value;
    PPOn = true;
}

void MyPrimaryGenerator::SetPPSigma(G4double sigma){
    if (sigma < 0){
        G4cerr << "PPSigma " << sigma << " is negative, ignored" << G4endl;
        return;
    }
    PPSigma = sigma;
}

void MyPrimaryGenerator::PrintPPValue(){
    G4cout << "PPValue = " << PPValue << ", PPSigma = " << PPSigma
           << ", on = " << PPOn << G4endl;
}