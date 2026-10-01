#include "3physics.hh"
#include "G4ChargeExchangePhysics.hh"

#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4BosonConstructor.hh"

#include "G4HadronElasticPhysics.hh"
#include "G4HadronicParameters.hh"
#include "G4ChipsElasticModel.hh"

#include "G4ChargeExchangeNP.hh"

MyPhysicsList::MyPhysicsList()
{
    // only use physics u need, otherwise too laggy
    RegisterPhysics (new G4EmStandardPhysics());
    //RegisterPhysics (new G4OpticalPhysics());
    RegisterPhysics (new G4EmExtraPhysics());
    RegisterPhysics (new G4DecayPhysics());
    RegisterPhysics (new G4IonPhysics());
    RegisterPhysics (new G4StoppingPhysics());
    RegisterPhysics (new G4NeutronTrackingCut());
    RegisterPhysics (new G4HadronPhysicsFTFP_BERT_ATL());
    RegisterPhysics (new G4HadronElasticPhysics());


    auto* ChargeExchangePhysics = new G4ChargeExchangePhysics();
    //ChargeExchangePhysics->SetCrossSectionFactor(3e20);

    RegisterPhysics (new G4ChargeExchangePhysics);
};

MyPhysicsList::~MyPhysicsList()
{};

void MyPhysicsList::ConstructParticle(){
    G4BosonConstructor bosons;
    bosons.ConstructParticle();
    
    G4LeptonConstructor leptons;
    leptons.ConstructParticle();

    G4MesonConstructor mesons;
    mesons.ConstructParticle();

    G4BaryonConstructor baryons;
    baryons.ConstructParticle();

    G4ShortLivedConstructor shortLived;
    shortLived.ConstructParticle();

    G4IonConstructor ions;
    ions.ConstructParticle();
}