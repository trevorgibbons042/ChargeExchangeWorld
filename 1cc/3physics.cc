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
#include "G4ChargeExchange.hh"
#include "G4ProcessTable.hh"
#include "G4HadronicProcess.hh"

MyPhysicsList::MyPhysicsList()
{
    fCEMessenger = new G4GenericMessenger(this, "/CE/", "Charge exchange settings");
    fCEMessenger->DeclareProperty("APScale", APScale, "Factor on the analyzing power, AP is clamped to -1.0 at the low end (1 = regular value)");

    // only use physics u need, otherwise too laggy
    RegisterPhysics (new G4EmStandardPhysics());
    RegisterPhysics (new G4EmExtraPhysics());
    RegisterPhysics (new G4DecayPhysics());
    RegisterPhysics (new G4IonPhysics());
    RegisterPhysics (new G4StoppingPhysics());
    RegisterPhysics (new G4NeutronTrackingCut());
    RegisterPhysics (new G4HadronPhysicsFTFP_BERT_ATL());
    RegisterPhysics (new G4HadronElasticPhysics());

    auto* ChargeExchangePhysics = new G4ChargeExchangePhysics();
    //ChargeExchangePhysics->SetCrossSectionFactor(3e20);

    RegisterPhysics (ChargeExchangePhysics);
};

MyPhysicsList::~MyPhysicsList()
{delete fCEMessenger;};

G4double MyPhysicsList::APScale = 1.0;

//the NP models are per thread and built at /run/initialize, before the macro sets /CE/APScale, so this runs every run
void MyPhysicsList::ApplyAPScale()
{
    auto* ceProcess = dynamic_cast<G4HadronicProcess*>(G4ProcessTable::GetProcessTable()->FindProcess("neutronChargeExNP", "neutron"));
    if (ceProcess != nullptr){
        for (auto* model : ceProcess->GetHadronicInteractionList()){
            if (auto* ceModel = dynamic_cast<G4ChargeExchange*>(model)){ceModel->SetAPScale(APScale);}
        }
    }
};

/*
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
*/