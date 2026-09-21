#include "8event.hh"

MyEventAction::MyEventAction(MyRunAction*)
{
    fEdep = 0;
};

MyEventAction::~MyEventAction(){
};

void MyEventAction::BeginOfEventAction(const G4Event*){
    fEdep = 0; //reset the value when new event starts
}

void MyEventAction::EndOfEventAction(const G4Event*){
    ////Energy Deposition (for 7run.cc, comes from adding up steps in 9stepping.cc)
    //G4cout << "Energy deposition: " << fEdep << G4endl;
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    //man->FillNtupleDColumn(2,0,fEdep);
    //man->AddNtupleRow(2);
}