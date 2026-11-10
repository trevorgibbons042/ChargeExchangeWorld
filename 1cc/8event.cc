#include "8event.hh"
#include "G4Event.hh"
#include "6detector.hh" //gEventTrackerHits, gEventTrackerTracks

MyEventAction::MyEventAction(MyRunAction* runAction) : fRunAction(runAction)
{};

MyEventAction::~MyEventAction(){};

void MyEventAction::BeginOfEventAction(const G4Event*)
{
    fNZdcAll = 0; fNZdcCharged = 0; fZdcKEAll = 0.; fZdcKECharged = 0.;
    gEventTrackerHits = 0; gEventTrackerTracks = 0;
    gEventHitZ = &fRunAction->hitZ; gEventHitKE = &fRunAction->hitKE;
    gEventHitZ->clear(); gEventHitKE->clear();
}

//one EventSummary row per event (ntuple 2, booked in 7run.cc): the particles that reached the ZDC, with no tracker or reconstruction requirement
void MyEventAction::EndOfEventAction(const G4Event* event)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->FillNtupleIColumn(2, 0, event->GetEventID());
    man->FillNtupleIColumn(2, 1, fNZdcAll);
    man->FillNtupleIColumn(2, 2, fNZdcCharged);
    man->FillNtupleDColumn(2, 3, fZdcKEAll);
    man->FillNtupleDColumn(2, 4, fZdcKECharged);
    man->FillNtupleIColumn(2, 5, gEventTrackerHits);
    man->FillNtupleIColumn(2, 6, gEventTrackerTracks);
    man->AddNtupleRow(2);
}