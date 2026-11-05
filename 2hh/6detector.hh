#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4TrackStatus.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//root stuff
#include "G4AnalysisManager.hh"
#include "G4VNtupleManager.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <cmath>
#include <set>
#include <vector>

//a track's ParticleInfo rows are written only if it hit at least this many distinct tracker planes, or it is a charge exchange proton
//(the CE protons are kept so CountCE and the CE-lost counts stay complete). Written to settings.csv by 7run.cc.
constexpr G4int kMinTrackerPlanesSaved = 2;

//one row of the ParticleInfo ntuple (booked in 7run.cc), held until the track ends
struct ParticleRow
{
    G4int boundaryType;
    G4ThreeVector pos;
    G4double KE;
    G4ThreeVector detectorCenter;
    G4double measuredKE;
};

//rows of the track being tracked on this worker thread, a thread tracks one track at a time
inline thread_local std::vector<ParticleRow> gTrackRows;

//per-event tracker counts of this worker thread, reset in BeginOfEventAction and written to EventSummary in EndOfEventAction (8event.cc)
//hits: every charged tracker entry of every particle, tracks: charged tracks that hit at least kMinTrackerPlanesSaved distinct planes
inline thread_local G4int gEventTrackerHits = 0;
inline thread_local G4int gEventTrackerTracks = 0;

//boundaryType: 0 tracker entry, 1 tracker exit, 2 ZDC, 3 track start
//the row is only buffered here, FlushParticleInfo (MyTrackingAction::PostUserTrackingAction) decides if the track is saved
//only charged particles are recorded: a neutral track gets no row and no tracker hit count, whatever the call site
inline void FillParticleInfo(const G4Track *track, G4int boundaryType, const G4ThreeVector &pos, G4double KE,
                             const G4ThreeVector &detectorCenter, G4double measuredKE)
{
    if (track->GetParticleDefinition()->GetPDGCharge() == 0.){return;}
    gTrackRows.push_back({boundaryType, pos, KE, detectorCenter, measuredKE});
    if (boundaryType == 0){gEventTrackerHits++;}
}

//writes the buffered rows of the finished track to ntuple 0 if it is worth saving, then empties the buffer
inline void FlushParticleInfo(const G4Track *track)
{
    //distinct planes hit, a plane is told apart by the z of the tracker center (1 um steps, safe against rounding)
    std::set<long> planes;
    for (const auto &row : gTrackRows){
        if (row.boundaryType == 0){planes.insert(std::lround(row.detectorCenter.z() / um));}
    }

    if (static_cast<G4int>(planes.size()) >= kMinTrackerPlanesSaved){gEventTrackerTracks++;}

    const G4VProcess *creator = track->GetCreatorProcess();
    const bool ceProton = track->GetParticleDefinition()->GetPDGEncoding() == 2212 && creator != nullptr &&
                          creator->GetProcessName() == "neutronChargeExNP";

    if (static_cast<G4int>(planes.size()) >= kMinTrackerPlanesSaved || ceProton){
        G4AnalysisManager *man = G4AnalysisManager::Instance();
        const G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
        for (const auto &row : gTrackRows){
            man->FillNtupleIColumn(0, 0, eventID);
            man->FillNtupleIColumn(0, 1, track->GetTrackID());
            man->FillNtupleIColumn(0, 2, track->GetParentID());
            man->FillNtupleIColumn(0, 3, track->GetParticleDefinition()->GetPDGEncoding());
            man->FillNtupleIColumn(0, 4, row.boundaryType);
            man->FillNtupleDColumn(0, 5, row.pos.x());
            man->FillNtupleDColumn(0, 6, row.pos.y());
            man->FillNtupleDColumn(0, 7, row.pos.z());
            man->FillNtupleDColumn(0, 8, row.KE);
            man->FillNtupleDColumn(0, 9, row.detectorCenter.x());
            man->FillNtupleDColumn(0, 10, row.detectorCenter.y());
            man->FillNtupleDColumn(0, 11, row.detectorCenter.z());
            man->FillNtupleDColumn(0, 12, row.measuredKE);
            man->AddNtupleRow(0);
        }
    }
    gTrackRows.clear();
}

class MySensitiveDetector : public G4VSensitiveDetector
{
public:
    MySensitiveDetector(G4String);
    ~MySensitiveDetector();
    void SetTrackVolume(const G4LogicalVolume *lv) {trackLV = lv;} //logicTrack, compared by pointer in ProcessHits
private:
    const G4LogicalVolume *trackLV = nullptr;
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *) override;
};

#endif
