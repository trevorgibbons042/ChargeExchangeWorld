#include "7run.hh"
#include <cassert>
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
int main() {
    G4RunManager runManager;
    // This isolated analysis lifecycle check does not transport particles.
    G4ParticleTable::GetParticleTable()->SetReadiness(true);
    MyRunAction action;
    auto* manager = G4AnalysisManager::Instance();
    for (int id = 0; id < 2; ++id) {
        G4Run run;
        run.SetRunID(id);
        action.BeginOfRunAction(&run);
        assert(manager->IsOpenFile());
        manager->FillNtupleIColumn(0, id);
        manager->FillNtupleDColumn(1, 1.0);
        manager->FillNtupleDColumn(2, 2.0);
        manager->FillNtupleDColumn(3, 3.0);
        manager->AddNtupleRow();
        action.EndOfRunAction(&run);
        assert(!manager->IsOpenFile());
    }
}
