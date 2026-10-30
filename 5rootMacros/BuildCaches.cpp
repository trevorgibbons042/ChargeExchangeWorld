#include <string>
#include <vector>

#include "CutCache.h"
#include "FitCore.h"
#include "Plot3DCols.h"

// Builds every cut cache the plot macros ask for (CutCache.h, spec A.7 C1 to C9) once, with all the threads, so no plot lane builds one or waits on another's lock.
// Run after CutTrees.cpp (needs cuts/cut_*.root from CutEvents.cpp, cuts/cut_angles_all.root from ReconAngleAll.cpp). A macro run alone still builds what it needs, under the lock.
void BuildCaches(){
    ROOT::EnableImplicitMT(NThreads());
    const std::string raw = "output0_t*.root", ce = "ChargeExchange_values", pi = "ParticleInfo";
    const std::string carbon = CarbonCeCut();
    GetCutTree(raw, ce, {carbon});                                              // C1: seed cut_CE
    GetCutTree(raw, ce, {"didChargeExchange>0"});                               // C2
    GetCutTree(raw, ce, {"didChargeExchange>0", "!(" + carbon + ")"});          // C3: from C2
    GetCutTree(raw, pi, {"@isParentCE", "boundaryType == 2 && pdg == 2212"});   // C4: from seed cut_ParentID
    GetCutTree(raw, pi, {"@isParentCE", "boundaryType == 0 && pdg == 2212"});   // C5
    GetCutTree(raw, pi, {"boundaryType == 2 && pdg == 2212"});                  // C6: seed cut_AllEvProtons
                                                                                // C7 to C9 load the maps themselves when they build (Plot3DColumns)
    Plot3DCE();                                                                 // C7: C4 + p3d_v1
    Plot3DAllEv();                                                              // C8: C6 + p3d_v1
    Plot3DRaw();                                                                // C9: AnglesRaw + p3d_v1
}
