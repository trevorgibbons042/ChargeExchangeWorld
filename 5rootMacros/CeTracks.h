#pragma once
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>

#include <map>
#include <set>
#include <vector>

#include "FitCore.h"

// CE neutrons (the selection of CutPlots.cpp (events)): event -> the track IDs of the CE neutrons. The carbon-block CE (CarbonCeCut, FitCore.h) are read straight from the sim output
// (ChargeExchange_values of output0_t*.root), so the recon does not wait for cuts/cut_CE.root. A tracker hit is a CE child proton when its ParentID is in the set of its event.
// blockZ: the carbon block centres (mm), halfBlockMm: half the block thickness, zdcFrontMm: the ZDC front face (StackZ(s).zdcFront). Needs ROOT::EnableImplicitMT first. False (and an Error) with no ChargeExchange_values.
inline bool CarbonCeTracks(const std::vector<double>& blockZ, double halfBlockMm, double zdcFrontMm, std::map<int, std::set<int>>& ceTracks, const char* caller){
    TChain ceChain("ChargeExchange_values");
    ceChain.Add("output0_t*.root");
    if (ceChain.LoadTree(0) < 0){Error(caller, "no ChargeExchange_values in output0_t*.root"); return false;}
    ROOT::RDataFrame dfCEAll(ceChain);
    // CarbonCeCut() as a typed lambda (a Filter string is JIT compiled, which fails on the Mac): a CE with fZ in front of the ZDC (fZ < zdcFrontMm) and inside one carbon block
    auto dfCE = dfCEAll.Filter([blockZ, halfBlockMm, zdcFrontMm](double did, double z){
        if (!(did > 0) || !(z < zdcFrontMm)) return false;
        for (const double zc : blockZ) if (z > zc - halfBlockMm && z < zc + halfBlockMm) return true;
        return false;
    }, {"didChargeExchange", "fZ"});
    std::vector<std::map<int, std::set<int>>> slotCeTracks(dfCEAll.GetNSlots());
    dfCE.ForeachSlot([&](unsigned slot, int eventID, int trackID){
        slotCeTracks[slot][eventID].insert(trackID);
    }, {"fEvent", "TrackID"});
    for (auto &m : slotCeTracks) {
        for (auto &[ev, ids] : m) ceTracks[ev].insert(ids.begin(), ids.end());
    }
    return true;
}
