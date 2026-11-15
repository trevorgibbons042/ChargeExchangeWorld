#pragma once
// The inputs of the signal vs background classifier of TmvaTrain.cpp and TmvaApply.cpp, in one place so training and apply cannot disagree.
// Every one is a column of AnglesRaw (ReconAngleAll.cpp) that a real detector also measures: the pixelated (D) positions and fit, the smeared ZDC KE and hit counts. No truth column and no label (inCE, isFromCE, evClass, nCEUsed...).
// Signal = a row whose fit used a CE child hit (isFromCE == 1), background = isFromCE == 0.
#include <cmath>
#include <string>
#include <vector>
#include "FitCore.h"

struct TmvaVar {std::string name; char type; std::string unit;}; // type 'D' double or 'I' int, as in the AnglesRaw branch

// The variables of a run: nChargedPlane<i>, the charged hits of the event on tracker plane i, has one column per plane of the run (StackZ, FitCore.h), so a run's block count sets the list.
// A classifier is trained and applied on runs with the same block count.
inline std::vector<TmvaVar> TmvaVariables(){
    std::vector<TmvaVar> v = {
        {"zdcKEAllD", 'D', "MeV"}, {"outDTheta", 'D', "rad"}, {"outDPhi", 'D', "rad"}, {"resDRms", 'D', "mm"}, {"impactD", 'D', "mm"},
        {"firstXD", 'D', "mm"}, {"firstYD", 'D', "mm"}, {"lastXD", 'D', "mm"}, {"lastYD", 'D', "mm"},
        {"nProtonsEv", 'I', ""}};
    const std::size_t nPlanes = StackZ(ReadSettings()).planeZ.size();
    // With 2 planes every fit uses both: nPlanesUsed and nHits are 2 on every row, r2D needs 3 hits and resHalfRatioD more than 2, so they are missing on every row. TMVA stops on a variable with no spread, so they join from 3 planes on.
    if (nPlanes >= 3){for (const char* n : {"r2D", "resHalfRatioD"}) v.push_back({n, 'D', ""}); for (const char* n : {"nPlanesUsed", "nHits"}) v.push_back({n, 'I', ""});}
    for (std::size_t i = 0; i < nPlanes; i++) v.push_back({"nChargedPlane" + std::to_string(i), 'I', ""});
    return v;
}

// TMVA cannot take a NaN, so a NaN input goes in as kTmvaMissing and every row is kept (r2D is NaN with fewer than 3 hits, resHalfRatioD with 2 or fewer, zdcKEAllD without a ZDC total).
// Dropping those rows would throw away every short track and bias the sample; the BDT learns from the missing value instead. TmvaSplit.cpp and TmvaApply.cpp both use this value.
static constexpr double kTmvaMissing = -999.0;
inline double TmvaValue(double x){return std::isfinite(x) ? x : kTmvaMissing;}
