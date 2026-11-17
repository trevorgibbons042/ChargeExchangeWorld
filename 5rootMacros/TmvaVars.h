#pragma once
// The inputs of the signal vs background classifier of TmvaTrain.cpp and TmvaApply.cpp, in one place so training and apply cannot disagree.
// Every one is a column of AnglesRaw (ReconAngleAll.cpp) that a real detector also measures: the pixelated (D) positions and fit, the smeared ZDC KE and hit counts. No truth column and no label (inCE, isFromCE, evClass, nCEUsed...).
// Signal = a row whose fit used a CE child hit (isFromCE == 1), background = isFromCE == 0.
#include <cmath>
#include <string>
#include <vector>
#include "FitCore.h"

struct TmvaVar {std::string name; char type; std::string unit;}; // type 'D' double or 'I' int, as in the AnglesRaw branch

// The variables of a run: the block count (StackZ, FitCore.h) decides whether nPlanesUsed is in the list, so a classifier is trained and applied on runs with the same block count.
// Dropped on purpose: outDPhi, resDRms, impactD, firstXD, firstYD, r2D, resHalfRatioD, nHits and nChargedPlane<i>.
inline std::vector<TmvaVar> TmvaVariables(){
    std::vector<TmvaVar> v = {
        {"zdcKEAllD", 'D', "MeV"}, {"outDTheta", 'D', "rad"}, {"lastXD", 'D', "mm"}, {"lastYD", 'D', "mm"}, {"nChargedParticle", 'I', ""}};
    // With 2 planes every fit uses both, so nPlanesUsed is 2 on every row and TMVA stops on a variable with no spread: it joins from 3 planes on.
    if (StackZ(ReadSettings()).planeZ.size() >= 3) v.push_back({"nPlanesUsed", 'I', ""});
    return v;
}

// TMVA cannot take a NaN, so a NaN input goes in as kTmvaMissing and every row is kept (r2D is NaN with fewer than 3 hits, resHalfRatioD with 2 or fewer, zdcKEAllD without a ZDC total).
// Dropping those rows would throw away every short track and bias the sample; the BDT learns from the missing value instead. TmvaSplit.cpp and TmvaApply.cpp both use this value.
static constexpr double kTmvaMissing = -999.0;
inline double TmvaValue(double x){return std::isfinite(x) ? x : kTmvaMissing;}
