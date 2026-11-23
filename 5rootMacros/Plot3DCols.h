#pragma once
// The Plot3D values of F1_2_Plot3D.cpp, by name, from the columns of AnglesRaw (the event level ones are repeated on each track row, ReconAngleAll.cpp): no file and no map of their own.
// The cut strings and expressions of F1_2_Plot3D.cpp say isReco(fEvent), rawVal(fEvent, k) and rawOk(fEvent, k); P3() rewrites them into the formula of slot k.
// A row of AnglesRaw is a track of an event, so every slot is a value of the row's own track (event level ones are repeated on every track of the event).
// Every row is a reconstructed fit, so isReco is 1.
// Slots: 0 resRms, 1 resDRms, 2 impact, 3 impactD, 4 r2, 5 r2D, 6 nChargedParticle, 7 nZdcCharged (NaN without it), 8 resHalfRatioD, 9 thetaRelN (deg), 10 nHits, 11 thetaRelND (deg), 12 firstXD, 13 firstYD,
// 14 and 15 RECONPHI_ESTNEU and RECONPHIBL_ESTNEU (phiRelN, phiRelND in degrees, 0 to 360), 16 recon vertex z (mm, from the block position line fit), 17 nTrackerHitsEv and 18 nTrackerTracksEv (NaN without them),
// 19 nMaxChargedPerPlaneEv (the most charged tracks that entered one tracker plane in the event), 20 and 21 thetaRelNDMin / thetaRelNDMax (deg) and 22 phiRelNDHalf (rad): the est neutron range of the block recon proton,
// 23 and 24 thetaRelNDBMin / thetaRelNDBMax (deg) and 25 phiRelNDBHalf (rad): the block grid range alone, 26 the average charged tracker plane hits per hit plane, nTrackerHitsEv / nPlanesHitEv (NaN without a hit plane or without the columns).
#include <cstdio>
#include <regex>
#include <string>
#include <unordered_map>


// The formula of slot k, over the columns of AnglesRaw and Events. TTree::Draw and RDataFrame strings both read it (TMath::QuietNaN() and TMath::RadToDeg() are in both).
inline std::string P3Formula(int k){
    static constexpr const char* kDeg = "TMath::RadToDeg()";
    static constexpr const char* kNaN = "TMath::QuietNaN()";
    char b[256];
    switch (k){
        case 0: return "resRms"; case 1: return "resDRms"; case 2: return "impact"; case 3: return "impactD"; case 4: return "r2"; case 5: return "r2D";
        case 6: return "nChargedParticle";
        case 7: snprintf(b, sizeof b, "(nZdcCharged>=0?nZdcCharged:%s)", kNaN); return b;
        case 8: return "resHalfRatioD";
        case 9: snprintf(b, sizeof b, "(thetaRelN*%s)", kDeg); return b;
        case 10: return "nHits";
        case 11: snprintf(b, sizeof b, "(thetaRelND*%s)", kDeg); return b;
        case 12: return "firstXD"; case 13: return "firstYD";
        case 14: snprintf(b, sizeof b, "((phiRelN<0?phiRelN+TMath::TwoPi():phiRelN)*%s)", kDeg); return b;
        case 15: snprintf(b, sizeof b, "((phiRelND<0?phiRelND+TMath::TwoPi():phiRelND)*%s)", kDeg); return b;
        case 16: return "z_bx_reconVertexOut";
        case 17: snprintf(b, sizeof b, "(nTrackerHitsEv>=0?nTrackerHitsEv:%s)", kNaN); return b;
        case 18: snprintf(b, sizeof b, "(nTrackerTracksEv>=0?nTrackerTracksEv:%s)", kNaN); return b;
        case 19: snprintf(b, sizeof b, "(nMaxChargedPerPlaneEv>=0?nMaxChargedPerPlaneEv:%s)", kNaN); return b;
        case 20: snprintf(b, sizeof b, "(thetaRelNDMin*%s)", kDeg); return b;
        case 21: snprintf(b, sizeof b, "(thetaRelNDMax*%s)", kDeg); return b;
        case 22: return "phiRelNDHalf";
        case 23: snprintf(b, sizeof b, "(thetaRelNDBMin*%s)", kDeg); return b;
        case 24: snprintf(b, sizeof b, "(thetaRelNDBMax*%s)", kDeg); return b;
        case 25: return "phiRelNDBHalf";
        case 26: snprintf(b, sizeof b, "((nPlanesHitEv>0&&nTrackerHitsEv>=0)?nTrackerHitsEv*1.0/nPlanesHitEv:%s)", kNaN); return b;
        default: return "TMath::QuietNaN()";
    }
}

// Rewrites the isReco/rawVal/rawOk calls of a TTree::Draw expression or cut into the formulas.
inline std::string P3(const std::string& s){
    thread_local std::unordered_map<std::string, std::string> cache; // the same expressions and cuts come back in many calls
    if (const auto hit = cache.find(s); hit != cache.end()) return hit->second;
    static const std::regex fn(R"(\bisReco\s*\(\s*fEvent\s*\))");
    static const std::regex call(R"(\braw(Val|Ok)\s*\(\s*fEvent\s*,\s*(\d+)\s*\))");
    std::string o = std::regex_replace(s, fn, "(1==1)"); // every row of AnglesRaw is a reconstructed fit
    std::string out;
    auto last = o.cbegin();
    for (std::sregex_iterator it(o.begin(), o.end(), call), end; it != end; ++it){
        const auto& m = *it;
        out.append(last, m[0].first);
        const std::string f = P3Formula(std::stoi(m[2].str()));
        out += m[1].str() == "Val" ? f : "(TMath::Finite(" + f + ")==1)";
        last = m[0].second;
    }
    out.append(last, o.cend());
    cache.emplace(s, out);
    return out;
}
