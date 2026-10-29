#pragma once
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <TMath.h>

// Theta limit of each carbon block (deg): a proton leaving the block center on the beam axis still lands in the ZDC (its corner, the largest radius) below it.
// Geometry as in 2construction.cc: elements along z at dz spacing around zCenter, nBlocksTracking tracker planes before each carbon block and after the last one.
// Shared by Plot3D.cpp (the 10_thetacut folder of 6_cuts) and SecondaryCE.cpp (the same cut in 7_secondaryCE).
struct BlockLimit {double zMm, thetaMaxDeg;};
inline std::vector<BlockLimit> gBlockLimits;
inline void FillBlockLimits(const std::map<std::string, double>& st){
    const int nBlocks = (int)st.at("nBlocks"), nTrk = (int)st.at("nBlocksTracking");
    const double dz = st.at("dz"), zCenter = 4.0, zZdc = 8.0, zdcHalf = 0.2; // m, the ZDC (2construction.cc) is 0.4 m x 0.4 m at z = 8 m
    const int nElements = (nBlocks + 1) * nTrk + nBlocks;
    const double zStart = zCenter - (nElements - 1) * dz / 2.0, corner = zdcHalf * std::sqrt(2.0);
    gBlockLimits.clear();
    for (int i = 0; i < nBlocks; i++){
        const double z = zStart + ((i + 1) * nTrk + i) * dz;
        gBlockLimits.push_back({z * 1000.0, std::atan(corner / (zZdc - z)) * TMath::RadToDeg()});
    }
}
// The block is the carbon block nearest the recon vertex z (mm); true if the recon block theta (deg) is within that block's limit.
inline bool BlockThetaOk(double thetaBlockDeg, double zVertexMm){
    if (gBlockLimits.empty() || !std::isfinite(thetaBlockDeg) || !std::isfinite(zVertexMm)) return false;
    size_t best = 0;
    for (size_t b = 1; b < gBlockLimits.size(); b++) if (std::fabs(gBlockLimits[b].zMm - zVertexMm) < std::fabs(gBlockLimits[best].zMm - zVertexMm)) best = b;
    return thetaBlockDeg <= gBlockLimits[best].thetaMaxDeg;
}
