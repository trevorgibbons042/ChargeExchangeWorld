#pragma once
#include <cmath>
#include <cstdlib>
#include <string>

// Block residual cut. Scattered protons, and CE protons that scatter afterwards, do not follow a straight line from the beam axis:
//   residual (mm): RMS transverse distance of the track's hits to its fitted line (detector block position)
// There is no cut on the fitted tracks as a whole: AnglesRaw keeps every fitted proton and the Angles view selects all of them (CutFiles.h).
// The block residual cut only applies to block quantities: every tree keeps
// resDRms (ReconAngleAll.cpp, CalcAPandT.cpp), and CutPlots.cpp (trees) cuts it with kBlockCut into the block entry lists the plots read (CutFiles.h),
// so a new threshold needs no rerun of ReconAngleAll.cpp or CalcAPandT.cpp. (The passResD branch holds the default threshold only.)
// The env var RES_CUT_MM (set by res_scan.sh) overrides the block residual threshold; "none" switches it off.
static double ResidualCut(double fallback)
{
    const char* env = std::getenv("RES_CUT_MM");
    if (!env || !*env) return fallback;
    return std::string(env) == "none" ? 1e9 : std::atof(env);
}
static const double kMaxResidualDMm = ResidualCut(2.0);
static const std::string kBlockCut = "resDRms <= " + std::to_string(kMaxResidualDMm); // TTree::Draw cut for every block position quantity
// The 2 mm block residual limit of the 3_ResCuts cut (CutPlots.cpp (samples), Plot3D.cpp).
// resDRms is a branch of AnglesRaw (ReconAngleAll.cpp) and of every cut_APandT.root tree (CalcAPandT.cpp).
static const double kRes2Mm = 2.0;

inline bool PassesBlockResidual(double resD){return resD <= kMaxResidualDMm;}
