#pragma once
#include <cmath>
#include <cstdlib>
#include <string>

// Track quality cuts. Scattered protons, and CE protons that scattered afterwards, do not follow a straight line from the beam axis:
//   residual (mm): RMS transverse distance of the track's hits to its fitted line (true position, and the detector block position)
//   impact (mm):   closest approach of the fitted line to the beam axis (a CE proton starts on the axis)
// Thresholds from QualityStudy.cpp on the Physics run (the four together keep 95% of the real CE protons). Set a threshold to 1e9 to switch a cut off.
// PassesQuality (the impact cuts; the true-position residual cut is off) is applied in ReconAngleAll.cpp before a track is written to Angles,
// AnglesCE, Angles_DE and AnglesCE_DE, so it holds for every plot. The block residual cut only applies to block quantities: every tree keeps
// resDRms (ReconAngleAll.cpp, CalcAPandT.cpp), and CutTrees.cpp cuts it with kBlockCut into the block entry lists the plots read (CutFiles.h),
// so a new threshold needs no rerun of ReconAngleAll.cpp or CalcAPandT.cpp. (The passResD branch holds the default threshold only.)
// The env var RES_CUT_MM (set by res_scan.sh) overrides the block residual threshold; "none" switches it off.
static double ResidualCut(double fallback)
{
    const char* env = std::getenv("RES_CUT_MM");
    if (!env || !*env) return fallback;
    return std::string(env) == "none" ? 1e9 : std::atof(env);
}
static const double kMaxResidualMm = 1e9, kMaxResidualDMm = ResidualCut(2.0);
static const double kMaxImpactMm = 16.2162, kMaxImpactDMm = 16.5846;
static const std::string kBlockCut = "resDRms <= " + std::to_string(kMaxResidualDMm); // TTree::Draw cut for every block position quantity
// Every KE cut also comes as a twin with this tighter block residual cut on top (folder <ke>_res2mm), on every plot that has KE cuts.
// resDRms is a branch of the Angles trees (ReconAngleAll.cpp) and of every cut_APandT.root tree (CalcAPandT.cpp).
static const double kRes2Mm = 2.0;
static const char* const kRes2Cut = "resDRms < 2";
static const char* const kRes2Dir = "_res2mm";
static const char* const kRes2Label = ", block residual < 2 mm";

inline bool PassesQuality(double res, double impact, double impactD)
{
    return res <= kMaxResidualMm && impact <= kMaxImpactMm && impactD <= kMaxImpactDMm;
}
inline bool PassesBlockResidual(double resD){return resD <= kMaxResidualDMm;}
