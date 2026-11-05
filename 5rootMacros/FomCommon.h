#pragma once
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TObject.h>
#include <TStyle.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Lists, labels and helpers shared by the plots under plots/5_fomCalc (PlotFomCumulative.cpp, PlotCEFraction.cpp).
// Functions with a function-local static, not globals with a runtime initializer: ACLiC cannot load those on LCRC (Failed to materialize symbols).
inline const std::vector<std::string>& Sets(){static const std::vector<std::string> v = {"1_reconTheta_estneu", "2_reconBlockTheta_estneu"}; return v;}
// the KE cut (CutFiles.h); its csv file comes from PlotAngle.cpp
inline const std::vector<std::string>& Kes(){static const std::vector<std::string> v = {"ke_gt80"}; return v;} // the FOM plots have the 80% KE cut and no residual cut folders
static const int kColors[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta, kOrange + 7, kCyan + 2, kGray + 2};
inline const std::map<std::string, std::string>& SetNames(){static const std::map<std::string, std::string> m = {{"1_reconTheta_estneu", "true position"}, {"2_reconBlockTheta_estneu", "block position"}}; return m;}

// Base folder of the FOM and polarization plots, under it <theta set>/<ke>/... (PlotFomPolarization.cpp)
static std::string FomRoot(){return "plots/5_fomCalc";}
// true when only the block position set is wanted (the position set has no fom_polarization plots)
static bool FomBlockOnly(){return false;}
// Subfolder of <theta set>/<ke>/ for the P maps and fits of PlotFomPolarization.cpp
static constexpr const char* kFitDir = ""; // the phi fits go straight into the theta set folder

// Every theta axis runs 0 to the last csv bin edge (the total row's theta_hi of ke_gt80.csv), the same for all plots.
static double ReadThetaMax(){
    std::ifstream in("plots/5_fomCalc/1_reconTheta_estneu/1_csv/ke_gt80.csv");
    std::string line, last;
    while (std::getline(in, line)) if (!line.empty()) last = line;
    std::stringstream ss(last); std::string f; std::getline(ss, f, ','); std::getline(ss, f, ',');
    return f.empty() ? 0.0 : std::stod(f);
}

