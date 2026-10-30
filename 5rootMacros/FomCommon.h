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
static const std::vector<std::string> kSets = {"1_reconTheta", "2_reconBlockTheta"};
// the KE cut (CutFiles.h); its csv file comes from PlotAngle.cpp
static const std::vector<std::string> kKes = {"ke_gt80"}; // the FOM plots have the 80% KE cut and no residual cut folders
static const int kColors[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta, kOrange + 7, kCyan + 2, kGray + 2};
static const std::map<std::string, std::string> kSetName = {{"1_reconTheta", "true position"}, {"2_reconBlockTheta", "block position"}};
static std::string KeName(const std::string& ke){
    const bool res2 = ke.size() > 7 && ke.compare(ke.size() - 7, 7, "_res2mm") == 0;
    const std::string base = res2 ? ke.substr(0, ke.size() - 7) : ke;
    return ("KE cut " + base.substr(5) + "%") + (res2 ? ", res < 2 mm" : "");
}

// Base folder of the FOM and polarization plots, under it <theta set>/<ke>/... (PlotFomPolarization.cpp)
static std::string FomRoot(){return "plots/5_fomCalc";}
// true when only the block position set is wanted (the position set has no fom_polarization plots)
static bool FomBlockOnly(){return false;}
// Subfolder of <theta set>/<ke>/ for the P maps and fits of PlotFomPolarization.cpp
static const std::string kFitDir = ""; // the phi fits go straight into the theta set folder

// Every theta axis runs 0 to the last csv bin edge (the total row's theta_hi of ke_gt80.csv), the same for all plots.
static double ReadThetaMax(){
    std::ifstream in("plots/5_fomCalc/1_reconTheta/1_csv/ke_gt80.csv");
    std::string line, last;
    while (std::getline(in, line)) if (!line.empty()) last = line;
    std::stringstream ss(last); std::string f; std::getline(ss, f, ','); std::getline(ss, f, ',');
    return f.empty() ? 0.0 : std::stod(f);
}

