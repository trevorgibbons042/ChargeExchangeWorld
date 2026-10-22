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

// Lists, labels and helpers shared by the plots under plots/6_fomCalc (PlotFomCumulative.cpp, PlotCEFraction.cpp).
static const std::vector<std::string> kEnergies = {"1_KEtrue", "2_KEDetector"}, kSets = {"FOM_calc/1_reconTheta", "FOM_calc/2_reconBlockTheta"};
static const std::vector<std::string> kKes = {"ke_all", "ke_gt50", "ke_gt60", "ke_gt70", "ke_gt80", "ke_gt90"};
static const int kColors[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta, kOrange + 7};
static const std::map<std::string, std::string> kEnergyName = {{"1_KEtrue", "true KE"}, {"2_KEDetector", "detector KE"}},
                                                kSetName = {{"FOM_calc/1_reconTheta", "true position"}, {"FOM_calc/2_reconBlockTheta", "block position"}};
static std::string KeName(const std::string& ke){return ke == "ke_all" ? std::string("no KE cut") : "KE cut " + ke.substr(5) + "%";}

// Every theta axis runs 0 to the last csv bin edge (the total row's theta_hi of ke_all.csv), the same for all plots.
static double ReadThetaMax(){
    std::ifstream in("plots/6_fomCalc/1_KEtrue/FOM_calc/1_reconTheta/1_csv/ke_all.csv");
    std::string line, last;
    while (std::getline(in, line)) if (!line.empty()) last = line;
    std::stringstream ss(last); std::string f; std::getline(ss, f, ','); std::getline(ss, f, ',');
    return f.empty() ? 0.0 : std::stod(f);
}

// The plot drawn on c (mg already drawn) again on a log y axis, to logPath. The axis starts at half the smallest positive point.
static void SaveLogCopy(TCanvas& c, TMultiGraph& mg, const std::string& logPath){
    double minPos = 1e300;
    for (TObject* o : *mg.GetListOfGraphs()){
        TGraph* g = static_cast<TGraph*>(o);
        for (int i = 0; i < g->GetN(); i++) if (g->GetY()[i] > 0) minPos = std::min(minPos, g->GetY()[i]);
    }
    if (minPos >= 1e300) return;
    mg.SetMinimum(minPos / 2);
    c.SetLogy(1); c.Modified(); c.Update();
    c.SaveAs(logPath.c_str());
    c.SetLogy(0);
}
