#include <ROOT/RDataFrame.hxx>
#include <TString.h>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
void Drop(double BEAMKE, int MAXP){
    ROOT::EnableImplicitMT(32);
    ROOT::RDataFrame d0("AnglesRaw", "cuts/cut_angles_all.root");
    auto df = d0.Define("thDeg", "thetaRelND * 180.0 / 3.14159265358979");
    std::vector<std::pair<std::string, std::string>> cuts = {
        {"KE>90%", Form("zdcKEAllD > 0.9 * %f", BEAMKE)},
        {"even planes", Form("nPlanesUsed %% 2 == 0 && nPlanesUsed <= %d", MAXP)},
        {"theta<=1deg", "thDeg >= 0 && thDeg <= 1.0"},
        {"avg hits/plane=1", "nPlanesHitEv > 0 && nTrackerHitsEv == nPlanesHitEv"},
        {"single track", "trackNo == 0 && nTrackerTracksEv == 1"},
        {"planes>=3", "nPlanesUsed >= 3"}};
    struct R { std::string name; ROOT::RDF::RResultPtr<ULong64_t> all, ce; };
    std::vector<R> res;
    for (int skip = -2; skip < (int)cuts.size(); skip++){
        std::string f = "true"; // skip -2: no cuts at all, skip -1: all cuts, else drop cut number skip
        if (skip != -2) for (int i = 0; i < (int)cuts.size(); i++) if (i != skip) f += " && (" + cuts[i].second + ")";
        auto s = df.Filter(f);
        res.push_back({skip == -2 ? "no cuts" : skip == -1 ? "all cuts" : "drop " + cuts[skip].first, s.Count(), s.Filter("evClass == 0").Count()});
    }
    for (auto& r : res) printf("%-24s CE %8llu all %10llu purity %.3f\n", r.name.c_str(), (unsigned long long)*r.ce, (unsigned long long)*r.all, *r.all ? (double)*r.ce / *r.all : 0.0);
}
