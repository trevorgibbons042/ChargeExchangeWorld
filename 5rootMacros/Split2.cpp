#include <ROOT/RDataFrame.hxx>
#include <cmath>
#include <cstdio>
void Split2(){
    ROOT::EnableImplicitMT(32);
    ROOT::RDataFrame df("AnglesRaw", "cuts/cut_angles_all.root");
    const double rad = 180.0 / M_PI;
    auto sel = df.Filter("zdcKEAllD > 0.8 * 50000.0 && nPlanesUsed % 2 == 0 && nPlanesUsed <= 12 && nPlanesHitEv > 0 && nTrackerHitsEv == nPlanesHitEv").Filter([=](double th){return th * rad >= 0.0 && th * rad <= 1.0;}, {"thetaRelND"});
    auto n = sel.Count();
    auto t0 = sel.Filter("trackNo == 0").Count();
    auto t1 = sel.Filter("trackNo > 0").Count();
    auto t1ce = sel.Filter("trackNo > 0 && evClass == 0").Count();
    auto t1ot = sel.Filter("trackNo > 0 && evClass > 1").Count();
    auto t0ce = sel.Filter("trackNo == 0 && evClass == 0").Count();
    auto t0ot = sel.Filter("trackNo == 0 && evClass > 1").Count();
    printf("selected %lld trackNo0 %lld (CE %lld other %lld) trackNo>0 %lld (CE %lld other %lld)\n", (long long)*n, (long long)*t0, (long long)*t0ce, (long long)*t0ot, (long long)*t1, (long long)*t1ce, (long long)*t1ot);
}
