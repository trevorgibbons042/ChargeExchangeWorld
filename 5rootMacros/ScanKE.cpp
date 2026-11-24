#include <ROOT/RDataFrame.hxx>
#include <TString.h>
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>
void ScanKE(double BEAMKE, int MAXP){
    ROOT::EnableImplicitMT(32);
    ROOT::RDataFrame df("AnglesRaw", "cuts/cut_angles_all.root");
    const double rad = 180.0 / M_PI;
    auto sel = df.Filter(Form("zdcKEAllD > 0.8 * %f && nPlanesUsed %% 2 == 0 && nPlanesUsed <= %d && nPlanesHitEv > 0 && nTrackerHitsEv == nPlanesHitEv && trackNo == 0 && nTrackerTracksEv == 1", BEAMKE, MAXP))
        .Define("thDeg", [=](double t){return t * rad;}, {"thetaRelND"}).Filter("thDeg >= 0 && thDeg <= 1.0");
    auto ke = sel.Take<double>("zdcKEAllD"); auto np = sel.Take<int>("nPlanesUsed"); auto r2 = sel.Take<double>("r2D");
    auto rs = sel.Take<double>("resDRms"); auto im = sel.Take<double>("impactD"); auto th = sel.Take<double>("thDeg"); auto cl = sel.Take<int>("evClass");
    const auto &K = *ke, &R2 = *r2, &RS = *rs, &IM = *im, &TH = *th; const auto &NP = *np, &CL = *cl;
    const size_t n = K.size(); size_t nce = 0; for (size_t i = 0; i < n; i++) nce += CL[i] == 0;
    printf("base rows %zu CE %zu purity %.3f\n", n, nce, n ? (double)nce / n : 0.0);
    for (double kf : {0.80, 0.85, 0.90, 0.92, 0.95}){
        size_t a = 0, c = 0;
        for (size_t i = 0; i < n; i++){ if (K[i] <= kf * BEAMKE) continue; a++; c += CL[i] == 0; }
        printf("KE>%.2f: CE %zu (%.1f%% of base CE) all %zu purity %.3f\n", kf, c, nce ? 100.0 * c / nce : 0.0, a, a ? (double)c / a : 0.0);
    }
}
