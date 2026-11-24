#include <ROOT/RDataFrame.hxx>
#include <TString.h>
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>
void Scan(double BEAMKE, int MAXP){
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
    double bestCE = 0, bestP = 0; double bk = 0, br = 0, brs = 0, bim = 0, bth = 0; int bnp = 0; size_t bAll = 0;
    for (double kf = 0.80; kf <= 1.0001; kf += 0.01)
    for (int mp = 2; mp <= 8; mp += 2)
    for (double r2c = 0.0; r2c <= 0.9991; r2c += 0.05)
    for (double rsc : {1e9, 8.0, 6.0, 4.0, 3.0, 2.0, 1.0})
    for (double imc : {1e9, 40.0, 30.0, 20.0, 15.0, 10.0, 5.0})
    for (double thc : {1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2}){
        size_t a = 0, c = 0;
        for (size_t i = 0; i < n; i++){
            if (K[i] <= kf * BEAMKE || NP[i] < mp || R2[i] < r2c || RS[i] > rsc || IM[i] > imc || TH[i] > thc) continue;
            a++; c += CL[i] == 0;
        }
        if (a && (double)c / a >= 0.8 && c > bestCE){bestCE = c; bestP = (double)c / a; bk = kf; bnp = mp; br = r2c; brs = rsc; bim = imc; bth = thc; bAll = a;}
    }
    printf("best CE %.0f of %zu (%.1f%%) purity %.3f all %zu | KE>%.2f beam planes>=%d r2>=%.2f resRms<=%.1f impact<=%.1f theta<=%.1f\n", bestCE, nce, nce ? 100.0 * bestCE / nce : 0.0, bestP, bAll, bk, bnp, br, brs, bim, bth);
}
