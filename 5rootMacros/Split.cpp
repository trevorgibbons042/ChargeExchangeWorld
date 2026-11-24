#include <ROOT/RDataFrame.hxx>
#include <cmath>
#include <cstdio>
void Split(){
    ROOT::EnableImplicitMT(32);
    ROOT::RDataFrame df("AnglesRaw", "cuts/cut_angles_all.root");
    const double rad = 180.0 / M_PI;
    auto sel = df.Filter("zdcKEAllD > 0.8 * 50000.0 && nPlanesUsed % 2 == 0 && nPlanesUsed <= 12 && nPlanesHitEv > 0 && nTrackerHitsEv == nPlanesHitEv").Filter([=](double th){return th * rad >= 0.0 && th * rad <= 1.0;}, {"thetaRelND"});
    auto n = sel.Count();
    auto cls = sel.Histo1D({"c", "", 6, -0.5, 5.5}, "evClass");
    auto oth = sel.Filter("evClass > 1");
    auto nOth = oth.Count();
    auto pi = oth.Filter("nPionUsed > 0").Count();
    auto ka = oth.Filter("nKaonUsed > 0").Count();
    auto nu = oth.Filter("nNucleusUsed > 0").Count();
    auto pr = oth.Filter("nProtonUsed > 0 && nPionUsed == 0 && nKaonUsed == 0 && nNucleusUsed == 0").Count();
    auto ceUsed = oth.Filter("nCEUsed > 0").Count();
    auto ceEvPart = oth.Filter("nCEPartEv > 0").Count();
    auto noCEPart = oth.Filter("nCEPartEv == 0").Count();
    auto first = oth.Filter("trackNo == 0").Count();
    auto multi = oth.Filter("nTrackerTracksEv > 1").Count();
    auto planes = oth.Histo1D({"p", "", 14, -0.5, 13.5}, "nPlanesUsed");
    printf("selected %lld\n", (long long)*n);
    for (int b = 1; b <= 6; b++) printf("evClass %d: %d\n", b - 1, (int)cls->GetBinContent(b));
    printf("evClass>1: %lld  withPion %lld withKaon %lld withNucleus %lld protonOnly %lld | CEused>0 %lld | CE particle in event %lld noCEpartInEvent %lld | trackNo0 %lld | multiTrackEv %lld\n",
        (long long)*nOth, (long long)*pi, (long long)*ka, (long long)*nu, (long long)*pr, (long long)*ceUsed, (long long)*ceEvPart, (long long)*noCEPart, (long long)*first, (long long)*multi);
    for (int b = 1; b <= 14; b++) if (planes->GetBinContent(b) > 0) printf("other planes %d: %d\n", b - 1, (int)planes->GetBinContent(b));
}
