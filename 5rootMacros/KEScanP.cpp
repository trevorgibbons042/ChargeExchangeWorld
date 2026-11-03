#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>

#include <cstdio>
#include <string>
#include <vector>

#include "FitCore.h"
#include "Plot3DCols.h"

// P of the phi cos fit at recon block theta start 0.15 deg (2_ReconD/2_phiReD_P_vs_start_theta of Plot3D.cpp part 4) for the 9_megaCut3/1_1particlecut sample with the ZDC KE cut at 80, 85, 90 and 95%,
// for all recon particles, all recon protons and CE recon protons, with the purity (CE recon protons / all recon protons (particles) in the cut) and the share of the CE recon protons
// that stay (CE recon protons in the cut / CE recon protons of cut_angles_all.root with theta 0 to 1 deg). All three counts and the fit use the block theta 0.15 to 1 deg and block phi 0.5 to 359.5 deg,
// like the plot. Cut files: cuts/cut_megacut3_1.root (80%), cuts/cut_megacut3_1_ke85/90/95.root (CutPlots.cpp (samples)). Writes cuts/ke_scan_P.csv. Run: root -l -b -q KEScanP.cpp
namespace {
constexpr int kStartBin = 15; // start theta 0.15 deg: ScanStartTheta bucket 15 (0.01 deg each)

struct Rows {std::vector<double> th, ph;}; // block theta and phi (deg) of the rows that pass

Rows Load(const std::string& file, const std::string& extra){
    ROOT::RDataFrame df("AnglesRaw", file);
    auto d = df.Define("thR", "outTheta * TMath::RadToDeg()").Define("th", "outDTheta * TMath::RadToDeg()").Define("ph", "outDPhi * TMath::RadToDeg()")
               .Filter("thR >= 0 && thR <= 1 && th >= 0.15 && th <= 1 && ph >= 0.5 && ph <= 359.5" + extra); // NaN fails every comparison
    auto th = d.Take<double>("th"), ph = d.Take<double>("ph");
    return {*th, *ph};
}

struct Fit {double p = 0, err = 0, aMean = 0; long long n = 0;};

Fit FitP(const Rows& r, double apScale){
    Fit f;
    f.n = (long long)r.th.size();
    if (f.n < 10) return f;
    const StartScan s = ScanStartTheta(r.th, r.ph, apScale);
    TH1D h("hKEScan", "", 50, 0, 360);
    for (int j = 1; j <= 50; j++) h.SetBinContent(j, s.phi[kStartBin][j]);
    h.SetEntries((double)f.n);
    f.aMean = s.aSum[kStartBin] / f.n;
    TF1 fit("fKEScan", "[0] * (1 + [1] * [2] * cos(x * TMath::Pi() / 180))", 0, 360);
    fit.FixParameter(0, f.n / 50.0); fit.FixParameter(2, f.aMean); fit.SetParameter(1, 0.5); fit.SetParLimits(1, -1, 1);
    h.Fit(&fit, "Q0");
    f.p = fit.GetParameter(1); f.err = fit.GetParError(1);
    return f;
}
}

void KEScanP(){
    ROOT::EnableImplicitMT(NThreads());
    const double apScale = ReadSettings().at("APScale");
    const double nCEtotal = (double)Load("cuts/cut_angles_all.root", " && isFromCE == 1").th.size();
    printf("CE recon protons before any cut (block theta 0.15 to 1 deg): %.0f\n", nCEtotal);
    FILE* cf = fopen("cuts/ke_scan_P.csv", "w");
    fprintf(cf, "ZDC_KE_cut_percent,sample,rows,P,P_err,CE_recon_protons,all_recon_protons,all_recon_particles,purity_CE_over_all_protons,purity_CE_over_all_particles,CE_in_cut_over_CE_total\n");
    printf("%-5s %-22s %9s %8s %8s | %8s %8s %8s | %8s %8s %8s\n", "KE%", "sample", "rows", "P", "err", "CE", "protons", "particles", "CE/prot", "CE/part", "CE/CEtot");
    for (const int pct : {80, 85, 90, 95}){
        const std::string file = pct == 80 ? "cuts/cut_megacut3_1.root" : Form("cuts/cut_megacut3_1_ke%d.root", pct);
        const Rows part = Load(file, ""), prot = Load(file, " && PDG == 2212"), ce = Load(file, " && isFromCE == 1");
        const double nPart = (double)part.th.size(), nProt = (double)prot.th.size(), nCE = (double)ce.th.size();
        const double purProt = nProt > 0 ? nCE / nProt : 0, purPart = nPart > 0 ? nCE / nPart : 0, eff = nCEtotal > 0 ? nCE / nCEtotal : 0;
        const std::pair<const char*, const Rows*> samples[3] = {{"all recon particles", &part}, {"all recon protons", &prot}, {"CE recon protons", &ce}};
        for (const auto& [name, rows] : samples){
            const Fit f = FitP(*rows, apScale);
            printf("%-5d %-22s %9lld %8.4f %8.4f | %8.0f %8.0f %8.0f | %8.4f %8.4f %8.4f\n", pct, name, f.n, f.p, f.err, nCE, nProt, nPart, purProt, purPart, eff);
            fprintf(cf, "%d,%s,%lld,%.5f,%.5f,%.0f,%.0f,%.0f,%.5f,%.5f,%.5f\n", pct, name, f.n, f.p, f.err, nCE, nProt, nPart, purProt, purPart, eff);
        }
    }
    fclose(cf);
}
