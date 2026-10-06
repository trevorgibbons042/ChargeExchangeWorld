#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <ROOT/RDataFrame.hxx>

#include <map>
#include <tuple>
#include <mutex>
#include <limits>
#include <cmath>
#include <iostream>
#include <sstream>
#include <filesystem>


struct EnergyPair {
    double KE3 = NAN;
    double KE2 = NAN;

    ULong64_t entry3 = std::numeric_limits<ULong64_t>::max();
    ULong64_t entry2 = std::numeric_limits<ULong64_t>::max();
};


void EnergyDifference()
{
    std::string folder = "pngs_EnergyDiff";
    std::filesystem::create_directories(folder);
    
    ROOT::EnableImplicitMT();
    bool Verbose = true;

    ROOT::RDataFrame df("ParticleInfo", "cut_EOL.root");
    using Key = std::tuple<int, int, int>;
    std::map<Key, EnergyPair> energies;
    std::mutex m;


    // Collect boundary 3 and boundary 2 energies
    df.Foreach([&](
        int eventID, int trackID, int pdg, int boundary, double KE,ULong64_t entry){
            if (boundary != 3 && boundary != 2) {std::cout << "heyyy" << std::endl; return;}

            std::lock_guard<std::mutex> lock(m);

            auto &e = energies[std::make_tuple(eventID, trackID, pdg)];

            if (boundary == 3 && entry < e.entry3) {
                e.KE3 = KE;
                e.entry3 = entry;
            }
            if (boundary == 2 && entry < e.entry2) {
                e.KE2 = KE;
                e.entry2 = entry;
            }

        },
        {"fEvent","TrackID","pdg","boundaryType","fKEnergy","rdfentry_"}
    );


    // KE3 - KE2
    TH1D *hEnergyDiff = new TH1D(
        "hEnergyDiff",
        "Energy Difference Between Boundary 3 and 2;"
        "KE_{3} - KE_{2} [MeV];Counts",
        100, -10, 10);

    int matchedTracks = 0;
    for (const auto &[key, e] : energies) {

        if (!std::isfinite(e.KE3) ||!std::isfinite(e.KE2)) {continue;}

        double deltaKE = e.KE3 - e.KE2;
        hEnergyDiff->Fill(deltaKE);
        matchedTracks++;

        if (Verbose == true) {
            std::cout
                << "Event: " << std::get<0>(key)
                << " TrackID: " << std::get<1>(key)
                << " PDG: " << std::get<2>(key)
                << " KE3: " << e.KE3
                << " KE2: " << e.KE2
                << " deltaKE: " << deltaKE
                << std::endl;
            
            std::cout << "Tracks with boundary 3 and boundary 2: " << matchedTracks << std::endl;
        }
    }

    TCanvas *canvas = new TCanvas("canvas", "Energy Difference",900,700);

    hEnergyDiff->Draw();
    canvas->Update();
    canvas->SaveAs((folder + "/energy_difference_0_2.png").c_str());
}
