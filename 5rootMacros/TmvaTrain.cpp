// Trains a BDT that separates signal (the fit used a CE child hit, isFromCE == 1) from background (isFromCE == 0), from the inputs of TmvaVars.h.
// Reads the two files TmvaSplit.cpp wrote (every row, a NaN input as -999, prescaled there). Run in the run folder (build/root/<run>/): root -l -b -q 'TmvaTrain.cpp'
// On the cluster a PBS job needs `module load openblas/0.3.26-gcc-12.3.0` after `module load root`: libTMVA.so needs libopenblas.so.0, which the compute nodes do not have in /usr/lib64.
// Writes resultTMVAtraining_SigBkg.root (TMVA plots) and dataset/weights/SigBkg_BDT.weights.xml (read by TmvaApply.cpp). Train and apply on runs with the same block count (TmvaVars.h).
// Adapted from TMVA_training_for_Grace.C: no std::filesystem (cling on the cluster), one method.
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "TCut.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/DataLoader.h"
#include "TMVA/Config.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

#include "TmvaVars.h"

void TmvaTrain(const char* signalFile = "cuts/cut_tmva_signal.root", const char* backgroundFile = "cuts/cut_tmva_background.root"){
    for (const char* f : {signalFile, backgroundFile}) if (gSystem->AccessPathName(f)){Error("TmvaTrain", "no %s (run TmvaSplit.cpp first)", f); return;}
    TFile* inSignal = TFile::Open(signalFile);
    TFile* inBackground = TFile::Open(backgroundFile);
    TTree* signalTree = inSignal ? (TTree*)inSignal->Get("TmvaSignal") : nullptr;
    TTree* backgroundTree = inBackground ? (TTree*)inBackground->Get("TmvaBackground") : nullptr;
    if (!signalTree || !backgroundTree){Error("TmvaTrain", "no TmvaSignal in %s or no TmvaBackground in %s", signalFile, backgroundFile); return;}

    const std::vector<TmvaVar> vars = TmvaVariables();
    for (const auto& x : vars) if (!signalTree->GetBranch(x.name.c_str()) || !backgroundTree->GetBranch(x.name.c_str())){Error("TmvaTrain", "no %s in the split files (run TmvaSplit.cpp of this version)", x.name.c_str()); return;}

    TMVA::Tools::Instance();
    TFile* out = TFile::Open("resultTMVAtraining_SigBkg.root", "RECREATE");
    // Transformations=I only: Grace's I;D;P;G,D builds decorrelated, PCA and Gaussianized copies of the inputs for control plots, the BDT does not use them, and the decorrelation fails on these inputs (integer counts, -999 for missing values) with a NaN table and a FATAL in GetSeparation
    TMVA::Factory factory("SigBkg", out, "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification");
    TMVA::DataLoader dataloader("dataset");
    for (const auto& x : vars) dataloader.AddVariable(x.name.c_str(), x.name.c_str(), x.unit.c_str(), x.type);

    // each class has its own file and tree, no cut here (TmvaSplit.cpp made them)
    dataloader.AddSignalTree(signalTree, 1.0);
    dataloader.AddBackgroundTree(backgroundTree, 1.0);
    const TCut noCut = "";
    // half of each class trains, the other half tests
    dataloader.PrepareTrainingAndTestTree(noCut, noCut, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

    factory.BookMethod(&dataloader, TMVA::Types::kBDT, "BDT",
                       "!H:!V:NTrees=850:MinNodeSize=10%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
    // threads for the parts of TMVA that can use them (BDT training itself is mostly single-threaded, so the time below shows how much it gains)
    TMVA::gConfig().EnableMT(NThreads());
    const auto tTrain = std::chrono::steady_clock::now();
    factory.TrainAllMethods();
    std::cout << "==> training took " << std::chrono::duration<double>(std::chrono::steady_clock::now() - tTrain).count() << " s on " << NThreads() << " threads" << std::endl;
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    out->Close();
    inSignal->Close();
    inBackground->Close();
    std::cout << "==> Wrote resultTMVAtraining_SigBkg.root and dataset/weights/SigBkg_BDT.weights.xml" << std::endl;
}
