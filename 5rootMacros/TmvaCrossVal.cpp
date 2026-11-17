// K-fold cross validation (TMVA::CrossValidation, ROOT 6.14+) of the signal vs background BDT of TmvaTrain.cpp, on a subsample of the files TmvaSplit.cpp wrote, with the same inputs (TmvaVars.h) and BDT options.
// Run in the run folder (build/root/<run>/): root -l -b -q 'TmvaCrossVal.cpp'  or  'TmvaCrossVal.cpp(5, 200000)'
// Every signal row is kept up to maxSignal, the background is thinned to about maxBackground rows by event number (all the tracks of an event stay together, so a fold never holds half an event of another fold).
// With dropVar = an input name the cross validation trains without that input (the drop-one test: the ROC loss of each input; "none" is the baseline), writes only crossval_drop/<dropVar>.txt, and reuses the subsample files if they exist (many runs side by side).
// With fromTmva = a resultTMVAtraining_SigBkg.root (TmvaTrain.cpp) the rows are the training and test events stored in it (dataset/TrainTree and TestTree, all of their inputs, thinned by entry), written to crossval_oldfile/ only,
// for a trained model whose split files are gone: root -l -b -q 'TmvaCrossVal.cpp(5, 200000, 20000, "resultTMVAtraining_SigBkg.root")'
// The subsample files cuts/cut_tmva_cv_signal.root and cut_tmva_cv_background.root are rebuilt each run. A cross validation trains folds times, so the subsample keeps it to minutes instead of the hours the full split files would take.
// Writes resultTMVACrossValidation.root and, for the plots of plots/10_BDTcuts/TMVA/: 6_crossvalidation_ROC.pdf (the ROC of every fold and the average) and crossvalidation_summary.txt (ROC integral and separation per fold, mean and standard deviation).
// The folds train in parallel, one process each (NumWorkerProcs = folds), so the time is about that of one fold.
// Like TmvaTrain.cpp, a PBS job needs `module load openblas/0.3.26-gcc-12.3.0` after `module load root`.
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"

#include "TMVA/CrossValidation.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include "TmvaVars.h"

namespace {
// a copy of the rows of `tree` with fEvent % every == 0 (every one for every <= 1), single threaded: it is a pass over two columns
ULong64_t Thin(const char* inFile, const char* treeName, const char* outFile, ULong64_t every, const std::vector<std::string>& columns){
    ROOT::RDataFrame df(treeName, inFile);
    auto kept = every > 1 ? df.Filter(Form("fEvent %% %llu == 0", every)) : df.Filter("true");
    auto snap = kept.Snapshot(treeName, outFile, columns);
    return *snap->Count();
}
}

void TmvaCrossVal(int folds = 5, ULong64_t maxBackground = 200000, ULong64_t maxSignal = 20000, const char* fromTmva = "", const char* dropVar = "", const char* bdtLabel = "", const char* bdtOpts = ""){
    const bool old = fromTmva && fromTmva[0] != 0;
    const std::string drop = dropVar ? dropVar : ""; // an input left out of the training (drop-one test); "none" leaves none out, the baseline of the test
    const std::string label = bdtLabel ? bdtLabel : ""; // bdtLabel names a run with other BDT options (bdtOpts), written to crossval_opts/<label>.txt only
    const std::string tag = (old ? "Old" : "") + (drop.empty() ? std::string("") : "_drop_" + drop) + (label.empty() ? std::string("") : "_opt_" + label); // the names of a run on an old file or of a drop-one run never collide with another run
    const char* signalFile = "cuts/cut_tmva_signal.root";
    const char* backgroundFile = "cuts/cut_tmva_background.root";
    if (folds < 2) folds = 2;
    std::vector<TmvaVar> vars;
    std::vector<std::string> columns;
    const std::string cvSig = "cuts/cut_tmva_cv" + (old ? std::string("old") : std::string("")) + "_signal.root", cvBkg = "cuts/cut_tmva_cv" + (old ? std::string("old") : std::string("")) + "_background.root";
    gSystem->mkdir("cuts", true);
    ULong64_t nSigUsed = 0, nBkgUsed = 0;
    if (old){
        // the rows of the TMVA file: the training and the test tree, classID 0 signal and 1 background; every input is a float there
        const bool reuse = !drop.empty() && !gSystem->AccessPathName(cvSig.c_str()) && !gSystem->AccessPathName(cvBkg.c_str()); // a drop-one run uses the subsample of the runs before it
        if (!reuse && gSystem->AccessPathName(fromTmva)){Error("TmvaCrossVal", "no %s", fromTmva); return;}
        if (reuse){
            TFile fr(cvSig.c_str()); TTree* tr = (TTree*)fr.Get("TmvaSignal");
            if (!tr){Error("TmvaCrossVal", "no TmvaSignal in %s", cvSig.c_str()); return;}
            for (auto* b : *tr->GetListOfBranches()){vars.push_back({b->GetName(), 'F', ""}); columns.push_back(b->GetName());}
            nSigUsed = tr->GetEntries();
            TFile fb(cvBkg.c_str()); TTree* tb = (TTree*)fb.Get("TmvaBackground");
            nBkgUsed = tb ? tb->GetEntries() : 0;
        } else {
        TChain chain("dataset/TrainTree"); chain.Add(fromTmva); chain.AddFile(fromTmva, TTree::kMaxEntries, "dataset/TestTree");
        for (auto* b : *chain.GetListOfBranches()){
            const std::string n = b->GetName();
            if (n == "classID" || n == "className" || n == "weight" || n == "BDT") continue;
            vars.push_back({n, 'F', ""});
            columns.push_back(n);
        }
        const ULong64_t nSig = chain.GetEntries("classID == 0"), nBkg = chain.GetEntries("classID == 1");
        // every k-th row of a class by entry number (no event number in the TMVA file), single threaded
        const ULong64_t kSig = nSig > maxSignal ? nSig / maxSignal : 1, kBkg = nBkg > maxBackground ? nBkg / maxBackground : 1;
        ROOT::RDataFrame df(chain);
        auto sig = df.Filter(Form("classID == 0 && rdfentry_ %% %llu == 0", kSig)).Snapshot("TmvaSignal", cvSig, columns);
        auto bkg = df.Filter(Form("classID == 1 && rdfentry_ %% %llu == 0", kBkg)).Snapshot("TmvaBackground", cvBkg, columns);
        nSigUsed = *sig->Count(); nBkgUsed = *bkg->Count();
        std::cout << "==> rows of " << fromTmva << ": " << nSigUsed << " signal (every " << kSig << "th of " << nSig << ") and " << nBkgUsed << " background (every " << kBkg << "th of " << nBkg << "), " << folds << " folds, inputs " << columns.size() << std::endl;
        }
    } else {
        for (const char* f : {signalFile, backgroundFile}) if (gSystem->AccessPathName(f)){Error("TmvaCrossVal", "no %s (run TmvaSplit.cpp first)", f); return;}
        vars = TmvaVariables();
        columns = {"fEvent", "trackNo"};
        for (const auto& x : vars) columns.push_back(x.name);
        ULong64_t nSig = 0, nBkg = 0;
        {
            TFile fs(signalFile), fb(backgroundFile);
            TTree* ts = (TTree*)fs.Get("TmvaSignal"); TTree* tb = (TTree*)fb.Get("TmvaBackground");
            if (!ts || !tb){Error("TmvaCrossVal", "no TmvaSignal or TmvaBackground in the split files"); return;}
            for (const auto& x : vars) if (!ts->GetBranch(x.name.c_str()) || !tb->GetBranch(x.name.c_str())){Error("TmvaCrossVal", "no %s in the split files (run TmvaSplit.cpp of this version)", x.name.c_str()); return;}
            nSig = ts->GetEntries(); nBkg = tb->GetEntries();
        }
        // the step of the event thinning, odd so it does not line up with the even events TmvaSplit.cpp kept (prescale 2)
        const auto step = [](ULong64_t n, ULong64_t most){ULong64_t k = most > 0 ? n / most : 1; if (k > 1 && k % 2 == 0) k++; return k < 1 ? (ULong64_t)1 : k;};
        const ULong64_t kSig = step(nSig, maxSignal), kBkg = step(nBkg, maxBackground);
        nSigUsed = Thin(signalFile, "TmvaSignal", cvSig.c_str(), kSig, columns);
        nBkgUsed = Thin(backgroundFile, "TmvaBackground", cvBkg.c_str(), kBkg, columns);
        std::cout << "==> subsample: " << nSigUsed << " signal rows (every " << kSig << "th event) and " << nBkgUsed << " background rows (every " << kBkg << "th event), " << folds << " folds" << std::endl;
    }
    if (nSigUsed < (ULong64_t)folds * 10 || nBkgUsed < (ULong64_t)folds * 10){Error("TmvaCrossVal", "too few rows in the subsample for %d folds", folds); return;}

    TFile* inSignal = TFile::Open(cvSig.c_str());
    TFile* inBackground = TFile::Open(cvBkg.c_str());
    TTree* signalTree = (TTree*)inSignal->Get("TmvaSignal");
    TTree* backgroundTree = (TTree*)inBackground->Get("TmvaBackground");

    TMVA::Tools::Instance();
    TFile* out = TFile::Open(("resultTMVACrossValidation" + tag + ".root").c_str(), "RECREATE");
    TMVA::DataLoader dataloader(("datasetCV" + tag).c_str());
    std::vector<TmvaVar> used;
    for (const auto& x : vars) if (x.name != drop){dataloader.AddVariable(x.name.c_str(), x.name.c_str(), x.unit.c_str(), x.type); used.push_back(x);}
    if (!drop.empty() && used.size() == vars.size() && drop != "none"){Error("TmvaCrossVal", "no input %s to drop", drop.c_str()); return;}
    dataloader.AddSignalTree(signalTree, 1.0);
    dataloader.AddBackgroundTree(backgroundTree, 1.0);
    const TCut noCut = "";
    dataloader.PrepareTrainingAndTestTree(noCut, noCut, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

    TMVA::gConfig().EnableMT(NThreads());
    // the same options as TmvaTrain.cpp, one fold's training set is (folds - 1) / folds of the subsample
    TMVA::CrossValidation cv(("SigBkgCV" + tag).c_str(), &dataloader, out, Form("!V:!Silent:Transformations=I:AnalysisType=Classification:NumFolds=%d:NumWorkerProcs=%d:SplitType=Random:FoldFileOutput=False", folds, folds)); // NumWorkerProcs: one process per fold, the folds train at the same time (one BDT training is single-threaded)
    cv.BookMethod(TMVA::Types::kBDT, "BDT",
                  (bdtOpts && bdtOpts[0]) ? bdtOpts : "!H:!V:NTrees=850:MinNodeSize=10%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
    cv.Evaluate();

    // the pdf and the summary: crossval/ (plots/10_BDTcuts/TMVA/ too for the split files; F9_PlotReconNew.cpp's job removes plots/10_BDTcuts when it starts, so crossval/ keeps a copy); crossval_oldfile/ only for a TMVA file
    const std::string outDir = !label.empty() ? "crossval_opts" : !drop.empty() ? "crossval_drop" : old ? "crossval_oldfile" : "crossval", tmvaDir = "plots/10_BDTcuts/TMVA";
    gSystem->mkdir(outDir.c_str(), true);
    if (!old && drop.empty() && label.empty()) gSystem->mkdir(tmvaDir.c_str(), true);
    const auto& results = cv.GetResults();
    {
        std::ofstream txt(outDir + (!label.empty() ? "/" + label + ".txt" : drop.empty() ? std::string("/crossvalidation_summary.txt") : "/" + drop + ".txt"));
        txt << "TMVA::CrossValidation of the BDT, " << folds << " folds, " << nSigUsed << " signal and " << nBkgUsed << " background rows (" << (old ? std::string("training and test rows of ") + fromTmva : std::string("a subsample of cuts/cut_tmva_*.root")) << "), inputs:";
        for (const auto& x : used) txt << " " << x.name;
        txt << "\ndropped: " << (drop.empty() ? std::string("none") : drop) << "\n";
        for (const auto& r : results){
            txt << "ROC integral per fold:";
            for (const auto& kv : r.GetROCValues()) txt << " " << kv.second; // fold number -> ROC integral
            txt << "\nmean " << r.GetROCAverage() << ", standard deviation " << r.GetROCStandardDeviation() << "\n";
        }
    }
    if (!results.empty()){
        TCanvas* c = results[0].Draw(); // the ROC curve of every fold
        if (c && drop.empty() && label.empty()){c->SaveAs((outDir + "/6_crossvalidation_ROC.pdf").c_str()); if (!old) c->SaveAs((tmvaDir + "/6_crossvalidation_ROC.pdf").c_str());}
    }
    if (!old && drop.empty() && label.empty()) gSystem->CopyFile((outDir + "/crossvalidation_summary.txt").c_str(), (tmvaDir + "/crossvalidation_summary.txt").c_str(), true);
    cv.Print();
    out->Close();
    inSignal->Close();
    inBackground->Close();
    std::cout << "==> Wrote resultTMVACrossValidation" << tag << ".root and " << outDir << "/6_crossvalidation_ROC.pdf, crossvalidation_summary.txt" << std::endl;
}
