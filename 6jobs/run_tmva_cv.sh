# K-fold cross validation of the signal vs background BDT on one run folder (5rootMacros/TmvaCrossVal.cpp). Run by Analysis.pbs with -v RUN=<run>,SCRIPT=run_tmva_cv.sh
# Reuses cuts/cut_tmva_signal.root and cuts/cut_tmva_background.root of an earlier run_tmva.sh. Writes tmva_crossval.log and resultTMVAcrossval_SigBkg.root.

source ./run_lib.sh # run <macro>: times it, ACLiC or interpreted

run 'TmvaCrossVal.cpp' 2>&1 | tee tmva_crossval.log
