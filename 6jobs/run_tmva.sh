# Signal vs background BDT on one run folder: TmvaSplit, TmvaTrain, TmvaApply (5rootMacros/Tmva*.cpp). Run by Analysis.pbs with -v RUN=<run>,SCRIPT=run_tmva.sh[,PRESCALE=2][,FORCE_RECON=1]
# ReconAngleAll runs only when cuts/cut_angles_all.root is missing (FORCE_RECON=1 reruns it). Every step writes its full output to its own log in the run folder, so `tail -f tmva_train.log` shows the progress.
# Outputs: cuts/cut_tmva_signal.root, cuts/cut_tmva_background.root, resultTMVAtraining_SigBkg.root (TMVA::TMVAGui), dataset/weights/SigBkg_BDT.weights.xml, cuts/tmva_scores.root.

source ./run_lib.sh # run <macro>: times it, ACLiC or interpreted

NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
export ANA_NTHREADS=$NCPU # ReconAngleAll's pool; the TMVA macros read NCPUS through NThreads() (FitCore.h)
mkdir -p cuts
rm -rf cuts/cut_tmva_signal.root cuts/cut_tmva_background.root cuts/tmva_scores.root dataset resultTMVAtraining_SigBkg.root # an old training must not be read by this run

if [ -n "$FORCE_RECON" ] || [ ! -f cuts/cut_angles_all.root ]; then
    run 'ReconAngleAll.cpp' 2>&1 | tee tmva_recon.log
fi
run "TmvaSplit.cpp(\"cuts/cut_angles_all.root\", ${PRESCALE:-2})" 2>&1 | tee tmva_split.log
run 'TmvaTrain.cpp' 2>&1 | tee tmva_train.log
run 'TmvaApply.cpp' 2>&1 | tee tmva_apply.log
unset ANA_NTHREADS
