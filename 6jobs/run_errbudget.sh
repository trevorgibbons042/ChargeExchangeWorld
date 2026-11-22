#source ~/.zshrc
# The error budget chain of one run folder: the recon and the cut files, then PErrorBudget.cpp (plots/11_errorBudget/*.txt). No plots, no pdfs.
# Submit: qsub -v RUN=<run folder in build/root>,SCRIPT=run_errbudget.sh Analysis.pbs
source ./run_lib.sh

NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
sh preflight.sh || exit 1
export PREFLIGHT_DONE=1

mkdir -p cuts
export ANA_NTHREADS=$NCPU
run 'ReconAngleAll.cpp'
run 'CutPlots.cpp("events")'
run 'F8_TrackerFits.cpp' # plots/8_trackerfits/track_purity.csv, the counts of Dtrack
run 'PErrorBudget.cpp'
unset ANA_NTHREADS
