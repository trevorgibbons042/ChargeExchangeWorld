# PErrorBudget.cpp only (plots/11_errorBudget/*.txt and P_table.md), on a run folder that already has cuts/cut_angles_all.root and cuts/cut_CE.root. Submit: qsub -v RUN=<run folder in build/root>,SCRIPT=run_peb.sh Analysis.pbs
source ./run_lib.sh
export ANA_NTHREADS=${NCPUS:-$(nproc)}
run 'PErrorBudget.cpp'
