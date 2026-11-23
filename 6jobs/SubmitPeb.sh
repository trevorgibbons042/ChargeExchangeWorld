#!/bin/bash
# PErrorBudget.cpp (plots/11_errorBudget/*.txt and P_table.md) for run folders that have cuts/cut_angles_all.root and cuts/cut_CE.root, one job per run folder.
# Usage: bash SubmitPeb.sh <run folder in build/root> [more run folders ...]
[ $# -gt 0 ] || { echo "usage: bash SubmitPeb.sh <run folder in build/root> [more run folders ...]"; exit 1; }
cd "$(dirname "$0")" || exit 1
for RUN in "$@"; do
    RUN=${RUN%/}
    [ -f "../build/root/$RUN/cuts/cut_CE.root" ] || { echo "no cuts/cut_CE.root in $RUN, skipped"; continue; }
    rsync -a --exclude __pycache__ --exclude '*_cpp.*' --exclude .DS_Store --exclude CLAUDE.md --exclude debug ../5rootMacros/ "../build/root/$RUN/" || exit 1
    rsync -a ../6jobs/run_peb.sh ../6jobs/run_lib.sh ../6jobs/preflight.sh "../build/root/$RUN/" || exit 1
    echo "$RUN $(qsub -v RUN=$RUN,SCRIPT=run_peb.sh Analysis.pbs)"
done
