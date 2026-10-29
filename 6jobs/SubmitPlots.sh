#!/bin/bash
# Submits the plotting chain for one run folder as PBS jobs (Plot.pbs), the lanes of 5rootMacros/plot.sh in parallel, instead of running plot.sh in one go.
# Usage: bash SubmitPlots.sh <folder in build/root>, e.g. bash SubmitPlots.sh Physics_xs1e10
# Needs cut_*.root and settings.csv already in the run folder (CutEvents.cpp, ReconAngleAll.cpp, CalcAPandT.cpp, CountCE.cpp done).
RUN=$1
[ -n "$RUN" ] || { echo "usage: bash SubmitPlots.sh <run folder>"; exit 1; }
cd "$(dirname "$0")" || exit 1
rsync -a --exclude __pycache__ ../5rootMacros/ "../build/root/$RUN/" || exit 1 # the jobs run the run folder's copy of the macros

# submit <job name> <step> [job ids it waits for]: prints the new job id
submit() {
    local name=$1 step=$2 after=$3
    local dep=()
    [ -n "$after" ] && dep=(-W depend=afterok:$after)
    qsub -N "$name" -v RUN=$RUN,STEP=$step "${dep[@]}" Plot.pbs
}

cut=$(submit CutTrees cuttrees)
# the lanes, after CutTrees writes cuts/<ke>.root
cuts=$(submit Cuts cuts "$cut")
fom=$(submit Fom fom "$cut")
ce=$(submit CEblock ceblock "$cut")
sec=$(submit SecCE secce "$cut")
# tidy, side box and crop once every lane is done
submit Pdfs pdfs "$cuts:$fom:$ce:$sec"
