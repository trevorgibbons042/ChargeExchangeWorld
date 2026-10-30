#!/bin/bash
# Submits the plotting chain for one run folder as PBS jobs (Plot.pbs), the lanes of 5rootMacros/plot.sh in parallel, instead of running plot.sh in one go.
# Usage: bash SubmitPlots.sh <folder in build/root>, e.g. bash SubmitPlots.sh Physics_xs1e10
# Needs cut_*.root and settings.csv already in the run folder (CutEvents.cpp, ReconAngleAll.cpp, CalcAPandT.cpp, CountCE.cpp done).
RUN=$1
[ -n "$RUN" ] || { echo "usage: bash SubmitPlots.sh <run folder>"; exit 1; }
cd "$(dirname "$0")" || exit 1
rsync -a --exclude __pycache__ ../5rootMacros/ "../build/root/$RUN/" || exit 1 # the jobs run the run folder's copy of the macros

# submit <job name> <step> [job ids it waits for] [ncpus, default the Plot.pbs line]: prints the new job id
submit() {
    local name=$1 step=$2 after=$3 ncpus=$4
    local dep=() res=()
    [ -n "$after" ] && dep=(-W depend=afterok:$after)
    [ -n "$ncpus" ] && res=(-l select=1:ncpus=$ncpus)
    qsub -N "$name" -v RUN=$RUN,STEP=$step "${dep[@]}" "${res[@]}" Plot.pbs
}

# CutTrees and the caches are on the critical path, ahead of every lane: all the cores
cut=$(submit CutTrees cuttrees "" 128)
# the cut caches, after CutTrees (a lane that asked for one that is missing would build it itself, under a lock)
caches=$(submit Caches caches "$cut" 128)
# the lanes, after CutTrees writes cuts/<ke>.root and the caches are built
cuts=$(submit Cuts cuts "$caches")
fom=$(submit Fom fom "$caches")
ce=$(submit CEblock ceblock "$caches")
sec=$(submit SecCE secce "$caches")
# tidy, side box and crop once every lane is done
submit Pdfs pdfs "$cuts:$fom:$ce:$sec"
