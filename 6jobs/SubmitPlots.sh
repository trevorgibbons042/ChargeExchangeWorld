#!/bin/bash
# Submits the plotting chain for one run folder as PBS jobs (Plot.pbs), the lanes of 5rootMacros/plot.sh in parallel, instead of running plot.sh in one go.
# Usage: bash SubmitPlots.sh <folder in build/root> [job ids CutTrees waits for, colon separated], e.g. bash SubmitPlots.sh Physics_xs1e10
# Needs cut_*.root and settings.csv already in the run folder (CutEvents.cpp, ReconAngleAll.cpp, CalcAPandT.cpp, CountCE.cpp done), or the job ids of the jobs that write them (SubmitEvents.sh).
RUN=$1
[ -n "$RUN" ] || { echo "usage: bash SubmitPlots.sh <run folder> [job ids to wait for]"; exit 1; }
cd "$(dirname "$0")" || exit 1
source ./submit_lib.sh
rsync -a --exclude __pycache__ ../5rootMacros/ "../build/root/$RUN/" || exit 1 # the jobs run the run folder's copy of the macros

# CutTrees and ZdcColumns are on the critical path, ahead of every lane: all the cores
cut=$(submit CutTrees cuttrees "$2" 128)
# the ZDC proton cut files with the Plot3D columns (cut_ZdcCE.root, cut_ZdcAllEv.root), after CutTrees
zdc=$(submit ZdcCols zdc "$cut" 128)
# the lanes, after CutTrees writes cuts/<ke>.root and ZdcColumns writes its cut files
cuts=$(submit Cuts cuts "$zdc")
fom=$(submit Fom fom "$zdc")
ce=$(submit CEblock ceblock "$zdc")
sec=$(submit SecCE secce "$zdc")
# tidy, side box and crop once every lane is done
submit Pdfs pdfs "$cuts:$fom:$ce:$sec"
