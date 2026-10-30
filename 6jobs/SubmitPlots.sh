#!/bin/bash
# Submits the plotting chain for one run folder as PBS jobs (Plot.pbs), the lanes of 5rootMacros/plot.sh in parallel, instead of running plot.sh in one go.
# Usage: bash SubmitPlots.sh <folder in build/root> [job ids CutPlots waits for, colon separated], e.g. bash SubmitPlots.sh Physics_xs1e10
# Needs cut_*.root and settings.csv already in the run folder (CutEvents.cpp, ReconAngleAll.cpp, CalcAPandT.cpp, CountCE.cpp done), or the job ids of the jobs that write them (SubmitEvents.sh).
RUN=$1
[ -n "$RUN" ] || { echo "usage: bash SubmitPlots.sh <run folder> [job ids to wait for]"; exit 1; }
cd "$(dirname "$0")" || exit 1
source ./submit_lib.sh
rsync -a --exclude __pycache__ ../5rootMacros/ "../build/root/$RUN/" || exit 1 # the jobs run the run folder's copy of the macros

# CutPlots (KE cut lists, ZDC proton cut files, sample cut files) is on the critical path, ahead of every lane: all the cores
cut=$(submit CutPlots cuttrees "$2" 128)
# the lanes start after CutPlots has written the cuts/ files
cuts=$(submit Cuts cuts "$cut")
fom=$(submit Fom fom "$cut")
ce=$(submit CEblock ceblock "$cut")
sec=$(submit SecCE secce "$cut")
# tidy, side box and crop once every lane has ended, so the pdfs of the lanes that worked are still finished when one failed
DEPKIND=afterany submit Pdfs pdfs "$cuts:$fom:$ce:$sec"
