#!/bin/bash
# Rebuilds plots/6_cuts of one run folder as PBS jobs, doing only what that folder reads (the other plots/ folders stay):
#   Analysis6 (CutEvents, ReconAngleAll, CalcAPandT true KE only, CutPlots zdc + samples; Plot.pbs analysis6) -> Plot3D parts of 6_cuts and QualityStudy side by side -> Pdfs (SubmitPlots.sh).
# Skipped against SubmitEvents.sh: CountCE, the detector-KE pass of CalcAPandT, the CutPlots trees part, and every plot lane outside 6_cuts.
# Usage: bash Submit6Cuts.sh <folder in build/root> [CutPlots | Plots], e.g. bash Submit6Cuts.sh Physics_801
# CutPlots: the cuts/ files of CutEvents, ReconAngleAll and CalcAPandT are already there, start at CutPlots. Plots: all cuts/ files are there, only the 6_cuts plot jobs.
# Needs output0_t*.root and settings.csv in the run folder.
RUN=$1
[ -n "$RUN" ] || { echo "usage: bash Submit6Cuts.sh <run folder> [CutPlots | Plots]"; exit 1; }
MODE=$2
case "$MODE" in ""|CutPlots|Plots) ;; *) echo "unknown option $MODE (CutPlots or Plots)"; exit 1;; esac
cd "$(dirname "$0")" || exit 1
source ./submit_lib.sh
mkdir -p "../build/root/$RUN/cuts" "../build/root/$RUN/plots/6_cuts" # SubmitPlots.sh redoes a folder that exists: it deletes plots/6_cuts and submits its jobs
rsync -a --exclude __pycache__ ../5rootMacros/ "../build/root/$RUN/" || exit 1 # the jobs run the run folder's copy of the macros
rm -f ../build/root/$RUN/Plot3D_part*.log

ana=
case "$MODE" in
    "") ana=$(submit Analysis6 analysis6 "" 128) ;;
    CutPlots) ana=$(submit CutPlots6 cutplots6 "" 128) ;;
esac
PDFS=1 bash SubmitPlots.sh "$RUN" "$ana" CutPlots 6_cuts
