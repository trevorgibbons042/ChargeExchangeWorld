#!/bin/bash
# Submits the analysis chain of 5rootMacros/run.sh for one run folder as PBS jobs (Plot.pbs), then the plot chain (SubmitPlots.sh), each job waiting for the ones it needs:
#   CutEvents -> ReconAngleAll (or ConvertAngles) -> CalcAPandT, and CountCE after CutEvents, then CutPlots and the plot jobs after both.
# Usage: bash SubmitEvents.sh <folder in build/root> [cutdone], e.g. bash SubmitEvents.sh Physics_801 cutdone
# cutdone: cuts/cut_angles_all.root is already there from a finished ReconAngleAll.cpp (the old 5-tree file or the new one). ConvertAngles.cpp then puts it in the one-tree
# format instead of ReconAngleAll.cpp fitting every proton again (about an hour at 100M events); the old file stays as cut_angles_all_old5.root.
# Needs output0_t*.root and settings.csv in the run folder.
# A plot job name (CutPlots, Cuts, Fom, CEblock, SecCE, Pdfs) instead of cutdone submits no analysis job and starts the plot chain at that job (SubmitPlots.sh), e.g. bash SubmitEvents.sh Physics_801 Cuts.
RUN=$1
[ -n "$RUN" ] || { echo "usage: bash SubmitEvents.sh <run folder> [cutdone | CutPlots|Cuts|Fom|CEblock|SecCE|Pdfs]"; exit 1; }
MODE=$2
cd "$(dirname "$0")" || exit 1
case "$MODE" in
    ""|cutdone) ;;
    CutPlots|Cuts|Fom|CEblock|SecCE|Pdfs) exec bash SubmitPlots.sh "$RUN" "" "$MODE";;
    *) echo "unknown option $MODE (cutdone, or a plot job: CutPlots Cuts Fom CEblock SecCE Pdfs)"; exit 1;;
esac
source ./submit_lib.sh
mkdir -p "../build/root/$RUN/cuts"
rsync -a --exclude __pycache__ ../5rootMacros/ "../build/root/$RUN/" || exit 1 # the jobs run the run folder's copy of the macros

# the analysis jobs take a whole node, as in run.sh; ConvertAngles.cpp runs one block of rows per thread too
cutev=$(submit CutEvents cutevents "" 128)
if [ "$MODE" = cutdone ]; then
    angles=$(submit Convert convert "$cutev" 128)
else
    angles=$(submit Recon recon "$cutev" 128)
fi
calc=$(submit Calc calc "$angles" 128)
count=$(submit CountCE countce "$cutev" 128)
# CutPlots and every plot job after CalcAPandT and CountCE
bash SubmitPlots.sh "$RUN" "$calc:$count"
