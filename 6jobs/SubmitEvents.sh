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

# CutEvents, Recon (or Convert), Calc, CountCE and CutPlots run as one job on a whole node, one after the other (Plot.pbs)
if [ "$MODE" = cutdone ]; then step=analysisconvert; else step=analysis; fi
ana=$(submit Analysis $step "" 128)
# the plot jobs after it start at Cuts (CutPlots is in the analysis job); SubmitPlots.sh only clears plots/ when it starts at CutPlots, so do it here
rm -rf "../build/root/$RUN/plots"; rm -f ../build/root/$RUN/Plot3D_part*.log
bash SubmitPlots.sh "$RUN" "$ana" Cuts
