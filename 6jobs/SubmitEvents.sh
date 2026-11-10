#!/bin/bash
# Submits the analysis chain of 5rootMacros/run.sh for one run folder as PBS jobs (Plot.pbs), then the plot chain (SubmitPlots.sh), each job waiting for the ones it needs:
#   CutEvents -> ReconAngleAll -> CalcAPandT, and CountCE after CutEvents, then CutPlots and the plot jobs after both.
# Usage: bash SubmitEvents.sh <folder in build/root>, e.g. bash SubmitEvents.sh Physics_801
# Needs output0_t*.root and settings.csv in the run folder.
# A plot job name (CutPlots, Cuts, Fom, CEblock, SecCE, Smear, Pdfs) as the second argument submits no analysis job and starts the plot chain at that job (SubmitPlots.sh), e.g. bash SubmitEvents.sh Physics_801 Cuts.
RUN=$1
[ -n "$RUN" ] || { echo "usage: bash SubmitEvents.sh <run folder> [CutPlots|Cuts|Fom|CEblock|SecCE|Smear|Pdfs]"; exit 1; }
MODE=$2
cd "$(dirname "$0")" || exit 1
case "$MODE" in
    "") ;;
    CutPlots|Cuts|Fom|CEblock|SecCE|Smear|Pdfs) exec bash SubmitPlots.sh "$RUN" "" "$MODE";;
    *) echo "unknown option $MODE (a plot job: CutPlots Cuts Fom CEblock SecCE Smear Pdfs)"; exit 1;;
esac
source ./submit_lib.sh
mkdir -p "../build/root/$RUN/cuts"
rsync -a --exclude __pycache__ ../5rootMacros/ "../build/root/$RUN/" || exit 1 # the jobs run the run folder's copy of the macros

# CutEvents, Recon, Calc, CountCE and CutPlots run as one job on a whole node, one after the other (Plot.pbs)
# ANA_AFTER (optional) is a job id this Analysis job waits for, so several runs read the disk one after another (ScanBlockAnalysis.sh); afterany so one failed run does not stall the rest. ANA_ID_FILE (optional) gets this job's id.
ana=$(DEPKIND=afterany submit Analysis analysis "$ANA_AFTER" 128)
[ -n "$ANA_ID_FILE" ] && echo "$ana" > "$ANA_ID_FILE"
# the plot jobs after it start at Cuts (CutPlots is in the analysis job); SubmitPlots.sh only clears plots/ when it starts at CutPlots, so do it here
rm -rf "../build/root/$RUN/plots"; rm -f ../build/root/$RUN/Plot3D_part*.log
bash SubmitPlots.sh "$RUN" "$ana" Cuts
