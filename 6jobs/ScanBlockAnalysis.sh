#!/bin/bash
# Analysis for a whole block scan: runs SubmitEvents.sh on every run folder of a ScanBlocks.sh scan.
# Usage: bash ScanBlockAnalysis.sh <scan folder in build/root> [plot job name], e.g. bash ScanBlockAnalysis.sh Physics_801_0.4m
# The scan folder holds <E>GeV/<n>block run folders (ScanBlocks.sh names them so); a run folder is any folder two levels down with a settings.csv.
# The second argument is passed on to every SubmitEvents.sh call.
SCAN=$1
MODE=$2
[ -n "$SCAN" ] || { echo "usage: bash ScanBlockAnalysis.sh <scan folder in build/root> [CutPlots|Cuts|Fom|CEblock|SecCE|Pdfs]"; exit 1; }
cd "$(dirname "$0")" || exit 1
[ -d "../build/root/$SCAN" ] || { echo "no ../build/root/$SCAN"; exit 1; }
found=0
prev="" # the Analysis jobs run one after another (each waits for the previous one), so the runs do not all read the disk at once
idfile=$(mktemp)
for dir in ../build/root/"$SCAN"/*/*/; do
    [ -f "${dir}settings.csv" ] || continue
    run=${dir#../build/root/}
    run=${run%/}
    echo "$run"
    ANA_AFTER=$prev ANA_ID_FILE=$idfile bash SubmitEvents.sh "$run" $MODE
    [ -s "$idfile" ] && prev=$(cat "$idfile") # empty for a plot-job mode: no Analysis job to chain
    found=$((found + 1))
done
rm -f "$idfile"
[ "$found" -gt 0 ] || { echo "no run folder with settings.csv under ../build/root/$SCAN/*/*"; exit 1; }
echo "$found run folders submitted"
