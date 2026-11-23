#!/bin/bash
# The runs the first ScanGrid.sh left broken (the eic project quota was full), grouped in batches, smallest total carbon first. A batch starts when the previous one has ended (afterany), so only a few runs hold disk at a time.
# Each entry is <thickness>:<block count>. PREV (job ids a:b:c) makes the first batch wait for those jobs. Usage: bash ScanGridBatches.sh (DRY=1 prints the batches only)
cd "$(dirname "$0")" || exit 1
export PATH=$PATH:/opt/pbs/bin
export FULL=1 FULLPBS=FullKeep.pbs ENERGIES="50:1.8" ADDTRK="1:10 2:4 3:2 4:0 5:0 6:0" EFF=0.95
BATCHES=("0.7:1 0.3:4" "0.7:2 0.3:5 0.5:3" "0.3:6 0.5:4" "0.7:3 0.5:5" "0.7:4 0.5:6 0.7:5" "0.7:6")
if [ -n "$ONLY_AFTER" ]; then BATCHES=("${BATCHES[@]:$ONLY_AFTER}"); fi # skip the first batches already submitted
prev="${PREV:-}"
for b in "${BATCHES[@]}"; do
    ids=""
    for run in $b; do
        t=${run%%:*}; n=${run##*:}
        if [ -n "$DRY" ]; then echo "batch [$b] run t=$t n=$n after=${prev:-none}"; continue; fi
        out=$(QSUBOPTS="${prev:+-W depend=afterany:$prev}" COUNTS="$n" bash ScanBlocksEach.sh run_100m.mac "$t")
        echo "$out" | grep -v '^$'
        id=$(echo "$out" | grep -o '[0-9]*\.imgt1' | head -1)
        ids="$ids:$id"
    done
    prev=${ids#:}
done
