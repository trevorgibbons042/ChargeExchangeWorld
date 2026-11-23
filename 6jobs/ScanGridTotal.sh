#!/bin/bash
# The 50 GeV grid by TOTAL carbon thickness: 1 to 6 CE blocks, total 0.1, 0.3, 0.5, 0.7, 0.9, 1.1 m split evenly over the blocks (n blocks of T/n each), 95% tracker efficiency, 100M neutrons.
# nAdditionalTrackers per gap by block count: 1:10 2:4 3:2 4:0 5:0 6:0 (trackers in all 12 12 12 8 10 12). One FullKeep.pbs job per run, runs whose folder already exists are skipped.
# Sliding window: run i waits (afterany) for run i-WINDOW, so at most WINDOW runs hold disk at once (the eic quota is shared). Smallest total, then fewest blocks, first.
# Usage: bash ScanGridTotal.sh (DRY=1 prints the list, WINDOW=6 by default)
cd "$(dirname "$0")" || exit 1
export PATH=$PATH:/opt/pbs/bin
export FULL=1 FULLPBS=FullKeep.pbs ENERGIES="50:1.8" ADDTRK="1:10 2:4 3:2 4:0 5:0 6:0" EFF=0.95
WINDOW=${WINDOW:-6}
ids=($PREVIDS) # job ids already running count toward the window (PREVIDS="id1 id2 ...")
i=${#ids[@]}
for T in 0.1 0.3 0.5 0.7; do for n in 1 2 3 4 5 6; do
    [ -d "../build/root/newPhysics_801_100m_${T}m/50GeV/${n}block" ] && { echo "skip T=$T n=$n (folder exists)"; continue; }
    dep=""; [ $i -ge $WINDOW ] && dep="-W depend=afterany:${ids[$((i-WINDOW))]}"
    if [ -n "$DRY" ]; then echo "run $i T=$T n=$n ${dep:+after run $((i-WINDOW))}"; ids+=("x"); i=$((i+1)); continue; fi
    out=$(QSUBOPTS="$dep" COUNTS="$n" bash ScanBlocks.sh run_100m.mac "$T")
    echo "$out" | grep -v '^$'
    ids+=("$(echo "$out" | grep -o '[0-9]*\.imgt1' | head -1)"); i=$((i+1))
done; done
