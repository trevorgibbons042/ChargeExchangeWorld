#!/bin/bash
# The 50 GeV grid: 1 to 6 CE blocks, every block 0.1, 0.3, 0.5, 0.7, 0.9 m thick (the total is n x thickness), 95% tracker efficiency, 100M neutrons each, one FullKeep.pbs job per run.
# nAdditionalTrackers per gap by block count: 1:10 2:4 3:2 4:0 5:0 6:0, so the trackers in all are 12 12 12 8 10 12 (2n + n x A). No TMVA, no USE_BDT.
# The jobs go in by total carbon (n x thickness, smallest first) so the fast runs start first and the slow thick ones queue behind them (the server runs only 15 of your jobs at once).
# Usage: bash ScanGrid.sh  (submits 30 jobs, DRY=1 only prints what it would run)
cd "$(dirname "$0")" || exit 1
export FULL=1 FULLPBS=FullKeep.pbs ENERGIES="50:1.8" ADDTRK="1:10 2:4 3:2 4:0 5:0 6:0" EFF=0.95
for t in 0.1 0.3 0.5 0.7; do for n in 1 2 3 4 5 6; do echo "$(awk -v t=$t -v n=$n 'BEGIN{printf "%.2f", t*n}') $n $t"; done; done | sort -n -k1,1 -k2,2 | while read -r tot n t; do
    if [ -n "$DRY" ]; then echo "total=$tot n=$n t=$t"; else COUNTS="$n" bash ScanBlocksEach.sh run_100m.mac "$t"; fi
done
