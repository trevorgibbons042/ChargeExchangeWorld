#!/bin/bash
# Resubmits the runs of a ScanThickness.sh scan that failed: a run folder that is still in build/root while none of our Full jobs has it (a run that finished was moved to DEST by Full.pbs).
# A run with all its sim output (128 output0_t*.root) is analysed again (Full.pbs with SKIPSIM=1, which also wipes its cuts/ and plots/); one with less is printed and left, it needs its sim again
# (qsub -v MAC=<macro>,RUN=<run>,DEST=<DEST> Full.pbs). A run is resubmitted at most twice (the count is in sweep_resubmitted.txt next to this script), so a run that fails for another reason stays visible.
# Usage: bash SweepFailed.sh          (only prints)
#        SUBMIT=1 bash SweepFailed.sh (submits)
cd "$(dirname "$0")" || exit 1
HERE=$PWD
DEST=${DEST:-/lcrc/group/eic/data/NPol/0.00000_Trevor}
NOUT=${NOUT:-128}
cd ../build/root || exit 1
active=$(for i in $(qstat -u "$USER" | awk 'NR>5 && $4 ~ /^Full_/ {print $1}'); do qstat -f "$i" | tr -d '\n\t ' | grep -o 'RUN=[^,]*' | head -1 | cut -d= -f2; done)
touch "$HERE/sweep_resubmitted.txt"
for d in newPhysics_801_100m_*m/*/*block; do
    [ -d "$d" ] || continue
    echo "$active" | grep -qx "$d" && continue
    n=$(ls "$d"/output0_t*.root 2>/dev/null | wc -l)
    tries=$(grep -cx "$d" "$HERE/sweep_resubmitted.txt")
    if [ "$n" -lt "$NOUT" ]; then echo "SIM INCOMPLETE ($n output files), needs its sim again: $d"; continue; fi
    if [ "$tries" -ge 2 ]; then echo "GAVE UP after $tries resubmits: $d"; continue; fi
    e=$(echo "$d" | cut -d/ -f2); b=$(echo "$d" | cut -d/ -f3 | sed 's/block//')
    if [ -n "$SUBMIT" ]; then
        id=$(cd "$HERE" && qsub -N "Full_${e}_${b}block" -m n -v RUN="$d",SKIPSIM=1,DEST="$DEST" Full.pbs)
        echo "$d" >> "$HERE/sweep_resubmitted.txt"
        echo "resubmitted $d as $id"
    else
        echo "would resubmit $d"
    fi
done
