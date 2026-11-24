#!/bin/bash
# Thickness scan in batches: one batch is one total thickness, 25 and 50 GeV x 1 to 6 blocks = 12 Full.pbs jobs (the sim, the whole analysis chain, then the output*_t*.root go and the run folder
# moves to DEST/<folder>_<t>m/<E>GeV/<n>block, Full.pbs with DEST). Each batch submits a 1 cpu job that waits for its 12 jobs (afterany: one failed run does not stop the scan) and starts the next one,
# so no more than 12 jobs and one waiting job exist at a time (LCRC allows about 70).
# Usage: bash ScanThickness.sh 0.1 0.2 0.3 ... 1.1   (the first thickness now, the others one after the other)
# Extra trackers per block count: 12 planes on 1 to 4 and 6 blocks, 10 on 5 (2 n planes of the stack + nAdditionalTrackers per gap, FitCore.h StackZ).
# Environment: TEMPLATE (macro in build/, default run_100m.mac), DEST (default /lcrc/group/eic/data/NPol/0.00000_Trevor), ENERGIES, COUNTS, ADDTRK override the defaults. PREV is set by the waiting job to the
# batch before: its 12 moved runs are checked first (plots/11_errorBudget/P_result.csv at DEST) and the missing ones are printed.
cd "$(dirname "$0")" || exit 1
SCAN_DIR=$PWD
T=$1; shift
[ -n "$T" ] || { echo "usage: bash ScanThickness.sh <thickness m> [next thicknesses ...]"; exit 1; }
export TEMPLATE=${TEMPLATE:-run_100m.mac}
export DEST=${DEST:-/lcrc/group/eic/data/NPol/0.00000_Trevor}
export COUNTS=${COUNTS:-1 2 3 4 5 6}
export ENERGIES=${ENERGIES:-25:0.92 50:1.8}
export ADDTRK=${ADDTRK:-1:10 2:4 3:2 4:1 5:0 6:0}
FOLDER=$(awk '/^\/output\/folder /{print $2}' "../build/$TEMPLATE")

if [ -n "$PREV" ]; then
    missing=0
    for es in $ENERGIES; do for n in $COUNTS; do
        f="$DEST/${FOLDER}_${PREV}m/${es%%:*}GeV/${n}block/plots/11_errorBudget/P_result.csv"
        [ -f "$f" ] || { echo "!!! MISSING $f"; missing=$((missing + 1)); }
    done; done
    echo ">>> batch ${PREV} m: $missing of $(echo $ENERGIES | wc -w)x$(echo $COUNTS | wc -w) runs missing at $DEST"
fi

echo ">>> batch $T m: $(date)"
ids=$(FULL=1 QSUBOPTS="-m n" bash ScanBlocks.sh "$TEMPLATE" "$T" | tee /dev/stderr | grep -E '^[0-9]+\.' | paste -sd: -)
[ -n "$ids" ] || { echo "no job was submitted for $T m"; exit 1; }
echo ">>> jobs $ids"

if [ $# -gt 0 ]; then
    # the waiting job: starts the next batch when every job of this one has ended
    # the waiting job does not see this shell's environment, so the settings of this scan (ENERGIES, COUNTS, ADDTRK, TEMPLATE, DEST) are written into its script: a lane of one energy stays a lane of that energy
    printf '#!/bin/bash -l\nPREV=%s ENERGIES="%s" COUNTS="%s" ADDTRK="%s" TEMPLATE=%s DEST=%s bash %s/ScanThickness.sh %s\n' "$T" "$ENERGIES" "$COUNTS" "$ADDTRK" "$TEMPLATE" "$DEST" "$SCAN_DIR" "$*" \
        | qsub -N "ScanNext_$1" -q backfill -A eic -l select=1:ncpus=1 -l walltime=01:00:00 -j oe -m n -o "$SCAN_DIR/ScanThickness_$1_${ENERGIES%%:*}.log" -W "depend=afterany:$ids"
fi
