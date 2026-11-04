#!/bin/bash
# Block scan with a fixed thickness per block: one sim run per block count (1, 3, 5, 7, 9) and beam energy (25, 50, 137.5 GeV), every block <thickness> m thick (total CE thickness = count x thickness).
# Writes one macro per count into build/ from a template macro and submits each with Run.sh.
# Usage: bash ScanBlocksEach.sh <template macro in build/> <thickness of each block in m> [dz in m, default 0.025], e.g. bash ScanBlocksEach.sh run_0.mac 0.4
# Each macro is the template with /detector/nBlocks, /detector/CarbonThickness (the block thickness) and /detector/dz set before /run/reinitializeGeometry,
# /randGun/Energy and /randGun/EnergySD set to the beam energy and its sigma, and /output/folder continued with _<thickness>mEach/<E>GeV/<n>block, so the runs sit in build/root/<folder>_<thickness>mEach/<E>GeV/<n>block (analyse them all with ScanBlockAnalysis.sh). The analysis reads output0_t*.root, so every count is its own run and run folder.
TEMPLATE=$1
THICK=$2
DZ=${3:-0.025}
[ -n "$TEMPLATE" ] && [ -n "$THICK" ] || { echo "usage: bash ScanBlocksEach.sh <template macro in build/> <thickness of each block m> [dz m]"; exit 1; }
cd "$(dirname "$0")" || exit 1
[ -f "../build/$TEMPLATE" ] || { echo "no ../build/$TEMPLATE"; exit 1; }
grep -q '^/run/reinitializeGeometry' "../build/$TEMPLATE" || { echo "$TEMPLATE has no /run/reinitializeGeometry line"; exit 1; }
grep -q '^/output/folder ' "../build/$TEMPLATE" || { echo "$TEMPLATE has no /output/folder line"; exit 1; }

COUNTS="1 3 5 7 9"
ENERGIES="25:0.92 50:1.8 137.5:5" # beam kinetic energy:Gaussian sigma, both in GeV
TRACKER_HALF=0.00025 # tracker planes are 500 um thick (2construction.cc)
ZDC_FRONT=7.7 # m, front face of the ZDC: the blocks sit in evenly spaced slots from z = 0 to here

# check every count before writing or submitting anything
for n in $COUNTS; do
    # dz is the gap between a block edge and its tracker edge. The blocks sit in slots 7.7/(n+1) apart; the tracker behind one block and the one in front of the next
    # (block + 2 gaps + 2 trackers wide) must not meet, which also covers the last block and the ZDC tracker.
    awk -v t="$THICK" -v dz="$DZ" -v n="$n" -v l="$ZDC_FRONT" -v th="$TRACKER_HALF" 'BEGIN{exit !(l/(n+1) <= t + 2*dz + 4*th)}' \
        && { echo "$n blocks: slots are $(awk -v n="$n" -v l="$ZDC_FRONT" 'BEGIN{printf "%.4g", l/(n+1)}') m apart and a block is $THICK m thick, so dz=$DZ m puts neighbouring trackers on top of each other, lower dz below $(awk -v t="$THICK" -v n="$n" -v l="$ZDC_FRONT" -v th="$TRACKER_HALF" 'BEGIN{printf "%.4g", (l/(n+1) - t - 4*th)/2}') m (or use thinner blocks)"; exit 1; }
done

BASE=${TEMPLATE%.mac}
for es in $ENERGIES; do
    e=${es%%:*}
    esd=${es##*:}
    for n in $COUNTS; do
        tag="${e}GeV_${n}block_${THICK}mEach"
        mac="${BASE}_${tag}.mac"
        awk -v n="$n" -v t="$THICK" -v dz="$DZ" -v e="$e" -v esd="$esd" -v suffix="_${THICK}mEach/${e}GeV/${n}block" '
            /^#?\/detector\/(nBlocks|CarbonThickness|dz) / { next }
            /^#?\/randGun\/(Energy|EnergySD) / { next }
            /^\/run\/reinitializeGeometry/ { print "/detector/nBlocks " n; print "/detector/CarbonThickness " t; print "/detector/dz " dz; print "/randGun/Energy " e; print "/randGun/EnergySD " esd }
            /^\/output\/folder / { print $0 suffix; next }
            { print }
        ' "../build/$TEMPLATE" > "../build/$mac"
        echo "$mac: $e GeV (sigma $esd), $n blocks x $THICK m (total $(awk -v t="$THICK" -v n="$n" 'BEGIN{printf "%.6g", t*n}') m)"
        bash Run.sh "$mac"
    done
done
