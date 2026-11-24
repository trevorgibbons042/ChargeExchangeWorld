#!/bin/bash
# Resim of the thickness scan in groups of 5 runs, three groups at a time (15 jobs). resim_order.txt has one line per run (thickness energy blocks, thin first, 25 GeV before 50 GeV, few blocks first).
# Usage: bash ResimGroup.sh <group number>: submits the 5 runs of that group (lines 5G+1 to 5G+5 of resim_order.txt, Full.pbs with DEST), then a 1 cpu waiting job (ScanNext_g<G+3>)
# that starts group G+3 when all 5 jobs of this group have ended (afterany). Start groups 0, 1 and 2 by hand (or from waiting jobs) to run three chains side by side.
cd "$(dirname "$0")" || exit 1
HERE=$PWD
G=$1
[ -n "$G" ] || { echo "usage: bash ResimGroup.sh <group number>"; exit 1; }
SIZE=5; LANES=3
LIST=$HERE/resim_order.txt
total=$(wc -l < "$LIST")
export DEST=${DEST:-/lcrc/group/eic/data/NPol/0.00000_Trevor}
export ADDTRK=${ADDTRK:-"1:10 2:4 3:2 4:1 5:0 6:0"}
export FULL=1 QSUBOPTS="-m n"
ids=""
for ((i = G * SIZE + 1; i <= G * SIZE + SIZE && i <= total; i++)); do
    read -r t e n < <(sed -n "${i}p" "$LIST")
    case $e in 25) sig=0.92 ;; 50) sig=1.8 ;; *) echo "bad energy $e"; exit 1 ;; esac
    id=$(ENERGIES="$e:$sig" COUNTS="$n" bash "$HERE/ScanBlocks.sh" run_100m.mac "$t" | grep -E '^[0-9]+\.' | head -1)
    echo "$(date) group $G: $t m $e GeV $n blocks as $id (run $i of $total)"
    ids="${ids:+$ids:}$id"
done
[ -n "$ids" ] || { echo "group $G is empty"; exit 0; }
next=$((G + LANES))
if [ $((next * SIZE)) -lt "$total" ]; then
    printf '#!/bin/bash -l\nbash %s/ResimGroup.sh %s\n' "$HERE" "$next" \
        | qsub -N "ScanNext_g$next" -q backfill -A eic -l select=1:ncpus=1 -l walltime=01:00:00 -j oe -m n -o "$HERE/ResimGroup_$next.log" -W "depend=afterany:$ids"
fi
