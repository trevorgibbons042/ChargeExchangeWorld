#!/bin/bash
# Pre-submits the first 45 runs (groups 0 to 8) of resim_order.txt as held Full jobs, and one 1 cpu waiting job (ScanNext_g<G>) per group that releases (qrls) the 5 held jobs of its group when the jobs
# 3 groups before have ended (afterany), so 15 jobs run at a time. Groups 0, 1 and 2 wait for the 0.1 m runs. Groups 9, 10 and 11 are started by waiting jobs (ResimGroup.sh), which then chain on by themselves.
cd "$(dirname "$0")" || exit 1
HERE=$PWD
LIST=$HERE/resim_order.txt
export DEST=${DEST:-/lcrc/group/eic/data/NPol/0.00000_Trevor}
export ADDTRK=${ADDTRK:-"1:10 2:4 3:2 4:1 5:0 6:0"}
export FULL=1 QSUBOPTS="-m n -h"
dep() { echo "$@" | sed 's/[0-9]\+/&.imgt1/g' | tr ' ' ':'; }
declare -A G
for g in 0 1 2 3 4 5 6 7 8; do
    for ((i = g * 5 + 1; i <= g * 5 + 5; i++)); do
        read -r t e n < <(sed -n "${i}p" "$LIST")
        case $e in 25) sig=0.92 ;; 50) sig=1.8 ;; esac
        id=$(ENERGIES="$e:$sig" COUNTS="$n" bash "$HERE/ScanBlocks.sh" run_100m.mac "$t" | grep -E '^[0-9]+\.' | head -1)
        G[$g]="${G[$g]:+${G[$g]} }$id"
    done
    echo "group $g held: ${G[$g]}"
done
wait_for() { # $1 group, $2 ids (space separated) it waits for, a ScanNext job that releases the group
    local ids; ids=$(echo $2 | tr ' ' ':')
    printf '#!/bin/bash -l\nqrls %s\n' "${G[$1]}" | qsub -N "ScanNext_g$1" -q backfill -A eic -l select=1:ncpus=1 -l walltime=01:00:00 -j oe -m n -o "$HERE/ResimGroup_$1.log" -W "depend=afterany:$ids"
}
wait_for 0 "7877724.imgt1 7877725.imgt1 7877726.imgt1 7877727.imgt1 7877728.imgt1 7877729.imgt1"
wait_for 1 "7877673.imgt1 7877674.imgt1 7877675.imgt1 7877676.imgt1 7877723.imgt1 7877730.imgt1"
wait_for 2 "7877724.imgt1 7877725.imgt1"
for g in 3 4 5 6 7 8; do wait_for $g "${G[$((g - 3))]}"; done
for g in 9 10 11; do
    ids=$(echo ${G[$((g - 3))]} | tr ' ' ':')
    printf '#!/bin/bash -l\nbash %s/ResimGroup.sh %s\n' "$HERE" "$g" | qsub -N "ScanNext_g$g" -q backfill -A eic -l select=1:ncpus=1 -l walltime=01:00:00 -j oe -m n -o "$HERE/ResimGroup_$g.log" -W "depend=afterany:$ids"
done
