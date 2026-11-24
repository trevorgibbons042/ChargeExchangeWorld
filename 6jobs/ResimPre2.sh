#!/bin/bash
# Pre-submits groups 12 to 17 (runs 61 to 90 of resim_order.txt) as held Full jobs with one waiting ScanNext_g<G> job per group that releases it when the 5 jobs 3 groups before have ended.
# The waiting jobs ScanNext_g12 to g14 of ResimGroup.sh are replaced (their dependency lists are the job ids of groups 9 to 11); groups 18 to 20 are started by new ResimGroup.sh waiting jobs.
cd "$(dirname "$0")" || exit 1
HERE=$PWD
LIST=$HERE/resim_order.txt
export DEST=${DEST:-/lcrc/group/eic/data/NPol/0.00000_Trevor}
export ADDTRK=${ADDTRK:-"1:10 2:4 3:2 4:1 5:0 6:0"}
export FULL=1 QSUBOPTS="-m n -h"
declare -A G
# job ids of groups 9 to 11 = the dependency list of the waiting jobs g12 to g14
for g in 12 13 14; do
    id=$(qstat -u "$USER" | awk 'NR>5 && $4 ~ /^ScanNext/ {print $1}' | while read -r j; do qstat -f "$j" | grep -q "Job_Name = ScanNext_g$g\$" && echo "$j"; done | head -1)
    [ -n "$id" ] || { echo "no ScanNext_g$g"; exit 1; }
    dep=$(qstat -f "$id" | tr -d '\n\t ' | grep -o 'depend=afterany:[^,]*' | head -1 | sed 's/depend=afterany://; s/@imgt1//g')
    G[$((g - 3))]=$(echo "$dep" | tr ':' ' ')
    echo "group $((g - 3)) ids: ${G[$((g - 3))]}"
    qdel "$id"
done
for g in 12 13 14 15 16 17; do
    for ((i = g * 5 + 1; i <= g * 5 + 5; i++)); do
        read -r t e n < <(sed -n "${i}p" "$LIST")
        case $e in 25) sig=0.92 ;; 50) sig=1.8 ;; esac
        id=$(ENERGIES="$e:$sig" COUNTS="$n" bash "$HERE/ScanBlocks.sh" run_100m.mac "$t" | grep -E '^[0-9]+\.' | head -1)
        G[$g]="${G[$g]:+${G[$g]} }$id"
    done
    echo "group $g held: ${G[$g]}"
done
for g in 12 13 14 15 16 17; do
    ids=$(echo ${G[$((g - 3))]} | tr ' ' ':')
    printf '#!/bin/bash -l\nqrls %s\n' "${G[$g]}" | qsub -N "ScanNext_g$g" -q backfill -A eic -l select=1:ncpus=1 -l walltime=01:00:00 -j oe -m n -o "$HERE/ResimGroup_$g.log" -W "depend=afterany:$ids"
done
for g in 18 19 20; do
    ids=$(echo ${G[$((g - 3))]} | tr ' ' ':')
    printf '#!/bin/bash -l\nbash %s/ResimGroup.sh %s\n' "$HERE" "$g" | qsub -N "ScanNext_g$g" -q backfill -A eic -l select=1:ncpus=1 -l walltime=01:00:00 -j oe -m n -o "$HERE/ResimGroup_$g.log" -W "depend=afterany:$ids"
done
