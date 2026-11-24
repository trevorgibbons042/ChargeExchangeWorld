#!/bin/bash
# Feeds the thickness scan resim in groups: when LOW or fewer jobs (Full_* and bgrej*) are running or queued, the next GROUP runs of resim_order.txt are submitted (LOW + GROUP = 14 jobs at most).
# resim_order.txt has one line per run: thickness energy blocks (thin first, 25 GeV before 50 GeV, few blocks first). resim_pos.txt holds how many lines are submitted.
# Runs as a PBS job (ResimFeeder.pbs) and submits itself again before its walltime ends.
# Usage: qsub -N ResimFeeder -q backfill -A eic -l select=1:ncpus=1 -l walltime=12:00:00 -j oe -m n -o $PWD/ResimFeeder.log ResimFeeder.pbs
cd "$(dirname "$0")" || exit 1
HERE=$PWD
LIST=$HERE/resim_order.txt
POS=$HERE/resim_pos.txt
GROUP=${GROUP:-5}
LOW=${LOW:-9}
export DEST=${DEST:-/lcrc/group/eic/data/NPol/0.00000_Trevor}
export ADDTRK=${ADDTRK:-"1:10 2:4 3:2 4:1 5:0 6:0"}
export FULL=1 QSUBOPTS="-m n"
total=$(wc -l < "$LIST")
end=$((SECONDS + 11 * 3600))
while [ $SECONDS -lt $end ]; do
    pos=$(cat "$POS" 2>/dev/null || echo 0)
    [ "$pos" -ge "$total" ] && { echo "$(date) all $total runs submitted"; exit 0; }
    cur=$(qstat -u "$USER" | awk 'NR>5 && ($10=="R" || $10=="Q") && ($4 ~ /^Full_/ || $4 ~ /^bgrej/) {c++} END{print c+0}')
    k=0
    [ "$cur" -le "$LOW" ] || k=$GROUP # more than LOW jobs left: no group yet
    while [ "$k" -lt "$GROUP" ] && [ "$pos" -lt "$total" ]; do
        read -r t e n < <(sed -n "$((pos + 1))p" "$LIST")
        case $e in 25) sig=0.92 ;; 50) sig=1.8 ;; *) echo "bad energy $e"; exit 1 ;; esac
        id=$(ENERGIES="$e:$sig" COUNTS="$n" bash "$HERE/ScanBlocks.sh" run_100m.mac "$t" | grep -E '^[0-9]+\.' | head -1)
        echo "$(date) submitted $t m $e GeV $n blocks as $id ($((pos + 1)) of $total, $cur jobs before the group)"
        pos=$((pos + 1)); k=$((k + 1))
        echo "$pos" > "$POS"
    done
    sleep 120
done
# walltime nearly over and runs left: start the next feeder
cd "$HERE" && qsub -N ResimFeeder -q backfill -A eic -l select=1:ncpus=1 -l walltime=12:00:00 -j oe -m n -o "$HERE/ResimFeeder.log" ResimFeeder.pbs
