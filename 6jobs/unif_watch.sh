#!/bin/bash
# qdel each unif_all job once it starts a run that unifx or unify already covers
cd /lcrc/project/eic/tgibbons/G4ChargeExchange/6jobs
declare -A JOB=([0]=7878637 [1]=7878638 [2]=7878639)
declare -A MAX=([0]=25 [1]=24 [2]=26)
left="0 1 2"
while [ -n "$left" ]; do
  new=""
  for i in $left; do
    n=$(grep -c "^===" unif_all_$i.log)
    if [ "$n" -gt "${MAX[$i]}" ] || grep -q ALLDONE unif_all_$i.log; then
      /opt/pbs/bin/qdel ${JOB[$i]}.imgt1 2>&1 || qdel ${JOB[$i]}.imgt1
      echo "$(date +%T) qdel unif_all$i at $n"
    else
      new="$new $i"
    fi
  done
  left=$(echo $new)
  sleep 60
done
echo WATCHDONE
