#!/bin/bash
# wait until the 0.3 m 25 GeV sim outputs are deleted (1block and 3block were the last), then release the held 0.7 m jobs
export PATH=$PATH:/opt/pbs/bin
B=/lcrc/project/eic/tgibbons/G4ChargeExchange/build/root/newPhysics_801_100m_0.3m/25GeV
while true; do
  n=0
  for b in 1 2 3 4 5 6; do n=$((n + $(ls $B/${b}block/output*.root 2>/dev/null | wc -l))); done
  [ $n -eq 0 ] && break
  sleep 60
done
for j in 7869615 7869616 7869617 7869618 7869619 7869620; do
  qalter -W depend=afterany:7869184 $j.imgt1
done
echo "$(date +%m-%d_%H:%M) released 0.7 m"
