#!/bin/bash -l
cd /lcrc/project/eic/tgibbons
B=/lcrc/project/eic/tgibbons/G4ChargeExchange
for d in $B/build/root/newPhysics_801_100m_*m/50GeV/*block; do
  [ -s $d/cuts/cut_angles_all.root ] || continue
  t=$(echo $d | sed 's#.*_\([0-9.]*\)m/50GeV/.*#\1#'); b=$(basename $d | sed 's/block//')
  qsub -N F9Mg_${t}_$b -o /lcrc/project/eic/tgibbons/f9mg_${t}_$b.log -v D=$d f9one.pbs
done
