#!/bin/bash -l
cd /lcrc/project/eic/tgibbons
qdel 7869184
sleep 3
B=/lcrc/project/eic/tgibbons/G4ChargeExchange
for d in $B/build/root/newPhysics_801_100m_*m/50GeV/*block; do
  [ -s $d/plots/6_ReconNew_1/6_megaCut/scan_stats.csv ] && continue
  [ -s $d/cuts/cut_angles_all.root ] || continue
  rm -rf $d/plots/6_ReconNew_1
  t=$(echo $d | sed 's#.*_\([0-9.]*\)m/50GeV/.*#\1#'); b=$(basename $d | sed 's/block//')
  qsub -N F9Fx_${t}_$b -o /lcrc/project/eic/tgibbons/f9fx_${t}_$b.log -v D=$d f9one.pbs
done
