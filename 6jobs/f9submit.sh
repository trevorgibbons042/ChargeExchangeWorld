#!/bin/bash -l
cd /lcrc/project/eic/tgibbons
qdel 7869237
B=/lcrc/project/eic/tgibbons/G4ChargeExchange
for d in $B/build/root/newPhysics_801_100m_*m/50GeV/*block; do
  [ -s $d/plots/6_ReconNew_1/6_megaCut/scan_stats.csv ] || continue
  [ $d/plots/6_ReconNew_1/1_noCuts/1_KEEnd.pdf -nt $B/5rootMacros/F9_PlotReconNew.cpp ] && continue
  t=$(echo $d | sed 's#.*_\([0-9.]*\)m/50GeV/.*#\1#'); b=$(basename $d | sed 's/block//')
  qsub -N F9Ke_${t}_$b -o /lcrc/project/eic/tgibbons/f9ke_${t}_$b.log -v D=$d f9one.pbs
done
