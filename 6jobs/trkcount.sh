#!/bin/bash -l
#PBS -q backfill
#PBS -A eic
#PBS -N TrkCount
#PBS -l select=1:ncpus=32
#PBS -l walltime=01:00:00
#PBS -j oe

# tracks in cuts/cut_angles_all.root (AnglesRaw, one row per recon track) of the 2-block tracker-test runs, before any cut
module load root
export LD_LIBRARY_PATH=/lcrc/project/eic/tgibbons/xxhash-install/lib:/lcrc/project/eic/tgibbons/gif-lib:$LD_LIBRARY_PATH
cat > /lcrc/project/eic/tgibbons/trkcount.C <<'EOF'
void trkcount(const char* d){
  ROOT::EnableImplicitMT(32);
  ROOT::RDataFrame df("AnglesRaw", TString(d)+"/cuts/cut_angles_all.root");
  auto all   = df.Count();
  auto first = df.Filter("trackNo==0").Count();
  auto more  = df.Filter("trackNo>0").Count();
  auto nonCEfirst = df.Filter("trackNo==0 && evClass==2").Count();
  auto nonCEall   = df.Filter("evClass==2").Count();
  auto ceFirst    = df.Filter("trackNo==0 && evClass==0").Count();
  printf("%s rows=%llu events_with_track=%llu extra_tracks=%llu (%.1f%% of rows) nonCE_rows=%.1f%% nonCE_among_first=%.1f%% CE_first=%llu\n", d,
    *all, *first, *more, 100.*(*more)/(*all), 100.*(*nonCEall)/(*all), 100.*(*nonCEfirst)/(*first), *ceFirst);
}
EOF
cd /lcrc/project/eic/tgibbons/G4ChargeExchange/build/root
for a in 0 1 2 3; do
  root -l -b -q "/lcrc/project/eic/tgibbons/trkcount.C(\"newPhysics_801_100m_trk${a}_0.3m/25GeV/2block\")" 2>&1 | grep -a -E "rows=|rror"
done > /lcrc/project/eic/tgibbons/trkcount.txt 2>&1
