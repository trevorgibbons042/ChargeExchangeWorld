#!/bin/bash
# Rewrites plots/6_ReconNew_1 (and its scan_stats.csv files) for run folders that already have the recon, one job per run folder.
# Usage: bash SubmitPlotNew.sh <run folder in build/root> [more run folders ...], e.g. bash SubmitPlotNew.sh newPhysics_801_250m_0.4m/25GeV/8block newPhysics_801_250m_0.4m/25GeV/9block
# A run folder with cuts/cut_angles_all.root runs only F9_PlotReconNew.cpp (step plotnew). One without it runs the recon first (step reconnew).
# Only the csv files: SCAN_ONLY=1 makes F9_PlotReconNew.cpp write cuts.txt and scan_stats.csv of every folder and no pdf. No Pdfs job either.
export SCAN_ONLY=1
[ $# -gt 0 ] || { echo "usage: bash SubmitPlotNew.sh <run folder in build/root> [more run folders ...]"; exit 1; }
cd "$(dirname "$0")" || exit 1
source ./submit_lib.sh
for RUN in "$@"; do
    RUN=${RUN%/}
    [ -f "../build/root/$RUN/settings.csv" ] || { echo "no ../build/root/$RUN/settings.csv, skipped"; continue; }
    rsync -a --exclude __pycache__ --exclude '*_cpp.*' --exclude .DS_Store --exclude CLAUDE.md --exclude debug ../5rootMacros/ "../build/root/$RUN/" || exit 1 # the job runs the run folder's copy of the macros
    rsync -a ../6jobs/run_all.sh ../6jobs/run_norecon.sh ../6jobs/plot.sh ../6jobs/run_lib.sh ../6jobs/preflight.sh ../6jobs/to_png.sh "../build/root/$RUN/" || exit 1 # the scripts of the chain live in 6jobs/
    if [ -f "../build/root/$RUN/cuts/cut_angles_all.root" ]; then step=plotnew; else step=reconnew; fi
    echo "$RUN $step $(submit PlotNew $step "" 128)"
done
