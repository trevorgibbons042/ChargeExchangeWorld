#!/bin/bash
# Submits the plotting chain for one run folder as PBS jobs (Plot.pbs), the lanes of 6jobs/plot.sh in parallel, instead of running plot.sh in one go.
# Usage: bash SubmitPlots.sh <folder in build/root> [job ids the first job waits for, colon separated] [first job] [plots/ folder to redo], e.g. bash SubmitPlots.sh Physics_xs1e10
# The Pdfs job is only submitted with PDFS=1 in the environment or with Pdfs as the first job (the examples below that mention it need one of them).
# first job: Plot3D (default), CEblock, SecCE or Pdfs. The jobs before it are not submitted, so e.g. "CEblock" runs CEblock, SecCE and Pdfs: bash SubmitPlots.sh Physics_801 "" CEblock
# CutPlots is still accepted as the first job and means Plot3D (the CutPlots zdc job wrote nothing and is gone). With no first job, or CutPlots, plots/ is wiped first.
# plots/ folder to redo: one folder of plots/ (a leading plots/ is dropped), e.g. bash SubmitPlots.sh Physics_801 "" Plot3D 4_CEblock. It is deleted, then only the job that writes it runs, then Pdfs.
#   1_position, 2_energy  the Plot3D job      4_CEblock  the CEblock lane      7_secondaryCE  the SecCE lane
# Needs cut_*.root and settings.csv already in the run folder (CutPlots.cpp (events), ReconAngleAll.cpp, CalcAPandT.cpp, CutPlots.cpp (count) done), or the job ids of the jobs that write them (SubmitEvents.sh).
RUN=$1
ORDER="Plot3D CEblock SecCE Pdfs"
# the Pdfs job (tidy_plots.py, crop_pdfs.py) is not submitted unless it is asked for: PDFS=1 in the environment, or Pdfs as the first job
[ "${3:-}" = Pdfs ] && PDFS=1
WIPE=; case "${3:-}" in ""|CutPlots) WIPE=1;; esac # a run from the start rebuilds every plot
FROM=${3:-Plot3D}; [ "$FROM" = CutPlots ] && FROM=Plot3D # the old name of the first job
TARGET=$4
[ -n "$RUN" ] || { echo "usage: bash SubmitPlots.sh <run folder> [job ids to wait for] [first job: $ORDER] [plots/ folder to redo]"; exit 1; }
case " $ORDER " in *" $FROM "*) ;; *) echo "unknown first job $FROM (one of: $ORDER)"; exit 1;; esac
cd "$(dirname "$0")" || exit 1
source ./submit_lib.sh

LANE=
if [ -n "$TARGET" ]; then
    TARGET=${TARGET#plots/}; TARGET=${TARGET%/}
    case "$TARGET" in ""|.*|*..*|/*) echo "bad folder: $4"; exit 1;; esac
    top=${TARGET%%/*}
    case $top in
        1_position|2_energy) LANE=plot3d ;;
        4_CEblock) LANE=ceblock ;;
        7_secondaryCE) LANE=secce ;;
        *) echo "unknown plots folder $top (one of: 1_position 2_energy 4_CEblock 7_secondaryCE)"; exit 1 ;;
    esac
    [ -d "../build/root/$RUN/plots/$TARGET" ] || { echo "no folder ../build/root/$RUN/plots/$TARGET"; exit 1; }
    rm -rf "../build/root/$RUN/plots/$TARGET" || exit 1 # the old plots go, so nothing stale stays beside the new ones
fi

rsync -a --exclude __pycache__ --exclude '*_cpp.*' --exclude .DS_Store --exclude CLAUDE.md --exclude debug ../5rootMacros/ "../build/root/$RUN/" || exit 1 # the jobs run the run folder's copy of the macros
rsync -a ../6jobs/run_all.sh ../6jobs/run_norecon.sh ../6jobs/plot.sh ../6jobs/run_lib.sh ../6jobs/preflight.sh ../6jobs/to_png.sh "../build/root/$RUN/" || exit 1 # the scripts of the chain live in 6jobs/
# a run from the start (no folder, no first job or CutPlots) rebuilds every plot: the old plots/ and the Plot3D logs go first, so nothing of an older macro version stays beside the new files. cuts/ stays.
if [ -z "$TARGET" ] && [ -n "$WIPE" ]; then
    rm -rf "../build/root/$RUN/plots"; rm -f ../build/root/$RUN/Plot3D_sub*.log
fi

if [ -n "$TARGET" ]; then # one folder: its job only (the cuts/ files are already there), then the python steps
    ids=$(submit Redo $LANE "$2")
    if [ -n "$PDFS" ]; then DEPKIND=afterany submit Pdfs pdfs "$ids"; fi
    exit
fi

# on <job>: true when the job is FROM or comes after it in ORDER
on() { [ "$1" = "$FROM" ] || case " $ORDER " in *" $FROM "*"$1 "*) return 0;; *) return 1;; esac; }

# the lanes wait for the ids in $2 (the analysis job of SubmitEvents.sh), if any
w=$2
on Plot3D && p3=$(submit Plot3D plot3d "$w" 128) # Plot3D runs three processes side by side (Plot.pbs)
on CEblock && ce=$(submit CEblock ceblock "$w")
on SecCE && sec=$(submit SecCE secce "$w")
# tidy and crop once every lane has ended, so the pdfs of the lanes that worked are still finished when one failed
ids=$(echo $p3 $ce $sec | tr ' ' ':')
if [ -n "$PDFS" ]; then DEPKIND=afterany submit Pdfs pdfs "${ids:-$w}"; fi
