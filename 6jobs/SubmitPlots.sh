#!/bin/bash
# Submits the plotting chain for one run folder as PBS jobs (Plot.pbs), the lanes of 5rootMacros/plot.sh in parallel, instead of running plot.sh in one go.
# Usage: bash SubmitPlots.sh <folder in build/root> [job ids CutPlots waits for, colon separated] [first job], e.g. bash SubmitPlots.sh Physics_xs1e10
# first job: CutPlots (default), Cuts, Fom, CEblock, SecCE or Pdfs. The jobs before it are not submitted (their files are already in cuts/), so e.g. "Cuts" runs Cuts, Fom, CEblock, SecCE and Pdfs: bash SubmitPlots.sh Physics_801 "" Cuts
# Needs cut_*.root and settings.csv already in the run folder (CutEvents.cpp, ReconAngleAll.cpp, CalcAPandT.cpp, CountCE.cpp done), or the job ids of the jobs that write them (SubmitEvents.sh).
RUN=$1
ORDER="CutPlots Cuts Fom CEblock SecCE Pdfs"
FROM=${3:-CutPlots}
[ -n "$RUN" ] || { echo "usage: bash SubmitPlots.sh <run folder> [job ids to wait for] [first job: $ORDER]"; exit 1; }
case " $ORDER " in *" $FROM "*) ;; *) echo "unknown first job $FROM (one of: $ORDER)"; exit 1;; esac
cd "$(dirname "$0")" || exit 1
source ./submit_lib.sh
rsync -a --exclude __pycache__ ../5rootMacros/ "../build/root/$RUN/" || exit 1 # the jobs run the run folder's copy of the macros

# on <job>: true when the job is FROM or comes after it in ORDER
on() { [ "$1" = "$FROM" ] || case " $ORDER " in *" $FROM "*"$1 "*) return 0;; *) return 1;; esac; }

# CutPlots (KE cut lists, ZDC proton cut files, sample cut files) is on the critical path, ahead of every lane: all the cores
on CutPlots && cut=$(submit CutPlots cuttrees "$2" 128)
# the lanes start after CutPlots has written the cuts/ files; without CutPlots they wait for the ids in $2, if any
w=${cut:-$2}
on Cuts && cuts=$(submit Cuts cuts "$w" 128) # Plot3D runs nine processes side by side (Plot.pbs): a whole node
on Fom && fom=$(submit Fom fom "$w")
on CEblock && ce=$(submit CEblock ceblock "$w")
on SecCE && sec=$(submit SecCE secce "$w")
# tidy, side box and crop once every lane has ended, so the pdfs of the lanes that worked are still finished when one failed
ids=$(echo $cuts $fom $ce $sec | tr ' ' ':')
DEPKIND=afterany submit Pdfs pdfs "${ids:-$w}"
