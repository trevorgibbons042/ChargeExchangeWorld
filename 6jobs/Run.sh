#!/bin/bash
# Submits Run.pbs, the sim, with one macro from build/.
# Usage: bash Run.sh <macro in build/>, e.g. bash Run.sh run_0.mac
MAC=$1
[ -n "$MAC" ] || { echo "usage: bash Run.sh <macro in build/>"; exit 1; }
cd "$(dirname "$0")" || exit 1
[ -f "../build/$MAC" ] || { echo "no ../build/$MAC"; exit 1; }
qsub -N "Run_${MAC%.mac}" -v MAC="$MAC" Run.pbs
