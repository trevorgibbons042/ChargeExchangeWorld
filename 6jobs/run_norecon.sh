# run_all.sh without ReconAngleAll.cpp: reuses the cuts/cut_angles_all.root already in the run folder (run it once with run_all.sh first)

source ./run_lib.sh # run <macro>: times it, ACLiC or interpreted, retries in the other mode on Failed to materialize symbols

if [ ! -f cuts/cut_angles_all.root ]; then
    echo "cuts/cut_angles_all.root missing: run ReconAngleAll.cpp (run_all.sh) first"
    exit 1
fi

NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}

sh preflight.sh || exit 1 # ROOT features, compile test, writable cuts/ plots/, disk; PREFLIGHT=0 skips
export PREFLIGHT_DONE=1 # plot.sh below does not repeat it

export ANA_NTHREADS=$NCPU
run 'CutPlots.cpp("events")' # rewrites the cut files cuts/cut_*.root (CutFiles.h)

# CutPlots count needs only the CutPlots events output, CalcAPandT only the angles: two chains side by side, half of the cores each
export ANA_NTHREADS=$(( NCPU / 2 > 0 ? NCPU / 2 : 1 ))
run 'CutPlots.cpp("count")' &
{
    run 'CalcAPandT.cpp' # the detector KE trees of cut_APandT.root
} &
wait
unset ANA_NTHREADS
zsh plot.sh
