# run.sh without ReconAngleAll.cpp: reuses the cuts/cut_angles_all.root already in the run folder (run it once with run.sh first)

# runs one macro and prints how long it took, compiled with ACLiC (+O) unless ACLIC= (empty), see plot.sh
run(){
    local start=$SECONDS
    root -l -b -q "${1/.cpp/.cpp${ACLIC-+O}}"
    echo ">>> $1 took $((SECONDS - start)) s"
}

if [ ! -f cuts/cut_angles_all.root ]; then
    echo "cuts/cut_angles_all.root missing: run ReconAngleAll.cpp (run.sh) first"
    exit 1
fi

NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}

sh preflight.sh || exit 1 # ROOT features, compile test, writable cuts/ plots/, disk; PREFLIGHT=0 skips
export PREFLIGHT_DONE=1 # plot.sh below does not repeat it

export ANA_NTHREADS=$NCPU
run 'CutEvents.cpp' # rewrites the cut files cuts/cut_*.root (CutOpen.h)

# CountCE needs only the CutEvents output, CalcAPandT only the angles: two chains side by side, half of the cores each
export ANA_NTHREADS=$(( NCPU / 2 > 0 ? NCPU / 2 : 1 ))
run 'CountCE.cpp' &
{
    run 'CalcAPandT.cpp'
    run 'CalcAPandT.cpp("fDEnergy")'
} &
wait
unset ANA_NTHREADS
zsh plot.sh
