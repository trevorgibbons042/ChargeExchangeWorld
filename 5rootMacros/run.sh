#source ~/.zshrc

# runs one macro and prints how long it took, compiled with ACLiC (+O) unless ACLIC= (empty), see plot.sh
run(){
    local start=$SECONDS
    root -l -b -q "${1/.cpp/.cpp${ACLIC-+O}}"
    echo ">>> $1 took $((SECONDS - start)) s"
}

NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}

sh preflight.sh || exit 1 # ROOT features, compile test, writable cuts/ plots/, disk; PREFLIGHT=0 skips
export PREFLIGHT_DONE=1 # plot.sh below does not repeat it

mkdir -p cuts # every cut_*.root and the KE cut lists
rm -rf cuts/cache # the cut caches (CutCache.h) are rebuilt from the new cut_*.root files
export ANA_NTHREADS=$NCPU
run 'CutEvents.cpp' # writes cuts/cut_*.root and registers them as cut cache seeds

# CountCE needs only the CutEvents output, ReconAngleAll and CalcAPandT only theirs: two chains side by side, half of the cores each
export ANA_NTHREADS=$(( NCPU / 2 > 0 ? NCPU / 2 : 1 ))
run 'CountCE.cpp' &
{
    run 'ReconAngleAll.cpp'
    run 'CalcAPandT.cpp'
    run 'CalcAPandT.cpp("fDEnergy")'
} &
wait
unset ANA_NTHREADS
zsh plot.sh
