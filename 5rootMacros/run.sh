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
export ANA_NTHREADS=$NCPU
run 'CutEvents.cpp' # writes the cut files cuts/cut_*.root (CutOpen.h)

# One after the other, all the cores each: ReconAngleAll alone on 8 threads took 39 s and CountCE (it reads only cut_CE.root) 4 s, against 61 s for the two side by side on half of the cores each
export ANA_NTHREADS=$NCPU
run 'ReconAngleAll.cpp'
run 'CountCE.cpp'
run 'CalcAPandT.cpp' # true KE, then detector KE
unset ANA_NTHREADS
zsh plot.sh
