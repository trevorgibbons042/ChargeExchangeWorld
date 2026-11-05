# Sourced by run.sh, run_norecon.sh, plot.sh and plot_cut.sh: the run function.
# runs one macro and prints how long it took. Compiled with ACLiC (+O: optimized, rebuilt only when the macro or a header changed);
# ACLIC= (empty) runs it interpreted instead, for a ROOT whose compiled macros do not load (the Mac one: runStaticInitializersOnce)
run(){
    local start=$SECONDS log=$(mktemp) mode=${ACLIC-+O} other
    nice -n 10 root -l -b -q "${1/.cpp/.cpp$mode}" 2>&1 | tee "$log"
    if grep -q "Failed to materialize symbols" "$log"; then # this ROOT cannot load the macro in that mode: rerun it in the other one (compiled <-> interpreted)
        [ -n "$mode" ] && other= || other=+O
        echo ">>> $1: Failed to materialize symbols, rerunning ${other:+compiled (ACLiC)}${other:-interpreted}"
        nice -n 10 root -l -b -q "${1/.cpp/.cpp$other}"
    fi
    rm -f "$log"
    echo ">>> $1 took $((SECONDS - start)) s"
}

