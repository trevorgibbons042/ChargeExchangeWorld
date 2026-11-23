# Drop-one test of the signal vs background BDT inputs on one run folder: TmvaCrossVal.cpp once per input left out, plus "none" as the baseline. Run by Analysis.pbs with -v RUN=<run>,SCRIPT=run_tmva_cv_drop.sh
# Reuses cuts/cut_tmva_signal.root and cut_tmva_background.root of an earlier run_tmva.sh. One run after the other: every run rebuilds the same subsample files (the same rows each time, so the folds are comparable).
# Writes crossval_drop/<input>.txt per run and the means to tmva_crossval_drop.log.

source ./run_lib.sh # run <macro>: times it, ACLiC or interpreted

for v in none lastXD lastYD nPlanesUsed nChargedParticle zdcKEAllD outDTheta; do
    run "TmvaCrossVal.cpp(5, 200000, 20000, \"\", \"$v\")" 2>&1 | tee tmva_crossval_drop_$v.log
done
echo "=== drop-one summary"
for f in crossval_drop/*.txt; do echo "$f: $(grep -E '^mean' $f)"; done | tee tmva_crossval_drop.log
