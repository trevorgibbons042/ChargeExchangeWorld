# BDT options test: TmvaCrossVal.cpp once per option set, the same subsample rows each time. Run by Analysis.pbs with -v RUN=<run>,SCRIPT=run_tmva_cv_opts.sh
# Writes crossval_opts/<label>.txt per run, the ROC integral per fold. Compare to crossval/crossvalidation_summary.txt (MinNodeSize=10%, MaxDepth=3).

source ./run_lib.sh # run <macro>: times it, ACLiC or interpreted

base='!H:!V:NTrees=850:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20'
for cfg in m1d3:1:3 m1d4:1:4 m1d5:1:5; do
    label=${cfg%%:*}; rest=${cfg#*:}; mns=${rest%%:*}; depth=${rest#*:}
    run "TmvaCrossVal.cpp(5, 200000, 20000, \"\", \"\", \"$label\", \"$base:MinNodeSize=${mns}%:MaxDepth=$depth\")" 2>&1 | tee tmva_crossval_opts_$label.log
done
echo "=== BDT options summary"
for f in crossval_opts/*.txt; do echo "$f: $(grep -E '^mean' $f)"; done | tee tmva_crossval_opts.log
