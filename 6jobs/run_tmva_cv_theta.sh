# CV of the BDT with the current TmvaVars.h inputs and MinNodeSize=1%, MaxDepth=3. Run by Analysis.pbs with -v RUN=<run>,SCRIPT=run_tmva_cv_theta.sh
# TmvaSplit first, because the split files need a column for every input of TmvaVars.h (a new input is not in the old ones). Writes crossval_opts/theta_m1d3.txt.

source ./run_lib.sh # run <macro>: times it, ACLiC or interpreted

run "TmvaSplit.cpp(\"cuts/cut_angles_all.root\", 2)" 2>&1 | tee tmva_split.log
opts='!H:!V:NTrees=850:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:MinNodeSize=1%:MaxDepth=3'
run "TmvaCrossVal.cpp(5, 200000, 20000, \"\", \"\", \"theta_m1d3\", \"$opts\")" 2>&1 | tee tmva_crossval_opts_theta_m1d3.log
echo "=== theta_m1d3"; grep -E '^(ROC|mean|TMVA)' crossval_opts/theta_m1d3.txt | tee tmva_crossval_theta.log
