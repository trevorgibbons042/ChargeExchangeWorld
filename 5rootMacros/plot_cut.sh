source ~/.zshrc

# Plots of one 6_cuts folder without the rest of Plot3D: zsh plot_cut.sh <i>, with i the P3_CUT number of Plot3D.cpp
# (0: 1_hitZDC and 2_KECUT80, 1 to 9: sampleCuts in order, 9 is the mega cut, 11_megaCut/1_megaCut). Needs the cuts/ files of run.sh (ZdcColumns,
# CutSamples, CutTrees) and, for the rank_*.csv files, cuts/rank_rows_0 to 8.tsv of an earlier full plot.sh. Same ACLIC switch as plot.sh.
i=${1:?usage: zsh plot_cut.sh <cut number, 0 to 9>}
export ACLIC=${ACLIC-} # interpreted by default, like Plot.pbs: ACLiC-compiled macros fail to load on the cluster ROOT (runStaticInitializersOnce), and three processes would build the same .so. ACLIC=+O compiles.

run(){
    local start=$SECONDS
    root -l -b -q "${1/.cpp/.cpp${ACLIC-+O}}"
    echo ">>> $1 took $((SECONDS - start)) s"
}

sh preflight.sh || exit 1
export PLOT_PIPELINE=1
NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
export ANA_NTHREADS=$NCPU
total=$SECONDS

run 'CutSamples.cpp' # rewrites every cuts/cut_*.root of 6_cuts in one pass, cut_megacut.root included

export ANA_NTHREADS=$(( NCPU / 3 > 0 ? NCPU / 3 : 1 ))
for part in 2 3 4; do (export P3_PART=$part P3_CUT=$i; run 'Plot3D.cpp') & done # cut_ratios.csv and the rank row, 1_additionalinfo, 2_phi_cosfit
wait
(export P3_PART=5 P3_CUT=; run 'Plot3D.cpp') # rank_*.csv from cuts/rank_rows_0 to 9.tsv

python3 tidy_plots.py plots
python3 cut_side_box.py plots # side box with the cut numbers on every 6_cuts pdf, after tidy_plots.py
python3 crop_pdfs.py plots # last of the pdf edits
echo ">>> cut $i took $((SECONDS - total)) s"
