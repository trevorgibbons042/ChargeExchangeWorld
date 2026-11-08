source ~/.zshrc

# Plots of one 6_cuts folder without the rest of Plot3D: zsh plot_cut.sh <i>, with i the P3_CUT number of Plot3D.cpp
# (0: 1_hitZDC and 2_KECUT/KECUT80 to KECUT95, 1 to 8: sampleCuts in order, 5 is 8_nZdcParticles/1_1zdcparticlecut, 6 is 9_megaCut3/1_1particlecut). Needs the cuts/ files of run.sh (CutPlots.cpp:
# trees, zdc, samples) and, for the rank_*.csv files, cuts/rank_rows_0 to 8.tsv of an earlier full plot.sh. Same ACLIC switch as plot.sh.
i=${1:?usage: zsh plot_cut.sh <cut number, 0 to 8>}
export ACLIC=${ACLIC-} # interpreted by default, like Plot.pbs: ACLiC-compiled macros fail to load on the cluster ROOT (runStaticInitializersOnce), and three processes would build the same .so. ACLIC=+O compiles.

source ./run_lib.sh # run <macro>: times it, ACLiC or interpreted, retries in the other mode on Failed to materialize symbols

sh preflight.sh || exit 1
export PLOT_PIPELINE=1
NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
export ANA_NTHREADS=$NCPU
total=$SECONDS

run 'CutPlots.cpp("samples")' # rewrites every cuts/cut_*.root of 6_cuts in one pass

export ANA_NTHREADS=$(( NCPU / 3 > 0 ? NCPU / 3 : 1 ))
for part in 2 3 4; do (export P3_PART=$part P3_CUT=$i; run 'Plot3D.cpp') & done # cut_ratios.csv and the rank row, 1_additionalinfo, 2_phi_cosfit
wait
(export P3_PART=5 P3_CUT=; run 'Plot3D.cpp') # rank_*.csv from cuts/rank_rows_0 to 8.tsv
case $i in # the folders of this cut number (Plot3D.cpp sampleCuts order): their SmearMatrix pages 2_matrix2 and 3_matrix3 in 1_additionalinfo, after Plot3D wrote their cuts.txt
    0) smearDirs=(1_hitZDC 2_KECUT/KECUT80 2_KECUT/KECUT85 2_KECUT/KECUT90 2_KECUT/KECUT95) ;;
    1) smearDirs=(3_ResCuts) ;;
    2) smearDirs=(4_nParticles/1_particleCut) ;;
    3) smearDirs=(4_nParticles/2_particleCut) ;;
    4) smearDirs=(5_nHits) ;;
    5) smearDirs=(8_nZdcParticles/1_1zdcparticlecut) ;;
    6) smearDirs=(9_megaCut3/1_1particlecut 9_megaCut3/1_1particlecut_KE85 9_megaCut3/1_1particlecut_KE90 9_megaCut3/1_1particlecut_KE95) ;;
    7) smearDirs=(10_nChargedTrackerHits/1_noExtraHits) ;;
    8) smearDirs=(10_nChargedTrackerHits/2_maxOneExtraHit) ;;
esac
zsh smear_cuts.sh "${smearDirs[@]}"

python3 tidy_plots.py plots
python3 cut_side_box.py plots # side box with the cut numbers on every 6_cuts pdf, after tidy_plots.py
python3 crop_pdfs.py plots # last of the pdf edits
echo ">>> cut $i took $((SECONDS - total)) s"
