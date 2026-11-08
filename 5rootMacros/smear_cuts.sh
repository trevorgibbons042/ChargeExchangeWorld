# The response matrices of SmearMatrix.cpp for the sample of each 6_cuts folder (SMEAR_CUT, SmearMatrix.cpp): pages 2_matrix2 and 3_matrix3 into plots/6_cuts/<folder>/1_additionalinfo/.
# zsh smear_cuts.sh [folder under 6_cuts ...], no argument: every folder with a 1_additionalinfo. Run after Plot3D.cpp (cuts.txt of each folder names the cut in the titles) and the cuts/ files of
# run.sh, before tidy_plots.py. SMEAR_JOBS (default 6) processes at a time, each with its share of the cores. Runs in zsh (plot.sh, plot_cut.sh) and bash (6jobs/Plot.pbs, step smearcuts); the environment
# (modules, ACLIC) is the caller's: ACLIC unset compiles with ACLiC (+O), ACLIC= runs interpreted.
source ./run_lib.sh # run <macro>: times it, ACLiC or interpreted, retries in the other mode on Failed to materialize symbols

dirs=("$@")
(( $# )) || dirs=(1_hitZDC 2_KECUT/KECUT80 2_KECUT/KECUT85 2_KECUT/KECUT90 2_KECUT/KECUT95 3_ResCuts 4_nParticles/1_particleCut 4_nParticles/2_particleCut 5_nHits 8_nZdcParticles/1_1zdcparticlecut 9_megaCut3/1_1particlecut 9_megaCut3/1_1particlecut_KE85 9_megaCut3/1_1particlecut_KE90 9_megaCut3/1_1particlecut_KE95 10_nChargedTrackerHits/1_noExtraHits 10_nChargedTrackerHits/2_maxOneExtraHit)

export PLOT_PIPELINE=1
NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
nj=${SMEAR_JOBS:-6}
export ANA_NTHREADS=$(( NCPU / nj > 0 ? NCPU / nj : 1 ))
total=$SECONDS
# build the macro once before the processes (ACLiC): two processes never build the same .so. Not with ACLIC= (interpreted).
if [ "${ACLIC-+O}" != "" ]; then nice -n 10 root -l -b -q -e ".L SmearMatrix.cpp${ACLIC-+O}" > /dev/null 2>&1; fi
n=0
for d in "${dirs[@]}"; do
    (export SMEAR_CUT=$d; run 'SmearMatrix.cpp') &
    (( ++n % nj == 0 )) && wait
done
wait
# the unfolding of the megaCut3 phi cosfits (and of the full sample matrix) needs the matrices written above: after them, one process with all the cores for its toys
case " ${dirs[*]} " in *" 9_megaCut3/1_1particlecut "*) ANA_NTHREADS=$NCPU run 'UnfoldPhi.cpp' ;; esac
echo ">>> smear_cuts ($n folders) took $((SECONDS - total)) s"
