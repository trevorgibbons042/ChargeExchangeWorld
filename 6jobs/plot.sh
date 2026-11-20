source ~/.zshrc

source ./run_lib.sh # run <macro>: times it, ACLiC or interpreted, retries in the other mode on Failed to materialize symbols

sh preflight.sh || exit 1 # skipped when run_all.sh already ran it (PREFLIGHT_DONE) or with PREFLIGHT=0
# Cores: CutPlots runs with all of them, each of the 4 lanes with a quarter (the macros read ANA_NTHREADS through NThreads(), FitCore.h).
# PLOT_PIPELINE=1: the macros skip BoxAndCropPlots (FitCore.h), the python steps below run once after every lane.
export PLOT_PIPELINE=1
NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
export ANA_NTHREADS=$NCPU
# Most root processes started at once by the compile step. The cap applies on the Mac only; elsewhere (LCRC) 1000 means no cap in practice. PLOT_MAXJOBS overrides it.
if [[ $(uname) == Darwin ]]; then MAXJOBS=${PLOT_MAXJOBS:-$(( NCPU / 2 > 0 ? NCPU / 2 : 1 ))}; else MAXJOBS=${PLOT_MAXJOBS:-1000}; fi

total=$SECONDS
# Compile every macro once, at the same time, before the lanes (ACLiC): two lanes never build the same .so, and the lanes skip the dependency check. Not with ACLIC= (interpreted).
if [ "${ACLIC-+O}" != "" ]; then
    nc=0
    for m in ${ALL_PLOTS:+CutPlots F1_2_Plot3D F3_ReconCalcVals F4_CEBlock F5_FomCalc F7_SecondaryCE F8_TrackerFits} F9_PlotReconNew ${USE_BDT:+F10_BDTCuts}; do
        nice -n 10 root -l -b -q -e ".L $m.cpp${ACLIC-+O}" > /dev/null 2>&1 # one after the other: compiling them side by side raced in ACLiC (F9 and F4 failed silently with no member RResultHandle), the errors went to /dev/null
    done
    wait
    echo ">>> compile step took $((SECONDS - total)) s"
fi
[ -n "${ALL_PLOTS:-}" ] && run 'CutPlots.cpp("zdc")' # cut_Plot3DCE.root, cut_Plot3DAllEv.root and cut_Plot3DAllPart.root that Plot3D reads, with all the threads

export ANA_NTHREADS=$(( NCPU / 4 > 0 ? NCPU / 4 : 1 )) # 5 lanes below, the cores are still split by 4

# Lanes run at the same time. A lane holds the macros that write into the same plots/ folders (stats.csv and fit_values.csv are
# truncated per folder) or read another's csv files, in order; no two lanes share a folder.
# Default: only F9_PlotReconNew (plots/6_ReconNew_1). The lanes of the other plot macros run with ALL_PLOTS=1.
if [ -n "${ALL_PLOTS:-}" ]; then
{ # 1_position, 2_energy: the three groups of Plot3D (P3_SUB 0 to 2, F1_2_Plot3D.cpp) side by side
    for sub in 0 1 2; do (export P3_SUB=$sub; run 'F1_2_Plot3D.cpp') & done
    wait
} &
{ # 4_CEblock
    run 'F4_CEBlock.cpp(50, 180, "fDEnergy")' # the block position and CE plots of plots/4_CEblock
} &
{ # 3_ReconCalcVals (CE truth vs recon, cuts/cut_APandT.root) and 5_fomCalc (FOM and purity per theta bin, cuts/cut_angles_all.root), no cuts
    run 'F3_ReconCalcVals.cpp("fDEnergy")'
    run 'F5_FomCalc.cpp'
} &
run 'F7_SecondaryCE.cpp' & # CE of the beam neutron vs the secondary neutrons by creator, plots/7_secondaryCE
run 'F8_TrackerFits.cpp' & # seed and extend tracks vs the track IDs (pure, mixed, missing), from the sim output, plots/8_trackerfits
fi
run 'F9_PlotReconNew.cpp' & # the event recon (ReconAngleAll.cpp, cuts/cut_angles_all.root), plots/6_ReconNew_1; skips itself when that file is missing. With cuts/tmva_scores.root (run_tmva.sh) it also writes the BDT cut folders of plots/10_BDTcuts
wait
[ -n "${USE_BDT:-}" ] && run 'F10_BDTCuts.cpp' # only with USE_BDT=1; plots/10_BDTcuts: P and efficiencies vs the BDT cut (the csv files F9 just wrote) and the TMVA plots of resultTMVAtraining_SigBkg.root; each part skips itself when its input is missing

# Folder and csv layout of plots/, see tidy_plots.py. Done last: the macros read and append to their csv files while they run.
python3 python/tidy_plots.py plots
python3 python/crop_pdfs.py plots # cuts the blank rest of ROOT's A4 page off every pdf
# to_png.sh (pdf to png) is not run here: run it by hand on the machine you view the plots on
cp settings.csv plots/ # the sim settings of this run beside its plots
# F9_PlotReconNew.cpp can die with no message (a 100 GB+ cut_angles_all.root, out of memory) and leave plots/6_ReconNew_1 without scan_stats.csv: say so loudly and exit 1, so the job and the run folder do not look finished
missing=0
[ -f plots/settings.csv ] || { echo "!!! MISSING plots/settings.csv"; missing=1; }
for d in 1_noCuts 2_KECUT80 3_nPlanesUsedEven 4_Theta0_1 5_avgTrackerHitsPerPlane1 6_megaCut; do
    [ -f plots/6_ReconNew_1/$d/scan_stats.csv ] || { echo "!!! MISSING plots/6_ReconNew_1/$d/scan_stats.csv (F9_PlotReconNew.cpp did not finish)"; missing=1; }
done
echo ">>> all plots took $((SECONDS - total)) s"
exit $missing
