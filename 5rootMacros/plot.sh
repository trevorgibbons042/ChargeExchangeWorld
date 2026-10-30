source ~/.zshrc

# runs one macro and prints how long it took. Compiled with ACLiC (+O: optimized, rebuilt only when the macro or a header changed);
# ACLIC= (empty) runs it interpreted instead, for a ROOT whose compiled macros do not load (the Mac one: runStaticInitializersOnce)
run(){
    local start=$SECONDS
    root -l -b -q "${1/.cpp/.cpp${ACLIC-+O}}"
    echo ">>> $1 took $((SECONDS - start)) s"
}

sh preflight.sh || exit 1 # skipped when run.sh already ran it (PREFLIGHT_DONE) or with PREFLIGHT=0
# Cores: CutTrees runs with all of them, each of the 4 lanes with a quarter (the macros read ANA_NTHREADS through NThreads(), FitCore.h).
# PLOT_PIPELINE=1: the macros skip BoxAndCropPlots (FitCore.h), the python steps below run once after every lane.
export PLOT_PIPELINE=1
NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
export ANA_NTHREADS=$NCPU

total=$SECONDS
run 'CutTrees.cpp' # cuts/ke_gt80.root: the shared 80% KE cut the plot macros read (CutFiles.h)

# Compile every macro once, at the same time, before the lanes (ACLiC): two lanes never build the same .so, and the lanes skip the dependency check. Not with ACLIC= (interpreted).
if [ "${ACLIC-+O}" != "" ]; then
    for m in ZdcColumns CutSamples Plot3D QualityStudy PlotAngle PlotCEFraction PlotPolarization PlotFomPolarization PlotDiffPhiCuts GraphBlocks CETCompare SecondaryCE; do
        root -l -b -q -e ".L $m.cpp${ACLIC-+O}" > /dev/null 2>&1 &
    done
    wait
fi
run 'ZdcColumns.cpp' # cuts/cut_ZdcCE.root and cut_ZdcAllEv.root (the ZDC protons with the Plot3D columns) that Plot3D reads, with all the threads
run 'CutSamples.cpp' # one cut file per cut of 6_cuts (cut_ResCuts.root, cut_zOrder.root, cut_KECUT80_CE.root, ...) that Plot3D reads, with all the threads

export ANA_NTHREADS=$(( NCPU / 4 > 0 ? NCPU / 4 : 1 ))

# Lanes run at the same time. A lane holds the macros that write into the same plots/ folders (stats.csv and fit_values.csv are
# truncated per folder) or read another's csv files, in order; no two lanes share a folder.
{ # 1_position, 2_energy, 6_cuts: the four independent parts of Plot3D (P3_PART, Plot3D.cpp) side by side, then QualityStudy
    for part in 1 2 3 4; do (export P3_PART=$part; run 'Plot3D.cpp') & done
    wait
    run 'QualityStudy.cpp'
} &
{ # 3_ReconCalcVals, 5_fomCalc: PlotCEFraction and PlotFomPolarization read the FOM tables PlotAngle writes (FomCommon.h)
    run 'PlotAngle.cpp("fDEnergy")'
    run 'PlotCEFraction.cpp'
    run 'PlotPolarization.cpp("fDEnergy")'
    run 'PlotFomPolarization.cpp("fDEnergy")'
} &
{ # 4_CEblock
    run 'PlotDiffPhiCuts.cpp(-180, -50, "fDEnergy")'
    run 'GraphBlocks.cpp'
    run 'CETCompare.cpp' # t and z of the CE in the carbon blocks vs the silicon trackers, plots/4_CEblock/1_CarbonVsSi
} &
run 'SecondaryCE.cpp' & # CE of the beam neutron vs the secondary neutrons by creator, plots/7_secondaryCE
wait

# Folder and csv layout of plots/, see tidy_plots.py. Done last: the macros read and append to their csv files while they run.
python3 tidy_plots.py plots
python3 cut_side_box.py plots # side box with the cut numbers on every 6_cuts pdf, after tidy_plots.py (the pdfs are final)
python3 crop_pdfs.py plots # cuts the blank rest of ROOT's A4 page off every pdf, the side box included
bash to_png.sh plots # last: every final pdf becomes a png
echo ">>> all plots took $((SECONDS - total)) s"
