source ~/.zshrc

# runs one macro and prints how long it took. Compiled with ACLiC (+O: optimized, rebuilt only when the macro or a header changed);
# ACLIC= (empty) runs it interpreted instead, for a ROOT whose compiled macros do not load (the Mac one: runStaticInitializersOnce)
run(){
    local start=$SECONDS log=$(mktemp) mode=${ACLIC-+O} other
    root -l -b -q "${1/.cpp/.cpp$mode}" 2>&1 | tee "$log"
    if grep -q "Failed to materialize symbols" "$log"; then # this ROOT cannot load the macro in that mode: rerun it in the other one (compiled <-> interpreted)
        [ -n "$mode" ] && other= || other=+O
        echo ">>> $1: Failed to materialize symbols, rerunning ${other:+compiled (ACLiC)}${other:-interpreted}"
        root -l -b -q "${1/.cpp/.cpp$other}"
    fi
    rm -f "$log"
    echo ">>> $1 took $((SECONDS - start)) s"
}

sh preflight.sh || exit 1 # skipped when run.sh already ran it (PREFLIGHT_DONE) or with PREFLIGHT=0
# Cores: CutPlots runs with all of them, each of the 4 lanes with a quarter (the macros read ANA_NTHREADS through NThreads(), FitCore.h).
# PLOT_PIPELINE=1: the macros skip BoxAndCropPlots (FitCore.h), the python steps below run once after every lane.
export PLOT_PIPELINE=1
NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
export ANA_NTHREADS=$NCPU

total=$SECONDS
# Compile every macro once, at the same time, before the lanes (ACLiC): two lanes never build the same .so, and the lanes skip the dependency check. Not with ACLIC= (interpreted).
if [ "${ACLIC-+O}" != "" ]; then
    for m in CutPlots Plot3D QualityStudy PlotAngle PlotCEFraction PlotPolarization PlotFomPolarization PlotDiffPhiCuts GraphBlocks CETCompare SecondaryCE; do
        root -l -b -q -e ".L $m.cpp${ACLIC-+O}" > /dev/null 2>&1 &
    done
    wait
fi
run 'CutPlots.cpp' # trees (cuts/ke_gt80.root, ...: the KE cut lists), zdc (cut_ZdcCE.root, cut_ZdcAllEv.root) and samples (one cut file per cut of 6_cuts) that the plot macros read, with all the threads

export ANA_NTHREADS=$(( NCPU / 4 > 0 ? NCPU / 4 : 1 ))

# Lanes run at the same time. A lane holds the macros that write into the same plots/ folders (stats.csv and fit_values.csv are
# truncated per folder) or read another's csv files, in order; no two lanes share a folder.
{ # 1_position, 2_energy, 6_cuts: the four independent parts of Plot3D (P3_PART, Plot3D.cpp) side by side, then QualityStudy
    for part in 1 3 4 2:0 2:1 2:2 2:3 2:4 2:5; do (export P3_PART=${part%%:*} P3_CUT=; [[ $part == *:* ]] && export P3_CUT=${part#*:}; run 'Plot3D.cpp') & done
    wait
    (export P3_PART=5 P3_CUT=; run 'Plot3D.cpp') # part 5: the rank_*.csv files from the rows of the part 2 processes
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
