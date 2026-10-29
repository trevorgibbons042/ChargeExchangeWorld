source ~/.zshrc

# runs one macro and prints how long it took. Compiled with ACLiC (+O: optimized, rebuilt only when the macro or a header changed);
# ACLIC= (empty) runs it interpreted instead, for a ROOT whose compiled macros do not load (the Mac one: runStaticInitializersOnce)
run(){
    local start=$SECONDS
    root -l -b -q "${1/.cpp/.cpp${ACLIC-+O}}"
    echo ">>> $1 took $((SECONDS - start)) s"
}

total=$SECONDS
run 'CutTrees.cpp' # cuts/<ke>.root: the shared KE and block residual cuts the plot macros read (CutFiles.h)

# Lanes run at the same time. A lane holds the macros that write into the same plots/ folders (stats.csv and fit_values.csv are
# truncated per folder) or read another's csv files, in order; no two lanes share a folder.
{ # 1_position, 2_energy, 6_cuts
    run 'Plot3D.cpp'
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
