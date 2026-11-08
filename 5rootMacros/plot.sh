source ~/.zshrc

source ./run_lib.sh # run <macro>: times it, ACLiC or interpreted, retries in the other mode on Failed to materialize symbols

sh preflight.sh || exit 1 # skipped when run.sh already ran it (PREFLIGHT_DONE) or with PREFLIGHT=0
# Cores: CutPlots runs with all of them, each of the 5 lanes with a fifth (the macros read ANA_NTHREADS through NThreads(), FitCore.h).
# PLOT_PIPELINE=1: the macros skip BoxAndCropPlots (FitCore.h), the python steps below run once after every lane.
export PLOT_PIPELINE=1
NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
export ANA_NTHREADS=$NCPU
# Most root processes started at once by the compile step and the Plot3D parts. Each Plot3D part loads the full trees, so starting all of them together
# overloads the CPU and fills the memory of a laptop. The cap applies on the Mac only; elsewhere (LCRC) 1000 means no cap in practice. PLOT_MAXJOBS overrides it.
if [[ $(uname) == Darwin ]]; then MAXJOBS=${PLOT_MAXJOBS:-$(( NCPU / 2 > 0 ? NCPU / 2 : 1 ))}; else MAXJOBS=${PLOT_MAXJOBS:-1000}; fi

total=$SECONDS
# Compile every macro once, at the same time, before the lanes (ACLiC): two lanes never build the same .so, and the lanes skip the dependency check. Not with ACLIC= (interpreted).
if [ "${ACLIC-+O}" != "" ]; then
    nc=0
    for m in CutPlots Plot3D QualityStudy PlotAngle PlotCEFraction PlotPolarization PlotFomPolarization PlotDiffPhiCuts GraphBlocks CETCompare SecondaryCE; do
        nice -n 10 root -l -b -q -e ".L $m.cpp${ACLIC-+O}" > /dev/null 2>&1 &
        (( ++nc % MAXJOBS == 0 )) && wait
    done
    wait
fi
if (( NCPU >= 32 )); then
    # the same split as 6jobs/Plot.pbs (cp_run): trees, zdc and the samples without a ZDC KE cut side by side, then the samples that read cut_ZdcAllEv.root (megacut3, the KE scan 85 90 95, hitZDC) as five processes, each with its share of the cores.
    # One process for all of them took 590 s of the 100M event 25 GeV 5block run, the split ones 26 to 160 s each. Fewer cores (the Mac): the single process below.
    cpStart=$SECONDS
    cp_run(){ # name threads part CUT_ONLY
        ( ANA_NTHREADS=$2 CUT_ONLY=$4 nice -n 10 root -l -b -q "CutPlots.cpp${ACLIC-+O}(\"$3\")" > CutPlots_$1.log 2>&1; echo ">>> CutPlots $1 exit $? after $((SECONDS - cpStart)) s" ) &
    }
    cp_run trees $(( NCPU / 32 )) trees ""
    cp_run zdc $(( NCPU / 2 )) zdc ""
    cp_run samplesA $(( NCPU / 2 )) samples "ResCuts,nParticles,nZdc1,nHits,nChHits"
    wait
    cp_run samplesB1 $(( NCPU / 5 )) samples "cut_hitZDC"
    cp_run samplesB2 $(( NCPU / 5 )) samples "cut_megacut3_1.root"
    cp_run samplesB3 $(( NCPU / 5 )) samples "_ke85"
    cp_run samplesB4 $(( NCPU / 5 )) samples "_ke90"
    cp_run samplesB5 $(( NCPU / 5 )) samples "_ke95"
    wait
    tail -qn 1 CutPlots_*.log
    echo ">>> CutPlots.cpp took $((SECONDS - cpStart)) s"
else
    run 'CutPlots.cpp' # trees (cuts/ke_gt80.root, ...: the KE cut lists), zdc (cut_ZdcCE.root, cut_ZdcAllEv.root) and samples (one cut file per cut of 6_cuts) that the plot macros read, with all the threads
fi

export ANA_NTHREADS=$(( NCPU / 5 > 0 ? NCPU / 5 : 1 )) # 5 lanes below

# Lanes run at the same time. A lane holds the macros that write into the same plots/ folders (stats.csv and fit_values.csv are
# truncated per folder) or read another's csv files, in order; no two lanes share a folder.
{ # 1_position, 2_energy, 6_cuts: the independent parts of Plot3D (P3_PART, P3_CUT, Plot3D.cpp) side by side, parts 2, 3 and 4 once per cut (4 also for the KE scan sets, 9 to 11; 3 and 4 of cut 0 once per folder, P3_SUB), then QualityStudy
    (export P3_PART=0 P3_CUT= P3_SUB=; run 'Plot3D.cpp') # cuts/plot3d_totals.csv once (Plot3D.cpp RawTotals), then the parts read it
    n=0
    for part in 1::0 1::1 1::2 1::3 3:0:0 3:0:1 3:0:2 3:0:3 3:1 3:2 3:3 3:4 3:5 3:6 3:7 3:8 4:6 4:9 4:10 4:11 2:0 2:1 2:2 2:3 2:4 2:5 2:6 2:7 2:8; do (export P3_PART=${part%%:*} P3_CUT= P3_SUB=; [[ $part == *:* ]] && { export P3_CUT=${part#*:}; [[ $P3_CUT == *:* ]] && { export P3_SUB=${P3_CUT#*:}; export P3_CUT=${P3_CUT%%:*}; }; }; run 'Plot3D.cpp') & (( ++n % MAXJOBS == 0 )) && wait; done
    wait
    (export P3_PART=5 P3_CUT=; run 'Plot3D.cpp') # part 5: the rank_*.csv files from the rows of the part 2 processes
    run 'PlotKEOverlay.cpp' # the KE overlays of the megaCut folder (6_cuts/9_megaCut3/KE_overlay), from the csv files of part 4, before tidy_plots.py moves them
    run 'QualityStudy.cpp'
} &
{ # 3_ReconCalcVals, 5_fomCalc: PlotCEFraction and PlotFomPolarization read the FOM tables PlotAngle writes (FomCommon.h)
    run 'PlotAngle.cpp("fDEnergy")'
    run 'PlotCEFraction.cpp'
    run 'PlotPolarization.cpp("fDEnergy")'
    run 'PlotFomPolarization.cpp("fDEnergy")'
} &
{ # 4_CEblock
    run 'PlotDiffPhiCuts.cpp(50, 180, "fDEnergy")'
    run 'GraphBlocks.cpp'
    run 'CETCompare.cpp' # t and z of the CE in the carbon blocks vs the silicon trackers, plots/4_CEblock/1_CarbonVsSi
} &
run 'SmearMatrix.cpp' & # response matrices of RECONTHETABL_ESTNEU and phi vs the truth relative to the neutron, plots/8_smear: a lane of its own, it reads only cut_angles_all.root and cut_APandT.root and writes only plots/8_smear
run 'SecondaryCE.cpp' & # CE of the beam neutron vs the secondary neutrons by creator, plots/7_secondaryCE
wait
zsh smear_cuts.sh # the SmearMatrix pages 2_matrix2 and 3_matrix3 of every 6_cuts folder's sample (1_additionalinfo); needs the cuts.txt of Plot3D part 2, so after the lanes

# Folder and csv layout of plots/, see tidy_plots.py. Done last: the macros read and append to their csv files while they run.
python3 tidy_plots.py plots
python3 cut_side_box.py plots # side box with the cut numbers on every 6_cuts pdf, after tidy_plots.py (the pdfs are final)
python3 crop_pdfs.py plots # cuts the blank rest of ROOT's A4 page off every pdf, the side box included
# to_png.sh (pdf to png) is not run here: run it by hand on the machine you view the plots on
cp settings.csv plots/ # the sim settings of this run beside its plots
echo ">>> all plots took $((SECONDS - total)) s"
