#source ~/.zshrc

source ./run_lib.sh # run <macro>: times it, ACLiC or interpreted, retries in the other mode on Failed to materialize symbols

NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}

sh preflight.sh || exit 1 # ROOT features, compile test, writable cuts/ plots/, disk; PREFLIGHT=0 skips
export PREFLIGHT_DONE=1 # plot.sh below does not repeat it

mkdir -p cuts # every cut_*.root and the KE cut lists
export ANA_NTHREADS=$NCPU
# One after the other, all the cores each: ReconAngleAll alone on 8 threads took 39 s and CutPlots count (it reads only cut_CE.root) 4 s, against 61 s for the two side by side on half of the cores each
# ReconAngleAll reads only the sim output, so it goes first; CutPlots count needs the cut_CE.root of CutPlots events
run 'ReconAngleAll.cpp'
run 'CutPlots.cpp("events")' # writes the cut files cuts/cut_*.root (CutFiles.h)
run 'CutPlots.cpp("count")'
run 'CalcAPandT.cpp' # the detector KE trees of cut_APandT.root
unset ANA_NTHREADS
zsh plot.sh
# the all-recon-particles cut file is only read by the 3D plots (plot.sh just drew them), the biggest file of the run after cut_angles_all.root: dropped by default
rm -f cuts/cut_Plot3DAllPart.root
