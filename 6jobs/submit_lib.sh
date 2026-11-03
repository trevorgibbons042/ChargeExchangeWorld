# Sourced by SubmitEvents.sh and SubmitPlots.sh, from the 6jobs folder, with RUN set (the run folder in build/root).
# submit <job name> <step> [job ids it waits for, colon separated] [ncpus, default the Plot.pbs line]: prints the new job id.
# The ncpus also go to the job as NCPUS, its thread count (Plot.pbs: ANA_NTHREADS), so a 128 cpu job runs 128 threads.
# CUT_ONLY (environment of the call, e.g. megacut3) goes to the job too: CutPlots.cpp then writes only the cut files whose name holds it and skips the trees and zdc parts.
# PARTS (set in the environment of the call, parts joined with +, e.g. 2:1+3:1+4:1) goes to the job as PARTS: the Plot3D parts of the plot3d step.
submit() {
    local name=$1 step=$2 after=$3 ncpus=$4
    local dep=() res=()
    [ -n "$after" ] && dep=(-W depend=${DEPKIND:-afterok}:$after) # DEPKIND=afterany: start when the jobs ended, whatever their exit code
    [ -n "$ncpus" ] && res=(-l select=1:ncpus=$ncpus)
    qsub -N "$name" -v RUN=$RUN,STEP=$step,NCPUS=${ncpus:-16}${PARTS:+,PARTS=$PARTS}${CUT_ONLY:+,CUT_ONLY=$CUT_ONLY} "${dep[@]}" "${res[@]}" Plot.pbs
}
