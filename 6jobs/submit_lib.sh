# Sourced by SubmitEvents.sh and SubmitPlots.sh, from the 6jobs folder, with RUN set (the run folder in build/root).
# submit <job name> <step> [job ids it waits for, colon separated] [ncpus, default the Plot.pbs line]: prints the new job id.
# The ncpus also go to the job as NCPUS, its thread count (Plot.pbs: ANA_NTHREADS), so a 128 cpu job runs 128 threads.
submit() {
    local name=$1 step=$2 after=$3 ncpus=$4
    local dep=() res=()
    [ -n "$after" ] && dep=(-W depend=afterok:$after)
    [ -n "$ncpus" ] && res=(-l select=1:ncpus=$ncpus)
    qsub -N "$name" -v RUN=$RUN,STEP=$step,NCPUS=${ncpus:-16} "${dep[@]}" "${res[@]}" Plot.pbs
}
