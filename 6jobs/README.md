# 6jobs: running on LCRC

PBS scripts for the simulation and the analysis chain. Queue `backfill`, account `eic`. Run everything from this folder on an LCRC login node:

```
cd /lcrc/project/eic/tgibbons/G4ChargeExchange/6jobs
```

`<run>` is a folder name inside `build/root`, for example `Physics_801`.

## Simulation

```
bash Run.sh <macro>
```

Runs `./sim <macro>` in `build/` on one 128-cpu node. The macro is a file in `build/`, for example `bash Run.sh run_0.mac`. The job is named `Run_<macro name>`. Plain `qsub Run.pbs` still works and runs `run_0.mac`. The job sets `G4FORCENUMBEROFTHREADS=128`, so the sim uses 128 threads whatever the macro's `/run/numberOfThreads` says.

## Analysis chain, one job per step

```
bash SubmitEvents.sh <run>
```

Submits CutEvents, ReconAngleAll, CalcAPandT and CountCE as dependent jobs, then the plot jobs from `SubmitPlots.sh`. It copies `5rootMacros/` into the run folder first. The run folder needs `output0_t*.root` and `settings.csv`.

If `cuts/cut_angles_all.root` already exists from an earlier ReconAngleAll run, skip the refit:

```
bash SubmitEvents.sh <run> cutdone
```

ConvertAngles then rewrites the file into the one-tree format. The old file is kept as `cut_angles_all_old5.root`.

## Plots only

```
bash SubmitPlots.sh <run> [jobids]
```

Needs `cut_*.root` and `settings.csv` in the run folder. `jobids` is an optional colon-separated list of jobs that CutPlots should wait for.

## Whole chain in a single job

```
rsync -a --exclude __pycache__ ../5rootMacros/ ../build/root/<run>/
qsub -v RUN=<run> Analysis.pbs
```

Runs `run.sh` on one 128-cpu node. Copy the macros in first, as above.

Without ReconAngleAll (reuses `cuts/cut_angles_all.root`):

```
qsub -v RUN=<run> AnalysisNoRecon.pbs
```

This one copies the macros itself.

## One step by hand

```
qsub -v RUN=<run>,STEP=<step>,NCPUS=16 Plot.pbs
```

Steps: `cutevents`, `recon`, `convert`, `calc`, `countce`, `cuttrees`, `cuts`, `fom`, `ceblock`, `secce`, `pdfs`. Copy the macros into the run folder first.

## Options

- `-v ACLIC=+O` compiles the macros with ACLiC. The default is interpreted, because ROOT 6.32 on the cluster fails to load compiled macros.
- `DEPKIND=afterany` makes a job start after its dependencies end, whatever their exit code. `SubmitPlots.sh` uses it for the final Pdfs job.
- `BASE=<path>` points the scripts at another copy of the project.

## Checking jobs

```
qstat -u $USER
qdel <jobid>
```

Logs land in this folder as `<jobname>.o<jobid>`.

## Notes

`CutEvent.pbs` only loads ROOT and prints host info. It runs no macro.
