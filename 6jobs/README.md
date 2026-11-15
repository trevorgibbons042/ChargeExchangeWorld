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

Runs `./sim <macro>` in `build/` on one 128-cpu node. The macro is a file in `build/`, for example `bash Run.sh run_0.mac`. The job is named `Run_<macro name>`. Plain `qsub Run.pbs` still works and runs `run_0.mac`. The job sets `G4FORCENUMBEROFTHREADS=128`, so the sim uses 128 threads whatever the macro's `/run/numberOfThreads` says. When the sim finishes cleanly, `stash_mac.sh` moves the macro into its run folder (`build/root/<run>/`), so `build/` keeps only the templates. `run_0.mac` and `run_nophysics.mac` always stay, and so does the macro of a failed run. `Full.pbs` does the same after its sim step.

```
bash ScanBlocks.sh <template macro> <total thickness m> [dz m]
```

Block scan. Writes one macro per block count (1, 3, 5, 7, 9) and beam energy from the template in `build/`, each with `/detector/CarbonThickness` set to the total divided by the count, and submits each with `Run.sh`. The energies are 25 GeV with a 0.92 GeV sigma, 50 GeV with 1.8 GeV and 137.5 GeV with 5 GeV, set through `/randGun/Energy` and `/randGun/EnergySD` and replacing any energy lines in the template, so one call submits 15 jobs. The run folder becomes `<folder>_<total>m/<E>GeV/<n>block` in `build/root/`, and `ScanBlockAnalysis.sh` analyses them all. dz (default 0.025 m) is the gap between the edge of a block and the edge of its tracker. It refuses if the trackers of neighbouring blocks would meet: the slot spacing 7.7 m / (n + 1) must exceed the block thickness plus 2 dz. Example: `bash ScanBlocks.sh run_0.mac 0.4`.

`COUNTS="1 2 3"` and `ENERGIES="25:0.92"` (environment, energy:sigma in GeV) replace the block counts and energies. `FULL=1` submits every run as one `Full.pbs` job, named `Full_<E>GeV_<n>block`, that runs the sim and then the whole analysis chain (`run_all.sh`: both recons, all plots) in the run folder, so a scan is one job per block count with nothing left to submit. Example, 25 GeV with blocks 1 to 11: `COUNTS="1 2 3 4 5 6 7 8 9 10 11" ENERGIES="25:0.92" FULL=1 bash ScanBlocks.sh run_801_250m.mac 0.4`.

```
bash ScanBlocksEach.sh <template macro> <thickness of each block m> [dz m]
```

Same scan as `ScanBlocks.sh`, but every block gets the given thickness, so the total grows with the count (n x thickness). The run folder becomes `<folder>_<thickness>mEach/<E>GeV/<n>block`, and `ScanBlockAnalysis.sh <folder>_<thickness>mEach` analyses them. Example: `bash ScanBlocksEach.sh run_0.mac 0.4`.

```
bash ScanBlockAnalysis.sh <scan folder in build/root> [plot job]
```

Runs `SubmitEvents.sh` on every run folder two levels below the scan folder (`<scan folder>/<E>GeV/<n>block`), found by their `settings.csv`. The second argument goes to each call. Example: `bash ScanBlockAnalysis.sh Physics_801_0.4m`.

## Analysis chain, one job per step

```
bash SubmitEvents.sh <run>
```

Submits one Analysis job that runs CutPlots events, ReconAngleAll, CalcAPandT, CutPlots count and CutPlots one after the other, then the plot jobs from `SubmitPlots.sh` (starting at Cuts). It copies `5rootMacros/` into the run folder first. The run folder needs `output0_t*.root` and `settings.csv`.

To start the chain at a plot job, give its name (`CutPlots`, `Cuts`, `CEblock`, `SecCE` or `Pdfs`) as the second argument. No analysis job is submitted, and that job and the ones after it in the list run:

```
bash SubmitEvents.sh <run> Cuts
```

This submits Cuts, CEblock, SecCE and Pdfs, and needs the `cuts/` files from CutPlots to be there already.

## Plots only

```
bash SubmitPlots.sh <run> [jobids] [first job]
```

Needs `cut_*.root` and `settings.csv` in the run folder. `jobids` is an optional colon-separated list of jobs that CutPlots should wait for (use `""` for none). `first job` is one of the plot job names above; the default is `CutPlots`, so the whole plot chain runs. A run from `CutPlots` deletes the old `plots/` and the `Plot3D_sub*.log` files of the run folder first (`cuts/` stays). A later first job keeps `plots/`.

### Redo one folder of plots/

```
bash SubmitPlots.sh <run> "" CutPlots 4_CEblock
```

The 4th argument is a folder of `plots/` (a leading `plots/` is dropped). The script deletes that folder, submits only the jobs that write it (no CutPlots; the `cuts/` files must already be there), then Pdfs. Which jobs run:

- `1_position`, `2_energy`: Plot3D.
`4_CEblock`: CEblock. `7_secondaryCE`: SecCE. These write csv files the next macro reads, so the whole lane reruns and only the named folder is deleted first.

An unknown folder is refused before anything is deleted. The new steps `plot3d` (Plot3D parts in `PARTS`, joined with `+`, for example `PARTS=2:1+3:1+4:1`) and `quality` can also be given to `Plot.pbs` by hand.

## Whole chain in a single job

```
qsub -v RUN=<run> Analysis.pbs
```

Runs `run_all.sh` on one 128-cpu node; the job copies the macros in itself.

Without ReconAngleAll (reuses `cuts/cut_angles_all.root`):

```
qsub -v RUN=<run>,SCRIPT=run_norecon.sh Analysis.pbs
```

Signal vs background BDT (`run_tmva.sh`: TmvaSplit, TmvaTrain, TmvaApply; ReconAngleAll only when `cuts/cut_angles_all.root` is missing, `FORCE_RECON=1` reruns it). Each step logs to `tmva_*.log` in the run folder; `PRESCALE` (default 20) keeps every 20th event:

```
qsub -v RUN=<run>,SCRIPT=run_tmva.sh Analysis.pbs
```

## One step by hand

```
qsub -v RUN=<run>,STEP=<step>,NCPUS=16 Plot.pbs
```

Steps: `cutevents`, `analysis`, `recon`, `reconnew`, `plotnew`, `calc`, `countce`, `cutzdc`, `plot3d`, `ceblock`, `calcplots` (3_ReconCalcVals and 5_fomCalc), `secce`, `pdfs`. Copy the macros into the run folder first.

## Options

- The macros are compiled with ACLiC by default. `-v ACLIC=` runs them interpreted instead.
- `DEPKIND=afterany` makes a job start after its dependencies end, whatever their exit code. `SubmitPlots.sh` uses it for the final Pdfs job.
- `BASE=<path>` points the scripts at another copy of the project.

## Checking jobs

```
qstat -u $USER
qdel <jobid>
```

Logs land in this folder as `<jobname>.o<jobid>`.
