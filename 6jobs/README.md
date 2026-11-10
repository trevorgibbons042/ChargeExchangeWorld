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

```
bash ScanBlocks.sh <template macro> <total thickness m> [dz m]
```

Block scan. Writes one macro per block count (1, 3, 5, 7, 9) and beam energy from the template in `build/`, each with `/detector/CarbonThickness` set to the total divided by the count, and submits each with `Run.sh`. The energies are 25 GeV with a 0.92 GeV sigma, 50 GeV with 1.8 GeV and 137.5 GeV with 5 GeV, set through `/randGun/Energy` and `/randGun/EnergySD` and replacing any energy lines in the template, so one call submits 15 jobs. The run folder becomes `<folder>_<total>m/<E>GeV/<n>block` in `build/root/`, and `ScanBlockAnalysis.sh` analyses them all. dz (default 0.025 m) is the gap between the edge of a block and the edge of its tracker. It refuses if the trackers of neighbouring blocks would meet: the slot spacing 7.7 m / (n + 1) must exceed the block thickness plus 2 dz. Example: `bash ScanBlocks.sh run_0.mac 0.4`.

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

Submits one Analysis job that runs CutEvents, ReconAngleAll, CalcAPandT, CountCE and CutPlots one after the other, then the plot jobs from `SubmitPlots.sh` (starting at Cuts). It copies `5rootMacros/` into the run folder first. The run folder needs `output0_t*.root` and `settings.csv`.

To start the chain at a plot job, give its name (`CutPlots`, `Cuts`, `Fom`, `CEblock`, `SecCE`, `Smear` or `Pdfs`) as the second argument. No analysis job is submitted, and that job and the ones after it in the list run:

```
bash SubmitEvents.sh <run> Cuts
```

This submits Cuts, Fom, CEblock, SecCE, Smear and Pdfs, and needs the `cuts/` files from CutPlots to be there already.

## Plots only

```
bash SubmitPlots.sh <run> [jobids] [first job]
```

Needs `cut_*.root` and `settings.csv` in the run folder. `jobids` is an optional colon-separated list of jobs that CutPlots should wait for (use `""` for none). `first job` is one of the plot job names above; the default is `CutPlots`, so the whole plot chain runs. A run from `CutPlots` deletes the old `plots/` and the `Plot3D_part*.log` files of the run folder first (`cuts/` stays). A later first job keeps `plots/`.

### Redo one folder of plots/

```
bash SubmitPlots.sh <run> "" CutPlots 6_cuts/3_ResCuts
```

The 4th argument is a folder of `plots/` (a leading `plots/` is dropped). The script deletes that folder, submits only the jobs that write it (no CutPlots; the `cuts/` files must already be there), then Pdfs. Which jobs run:

- `1_position`, `2_energy`: Plot3D part 1.
- `6_cuts`: Plot3D parts 2, 3 and 4 for every cut, part 5 (the `rank_*.csv` files), QualityStudy.
- `6_cuts/<cut folder>` (for example `6_cuts/6_megaCut3` or `6_cuts/6_megaCut3/1_1particlecut`): parts 2, 3 and 4 of that cut only, then part 5. QualityStudy also runs for `3_ResCuts`. The cut numbers are the `P3_CUT` numbers of `Plot3D.cpp`, listed in `cut_indexes` in `SubmitPlots.sh`.
- `3_ReconCalcVals`, `5_fomCalc`: the Fom jobs. `4_CEblock`: CEblock. `7_secondaryCE`: SecCE. `8_smear`: Smear. These write csv files the next macro reads, so the whole lane reruns and only the named folder is deleted first.

An unknown folder is refused before anything is deleted. The new steps `plot3d` (Plot3D parts in `PARTS`, joined with `+`, for example `PARTS=2:1+3:1+4:1`) and `quality` can also be given to `Plot.pbs` by hand.

## Whole chain in a single job

```
qsub -v RUN=<run> Analysis.pbs
```

Runs `run.sh` on one 128-cpu node; the job copies the macros in itself.

Without ReconAngleAll (reuses `cuts/cut_angles_all.root`):

```
qsub -v RUN=<run>,SCRIPT=run_norecon.sh Analysis.pbs
```

## One step by hand

```
qsub -v RUN=<run>,STEP=<step>,NCPUS=16 Plot.pbs
```

Steps: `cutevents`, `recon`, `calc`, `countce`, `cuttrees`, `cuts`, `fom`, `ceblock`, `secce`, `smear`, `pdfs`. Copy the macros into the run folder first.

## 6_cuts only

```
bash Submit6Cuts.sh <run> [CutPlots | Plots]
```

Rebuilds `plots/6_cuts` only, from the sim output, without the rest of the chain. Jobs: Analysis6 (CutEvents, ReconAngleAll, CalcAPandT true KE only, CutPlots zdc and samples), then the Plot3D parts of `6_cuts` and QualityStudy (side by side, both after Analysis6), then Pdfs. It skips CountCE, the detector-KE pass of CalcAPandT, the CutPlots `trees` part, and every plot lane outside `6_cuts`. Other folders of `plots/` stay as they are; `plots/6_cuts` is deleted first.

- `CutPlots`: skips CutEvents, ReconAngleAll and CalcAPandT (their `cuts/` files must be there, including `cut_APandT.root`) and starts at CutPlots.
- `Plots`: no analysis at all, only the 6_cuts plot jobs (all the `cuts/` files must be there).

A `cut_APandT.root` written by `Analysis6` has no `_DE` trees, so the Fom lane needs a full `calc` step before it runs.

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
