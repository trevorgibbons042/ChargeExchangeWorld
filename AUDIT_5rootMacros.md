# 5rootMacros audit

Read-only. Nothing in `5rootMacros/` was changed for this. "Checked" means I read the code or ran a count; "suspected" means I did not confirm it. No timing was measured, so the compute items are ranked by how much work the code does, not by seconds.

## Inconsistencies

1. **Start-theta plots say "to 5 deg".** The P-vs-start-theta titles in `Plot3D.cpp` and the scan cap (`th >= 5.0`) predate the 0 to 1 degree baseline. Every sample is now at most 1 degree, so the scan ends near 0.5 to 1 degree. Checked. Fix: say "to 1 deg" and cap at 1.
2. **The 0 to 1 degree theta baseline only covers `6_cuts`.** `3_ReconCalcVals`, `4_CEblock`, `5_fomCalc` and `7_secondaryCE` still use all theta. The rank FOM in `6_cuts` still integrates theta from 0 to 15 degrees (`fomStats`) over samples that stop at 1 degree. Checked. Decide whether that is intended.
3. **Each cut is defined twice.** The `sampleCuts` expressions in `Plot3D.cpp` (used for the 80% KE rows on the ZDC trees) and the filters in `CutSamples.cpp` (used for the cut files) are separate copies of the same cut. The nHits 4/8/12/16 test and the 0 to 1 degree test appear in both, and again inside the mega cut. They have to be edited by hand together. Checked.
4. **Block phi range is not applied everywhere.** `kBlockPhiLo/Hi` (`FitCore.h`) is used by `Plot3D.cpp`, `PlotAngle.cpp` and `PlotFomPolarization.cpp`. The block slice fits in `PlotPolarization.cpp` (`FitPhi`, csv only) still include the phi = 0 spike, so their P values differ from the plotted fits. Checked.
5. **Removed variables are still carried around.** `impact`, `impactD`, `r2`, `r2D`, `resHalfRatioD`, `firstXD`, `firstYD` are still in `CutSamples.cpp` `rawCols`, `Plot3DCols.h` (17 positional values and the alias list), the `QualityStudy.cpp` column list, and every ZDC cut file. `Plot3D.cpp` and `QualityStudy.cpp` no longer plot them. The positional `rawVal(..., k)` indices make removal risky, so this needs a careful pass. Checked.
6. **`QualityStudy.cpp` still computes four variables, plots one.** The joint threshold search, `names`/`labels`, the impact equation text and the "all four together" printout still cover impact. `thresholds.csv` now has one row. Checked.
7. **`plot.sh`, `plot_cut.sh` and `6jobs/Plot.pbs` repeat the lane list and `run()`.** The retry on `Failed to materialize symbols` is in the first two only. `Plot.pbs` calls `root` directly and has no retry. Checked.
8. **Missing inputs can give empty plots with no message.** `PlotCEFraction.cpp` now reports a missing FOM csv, but the same read in `FomCommon.h:36` and the `cut_APandT.root` readers do not. This is how the empty `5_fomCalc` graphs happened. Checked.

## Wasted computing

1. **Every Plot3D process repeats full passes.** `Plot3D.cpp` runs as 9 processes in `plot.sh`. Each one counts `nCERecoAll`/`nRecoAll` (all ZDC protons) and `nCEreconTotal` (all 1.8M fitted protons) before the part check, even when the part does not use them. Checked in the code; the cost is not measured. Fix: compute inside the parts that need them, or write them once to a small csv.
2. **Wide reads for few columns.** `QualityStudy.cpp` reads about 12 columns of `cut_angles_all.root` (368 MB) and uses about 7. `CutSamples.cpp` writes 23 columns into every cut file. Checked.
3. **`cut_side_box.py` and `crop_pdfs.py` each render every pdf through `gs`.** That is two `gs` runs per pdf per pass, repeated on every rerun even when one macro changed. Suspected to dominate post-processing time. Fix: skip pdfs whose mtime did not change.
4. **`run.sh` runs `CalcAPandT.cpp` twice (plain and `"fDEnergy"`).** `plot.sh` only runs the `fDEnergy` variants of `PlotAngle`, `PlotPolarization` and `PlotFomPolarization`. I did not confirm that nothing else reads the plain trees (`CalcVals*` without `_DE`). Suspected.
5. **Small items.** `Plot3D.cpp` opens `cut_CE.root` only to test whether `fDEnergy` exists. `TestDrawPlan.cpp` and `DrawPlan.h` tests are copied into every run folder. `CutTrees.cpp` is rerun on every `plot.sh`.

## Not a problem

`ReadSettings()` is memoised by file stamp (`FitCore.h`), so calls inside loops do not reread the file.

## Run parameters: beam energy, EnergySD, beam angle and vertex, nBlocks

Checked by reading the code. The load check (`.L` of Plot3D, PlotAngle, SecondaryCE, PlotDiffPhiCuts, PlotPolarization) gives 0 errors, and `KeAxisMax()` and `ThetaAxisMax()` give 30000 MeV and 2.3 deg on a 25 GeV/c run, the old fixed values. No plot was rerun at 50 or 137.5 GeV yet.

### Already follows the settings
- Beam momentum, KE and momentum SD come from `settings.csv` (`kBeamMomentum`, `kBeamKE`, `kBeamMomentumSD`). t, A_P and the beam error fractions follow them. The 80% KE cut uses `kBeamKE` (`CutFiles.h`).
- Block and tracker z come from `StackZ` (`FitCore.h`), used by `CarbonCeCut`, `CountCE.cpp`, `BlockLimits.h` and `CETCompare.cpp`. The even-hit cuts use `MaxTrackerHits()` (2 * nBlocks), and the `GraphBlocks.cpp` grid is sized from the number of `Block` trees.
- The ZDC cell bins read `nCols_detector` and `nRows_detector`.

### Fixed
- **KE axes fixed at 0 to 30000 MeV** (`Plot3D.cpp`, `SecondaryCE.cpp`, `PlotDiffPhiCuts.cpp`). At 50 and 137.5 GeV the beam peak was off the axis. They now run to `KeAxisMax()` (1.2 * beam KE, rounded up to 5 GeV). The `_0to30k` and `_10to30k` file tags are built from that range, so they are unchanged at 25 GeV/c.
- **Theta axes fixed at 0 to 2.3 deg** (`PlotAngle.cpp`, `PlotDiffPhiCuts.cpp`, `Plot3D.cpp` thetaBlock). Theta at a given t goes as 1/p, so at 137.5 GeV the CE protons sat in the first few bins. The axes use `ThetaAxisMax()` (2.3 deg * 25000 / p).
- **`PlotPolarization.cpp` theta range and total-fit windows** (`kThetaMax` = 5 deg, 0.5 to 2 deg) scale by the same factor. The pdf names are built from the windows and stay `0p5to2deg` and `0to5deg` at 25 GeV/c.
- **`nChargedTrackerHitsEv` axis** was fixed at 0 to 40. It is now 5 * `MaxTrackerHits()`, which is the same 40 for 4 blocks.
- **Silent fallback to the 25 GeV/c beam** when `settings.csv` has no `BeamMomentum` now prints a warning.

### Not changed
- **t axes (0 to 0.05) stay.** t does not scale with p, so these axes are fine at every energy. An earlier draft of this audit said otherwise; that was wrong.
- **A_P is taken at the middle of a 0.05 deg theta bin** (`Plot3D.cpp`, `Plot3DCols.h`), and the start-theta scan starts at 0.15 deg (`KEScanP.cpp`). At 137.5 GeV the CE angles are about 5 times smaller, so a 0.05 deg bin is coarse and a 0.15 deg start may remove most CE protons. This is a physics choice, so I left it for you to decide.
- **EnergySD.** Each event uses the mean p for t and A_P, and the SD only enters as the error (`kBeamTFrac`, `kBeamAFrac`). That is fine for a small SD. A large SD also puts the fixed 80% KE cut inside the low tail of the beam.
- **Beam angle and vertex spread (`randPAngle`, `randPos`).** `7run.cc` does not write them to `settings.csv`, so the macros cannot tell a run had them. `ReconAngleAll.cpp` measures theta and phi from the z axis and `ImpactToAxis` assumes the beam at x = y = 0, while the true `LastTheta` is relative to the neutron direction. A beam angle of up to 5 mrad (0.29 deg) is as large as the CE angle at 50 to 137.5 GeV, so recon theta, t, A_P and phi are biased on such a run. The vertex sigmas (5.5 um, 3 um, 1.5 cm in z) are harmless. Fixing the angle needs the beam direction per event in the ntuple and a change to the sim.
- `nParticlesEv` (0 to 19) and `nZdcCharged` (0 to 9) stay fixed; they are fine for the current scans.
