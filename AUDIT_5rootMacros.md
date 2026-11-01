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
