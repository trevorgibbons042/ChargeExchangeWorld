# The cuts in the analysis

## Event and track cuts (CutPlots.cpp (events), one `cuts/cut_<name>.root` each)

- **cut_CE**: events with a carbon charge exchange. `CarbonCeCut()` selects `ChargeExchange_values` rows with `didChargeExchange > 0`, `fZ` in front of the ZDC face (`zdcFrontZ`), and fZ inside one of the carbon blocks. The block windows come from `StackZ` and the run's `settings.csv` (nBlocks, CarbonThickness, firstBlockZ, zdcFrontZ, blockRegion), so they follow the geometry of the run.
- **cut_ParentID**: tracks whose parent is a CE track of the same event.
- **cut_AllEvProtons**: every proton that reaches the ZDC (`boundaryType == 2 && pdg == 2212`), for all events.
- **cut_ParentZdc, cut_ParentTracker, cut_CEtracker, cut_CEall**: the ZDC and first-tracker protons of the CE children, and the CE outside the carbon blocks.
- **Plot3D sets** (no cut files): `F1_2_Plot3D.cpp` opens `AnglesRaw` once per set with a row filter (`SetRowFilter`), evClass 0, 0 or 1, and any, and `P3()` (Plot3DCols.h) adds the Plot3D columns.

## Quality cut

None. `ReconAngleAll.cpp` keeps every track `SeedExtend` finds (`SeedExtend.h`) and fits. `impact` and `impactD` (closest approach of the fitted line to the beam axis, true and block position, mm) stay branches, and no cut uses them.

`inCE` marks the protons of an event in cut_CE. The view `Angles` is every row of `AnglesRaw`, and `AnglesCE` is `inCE == 1` (`AngleSelection`, CutFiles.h).

## Residual cut

None. No plot, entry list or column applies a residual cut. `resDRms` (RMS distance of the block hits to the fitted line, mm) stays a branch.

## KE cut (F9_PlotReconNew.cpp, the `2_KECUT80` and `6_megaCut` folders)

- The smeared total ZDC KE of the event (`zdcKEAllD`, ReconAngleAll.cpp: `zdcKEAll`, the KE of every particle at the ZDC entry, charged and neutral, smeared once with `/zdc/resStochastic` and `/zdc/resConstant`) above 80% of the beam KE. It is a whole-event cut: an event below it is dropped with all its particles. The folder titles name it.
