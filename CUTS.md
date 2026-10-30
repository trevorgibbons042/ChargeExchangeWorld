# The cuts in the analysis

## Event and track cuts (CutEvents.cpp, one `cuts/cut_<name>.root` each)

- **cut_CE**: events with a carbon charge exchange. `CarbonCeCut()` selects `ChargeExchange_values` rows with `didChargeExchange > 0`, `fZ < 7000`, and fZ inside one of the four carbon blocks (2775-2975, 3525-3725, 4275-4475, 5025-5225).
- **cut_ParentID**: tracks whose parent is a CE track of the same event.
- **cut_FakeID**: tracks in an event that has a CE track which is not their parent.
- **cut_AllEvProtons**: every proton that reaches the ZDC (`boundaryType == 2 && pdg == 2212`), for all events.
- **cut_ParentZdc, cut_ParentTracker, cut_CEtracker, cut_CEall**: the ZDC and first-tracker protons of the CE children, and the CE outside the carbon blocks.
- **cut_ZdcCE, cut_ZdcAllEv** (ZdcColumns.cpp): the ZDC protons above, with the Plot3D columns added.

## Quality cut (ReconAngleAll.cpp, the `passQuality` branch)

A proton passes if its z order is right (`zOrderOK`) and `PassesQuality` holds (QualityCut.h):
- impact (closest approach of the fitted line to the beam axis, true position) at most 16.2162 mm.
- impactD (same, block position) at most 16.5846 mm.
- The true-position residual cut is switched off (limit 1e9).

`inCE` marks the protons of an event in cut_CE. The view `Angles` is `passQuality == 1`, and `AnglesCE` is `passQuality == 1 && inCE == 1`.

## Block residual cut (CutTrees.cpp, the `_block` lists)

`resDRms <= 2.0` mm, the RMS distance of the block hits to the fitted line. It applies only to quantities from block positions. The env var `RES_CUT_MM` overrides it.

## KE cuts (CutFiles.h, one list file `cuts/<dir>.root` each)

- `ke_all`: no cut.
- `ke_gt50`, `ke_gt60`, `ke_gt70`, `ke_gt80`, `ke_gt90`, `ke_gt95`, `ke_gt99`: `finalKE` above that percentage of the beam KE (`finalKED` for the detector-energy plots).
- Each KE cut also comes with a twin folder `<ke>_res2mm` that adds `resDRms < 2`.
