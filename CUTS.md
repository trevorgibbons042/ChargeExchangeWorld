# The cuts in the analysis

## Event and track cuts (CutEvents.cpp, one `cuts/cut_<name>.root` each)

- **cut_CE**: events with a carbon charge exchange. `CarbonCeCut()` selects `ChargeExchange_values` rows with `didChargeExchange > 0`, `fZ < 7000`, and fZ inside one of the four carbon blocks (2775-2975, 3525-3725, 4275-4475, 5025-5225).
- **cut_ParentID**: tracks whose parent is a CE track of the same event.
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

## KE cut (CutFiles.h, the list file `cuts/ke_gt80.root`)

- `ke_gt80`: `finalKE` above 80% of the beam KE (`finalKED` for the detector-energy plots). It is the only KE cut, and every plot uses it.

## One cut file per cut of 6_cuts (CutSamples.cpp)

- AnglesRaw cuts: `cut_ResCuts` (`passQuality` and `resDRms < 2`), `cut_nProtons1`, `cut_nProtons2`, `cut_zOrder`, `cut_resHalf`, `cut_nHits` (3 to 16), `cut_firstHit` (within 5 mm), `cut_thetacut` (below the ZDC corner limit).
- ZDC cuts: `cut_hitZDC_CE`, `cut_hitZDC_AllEv` (reconstructed track), `cut_KECUT80_CE`, `cut_KECUT80_AllEv` (also above 80% of the beam KE).
