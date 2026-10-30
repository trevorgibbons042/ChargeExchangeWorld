# What each step of the analysis does

Run with `bash 6jobs/SubmitEvents.sh <run> [cutdone]`. Each step is one PBS job (`6jobs/Plot.pbs`, STEP name in brackets) and waits for the jobs it needs.

## Analysis chain (SubmitEvents.sh)

1. **CutEvents** [cutevents] (`CutEvents.cpp`). Reads the sim `output*.root` files and writes one `cuts/cut_<name>.root` per cut: the CE events (`cut_CE`), the CE events in the tracker (`cut_CEtracker`, `cut_CEall`), the parent-track cuts, and all-event protons. Each file holds only the rows that pass.
2. **ReconAngleAll** [recon] (`ReconAngleAll.cpp`). Merges the hits of every proton, fits its track, and writes `cuts/cut_angles_all.root`. That file holds one tree, `AnglesRaw`, with `passQuality` (z order and residual cut), `inCE` (event has a carbon CE) and `finalKED` (detector KE). This is the slow fit, about 54 min at 100M events.
3. **ConvertAngles** [convert] (`ConvertAngles.cpp`), used instead of step 2 with `cutdone`. Turns an old 5-tree `cut_angles_all.root` into the one-tree format without refitting. Each thread converts one block of rows, and the blocks are merged in order. The old file is kept as `cut_angles_all_old5.root`.
4. **CalcAPandT** [calc] (`CalcAPandT.cpp`, run twice: plain and `"fDEnergy"`). Computes the analyzing power A and the theta of each CE proton, and writes `cuts/cut_APandT.root`.
5. **CountCE** [countce] (`CountCE.cpp`). Counts the CE events per block and cut. It only needs step 1.

## Plot chain (SubmitPlots.sh)

6. **CutTrees** [cuttrees] (`CutTrees.cpp`). Builds the KE and residual cut lists on `AnglesRaw`, so the plot macros select rows without copying them.
7. **ZdcColumns** [zdc] (`ZdcColumns.cpp`). Adds the Plot3D columns to the ZDC proton cut files and writes `cut_ZdcCE.root` and `cut_ZdcAllEv.root`. Only Plot3D reads them.
8. **Lanes**, all in parallel after step 6 (the cuts lane after step 7):
   - **cuts** [cuts]: `Plot3D.cpp` and `QualityStudy.cpp` write `1_position`, `2_energy` and `6_cuts`.
   - **fom** [fom]: `PlotAngle`, `PlotCEFraction`, `PlotPolarization` and `PlotFomPolarization` write `3_ReconCalcVals` and `5_fomCalc`.
   - **ceblock** [ceblock]: `PlotDiffPhiCuts`, `GraphBlocks` and `CETCompare` write `4_CEblock`.
   - **secce** [secce]: `SecondaryCE.cpp` writes `7_secondaryCE`.
9. **Pdfs** [pdfs]. Runs after every lane has ended. `tidy_plots.py` merges the per-cut csv files, `cut_side_box.py` stamps the cut numbers into a side box, `crop_pdfs.py` crops the pdfs, and `to_png.sh` makes the pngs.
