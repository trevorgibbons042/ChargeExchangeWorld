# Scan B: 5rootMacros (10 files)

Caveat: the Read tool returned function bodies elided ("caveman: body elided") for 9 of the 10 files. Only Plot3D.cpp was read in full. The other 9 were reconstructed from the visible signatures and comments, plus Grep of key lines. Line numbers below come from the on-disk files via Grep, not from the elided view. PlotCEFraction.cpp, CETCompare.cpp and GraphBlocks.cpp are only partly covered.

## Per file

### Plot3D.cpp (674 lines)
- Purpose: ZDC position, energy and KE-cut plots, the 6_cuts cut-folder study, ranking csv files, phi cos fits and P vs start theta.
- Single function `Plot3D()` at :62. Helper lambdas: `drawPlot` :144, `drawXY` :170, `fomStats` :341, `phiFit` :359.
- Reads (TFile::Open):
  - `cuts/cut_angles_all.root` at :77 (Angles, AnglesRaw) and again at :332 (AnglesRaw).
  - `cuts/cut_CE.root` :113, `cuts/cut_ParentID.root` :114, `cuts/cut_FakeID.root` :115, `cuts/cut_AllEvProtons.root` :127.
  - `TChain ceAll("ChargeExchange_values")` on `output0_t*.root`, rebuilt inside the cut loops at :398 and :411.
- Serial per-entry loops:
  - :84-89 Angles tree with SetBranchAddress/GetEntry, filling 3 std::set<pair<int,int>> (gRecoTracks, gRecoTracksD, gRes2Tracks).
  - :100-106 AnglesRaw, filling `std::map<pair,array<double,17>> gRawVals`.
  - Every Draw cut then calls the C++ functions `isReco`, `isRecoD`, `isRes2`, `rawOk`, `rawVal`, `thetaOk` (:41-58). Each one is a set or map lookup per entry, run through the TTreeFormula interpreter.
- Draw patterns:
  - `drawPlot` does Clear, `t->Draw(expr, cut, opt)`, GetHistogram, `FitCore` or `FitBlocks` in the tweak, `events2->GetEntries(cut)` for the CE percent box, WritePlotNumbers, SaveAs. Each plot makes about 2 tree passes (Draw plus GetEntries).
  - Overlay plots (:237, :294, :429, :549) do a second `Draw(..., "same")` on another tree.
  - 1_additionalinfo (:579-636): for each of 11 infoCuts (hitZDC, KECUT80, 8 sampleCuts) and each of 14 kInfoVars, two `Draw(..., "goff")` calls, one on tCE and one on tAll, plus a quantile sort of all of gRawVals (:582) when lo >= hi.
  - `fomStats` and `phiFit` use `SetEstimate(GetEntries()+1)`, then `Draw("rawVal(...)","goff")` and copy `GetV1()` into a C++ loop.
  - phiFit does 1 Draw for theta, then 1 Draw for the phi histogram, then a TF1 fit. The `fPhiCosCut` fit uses N0 fixed, A fixed, P fitted.
  - Per sampleCut x ke80: `GetEntries` x4, `fomStats` (2 Draws), `phiFit` x2 (4 Draws), a fresh TChain, then 1 drawPlot with 2 Draws.
  - phiSets loop (:478-543): 10 sets x 2 (j): phiFit, a Draw of theta:phi, then a 100-step P-vs-start-theta loop. Each step refills a TH1D over n entries and runs a TF1 "QN" fit, up to 100 steps, plus 2 graphs saved per set.
- Cuts (verbatim):
  - `boundaryType == 2 && pdg == 2212`, plus `&& isReco(fEvent, TrackID)` (position plots) or `isRecoD` (block).
  - `fabs(fY) < 2`, `fabs(fDy) < 2`.
  - KE: `fDEnergy < 0.8*beam` and `fDEnergy > 0.8*beam`, with `kBeamKE` from FitCore.h; pct loop is `{80}`.
  - sampleCuts (:313-321):
    - `isRes2(fEvent, TrackID)`
    - `rawVal(.,.,6) <= 1`
    - `rawVal(.,.,6) <= 2`
    - `rawVal(.,.,7) > 0.5`
    - `rawVal(.,.,8) <= 2`
    - `fabs(rawVal(.,.,12)) <= 5 && fabs(rawVal(.,.,13)) <= 5`
    - `rawVal(.,.,10) >= 3 && rawVal(.,.,10) <= 16`
    - `thetaOk(fEvent, TrackID)`
  - `isFromCE==1` (rawCE). `rawOk(fEvent, TrackID, k)` is added to nearly every Draw.
  - `kRes2Mm`, `PassesBlockResidual` come from QualityCut.h. The theta limit comes from BlockLimits.h (`FillBlockLimits(ReadSettings())`).
- Writes:
  - `plots/1_position/1_true/`, `plots/1_position/2_block/`, `plots/2_energy/`.
  - `plots/6_cuts/{1_hitZDC,2_KECUT80,3_ResCuts,4_nProtons,5_zOrder,6_resHalf,8_nHits,9_firstHit,10_thetacut}/` with `1_additionalinfo/` and `2_phi_cosfit/`.
  - `cut_ratios.csv` per cut folder and `block_theta_limits.csv`.
  - `rank_*.csv` x7 in `6_cuts/`: CE_recon_over_recon_proton_in_cut, CE_recon_over_CE_recon_before_cut, nonCE_recon_over_recon_proton_in_cut, FOM_theta, FOM_thetaBlock, P_phiRe, P_phiReD.
  - stats.csv and fit_values.csv via WritePlotNumbers, P_vs_start_theta.csv via WriteStatsRow.
  - Ends with `BoxAndCropPlots()` :672.
- Notes:
  - The `doTheta` block :651 is false, so the theta/phi angle plots are dead code.
  - `ReadSettings()` (settings.csv from disk) is called at :111, :340 (twice), :366, :490 and :510. Inside phiFit it is hoisted out of the event loop.
  - The TChain `ceAll` over `output0_t*.root` is rebuilt about 3 times per cut and per ke80, each time only to call `GetEntries()`.

### PlotAngle.cpp (386 lines)
- Purpose: CE and recon theta, phi, t and AP distributions with theta/KE cuts, FOM csv tables, and phi cos fits. Argument `energy` is "fKEnergy" or "fDEnergy" (suffix `_DE`).
- Functions: `DrawGroup` :51, `ThetaEdges` :120, `WriteFomCsv` :133, `PhiCosFit` :179, `PlotAngle` :226.
- Reads:
  - `cuts/cut_APandT.root` :233, trees CE, CalcValsRe, CalcValsReD, CalcValsReDIFF, CalcValsReDDIFF (+`_DE`). Friend trees are joined by entry index: `AddFriend(events2)` for 4 trees :243, `AddFriend(events)` :245.
  - `cuts/cut_angles_all.root` :335, trees Angles and AnglesCE.
  - `TChain rawCE("ChargeExchange_values")` on `output0_t*.root` :369.
- Cuts and cut mechanism: `ApplyCut(tree, ke.dir, block)` sets a TEntryList built by `CutTrees.cpp` (CutFiles.h `CutList` :37, `ApplyCut` :52). The KE cut set is `KeCuts()`, no cut plus 50 to 99% of `kBeamKE`, each with a `_res2mm` twin. Block trees (events3, events5) get `block=true`. A Draw cut string such as `(thetaRe<op>)` for the theta cuts. `CarbonCeCut()` for rawCE (FitCore.h).
- Loops and Draw/Fit:
  - `DrawGroup` is nested keDirs x thetas x specs. Per spec it does: `GetEntries("isFromCE == 1")` (:83), an optional apTree Draw (:88), `Draw` (:95), `FitCore` if s.fit, `FitCosPhi` if s.apTree, `GetEntries(cut && isFromCE==1)` (:103), then SaveAs. About 3 to 4 tree passes per plot. The DrawGroup call (:361) uses only `keAll.front()` (ke_all), with theta cuts on the 2_phi plots only.
  - The 37 specsCE plots plus AP-diff plots (specsCE.push_back loop :313-325).
  - FOM: `WriteFomCsv` x2 (:347, :348) with `thetaEdges = ThetaEdges(maxTheta, 0.001)`, so roughly maxTheta/0.001 bins. Per KE cut (only one, `keFom`): `Draw("theta<sfx>:ap<sfx>","","goff")` over the whole tree, a per-entry loop with `upper_bound` (:147-152), then the same for `all` (:155-160). This is the heaviest per-entry loop in the file and scales with all recon protons.
  - `PhiCosFit` x3 (:374, :378, :382): the weighted version Draws twice and the unweighted once, each followed by `h->Fit`.
  - `all->GetMaximum("outTheta")` and `GetMaximum("outDTheta")` :341, 2 tree scans.
  - The true-KE pass (`energy != "fDEnergy"`) returns right after the FOM tables, :359.
- Writes: `plots/3_ReconCalcVals/<ke>/...`, `plots/5_fomCalc/{1_reconTheta,2_reconBlockTheta}/1_csv/<ke>.csv`, stats.csv and fit_values.csv via WritePlotNumbers.
- Headers: FitCore.h, CutFiles.h. RDataFrame is included at :3 but not used.
- Notes: `ReadSettings()` is called at global scope :115. The time-stamp printf lines (TDatime) at :339, :344, :350, :384.

### PlotPolarization.cpp (~330 lines)
- Purpose: P from phi modulation in thin theta slices (kThetaStep 0.1 deg, to kThetaMax 5.0), for CE and all protons, true and block position, with AP from CE or recon. Argument `energy`.
- Functions: `Gather` :44 (static cache keyed by tree+expr+KE dir+block), `SliceOf` :60, `FitPhi` :64, `FitSlice` :76, `DrawSliceGrid` :105, `PlotPolarization` :137.
- Reads: `cuts/cut_APandT.root` :143, `cuts/cut_angles_all.root` :144. AP source tree via friend.
- Loops and Draw/Fit:
  - `Gather` does `SetEstimate(GetEntries()+1)`, `ApplyCut(t, keDir, block)` :50, one `Draw(expr,"","goff")` per unique key, then copies the columns into vectors.
  - Slice-width and AP-source loops reuse the cached columns (:194-195), so tree reads happen once per column set.
  - Per-slice fitting is in memory: `FitPhi` fills a TH1D (kPhiBins 100) and `h.Fit(&f, opt)` runs a TF1 "LQ" fit per slice. This is Minuit, serial, one fit per slice per pass.
  - Slice and pass counts are not visible; a `DrawSliceGrid` 3x3 page is written per pass.
  - `DrawCEPercent` :298 and `ce->GetEntries()` / `all->GetEntries()` :245.
- KE cut: only `KeCuts(keLabel).front()`, so ke_all only. Block set uses `kBlockCut` (QualityCut.h).
- Writes: `plots/5_fomCalc/<set.dir>/<ke>/` (slices_<tag>.csv, 1_P_vs_theta_<tag>.pdf, 2_pct_error_vs_theta_<tag>.pdf) and `plots/5_fomCalc/all_slice_fits.csv` (fopen "w", :168).

### PlotFomPolarization.cpp (~330 lines)
- Purpose: 2D maps of P_fit and % error over theta ranges [lo, hi), 0 to 5 deg in 0.05 deg steps (100 x 100 bins, about 5000 fits per map).
- Functions: `FitRange` :50 (direct Newton solve for the binned Poisson likelihood, replaces the Minuit fit that took about 80 s per map), `MakeMap` :84, `SaveMap` :92, `BinPhi` :165, `AddProton` :179, `PlotFomPolarization` :184.
- Reads: `cuts/cut_APandT.root` :196, `cuts/cut_angles_all.root` :198.
  - `all` tree: manual `SetBranchAddress` + `GetEntry` loop over all entries (:220-226), binning into per-slice phi histograms (PhiAP).
  - `ce` tree: `AddFriend` (:236), `ApplyCut(ce, ke.dir, set.block)` (:240), a single `Draw("thetaCE:phiCE:<ap>:apCE:theta:phi","","goff")` (:243), then per-entry `AddProton` loops (:247-251).
- Fits: the map fits use `FitRange` (no Minuit). Only `hd->Fit(&fFit,"Q")` at :306 is a TF1 fit, used for the wide fit plots.
- KE cut: ke_all only (`.front()`, :190). Block set uses `kBlockCut`.
- Writes: `plots/5_fomCalc/<1_reconTheta,2_reconBlockTheta>/<ke>/` with `5_P_fit_2d_<tag>.pdf/.csv` and `6_P_pct_error_2d_<tag>.pdf`. It writes the same 5_fomCalc folders as PlotPolarization.cpp.

### PlotDiffPhiCuts.cpp (~210 lines)
- Purpose: recon protons whose block phi differs from the CE phi by `[phiLo, phiHi]` (default -180 to -50): ZDC positions, KE, block residual, R^2 and impact for the window and for all recon protons.
- Function `PlotDiffPhiCuts` :25.
- Reads:
  - `cuts/cut_APandT.root` (CalcValsReDDIFF) :31, `cuts/cut_ParentID.root` :33, `cuts/cut_angles_all.root` :51, `cuts/cut_AllEvProtons.root` :80.
  - Draws to `goff` and loops in C++: :40 (`phiReD >= phiLo && phiReD < phiHi`), :55 (9-column Angles pull), :72 (`boundaryType == 2 && pdg == 2212`, 7 columns), :84 (same cut on AllEv), :176 (AnglesRaw columns).
  - Matching is by (fEvent, ParentID) or (fEvent, TrackID) in std::map lookups, in loops at :44, :62, :74, :85, :146, :163, :194.
- Cut: `ApplyCut(diff, "ke_all", true)` :39, the block residual cut.
- Writes: `plots/4_CEblock/3_diffphiBlock/`.

### PlotCEFraction.cpp (~150 lines)
- Purpose: cumulative CE percent and FOM-per-bin summary plots from csv, not from trees. No tree reads.
- Function `PlotCEFraction` :24. Reads `plots/5_fomCalc/<set>/ke_all/stats.csv` :50 (for each set) and the 1_csv ke csv files written by PlotAngle.cpp.
- Per theta set it draws TGraphs, binW 0.05, with `ReadSettings().at("TotalNeutrons")` :69 and `ReadSettings().at("APScale")` :88.
- Writes: `plots/5_fomCalc/<set>/ke_all/` (1_ce_percent_of_recon_per_bin, 2_nonce_percent_of_recon_per_bin, 3_ce_percent_of_neutrons_per_bin, 4_ce_percent_x_A2_per_bin, 5_integral_ce_x_A2_dtheta, 1_cumulative_ce_percent_log).
- Cheap, csv-only.

### SecondaryCE.cpp (~400 lines)
- Purpose: CE split into the beam neutron (primary) and secondary neutrons by `CECreator`: CEPolY, counts, KE, x, vertex z, phi fit, first tracker x, and event count by group.
- Function `SecondaryCE` :28.
- Reads: `TChain ce("ChargeExchange_values")` :30 on `output0_t*.root`, `TChain pi("ParticleInfo")` :146.
  - `ce.Draw(">>ceList", CarbonCeCut())` :38 makes a TEventList (CarbonCeCut is from FitCore.h).
  - `pi.Draw(">>zdcList","boundaryType == 2 && pdg == 2212")` :152. `pi.Draw(">>trkList","boundaryType == 0 && pdg == 2212")` :329.
  - Manual `SetBranchAddress` + `GetEntry(list->GetEntry(k))` loops over those lists at :52, :154, :195, :243, :295, :326, :331, :369.
  - So the CE chain is read several times (at least 5 passes over the list). Each pass uses the same list but different branch sets.
- Fits: TF1 "Q0" phi fit per group (4 groups) :260.
- Uses `ReadSettings().at("PPValue")` :127, :136.
- Writes: `plots/7_secondaryCE/`. Numbered pdfs include 2_CE_count_by_creator, 3_CEpolY_distribution, 14_first_tracker_x_by_creator, 15_..._cut80_by_creator.
- Headers: BlockLimits.h, FitCore.h, QualityCut.h.

### QualityStudy.cpp (~320 lines)
- Purpose: choose the thresholds in QualityCut.h as a common CE-percentile (default `keepCE = 0.95`) for residual and impact, true and block, plus R^2 and the other candidate variables.
- Functions: `Quantile` :27 (copies the vector by value), `QualityStudy` :27.
- Reads `cuts/cut_angles_all.root` :35 (AnglesRaw). The manual loops are `for (Long64_t e = 0; e < t->GetEntries(); e++){ t->GetEntry(e); ...}` at :46, :148, :204, :253, :268. That is about 5 serial passes over AnglesRaw, each with SetBranchAddress, and the passes at :253 and :268 re-scan the tree per variable.
- Fits: none. Quantile sorts of each variable vector.
- Writes: `plots/6_cuts/3_ResCuts/` (1_<name>.pdf to 4_, 5_r2 pdfs, thresholds.csv) and `plots/6_cuts/<folder>/` for the other candidate variables :215.
- Same output folders as Plot3D.cpp (6_cuts/3_ResCuts, and other 6_cuts subfolders).

### GraphBlocks.cpp (~160 lines)
- Purpose: CE per carbon block, z and xy position per block, and recon-only versions.
- Function `GraphBlocks` :16.
- Reads `cuts/cut_blocks.root` :23 (one tree per block, `blocks` vector), `cuts/cut_angles_all.root` :124 (AnglesCE).
- Draws: per-block z histogram :76, xy colz plots :104, grid :87.
- `GetEntries("reachedDetector==1")` :55 and `GetEntries("reachedNextTracker==1")` :57.
- Loop at :130: AnglesCE builds the `reconKeys` set. Then per block (:139-154): `Draw("fEvent:TrackID:reachedDetector","","goff")`, a loop over `GetSelectedRows` building a `TEntryList`, `SetEntryList`, `Draw` z histogram.
- Uses `ReadSettings().at("nCols_blocks")` and `nCols_detector` :91.
- Writes: `plots/4_CEblock/`, including `CE_per_block.pdf` and its stats.

### CETCompare.cpp (~60 lines)
- Purpose: t of the CE in carbon blocks vs the CE in the silicon trackers, plus a z-position plot.
- Function `CETCompare` :14.
- Reads `TChain ce("ChargeExchange_values")` on `output0_t*.root` :15 (the whole chain).
- Cuts: `carbon = CarbonCeCut()`, `tracker = "didChargeExchange>0 && !(" + carbon + ")"`.
- About 5 Draws "goff" (`fZ` x3 :28, `LastT` x2 :43) and 4 GetEntries (:29), so about 9 passes over the chain.
- Writes: `plots/4_CEblock/1_CarbonVsSi/` (z_CE_blocks_vs_trackers.pdf, t_carbon_vs_tracker.pdf, t_carbon_vs_tracker_norm.pdf).

## Cross-file

- Tree and file reuse:
  - `cuts/cut_angles_all.root` is opened by Plot3D (twice), PlotAngle, PlotPolarization, PlotFomPolarization, PlotDiffPhiCuts, QualityStudy and GraphBlocks.
  - `cuts/cut_APandT.root` is opened by PlotAngle, PlotPolarization, PlotFomPolarization and PlotDiffPhiCuts.
  - `cuts/cut_AllEvProtons.root` is opened by Plot3D and PlotDiffPhiCuts.
  - `cuts/cut_ParentID.root` is opened by Plot3D and PlotDiffPhiCuts.
  - `output0_t*.root` TChain on `ChargeExchange_values` is opened by Plot3D (about 3 per cut), PlotAngle, SecondaryCE and CETCompare. SecondaryCE also chains `ParticleInfo`.
- Repeated cuts:
  - `boundaryType == 2 && pdg == 2212` appears in Plot3D (as `zdcProtons`), PlotDiffPhiCuts (:72, :84) and SecondaryCE (:152).
  - `CarbonCeCut()` appears in PlotAngle, SecondaryCE and CETCompare.
  - The block residual cut (`kBlockCut` / `PassesBlockResidual`, `kRes2Mm`) appears in Plot3D, PlotAngle (via the block flag), PlotPolarization, PlotFomPolarization, PlotDiffPhiCuts and QualityStudy.
  - `fDEnergy > 0.8*beam` (80% KE) appears in Plot3D and in the KeCuts lists of PlotAngle.
  - `isFromCE==1` appears in Plot3D and PlotAngle.
- Repeated work:
  - The phi cos fit `N0 (1 + P*A*cos(phi))` is implemented separately in Plot3D (`phiFit`, `fPhiCosCut`, `fPstart`), PlotAngle (`PhiCosFit`, `fPhiCos`), PlotPolarization (`FitPhi`), PlotFomPolarization (`FitRange`), and SecondaryCE (phiFit[g]).
  - The theta-binned model A, `APFromT(TFromAngle(...))` at 0.05 deg bin centres, is computed per entry in Plot3D (:350, :368, :490) and PlotFomPolarization.
  - The FOM formula (sum of A^2 * 0.05 / N_neutrons) is in Plot3D `fomStats`, PlotCEFraction (`5_integral_ce_x_A2_dtheta`) and the PlotAngle csv.
  - Plot3D, PlotPolarization and PlotFomPolarization each pull the same CE and recon theta/phi columns from the same trees into vectors.
- Serial hot spots:
  - Plot3D: an interpreter-evaluated `rawVal`/`rawOk` per entry in nearly every Draw. Many Draws on the same events2/eventsAllEv/rawT trees (cuts x ke80 x variables). The 10 phiSets x 2 x up to 100 fit steps.
  - QualityStudy: 5 serial GetEntry loops. SecondaryCE: about 8 passes over lists.
  - PlotAngle `WriteFomCsv` at 0.001 deg bins.
  - Fit counts: PlotPolarization slice fits are Minuit TF1; PlotFomPolarization uses a direct Newton solver instead.
  - Everything is serial inside a single root process. No ROOT::EnableImplicitMT or RDataFrame use was seen in these 10 files, apart from the unused include in PlotAngle.
- Output folder sharing:
  - `plots/5_fomCalc/<set>/`: PlotAngle (1_csv), PlotPolarization, PlotFomPolarization and PlotCEFraction all write here. PlotCEFraction reads the stats.csv that the others write.
  - `plots/6_cuts/`: Plot3D and QualityStudy both write 3_ResCuts.
  - `plots/4_CEblock/`: GraphBlocks (root), CETCompare (1_CarbonVsSi), PlotDiffPhiCuts (3_diffphiBlock).
  - `plots/3_ReconCalcVals/`: PlotAngle only. `plots/7_secondaryCE/`: SecondaryCE only. `plots/1_position` and `plots/2_energy`: Plot3D only.
- Shared headers: FitCore.h (in all 10, ReadSettings, FitCore, WritePlotNumbers, DrawCEPercent, CarbonCeCut), CutFiles.h (PlotAngle, PlotPolarization, PlotFomPolarization, PlotDiffPhiCuts), QualityCut.h (Plot3D, SecondaryCE, QualityStudy), BlockLimits.h (Plot3D, SecondaryCE), FomCommon.h (PlotAngle? no, PlotFomPolarization, PlotCEFraction).
