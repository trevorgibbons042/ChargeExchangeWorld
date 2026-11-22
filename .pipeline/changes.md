# Changes: faster ROOT analysis (spec .pipeline/spec.md)

Nothing committed. Edits only in `5rootMacros/` and `6jobs/`. Run folders under `build/` untouched. All tests ran on scratch copies of
`build/root/Physics_xs1000_noFTF_labphi_4block_801_worldphi_100k` in the session scratchpad
(`.../scratchpad/T`, `T2`, reference `Tref` = the HEAD macros).

User decisions applied: OQ1 rule (same png, csv within 1e-7 relative, counts exact), CutTrees builds only ke_all unless `CUTTREES_ALL=1`,
no `rdfentry_` assumption (DrawPlan checks it at run time and falls back to `TTree::Draw`), `TEntryList::EnterRange` not used.

## Status by phase

| phase | status |
|---|---|
| P0 quick wins | done, tested |
| P1 cut cache + consumers | done, tested (Plot3D, PlotAngle, CETCompare, SecondaryCE, PlotDiffPhiCuts, CountCE, ReconAngleAll, CutEvents seeds) |
| P2 DrawPlan | `DrawPlan.h` + `TestDrawPlan.cpp` done and tested. Only CETCompare and PlotDiffPhiCuts were converted to it. NOT converted: Plot3D, PlotAngle, GraphBlocks, QualityStudy, SecondaryCE, PlotPolarization, PlotFomPolarization (they use the P1 caches and plain TTree::Draw / GetEntry, plus `EnableImplicitMT(NThreads())`) |
| P3 front-end | CountCE (C1 tree, z-window CopyTree) and ReconAngleAll CE map (C1 tree) done. NOT done: the parallel per-track loop in ReconAngleAll (`TThreadExecutor`), the CountCE `GetLeaf` hoist, `SetBranchStatus` in the remaining GetEntry loops (GraphBlocks :127-130, Plot3D, CountCE :206-209) |

## Files changed

### P0
- `5rootMacros/FitCore.h`: `NThreads()` (env `ANA_NTHREADS`, else hardware threads). `ReadSettings` memoised by path + mtime + size (returns a copy). `CarbonCeCut()` memoised (rebuilt when settings change). `BoxAndCropPlots()` returns at once when `PLOT_PIPELINE` is set. Added includes `<thread> <cstdlib> <sys/stat.h>`.
- `5rootMacros/CutTrees.cpp`: only the `ke_all` lists by default (`CUTTREES_ALL=1` builds all 8). "No cut" lists are built with `TEntryList::Enter(i)` for every entry, no scan (title stays "no cut"). Block lists still use `Draw(">>", kBlockCut, "entrylist")`.
- `5rootMacros/CutFiles.h`: `ListCut(dir, block)` returns the string CutTrees used for a list (DrawPlan `listCut`).
- `5rootMacros/crop_pdfs.py`, `cut_side_box.py`: per-pdf jobs run in `ProcessPoolExecutor(max_workers=os.cpu_count())`, same printed counts. `unbox_and_stamp` is a module-level function (picklable).
- `5rootMacros/plot.sh`: `export PLOT_PIPELINE=1`, `NCPU` (`NCPUS`, `nproc`, `sysctl`), `ANA_NTHREADS=NCPU` for CutTrees and BuildCaches, `NCPU/4` for the lanes. Pre-compile of all macros in parallel before the lanes (skipped with `ACLIC=`), then `run 'BuildCaches.cpp'` after CutTrees. Lane contents and the post-processing order are unchanged.
- `6jobs/Plot.pbs`: `export ANA_NTHREADS=${NCPUS:-16} PLOT_PIPELINE=1`, new `caches` step.
- `6jobs/SubmitPlots.sh`: `caches` job after `cuttrees`, the lanes depend on `caches`.
- `6jobs/Analysis.pbs` (new): 128 cpus, 12 h, same module and LD_LIBRARY_PATH as the other jobs, `cd $BASE/build/root/$RUN && zsh run.sh`, `NCPUS=128`.
- `ReconAngleAll.cpp`, `CalcAPandT.cpp`, `CountCE.cpp`, `CutEvents.cpp`: `EnableImplicitMT(128)` replaced by `NThreads()`.

### P1
- `5rootMacros/CutCache.h` (new): spec A. `SplitAtoms`, `AtomKey`, `GetCutTree`, `SeedCutCaches`, `ClearCutCache`, `CacheTree`, `ColumnAdder`. FNV1a64 signature, sidecar written last, `.lock` via `O_EXCL`, tmp + rename, closest-subset lookup, log line per call.
  Differences from the spec text (read these):
  - Candidate filter also requires that every opaque atom (`@...`) of the request is in the candidate. Without it C4 would pick the seed `cut_AllEvProtons` (2 atoms) over `cut_ParentID` (1 atom) and then die on the missing `@isParentCE`.
  - Atoms are wrapped in parentheses when joined into the Filter string (an atom such as `(a) || (b)` has its outer parentheses stripped, so joining without parentheses would break the precedence).
  - A lock whose holder is a dead process on the same host (`kill(pid,0)` ESRCH) is stale at once, besides the `CACHE_LOCK_STALE_S` age rule. A crashed build otherwise blocks everyone for 6 h (this happened in testing).
  - The stamp of a built cache also lists its base file, so a rebuilt base makes it stale.
- `5rootMacros/Plot3DCols.h` (new): the Plot3D maps and `isReco/isRecoD/isRes2/rawOk/rawVal/thetaOk` (now `inline`), `LoadPlot3DMaps()` (also fills the block limits), `Plot3DColumns()` (tag `p3d_v1`: `isReco isRecoD isRes2 thetaOk raw_0..16 ok_0..16`, one cached map lookup per row and thread), `P3()` regex rewrite, `Plot3DCE()` (C7), `Plot3DAllEv()` (C8), `Plot3DRaw()` (C9). No `unordered_map` change.
- `5rootMacros/BuildCaches.cpp` (new): builds C1 to C9 with all threads.
- `5rootMacros/CutEvents.cpp`: `SeedCutCaches()` after the merge loop.
- `5rootMacros/Plot3D.cpp`: maps moved to the header; events2 = C7, eventsAllEv = C8, rawT = C9; `P3()` on every expr and cut at every Draw/GetEntries; one `TChain ceAllChain` count (`nCEallZdcAll`) replaces the 17 chain rebuilds; `ROOT::EnableImplicitMT(NThreads())`.
- `5rootMacros/PlotAngle.cpp`: rawCE is the C1 tree with an empty cut.
- `5rootMacros/SecondaryCE.cpp`: `ce` = C1, `pi` = C4, new `pt` = C5 (first tracker protons); the TEventLists are now "all rows". Loop bodies unchanged.
- `5rootMacros/PlotDiffPhiCuts.cpp`: `pos` = C4, `allEv` = C6 (the seed, exact hit).
- `5rootMacros/CETCompare.cpp`: trees C1/C2/C3 with no selection.
- `5rootMacros/CountCE.cpp`: pass 1 on the C1 tree with no filter; `CopyTree` on the C1 tree with only the z window. Pass 2 stays on raw.
- `5rootMacros/ReconAngleAll.cpp`: the CE map reads the C1 tree.
- `5rootMacros/run.sh`: `rm -rf cuts/cache`; CutEvents (full cores); then CountCE next to the ReconAngleAll -> CalcAPandT -> CalcAPandT("fDEnergy") chain with `ANA_NTHREADS=NCPU/2` each; `wait`; plot.sh.
- Every plot macro that was not otherwise touched (GraphBlocks, QualityStudy, PlotPolarization, PlotFomPolarization, PlotCEFraction, and PlotAngle, SecondaryCE, Plot3D) starts with `ROOT::EnableImplicitMT(NThreads())`.

### P2
- `5rootMacros/DrawPlan.h` (new): `H`, `N`, `C`, `Run`, `Default`. One RDataFrame pass per tree under IMT; rows `(rdfentry_, values...)` collected with `Aggregate` into one flat vector per slot, sorted by entry, then filled on the main thread in entry order. Dedup by tree + effective selection (`(sel) && (listCut)`). The tree's TEntryList is removed for the RDF pass and restored (listCut stands for it). `rdfentry_` is checked once per tree against `TTree::GetEntry` on about 1000 sampled entries; if it fails, or booking or the event loop throws, the request runs as plain `TTree::Draw` (one warning). `DRAWPLAN_FALLBACK=1` forces the fallback. Not done: `DRAWPLAN_MEM_GB` batching, auto-binned varexp (needs explicit bins or a model), friend columns with a dot (untested, would fall back).
- `5rootMacros/TestDrawPlan.cpp` (new): T2.
- `5rootMacros/CETCompare.cpp`: the 5 fills go through DrawPlan (models `hAll hZC hZT hC hT`).
- `5rootMacros/PlotDiffPhiCuts.cpp`: the 5 column pulls are `C()` requests in one `Run()`; the diff request uses `ListCut("ke_all", true)` instead of `ApplyCut`.

## What the Tester should check

Results on the 100k scratch run (reference = HEAD macros, interpreted, `ACLIC=`):
- T4 (Plot3D columns): C7/C8/C9, every row, `isReco/isRecoD/isRes2/thetaOk` and all 17 `raw_k/ok_k` equal the old functions. 0 mismatches. `CECreator` is present in a built cache (C2).
- T1 (cache): C1 and C6 hit the seeds; C2 builds (+1 atom from raw); C3 builds from C2 (+1 atom), entries 212 = raw `TChain::GetEntries`; C4 builds from the cut_ParentID seed (+2 atoms), 10625 rows = the ZDC protons of cut_ParentID. Touching settings.csv made all stale and rebuilt. Two concurrent `root` processes asking for the same cache gave one build and one wait.
- T2 (`TestDrawPlan.cpp`): 1D, 1D cut, 2D, model (Sumw2), weighted, listCut, tree entry list set: class, bins, bin contents, entries, mean and stddev bit-equal (`%a`), attributes, title, Cols in order, N: all identical. The only reported "mismatch" is a test artefact (the reference model histogram is named `hM2`). Fallback mode matches except where the test deliberately put an entry list on the tree.
- T5 (pipeline, equivalent chain): full plot chain with the new macros vs the HEAD macros: same 487 files, 94 csv files compared, 389 png files compared byte for byte. One difference: `3_ReconCalcVals/4_ap/5_apBlockDiff.png` and its `fit_values.csv` (the 5th decimal of a fit constant, and the near-zero background parameter). The inputs are unchanged, so this is a Minuit fit difference under OQ1. Please look at it by eye.
- `run.sh` chain (CutEvents, CountCE, ReconAngleAll, CalcAPandT x2) vs the original cut files: all trees have the same entry counts and identical sorted value lists for every scalar branch (cut_blocks, cut_APandT, cut_angles_all, cut_CE, cut_ParentID, cut_AllEvProtons). Seeds written by CutEvents.
- Timings on this 100k run are too small to mean anything. T6 (100M on LCRC through `Analysis.pbs`) was not run.

Please also run on LCRC:
- `TestDrawPlan.cpp` there (OQ3): it prints whether `rdfentry_` is the entry number; with a warning in the log the requests fell back to `TTree::Draw`.
- `CutCache` concurrency with 5 lanes (the log lines `[CutCache] ... built|hit`). BuildCaches should be the only builder.
- plot.sh's pre-compile loop, which was not run (Mac has no working ACLiC here).

## Environment caveats (important)

- ROOT 6.38.04 exists locally. Every changed macro passes `root -l -b -q -e '.L File.cpp'` except `PlotAngle.cpp`, which throws at load because of its pre-existing file-scope `ReadSettings()` and has no settings.csv in `5rootMacros/`. It loads in a run folder.
- This Mac's cling is flaky ("Failed to materialize symbols ... [abi:nqe220106]"):
  - Interpreted mode (`ACLIC=`) now fails for `Plot3D.cpp` and BuildCaches, because the header code (RDataFrame Snapshot and the `Define` lambdas in `Plot3DCols.h`) is JIT-instantiated. The HEAD Plot3D ran interpreted. Small macros (CutCache users, DrawPlan users, CETCompare, PlotDiffPhiCuts, TestDrawPlan) run interpreted fine.
  - ACLiC (`+O`) needs `vdt` on the include path here (`gSystem->AddIncludePath("-I<spack vdt>/include")`), and then jitted strings (RDF `Filter("...")`, `Define` strings, string-literal macro arguments) fail at run time.
  - I therefore verified with a mixed driver: caches built interpreted, Plot3D/others compiled with ACLiC, CETCompare/PlotDiffPhiCuts interpreted. Compiled mode on Linux (LCRC) is what the spec targets and is the untested part.
  - If the user still runs `ACLIC=` on the Mac, Plot3D needs the caches prebuilt and probably a compiled run. A possible fix: build code only in BuildCaches.cpp behind `#ifdef CUTCACHE_BUILDER`, so plot macros do not instantiate RDF. Not done, since it changes the spec's API (ColumnAdder).
- The cut cache needs caches built by `BuildCaches.cpp` or lazily by the first macro. Run alone, a macro builds what it needs under the lock (works on Linux).

## Not verified
- 100M-event behaviour, memory, LCRC filesystem lock semantics, PBS jobs (only `bash -n`).
- `plot.sh` end to end (shell syntax checked with `zsh -n`; the lanes were run by hand in a driver).
