VERDICT: NEEDS WORK

Scope: spec.md, changes.md, test-results.md, `git diff`, new files CutCache.h, Plot3DCols.h, DrawPlan.h, BuildCaches.cpp, TestDrawPlan.cpp, 6jobs/Analysis.pbs.
The design is sound. The thread split is correct: all TH1 fills, fits, canvases and gStyle use stay on the main thread, and the MT lambdas only read globals filled before the loop.
Two defects must be fixed before a 100M run: a silent wrong result on the DrawPlan fallback path, and DrawPlan memory.

## HIGH (must fix)

1. DrawPlan fallback drops listCut, so PlotDiffPhiCuts gives wrong counts in exactly the case OQ3 worries about.
   - DrawPlan.h:242 (FallbackHist), :340 (counts), :349 (cols) all run `r->sel` only. The listCut is not added.
   - PlotDiffPhiCuts.cpp:44 passes `ListCut("ke_all", true)` and does NOT set the ke_all block list on `diff`.
   - If `rdfentry_` fails the check on LCRC, if booking or the loop throws, or with DRAWPLAN_FALLBACK=1, cDiff loses the `resDRms <= 2` cut. The phi window then holds unfiltered protons, and keys, hKE and every 4_CEblock window plot go wrong. Only a Warning is printed.
   - TestDrawPlan h6 ("listCut without a list set") is this exact case. The fallback "mismatch" that changes.md waves off is this bug, not a test artefact.
   - Fix: in every fallback, draw with `Effective(sel, listCut)`. The list and the cut are the same set, so this is also correct when the list is set. Then fix the header comment at DrawPlan.h:7.
2. DrawPlan memory at 100M (spec item DRAWPLAN_MEM_GB not done).
   - This 100k run has Angles = 150k rows and AnglesRaw = 213k rows. If the rows scale with events, 100M gives Angles about 150M rows.
   - PlotDiffPhiCuts puts cAng and cEx in one group (same tree, sel ""), so stride is 19 doubles, about 23 GB per copy.
   - Copies alive at once:
     - the per-slot vectors. DrawPlan.h:289 calls `slots[i].clear()`, which keeps the capacity.
     - slots[0] growing to hold every row (reallocations double it).
     - `order` (16 B per row).
     - `g->rows` (:317, kept until `plan` dies).
     - `r->v` (:356-359).
     - `angCols`/`exCols` in PlotDiffPhiCuts.cpp:71, :200.
   - Estimated peak is about 90 GB, against about 25-30 GB for the old SetEstimate/GetVal path. Plot3D (next item) runs in the next lane at the same time.
   - Fix:
     - Free the per-slot buffers in the merger (`std::vector<double>().swap(slots[i])`).
     - Release `g->rows` after the results are extracted (clear `flat` too).
     - Let PlotDiffPhiCuts read `V()` directly, without the assign copies.
     - Or implement DRAWPLAN_MEM_GB batching.
   - Measure with `/usr/bin/time -v` on LCRC.

## MEDIUM

3. Plot3D always calls LoadPlot3DMaps (Plot3D.cpp:57), even when C7/C8/C9 are cache hits.
   - At 100M that is a serial GetEntry loop over about 213M AnglesRaw rows into a `std::map<pair,array<double,17>>`, about 40 GB, plus 3 sets of about 150M nodes.
   - The maps are only needed to build the caches, plus `gRawVals.empty()` at :521. This is minutes and tens of GB with nothing gained, and it goes against "speed first".
   - Fix: load lazily inside the ColumnAdder or BuildCaches only. Keep a cheap flag in place of the `gRawVals.empty()` test.
4. Cache key misses code dependencies (Plot3DCols.h:88).
   - The p3d_v1 columns depend on kMaxResidualDMm/kRes2Mm (QualityCut.h) and on BlockThetaOk (BlockLimits.h).
   - The stamp only covers cut_angles_all.root and settings.csv. After the user edits the 2 mm cut and rsyncs it, isRecoD/isRes2 stay at the old values with a valid stamp, so counts are silently wrong.
   - Fix: add "QualityCut.h", "BlockLimits.h", "Plot3DCols.h" to kPlot3DTagDeps (rsync updates their mtime), or hash kBlockCut into the tag.
5. The rdfentry_ check is weak (DrawPlan.h:199-205).
   - It takes the first scalar leaf: TrackID in Angles, fEvent in ParticleInfo/ChargeExchange_values. Both repeat heavily, and neighbouring entries share fEvent.
   - A small entry offset under MT can pass the check and give a wrong fill order. The counts stay right but the stats sums drift.
   - Fix: prefer a Double_t leaf with continuous values (finalKE, fX, fZ), or compare 2-3 leaves.
6. Stale lock across hosts (CutCache.h:255-275, 428-437).
   - The dead-pid takeover works only on the same host. PBS jobs (caches job and lanes, or reruns) land on other nodes.
   - A caches job killed mid-build (walltime, OOM) leaves a lock that makes any later builder sleep up to 6 h with no message.
   - Fix: print a line while waiting (holder host/pid, age). Write the PBS_JOBID into the lock. Consider a lower CACHE_LOCK_STALE_S in Plot.pbs/Analysis.pbs.
7. PBS caches step gets ncpus=16 (Plot.pbs:5) and so builds C7-C9 (about 213M rows plus 38 columns) with 16 threads. It sits on the critical path ahead of all 4 lanes. Fix: give it its own select (e.g. 64-128 cpus) or a separate pbs.

## LOW

8. Cache reuse does not check the file set (CutCache.h:380): it matches on the glob string, and StampValid only checks files that were listed. If a new output0_t<N>.root is added and no old file changes, a stale cache is reused. Fix: compare Resolve(source) with the files in the stamp.
9. Row order of the caches is nondeterministic (MT Snapshot). Order-dependent code:
   - SecondaryCE.cpp:335 `zdcDE[...] = piDE` (last row wins). This only matters if a proton has more than 1 ZDC row.
   - The FirstHit tie on equal z.
   - Counts could then differ from HEAD. A check on LCRC is enough.
10. Missing checks:
    - PlotAngle.cpp:374-377 has no null check on `rawCE` before Draw (crash if the cache fails).
    - plot.sh:23 sends pre-compile errors to /dev/null.
    - The plot.sh:20 comment says the lanes "skip the dependency check". They do not, and no harm comes of it.
11. crop_pdfs.py and cut_side_box.py use `os.cpu_count()`. On a PBS node that is every core of the node, not NCPUS. Use `len(os.sched_getaffinity(0))` or NCPUS.
12. GetCutTree Enable/DisableImplicitMT when the macro had no IMT (CETCompare calls GetCutTree before EnableImplicitMT). This re-creates the TBB pool and is harmless. Better: enable IMT first in the macro.

## Checked, OK

- **CECreator:** in C1. cuts/cut_CE.root ChargeExchange_values has it, checked on the 100k run. C4 is ParticleInfo, which never had it. The only reader is SecondaryCE, from C1. Not a concern.
- **AddFriend joins:** PlotAngle.cpp:246-248 and PlotPolarization/PlotFomPolarization friend only trees of cut_APandT.root (the CalcAPandT loop, unchanged). No cache tree is friended. OK.
- **Cache subsets:**
  - Each seed's atoms match what CutEvents filters (carbon, @inCEEvent, @isParentCE, @isFakeCE, bt==2 && pdg==2212).
  - The subset rule and the rule that opaque atoms must be in the candidate are correct.
  - C3 builds from C2, and C4/C5 build from the ParentID seed.
  - SecondaryCE's switch to @isParentCE rows matches its groupOf keys (carbon CE).
- **Cache writes:** tmp files are unique per host and pid. The rename is atomic and the sidecar is written last. O_EXCL works on GPFS/Lustre.
- **Order of steps:**
  - plot.sh: tidy, side box, crop, to_png.
  - Plot.pbs: the pdfs step afterok all lanes, then caches afterok cuttrees.
  - PLOT_PIPELINE stops BoxAndCropPlots in the lanes.
  - run.sh runs CountCE in parallel with ReconAngleAll and CalcAPandT. There is no shared output, and both only hit C1.
- **Threads:** lanes 4 x NCPU/4 = NCPU, no oversubscription. DrawPlan fills on the main thread in entry order. Minuit runs in the main thread only.

## The 4_ap "Background" diffs: not a real concern

- FitCore.h:303/306 limits Background to [0, b0]. The fit puts it on the limit (about 1e-6, error 0.386).
- The change is 4.6e-9 absolute, about 1e-8 of its error. Constant/Mean2/Sigma2 move by less than 1e-7 relative. chi2, ndf and entries are the same.
- This is Minuit converging near a bound, not a physics change. OQ1 should allow fit parameters to match within 1e-6 of their error.
- The png differs only in the printed Background digits. Confirm by eye with gs.
- To rule out real nondeterminism, run PlotAngle("fDEnergy") with ANA_NTHREADS=1 and with 32 on the same inputs. If the two runs differ, the cause is IMT and not the caches.

## Run on LCRC first (in this order, on a 100k copy, then 100M)

1. `root -l -b -q 'TestDrawPlan.cpp+O'`, normal and with DRAWPLAN_FALLBACK=1. Note any rdfentry_ warning. After fix 1, h6 must match in fallback mode.
2. `rm -rf cuts/cache; /usr/bin/time -v root -l -b -q 'BuildCaches.cpp+O'`. Record wall time, peak RSS, cuts/cache size and the [CutCache] built/hit lines.
3. Full `zsh plot.sh` compiled, against HEAD macros on the same inputs. Plot3D output has never been compared in compiled mode. Check the pre-compile loop.
4. `/usr/bin/time -v` on Plot3D and PlotDiffPhiCuts at 100M: peak RSS against node memory, with all 4 lanes running.
5. Lock test: 2 processes asking for an uncached request (C5) on 2 nodes, then `kill -9` the builder and time the takeover.
6. Analysis.pbs end to end at 100M (T6), with the per-step `>>> took` times.
