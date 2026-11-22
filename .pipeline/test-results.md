# Test results (tester)

Scratch copy: <scratchpad>/V (copy of 100k run folder, new macros, cuts/cache wiped). Nothing under build/ touched.

| # | check | result |
|---|---|---|
| 1a | bash -n 6jobs/*.pbs, SubmitPlots.sh (Analysis, CutEvent, Plot, Run, SubmitPlots) | PASS |
| 1b | zsh -n 5rootMacros/plot.sh run.sh to_png.sh | PASS |
| 1c | python3 -I -m py_compile crop_pdfs.py cut_side_box.py | PASS |
| 2 | .L load, 15 of 17 cpp files (CalcAPandT CETCompare CountCE CutEvents CutTrees GraphBlocks PlotAngle(in run folder) PlotCEFraction PlotDiffPhiCuts PlotFomPolarization PlotPolarization QualityStudy ReconAngleAll SecondaryCE TestDrawPlan) | PASS |
| 2 | .L BuildCaches.cpp, Plot3D.cpp: "Failed to materialize symbols" | UNTESTABLE on Mac (known cling limit) |
| 3a | TestDrawPlan.cpp (ACLIC=, 28116-entry tree): 1D, 1D cut, 2D, model, weighted, listCut, tree entry list, N, Cols; mean/sd bit-equal | PASS (1 reported MISMATCH is the hM/hM2 test-artefact name, as the coder said) |
| 3b | Cache test (TestCutCache.cpp): C1 hit seed 36904; C2 built 37116 = raw; C3 built 212 = raw TChain count; C4 built 10625 = ParentID ZDC protons; C6 hit seed 65438 | PASS |
| 3c | Cache sidecars written (.cuts, with .root for built ones), no stale .lock left | PASS |
| 4 | Known diff, see below | PASS under OQ1 for 5 of 6 cells; 1 cell outside 1e-7 (see below) |
| 5 | git status: no change under build/, no commit (HEAD still 001c9f5) | PASS |

## 4. Known diff (T vs Tref fit_values.csv in 3_ReconCalcVals/4_ap/csv)
Same 24 rows both. Entries, chi2, ndf are identical. Differing value cells (new vs HEAD, relative):
- 5_apBlockDiff Constant 1224.5728 vs 1224.5727: 8.2e-8
- 5_apBlockDiff Mean2 -0.041055466 vs -0.041055467: 2.4e-8
- 5_apBlockDiff Sigma2 0.057505733 vs 0.057505736: 5.2e-8
- 5_apBlockDiff Background 1.193958e-06 vs 1.1893834e-06: 3.9e-3 rel, abs 4.6e-9 (parameter is about zero, error 0.386, so 1e-8 of its error)
- 5_apBlockDiff Background error 0.3861 vs 0.386: 2.6e-4 (print rounding)
- 6_apPct Background 5.0757008e-07 vs 5.0757024e-07: 3.2e-7 abs 1.6e-13 (near-zero parameter). NOT listed in the coder's notes: a second file differs.
5_apBlockDiff.png: 63221 vs 63126 bytes (0.15%). I did not view the png by eye (gs render not done); Minuit noise on a near-zero background, no physics change.
Strict reading of OQ1 (1e-7 relative on float cells): Background value cells fail it, but they are fits to zero. Needs your call.

## Informational
- CECreator leaf is absent in C4 (built from cut_ParentID source). The coder's note said it is present in C2 (built from raw); I did not test C2, so check whether any consumer needs CECreator from C4.
- Fresh runs print the ROOT 6.38 Snapshot compression-default Info line (zstd 505 vs 101): cache files are zstd now.

## UNTESTABLE here: run on LCRC
- BuildCaches.cpp and Plot3D.cpp end to end (RDF JIT fails in Mac cling), including Plot3D output vs HEAD.
- TestDrawPlan.cpp: look for the rdfentry_ fallback warning (OQ3).
- 5-lane CutCache concurrency: [CutCache] built|hit lines, BuildCaches the only builder.
- plot.sh precompile loop and full plot.sh, PBS jobs (Analysis.pbs, SubmitPlots.sh), 100M run timing.
