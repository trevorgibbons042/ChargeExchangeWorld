# Scan A: ROOT chain front end (headers, CutEvents/ReconAngleAll/CalcAPandT/CountCE/CutTrees, shell, PBS)

All paths relative to /Users/trevorg04/G4ChargeExchange. Note: Read on the headers showed bodies elided except where I re-read by offset; header facts below come from signatures, comments and the ranges I re-read (CarbonCeCut, KeCuts, CutList).

## Shared headers (5rootMacros/)

### FitCore.h
- Purpose: fits, stats boxes, CE side boxes, csv writers, settings reader, shared beam constants and CE selection.
- Constants: kNeutronMass = 939.565, kBeamMomentum = 25000.0 (MeV/c), kBeamKE = sqrt(p^2+m^2)-m (24078 MeV).
- Key funcs: TFromAngle (l.41), APFromT (l.48), ReadSettings(path="settings.csv") (l.55, reads settings.csv, throws if missing), CarbonCeCut (l.73), PlotTitle (l.90), DrawFitTails (l.117), DrawCEPercent (l.140), DrawCELost (l.154), AddFitEquation (l.178), FitCosPhi (l.196, N0(1+AP P cos phi), AP fixed, P limited -1..1), FitBlocks (l.211, "I" option fit), FitCore (l.211), StatsTitle (l.310), WriteFitValues (l.317), WriteHistStats (l.335), WriteStatsRow (l.349), WritePlotNumbers (l.359), BoxAndCropPlots (l.366, calls cut_side_box.py / crop_pdfs.py).
- CarbonCeCut() calls ReadSettings() on every call (file read, no cache). It builds, from settings nBlocks, nBlocksTracking, dz, CarbonThickness: zStart = 4000 - (nElements-1)*dzMm/2, nElements = (nBlocks+1)*nTracking+nBlocks. Returns `"didChargeExchange>0 && fZ<7000 && (" + "(fZ>lo && fZ<hi) || ..." + ")"`, one window per carbon block (string cut, JIT by TTree::Draw/RDF Filter).
- Static initializer kPlainFonts sets gStyle fonts at include time.
- Used by: all macros in this group (plus plot macros outside the group).

### CutFiles.h
- Purpose: KE cut list and TEntryList helpers shared by plot macros. Includes FitCore.h, QualityCut.h.
- KeCuts(keLabel) (l.21): returns {"ke_all", ""} then for pct in {50, 60, 70, 80, 90, 95, 99}: dir "ke_gt<pct>", cut `"finalKE > " + (int)(kBeamKE*pct/100)`. Total 8 KE cuts. (Comment mentions _res2mm twins; the vector in this file holds only these 8.)
- kApTrees = {"CE","CalcValsRe","CalcValsReD","CalcValsReDIFF","CalcValsReDDIFF"} (cut_APandT.root); kAngleTrees = {"Angles","AnglesCE"} (cut_angles_all.root); kBlockTrees = {"CalcValsReD","CalcValsReDDIFF","Angles","AnglesCE"}.
- CutList(dir, t, block) (l.37): opens cuts/<dir>.root READ once (static map), gets TEntryList named <tree>[_block]. ApplyCut (l.52): t->SetEntryList(...).

### QualityCut.h
- ResidualCut(fallback) (l.15): env RES_CUT_MM overrides, "none" disables. kMaxResidualMm = 1e9, kMaxResidualDMm = ResidualCut(2.0), kMaxImpactMm = 16.2162, kMaxImpactDMm = 16.5846.
- kBlockCut = `"resDRms <= " + std::to_string(kMaxResidualDMm)`; kRes2Cut = "resDRms < 2"; kRes2Dir = "_res2mm"; kRes2Label = ", block residual < 2 mm".
- PassesQuality(res, impact, impactD) (l.31, impact cuts; true-position residual cut off); PassesBlockResidual(resD) (l.32).

### BlockLimits.h
- gBlockLimits, FillBlockLimits(settings) (l.14), BlockThetaOk(thetaBlockDeg, zVertexMm) (l.26). Per-carbon-block theta max (deg) so proton lands in ZDC. Used by Plot3D.cpp, SecondaryCE.cpp (outside this group).

### FomCommon.h
- Plot-side constants/helpers for 5_fomCalc: kSets {"1_reconTheta","2_reconBlockTheta"}, kKes {"ke_all"}, kColors, KeName, FomRoot (l.28), FomBlockOnly (l.30), kFitDir "fom_polarization/", ReadThetaMax (l.35, reads ke_all.csv), SaveLogCopy (l.44). No tree I/O.

## Macros (5rootMacros/)

### CutEvents.cpp, CutEvents() l.19
- Purpose: split raw ParticleInfo and ChargeExchange_values into CE-event-selected cut files.
- Reads: TChain "ParticleInfo" and "ChargeExchange_values" over `output0_t*.root` (chain level, all thread files).
- Writes: cuts/cut_CE.root, cuts/cut_ParentID.root, cuts/cut_FakeID.root (each = both trees merged via TFileMerger from cuts/tmp_<tree>_<name>.root), cuts/cut_AllEvProtons.root (`boundaryType == 2 && pdg == 2212` on ParticleInfo).
- Parallelism: ROOT::EnableImplicitMT(128) (l.20). 7 lazy Snapshots (RSnapshotOptions fLazy=true) run together via RunGraphs (l.81), then DisableImplicitMT (merge crashed with MT).
- Serial steps: (1) `ChargeExchange->Draw("fEvent:TrackID", CarbonCeCut().c_str(), "goff")` l.29 (TChain::Draw, single thread, string cut JIT; pulls all CE rows into memory), (2) loop l.35-40 builds unordered_map ceTracksInEvent and unordered_set ceEventParent, (3) TFileMerger merge loop l.84-92 (3 merges, serial), tmp file removal.
- Cuts: Filter lambdas inEvent(fEvent), isParent(fEvent, ParentID), isFake(fEvent, ParentID) on ParticleInfo (lookup tables, captured by reference, read-only in threads); dfCE = Filter(CarbonCeCut()) on ChargeExchange_values (snapshotted 3 times, once per cut file).
- Per-event: RDF event loop over ParticleInfo (one pass for 4 snapshots) and ChargeExchange_values (one pass for 3 snapshots). Snapshot copies all branches.
- Headers: FitCore.h.

### ReconAngleAll.cpp, ReconAngleAll() l.147
- Purpose: per-proton straight-line fit (true and block positions), quality variables; writes angle trees.
- Reads: ChargeExchange_values via RDF over `output0_t*.root` (l.151) for CE track map; then per raw file (TChain file list l.227-229, `ROOT::RDataFrame df("ParticleInfo", rawFile->GetTitle())` l.230) ParticleInfo branches fEvent, TrackID, pdg, fX, fY, fZ, fDx, fDy, fDz, fKEnergy, fDEnergy, ParentID, rdfentry_.
- Writes cuts/cut_angles_all.root: trees Angles, AnglesCE, Angles_DE, AnglesCE_DE, AnglesRaw (same ~37 branches; finalKE branch points to finalKED for the _DE trees).
- Functions: fitTrack (l.40), ResidualRms (l.109), RSquared (l.122), ImpactToAxis (l.141); struct Track (struct of arrays, l.22), FitResult l.18.
- Parallelism: EnableImplicitMT(128) l.148. Per-file RDF `Filter("pdg == 2212 && boundaryType != 2").ForeachSlot` l.234 pushes points into slotTracks[slot][{event,track,pdg}] (std::map per slot, one slot vector per raw file), then serial merge l.241-243 (Track::append).
- Serial hot spot: the whole per-track loop l.248-349 is single threaded (sort hits by entry, sort points by z, unique, fitTrack, ResidualRms x4 with Track::select copies, RSquared, std::map lookups ceTracks.find, 5 TTree::Fill into one TFile). Also serial: raw-file loop (files processed one after another; memory held per file only), nProtonsInEvent map l.245-246, ceTracks merge l.157-160 (map<int,set<int>>).
- Cuts applied: `pdg == 2212 && boundaryType != 2` (RDF filter), `t.size() < 2` skip, non-finite bx/by/bxD/byD skip, `if (!zOrderOK) continue;` and `if (!PassesQuality(resRms, impact, impactD)) continue;` before outTree/outTreeDE Fill; outTreeCE/outTreeCEDE filled when event is in ceTracks. AnglesRaw filled before quality cuts. isFromCE = ceTracks[event] contains ParentID.
- CE selection: CarbonCeCut() (same as CutEvents).
- Headers: QualityCut.h, FitCore.h.

### CalcAPandT.cpp, CalcAPandT(std::string energy = "fKEnergy") l.37
- Purpose: pair CE truth with recon, compute t and AP. Run twice by run.sh (plain, "fDEnergy").
- Reads: ChargeExchange_values from cuts/cut_CE.root (all branches used: fEvent, TrackID, LastPhiLab, LastTheta, LastMomentum, LastT, LastAP, fZ, rdfentry_); tree "Angles"+suffix from cuts/cut_angles_all.root filtered `isFromCE==1` (branches fEvent, ParentID, outTheta, outPhi, outDTheta, outDPhi, finalKE, PDG, rdfentry_, passResD, resDRms). Reads settings.csv APScale.
- Writes cuts/cut_APandT.root: trees CE, CalcValsRe, CalcValsReD, CalcValsReDIFF, CalcValsReDDIFF (+ "_DE" suffix on the second pass; file opened UPDATE for fDEnergy, RECREATE otherwise).
- Parallelism: EnableImplicitMT(128) l.40; two ForeachSlot passes l.180 and l.189 fill per-slot std::map<tuple,Track>; mergeSlots lambda l.164 (serial merge + std::stable_sort of points by entry).
- Serial hot spot: final join loop l.222-287 over tracksCE (map find in tracksAngles x2, 5 TTree::Fill per CE). All CE data and all CE-event Angles points held in memory as std::map.
- Cuts: `Filter("isFromCE==1")`, `if (pdg != 2212) {return;}`, skip if key not in both Angles maps. Uses last point of each track.
- Headers: FitCore.h (pBeam = kBeamMomentum, TFromAngle, APFromT).
- Same Angles tree is read once per pass (two passes = true KE then detector KE read the same input structure again).

### CountCE.cpp, CountCE() l.20
- Purpose: find carbon blocks from CE z positions, count CEs reaching next tracker and ZDC; writes block trees.
- Reads: TChains ChargeExchange_values and ParticleInfo over `output0_t*.root`. `ce->GetEntries()`, `info->GetEntries()` l.31 (entry count on chains).
- Writes cuts/cut_blocks.root: trees Block1..BlockN (CopyTree of CE rows + added branches nextTrackerZ, reachedNextTracker, trackerX/Y, trackerDX/DY, detectorDX/DY, detectorX/Y, reachedDetector).
- Cut: `cut = CarbonCeCut()` (third separate computation of the same string); CopyTree selection `Form("(%s) && fZ>=%.6f && fZ<=%.6f", cut.Data(), lo - 1.0, hi + 1.0)` l.172. gapMm = 300.0, binMm = 10.0, thickness filter keeps runs >= 0.5*longest.
- Parallelism: EnableImplicitMT(128) l.33. RDF pass 1 over ChargeExchange_values with `Filter(cut)` + ForeachSlot (slotZ, slotKeys) l.43; RDF pass 2 over ParticleInfo l.104 ForeachSlot into per-slot Local maps (trackerZ, protonHits, protonDetectorHits), merged serially l.113-117.
- Serial hot spots: CopyTree per block (l.175, single thread, string cut JIT, plus rewrites all CE branches); per-entry loop l.206-234 with `t->GetLeaf("fEvent")->GetValue()` and `GetLeaf("TrackID")->GetValue()` string lookups per entry, std::map lower_bound, 11 branch->Fill() calls per entry; std::sort of z, std::set/map merges.
- Whole pass 2 reads every ParticleInfo row (all pdg) of all files.
- Headers: FitCore.h.

### CutTrees.cpp, CutTrees() l.14
- Purpose: turn KE and block residual cuts into TEntryLists once so plot macros do not rerun Draw cuts.
- Reads: cuts/cut_APandT.root (trees kApTrees x {"", "_DE"}) and cuts/cut_angles_all.root (kAngleTrees x {"", "_DE"}) = 14 trees.
- Writes: cuts/<ke dir>.root for each of 8 KeCuts (ke_all, ke_gt50, 60, 70, 80, 90, 95, 99), one TEntryList per tree and "<tree>_block" for the 4 kBlockTrees (x2 suffix = 8 block lists). Per KE: 22 lists. 8 x 22 = 176 lists.
- Loop l.26-44: for ke in KeCuts, for tree in trees, for block in {false,true}: `t->Draw(">>name", cut, "entrylist")` l.37 (serial, string cut JIT, full-tree scan each; Angles/AnglesCE and CalcVals trees are the large ones). Cut = ke.cut, block: `ke.cut + " && " + kBlockCut` (or kBlockCut alone for ke_all).
- Fully serial macro; no MT, no RDF.
- Headers: CutFiles.h (FitCore.h, QualityCut.h).

## Shell and PBS (5rootMacros/, 6jobs/)

### run.sh
- zsh, `#source ~/.zshrc`; run() helper compiles with ACLiC `root -l -b -q "${1/.cpp/.cpp${ACLIC-+O}}"` and prints `>>> <macro> took N s`.
- Sequence (all serial): mkdir -p cuts; CutEvents.cpp; ReconAngleAll.cpp; CalcAPandT.cpp; CalcAPandT.cpp("fDEnergy"); CountCE.cpp; then `zsh plot.sh`.
- CountCE does not depend on ReconAngleAll/CalcAPandT outputs (reads raw output0_t*.root); CutEvents outputs (cut_CE.root) are read by CalcAPandT.

### plot.sh
- `source ~/.zshrc`; same run() helper. First `run 'CutTrees.cpp'` serial, then four parallel lanes (shell `&`), then `wait`:
  - lane 1: Plot3D.cpp, QualityStudy.cpp (1_position, 2_energy, 6_cuts)
  - lane 2: PlotAngle.cpp("fDEnergy"), PlotCEFraction.cpp, PlotPolarization.cpp("fDEnergy"), PlotFomPolarization.cpp("fDEnergy") (3_ReconCalcVals, 5_fomCalc)
  - lane 3: PlotDiffPhiCuts.cpp(-180, -50, "fDEnergy"), GraphBlocks.cpp, CETCompare.cpp (4_CEblock)
  - lane 4: SecondaryCE.cpp (7_secondaryCE)
- Then serial: python3 tidy_plots.py plots; python3 cut_side_box.py plots; python3 crop_pdfs.py plots; bash to_png.sh plots; prints `>>> all plots took N s`.
- Lane design comment: macros writing same plots/ folders run in order in one lane; no two lanes share a folder.

### 6jobs/Run.pbs
- PBS -A eic, select=1:ncpus=128, walltime=12:00:00; cd build/; `./sim run_0.mac` (sim only, no analysis).

### 6jobs/CutEvent.pbs
- PBS -A eic -N CutEvent, select=1:ncpus=128, walltime=12:00:00, mail bea. Loads root, sets LD_LIBRARY_PATH (xxhash, gif-lib), prints job info. Contains no command that runs CutEvents (file ends after echo lines).

### 6jobs/Plot.pbs
- select=1:ncpus=16, walltime=12:00:00. Takes env RUN, STEP. cd $BASE/build/root/$RUN (BASE default /lcrc/project/eic/tgibbons/G4ChargeExchange).
- STEP cases: cuttrees (CutTrees.cpp), cuts (Plot3D.cpp, QualityStudy.cpp), fom (PlotAngle.cpp("fDEnergy"), PlotCEFraction.cpp, PlotPolarization.cpp("fDEnergy"), PlotFomPolarization.cpp("fDEnergy")), ceblock (PlotDiffPhiCuts.cpp(-180, -50, "fDEnergy"), GraphBlocks.cpp, CETCompare.cpp), secce (SecondaryCE.cpp), pdfs (tidy_plots.py, cut_side_box.py, crop_pdfs.py, to_png.sh). Each macro `time root -l -b -q` with +O ACLiC.

### 6jobs/SubmitPlots.sh
- rsync 5rootMacros/ into build/root/$RUN/, then qsub Plot.pbs: CutTrees (no dep); Cuts, Fom, CEblock, SecCE each `-W depend=afterok:<CutTrees>`; Pdfs `afterok:cuts:fom:ce:sec`. Assumes cut_*.root and settings.csv already exist (CutEvents, ReconAngleAll, CalcAPandT, CountCE not submitted by it).

## Cross-file observations

Duplicated work
- CarbonCeCut() string (FitCore.h l.73, re-reads settings.csv per call) is evaluated in CutEvents (Draw l.29 + Filter l.56), ReconAngleAll (Filter l.153), CountCE (Filter l.43, CopyTree l.172). Each applies it as a JIT-compiled string cut over the full ChargeExchange_values tree of all output0_t*.root.
- The raw ChargeExchange_values tree (all `output0_t*.root`) is scanned with the CE cut separately in CutEvents (twice: Draw and RDF snapshot passes), ReconAngleAll, CountCE.
- The raw ParticleInfo tree is read in full by CutEvents (snapshot pass), ReconAngleAll (per file, pdg==2212 filter), CountCE (all rows, ForeachSlot).
- ReconAngleAll builds event -> CE track-ID map (std::map<int,set<int>>); CutEvents builds an equivalent (unordered_map + unordered_set of event/track) from the same selection.
- CalcAPandT runs twice per run (plain, "fDEnergy"); each pass re-reads cut_CE.root ChargeExchange_values and re-reads Angles(_DE) with the same joins; the CE-side join is identical in both passes.
- CutEvents' cuts/cut_ParentID.root, cut_FakeID.root, cut_AllEvProtons.root are written here; consumers are outside this group (plot macros). CalcAPandT reads only cut_CE.root.
- CutTrees applies the same 8 KE cuts to 14 trees: 176 Draw-with-entrylist scans of trees already filtered upstream.
- Per-track map/key (event, track, type) accumulate-then-merge pattern (slot vectors of std::map then serial merge) appears in ReconAngleAll, CalcAPandT (mergeSlots), CountCE.

Serial hot spots (per macro)
- CutEvents: TChain::Draw l.29, unordered_map build l.35-40, TFileMerger loop l.84-92.
- ReconAngleAll: per-file loop plus per-track fit/sort/Fill loop l.248-349, one raw file at a time.
- CalcAPandT: mergeSlots stable_sort, join/Fill loop l.222-287.
- CountCE: CopyTree per block, per-entry GetLeaf loop with 11 branch Fills l.206-234, serial merge of slot containers.
- CutTrees: all 176 entrylist Draws serial, single process.
- run.sh: the 5 macro invocations plus plot.sh run strictly in sequence.
- plot.sh: CutTrees first (serial), then 4 lanes, then 4 serial post-processing scripts (tidy_plots.py, cut_side_box.py, crop_pdfs.py, to_png.sh).

Parallelism setting
- ROOT::EnableImplicitMT(128) hard-coded in CutEvents l.20, ReconAngleAll l.148, CalcAPandT l.40, CountCE l.33. Plot.pbs requests 16 cpus; Run.pbs and CutEvent.pbs request 128.
- Each run folder loads whole file set via wildcard `output0_t*.root` (one file per sim thread); RDF parallelizes by tree cluster/file.

Memory
- ReconAngleAll: holds one raw file's proton points at a time (comment: all files at once does not fit at 100M events). CalcAPandT: holds all CE points and all CE-event Angles points (std::map<tuple,Track>) and all CE keys; CutEvents holds all CE rows from Draw (GetV1/GetV2); CountCE holds protonHits map for all CE protons and a vector<double> z of all CEs.

Output folders shared
- cuts/: written by CutEvents (cut_CE, cut_ParentID, cut_FakeID, cut_AllEvProtons, tmp_*), ReconAngleAll (cut_angles_all.root), CalcAPandT (cut_APandT.root; second pass UPDATE appends _DE trees), CountCE (cut_blocks.root), CutTrees (cuts/ke_*.root).
- plots/ lanes (from plot.sh, plot macros outside this group): lane1 1_position, 2_energy, 6_cuts; lane2 3_ReconCalcVals, 5_fomCalc; lane3 4_CEblock; lane4 7_secondaryCE. Post-processing scripts touch all of plots/.
- Run folder inputs: output0_t*.root, settings.csv (read by ReadSettings in FitCore.h: CarbonCeCut, CalcAPandT APScale, FillBlockLimits).

Job structure
- Cluster jobs exist for: sim (Run.pbs, 128 cpu), plot lanes (Plot.pbs, 16 cpu each, dependent chain via SubmitPlots.sh). No PBS job runs CutEvents, ReconAngleAll, CalcAPandT, CountCE (CutEvent.pbs has no run command). CutTrees job is first in the dependency chain.
