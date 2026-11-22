# Spec: faster ROOT analysis (100M events, 1x xs, 128 cores, LCRC PBS)

All paths are relative to /Users/trevorg04/G4ChargeExchange. Edit macros only in `5rootMacros/` and PBS files only in `6jobs/`. Line numbers come from the files on disk today.

## OPEN QUESTIONS (answer before Phase 2 lands; Phase 0/1 do not depend on them)

1. **What "identical" means.** Today's run is already not bit-reproducible. CutEvents writes cut_CE/ParentID/FakeID/AllEvProtons with MT `Snapshot` (CutEvents.cpp:61-81), so entry order changes from run to run, and with it the float sums in TH1 stats. Proposed rule: same set of files; the same rendered png (gs -r100) or a png diff that is only in the last printed digit; csv cells equal after parsing, with integer columns exact and float columns within 1e-7 relative. Is this OK? Default if there is no answer: yes.
2. **Timings at 100M.** Do per-step `>>> X took N s` numbers from a 100M run exist (run.sh and plot.sh print them)? They decide the lane rebalancing in Idea list items 6 and 8. Default: ship Phases 0-2, then measure.
3. **ROOT version on LCRC (`module load root`).** It decides whether RDF reports `rdfentry_` as the global TTree entry number under IMT (DrawPlan relies on it, with a fallback) and whether `TEntryList::EnterRange` exists. Default: the Coder checks at runtime with `TestDrawPlan.cpp` and takes the fallbacks.
4. **Default for CutTrees.** Should it build only the ke_all lists by default (the only ones any macro reads, see A.6) and build the 7 ke_gt* files only when `CUTTREES_ALL=1`? Default: yes.

---

## 0. Goal and hard rules

- Speed comes first. Plot outputs (pdf, stats.csv, fit_values.csv, the other csv files) stay the same under the rule in OQ1.
- The plot.sh order stays: CutTrees, then the plot macros, then tidy_plots.py, cut_side_box.py, crop_pdfs.py (last of the three), then to_png.sh.
- All graphics, fits and file writes stay in the main thread of each process. Only event loops are multithreaded (section D).
- Phases, each one shippable on its own:
  - P0: quick wins.
  - P1: cut cache plus consumers (requirements 2 and 3).
  - P2: DrawPlan, multithreaded plotting (requirement 1).
  - P3: the front-end macros.

---

## A. Cut cache (requirements 2 and 3), new header `5rootMacros/CutCache.h`

### A.1 Concepts
- **Atom.** One boolean term of a cut. A cut string is split at top-level `&&` into atoms: outside parentheses and outside quotes. One pair of outer parentheses around an atom is removed if it wraps the whole atom.
  - `CarbonCeCut()` (FitCore.h:73-86) gives 3 atoms: `didChargeExchange>0`, `fZ<7000`, and `(fZ>a && fZ<b) || ...`.
  - Each atom has two forms. The original text is what gets applied. The key form has all whitespace removed outside quotes and is used only for comparing.
- **Opaque atom.** Starts with `@`, for example `@isParentCE`. It names a selection that a string cannot express (the lambdas at CutEvents.cpp:43-52). It can be part of a cached set. The builder can never apply one, so a request that holds an opaque atom can only be met by a cache or seed that already contains it.
- **Columns tag.** An optional string such as `p3d_v1`. It means "the base rows plus these derived columns", added by a compiled callback. It is part of the signature. A cache with tag T can be used only for requests with tag T, and a cache with no tag can serve as the base for a request with a tag.
- **Source.** A file path or glob (`output0_t*.root`, `cuts/cut_ParentID.root`) plus a tree name. When resolving, a glob is expanded and sorted.

### A.2 Signature and file naming
- `sig = FNV1a64(sourceCanon + "\n" + tree + "\n" + tag + "\n" + join(sorted unique atom keys, "\n"))`, printed as 16 hex digits. Use FNV-1a, not std::hash, because it is the same across builds and hosts.
- `sourceCanon` is the sorted list of resolved source paths, relative to the run folder, joined with `;`.
- Layout, inside the run folder:
  ```
  cuts/cache/<tree>/<srcTag>/<sig>.root     tree named <tree>, all source branches + tag columns
  cuts/cache/<tree>/<srcTag>/<sig>.cuts     sidecar, written LAST (its presence = the cache is complete)
  cuts/cache/<tree>/<srcTag>/<sig>.lock     present only while a build runs
  ```
  `srcTag` is the sanitised basename of the source with no extension or wildcard: `output0_t`, `cut_ParentID`, `cut_angles_all`.
- Sidecar format, plain text, one `key=value` per line, then the atoms:
  ```
  file=cuts/cache/ParticleInfo/output0_t/<sig>.root   (seeds point at cuts/cut_*.root)
  source=output0_t*.root
  tree=ParticleInfo
  tag=
  entries=<N>
  stamp=<path>:<size>:<mtime>;...   every resolved source file + settings.csv + tag dependency files
  atom=<original atom text>          one line per atom, sorted by key form
  ```

### A.3 API (CutCache.h)
```cpp
#pragma once
// needs FitCore.h (CarbonCeCut); uses ROOT::RDataFrame, TSystem, <fstream>, <sys/file.h> or O_EXCL open
struct CacheTree {
    TFile* file = nullptr;          // owned by the cache module (static map), opened READ
    TTree* tree = nullptr;
    std::string path;               // .root file in use (cache, seed, or the source itself when nothing is cached)
    std::vector<std::string> atoms; // atoms (original text) that every row of tree already passes
    bool built = false;             // true if this call built it
};
using ColumnAdder = std::function<ROOT::RDF::RNode(ROOT::RDF::RNode)>; // adds the tag columns (Define with compiled lambdas)

std::vector<std::string> SplitAtoms(const std::string& cut);                 // top-level && split, outer parens stripped, "" and "1" dropped
std::string AtomKey(const std::string& atom);                                // whitespace removed outside quotes
// Rows of source/tree that pass every atom of `cuts` (+ opaque atoms), with tag columns if tag non-empty.
CacheTree GetCutTree(const std::string& source, const std::string& tree,
                     const std::vector<std::string>& cuts,                   // cut strings, each split by SplitAtoms; "@..." kept whole
                     const std::string& tag = "", const ColumnAdder& addCols = nullptr,
                     const std::vector<std::string>& tagDeps = {});          // files the tag columns read (in the stamp)
void SeedCutCaches();          // writes the seed sidecars for cuts/cut_*.root if missing or stale (A.5); GetCutTree calls it once
void ClearCutCache();          // rm -rf cuts/cache (run.sh)
unsigned NThreads();           // env ANA_NTHREADS, else std::thread::hardware_concurrency(); put in FitCore.h so all macros share it
```

### A.4 Lookup: closest subset (requirement 3)
1. `need` = the union of the atoms of every element of `cuts`, deduplicated by key.
2. Scan `cuts/cache/<tree>/*/` for `*.cuts` (all srcTags, because seeds sit under the raw source's dir, see A.5). Keep a candidate when all of these hold:
   - its `source` equals the request source, and its `tree` matches;
   - the stamp is valid: every listed file exists with the same size and mtime;
   - its atom keys ⊆ the `need` keys;
   - its tag is empty or equal to the request tag.

   A stale sidecar (stamp mismatch) gets its `.root` and `.cuts` deleted. If the `.lock` of that sig is held, it is skipped instead.
3. An exact match (same atom set and same tag) is opened and returned.
4. Otherwise pick the candidate with the most atoms. The user's rule is "largest cached subset". On a tie, take the smaller `entries`, then the lexicographically smaller sig. With no candidate, the base is the source itself and it has 0 atoms.
5. `missing = need − base atoms`. If `missing` holds an opaque atom, call `Fatal` with the atom name: no cache or seed provides it.
6. Build, as in A.6, from the base: Filter on `missing`, joined with `&&` in sorted key order, then `addCols` if the tag is set and the base has none. Write the new cache under the request's full `need` set and tag.
- Correctness holds because each cache holds exactly the rows of the source that pass its atoms, and atoms are pure functions of a row. So Filter(missing) on any subset cache gives the same row set as Filter(need) on the source. Row order follows the base (see ASSUMPTIONS).

### A.5 Seeds: CutEvents outputs registered as caches
`SeedCutCaches()` writes sidecars with `file=` pointing at the existing files. The source is `output0_t*.root`. The stamp covers the raw files, settings.csv and the seed file itself. A seed is valid only if it is newer than every raw file.

| file | tree | atoms |
|---|---|---|
| cuts/cut_CE.root | ChargeExchange_values | SplitAtoms(CarbonCeCut()) |
| cuts/cut_CE.root | ParticleInfo | `@inCEEvent` |
| cuts/cut_ParentID.root | ParticleInfo | `@isParentCE` |
| cuts/cut_FakeID.root | ParticleInfo | `@isFakeCE` |
| cuts/cut_AllEvProtons.root | ParticleInfo | `boundaryType == 2`, `pdg == 2212` |

CutEvents.cpp calls `SeedCutCaches()` after the merge loop (after line 92), so a fresh run is seeded at once. Old run folders get seeded on the first GetCutTree call.

### A.6 Build (writer)
- Lock: `open(<sig>.lock, O_CREAT|O_EXCL)`, then write `hostname pid time` into it. Lustre and GPFS keep O_EXCL atomic.
  - If the lock exists, poll every 2 s until the sidecar appears.
  - A lock older than `CACHE_LOCK_STALE_S` (default 21600) is removed and taken over.
- Write path:
  1. `ROOT::RDataFrame df(base tree or TChain of the source)`, then `.Filter(missingJoined)`, then `addCols`.
  2. `.Snapshot(tree, "<sig>.root.tmp.<host>.<pid>")` with all columns. Keep the IMT that is on (`NThreads()`).
  3. `rename` to `<sig>.root`.
  4. Write `<sig>.cuts.tmp`, then `rename` it.
  5. Unlink the lock.
- Snapshot of `CECreator` (char leaf, read at SecondaryCE.cpp:44-47) must survive. CutEvents.cpp:78 already snapshots that tree, so it is proven to work. The Coder still checks `CECreator` in a built cache, in test T3.
- Log one line per call: `[CutCache] <tree> <need atoms> <- base <path> (+<k> atoms) built|hit N entries in S s`.

### A.7 Which repeated cuts get a cache, and the consumers changed

| # | cut (repeat sites) | request | base it finds |
|---|---|---|---|
| C1 | CarbonCeCut on raw ChargeExchange_values: CutEvents.cpp:29,56; ReconAngleAll.cpp:153; CountCE.cpp:43,172; PlotAngle.cpp:371-374; SecondaryCE.cpp:38; CETCompare.cpp:16,28,29,43 | `GetCutTree("output0_t*.root","ChargeExchange_values",{CarbonCeCut()})` | seed cut_CE (exact hit) |
| C2 | `didChargeExchange>0` (CETCompare.cpp:28,29) | `{ "didChargeExchange>0" }` | none, so it builds from raw (one MT pass) |
| C3 | tracker CE `didChargeExchange>0 && !(carbon)` (CETCompare.cpp:16,28,29,43) | `{ "didChargeExchange>0", "!(" + carbon + ")" }` | C2 (closest subset, +1 atom) |
| C4 | ZDC protons of the CE children: Plot3D.cpp:176-179 and every events2 Draw; PlotDiffPhiCuts.cpp:72; SecondaryCE.cpp:152 (only rows with a groupOf hit are used, and they are exactly `@isParentCE`) | `{ "@isParentCE", "boundaryType == 2 && pdg == 2212" }` | seed cut_ParentID (+2 atoms) |
| C5 | first-tracker protons, SecondaryCE.cpp:329 (rows used only on a groupOf hit, so `@isParentCE`) | `{ "@isParentCE", "boundaryType == 0 && pdg == 2212" }` | seed cut_ParentID (+2) |
| C6 | AllEv with a redundant ZDC cut, PlotDiffPhiCuts.cpp:84 | `{ "boundaryType == 2 && pdg == 2212" }` | seed cut_AllEvProtons (exact) |
| C7 | Plot3D join columns on C4: isReco, isRecoD, isRes2, rawVal, rawOk, thetaOk (Plot3D.cpp:41-58, called per entry through TTreeFormula in about 200 Draws) | C4 + tag `p3d_v1` | C4 (+0 atoms, +tag) |
| C8 | the same on AllEv (Plot3D.cpp:191-264 etc.) | C6 + tag `p3d_v1` | C6 |
| C9 | the same on AnglesRaw (Plot3D.cpp:332-335, sampleCuts :313-321) | `GetCutTree("cuts/cut_angles_all.root","AnglesRaw",{},"p3d_v1",…)` | source (adds columns only) |

- **Not cached as row copies: `kBlockCut` and the KE TEntryLists (CutFiles.h, CutTrees.cpp).** PlotAngle joins trees by entry index through `AddFriend` (PlotAngle.cpp:241-245), and a copied subset would break that alignment. They stay TEntryLists; see P0.3 for the speedup.
- **Tag `p3d_v1` columns.** Built by `Plot3DColumns()` in Plot3D.cpp, which captures the existing `gRecoTracks`, `gRecoTracksD`, `gRes2Tracks` and `gRawVals`. They are read-only, so sharing them across threads is safe. Columns added:
  - `isReco`, `isRecoD`, `isRes2`, `thetaOk`: Int_t 0/1, the same functions as today.
  - `raw_0` … `raw_16`: Double_t, `gRawVals.at(key)[k]`, or NaN if the key is missing.
  - `ok_0` … `ok_16`: Int_t, the same logic as `rawOk`.

  tagDeps = `cuts/cut_angles_all.root`. Bump the tag string whenever the column logic changes.
- **Plot3D expression rewrite.** `static std::string P3(const std::string&)` uses std::regex and is applied to every expr and cut before Draw/DrawPlan:
  - `isReco(fEvent, TrackID)` becomes `(isReco==1)`; the same for `isRecoD`, `isRes2` and `thetaOk`.
  - `rawVal(fEvent, TrackID, K)` becomes `raw_K`.
  - `rawOk(fEvent, TrackID, K)` becomes `(ok_K==1)`.

  The call sites keep their strings. The global functions stay, unused, so the old behaviour can be checked against the new.

---

## B. Multithreaded plotting (requirement 1), new header `5rootMacros/DrawPlan.h`

### B.1 Approach chosen
- **In process: "MT extract, ordered serial fill".**
  - One RDataFrame event loop per input tree runs under `ROOT::EnableImplicitMT(NThreads())`. It evaluates selections and expressions and collects `(entry, values...)` with `Take`. All booked requests run together through `ROOT::RDF::RunGraphs`.
  - The main thread then sorts each request by entry and fills the histogram in entry order. This is the same fill order as TTree::Draw, so bin contents and the stats sums match it.
  - The costly parts run in parallel: I/O, decompression, the jitted cut/expression evaluation, and the removed TTreeFormula-to-cling calls. Filling stays serial, at about 10 ns per entry.
- **Why not RDF Histo1D with per-slot merge.** Slot assignment is not deterministic, so stats sums would differ in the last bits from run to run. TTree::Draw also makes TH1F/TH2F, and per-slot TH1F merges lose counts above 2^24 per bin in a different way.
- **Why not threads that each draw a plot.** TCanvas, gPad, gStyle, TF1 fits (TMinuit global, `TVirtualFitter`) and `gDirectory` lookups by name are not thread-safe.
- **Between processes:** the lanes in plot.sh, with a split thread budget (C.1).

### B.2 API
```cpp
#pragma once
// TTree::Draw-compatible lazy requests on a TTree (friends allowed), run in one MT pass per tree.
class DrawPlan {
public:
    struct Hist  {std::shared_ptr<struct HistReq> r; TH1* Get() const;};   // Get() only after Run(): a fresh Clone in gDirectory named like Draw would name it (replaces an object of that name, like Draw)
    struct Count {std::shared_ptr<struct CountReq> r; Long64_t Get() const;};
    struct Cols  {std::shared_ptr<struct ColsReq> r; Long64_t N() const; const double* V(int i) const;}; // V(0) = first expr, like GetV1()
    // varexp: "e1[:e2]>>name(nx,lo,hi[,ny,lo,hi])", or ">>name" filling a model the caller passes (TH1D with Sumw2, PlotAngle.cpp:192-194)
    // selection: TTree::Draw selection; if weighted=true it is used as the weight (entries with weight 0 dropped, like TSelectorDraw)
    // listCut: the string the TEntryList set on t stands for (KeCut.cut / kBlockCut / ""), and-ed into the selection (CutTrees makes the lists from exactly this string)
    Hist  H(TTree* t, const std::string& varexp, const std::string& selection = "", const std::string& listCut = "", TH1* model = nullptr, bool weighted = false);
    Count N(TTree* t, const std::string& selection, const std::string& listCut = "");          // == t->GetEntries(selection) with the list set
    Cols  C(TTree* t, const std::vector<std::string>& exprs, const std::string& selection = "", const std::string& listCut = ""); // == SetEstimate(all)+Draw(goff)+GetV*
    void Run();                          // book every pending request (dedup by tree+expr+selection+listCut), RunGraphs, then ordered fills
    static DrawPlan& Default();          // one per process
};
```
- **Dedup.** Identical requests share one result. For example `events2->GetEntries(cut)` is repeated 3 times per sample cut (Plot3D.cpp:408,419) and is computed once.
- **Expression translation**, TTreeFormula to RDF C++:
  - `"e1:e2"` gives x = e2 and y = e1, as TTree::Draw does.
  - Each expr becomes `Define("__dpN", "(double)(" + e + ")")`. Each selection becomes `Filter("(" + sel + ")!=0")`, or a weight column when `weighted`.
  - Friend columns must use the qualified name (`CE_DE.thetaCE`, as PlotAngle.cpp:260,315 do already).
  - `TMath::`, `fabs` and `fmod` work as they are.
  - Rule for the Coder: no int/int division in an expression. TTreeFormula works in double and C++ does not. A scan of today's expressions found none.
- **Fallback, always correct.** If jitting fails, or `TestDrawPlan` shows `rdfentry_` is not the TTree entry number (OQ3), that request runs as plain `t->Draw` on the main thread during `Run()`, with a one-time warning.
- **Histogram fidelity.** The result must match TTree::Draw in all of these:
  - class TH1F/TH2F, or the model's class;
  - name and bin edges;
  - line, fill and marker attributes copied from the TTree (TSelectorDraw copies `t->GetLineColor()` etc.). The Coder confirms the attribute list in T2;
  - `SetDirectory(gDirectory)`, and `fEntries` including under- and overflow.
- **Memory.** Requests on the same tree and selection share the `entry` vector. `DRAWPLAN_MEM_GB` (default 64) splits RunGraphs into batches when the estimated total of (rows × columns × 8 B) is above it. Estimate rows with `t->GetEntries()`.

### B.3 How a macro uses it: book, run, render
Each macro section that does Draw/GetEntries is rewritten in 3 steps:
1. Book every request into local handles. The loops and cut strings stay the same.
2. `plan.Run()`.
3. Run the render code in the original order, with `t->Draw(...)` / `GetHistogram()` / `GetEntries(cut)` / `GetV1()` replaced by `h.Get()` / `n.Get()` / `c.V(0)`, plus `h->Draw(opt)` for non-goff.

Steps 1 and 3 must loop over the same lists in the same order. Store the handles in a `std::vector` in step 1 and walk it in step 3. Writing csv and pdf files only in step 3 keeps today's append order.

### B.4 Per-macro DrawPlan scope
- **Plot3D.cpp.** All `drawPlot` and `drawXY` calls (:144-174), their overlay Draws (:239,296,431,551), the `GetEntries` calls (:157,159,179,288-289,335,401,408,419), `fomStats` (:341-355), `phiFit` (:359-395), the P-vs-start Draw (:487), and the 1_additionalinfo Draws (:602,615,617).
  - Book them in one pass over the same loops; `drawPlot` gets a `book` and a `render` overload.
  - The data-dependent parts stay in render on in-memory columns: aMean, the 100-step start-theta loop (:493-508), the fits, and the rank csv (:444-470).
  - The `rawCE`/`rawT` requests use the C9 tree. The events2 requests use C7 and the eventsAllEv requests use C8, all after `P3()`.
- **PlotAngle.cpp.** `DrawGroup` (:51-113):
  - the per-spec `GetEntries` (:83, :103, :105), the `apTree` Draw (:88) and the main Draw (:95), with listCut = `""` or `kBlockCut` for `s.block`;
  - `WriteFomCsv` columns (:144, :155);
  - `PhiCosFit` (:187, :194, :197), where rawCE becomes the C1 tree (cut already applied, so the selection is "");
  - the `GetMaximum` calls (:341), which become `C()` with max taken in memory, or RDF `Max` booked in the same plan.
- **CETCompare.cpp.** The 5 Draws and 3 GetEntries (:28-29, :43) go to the C1/C2/C3 trees, with no selection, plus `N()`.
- **GraphBlocks.cpp.**
  - The per-block `GetMinimum`/`GetMaximum`/`GetEntries(sel)`/Draw calls (:51-57, :76, :87, :104) are booked for all blocks, then rendered.
  - The reconKeys loop (:127-130) gets `SetBranchStatus("*",0)` plus the 2 branches it needs.
- **QualityStudy.cpp.** The 5 serial GetEntry loops (:46, :148, :204, :253, :268) become one `C()` request per distinct branch set over AnglesRaw, with all of them in one Run. The loops then iterate the column arrays in entry order. Quantile results stay identical because the order is the same.
- **SecondaryCE.cpp.**
  - The `ce` list passes (:52, :195, :242, :295) become one `C()` over the C1 tree with fEvent, TrackID, CEPolY, LastAP, LastPhiLab, fX, fY, CEVertexZ (if present).
  - `CECreator` is a char leaf, so it is read once with a TTreeReader loop on C1 (CECreator branch only) into a `std::vector<std::string>` in entry order.
  - The zdcList passes (:154, :222, :301, :326, :368) become one `C()` over the C4 tree. The trkList pass (:331) becomes `C()` over C5.
  - The loop bodies stay as they are, run on arrays.
- **PlotDiffPhiCuts.cpp.** The `goff` column pulls (:40, :55, :72, :84, :176) become `C()` with one Run. `:72` uses the C4 tree and `:84` uses C6.
- **PlotPolarization.cpp `Gather` (:44-58) and PlotFomPolarization.cpp `:243`.** These become `C()` with listCut; the `all` loop at PlotFomPolarization.cpp:220-226 becomes `C()`. Fits stay serial; they are cheap per the scan.
- **PlotCEFraction.cpp.** No change. It is csv-only and cheap.

---

## C. Processes, threads, PBS

### C.1 Thread budget
- `NThreads()` (FitCore.h) replaces the hard-coded `EnableImplicitMT(128)` in CutEvents.cpp:20, ReconAngleAll.cpp:148, CalcAPandT.cpp:40 and CountCE.cpp:33. Every plot macro calls `ROOT::EnableImplicitMT(NThreads())` at the start. On a laptop 128 threads would oversubscribe.
- `plot.sh`: `NCPU=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}`; each lane runs with `ANA_NTHREADS=$((NCPU/4))`. CutTrees and the cache prebuild run with the full NCPU.
- `PLOT_PIPELINE=1` is exported by plot.sh and Plot.pbs. `BoxAndCropPlots()` (FitCore.h:366-369) returns at once when it is set. Today Plot3D.cpp:672, QualityStudy.cpp:326 and PlotDiffPhiCuts.cpp:203 run cut_side_box.py and crop_pdfs.py over all of plots/ while other lanes are still writing. That races and wastes time, and plot.sh runs both scripts after `wait` anyway, so the final pdfs do not change. crop is idempotent (crop_pdfs.py:45-46) and side box unboxes first (cut_side_box.py:221).

### C.2 Cache prebuild step
New `5rootMacros/BuildCaches.cpp`, `void BuildCaches()`. It calls GetCutTree for C1-C9, so no lane builds a cache and lanes never wait on each other's locks.
- C7-C9 need the Plot3D maps. Move `LoadPlot3DMaps()` (Plot3D.cpp:76-109) and `Plot3DColumns()` into a new header `5rootMacros/Plot3DCols.h`, included by Plot3D.cpp and BuildCaches.cpp.
- plot.sh runs it after CutTrees, at full NCPU. Run alone, a macro still works, because it builds lazily under the lock.

### C.3 PBS
- New `6jobs/Analysis.pbs`, `select=1:ncpus=128`, walltime 12:00:00, the same module and LD_LIBRARY_PATH as Plot.pbs:14-17. It does `cd $BASE/build/root/$RUN && zsh run.sh`. At present no job runs CutEvents, ReconAngleAll, CalcAPandT or CountCE (CutEvent.pbs has no command).
- `6jobs/Plot.pbs`: add `export ANA_NTHREADS=${NCPUS:-16} PLOT_PIPELINE=1`. Add a `caches` STEP that runs `BuildCaches.cpp`.
- `6jobs/SubmitPlots.sh`: the `caches` step comes after `cuttrees`, and the lanes depend on `caches`.

---

## D. Thread and race safety rules
- **Main thread only:** TCanvas, gPad, gStyle (FitCore.h:28-31 static init included), TLegend, TLatex, SaveAs, every TF1 fit (Minuit and `TVirtualFitter` are global), WritePlotNumbers/WriteStatsRow/fopen, and `gDirectory->Get` by name.
- **Worker threads (RDF) only:**
  - jitted selection and Define expressions;
  - compiled Define lambdas that read const maps (C7-C9);
  - `Take`, `Count`, `Max`/`Min`.

  No histogram filling, no TFile writes. Snapshot is the exception, because RDF manages it.
- **TFile:** RDF opens its own per-slot handles. Main-thread TFiles stay read-only. Never write into a TFile while RDF reads it.
- **Cache files across processes:** O_EXCL lock, then tmp, then rename, then the sidecar last (A.6). Readers trust only the sidecar plus the stamp.
- **Lanes:** the existing "no two lanes share a plots/ folder" rule (plot.sh:14-15) stays. Nothing in this spec adds a writer to a folder.
- **CutEvents merge:** keep `DisableImplicitMT` before TFileMerger (CutEvents.cpp:82). The seed write comes after it.
- **ACLiC:** the lanes compile different macros. Each must include CutCache.h and DrawPlan.h through relative `#include "..."`. Pre-compile the headers' users once per run folder before the lanes (plot.sh: `for m in ...; root -l -b -q -e '.L m+O' &; wait`), so that two jobs never compile the same `.so` at once.

---

## E. File-by-file change list

| file | phase | change |
|---|---|---|
| `5rootMacros/FitCore.h` | P0 | `NThreads()`. `ReadSettings` memoised: a static cache keyed by path and mtime, returning a copy (it is called per plot, e.g. Plot3D.cpp:340,366,490,510, GraphBlocks.cpp:91). `CarbonCeCut()` memoised too. `BoxAndCropPlots()` honours `PLOT_PIPELINE`. |
| `5rootMacros/CutCache.h` | P1 | new, section A. |
| `5rootMacros/Plot3DCols.h` | P1 | new: `LoadPlot3DMaps()` (moved from Plot3D.cpp:76-109, with `unordered_map` + reserve allowed only if the lookups stay equal), `Plot3DColumns()`, `P3()`. |
| `5rootMacros/DrawPlan.h` | P2 | new, section B. |
| `5rootMacros/BuildCaches.cpp` | P1 | new, C.2. |
| `5rootMacros/TestDrawPlan.cpp` | P2 | new, test T2. |
| `5rootMacros/CutEvents.cpp` | P1 | `NThreads()`. `SeedCutCaches()` after :92. No other logic change. |
| `5rootMacros/CutTrees.cpp` | P0 | Default: only `KeCuts().front()` (ke_all), unless `CUTTREES_ALL=1`. For ke_all non-block lists (cut ""), build the list without a scan: `TEntryList l(name,"no cut",t); for i<N: l.Enter(i)`. The title stays "no cut". Block lists stay `Draw(">>",kBlockCut,"entrylist")`. Run the 14 trees' block Draws in parallel processes? No: keep it serial, since only 8 cheap block Draws remain. |
| `5rootMacros/CutFiles.h` | P0 | No API change. Add `static std::string ListCut(const std::string& dir, bool block)` that returns the string CutTrees used (DrawPlan listCut). |
| `5rootMacros/Plot3D.cpp` | P1+P2 | C7/C8/C9 trees replace events2/eventsAllEv/rawT. `P3()` on every expr and cut. A single `nCEallZdc` replaces the `TChain ceAll` rebuilt at :398 and :411 (17 opens of every raw file), computed once with `TChain::GetEntries()` before the loops. DrawPlan per B.4. |
| `5rootMacros/PlotAngle.cpp` | P1+P2 | rawCE (:369-374) uses the C1 tree with an empty cut. DrawPlan per B.4. |
| `5rootMacros/CETCompare.cpp` | P1+P2 | C1/C2/C3. The hists stay TH1D models (`hZAll` etc., :27, :42) passed as `model`. |
| `5rootMacros/SecondaryCE.cpp` | P1+P2 | C1/C4/C5 replace the raw TChains (:30-31, :146-147). The `GetBranch` checks (:32, :290) run on the C1 tree. Columns per B.4. |
| `5rootMacros/PlotDiffPhiCuts.cpp` | P1+P2 | C4 for `pos` (:72), C6 for `allEv` (:84), DrawPlan columns. |
| `5rootMacros/GraphBlocks.cpp`, `QualityStudy.cpp`, `PlotPolarization.cpp`, `PlotFomPolarization.cpp` | P2 | B.4, plus `EnableImplicitMT(NThreads())`. Every manual GetEntry loop that stays gets `SetBranchStatus("*",0)` plus the branches it needs. |
| `5rootMacros/CountCE.cpp` | P3 | Pass 1 (:40-46) uses the C1 tree with no filter. CopyTree (:172-175) uses the C1 tree with only the z window `fZ>=lo-1 && fZ<=hi+1`; the carbon atoms are already applied. Pass 2 stays on raw: it needs `boundaryType==1` of all particles (:107). |
| `5rootMacros/ReconAngleAll.cpp` | P3 | The CE map (:151-160) reads the C1 tree. The per-track loop (:248-349) becomes two steps. Step one is a parallel compute into `std::vector<Row>`, one per track, in map order, with `ROOT::TThreadExecutor::Foreach` over index chunks. `nProtonsInEvent[...]` (:328) becomes `.at()`, because the operator[] insert is not thread-safe. Step two is a serial in-order loop that copies the Row into the branch variables and runs the same 5 Fill calls under the same conditions (:343-348). This keeps the output entries identical. |
| `5rootMacros/run.sh` | P1 | `rm -rf cuts/cache` before CutEvents. After CutEvents, run `CountCE.cpp` in the background next to `ReconAngleAll`, then CalcAPandT ×2, then `wait`. Each gets `ANA_NTHREADS=$((NCPU/2))` while both run. |
| `5rootMacros/plot.sh` | P0+P1 | `export PLOT_PIPELINE=1`, NCPU and ANA_NTHREADS as in C.1, pre-compile, `run 'BuildCaches.cpp'` after CutTrees. The lanes and the post-processing order stay the same. |
| `5rootMacros/crop_pdfs.py`, `cut_side_box.py` | P0 | Collect the per-pdf jobs from the existing `os.walk` loops (crop_pdfs.py:60-61, cut_side_box.py:208-229). Run `crop(path)` / `unbox(path)+stamp(...)` in `concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count())`. Sum the return values for the same printed counts. Each pdf is processed by exactly one worker. |
| `6jobs/Analysis.pbs`, `Plot.pbs`, `SubmitPlots.sh` | P0/P1 | C.3. |

---

## F. Test plan
Use a reference run folder R at about 1M events, e.g. the newest under build/root/.
- T0, baseline: copy R to R_ref and run today's `zsh run.sh` there. Save `plots/`, all csv files and the `>>> took` times.
- T1, CutCache unit test, macro `TestCutCache.cpp` (scratch, not committed). In R:
  - `GetCutTree(C1)` must hit the seed;
  - C2, then C3, must log the base as C2 with +1 atom;
  - the entries of C3 must equal `TChain::GetEntries(trackerCut)` on raw;
  - deleting `settings.csv`'s mtime match (touch) must make C1-C3 stale and rebuilt;
  - two concurrent `root` processes asking for C3 must give one build and one wait.
- T2, `TestDrawPlan.cpp`. For about 10 Plot3D and PlotAngle requests (1D, 2D, weighted, friend, listCut), compare against `TTree::Draw` on the same tree:
  - class, name, nbins and edges, all bin contents including under/overflow, `GetEntries`;
  - `GetMean`/`GetStdDev` bit-equal (print `%a`);
  - line, fill and marker attributes;
  - `Cols` arrays element-wise in order.

  Also check that `rdfentry_` equals the TTree entry under IMT (OQ3).
- T3: the C4/C5 caches have `CECreator`/`fDEnergy` etc.; SecondaryCE run on cache equals the raw version (csv check per T5).
- T4, the Plot3D column equivalence: for C7/C8/C9, over every row, `isReco == isReco(fEvent,TrackID)` etc. for all 17 `raw_k`/`ok_k` (NaN equals NaN).
- T5, whole pipeline: in R, `rsync` the macros, then `zsh run.sh`. Then `diff -r --brief` the file lists of plots/ against R_ref, and run a python compare script (scratch) over every csv: exact for integers, rel 1e-7 for floats. Render every pdf with `gs -q -sDEVICE=png16m -r100` and pixel-compare; list the differing pdfs and review them by eye. A difference must be explained by OQ1, or it is a bug.
- T6, timing: put the `>>> took` lines of T0 and T5 side by side. Then do one 100M run on LCRC through `Analysis.pbs` and report the per-step times.
- Per CLAUDE.md, render and look at a sample: one pdf per plots/ subfolder, plus `6_cuts/rank_FOM_theta.csv`.

---

## G. More speedup ideas, ranked by expected gain (citations are file:line)

1. **Remove TTreeFormula-to-cling calls per entry in Plot3D.** `isReco`/`rawVal`/`rawOk`/`thetaOk` (Plot3D.cpp:41-58) are called from nearly every Draw cut (:178-180, :313-321, :345, :361, :485-487, :594-602) through the interpreter. They become plain columns (C7-C9). Probably the largest single plotting gain.
2. **No raw ParticleInfo scans in plotting.** SecondaryCE scans the raw chain twice (SecondaryCE.cpp:146-152, :329) and reads all branches on each of 6 GetEntry passes (no SetBranchStatus). It moves to the C4/C5 caches, built from cut_ParentID.
3. **CutTrees does 176 full-tree Draws but only the ke_all lists are read** (CutTrees.cpp:26-44 vs consumers PlotAngle.cpp:346,361, PlotPolarization.cpp:148, PlotFomPolarization.cpp:190, PlotDiffPhiCuts.cpp:39). The fix (P0) leaves 8 block Draws plus scan-free "no cut" lists.
4. **Parallel per-track loop in ReconAngleAll** (ReconAngleAll.cpp:248-349, serial today, 128 cores idle while it runs) (P3).
5. **Python post-processing in parallel.** crop_pdfs.py:60-61 and cut_side_box.py:208-229 call `gs` once per pdf, serially, over hundreds of pdfs. to_png.sh:6 already uses `xargs -P`.
6. **Earlier starts in run.sh and plot.sh.** CountCE does not need ReconAngleAll (run.sh:12-15). After P1, SecondaryCE and CETCompare need only CutEvents outputs, so they can start right after CutEvents, overlapping ReconAngleAll and CalcAPandT. Needs a small DAG in run.sh (`wait $pid` per dependency).
7. **Skip BoxAndCropPlots inside the lanes** (FitCore.h:366-369, called at Plot3D.cpp:672, QualityStudy.cpp:326, PlotDiffPhiCuts.cpp:203). Today it runs gs over the whole plots/ tree 3 extra times (P0).
8. **Split lanes further.** Lane 3 macros write disjoint subfolders: PlotDiffPhiCuts writes `3_diffphiBlock`, CETCompare `1_CarbonVsSi` with no csv, GraphBlocks the 4_CEblock root and `Block*/`. They could be 3 lanes. Lane 1's Plot3D and QualityStudy write different file names in 6_cuts. Verify that no stats.csv or fit_values.csv path is shared (grep the WritePlotNumbers/WriteStatsRow paths) before splitting.
9. **CalcAPandT runs twice** (run.sh:13-14) with the same CE-side join (CalcAPandT.cpp:164-287). Add a `"both"` mode that reads cut_CE once and writes the plain and `_DE` trees in one process, with identical tree content.
10. **Hoist repeated TChain opens and settings reads.** Plot3D builds `TChain ceAll` 17× (:398, :411). ReadSettings is read from disk per call (FitCore.h:55, used inside loops at Plot3D.cpp:510 and GraphBlocks.cpp:91). CarbonCeCut also re-reads settings.csv (FitCore.h:74) (P0 memo).
11. **Branch status in the remaining GetEntry loops:** Plot3D.cpp:84-106, GraphBlocks.cpp:127-130, CountCE.cpp:206-209. CountCE also does a `GetLeaf("...")` name lookup per entry; cache the TLeaf* once before the loop.
12. **Duplicate CE trees in CutEvents.** The same `dfCE` is snapshotted 3× (CutEvents.cpp:78) into cut_CE, cut_ParentID and cut_FakeID. Plot3D reads only `ParticleInfo` from ParentID and FakeID (Plot3D.cpp:120-124 get the CE trees but never draw them). Write it once. The Coder must confirm with grep that no other reader uses those copies.
13. **ACLiC pre-compile** once per run folder before the lanes (D). This removes repeated dependency checks and parallel compile races in Plot.pbs jobs that share a folder.
14. **PlotPolarization's Minuit slice fits** (PlotPolarization.cpp:64-76) stay serial for identical output. If they show up in timings, they could run in a forked process per energy set. A Newton solver like FitRange (PlotFomPolarization.cpp:50) would change the numbers, so it is not proposed.

---

## ASSUMPTIONS
- In ChargeExchange_values and ParticleInfo, every branch the requests use is `double` or `int`, as the ForeachSlot signatures show (ReconAngleAll.cpp:235-238, CountCE.cpp:104-105). So jitted C++ gives the same values as TTreeFormula (doubles). T2 checks this. If cling uses FMA contraction, last-bit differences fall under OQ1.
- In AnglesRaw, `(fEvent, TrackID)` is unique, so `rawVal` on AnglesRaw rows returns that row's own values. The C9 columns use the same map lookup anyway, so the semantics hold even with duplicates.
- In SecondaryCE, the ParticleInfo rows used (zdcList, trkList) are exactly those with `(fEvent, ParentID)` in the carbon-CE set. This holds because `groupOf` keys come from the CarbonCeCut rows (SecondaryCE.cpp:59), and `isParent` in CutEvents uses the same set (CutEvents.cpp:39,45). Last-wins map inserts (SecondaryCE.cpp:326) can depend on row order only when there are duplicate ZDC rows per proton. Those are accepted under OQ1, since the seed order is already non-deterministic.
- Row order of a cache equals the order of its base. Seeds come from MT Snapshot and are not in raw order, the same as today's inputs to Plot3D/PlotDiffPhiCuts.
- `TChain("ChargeExchange_values").GetEntries()` on raw (Plot3D.cpp:401,412) counts all rows, not carbon CE only, and stays that way.
- The LCRC shared filesystem supports atomic `rename` and `O_EXCL` create, as Lustre and GPFS do.
- crop_pdfs.py is idempotent (crop_pdfs.py:45-46), and cut_side_box.py unboxes before it stamps (cut_side_box.py:221). So dropping the in-macro BoxAndCropPlots calls leaves the final pdfs unchanged.
- At 100M events the Plot3D maps (`gRawVals` etc.) fit in node memory. The macro already holds them today.
- Only `_DE` (fDEnergy) plot passes run in plot.sh (plot.sh:21-27). The true-KE CutTrees lists are built only for compatibility.
