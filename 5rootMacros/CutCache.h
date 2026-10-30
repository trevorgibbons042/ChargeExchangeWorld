#pragma once
// Cut cache (see .pipeline/spec.md section A): rows of a source tree that pass a set of cut atoms, built once per run folder and reused by every macro.
// Layout inside the run folder:
//   cuts/cache/<tree>/<srcTag>/<sig>.root   tree <tree>: every source branch + the tag columns
//   cuts/cache/<tree>/<srcTag>/<sig>.cuts   sidecar, written LAST (its presence = the cache is complete)
//   cuts/cache/<tree>/<srcTag>/<sig>.lock   only while a build runs
// The CutEvents outputs (cuts/cut_*.root) are registered as seed caches by SeedCutCaches(), the sidecar points at the existing file.
#include <TChain.h>
#include <TError.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>
#include <ROOT/RDFHelpers.hxx>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cerrno>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fcntl.h>
#include <fstream>
#include <functional>
#include <glob.h>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/time.h>
#include <thread>
#include <unistd.h>
#include <vector>

#include "FitCore.h" // CarbonCeCut, ReadSettings, NThreads

struct CacheTree {
    TFile* file = nullptr;          // owned by the cache module (static map), opened READ; null when the source is a glob (tree is then a TChain)
    TTree* tree = nullptr;
    std::string path;               // .root file in use (cache, seed, or the source itself when nothing is cached)
    std::vector<std::string> atoms; // atoms (original text) that every row of tree already passes
    bool built = false;             // true if this call built it
};
using ColumnAdder = std::function<ROOT::RDF::RNode(ROOT::RDF::RNode)>; // adds the tag columns (Define with compiled lambdas)

namespace CutCacheImpl {

// Directory helpers on TSystem: std::filesystem needs libstdc++fs on old GCC, which cling does not load (JIT "Failed to materialize symbols").
inline bool IsDir(const std::string& p){
    FileStat_t st;
    return gSystem->GetPathInfo(p.c_str(), st) == 0 && R_ISDIR(st.fMode);
}
inline std::vector<std::string> ListDir(const std::string& p){
    std::vector<std::string> o;
    void* d = gSystem->OpenDirectory(p.c_str());
    if (!d) return o;
    while (const char* e = gSystem->GetDirEntry(d)){
        const std::string n = e;
        if (n != "." && n != "..") o.push_back(p + "/" + n);
    }
    gSystem->FreeDirectory(d);
    return o;
}
inline void MakeDirs(const std::string& p){ gSystem->mkdir(p.c_str(), kTRUE); }
inline void RemoveAll(const std::string& p){
    if (IsDir(p)){
        for (const auto& c : ListDir(p)) RemoveAll(c);
        gSystem->Unlink(p.c_str());
    } else gSystem->Unlink(p.c_str());
}

inline std::string Trim(const std::string& s){
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) a++;
    while (b > a && std::isspace((unsigned char)s[b - 1])) b--;
    return s.substr(a, b - a);
}

// true if the '(' at index 0 closes at the last character (quotes respected)
inline bool OuterParensWrap(const std::string& a){
    if (a.size() < 2 || a.front() != '(' || a.back() != ')') return false;
    int depth = 0; char q = 0;
    for (size_t i = 0; i < a.size(); i++){
        const char c = a[i];
        if (q){ if (c == q && a[i - 1] != '\\') q = 0; continue; }
        if (c == '"' || c == '\''){ q = c; continue; }
        if (c == '(') depth++;
        else if (c == ')'){ depth--; if (depth == 0) return i == a.size() - 1; }
    }
    return false;
}

inline uint64_t Fnv1a64(const std::string& s){
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s){ h ^= c; h *= 1099511628211ULL; }
    return h;
}

inline std::string Hex16(uint64_t v){ char b[20]; std::snprintf(b, sizeof b, "%016llx", (unsigned long long)v); return b; }

inline bool IsOpaque(const std::string& atom){ return !atom.empty() && atom[0] == '@'; }

inline std::string HostName(){ char h[256] = "host"; ::gethostname(h, sizeof h - 1); h[sizeof h - 1] = 0; return h; }

inline long long FileSize(const std::string& p){ struct stat st; return ::stat(p.c_str(), &st) == 0 ? (long long)st.st_size : -1; }
inline long long FileMTime(const std::string& p){ struct stat st; return ::stat(p.c_str(), &st) == 0 ? (long long)st.st_mtime : -1; }
inline bool Exists(const std::string& p){ struct stat st; return ::stat(p.c_str(), &st) == 0; }

// sorted list of the files a source (path or glob) stands for
inline std::vector<std::string> Resolve(const std::string& source){
    std::vector<std::string> out;
    if (source.find_first_of("*?[") == std::string::npos){ if (Exists(source)) out.push_back(source); return out; }
    glob_t g{};
    if (::glob(source.c_str(), 0, nullptr, &g) == 0) for (size_t i = 0; i < g.gl_pathc; i++) out.push_back(g.gl_pathv[i]);
    ::globfree(&g);
    std::sort(out.begin(), out.end());
    return out;
}

inline std::string SrcTag(const std::string& source){
    std::string b = source.substr(source.find_last_of('/') == std::string::npos ? 0 : source.find_last_of('/') + 1);
    const size_t dot = b.find(".root");
    if (dot != std::string::npos) b.resize(dot);
    std::string o;
    for (char c : b) if (std::isalnum((unsigned char)c) || c == '_' || c == '-') o += c;
    return o.empty() ? "src" : o;
}

inline std::string CacheDir(const std::string& tree, const std::string& srcTag){ return "cuts/cache/" + tree + "/" + srcTag; }

inline std::string Join(const std::vector<std::string>& v, const std::string& sep){
    std::string o;
    for (size_t i = 0; i < v.size(); i++) o += (i ? sep : "") + v[i];
    return o;
}

inline std::string SourceCanon(const std::vector<std::string>& files){ return Join(files, ";"); }

inline std::string Signature(const std::string& canon, const std::string& tree, const std::string& tag, const std::vector<std::string>& sortedKeys){
    return Hex16(Fnv1a64(canon + "\n" + tree + "\n" + tag + "\n" + Join(sortedKeys, "\n")));
}

// "path:size:mtime" for every file that exists
inline std::vector<std::string> StampOf(const std::vector<std::string>& files){
    std::vector<std::string> o;
    for (const auto& f : files) if (Exists(f)) o.push_back(f + ":" + std::to_string(FileSize(f)) + ":" + std::to_string(FileMTime(f)));
    return o;
}

inline bool StampValid(const std::string& stamp){
    std::stringstream ss(stamp); std::string item;
    while (std::getline(ss, item, ';')){
        if (item.empty()) continue;
        const size_t b = item.rfind(':'); if (b == std::string::npos) return false;
        const size_t a = item.rfind(':', b - 1); if (a == std::string::npos) return false;
        const std::string path = item.substr(0, a);
        if (std::to_string(FileSize(path)) != item.substr(a + 1, b - a - 1) || std::to_string(FileMTime(path)) != item.substr(b + 1)) return false;
    }
    return true;
}

struct Sidecar {
    std::string sig, dir, cutsPath;
    std::map<std::string, std::string> kv;
    std::vector<std::string> atoms;      // original text
    std::vector<std::string> keys;       // AtomKey of each, sorted
    bool ok = false;
};

} // namespace CutCacheImpl

// Top-level && split, outer parentheses (one pair, only if they wrap the whole atom) stripped, "" and "1" dropped.
inline std::vector<std::string> SplitAtoms(const std::string& cut){
    using namespace CutCacheImpl;
    std::vector<std::string> out; std::string cur; int depth = 0; char q = 0;
    auto flush = [&](){
        std::string a = Trim(cur); cur.clear();
        if (OuterParensWrap(a)) a = Trim(a.substr(1, a.size() - 2));
        if (!a.empty() && a != "1") out.push_back(a);
    };
    for (size_t i = 0; i < cut.size(); i++){
        const char c = cut[i];
        if (q){ cur += c; if (c == q && cut[i - 1] != '\\') q = 0; continue; }
        if (c == '"' || c == '\''){ q = c; cur += c; continue; }
        if (c == '(') depth++; else if (c == ')') depth--;
        if (depth == 0 && c == '&' && i + 1 < cut.size() && cut[i + 1] == '&'){ flush(); i++; continue; }
        cur += c;
    }
    flush();
    return out;
}

// Whitespace removed outside quotes (comparison form only, the original text is what gets applied).
inline std::string AtomKey(const std::string& atom){
    std::string o; char q = 0;
    for (size_t i = 0; i < atom.size(); i++){
        const char c = atom[i];
        if (q){ o += c; if (c == q && atom[i - 1] != '\\') q = 0; continue; }
        if (c == '"' || c == '\''){ q = c; o += c; continue; }
        if (!std::isspace((unsigned char)c)) o += c;
    }
    return o;
}

namespace CutCacheImpl {

inline std::vector<std::string> Lines(const std::string& path){
    std::vector<std::string> o; std::ifstream in(path); std::string l;
    while (std::getline(in, l)) o.push_back(l);
    return o;
}

inline Sidecar ParseSidecar(const std::string& cutsPath){
    Sidecar sc; sc.cutsPath = cutsPath;
    const size_t slash = cutsPath.find_last_of('/');
    sc.dir = cutsPath.substr(0, slash);
    sc.sig = cutsPath.substr(slash + 1, cutsPath.size() - slash - 1 - 5); // minus ".cuts"
    for (const auto& l : Lines(cutsPath)){
        const size_t eq = l.find('=');
        if (eq == std::string::npos) continue;
        const std::string k = l.substr(0, eq), v = l.substr(eq + 1);
        if (k == "atom") sc.atoms.push_back(v); else sc.kv[k] = v;
    }
    std::sort(sc.atoms.begin(), sc.atoms.end(), [](const std::string& a, const std::string& b){ return AtomKey(a) < AtomKey(b); });
    for (const auto& a : sc.atoms) sc.keys.push_back(AtomKey(a));
    sc.ok = sc.kv.count("file") && sc.kv.count("source") && sc.kv.count("tree") && sc.kv.count("stamp") && sc.kv.count("entries");
    return sc;
}

inline void WriteSidecar(const std::string& cutsPath, const std::string& file, const std::string& source, const std::string& tree, const std::string& tag,
                         Long64_t entries, const std::vector<std::string>& stamp, const std::map<std::string, std::string>& atomsByKey){
    const std::string tmp = cutsPath + ".tmp." + HostName() + "." + std::to_string((long long)::getpid());
    {
        std::ofstream o(tmp);
        o << "file=" << file << "\nsource=" << source << "\ntree=" << tree << "\ntag=" << tag << "\nentries=" << entries << "\nstamp=" << Join(stamp, ";") << "\n";
        for (const auto& kv : atomsByKey) o << "atom=" << kv.second << "\n"; // map is sorted by key form
    }
    std::rename(tmp.c_str(), cutsPath.c_str());
}

// One TFile per path for the whole process, opened READ.
inline TFile* OpenRead(const std::string& path){
    static std::map<std::string, TFile*> files;
    auto it = files.find(path);
    if (it != files.end()) return it->second;
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (f && f->IsZombie()){ delete f; f = nullptr; }
    files[path] = f;
    return f;
}

inline std::vector<Sidecar> ScanSidecars(const std::string& tree){
    std::vector<Sidecar> o;
    const std::string root = "cuts/cache/" + tree;
    if (!IsDir(root)) return o;
    for (const auto& sub : ListDir(root)){
        if (!IsDir(sub)) continue;
        for (const auto& p : ListDir(sub)){
            if (p.size() > 5 && p.compare(p.size() - 5, 5, ".cuts") == 0){
                Sidecar sc = ParseSidecar(p);
                if (sc.ok) o.push_back(sc);
            }
        }
    }
    return o;
}

inline void Remove(const std::string& p){ std::remove(p.c_str()); }

inline bool LockHeld(const std::string& dir, const std::string& sig){ return Exists(dir + "/" + sig + ".lock"); }

// A lock whose file was not touched for this long (default 600 s) is stale: the builder beats every 30 s (LockBeat), so this also catches a dead job on another node.
inline long long LockStaleSeconds(){
    const char* e = std::getenv("CACHE_LOCK_STALE_S");
    return e ? std::atoll(e) : 600;
}

// Touches the lock file while a build holds it, from its own thread (a build can run long inside one RDataFrame call). Stops in the destructor, which every exit of the build passes.
struct LockBeat {
    std::string path; std::atomic<bool> stop{false}; bool skewWarned = false; std::thread th;
    explicit LockBeat(const std::string& p) : path(p), th([this]{
        const long long ticks = 2 * std::max(1LL, std::min(30LL, LockStaleSeconds() / 4)); // 0.5 s ticks between two touches
        long long n = 0;
        while (!stop){
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
            if (!stop && ++n >= ticks){
                ::utimes(path.c_str(), nullptr); n = 0;
                // the file system's time should be this node's time: a big gap (clock skew, attribute caching) makes the stale check of the waiting nodes wrong
                const long long gap = FileMTime(path) - (long long)std::time(nullptr);
                if (!skewWarned && (gap > 120 || gap < -120)){ skewWarned = true; Warning("CutCache", "lock %s mtime is %lld s off this node's clock: waiting builders may misjudge it as stale", path.c_str(), gap); }
            }
        }
    }) {}
    LockBeat(const LockBeat&) = delete; LockBeat& operator=(const LockBeat&) = delete;
    ~LockBeat(){ stop = true; th.join(); }
};

// Creates the lock with O_EXCL. A lock not touched for CACHE_LOCK_STALE_S (default 600 s) is removed first.
inline bool TryLock(const std::string& lockPath){
    const long long staleS = LockStaleSeconds();
    bool stale = Exists(lockPath) && (long long)std::time(nullptr) - FileMTime(lockPath) > staleS;
    if (Exists(lockPath) && !stale){ // a lock of a dead process on this host (a crashed build) is stale at once
        const auto info = Lines(lockPath);
        std::istringstream is(info.empty() ? "" : info[0]);
        std::string host; long long pid = 0;
        if (is >> host >> pid && host == HostName() && pid > 0 && ::kill((pid_t)pid, 0) != 0 && errno == ESRCH) stale = true;
    }
    if (stale){
        Warning("CutCache", "removing stale lock %s", lockPath.c_str());
        Remove(lockPath);
    }
    const int fd = ::open(lockPath.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0644);
    if (fd < 0) return false;
    const char* job = std::getenv("PBS_JOBID"); // who holds it, for the waiting log line of another builder
    const std::string info = HostName() + " " + std::to_string((long long)::getpid()) + " " + std::to_string((long long)std::time(nullptr)) + " " + (job ? job : "-") + "\n";
    if (::write(fd, info.data(), info.size()) < 0) {}
    ::close(fd);
    return true;
}

inline void DeleteStale(const Sidecar& sc){
    const std::string file = sc.kv.at("file");
    if (file.compare(0, 11, "cuts/cache/") == 0) Remove(file); // never delete a seed (cuts/cut_*.root)
    Remove(sc.cutsPath);
}

inline CacheTree OpenAsCacheTree(const std::string& path, const std::string& source, const std::string& tree, const std::vector<std::string>& atoms){
    CacheTree ct; ct.path = path; ct.atoms = atoms;
    if (path == source && source.find_first_of("*?[") != std::string::npos){ // the raw source itself: a chain
        TChain* c = new TChain(tree.c_str());
        c->Add(source.c_str());
        ct.tree = c;
        return ct;
    }
    ct.file = OpenRead(path);
    if (ct.file) ct.tree = (TTree*)ct.file->Get(tree.c_str());
    return ct;
}

constexpr const char* kRawSource = "output0_t*.root";

} // namespace CutCacheImpl

// Seed sidecars for the CutEvents outputs (cuts/cut_*.root), written if missing or stale. A seed is valid only if its file is not older than every raw file.
inline void SeedCutCaches(){
    using namespace CutCacheImpl;
    const auto raw = Resolve(kRawSource);
    if (raw.empty()) return;
    long long rawMax = 0;
    for (const auto& r : raw) rawMax = std::max(rawMax, FileMTime(r));
    std::vector<std::string> carbon;
    try { carbon = SplitAtoms(CarbonCeCut()); } catch (const std::exception&){ return; } // no settings.csv: no seeds
    struct Seed {std::string file, tree; std::vector<std::string> atoms;};
    const std::vector<Seed> seeds = {
        {"cuts/cut_CE.root", "ChargeExchange_values", carbon},
        {"cuts/cut_CE.root", "ParticleInfo", {"@inCEEvent"}},
        {"cuts/cut_ParentID.root", "ParticleInfo", {"@isParentCE"}},
        {"cuts/cut_FakeID.root", "ParticleInfo", {"@isFakeCE"}},
        {"cuts/cut_AllEvProtons.root", "ParticleInfo", {"boundaryType == 2", "pdg == 2212"}},
    };
    const std::string canon = SourceCanon(raw);
    for (const auto& s : seeds){
        if (!Exists(s.file) || FileMTime(s.file) < rawMax) continue;
        std::map<std::string, std::string> byKey; for (const auto& a : s.atoms) byKey[AtomKey(a)] = a;
        std::vector<std::string> keys; for (const auto& kv : byKey) keys.push_back(kv.first);
        const std::string sig = Signature(canon, s.tree, "", keys);
        const std::string dir = CacheDir(s.tree, SrcTag(kRawSource));
        const std::string cuts = dir + "/" + sig + ".cuts";
        if (Exists(cuts)){
            const Sidecar sc = ParseSidecar(cuts);
            if (sc.ok && StampValid(sc.kv.at("stamp"))) continue;
        }
        TFile* f = OpenRead(s.file);
        TTree* t = f ? (TTree*)f->Get(s.tree.c_str()) : nullptr;
        if (!t) continue;
        MakeDirs(dir);
        std::vector<std::string> files = raw; files.push_back("settings.csv"); files.push_back(s.file);
        WriteSidecar(cuts, s.file, kRawSource, s.tree, "", t->GetEntries(), StampOf(files), byKey);
    }
}

inline void ClearCutCache(){
    CutCacheImpl::RemoveAll("cuts/cache");
}

// Rows of source/tree that pass every atom of `cuts` (+ opaque atoms), with the tag columns if tag is non-empty. See spec A.4 for the lookup rules.
inline CacheTree GetCutTree(const std::string& source, const std::string& tree, const std::vector<std::string>& cuts,
                            const std::string& tag = "", const ColumnAdder& addCols = nullptr, const std::vector<std::string>& tagDeps = {}){
    using namespace CutCacheImpl;
    static bool seeded = false;
    if (!seeded){ seeded = true; SeedCutCaches(); }
    const auto t0 = std::chrono::steady_clock::now();

    const std::vector<std::string> files = Resolve(source);
    if (files.empty()){ Error("GetCutTree", "no file matches %s", source.c_str()); return CacheTree(); }
    const std::string canon = SourceCanon(files);

    std::map<std::string, std::string> need; // key -> original text
    for (const auto& c : cuts) for (const auto& a : SplitAtoms(c)) need.emplace(AtomKey(a), a);
    std::vector<std::string> needKeys, needAtoms;
    for (const auto& kv : need){ needKeys.push_back(kv.first); needAtoms.push_back(kv.second); }
    const std::set<std::string> needSet(needKeys.begin(), needKeys.end());
    const std::string sig = Signature(canon, tree, tag, needKeys);
    const std::string srcTag = SrcTag(source);
    const std::string dir = CacheDir(tree, srcTag);

    auto log = [&](const std::string& base, size_t added, bool built, Long64_t n){
        const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("[CutCache] %s {%s} <- base %s (+%zu atoms) %s %lld entries in %.1f s\n", tree.c_str(), Join(needAtoms, " ; ").c_str(), base.c_str(), added,
                    built ? "built" : "hit", (long long)n, s);
    };

    // nothing to cut and no columns: the source itself
    if (need.empty() && tag.empty()){
        CacheTree ct = OpenAsCacheTree(files.size() == 1 ? files[0] : source, files.size() == 1 ? files[0] : source, tree, {});
        log(ct.path, 0, false, ct.tree ? ct.tree->GetEntries() : 0);
        return ct;
    }

    for (int attempt = 0;; attempt++){
        // candidates: valid sidecars that are subsets of the request
        std::vector<Sidecar> cands;
        for (auto& sc : ScanSidecars(tree)){
            if (sc.kv.at("source") != (source) && sc.kv.at("source") != canon) continue;
            if (sc.kv.at("tree") != tree) continue;
            if (!StampValid(sc.kv.at("stamp"))){
                if (!LockHeld(sc.dir, sc.sig)) DeleteStale(sc);
                continue;
            }
            const std::string ctag = sc.kv.count("tag") ? sc.kv.at("tag") : "";
            if (!ctag.empty() && ctag != tag) continue;
            bool subset = true;
            for (const auto& k : sc.keys) if (!needSet.count(k)) { subset = false; break; }
            if (!subset) continue;
            // every opaque atom of the request must already be in the candidate: the builder can never apply one
            std::set<std::string> have(sc.keys.begin(), sc.keys.end());
            bool opaqueOk = true;
            for (const auto& k : needKeys) if (IsOpaque(k) && !have.count(k)) { opaqueOk = false; break; }
            if (!opaqueOk) continue;
            if (!Exists(sc.kv.at("file"))) { Remove(sc.cutsPath); continue; }
            cands.push_back(sc);
        }
        // exact hit
        for (const auto& sc : cands){
            const std::string ctag = sc.kv.count("tag") ? sc.kv.at("tag") : "";
            if (sc.keys.size() == needKeys.size() && ctag == tag){
                CacheTree ct = OpenAsCacheTree(sc.kv.at("file"), source, tree, sc.atoms);
                if (!ct.tree){ Remove(sc.cutsPath); continue; }
                log(ct.path, 0, false, ct.tree->GetEntries());
                return ct;
            }
        }
        // largest cached subset; tie: fewer entries, then smaller sig
        const Sidecar* best = nullptr;
        for (const auto& sc : cands){
            if (!best){ best = &sc; continue; }
            if (sc.keys.size() != best->keys.size()){ if (sc.keys.size() > best->keys.size()) best = &sc; continue; }
            const long long e1 = std::atoll(sc.kv.at("entries").c_str()), e2 = std::atoll(best->kv.at("entries").c_str());
            if (e1 != e2){ if (e1 < e2) best = &sc; continue; }
            if (sc.sig < best->sig) best = &sc;
        }
        std::vector<std::string> baseKeys;
        if (best) baseKeys = best->keys;
        const std::set<std::string> baseSet(baseKeys.begin(), baseKeys.end());
        std::vector<std::string> missingAtoms; // sorted key order
        for (const auto& kv : need) if (!baseSet.count(kv.first)) missingAtoms.push_back(kv.second);
        for (const auto& a : missingAtoms) if (IsOpaque(a)) Fatal("GetCutTree", "opaque atom %s is in no cache or seed of tree %s (run CutEvents.cpp)", a.c_str(), tree.c_str());

        // build under the lock
        MakeDirs(dir);
        const std::string lockPath = dir + "/" + sig + ".lock", outPath = dir + "/" + sig + ".root", cutsPath = dir + "/" + sig + ".cuts";
        if (!TryLock(lockPath)){ // someone else builds it: wait for the sidecar (or for the lock to vanish), then look up again
            long long waitedS = 0; bool skewWarned = false;
            while (Exists(lockPath) && !Exists(cutsPath)){
                long long age = (long long)std::time(nullptr) - FileMTime(lockPath); // seconds since the holder's last beat
                if (age < -120){ // the holder's clock is ahead of this node's: age is meaningless, wait on instead of taking the lock over
                    if (!skewWarned){ skewWarned = true; Warning("CutCache", "lock %s is dated %lld s in the future (clock skew between nodes): not treating it as stale", lockPath.c_str(), -age); }
                    age = 0;
                }
                if (age > LockStaleSeconds()) break;
                gSystem->Sleep(2000); waitedS += 2;
                const auto info = Lines(lockPath); // holder died on this host: stop waiting, the next TryLock takes the lock over
                if (waitedS % 60 == 0) Info("CutCache", "waiting %lld s for %s, held by %s (no beat for %lld s, stale after %lld s)", waitedS, lockPath.c_str(), info.empty() ? "?" : info[0].c_str(), age, LockStaleSeconds());
                std::istringstream is(info.empty() ? "" : info[0]);
                std::string host; long long pid = 0;
                if (is >> host >> pid && host == HostName() && pid > 0 && ::kill((pid_t)pid, 0) != 0 && errno == ESRCH) break;
            }
            continue;
        }
        LockBeat beat(lockPath); // keeps the lock fresh for the whole build; its destructor runs on every exit below
        if (Exists(cutsPath)){ Remove(lockPath); continue; } // built between the scan and the lock
        Remove(cutsPath); Remove(outPath);

        const std::string baseFile = best ? best->kv.at("file") : (files.size() == 1 ? files[0] : source);
        const std::string tmpOut = outPath + ".tmp." + HostName() + "." + std::to_string((long long)::getpid());
        std::string joined;
        for (size_t i = 0; i < missingAtoms.size(); i++) joined += (i ? " && (" : "(") + missingAtoms[i] + ")";
        const std::string btag = best && best->kv.count("tag") ? best->kv.at("tag") : "";
        const bool needAdd = !tag.empty() && btag != tag;
        // CUTCACHE_COPYTREE=1: single-threaded TTree::CopyTree instead of the RDataFrame Snapshot, for a ROOT whose Snapshot JIT does not run (the Mac cling).
        // It cannot add columns, so a cache that needs the ColumnAdder still takes the Snapshot.
        const char* copyEnv = std::getenv("CUTCACHE_COPYTREE");
        const bool copyTree = !needAdd && copyEnv && std::string(copyEnv) == "1";
        const bool wasMT = ROOT::IsImplicitMTEnabled();
        if (!copyTree && !wasMT) ROOT::EnableImplicitMT(NThreads());
        Long64_t entries = 0;
        if (copyTree){
            TChain chain(tree.c_str()); chain.Add(baseFile.c_str());
            TFile out(tmpOut.c_str(), "RECREATE");
            if (TTree* copy = chain.CopyTree(joined.c_str())) copy->Write();
            out.Close();
        } else {
            ROOT::RDataFrame df(tree, baseFile);
            ROOT::RDF::RNode node = df;
            if (!missingAtoms.empty()) node = node.Filter(joined);
            if (needAdd){
                if (!addCols){ Remove(lockPath); Fatal("GetCutTree", "tag %s requested without a ColumnAdder", tag.c_str()); }
                node = addCols(node);
            }
            node.Snapshot(tree, tmpOut); // all columns
        }
        if (!copyTree && !wasMT) ROOT::DisableImplicitMT();
        {
            TFile f(tmpOut.c_str(), "READ");
            TTree* t = f.IsZombie() ? nullptr : (TTree*)f.Get(tree.c_str());
            entries = t ? t->GetEntries() : 0;
        }
        std::rename(tmpOut.c_str(), outPath.c_str());
        std::vector<std::string> stampFiles = files; stampFiles.push_back("settings.csv");
        for (const auto& d : tagDeps) stampFiles.push_back(d);
        if (best) stampFiles.push_back(best->kv.at("file")); // a rebuilt base makes this cache stale
        WriteSidecar(cutsPath, outPath, source, tree, tag, entries, StampOf(stampFiles), need);
        Remove(lockPath);
        CacheTree ct = OpenAsCacheTree(outPath, source, tree, needAtoms);
        ct.built = true;
        log(best ? best->kv.at("file") : source, missingAtoms.size(), true, entries);
        return ct;
    }
}
