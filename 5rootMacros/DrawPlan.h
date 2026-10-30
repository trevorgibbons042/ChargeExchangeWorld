#pragma once
// TTree::Draw compatible lazy requests, run in one multithreaded RDataFrame pass per tree (spec section B).
//   book:   h = plan.H(tree, "e1:e2>>name(nx,lo,hi,ny,lo,hi)", selection, listCut), n = plan.N(tree, selection), c = plan.C(tree, {"e1", "e2"}, selection)
//   run:    plan.Run()   one event loop per tree (ROOT::EnableImplicitMT(NThreads())), values collected with their entry numbers
//   render: h.Get(), n.Get(), c.V(0) ... the histograms are filled on the main thread in entry order, the order TTree::Draw fills in, so the stats sums match it.
// All TH1/TH2 filling, fits and graphics stay on the main thread; only the event loop is multithreaded. A request that cannot be translated (or if rdfentry_ is not the TTree entry number)
// runs as the plain TTree::Draw on the main thread during Run(), with a one-time warning. The fallback draws with the effective selection (selection and listCut), so the caller does not
// need to set the TEntryList on the tree; listCut is the string the list stands for.
#include <TError.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TEntryList.h>
#include <TLeaf.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RDFHelpers.hxx>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "FitCore.h" // NThreads

struct HistReq; struct CountReq; struct ColsReq;

class DrawPlan {
public:
    struct Hist  {std::shared_ptr<HistReq> r; TH1* Get() const;};
    struct Count {std::shared_ptr<CountReq> r; Long64_t Get() const;};
    struct Cols  {std::shared_ptr<ColsReq> r; Long64_t N() const; const double* V(int i) const;};

    Hist  H(TTree* t, const std::string& varexp, const std::string& selection = "", const std::string& listCut = "", TH1* model = nullptr, bool weighted = false);
    Count N(TTree* t, const std::string& selection, const std::string& listCut = "");
    Cols  C(TTree* t, const std::vector<std::string>& exprs, const std::string& selection = "", const std::string& listCut = "");
    void Run();
    static DrawPlan& Default(){static DrawPlan p; return p;}
    static bool ForceFallback(){const char* e = std::getenv("DRAWPLAN_FALLBACK"); return e && std::string(e) == "1";} // test switch: every request runs as plain TTree::Draw

    struct Group {                       // requests on one tree with the same effective selection: one Filter node
        TTree* tree = nullptr; std::string sel; bool weighted = false;
        std::vector<std::string> exprs;  // column expressions, deduplicated
        bool needRows = false;
        std::vector<std::vector<double>> col; // after Run: entry-sorted, one vector per expression, the weight last if weighted
        Long64_t n = 0;
        bool done = false, failed = false;
        int Index(const std::string& e){
            for (size_t i = 0; i < exprs.size(); i++) if (exprs[i] == e) return (int)i;
            exprs.push_back(e); return (int)exprs.size() - 1;
        }
        size_t stride() const {return 1 + exprs.size() + (weighted ? 1 : 0);}
    };
    std::vector<std::shared_ptr<Group>> groups;
    std::vector<std::shared_ptr<HistReq>> hists; std::vector<std::shared_ptr<CountReq>> counts; std::vector<std::shared_ptr<ColsReq>> cols;
    std::map<TTree*, int> entryOk; // 1: rdfentry_ is the TTree entry number (checked), 0: not
};

struct HistReq {
    TTree* t = nullptr; std::string varexp, sel, listCut, name; TH1* model = nullptr; bool weighted = false;
    std::shared_ptr<DrawPlan::Group> g; int ix = -1, iy = -1; // column index of x and y
    int dim = 1; std::vector<double> bins; // explicit bins, empty with a model
    TH1* result = nullptr; bool done = false;
};
struct CountReq {TTree* t = nullptr; std::string sel, listCut; std::shared_ptr<DrawPlan::Group> g; Long64_t n = 0; bool done = false;};
struct ColsReq  {TTree* t = nullptr; std::vector<std::string> exprs; std::string sel, listCut; std::shared_ptr<DrawPlan::Group> g; std::vector<int> ix;
                 std::vector<std::vector<double>> v; Long64_t n = 0; bool done = false;};

inline TH1* DrawPlan::Hist::Get() const {return r->result;}
inline Long64_t DrawPlan::Count::Get() const {return r->n;}
inline Long64_t DrawPlan::Cols::N() const {return r->n;}
inline const double* DrawPlan::Cols::V(int i) const {return r->v.empty() ? r->g->col.at(r->ix.at(i)).data() : r->v.at(i).data();} // v only holds the fallback copy

namespace DrawPlanImpl {

inline std::string Trim(const std::string& s){
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) a++;
    while (b > a && std::isspace((unsigned char)s[b - 1])) b--;
    return s.substr(a, b - a);
}

// "a:b" at depth 0, not "::", not inside quotes
inline std::vector<std::string> SplitColon(const std::string& e){
    std::vector<std::string> out; std::string cur; int depth = 0; char q = 0;
    for (size_t i = 0; i < e.size(); i++){
        const char c = e[i];
        if (q){cur += c; if (c == q) q = 0; continue;}
        if (c == '"' || c == '\''){q = c; cur += c; continue;}
        if (c == '(' || c == '[') depth++; else if (c == ')' || c == ']') depth--;
        if (c == ':' && depth == 0){
            if (i + 1 < e.size() && e[i + 1] == ':'){cur += "::"; i++; continue;}
            out.push_back(Trim(cur)); cur.clear(); continue;
        }
        cur += c;
    }
    out.push_back(Trim(cur));
    return out;
}

inline std::string Effective(const std::string& sel, const std::string& listCut){
    std::string a = Trim(sel), b = Trim(listCut);
    if (a.empty()) return b.empty() ? "" : "(" + b + ")";
    if (b.empty()) return "(" + a + ")";
    return "(" + a + ") && (" + b + ")";
}

inline bool ParseBins(const std::string& spec, std::vector<double>& bins){ // "(nx,lo,hi[,ny,lo,hi])"
    bins.clear();
    if (spec.size() < 2 || spec.front() != '(' || spec.back() != ')') return false;
    std::stringstream ss(spec.substr(1, spec.size() - 2)); std::string item;
    while (std::getline(ss, item, ',')){
        char* end = nullptr; const double v = std::strtod(item.c_str(), &end);
        if (end == item.c_str()) return false;
        bins.push_back(v);
    }
    return bins.size() == 3 || bins.size() == 6;
}

} // namespace DrawPlanImpl

inline DrawPlan::Hist DrawPlan::H(TTree* t, const std::string& varexp, const std::string& selection, const std::string& listCut, TH1* model, bool weighted){
    using namespace DrawPlanImpl;
    auto r = std::make_shared<HistReq>();
    r->t = t; r->varexp = varexp; r->sel = selection; r->listCut = listCut; r->model = model; r->weighted = weighted;
    const size_t arrow = varexp.find(">>");
    std::string left = arrow == std::string::npos ? varexp : varexp.substr(0, arrow), right = arrow == std::string::npos ? "" : Trim(varexp.substr(arrow + 2));
    const size_t par = right.find('(');
    r->name = Trim(par == std::string::npos ? right : right.substr(0, par));
    const auto parts = SplitColon(left);
    const bool binsOk = par == std::string::npos ? model != nullptr : ParseBins(right.substr(par), r->bins);
    if (arrow == std::string::npos || r->name.empty() || parts.size() < 1 || parts.size() > 2 || !binsOk || (par != std::string::npos && model)){
        Error("DrawPlan::H", "unsupported varexp '%s' (needs >>name(bins) or a model with >>name)", varexp.c_str());
        throw std::runtime_error("DrawPlan: unsupported varexp");
    }
    r->dim = (int)parts.size();
    if (!r->bins.empty() && (int)r->bins.size() != 3 * r->dim){Error("DrawPlan::H", "bins do not match the dimension: %s", varexp.c_str()); throw std::runtime_error("DrawPlan: bins");}
    auto g = std::make_shared<Group>();
    g->tree = t; g->weighted = weighted;
    // a weight is multiplied by the list cut, not and-ed with it
    g->sel = weighted && !Trim(selection).empty() && !Trim(listCut).empty() ? "(" + Trim(selection) + ")*(" + Trim(listCut) + ")" : Effective(selection, listCut);
    bool shared = false;
    for (auto& o : groups) if (!shared && o->tree == t && o->sel == g->sel && o->weighted == weighted){g = o; shared = true;}
    if (!shared) groups.push_back(g);
    g->needRows = true;
    r->g = g;
    // x is the last expression, y the first of two ("y:x"), like TTree::Draw
    r->ix = g->Index(parts.back());
    if (r->dim == 2) r->iy = g->Index(parts.front());
    hists.push_back(r);
    return {r};
}

inline DrawPlan::Count DrawPlan::N(TTree* t, const std::string& selection, const std::string& listCut){
    using namespace DrawPlanImpl;
    auto r = std::make_shared<CountReq>();
    r->t = t; r->sel = selection; r->listCut = listCut;
    const std::string eff = Effective(selection, listCut);
    std::shared_ptr<Group> g;
    for (auto& o : groups) if (o->tree == t && o->sel == eff && !o->weighted){g = o; break;}
    if (!g){g = std::make_shared<Group>(); g->tree = t; g->sel = eff; groups.push_back(g);}
    r->g = g; counts.push_back(r);
    return {r};
}

inline DrawPlan::Cols DrawPlan::C(TTree* t, const std::vector<std::string>& exprs, const std::string& selection, const std::string& listCut){
    using namespace DrawPlanImpl;
    auto r = std::make_shared<ColsReq>();
    r->t = t; r->exprs = exprs; r->sel = selection; r->listCut = listCut;
    const std::string eff = Effective(selection, listCut);
    std::shared_ptr<Group> g;
    for (auto& o : groups) if (o->tree == t && o->sel == eff && !o->weighted){g = o; break;}
    if (!g){g = std::make_shared<Group>(); g->tree = t; g->sel = eff; groups.push_back(g);}
    g->needRows = true;
    r->g = g;
    for (const auto& e : exprs) r->ix.push_back(g->Index(e));
    cols.push_back(r);
    return {r};
}

namespace DrawPlanImpl {

// Is rdfentry_ the TTree entry number of this tree under the current IMT? Checked once per tree on a sample of entries against TTree::GetEntry.
inline bool EntryNumberOk(DrawPlan& plan, TTree* t){
    auto it = plan.entryOk.find(t);
    if (it != plan.entryOk.end()) return it->second == 1;
    bool ok = false;
    try {
        // Up to 3 scalar leaves, continuous ones (Double_t, Float_t) first: a leaf like fEvent repeats over neighbouring entries, so a small entry offset would pass on it alone.
        TObjArray* leaves = t->GetListOfLeaves();
        std::vector<TLeaf*> chk;
        for (const bool continuous : {true, false})
            for (int i = 0; i < leaves->GetEntries() && chk.size() < 3; i++){
                TLeaf* l = (TLeaf*)leaves->At(i);
                const std::string tn = l->GetTypeName();
                const bool isCont = tn == "Double_t" || tn == "Float_t", isInt = tn == "Int_t" || tn == "Long64_t";
                if (l->GetLeafCount() == nullptr && l->GetLen() == 1 && (continuous ? isCont : isInt)) chk.push_back(l);
            }
        if (chk.empty() || t->GetEntries() == 0) throw std::runtime_error("no scalar leaf");
        const Long64_t n = t->GetEntries();
        const Long64_t step = std::max<Long64_t>(1, n / 997);
        ROOT::RDataFrame df(*t);
        ROOT::RDF::RNode sample = df.Filter("rdfentry_ % " + std::to_string(step) + " == 0");
        for (size_t c = 0; c < chk.size(); c++) sample = sample.Define("__chk" + std::to_string(c), "(double)(" + std::string(chk[c]->GetName()) + ")");
        auto es = sample.Take<ULong64_t>("rdfentry_");
        std::vector<ROOT::RDF::RResultPtr<std::vector<double>>> vs;
        for (size_t c = 0; c < chk.size(); c++) vs.push_back(sample.Take<double>("__chk" + std::to_string(c)));
        TEntryList* saved = t->GetEntryList(); t->SetEntryList(nullptr);
        {std::vector<ROOT::RDF::RResultHandle> hh{es}; for (auto& v : vs) hh.push_back(v); ROOT::RDF::RunGraphs(hh);}
        ok = !es->empty();
        for (auto& v : vs) ok = ok && es->size() == v->size();
        for (size_t i = 0; ok && i < es->size(); i++){
            t->GetEntry((*es)[i]);
            for (size_t c = 0; ok && c < chk.size(); c++){
                const double a = chk[c]->GetValue(), b = (*vs[c])[i];
                if (a != b && !(std::isnan(a) && std::isnan(b))) ok = false;
            }
        }
        t->SetEntryList(saved);
    } catch (const std::exception& e){ Warning("DrawPlan", "entry number check failed (%s): plain TTree::Draw", e.what()); ok = false; }
    plan.entryOk[t] = ok ? 1 : 0;
    if (!ok) Warning("DrawPlan", "rdfentry_ is not the TTree entry number of %s: its requests run as plain TTree::Draw", t->GetName());
    return ok;
}

// Fill like TTree::Draw: line, fill and marker attributes of the tree, in gDirectory, an old object of that name replaced.
inline TH1* MakeHist(HistReq& r){
    if (r.model) return r.model;
    if (TObject* old = gDirectory ? gDirectory->Get(r.name.c_str()) : nullptr) delete old;
    std::string title = r.varexp.substr(0, r.varexp.find(">>"));
    if (!r.sel.empty()) title += " {" + r.sel + "}";
    TH1* h;
    if (r.dim == 1) h = new TH1F(r.name.c_str(), title.c_str(), (int)r.bins[0], r.bins[1], r.bins[2]);
    else h = new TH2F(r.name.c_str(), title.c_str(), (int)r.bins[0], r.bins[1], r.bins[2], (int)r.bins[3], r.bins[4], r.bins[5]);
    h->SetDirectory(gDirectory);
    r.t->TAttLine::Copy(*h); r.t->TAttFill::Copy(*h); r.t->TAttMarker::Copy(*h);
    return h;
}

// The fallback ignores the tree's entry list like the event loop does (the list is in the selection as listCut).
struct NoEntryList {
    TTree* t; TEntryList* l;
    explicit NoEntryList(TTree* tree) : t(tree), l(tree->GetEntryList()) {t->SetEntryList(nullptr);}
    ~NoEntryList(){t->SetEntryList(l);}
};

inline void FallbackHist(HistReq& r){
    NoEntryList noList(r.t);
    // Draw fills the object of the >>name in the directory: put a model there for the call, or it would make a new histogram and leave the model empty
    TDirectory* modelDir = r.model ? r.model->GetDirectory() : nullptr;
    if (r.model) r.model->SetDirectory(gDirectory);
    r.t->Draw(r.varexp.c_str(), r.g->sel.c_str(), "goff"); // the effective selection: listCut is not on the tree
    if (r.model) r.model->SetDirectory(modelDir);
    r.result = r.model ? r.model : (TH1*)(gDirectory->Get(r.name.c_str()) ? gDirectory->Get(r.name.c_str()) : r.t->GetHistogram());
    r.done = true;
}

} // namespace DrawPlanImpl

inline void DrawPlan::Run(){
    using namespace DrawPlanImpl;
    static bool warned = false;
    const bool wasMT = ROOT::IsImplicitMTEnabled();
    if (!wasMT) ROOT::EnableImplicitMT(NThreads());

    // 1. book one filter node + one row collector (or Count) per group, all run together
    struct Booked {std::shared_ptr<Group> g; std::vector<std::vector<double>>* buf = nullptr; std::shared_ptr<ROOT::RDF::RResultPtr<ULong64_t>> cnt; std::shared_ptr<ROOT::RDF::RResultPtr<std::vector<std::vector<double>>>> rows;};
    std::vector<Booked> booked;
    std::vector<ROOT::RDF::RResultHandle> handles;
    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfs;
    std::map<TTree*, TEntryList*> savedLists;
    auto pending = groups;
    for (auto& g : pending){
        if (g->done) continue;
        bool want = false;
        for (auto& r : hists) if (r->g == g && !r->done) want = true;
        for (auto& r : counts) if (r->g == g && !r->done) want = true;
        for (auto& r : cols) if (r->g == g && !r->done) want = true;
        if (!want) {g->done = true; continue;}
        if (ForceFallback() || !EntryNumberOk(*this, g->tree)){g->failed = true; continue;}
        try {
            if (!savedLists.count(g->tree)){savedLists[g->tree] = g->tree->GetEntryList(); g->tree->SetEntryList(nullptr);} // the list is in the selection (listCut)
            dfs.emplace_back(new ROOT::RDataFrame(*g->tree));
            ROOT::RDF::RNode node = *dfs.back();
            std::string wexpr;
            if (g->weighted){ // the selection is the weight: entries with weight 0 are dropped
                node = node.Define("__w", "(double)(" + g->sel + ")").Filter("__w != 0");
            } else if (!g->sel.empty()) node = node.Filter(g->sel + " != 0");
            Booked b; b.g = g;
            if (!g->needRows){ // N only
                b.cnt = std::make_shared<ROOT::RDF::RResultPtr<ULong64_t>>(node.Count()); handles.push_back(*b.cnt);
            } else {
                std::string list = "ROOT::RVecD{(double)rdfentry_";
                for (size_t i = 0; i < g->exprs.size(); i++) list += ", (double)(" + g->exprs[i] + ")";
                if (g->weighted) list += ", __w";
                list += "}";
                node = node.Define("__v", list);
                // Columnar per slot: column 0 is the entry, then the values, the weight last. The merger appends column by column and frees each slot column at once,
                // and the sort below moves one column at a time, so the peak is about 1x the group, not the 2x of a flat buffer plus its sorted copy.
                if (const char* gb = std::getenv("DRAWPLAN_MEM_GB")){
                    const double upper = (double)g->tree->GetEntries() * g->stride() * sizeof(double) / 1e9;
                    if (upper > std::atof(gb)) Warning("DrawPlan", "group %s: up to %.1f GB if every entry passes (DRAWPLAN_MEM_GB=%s)", g->sel.c_str(), upper, gb);
                }
                using Cols2D = std::vector<std::vector<double>>;
                auto r = node.Aggregate([](Cols2D& acc, const ROOT::RVecD& v){
                                            if (acc.size() < v.size()) acc.resize(v.size());
                                            for (size_t c = 0; c < v.size(); c++) acc[c].push_back(v[c]);
                                        },
                                        [](std::vector<Cols2D>& slots){
                                            Cols2D& dst = slots[0];
                                            size_t ncol = 0; for (auto& s : slots) ncol = std::max(ncol, s.size());
                                            dst.resize(ncol);
                                            for (size_t c = 0; c < ncol; c++){
                                                size_t tot = 0; for (auto& s : slots) if (c < s.size()) tot += s[c].size();
                                                dst[c].reserve(tot);
                                                for (size_t i = 1; i < slots.size(); i++) if (c < slots[i].size()){
                                                    dst[c].insert(dst[c].end(), slots[i][c].begin(), slots[i][c].end());
                                                    std::vector<double>().swap(slots[i][c]); // clear() would keep the capacity
                                                }
                                            }
                                        },
                                        "__v", Cols2D());
                b.rows = std::make_shared<ROOT::RDF::RResultPtr<Cols2D>>(r); handles.push_back(r);
            }
            booked.push_back(b);
        } catch (const std::exception& e){
            if (!warned){warned = true; Warning("DrawPlan", "cannot book %s (%s): plain TTree::Draw", g->sel.c_str(), e.what());}
            g->failed = true;
        }
    }
    bool ranOk = true;
    if (!handles.empty()){
        try { ROOT::RDF::RunGraphs(handles); }
        catch (const std::exception& e){ranOk = false; Warning("DrawPlan", "event loop failed (%s): plain TTree::Draw for every request", e.what());}
    }
    for (auto& kv : savedLists) kv.first->SetEntryList(kv.second);
    if (!wasMT) ROOT::DisableImplicitMT();
    if (!ranOk) for (auto& b : booked) b.g->failed = true;

    // 2. sort every collected group by entry
    if (ranOk) for (auto& b : booked){
        if (!b.g->needRows){b.g->n = (Long64_t)b.cnt->GetValue(); b.g->done = true; continue;}
        const size_t stride = b.g->stride();
        std::vector<std::vector<double>> src; src.swap(**b.rows); // take the buffers, so each column can be freed as it is sorted
        src.resize(stride);
        const size_t total = src[0].size();
        std::vector<size_t> order(total); std::iota(order.begin(), order.end(), size_t(0));
        {
            const std::vector<double>& entry = src[0];
            std::sort(order.begin(), order.end(), [&](size_t x, size_t y){return entry[x] < entry[y];});
        }
        auto& col = b.g->col; col.assign(stride - 1, {});
        for (size_t c = 1; c < stride; c++){
            std::vector<double> sorted(total);
            for (size_t k = 0; k < total; k++) sorted[k] = src[c][order[k]];
            col[c - 1].swap(sorted);
            std::vector<double>().swap(src[c]);
        }
        std::vector<double>().swap(src[0]); std::vector<size_t>().swap(order);
        b.g->n = (Long64_t)total; b.g->done = true;
        Info("DrawPlan", "group %s: %zu rows x %zu columns = %.2f GB", b.g->sel.empty() ? "(no cut)" : b.g->sel.c_str(), total, stride - 1, (double)total * (stride - 1) * sizeof(double) / 1e9);
    }

    // 3. results, in the order of booking
    for (auto& r : hists){
        if (r->done) continue;
        if (!r->g->done || r->g->failed){FallbackHist(*r); continue;}
        TH1* h = MakeHist(*r);
        if (r->model) h->Reset();
        const Group& g = *r->g;
        const bool two = r->dim == 2;
        const double* xc = g.col[r->ix].data(); const double* yc = two ? g.col[r->iy].data() : nullptr;
        const double* wc = g.weighted ? g.col.back().data() : nullptr;
        for (Long64_t k = 0; k < g.n; k++){
            const double w = wc ? wc[k] : 1.0;
            if (two) ((TH2*)h)->Fill(xc[k], yc[k], w); else h->Fill(xc[k], w);
        }
        r->result = h; r->done = true;
    }
    for (auto& r : counts){
        if (r->done) continue;
        if (!r->g->done || r->g->failed){NoEntryList noList(r->t); r->n = r->t->GetEntries(r->g->sel.c_str()); r->done = true; continue;}
        r->n = r->g->n; r->done = true;
    }
    for (auto& r : cols){
        if (r->done) continue;
        if (!r->g->done || r->g->failed){
            NoEntryList noList(r->t);
            r->v.assign(r->exprs.size(), {});
            for (size_t i = 0; i < r->exprs.size(); i++){
                r->t->SetEstimate(r->t->GetEntries() + 1);
                const Long64_t n = r->t->Draw(r->exprs[i].c_str(), r->g->sel.c_str(), "goff");
                r->v[i].assign(r->t->GetV1(), r->t->GetV1() + n); r->n = n;
            }
            r->done = true; continue;
        }
        r->n = r->g->n; r->done = true; // V() reads the group columns, no copy
    }
    for (auto& g : pending){ // columns nobody reads through V() are done with
        bool used = false;
        for (auto& r : cols) if (r->g == g) used = true;
        if (!used) std::vector<std::vector<double>>().swap(g->col);
    }
    groups.clear(); hists.clear(); counts.clear(); cols.clear(); // handles stay valid (shared_ptr)
}

