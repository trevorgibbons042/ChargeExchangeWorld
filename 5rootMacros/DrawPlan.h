#pragma once
// TTree::Draw compatible lazy requests, run in one multithreaded RDataFrame pass per tree.
//   book:   h = plan.H(tree, "expr>>name", selection, &model) (1D, fills the given model histogram), c = plan.C(tree, {"e1", "e2"}, selection)
//   run:    plan.Run()   one event loop per tree (ROOT::EnableImplicitMT(NThreads())), values collected with their entry numbers
//   render: h.Get(), c.V(0) ... the histograms are filled on the main thread in entry order, the order TTree::Draw fills in, so the stats sums match it.
// All TH1 filling, fits and graphics stay on the main thread; only the event loop is multithreaded. The entry numbers come from WithEntryNumber (FitCore.h), not rdfentry_, which is wrong under
// implicit MT on a tree with several clusters. A request that cannot be translated
// runs as the plain TTree::Draw on the main thread during Run(), with a one-time warning.
#include <TError.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "FitCore.h" // NThreads

struct HistReq; struct ColsReq;

class DrawPlan {
public:
    struct Hist  {std::shared_ptr<HistReq> r; TH1* Get() const;};
    struct Cols  {std::shared_ptr<ColsReq> r; Long64_t N() const; const double* V(int i) const;};

    Hist  H(TTree* t, const std::string& varexp, const std::string& selection, TH1* model);
    Cols  C(TTree* t, const std::vector<std::string>& exprs, const std::string& selection = "");
    void Run();
    static bool ForceFallback(){const char* e = std::getenv("DRAWPLAN_FALLBACK"); return e && std::string(e) == "1";} // test switch: every request runs as plain TTree::Draw

    struct Group {                       // requests on one tree with the same selection: one Filter node
        TTree* tree = nullptr; std::string sel;
        std::vector<std::string> exprs;  // column expressions, deduplicated
        std::vector<std::vector<double>> col; // after Run: entry-sorted, one vector per expression
        Long64_t n = 0;
        bool done = false, failed = false;
        int Index(const std::string& e){
            for (size_t i = 0; i < exprs.size(); i++) if (exprs[i] == e) return (int)i;
            exprs.push_back(e); return (int)exprs.size() - 1;
        }
        size_t stride() const {return 1 + exprs.size();}
    };
    std::vector<std::shared_ptr<Group>> groups;
    std::vector<std::shared_ptr<HistReq>> hists; std::vector<std::shared_ptr<ColsReq>> cols;
    bool unordered = false; // set before Run(): the rows of a many-cluster tree may come back in any order (no rdfentry_), so it is read on all threads; the columns of a row stay together
};

struct HistReq {
    TTree* t = nullptr; std::string varexp, sel, name; TH1* model = nullptr;
    std::shared_ptr<DrawPlan::Group> g; int ix = -1; // column index of x
    TH1* result = nullptr; bool done = false;
};
struct ColsReq  {TTree* t = nullptr; std::vector<std::string> exprs; std::string sel; std::shared_ptr<DrawPlan::Group> g; std::vector<int> ix;
                 std::vector<std::vector<double>> v; Long64_t n = 0; bool done = false;};

inline TH1* DrawPlan::Hist::Get() const {return r->result;}
inline Long64_t DrawPlan::Cols::N() const {return r->n;}
inline const double* DrawPlan::Cols::V(int i) const {return r->v.empty() ? r->g->col.at(r->ix.at(i)).data() : r->v.at(i).data();} // v only holds the fallback copy

namespace DrawPlanImpl {

inline std::string Trim(const std::string& s){
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) a++;
    while (b > a && std::isspace((unsigned char)s[b - 1])) b--;
    return s.substr(a, b - a);
}

// the selection in parentheses, "" when it is empty
inline std::string Effective(const std::string& sel){
    const std::string a = Trim(sel);
    return a.empty() ? "" : "(" + a + ")";
}

// the group of a tree and selection, made when missing
inline std::shared_ptr<DrawPlan::Group> PlainGroup(std::vector<std::shared_ptr<DrawPlan::Group>>& groups, TTree* t, const std::string& eff){
    for (auto& o : groups) if (o->tree == t && o->sel == eff) return o;
    auto g = std::make_shared<DrawPlan::Group>(); g->tree = t; g->sel = eff; groups.push_back(g);
    return g;
}

} // namespace DrawPlanImpl

inline DrawPlan::Hist DrawPlan::H(TTree* t, const std::string& varexp, const std::string& selection, TH1* model){
    using namespace DrawPlanImpl;
    auto r = std::make_shared<HistReq>();
    r->t = t; r->varexp = varexp; r->sel = selection; r->model = model;
    const size_t arrow = varexp.find(">>");
    const std::string left = Trim(arrow == std::string::npos ? varexp : varexp.substr(0, arrow));
    r->name = arrow == std::string::npos ? "" : Trim(varexp.substr(arrow + 2));
    if (!model || r->name.empty() || r->name.find('(') != std::string::npos || left.empty()){
        Error("DrawPlan::H", "unsupported varexp '%s' (needs expr>>name and a model)", varexp.c_str());
        throw std::runtime_error("DrawPlan: unsupported varexp");
    }
    auto g = PlainGroup(groups, t, Effective(selection));
    r->g = g;
    r->ix = g->Index(left);
    hists.push_back(r);
    return {r};
}

inline DrawPlan::Cols DrawPlan::C(TTree* t, const std::vector<std::string>& exprs, const std::string& selection){
    using namespace DrawPlanImpl;
    auto r = std::make_shared<ColsReq>();
    r->t = t; r->exprs = exprs; r->sel = selection;
    auto g = PlainGroup(groups, t, Effective(selection));
    r->g = g;
    for (const auto& e : exprs) r->ix.push_back(g->Index(e));
    cols.push_back(r);
    return {r};
}

namespace DrawPlanImpl {

inline void FallbackHist(HistReq& r){
    // Draw fills the object of the >>name in the directory: put the model there for the call, or it would make a new histogram and leave the model empty
    TDirectory* modelDir = r.model->GetDirectory();
    r.model->SetDirectory(gDirectory);
    r.t->Draw(r.varexp.c_str(), r.g->sel.c_str(), "goff");
    r.model->SetDirectory(modelDir);
    r.result = r.model; r.done = true;
}

} // namespace DrawPlanImpl

inline void DrawPlan::Run(){
    using namespace DrawPlanImpl;
    static bool warned = false;
    const bool wasMT = ROOT::IsImplicitMTEnabled();
    if (!wasMT) ROOT::EnableImplicitMT(NThreads());

    // 1. book one filter node + one row collector per group, all run together
    struct Booked {std::shared_ptr<Group> g; std::shared_ptr<ROOT::RDF::RResultPtr<std::vector<std::vector<double>>>> rows; bool ordFree = false;};
    std::vector<Booked> booked;
    std::vector<ROOT::RDF::RResultHandle> handles;
    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfs;
    auto pending = groups;
    for (auto& g : pending){
        if (g->done) continue;
        bool want = false;
        for (auto& r : hists) if (r->g == g && !r->done) want = true;
        for (auto& r : cols) if (r->g == g && !r->done) want = true;
        if (!want) {g->done = true; continue;}
        const bool ordFree = unordered && MultiCluster(g->tree); // no entry number is read: the rows of the tree are not put back in entry order
        if (ForceFallback()){g->failed = true; continue;}
        try {
            dfs.emplace_back(new ROOT::RDataFrame(*g->tree)); // the number of slots is fixed when the RDataFrame is built
            auto& hv = handles;
            ROOT::RDF::RNode node = *dfs.back();
            if (!ordFree) node = WithEntryNumber(node, dfs.back()->GetNSlots()); // exact entry numbers, before any Filter
            if (!g->sel.empty()) node = node.Filter(g->sel + " != 0");
            Booked b; b.g = g; b.ordFree = ordFree;
            std::string list = ordFree ? "ROOT::RVecD{0.0" : "ROOT::RVecD{(double)entryNo"; // column 0 is not read when the order is free
            for (size_t i = 0; i < g->exprs.size(); i++) list += ", (double)(" + g->exprs[i] + ")";
            list += "}";
            node = node.Define("__v", list);
            // Columnar per slot: column 0 is the entry, then the values. The merger appends column by column and frees each slot column at once,
            // and the sort below moves one column at a time, so the peak is about 1x the group, not the 2x of a flat buffer plus its sorted copy.
            if (const char* gb = std::getenv("DRAWPLAN_MEM_GB")){
                const double upper = (double)g->tree->GetEntries() * g->stride() * sizeof(double) / 1e9;
                if (upper > std::atof(gb)) Warning("DrawPlan", "group %s: up to %.1f GB if every entry passes (DRAWPLAN_MEM_GB=%s)", g->sel.c_str(), upper, gb);
            }
            using Cols2D = RowCols; // FitCore.h
            auto r = node.Aggregate([](Cols2D& acc, const ROOT::RVecD& v){RowsAppend(acc, v);}, [](std::vector<Cols2D>& slots){RowsMerge(slots);}, "__v", Cols2D());
            b.rows = std::make_shared<ROOT::RDF::RResultPtr<Cols2D>>(r); hv.push_back(r);
            booked.push_back(b);
        } catch (const std::exception& e){
            if (!warned){warned = true; Warning("DrawPlan", "cannot book %s (%s): plain TTree::Draw", g->sel.c_str(), e.what());}
            g->failed = true;
        }
    }
    bool ranOk = true;
    if (!handles.empty()){
        try {
            ROOT::RDF::RunGraphs(handles);
        }
        catch (const std::exception& e){ranOk = false; Warning("DrawPlan", "event loop failed (%s): plain TTree::Draw for every request", e.what());}
    }
    if (!wasMT) ROOT::DisableImplicitMT();
    if (!ranOk) for (auto& b : booked) b.g->failed = true;

    // 2. sort every collected group by entry
    if (ranOk) for (auto& b : booked){
        const size_t stride = b.g->stride();
        std::vector<std::vector<double>> src; src.swap(**b.rows); // take the buffers, so each column can be freed as it is sorted
        src.resize(stride);
        const size_t total = src[0].size();
        std::vector<size_t> order(total); std::iota(order.begin(), order.end(), size_t(0));
        if (!b.ordFree){ // an order-free group keeps the rows as the threads delivered them
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
        TH1* h = r->model;
        h->Reset();
        const Group& g = *r->g;
        const double* xc = g.col[r->ix].data();
        for (Long64_t k = 0; k < g.n; k++) h->Fill(xc[k]);
        r->result = h; r->done = true;
    }
    for (auto& r : cols){
        if (r->done) continue;
        if (!r->g->done || r->g->failed){
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
    groups.clear(); hists.clear(); cols.clear(); // handles stay valid (shared_ptr)
}
