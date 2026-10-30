#pragma once
#include <TDirectory.h>
#include <TError.h>
#include <TFriendElement.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TObject.h>
#include <TNamed.h>
#include <TAttLine.h>
#include <TAttFill.h>
#include <TAttMarker.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <map>
#include <cstdlib>
#include <functional>
#include <regex>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// Multithreaded replacement for TTree::Draw(exprs, cut, "goff") and for GetEntry loops: one RDataFrame event loop on all the threads (ROOT::EnableImplicitMT(NThreads()) first, FitCore.h),
// on one thread for a tree with several clusters (rdfentry_ is wrong under implicit MT there, SerialScope in FitCore.h)
// copies expressions of the tree out as columns. Like TTree::Draw it sees the TEntryList set on t (ApplyCut, CutFiles.h) and the friends of t, written "Friend.branch" in the expressions.
// The rows come back in entry order, one vector per expression, so a fill or a sum over them runs in the same order as the serial Draw loop did.
// T is double or std::string (CECreator: "std::string(CECreator.data())"). filter is an optional row selection, a C++ expression like a TTree::Draw cut.
// Every call JIT compiles its expressions: ask for all the columns of a tree at once.
inline std::string ColumnsExpr(TTree* t, std::string expr, std::map<std::string, std::string>& aliases){
    // a TTree alias (SetAlias, e.g. raw_k of Plot3DRawAliases) is not a column of the RDataFrame: write its expansion in its place, again until no alias is left (ok_9 -> TMath::Finite(raw_9) -> ...)
    if (t->GetListOfAliases() && t->GetListOfAliases()->GetSize() > 0){
        for (int round = 0; round < 8; round++){
            bool changed = false;
            for (auto* o : *t->GetListOfAliases()){
                const std::string name = o->GetName();
                if (t->GetBranch(name.c_str()) || expr.find(name) == std::string::npos) continue;
                const std::string out = std::regex_replace(expr, std::regex("\\b" + name + "\\b"), "(" + std::string(((TNamed*)o)->GetTitle()) + ")");
                if (out != expr){expr = out; changed = true;}
            }
            if (!changed) break;
        }
        // TMath::Finite(int branch) is ambiguous (Double_t or Float_t overload) in the JIT: the argument as a double
        expr = std::regex_replace(expr, std::regex("TMath::Finite\\("), "TMath::Finite((double)");
    }
    // a friend's branch is the column "Friend.branch" of the RDataFrame, which an expression cannot spell: it goes through the alias Friend__branch
    if (t->GetListOfFriends()) for (auto* o : *t->GetListOfFriends()){
        const std::string name = ((TFriendElement*)o)->GetTree()->GetName();
        const std::regex re("\\b" + name + "\\.([A-Za-z_][A-Za-z0-9_]*)");
        for (std::smatch m; std::regex_search(expr, m, re);){
            aliases[name + "__" + m[1].str()] = name + "." + m[1].str();
            expr = m.prefix().str() + name + "__" + m[1].str() + m.suffix().str();
        }
    }
    return expr;
}

#include "FitCore.h" // SerialScope: a multi-cluster tree is read on one thread, rdfentry_ is wrong under implicit MT there

template <class T = double>
inline std::vector<std::vector<T>> Columns(TTree* t, const std::vector<std::string>& exprs, const std::string& filter = ""){
    static_assert(std::is_same<T, double>::value || std::is_same<T, std::string>::value, "Columns: T is double or std::string");
    SerialScope serial(t); // the RDataFrame and its lazy Takes below all run inside it
    const std::string type = std::is_same<T, double>::value ? "double" : "std::string";
    using Row = std::pair<ULong64_t, T>; // (entry, value): the entries put the rows of the threads back in order
    std::map<std::string, std::string> aliases;
    const std::string cut = filter.empty() ? "" : ColumnsExpr(t, filter, aliases);
    std::vector<std::string> defs;
    for (const auto& e : exprs) defs.push_back("std::make_pair((ULong64_t)rdfentry_, (" + type + ")(" + ColumnsExpr(t, e, aliases) + "))");
    ROOT::RDataFrame df0(*t);
    ROOT::RDF::RNode df = df0;
    for (const auto& [alias, column] : aliases) df = df.Alias(alias, column);
    if (!cut.empty()) df = df.Filter(cut, "columns filter");
    std::vector<ROOT::RDF::RResultPtr<std::vector<Row>>> taken;
    for (size_t i = 0; i < exprs.size(); i++){
        const std::string col = "tcol" + std::to_string(i);
        df = df.Define(col, defs[i]);
        taken.push_back(df.Take<Row>(col));
    }
    std::vector<std::vector<T>> cols(exprs.size());
    for (size_t i = 0; i < exprs.size(); i++){
        std::vector<Row>& rows = const_cast<std::vector<Row>&>(*taken[i]);
        std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b){return a.first < b.first;});
        cols[i].reserve(rows.size());
        for (auto& r : rows) cols[i].push_back(std::move(r.second));
    }
    return cols;
}

// The expressions and the histogram of a TTree::Draw string "y:x>>name(nx,lo,hi[,ny,lo,hi])": exprs in Draw order (y first), target the part after ">>".
struct DrawSpec {std::vector<std::string> exprs; std::string target;};
inline DrawSpec SplitDraw(const std::string& drawExpr){
    DrawSpec d;
    const size_t arrow = drawExpr.find(">>");
    const std::string e = drawExpr.substr(0, arrow);
    if (arrow != std::string::npos) d.target = drawExpr.substr(arrow + 2);
    int depth = 0; std::string cur;
    for (size_t i = 0; i < e.size(); i++){
        const char c = e[i];
        if (c == '(' || c == '[') depth++;
        if (c == ')' || c == ']') depth--;
        if (c == ':' && e[i + 1] == ':'){cur += "::"; i++; continue;} // TMath::Pi()
        if (c == ':' && depth == 0){d.exprs.push_back(cur); cur.clear(); continue;}
        cur += c;
    }
    d.exprs.push_back(cur);
    return d;
}

// The empty TH1F or TH2F TTree::Draw makes from the target "name(nx,lo,hi)" or "name(nx,lo,hi,ny,lo,hi)": in gDirectory, replacing a histogram of that name. nullptr if the target has no binning.
inline TH1* MakeHist(const std::string& target){
    const size_t open = target.find('(');
    if (open == std::string::npos || target.back() != ')'){Error("MakeHist", "no binning in %s", target.c_str()); return nullptr;}
    const std::string name = target.substr(0, open);
    std::vector<double> a;
    for (size_t p = open + 1; p < target.size();){
        size_t q = target.find_first_of(",)", p);
        a.push_back(std::atof(target.substr(p, q - p).c_str()));
        p = q + 1;
    }
    delete gDirectory->Get(name.c_str());
    if (a.size() == 3) return new TH1F(name.c_str(), name.c_str(), (int)a[0], a[1], a[2]);
    if (a.size() == 6) return new TH2F(name.c_str(), name.c_str(), (int)a[0], a[1], a[2], (int)a[3], a[4], a[5]);
    Error("MakeHist", "binning of %s is not 3 or 6 numbers", target.c_str());
    return nullptr;
}

// The entries of t that pass cut (a TTree::GetEntries(cut) with all the threads): the entry list set on t and the friends are honoured as in Columns. Empty cut or "1": every entry.
inline Long64_t CountMT(TTree* t, const std::string& cut = ""){
    std::map<std::string, std::string> aliases;
    const bool all = cut.empty() || cut == "1";
    const std::string c = all ? "" : ColumnsExpr(t, cut, aliases);
    ROOT::RDataFrame df0(*t);
    ROOT::RDF::RNode df = df0;
    for (const auto& [alias, column] : aliases) df = df.Alias(alias, column);
    if (!all) df = df.Filter(c, "count filter");
    return (Long64_t)df.Count().GetValue();
}

// TTree::Draw(drawExpr, cut, "goff") as one event loop on all the threads, no row vectors: drawExpr is "x>>name(nx,lo,hi)" or "y:x>>name(nx,lo,hi,ny,lo,hi)". Returns the TH1F / TH2F that MakeHist makes
// (in gDirectory, named name), filled with the bins and stats of the threads merged, with the line/fill/marker attributes of t as TTree::Draw gives them. It is not drawn: call h->Draw(opt) ("same" for an overlay).
// Bin contents are the same as the serial Draw; mean and RMS can differ in the last digits (the sums are added up per thread).
// Several HistMT of one tree in one event loop: reqs are (drawExpr, cut) pairs, the histograms come back in the same order (nullptr for a target without binning). The names in the
// drawExprs must differ, MakeHist replaces a histogram of the same name.
inline std::vector<TH1*> HistsMT(TTree* t, const std::vector<std::pair<std::string, std::string>>& reqs){
    std::map<std::string, std::string> aliases;
    std::vector<DrawSpec> specs;
    std::vector<TH1*> hs;
    std::vector<std::string> cuts;
    std::vector<std::vector<std::string>> defs;
    for (const auto& [drawExpr, cut] : reqs){
        specs.push_back(SplitDraw(drawExpr));
        hs.push_back(MakeHist(specs.back().target));
        cuts.push_back(cut.empty() || cut == "1" ? "" : ColumnsExpr(t, cut, aliases));
        defs.emplace_back();
        for (const auto& e : specs.back().exprs) defs.back().push_back("(double)(" + ColumnsExpr(t, e, aliases) + ")");
    }
    ROOT::RDataFrame df0(*t);
    ROOT::RDF::RNode df = df0;
    for (const auto& [alias, column] : aliases) df = df.Alias(alias, column);
    std::vector<std::function<TH1*()>> src(reqs.size()); // the first call runs the event loop of every booked histogram
    for (size_t j = 0; j < reqs.size(); j++){
        TH1* h = hs[j];
        if (!h) continue;
        ROOT::RDF::RNode d = df;
        if (!cuts[j].empty()) d = d.Filter(cuts[j], "hist filter");
        const std::string col = "hcol" + std::to_string(j) + "_";
        for (size_t i = 0; i < defs[j].size(); i++) d = d.Define(col + std::to_string(i), defs[j][i]);
        const std::string modelName = std::string(h->GetName()) + "_mt";
        if (specs[j].exprs.size() == 1 && h->GetDimension() == 1){
            auto r = d.Histo1D(ROOT::RDF::TH1DModel(modelName.c_str(), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()), col + "0");
            src[j] = [r]() mutable {return (TH1*)r->Clone();}; // mutable: operator-> of a result pointer is not const
        } else if (specs[j].exprs.size() == 2 && h->GetDimension() == 2){
            auto r = d.Histo2D(ROOT::RDF::TH2DModel(modelName.c_str(), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                                                    h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax()), col + "1", col + "0"); // x is the last expression of the Draw string
            src[j] = [r]() mutable {return (TH1*)r->Clone();};
        } else Error("HistMT", "%s: %zu expressions for a %dD histogram", reqs[j].first.c_str(), specs[j].exprs.size(), h->GetDimension());
    }
    for (size_t j = 0; j < reqs.size(); j++){
        TH1* h = hs[j];
        if (!h || !src[j]) continue;
        TH1* s = src[j]();
        s->SetDirectory(nullptr);
        for (int i = 0; i < h->GetNcells(); i++) h->SetBinContent(i, s->GetBinContent(i)); // SetBinContent, not Add: Add would turn Sumw2 on and the plot would show error bars
        double stats[13] = {0};
        s->GetStats(stats);
        h->PutStats(stats);
        h->SetEntries(s->GetEntries());
        delete s;
        t->TAttLine::Copy(*(TAttLine*)h); t->TAttFill::Copy(*(TAttFill*)h); t->TAttMarker::Copy(*(TAttMarker*)h);
    }
    return hs;
}
inline TH1* HistMT(TTree* t, const std::string& drawExpr, const std::string& cut = ""){return HistsMT(t, {{drawExpr, cut}})[0];}
