#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>

#include "CutFiles.h"
#include "DrawPlan.h"

// T2 (spec F): DrawPlan requests against TTree::Draw on the same tree: class, name, bins, all bin contents including under/overflow, GetEntries, GetMean/GetStdDev bit-equal (%a),
// line/fill/marker attributes, Cols element-wise, N, a model histogram, a weighted request, a listCut request. Needs a run folder with cuts/cut_APandT.root and cuts/ke_all.root.
// Run: root -l -b -q TestDrawPlan.cpp   (DRAWPLAN_FALLBACK=1 forces the plain Draw path, which must then match trivially)
static int gBad = 0;
static void Check(bool ok, const char* what){if (!ok){printf("  MISMATCH: %s\n", what); gBad++;}}
static bool SameHist(TH1* a, TH1* b, const char* tag){
    printf("%s: class %s/%s name %s/%s entries %.0f/%.0f mean %a/%a sd %a/%a\n", tag, a->ClassName(), b->ClassName(), a->GetName(), b->GetName(), a->GetEntries(), b->GetEntries(), a->GetMean(), b->GetMean(), a->GetStdDev(), b->GetStdDev());
    Check(std::string(a->ClassName()) == b->ClassName(), "class"); Check(std::string(a->GetName()).substr(1) == std::string(b->GetName()).substr(1), "name (first letter differs by design: h = DrawPlan, r = Draw)");
    Check(a->GetNbinsX() == b->GetNbinsX() && a->GetNbinsY() == b->GetNbinsY(), "nbins");
    Check(a->GetXaxis()->GetXmin() == b->GetXaxis()->GetXmin() && a->GetXaxis()->GetXmax() == b->GetXaxis()->GetXmax(), "x edges");
    bool same = true;
    for (int i = 0; i < a->GetNcells(); i++) if (a->GetBinContent(i) != b->GetBinContent(i)) same = false;
    Check(same, "bin contents incl. under/overflow");
    Check(a->GetEntries() == b->GetEntries(), "entries"); Check(a->GetMean() == b->GetMean(), "mean bit-equal"); Check(a->GetStdDev() == b->GetStdDev(), "stddev bit-equal");
    Check(a->GetLineColor() == b->GetLineColor() && a->GetLineWidth() == b->GetLineWidth() && a->GetFillColor() == b->GetFillColor() && a->GetMarkerStyle() == b->GetMarkerStyle(), "attributes");
    if (std::string(tag).find("listCut") == std::string::npos) Check(std::string(a->GetTitle()) == b->GetTitle(), "title"); // a listCut is not part of the title: the tree's list is not either
    return same;
}
void TestDrawPlan(){
    ROOT::EnableImplicitMT(NThreads());
    TFile* f = TFile::Open("cuts/cut_APandT.root", "READ");
    TTree* t = f ? (TTree*)f->Get("CalcValsReD_DE") : nullptr;
    if (!t){printf("needs cuts/cut_APandT.root, tree CalcValsReD_DE\\n"); return;}
    gROOT->cd();
    // find a few numeric branches
    std::vector<std::string> br;
    for (auto* l : *t->GetListOfLeaves()) br.push_back(l->GetName());
    printf("tree %s: %lld entries, %zu leaves, first: %s %s %s\n", t->GetName(), t->GetEntries(), br.size(), br[0].c_str(), br.size() > 1 ? br[1].c_str() : "", br.size() > 2 ? br[2].c_str() : "");
    const std::string a = "finalKE", b = "phiReD", c = "thetaReD";
    DrawPlan plan;
    // booked requests
    auto h1 = plan.H(t, a + ">>hA(100,0,30000)");
    auto h2 = plan.H(t, a + ">>hB(100,0,30000)", c + " > 0.5 && " + c + " < 4");
    auto h3 = plan.H(t, b + ":" + c + ">>hC(40,0,5,36,-180,180)", c + " < 5");
    TH1D model("hM", "", 50, 0, 360); model.Sumw2();
    auto h4 = plan.H(t, b + ">>hM", c + " > 0.2", "", &model);
    auto h5 = plan.H(t, a + ">>hW(100,0,30000)", "(" + c + " > 0.2)*0.5", "", nullptr, true);
    auto h6 = plan.H(t, a + ">>hL(100,0,30000)", "", kBlockCut); // listCut without a list set on the tree
    auto h7 = plan.H(t, a + ">>hE(100,0,30000)", c + " < 4", ListCut("ke_all", true)); // the tree has the ke_all block list set below (CutTrees.cpp)
    auto n1 = plan.N(t, c + " > 1");
    auto cs = plan.C(t, {a, b}, c + " > 1");
    ApplyCut(t, "ke_all", true); // the list stays on the tree while DrawPlan runs: Run() must ignore it, ListCut stands for it
    plan.Run();
    t->SetEntryList(nullptr);
    ApplyCut(t, "ke_all", true);
    t->Draw((a + ">>rE(100,0,30000)").c_str(), (c + " < 4").c_str(), "goff"); TH1* r7 = (TH1*)gDirectory->Get("rE");
    SameHist(h7.Get(), r7, "listCut vs the tree's entry list");
    printf("tree list restored after Run: %d\n", t->GetEntryList() != nullptr);
    t->SetEntryList(nullptr);
    // reference with TTree::Draw
    t->Draw((a + ">>rA(100,0,30000)").c_str(), "", "goff"); TH1* r1 = (TH1*)gDirectory->Get("rA");
    t->Draw((a + ">>rB(100,0,30000)").c_str(), (c + " > 0.5 && " + c + " < 4").c_str(), "goff"); TH1* r2 = (TH1*)gDirectory->Get("rB");
    t->Draw((b + ":" + c + ">>rC(40,0,5,36,-180,180)").c_str(), (c + " < 5").c_str(), "goff"); TH1* r3 = (TH1*)gDirectory->Get("rC");
    TH1D model2("hM2", "", 50, 0, 360); model2.Sumw2();
    t->Draw((b + ">>hM2").c_str(), (c + " > 0.2").c_str(), "goff");
    t->Draw((a + ">>rW(100,0,30000)").c_str(), ("(" + c + " > 0.2)*0.5").c_str(), "goff"); TH1* r5 = (TH1*)gDirectory->Get("rW");
    t->Draw((a + ">>rL(100,0,30000)").c_str(), kBlockCut.c_str(), "goff"); TH1* r6 = (TH1*)gDirectory->Get("rL");
    SameHist(h1.Get(), r1, "1D no cut"); SameHist(h2.Get(), r2, "1D cut"); SameHist(h3.Get(), r3, "2D"); SameHist(h4.Get(), &model2, "model");
    SameHist(h5.Get(), r5, "weighted (selection is the weight)");
    SameHist(h6.Get(), r6, "listCut as selection");
    printf("N: %lld vs %lld\n", n1.Get(), t->GetEntries((c + " > 1").c_str())); Check(n1.Get() == t->GetEntries((c + " > 1").c_str()), "N");
    t->SetEstimate(t->GetEntries() + 1);
    const Long64_t nd = t->Draw((a + ":" + b).c_str(), (c + " > 1").c_str(), "goff"); // v1 = a, v2 = b
    Check(nd == cs.N(), "cols n");
    bool eq = nd == cs.N();
    for (Long64_t i = 0; eq && i < nd; i++) if (t->GetV1()[i] != cs.V(0)[i] || t->GetV2()[i] != cs.V(1)[i]) eq = false;
    Check(eq, "cols element-wise in order"); printf("cols %lld rows\n", nd);
    printf(gBad ? "TestDrawPlan: %d MISMATCHES\n" : "TestDrawPlan: all identical\n", gBad);
}
