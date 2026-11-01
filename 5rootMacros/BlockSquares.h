#pragma once
// Block position maps (one bin per detector or tracker cell, bins on the cell grid) drawn as one filled square per occupied cell, shared by Plot3D.cpp (1_position/2_block),
// GraphBlocks.cpp (4_CEblock, the *_xy_detector plots) and PlotDiffPhiCuts.cpp (the block xy plots).
#include <TAxis.h>
#include <TColor.h>
#include <TH1.h>
#include <TPolyLine.h>

#include <algorithm>
#include <cmath>

// Call right after h->Draw("colz") on the pad that shows h. Each occupied bin becomes a TPolyLine filled with the color the colz bin would have (log or linear z, current palette).
// A TBox is not used: ROOT writes it to the pdf as a fill plus a 0.4 pt stroke ("re b*"), which grows the square past its cell and covers its neighbors; a filled polygon has no stroke.
// The edges come from one formula of the integer bin index, so a square and its neighbor share the identical edge value (GetBinUpEdge(i) and GetBinLowEdge(i + 1) differ in the last
// digits, and the pdf rounds each corner). The histogram then keeps only the axes, colorbar and stats: its own colz squares are emptied (stats are put back) so they do not show through.
// halfRange > 0: hard frame of +-halfRange, the axis is relabeled to it and the squares at the edge are clipped by it (h then has whole cells past +-halfRange).
inline void DrawBlockSquares(TH1* h, double halfRange = 0, bool logz = true){
    const double hr = halfRange > 0 ? halfRange : 1e9;
    const double zMax = std::max(2.0, h->GetMaximum());
    const double xLo = h->GetXaxis()->GetXmin(), xW = h->GetXaxis()->GetBinWidth(1), yLo = h->GetYaxis()->GetXmin(), yW = h->GetYaxis()->GetBinWidth(1);
    for (int ix = 1; ix <= h->GetNbinsX(); ix++) for (int iy = 1; iy <= h->GetNbinsY(); iy++){
        const double c = h->GetBinContent(ix, iy);
        if (c < (logz ? 1.0 : 1e-12)) continue;
        const double x1 = std::max(-hr, xLo + (ix - 1) * xW), x2 = std::min(hr, xLo + ix * xW);
        const double y1 = std::max(-hr, yLo + (iy - 1) * yW), y2 = std::min(hr, yLo + iy * yW);
        if (x1 >= x2 || y1 >= y2) continue;
        double bx[5] = {x1, x2, x2, x1, x1}, by[5] = {y1, y1, y2, y2, y1};
        TPolyLine* b = new TPolyLine(5, bx, by);
        const double f = logz ? std::log(c) / std::log(zMax) : c / zMax;
        b->SetFillColor(TColor::GetColorPalette(int(254 * std::min(1.0, std::max(0.0, f))))); b->SetFillStyle(1001);
        b->Draw("f");
    }
    double st[13]; h->GetStats(st); const double ent = h->GetEntries();
    for (int ix = 1; ix <= h->GetNbinsX(); ix++) for (int iy = 1; iy <= h->GetNbinsY(); iy++) h->SetBinContent(ix, iy, 0);
    if (halfRange > 0){h->GetXaxis()->SetLimits(-halfRange, halfRange); h->GetYaxis()->SetLimits(-halfRange, halfRange);}
    h->PutStats(st); h->SetEntries(ent);
    h->SetMinimum(logz ? 1 : 0); h->SetMaximum(zMax);
}
