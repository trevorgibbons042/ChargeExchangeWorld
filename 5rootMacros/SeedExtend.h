#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

// Seed and extend track finding on the tracker hits of one event (block positions x, y, z, any particle).
// planes: every tracker plane the event has a hit on (z in um, ascending, lround(z * 1000)), whatever particle, also the hits no track holds any more. A track runs through
// consecutive planes of that list: it cannot jump over a plane, since a charged particle leaves a hit on every plane it crosses.
// A seed is a pair of hits, one on the lowest plane that still has a hit and one on the next plane of the list (a plane with no hit left there gives no seed; the hits left alone below
// the seed plane are not tracks). Its line is followed plane by plane: on each next plane the hit nearest the line inside a road is taken (3 sigma grown by the extrapolation of
// the line fitted so far, added in quadrature to a scattering allowance of kSeedScatterRad times the distance to the last hit), the line is refitted (least squares, all hits so far) after every hit that is added, and the track stops on the first plane with no hit in the road. The candidate with the most planes, then
// the smallest summed squared distance in sigma, is the track. Returns the indices into x, y, z of its hits, one per plane, lowest z first; empty means no track.
// Positions only, no track ID. Hits of the same plane and cell are the same point, so which one is taken changes no fit (stable sort, the order the hits come in).
// Multiple scattering allowance of the road (rad): the next hit may sit this angle times the distance to the last hit off the line, on top of the position error. 5 mrad:
// the hits of the kaon and the pion of event 263272 of the xs1 run are 5 to 7 mm off the line of the hits before them at the next tracker, 1400 mm on.
constexpr double kSeedScatterRad = 0.005;

// maxSkip: 0 (default) a track cannot jump over a plane. 1 lets it jump over one plane in a row that has no hit in the road (a tracker with an efficiency below 1, /detector/trackerEfficiency in
// the sim, leaves no hit now and then): the seed pair may then sit on planes 2 apart, and a plane without a hit in the road is passed once, the next plane must have one. A track with a gap
// has fewer hits than the same track without, so the finding still prefers the track with no gap.
inline std::vector<std::size_t> SeedExtend(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<long>& planes, double sigma, int maxSkip = 0){
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < z.size(); i++) idx.push_back(i);
    std::stable_sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b){
        const long pa = std::lround(z[a] * 1000.0), pb = std::lround(z[b] * 1000.0);
        if (pa != pb) return pa < pb;
        if (x[a] != x[b]) return x[a] < x[b];
        return y[a] < y[b];
    });
    // the hits left on each plane of the list
    std::vector<std::vector<std::size_t>> on(planes.size());
    for (std::size_t i : idx) {
        const auto it = std::lower_bound(planes.begin(), planes.end(), std::lround(z[i] * 1000.0));
        if (it != planes.end() && *it == std::lround(z[i] * 1000.0)) on[it - planes.begin()].push_back(i);
    }
    std::vector<std::size_t> best, cur;
    double bestChi2 = INFINITY;
    for (std::size_t p0 = 0; p0 + 1 < planes.size() && best.empty(); p0++) {
        if (on[p0].empty()) continue; // the first plane that has a hit and a next plane with a hit gives the seeds
        for (std::size_t p1 = p0 + 1; p1 <= p0 + 1 + (std::size_t)maxSkip && p1 < planes.size(); p1++) for (const std::size_t a : on[p0]) for (const std::size_t b : on[p1]) {
            if (z[b] == z[a]) continue;
            cur.assign({a, b});
            // the running least squares line x = ax + bx z, y = ay + by z through the hits taken so far, refitted after every hit that is added; the sums are about the first hit
            double n = 2.0, sz = z[b] - z[a], sx = x[b] - x[a], sy = y[b] - y[a], szz = sz * sz, szx = sz * sx, szy = sz * sy; // the first hit is the origin (0, 0, 0)
            // the seed pair counts too: two candidates with as many planes are also told apart by how far apart their seed hits sit, in sigma
            double chi2 = ((x[b] - x[a]) * (x[b] - x[a]) + (y[b] - y[a]) * (y[b] - y[a])) / (sigma * sigma);
            bool skipped = false; // the last plane was passed with no hit
            for (std::size_t pl = p1 + 1; pl < planes.size(); pl++) {
                const double zc = (on[pl].empty() ? planes[pl] / 1000.0 : z[on[pl].front()]) - z[a];
                const double den = n * szz - sz * sz, bx = (n * szx - sz * sx) / den, by = (n * szy - sz * sy) / den;
                const double px = x[a] + (sx - bx * sz) / n + bx * zc, py = y[a] + (sy - by * sz) / n + by * zc;
                // the prediction error of a least squares line at zc, in units of the hit error: 1 + 1/n + (zc - mean z)^2 / sum of (z - mean z)^2
                const double zm = sz / n, spread = szz - n * zm * zm;
                const double posErr = 3.0 * sigma * std::sqrt(1.0 + 1.0 / n + (zc - zm) * (zc - zm) / spread);
                const double scatter = kSeedScatterRad * (zc - (z[cur.back()] - z[a])); // the track bends in the material, the line is off by an angle times the distance from the last hit
                const double road = std::sqrt(posErr * posErr + scatter * scatter);
                std::size_t pick = z.size(); double d2Best = road * road;
                for (const std::size_t h : on[pl]) {
                    const double d2 = (x[h] - px) * (x[h] - px) + (y[h] - py) * (y[h] - py);
                    if (d2 < d2Best) {d2Best = d2; pick = h;}
                }
                if (pick == z.size()) { // nothing in the road on this plane: the track ends here, or passes it once when skipping is on
                    if (maxSkip > 0 && !skipped) {skipped = true; continue;}
                    break;
                }
                skipped = false;
                cur.push_back(pick); chi2 += d2Best / (sigma * sigma);
                const double dz = z[pick] - z[a], dx = x[pick] - x[a], dy = y[pick] - y[a];
                n += 1.0; sz += dz; sx += dx; sy += dy; szz += dz * dz; szx += dz * dx; szy += dz * dy;
            }
            if (cur.size() > best.size() || (cur.size() == best.size() && chi2 < bestChi2)) {best = cur; bestChi2 = chi2;}
        }
    }
    return best;
}
