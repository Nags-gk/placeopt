#include "placeopt.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace placeopt {

namespace {
using Clock = std::chrono::steady_clock;

inline double pin_x(const Netlist& nl, const Placement& p, const Pin& q) {
    return q.owner >= 0 ? p.site[q.owner] * nl.site_w + q.dx : nl.term_x[-1 - q.owner];
}
inline double pin_y(const Netlist& nl, const Placement& p, const Pin& q) {
    return q.owner >= 0 ? p.row[q.owner] * nl.row_h + q.dy : nl.term_y[-1 - q.owner];
}

// Row occupancy grid: occ[r * sites + s] = cell id or -1.
struct Grid {
    int rows, sites;
    std::vector<int> occ;
    Grid(const Netlist& nl, const Placement& p) : rows(nl.rows), sites(nl.sites), occ(size_t(nl.rows) * nl.sites, -1) {
        for (int c = 0; c < nl.num_cells(); ++c)
            for (int k = 0; k < nl.width[c]; ++k) occ[size_t(p.row[c]) * sites + p.site[c] + k] = c;
    }
    int at(int r, int s) const { return occ[size_t(r) * sites + s]; }
    // true if [s, s+w) in row r holds only cells in {a, b} or is empty
    bool free_except(int r, int s, int w, int a, int b) const {
        if (s < 0 || s + w > sites) return false;
        for (int k = 0; k < w; ++k) {
            int o = at(r, s + k);
            if (o != -1 && o != a && o != b) return false;
        }
        return true;
    }
    void fill(int r, int s, int w, int v) {
        for (int k = 0; k < w; ++k) occ[size_t(r) * sites + s + k] = v;
    }
};
}  // namespace

void Netlist::build() {
    cell_nets.assign(width.size(), {});
    for (int n = 0; n < static_cast<int>(nets.size()); ++n)
        for (const Pin& q : nets[n])
            if (q.owner >= 0) {
                auto& v = cell_nets[q.owner];
                if (v.empty() || v.back() != n) v.push_back(n);
            }
    long long need = std::accumulate(width.begin(), width.end(), 0LL);
    if (need > 1LL * rows * sites) throw std::invalid_argument("cells do not fit in the die");
}

double net_hpwl(const Netlist& nl, const Placement& p, int net) {
    double x0 = std::numeric_limits<double>::max(), x1 = -x0, y0 = x0, y1 = -x0;
    for (const Pin& q : nl.nets[net]) {
        double x = pin_x(nl, p, q), y = pin_y(nl, p, q);
        x0 = std::min(x0, x); x1 = std::max(x1, x);
        y0 = std::min(y0, y); y1 = std::max(y1, y);
    }
    return (x1 - x0) + (y1 - y0);
}

double total_hpwl(const Netlist& nl, const Placement& p) {
    double s = 0;
    for (int n = 0; n < static_cast<int>(nl.nets.size()); ++n) s += net_hpwl(nl, p, n);
    return s;
}

bool is_legal(const Netlist& nl, const Placement& p, std::string* why) {
    std::vector<int> occ(size_t(nl.rows) * nl.sites, -1);
    for (int c = 0; c < nl.num_cells(); ++c) {
        if (p.row[c] < 0 || p.row[c] >= nl.rows || p.site[c] < 0 || p.site[c] + nl.width[c] > nl.sites) {
            if (why) *why = "cell " + std::to_string(c) + " out of die";
            return false;
        }
        for (int k = 0; k < nl.width[c]; ++k) {
            int& o = occ[size_t(p.row[c]) * nl.sites + p.site[c] + k];
            if (o != -1) {
                if (why) *why = "cells " + std::to_string(o) + " and " + std::to_string(c) + " overlap";
                return false;
            }
            o = c;
        }
    }
    return true;
}

Placement random_placement(const Netlist& nl, uint64_t seed) {
    // Shuffle cells, then pack them into rows chosen at random among those with room,
    // leaving random gaps so free space is spread across the die.
    std::mt19937_64 rng(seed);
    const int n = nl.num_cells();
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);
    long long used = std::accumulate(nl.width.begin(), nl.width.end(), 0LL);
    double slack = 1.0 - double(used) / (double(nl.rows) * nl.sites);
    std::vector<int> frontier(nl.rows, 0);
    Placement p{std::vector<int>(n), std::vector<int>(n)};
    std::uniform_real_distribution<double> u(0, 1);
    for (int c : order) {
        int w = nl.width[c];
        std::vector<int> ok;
        for (int r = 0; r < nl.rows; ++r)
            if (frontier[r] + w <= nl.sites) ok.push_back(r);
        if (ok.empty()) throw std::runtime_error("random placement failed: die too full");
        int r = ok[rng() % ok.size()];
        int gap = (u(rng) < slack) ? 1 : 0;
        if (frontier[r] + gap + w > nl.sites) gap = 0;
        p.row[c] = r;
        p.site[c] = frontier[r] + gap;
        frontier[r] += gap + w;
    }
    return p;
}

AnnealResult anneal(const Netlist& nl, const Placement& init, const AnnealParams& prm) {
    auto t_start = Clock::now();
    const int n = nl.num_cells();
    std::mt19937_64 rng(prm.seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    Placement p = init;
    Grid g(nl, p);

    std::vector<double> nh(nl.nets.size());
    double cost = 0;
    for (size_t k = 0; k < nl.nets.size(); ++k) cost += (nh[k] = net_hpwl(nl, p, int(k)));

    std::vector<int> touched_mark(nl.nets.size(), -1), touched;
    int stamp = 0;
    std::vector<double> new_h;

    // Evaluate the HPWL delta of moving cells (a -> ra,sa) and optionally (b -> rb,sb).
    auto delta_for = [&](int a, int ra, int sa, int b, int rb, int sb) {
        ++stamp;
        touched.clear();
        for (int net : nl.cell_nets[a]) if (touched_mark[net] != stamp) { touched_mark[net] = stamp; touched.push_back(net); }
        if (b >= 0) for (int net : nl.cell_nets[b]) if (touched_mark[net] != stamp) { touched_mark[net] = stamp; touched.push_back(net); }
        int oa_r = p.row[a], oa_s = p.site[a], ob_r = b >= 0 ? p.row[b] : 0, ob_s = b >= 0 ? p.site[b] : 0;
        p.row[a] = ra; p.site[a] = sa;
        if (b >= 0) { p.row[b] = rb; p.site[b] = sb; }
        double d = 0;
        new_h.resize(touched.size());
        for (size_t i = 0; i < touched.size(); ++i) {
            new_h[i] = net_hpwl(nl, p, touched[i]);
            d += new_h[i] - nh[touched[i]];
        }
        p.row[a] = oa_r; p.site[a] = oa_s;
        if (b >= 0) { p.row[b] = ob_r; p.site[b] = ob_s; }
        return d;
    };

    // Propose a move. Returns false if no legal move was found this try.
    struct Move { int a, ra, sa, b, rb, sb; };
    double win = prm.start_window_full ? 1.0 : prm.init_window;  // fraction of die
    auto propose = [&](Move& m) {
        int a = int(rng() % n);
        // Range-limited move; 5% of moves ignore the window so the search can still
        // escape a frozen region late in the anneal.
        double w = u01(rng) < 0.05 ? 1.0 : win;
        int wr = std::max(1, int(std::lround(w * nl.rows))), ws = std::max(3, int(std::lround(w * nl.sites)));
        int r = std::clamp(p.row[a] + int(rng() % (2 * wr + 1)) - wr, 0, nl.rows - 1);
        int s = std::clamp(p.site[a] + int(rng() % (2 * ws + 1)) - ws, 0, nl.sites - nl.width[a]);
        int b = g.at(r, s);
        if (b == a) return false;
        if (b == -1 || u01(rng) < 0.3) {  // displacement into (r, s) if free
            if (!g.free_except(r, s, nl.width[a], a, -1)) return false;
            m = {a, r, s, -1, 0, 0};
            return true;
        }
        // swap a with b: a takes b's left edge, b takes a's left edge
        int ra = p.row[b], sa = p.site[b], rb = p.row[a], sb = p.site[a];
        // shift a's target left if it would overrun the row end
        sa = std::min(sa, nl.sites - nl.width[a]);
        sb = std::min(sb, nl.sites - nl.width[b]);
        // check fits ignoring a and b themselves, and that a and b do not overlap each other
        if (!g.free_except(ra, sa, nl.width[a], a, b) || !g.free_except(rb, sb, nl.width[b], a, b)) return false;
        if (ra == rb && !(sa + nl.width[a] <= sb || sb + nl.width[b] <= sa)) return false;
        m = {a, ra, sa, b, rb, sb};
        return true;
    };
    auto apply = [&](const Move& m) {
        g.fill(p.row[m.a], p.site[m.a], nl.width[m.a], -1);
        if (m.b >= 0) g.fill(p.row[m.b], p.site[m.b], nl.width[m.b], -1);
        p.row[m.a] = m.ra; p.site[m.a] = m.sa;
        g.fill(m.ra, m.sa, nl.width[m.a], m.a);
        if (m.b >= 0) { p.row[m.b] = m.rb; p.site[m.b] = m.sb; g.fill(m.rb, m.sb, nl.width[m.b], m.b); }
        for (size_t i = 0; i < touched.size(); ++i) nh[touched[i]] = new_h[i];
    };

    // Initial temperature: 20 x std-dev of deltas over N random moves (VPR).
    double T = prm.init_temp;
    if (T < 0) {
        std::vector<double> ds;
        for (int i = 0, tries = 0; i < n && tries < 50 * n; ++tries) {
            Move m;
            if (!propose(m)) continue;
            ds.push_back(delta_for(m.a, m.ra, m.sa, m.b, m.rb, m.sb));
            ++i;
        }
        double mu = std::accumulate(ds.begin(), ds.end(), 0.0) / std::max<size_t>(1, ds.size()), var = 0;
        for (double d : ds) var += (d - mu) * (d - mu);
        T = 20.0 * std::sqrt(var / std::max<size_t>(1, ds.size())) + 1e-9;
    }

    // Keep the best placement seen: a late uphill move can otherwise freeze in.
    Placement best = p;
    double best_cost = cost;
    AnnealResult res;
    const long long moves_per_T = std::max<long long>(1, (long long)(prm.moves_per_cell * std::pow(n, 4.0 / 3.0)));
    long long total = 0;
    const double nets = std::max<size_t>(1, nl.nets.size());
    int snap_every_steps = 0;
    while (true) {
        long long acc = 0, tried = 0;
        for (long long k = 0; k < moves_per_T && total < prm.max_moves; ++k, ++total) {
            Move m;
            if (!propose(m)) continue;
            ++tried;
            double d = delta_for(m.a, m.ra, m.sa, m.b, m.rb, m.sb);
            if (d <= 0 || u01(rng) < std::exp(-d / T)) {
                apply(m);
                cost += d;
                ++acc;
                if (cost < best_cost - 1e-9) { best_cost = cost; best = p; }
            }
        }
        double ar = tried ? double(acc) / tried : 0.0;
        res.trace.push_back({T, cost, ar, win, total});
        // VPR adaptive schedule: cool slowly in the productive acceptance range.
        double alpha = ar > 0.96 ? 0.5 : ar > 0.8 ? 0.9 : ar > 0.15 ? 0.95 : 0.8;
        T *= alpha;
        // Range limit tracks a 0.44 acceptance target.
        win = std::clamp(win * (1.0 - 0.44 + ar), 1.0 / std::max(nl.rows, nl.sites), 1.0);
        if (prm.snapshots > 0 && res.trace.size() % 8 == 1 && int(res.snapshots.size()) < prm.snapshots)
            res.snapshots.push_back(p);
        (void)snap_every_steps;
        if (T < prm.exit_ratio * cost / nets || total >= prm.max_moves) break;
    }
    // Restart the final greedy (T = 0) pass from the best placement seen.
    if (best_cost < cost - 1e-9) {
        p = best;
        g = Grid(nl, p);
        cost = 0;
        for (size_t k = 0; k < nl.nets.size(); ++k) cost += (nh[k] = net_hpwl(nl, p, int(k)));
    }
    for (long long k = 0; k < moves_per_T && total < prm.max_moves; ++k, ++total) {
        Move m;
        if (!propose(m)) continue;
        double d = delta_for(m.a, m.ra, m.sa, m.b, m.rb, m.sb);
        if (d < -1e-12) { apply(m); cost += d; }
    }
    res.placement = p;
    res.hpwl = total_hpwl(nl, p);  // recomputed from scratch (guards against drift)
    res.moves = total;
    res.seconds = std::chrono::duration<double>(Clock::now() - t_start).count();
    if (prm.snapshots > 0) res.snapshots.push_back(p);
    return res;
}

// ---------------------------------------------------------------- analytical

namespace {
// Sparse symmetric matrix in coordinate-accumulated row form.
struct SpMat {
    int n;
    std::vector<std::vector<std::pair<int, double>>> rows;
    std::vector<double> diag;
    explicit SpMat(int n_) : n(n_), rows(n_), diag(n_, 0.0) {}
    void add_edge(int i, int j, double w) {  // between movable i and j
        diag[i] += w; diag[j] += w;
        rows[i].push_back({j, -w}); rows[j].push_back({i, -w});
    }
    void mul(const std::vector<double>& x, std::vector<double>& y) const {
        for (int i = 0; i < n; ++i) {
            double s = diag[i] * x[i];
            for (auto& [j, v] : rows[i]) s += v * x[j];
            y[i] = s;
        }
    }
};

// Jacobi-preconditioned conjugate gradient.
void pcg(const SpMat& A, const std::vector<double>& b, std::vector<double>& x, int iters = 500, double tol = 1e-6) {
    int n = A.n;
    std::vector<double> r(n), z(n), p(n), Ap(n);
    A.mul(x, Ap);
    for (int i = 0; i < n; ++i) { r[i] = b[i] - Ap[i]; z[i] = r[i] / A.diag[i]; p[i] = z[i]; }
    double rz = 0, bn = 0;
    for (int i = 0; i < n; ++i) { rz += r[i] * z[i]; bn += b[i] * b[i]; }
    bn = std::sqrt(bn) + 1e-12;
    for (int it = 0; it < iters; ++it) {
        A.mul(p, Ap);
        double pAp = 0;
        for (int i = 0; i < n; ++i) pAp += p[i] * Ap[i];
        double alpha = rz / (pAp + 1e-30), rn = 0;
        for (int i = 0; i < n; ++i) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; rn += r[i] * r[i]; }
        if (std::sqrt(rn) / bn < tol) break;
        double rz2 = 0;
        for (int i = 0; i < n; ++i) { z[i] = r[i] / A.diag[i]; rz2 += r[i] * z[i]; }
        double beta = rz2 / rz;
        rz = rz2;
        for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
    }
}
}  // namespace

GlobalResult global_place(const Netlist& nl, int iterations, double target_density) {
    auto t0 = Clock::now();
    const int n = nl.num_cells();
    const double W = nl.sites * nl.site_w, H = nl.rows * nl.row_h;
    GlobalResult res;
    res.x.assign(n, W / 2);
    res.y.assign(n, H / 2);
    std::vector<double> ax(n, W / 2), ay(n, H / 2);  // anchor targets from cell shifting
    double anchor_w = 0.0;
    std::vector<double> cw(n);
    for (int c = 0; c < n; ++c) cw[c] = nl.width[c] * nl.site_w;

    for (int it = 0; it <= iterations; ++it) {
        // Build x and y systems (clique net model, weight 1/(k-1)). Pin offsets are
        // ignored here; legalization and detailed placement restore pin accuracy.
        SpMat A(n);
        std::vector<double> bx(n, 0.0), by(n, 0.0);
        for (const auto& net : nl.nets) {
            int k = int(net.size());
            if (k < 2) continue;
            double w = 1.0 / (k - 1);
            for (int i = 0; i < k; ++i)
                for (int j = i + 1; j < k; ++j) {
                    const Pin &a = net[i], &b = net[j];
                    if (a.owner >= 0 && b.owner >= 0) {
                        if (a.owner != b.owner) A.add_edge(a.owner, b.owner, w);
                    } else if (a.owner >= 0 || b.owner >= 0) {
                        int c = a.owner >= 0 ? a.owner : b.owner;
                        int t = -1 - (a.owner >= 0 ? b.owner : a.owner);
                        A.diag[c] += w; bx[c] += w * nl.term_x[t]; by[c] += w * nl.term_y[t];
                    }
                }
        }
        for (int c = 0; c < n; ++c) {
            if (A.diag[c] == 0) A.diag[c] = 1e-6;  // isolated cell
            if (anchor_w > 0) { A.diag[c] += anchor_w; bx[c] += anchor_w * ax[c]; by[c] += anchor_w * ay[c]; }
        }
        pcg(A, bx, res.x);
        pcg(A, by, res.y);
        if (it == iterations) break;

        // Cell shifting (FastPlace): bin the die, then remap each bin row/column's
        // boundaries in proportion to utilization so dense bins expand.
        int nb = std::max(2, int(std::sqrt(n / 4.0)));
        double bw = W / nb, bh = H / nb;
        std::vector<double> util(size_t(nb) * nb, 0.0);
        for (int c = 0; c < n; ++c) {
            int i = std::clamp(int(res.x[c] / bw), 0, nb - 1), j = std::clamp(int(res.y[c] / bh), 0, nb - 1);
            util[size_t(j) * nb + i] += cw[c] * nl.row_h;
        }
        double cap = bw * bh * target_density, over = 0;
        for (double u : util) over += std::max(0.0, u - cap);
        res.overflow.push_back(over / std::max(1e-9, std::accumulate(cw.begin(), cw.end(), 0.0) * nl.row_h));

        auto shift = [&](bool horizontal) {
            for (int line = 0; line < nb; ++line) {
                // new bin boundaries along this row (horizontal) or column
                std::vector<double> u(nb);
                for (int k = 0; k < nb; ++k) u[k] = (horizontal ? util[size_t(line) * nb + k] : util[size_t(k) * nb + line]) / cap;
                const double L = horizontal ? bw : bh;
                std::vector<double> nbnd(nb + 1);
                nbnd[0] = 0; nbnd[nb] = nb * L;
                for (int k = 1; k < nb; ++k) {
                    double a = 1.0 + u[k - 1], b = 1.0 + u[k];  // FastPlace: delta = 1.5 (rounded to 1)
                    nbnd[k] = ((k - 1) * L * b + (k + 1) * L * a) / (a + b);
                }
                for (int c = 0; c < n; ++c) {
                    double& v = horizontal ? ax[c] : ay[c];
                    double pos = horizontal ? res.x[c] : res.y[c];
                    double other = horizontal ? res.y[c] : res.x[c];
                    int ol = std::clamp(int(other / (horizontal ? bh : bw)), 0, nb - 1);
                    if (ol != line) continue;
                    int k = std::clamp(int(pos / L), 0, nb - 1);
                    double f = (pos - k * L) / L;
                    v = nbnd[k] + f * (nbnd[k + 1] - nbnd[k]);
                }
            }
        };
        shift(true);
        shift(false);
        anchor_w = 0.01 * (it + 1) * (it + 1);  // stiffen anchors progressively
    }
    res.seconds = std::chrono::duration<double>(Clock::now() - t0).count();
    return res;
}

Placement legalize(const Netlist& nl, const std::vector<double>& x, const std::vector<double>& y) {
    // Greedy nearest-slot legalization in left-edge order (Tetris-like): each cell
    // takes the free run of sites, over all rows, closest to its global position.
    const int n = nl.num_cells();
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::vector<double> left(n);
    for (int c = 0; c < n; ++c) left[c] = x[c] - nl.width[c] * nl.site_w / 2;
    std::sort(order.begin(), order.end(), [&](int a, int b) { return left[a] < left[b]; });
    std::vector<char> occ(size_t(nl.rows) * nl.sites, 0);
    auto free_run = [&](int r, int s, int w) {
        if (s < 0 || s + w > nl.sites) return false;
        for (int k = 0; k < w; ++k) if (occ[size_t(r) * nl.sites + s + k]) return false;
        return true;
    };
    Placement p{std::vector<int>(n), std::vector<int>(n)};
    for (int c : order) {
        const int w = nl.width[c];
        const int want = std::clamp(int(std::lround(left[c] / nl.site_w)), 0, nl.sites - w);
        double best = std::numeric_limits<double>::max();
        int br = -1, bs = -1;
        for (int r = 0; r < nl.rows; ++r) {
            double dy = std::abs((r + 0.5) * nl.row_h - y[c]);
            if (dy >= best) continue;
            // nearest free run in this row, searching outward from `want`
            for (int d = 0; d < nl.sites; ++d) {
                double cost = d * nl.site_w + dy;
                if (cost >= best) break;
                if (free_run(r, want - d, w)) { best = cost; br = r; bs = want - d; break; }
                if (d && free_run(r, want + d, w)) { best = cost; br = r; bs = want + d; break; }
            }
        }
        if (br < 0) throw std::runtime_error("legalization failed: no free slot");
        p.row[c] = br; p.site[c] = bs;
        for (int k = 0; k < w; ++k) occ[size_t(br) * nl.sites + bs + k] = 1;
    }
    return p;
}

}  // namespace placeopt
