// PlaceOpt: row-based standard-cell placement (C++17).
//
// Model: cells sit on a grid of rows x sites. A cell occupies `width` contiguous
// sites in one row. Nets connect cell pins (offset from the cell origin) and
// fixed I/O terminals. Objective: total half-perimeter wirelength (HPWL).
#pragma once
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace placeopt {

struct Pin {
    int owner;      // >= 0: cell index; < 0: terminal index (-1 - t)
    double dx, dy;  // offset from cell origin (um); ignored for terminals
};

struct Netlist {
    int rows = 0, sites = 0;          // placement grid
    double site_w = 0.46, row_h = 2.72;
    std::vector<int> width;           // cell width in sites
    std::vector<double> term_x, term_y;
    std::vector<std::vector<Pin>> nets;

    // derived
    std::vector<std::vector<int>> cell_nets;
    void build();
    int num_cells() const { return static_cast<int>(width.size()); }
};

struct Placement {
    std::vector<int> row, site;       // per cell: row index and left site
};

struct AnnealParams {
    uint64_t seed = 1;
    double moves_per_cell = 20.0;     // moves per temperature = moves_per_cell * N^(4/3) / N (VPR style)
    double init_temp = -1.0;          // < 0: auto (20 x std-dev of random-move deltas)
    double exit_ratio = 0.005;        // stop when T < exit_ratio * cost / nets
    long long max_moves = 50'000'000;
    bool start_window_full = true;    // initial range limit = full die
    double init_window = 1.0;         // fraction of die if !start_window_full
    int snapshots = 0;                // number of placement snapshots to keep
};

struct TempStep {
    double temp, cost, accept_rate, window;
    long long moves;
};

struct AnnealResult {
    Placement placement;
    double hpwl = 0.0;
    long long moves = 0;
    double seconds = 0.0;
    std::vector<TempStep> trace;
    std::vector<Placement> snapshots;
};

double net_hpwl(const Netlist& nl, const Placement& p, int net);
double total_hpwl(const Netlist& nl, const Placement& p);
bool is_legal(const Netlist& nl, const Placement& p, std::string* why = nullptr);

Placement random_placement(const Netlist& nl, uint64_t seed);
AnnealResult anneal(const Netlist& nl, const Placement& init, const AnnealParams& prm);

// Analytical global placement: quadratic (clique model) solved with
// Jacobi-preconditioned conjugate gradient, plus iterative cell shifting with
// anchor pseudo-nets to spread cells (FastPlace-style).
struct GlobalResult {
    std::vector<double> x, y;     // cell centers (um)
    std::vector<double> overflow; // density overflow per iteration
    double seconds = 0.0;
};
GlobalResult global_place(const Netlist& nl, int iterations = 30, double target_density = 1.0);

// Greedy legalization: cells in left-edge order each take the nearest free run
// of sites (Manhattan displacement) over all rows.
Placement legalize(const Netlist& nl, const std::vector<double>& x, const std::vector<double>& y);

}  // namespace placeopt
