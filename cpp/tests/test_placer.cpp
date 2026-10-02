// Standalone C++ tests (registered with CTest).
#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>

#include "placeopt.hpp"

using namespace placeopt;

static Netlist chain(int n, int rows, int sites) {
    // n cells in a chain between two terminals on opposite die edges.
    Netlist nl;
    nl.rows = rows; nl.sites = sites;
    nl.width.assign(n, 3);
    nl.term_x = {0.0, sites * nl.site_w};
    nl.term_y = {rows * nl.row_h / 2, rows * nl.row_h / 2};
    nl.nets.push_back({{-1, 0, 0}, {0, 0.5, 1.0}});
    for (int i = 0; i + 1 < n; ++i) nl.nets.push_back({{i, 1.0, 1.0}, {i + 1, 0.3, 1.0}});
    nl.nets.push_back({{n - 1, 1.0, 1.0}, {-2, 0, 0}});
    nl.build();
    return nl;
}

int main() {
    // HPWL of a hand-placed 2-pin net
    {
        Netlist nl;
        nl.rows = 2; nl.sites = 10; nl.width = {2, 2};
        nl.nets = {{{0, 0.1, 0.2}, {1, 0.3, 0.4}}};
        nl.build();
        Placement p{{0, 1}, {0, 5}};
        double want = std::abs((5 * 0.46 + 0.3) - 0.1) + std::abs((2.72 + 0.4) - 0.2);
        assert(std::abs(total_hpwl(nl, p) - want) < 1e-9);
    }
    // Annealing a chain: legal and much shorter than random
    {
        Netlist nl = chain(30, 6, 40);
        Placement r = random_placement(nl, 7);
        assert(is_legal(nl, r));
        AnnealParams prm; prm.seed = 3;
        AnnealResult a = anneal(nl, r, prm);
        std::string why;
        assert(is_legal(nl, a.placement, &why));
        assert(std::abs(a.hpwl - total_hpwl(nl, a.placement)) < 1e-6);
        double lower = 40 * nl.site_w;   // chain must span the die horizontally
        std::printf("random %.1f  annealed %.1f  lower bound %.1f\n", total_hpwl(nl, r), a.hpwl, lower);
        assert(a.hpwl < 0.5 * total_hpwl(nl, r));
        assert(a.hpwl >= lower - 1e-6);
        // determinism
        AnnealResult b = anneal(nl, r, prm);
        assert(b.hpwl == a.hpwl);
    }
    // Global placement + legalization is legal
    {
        Netlist nl = chain(40, 8, 40);
        GlobalResult g = global_place(nl, 20);
        Placement p = legalize(nl, g.x, g.y);
        assert(is_legal(nl, p));
    }
    std::puts("all C++ tests passed");
    return 0;
}
