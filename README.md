# PlaceOpt

**Standard-cell placement treated as combinatorial optimization. A C++17 simulated-annealing placer and an analytical flow (global placement, then legalization, then detailed placement) run on Yosys-synthesized netlists built from real SkyWater SKY130 HD cells. On small instances, results are checked against exact optima from OR-Tools CP-SAT.**

Live dashboard: **https://nags-gk.github.io/placeopt/** (placement viewer with an SA animation, flow comparison, annealing schedule, and optimality gaps)

## Pipeline

```
RTL (bench/designs.v) --Yosys synth + ABC--> SKY130 HD gate netlist --> placement netlist
                                                (cells, nets, pin offsets from LEF, I/O terminals)
   |-- C++ simulated annealing (VPR-style adaptive schedule, range limiting, incremental HPWL)
   |-- C++ analytical: quadratic global placement (clique model, PCG solver)
   |                   + FastPlace-style cell shifting --> greedy legalization --> warm-started SA (detailed placement)
   `-- OR-Tools CP-SAT exact model (small instances) --> optimality gap
```

| Component | What it does |
|---|---|
| `scripts/build_cells.py` | Parses the official `sky130_fd_sc_hd` LEFs for cell footprints and signal-pin centers, then writes an area/function Liberty subset (18 cells, no timing tables) for Yosys/ABC |
| `scripts/synth.py` | Synthesizes 11 RTL designs (ISCAS c17, adders, ALUs, multipliers, CRC, counter) and exports the placement netlists |
| `cpp/src/placer.cpp` | C++17 core: incremental-HPWL simulated annealing with swap and displacement moves on a row/site occupancy grid; a quadratic global placer solved by Jacobi-preconditioned conjugate gradient with cell-shifting anchors; a nearest-free-slot legalizer; a legality checker |
| `cpp/src/bindings.cpp` | pybind11 module `placeopt._core` |
| `placeopt/exact.py` | CP-SAT model with the same grid, pin offsets and terminals: per-row `NoOverlap` intervals, and min/max per-net bounding boxes in integer nm |
| `scripts/run_benchmark.py` | Flow comparison over 5 seeds, plus the optimality study |

## Results

Median HPWL (µm) over 5 seeds, at 60% utilization on a 2-core VM:

| Design | Cells | Random | SA | SA time | Global+LG | Global+LG+DP | G+LG+DP time |
|---|---|---|---|---|---|---|---|
| c17 | 6 | 67 | 37 | 0.00 s | 44 | 37 | 0.00 s |
| adder8 | 24 | 682 | 369 | 0.02 s | 483 | 372 | 0.00 s |
| crc16 | 57 | 1,748 | 924 | 0.05 s | 986 | 921 | 0.02 s |
| counter_cmp | 69 | 2,103 | 874 | 0.06 s | 1,077 | 841 | 0.02 s |
| adder32 | 97 | 5,291 | 2,620 | 0.11 s | 3,115 | 2,607 | 0.03 s |
| mult6 | 141 | 5,280 | 1,525 | 0.18 s | 1,551 | 1,471 | 0.04 s |
| alu8 | 166 | 5,954 | 2,051 | 0.24 s | 2,357 | 2,026 | 0.06 s |
| mult8 | 291 | 14,384 | 3,375 | 0.50 s | 3,416 | 3,136 | 0.10 s |
| mult12 | 696 | 51,465 | 9,009 | 1.80 s | 10,838 | 8,363 | 0.37 s |
| alu32 | 1106 | 87,917 | 18,944 | 5.08 s | 37,633 | 18,999 | 1.79 s |
| mult16 | 1236 | 117,867 | 16,927 | 4.07 s | 28,156 | 16,777 | 1.16 s |

**Takeaway:** the analytical flow followed by a warm-started detailed anneal matches or beats pure simulated annealing on 9 of 11 designs (it loses on adder8 by 0.8% and alu32 by 0.3%). It runs about 3.5× faster on the largest design (mult16, 1,236 cells).

**Optimality vs. CP-SAT** (on 12 instances that CP-SAT solved to proven optimality, 6–10 cells):
- Best of 10 SA seeds hits the exact optimum on 4 of 12 instances. The median gap is 0.93% for best-of-10 and 4.63% for a single seed.
- With 10× more moves per temperature, the median best-of-10 gap drops to 0.00%.
- At 10–24 cells, CP-SAT often can't prove optimality within 60–300 s, so those rows are compared with its best solution. This is exactly the scaling wall that makes heuristics necessary.

## Debugging note

An early version gave *worse* results with more annealing effort on tiny instances. The trace showed why: the range window collapsed to a single site and the anneal froze just after an uphill move. The fix keeps a 3-site floor on the window, makes 5% of moves global, and returns the best placement seen. After the fix, quality improves monotonically with effort. The dashboard's annealing-schedule plots show the acceptance and window dynamics.

## Honest caveats
- The objective is HPWL only. There is no timing, density or congestion term, so the analytical flow packs cells toward the input-pin edge.
- Cell orientation and flipping are not modeled, and the Liberty subset has no timing, so ABC maps for area only.
- I/O pins are fixed, evenly spaced on the left edge (inputs) and right edge (outputs). Clock nets are excluded because clock-tree synthesis routes them.
- Designs top out at about 1,200 cells. Industrial placers handle millions; the asymptotic comparison is in the runtime plot.

## Build and run

```bash
sudo apt-get install yosys cmake g++ && pip install pybind11 ortools pytest
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
ctest --test-dir build && pytest -q
python scripts/synth.py                 # RTL -> SKY130 netlists
python scripts/run_benchmark.py         # flows, viz, exact (CP-SAT runs take ~25 min)
python scripts/sa_long.py && python scripts/make_report.py
```

To regenerate the cell data from the official LEFs, sparse-clone `google/skywater-pdk-libs-sky130_fd_sc_hd` into `.sc_src` and run `python scripts/build_cells.py`.

## License

Apache-2.0. Cell dimensions and pin locations come from SkyWater's SKY130 HD library (Apache-2.0).
