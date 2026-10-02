"""PlaceOpt benchmarks.

1. Flow comparison on synthesized SKY130 designs: random, SA, GP+LG, GP+LG+DP.
2. Optimality study: SA vs. exact CP-SAT optimum on small instances.
Writes results/flows.json, results/exact.json, results/viz.json.
"""
import json, random, statistics, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from placeopt import _core  # noqa: E402
from placeopt.design import Design, load  # noqa: E402
from placeopt.exact import solve  # noqa: E402
from placeopt.flows import run  # noqa: E402

DESIGNS = ["c17", "adder8", "crc16", "counter_cmp", "adder32", "mult6", "alu8", "mult8", "mult12", "alu32", "mult16"]
SEEDS = 5


def flows():
    out = {}
    for name in DESIGNS:
        d = load(name)
        rec = {"cells": d.nl.num_cells(), "nets": len(d.nl.nets), "rows": d.nl.rows, "sites": d.nl.sites,
               "die_um": d.die, "methods": {}}
        for m in ["random", "sa", "gp_lg", "gp_lg_dp"]:
            seeds = [0] if m == "gp_lg" else range(SEEDS)
            rs = [run(d, m, seed=s) for s in seeds]
            rec["methods"][m] = {"hpwl": [r["hpwl"] for r in rs], "seconds": [r["seconds"] for r in rs],
                                 "median_hpwl": statistics.median(r["hpwl"] for r in rs),
                                 "median_seconds": statistics.median(r["seconds"] for r in rs)}
            if m == "gp_lg_dp":
                rec["methods"][m]["hpwl_after_lg"] = rs[0]["hpwl_after_lg"]
        out[name] = rec
        ms = rec["methods"]
        print(f'{name:>12} {rec["cells"]:5d} cells | random {ms["random"]["median_hpwl"]:9.0f} | '
              f'SA {ms["sa"]["median_hpwl"]:8.0f} ({ms["sa"]["median_seconds"]:.2f}s) | '
              f'GP+LG {ms["gp_lg"]["median_hpwl"]:8.0f} | GP+LG+DP {ms["gp_lg_dp"]["median_hpwl"]:8.0f} '
              f'({ms["gp_lg_dp"]["median_seconds"]:.2f}s)', flush=True)
    return out


def tiny(n_cells: int, seed: int) -> Design:
    """Random small netlist built from real SKY130 cell footprints."""
    rng = random.Random(seed)
    lib = json.loads((ROOT / "data/cells.json").read_text())
    comb = [k for k in lib["cells"] if "dfxtp" not in k]
    types = [rng.choice(comb) for _ in range(n_cells)]
    widths = [round(lib["cells"][t]["w"] / lib["site_w"]) for t in types]
    nl = _core.Netlist()
    nl.rows = 2 if n_cells <= 8 else 3
    nl.sites = int(sum(widths) / (0.7 * nl.rows)) + 1
    nl.width = widths
    W, H = nl.sites * nl.site_w, nl.rows * nl.row_h
    nl.term_x, nl.term_y = [0.0, 0.0, W, W], [H * 0.25, H * 0.75, H * 0.3, H * 0.7]
    nets = []
    outs = {c: [p for p, v in lib["cells"][types[c]]["pins"].items() if v["dir"] == "output"][0] for c in range(n_cells)}
    ins = {c: [p for p, v in lib["cells"][types[c]]["pins"].items() if v["dir"] == "input"] for c in range(n_cells)}
    for c in range(n_cells):        # each output drives 1-3 random inputs or an output terminal
        sinks = rng.sample([x for x in range(n_cells) if x != c], k=min(n_cells - 1, rng.randint(1, 3)))
        pins = [_core.Pin(c, *[lib["cells"][types[c]]["pins"][outs[c]][k] for k in ("x", "y")])]
        for s in sinks:
            p = rng.choice(ins[s])
            pins.append(_core.Pin(s, lib["cells"][types[s]]["pins"][p]["x"], lib["cells"][types[s]]["pins"][p]["y"]))
        if rng.random() < 0.25:
            pins.append(_core.Pin(-3 - rng.randint(0, 1)))
        nets.append(pins)
    for t in (0, 1):                 # input terminals drive random cells
        s = rng.randrange(n_cells)
        p = rng.choice(ins[s])
        nets.append([_core.Pin(-1 - t), _core.Pin(s, lib["cells"][types[s]]["pins"][p]["x"], lib["cells"][types[s]]["pins"][p]["y"])])
    nl.nets = nets
    nl.build()
    return Design(f"tiny{n_cells}_{seed}", [{"name": f"u{i}", "type": t} for i, t in enumerate(types)], [], nl, lib["cells"], 0.7)


def exact_study():
    cases = [("c17", load("c17"))] + [(f"rand{n}_s{s}", tiny(n, s)) for n in (6, 8, 10, 12) for s in range(5)] \
        + [("adder8", load("adder8"))]
    out = []
    for name, d in cases:
        sa = [run(d, "sa", seed=s)["hpwl"] for s in range(10)]
        best_sa = min(sa)
        ex = solve(d, time_limit=60 if d.nl.num_cells() <= 12 else 300,
                   hint=run(d, "sa", seed=int(min(range(10), key=lambda s: sa[s])))["placement"])
        row = {"case": name, "cells": d.nl.num_cells(), "sa_best": best_sa, "sa_median": statistics.median(sa),
               "sa_hit_rate": sum(abs(x - ex.get("hpwl", -1)) < 1e-3 for x in sa) / len(sa), **{k: ex.get(k) for k in ("status", "hpwl", "bound", "seconds")}}
        if ex.get("hpwl"):
            row["gap_best_pct"] = 100 * (best_sa - ex["hpwl"]) / ex["hpwl"]
            row["gap_median_pct"] = 100 * (row["sa_median"] - ex["hpwl"]) / ex["hpwl"]
            row["gap_to_bound_pct"] = 100 * (best_sa - ex["bound"]) / ex["bound"]
        out.append(row)
        print(f'{name:>12} n={row["cells"]:3d} {row["status"]:>8} opt={row.get("hpwl")} bound={row.get("bound")} '
              f'SA best={best_sa:.3f} median={row["sa_median"]:.3f} hit={row["sa_hit_rate"]:.1f} ({row["seconds"]:.1f}s)', flush=True)
    return out


def viz():
    """Placements and an SA animation for the dashboard."""
    out = {}
    for name in ["mult8", "alu32"]:
        d = load(name)
        rec = {"rows": d.nl.rows, "sites": d.nl.sites, "site_w": d.nl.site_w, "row_h": d.nl.row_h,
               "width": list(d.nl.width), "types": [c["type"].replace("sky130_fd_sc_hd__", "") for c in d.cells],
               "term_x": list(d.nl.term_x), "term_y": list(d.nl.term_y),
               "nets": [[q.owner for q in net] for net in d.nl.nets], "methods": {}}
        for m in ["random", "sa", "gp_lg", "gp_lg_dp"]:
            r = run(d, m, seed=1, snapshots=14 if m == "sa" else 0)
            rec["methods"][m] = {"hpwl": r["hpwl"], "placement": r["placement"], "trace": r["trace"]}
            if m == "sa":
                rec["sa_snapshots"] = r["snapshots"]
            if "gp" in r:
                rec["gp_xy"] = {"x": r["gp"]["x"], "y": r["gp"]["y"]}
        out[name] = rec
    return out


if __name__ == "__main__":
    (ROOT / "results").mkdir(exist_ok=True)
    what = sys.argv[1:] or ["flows", "exact", "viz"]
    if "flows" in what:
        (ROOT / "results/flows.json").write_text(json.dumps(flows()))
    if "viz" in what:
        (ROOT / "results/viz.json").write_text(json.dumps(viz()))
    if "exact" in what:
        (ROOT / "results/exact.json").write_text(json.dumps(exact_study()))
