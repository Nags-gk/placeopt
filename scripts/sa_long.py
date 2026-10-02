"""Recompute the SA columns of results/exact.json against the stored CP-SAT results:
default-effort SA (10 seeds) and 10x-effort SA (10 seeds, 10x moves per temperature)."""
import json, statistics, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from placeopt import _core  # noqa: E402
from placeopt.design import load  # noqa: E402
from run_benchmark import tiny  # noqa: E402

rows = json.loads((ROOT / "results/exact.json").read_text())
for r in rows:
    name = r["case"]
    if name.startswith("rand"):
        n, s = name[4:].split("_s")
        d = tiny(int(n), int(s))
    else:
        d = load(name)
    def sa(mpc):
        out = []
        for seed in range(10):
            p = _core.AnnealParams(); p.seed = seed; p.moves_per_cell = mpc
            out.append(_core.anneal(d.nl, _core.random_placement(d.nl, seed), p).hpwl)
        return out
    std, long = sa(20.0), sa(200.0)
    r["sa_best"], r["sa_median"], r["sa_long_best"] = min(std), statistics.median(std), min(long)
    if r.get("hpwl"):
        opt = r["hpwl"]
        r["sa_hit_rate"] = sum(x <= opt + 1e-3 for x in std) / len(std)
        r["gap_best_pct"] = 100 * (r["sa_best"] - opt) / opt
        r["gap_median_pct"] = 100 * (r["sa_median"] - opt) / opt
        r["gap_to_bound_pct"] = 100 * (r["sa_best"] - r["bound"]) / r["bound"]
        r["sa_beats_cpsat"] = r["sa_best"] < opt - 1e-3
        print(f'{name:>12} {r["status"]:>8}: SA best {r["gap_best_pct"]:+.2f}%  median {r["gap_median_pct"]:+.2f}%  '
              f'long {100 * (r["sa_long_best"] - opt) / opt:+.2f}%')
(ROOT / "results/exact.json").write_text(json.dumps(rows))
