import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from placeopt import _core  # noqa: E402
from placeopt.design import load  # noqa: E402
from placeopt.exact import solve  # noqa: E402
from placeopt.flows import run  # noqa: E402


@pytest.mark.parametrize("method", ["random", "sa", "gp_lg", "gp_lg_dp"])
def test_flows_produce_legal_placements(method):
    d = load("alu8")
    r = run(d, method, seed=3)
    p = _core.Placement(r["placement"]["row"], r["placement"]["site"])
    ok, why = _core.is_legal(d.nl, p)
    assert ok, why
    assert r["hpwl"] == pytest.approx(_core.total_hpwl(d.nl, p))


def test_optimizers_beat_random():
    d = load("mult6")
    rnd = run(d, "random", seed=0)["hpwl"]
    assert run(d, "sa", seed=0)["hpwl"] < 0.5 * rnd
    assert run(d, "gp_lg_dp", seed=0)["hpwl"] < 0.5 * rnd


def test_sa_is_deterministic_per_seed():
    d = load("adder8")
    assert run(d, "sa", seed=5)["hpwl"] == run(d, "sa", seed=5)["hpwl"]


def test_exact_matches_or_beats_sa_on_c17():
    d = load("c17")
    ex = solve(d, time_limit=30)
    assert ex["status"] == "OPTIMAL"
    sa = min(run(d, "sa", seed=s)["hpwl"] for s in range(5))
    assert ex["hpwl"] <= sa + 1e-3
    p = _core.Placement(ex["placement"]["row"], ex["placement"]["site"])
    assert _core.is_legal(d.nl, p)[0]
    assert _core.total_hpwl(d.nl, p) == pytest.approx(ex["hpwl"], abs=0.01)  # nm rounding


def test_netlists_reference_only_library_cells():
    cells = json.loads(Path("data/cells.json").read_text())["cells"]
    for f in Path("data/netlists").glob("*.json"):
        nl = json.loads(f.read_text())
        assert all(c["type"] in cells for c in nl["cells"]), f.name
