"""Placement flows built from the C++ primitives."""
from __future__ import annotations

import time

from . import _core
from .design import Design


def _params(seed: int, **kw) -> "_core.AnnealParams":
    p = _core.AnnealParams()
    p.seed = seed
    for k, v in kw.items():
        setattr(p, k, v)
    return p


def run(d: Design, method: str, seed: int = 0, snapshots: int = 0) -> dict:
    """Run one flow. Returns HPWL (um), runtime (s), placement and trace."""
    nl = d.nl
    t = time.perf_counter()
    trace, snaps, extra = [], [], {}
    if method == "random":
        p = _core.random_placement(nl, seed)
    elif method == "sa":
        r = _core.anneal(nl, _core.random_placement(nl, seed), _params(seed, snapshots=snapshots))
        p, trace, snaps = r.placement, r.trace, r.snapshots
    elif method == "gp_lg":
        g = _core.global_place(nl, 30)
        p = _core.legalize(nl, g.x, g.y)
        extra["gp"] = {"x": list(g.x), "y": list(g.y), "overflow": list(g.overflow)}
    elif method == "gp_lg_dp":
        # Analytical global placement, Tetris legalization, then a cool, local SA
        # (detailed placement) seeded from the legal result.
        g = _core.global_place(nl, 30)
        p0 = _core.legalize(nl, g.x, g.y)
        extra["hpwl_after_lg"] = _core.total_hpwl(nl, p0)
        # Warm start: local window and a temperature scaled to the mean net length,
        # so the anneal refines the analytical solution instead of melting it.
        t0 = 0.2 * extra["hpwl_after_lg"] / max(1, len(nl.nets))
        r = _core.anneal(nl, p0, _params(seed, start_window_full=False, init_window=0.1, init_temp=t0,
                                         snapshots=snapshots))
        p, trace, snaps = r.placement, r.trace, r.snapshots
    else:
        raise ValueError(method)
    secs = time.perf_counter() - t
    ok, why = _core.is_legal(nl, p)
    if not ok:
        raise RuntimeError(f"{method} produced an illegal placement: {why}")
    return {"method": method, "seed": seed, "hpwl": _core.total_hpwl(nl, p), "seconds": secs,
            "placement": {"row": list(p.row), "site": list(p.site)},
            "trace": [{"temp": s.temp, "cost": s.cost, "accept": s.accept_rate, "window": s.window, "moves": s.moves}
                      for s in trace],
            "snapshots": [{"row": list(s.row), "site": list(s.site)} for s in snaps], **extra}
