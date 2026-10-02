"""Exact minimum-HPWL placement with OR-Tools CP-SAT (small instances only).

Same model as the C++ placer: cells on a rows x sites grid, pin offsets from
the LEF, fixed I/O terminals. Coordinates are in integer nanometres, so the
optimum is exact up to 1 nm rounding of pin offsets.
"""
from __future__ import annotations

import time

from ortools.sat.python import cp_model

from .design import Design


def solve(d: Design, time_limit: float = 120.0, workers: int = 2, hint=None) -> dict:
    nl = d.nl
    n, R, S = nl.num_cells(), nl.rows, nl.sites
    sw, rh = round(nl.site_w * 1000), round(nl.row_h * 1000)
    m = cp_model.CpModel()
    site = [m.NewIntVar(0, S - nl.width[c], f"s{c}") for c in range(n)]
    row = [m.NewIntVar(0, R - 1, f"r{c}") for c in range(n)]
    in_row = [[m.NewBoolVar(f"b{c}_{r}") for r in range(R)] for c in range(n)]
    for c in range(n):
        m.AddExactlyOne(in_row[c])
        m.Add(row[c] == sum(r * in_row[c][r] for r in range(R)))
    for r in range(R):
        ivs = [m.NewOptionalIntervalVar(site[c], nl.width[c], site[c] + nl.width[c], in_row[c][r], f"i{c}_{r}")
               for c in range(n)]
        m.AddNoOverlap(ivs)
    # Symmetry breaking is not added: cells of the same type still differ by connectivity.
    obj = []
    W, H = S * sw, R * rh
    for k, net in enumerate(nl.nets):
        xs, ys = [], []
        for q in net:
            if q.owner >= 0:
                xs.append(site[q.owner] * sw + round(q.dx * 1000))
                ys.append(row[q.owner] * rh + round(q.dy * 1000))
            else:
                xs.append(round(nl.term_x[-1 - q.owner] * 1000))
                ys.append(round(nl.term_y[-1 - q.owner] * 1000))
        for coords, L, tag in ((xs, W, "x"), (ys, H, "y")):
            hi, lo = m.NewIntVar(-10_000, L + 10_000, f"{tag}hi{k}"), m.NewIntVar(-10_000, L + 10_000, f"{tag}lo{k}")
            m.AddMaxEquality(hi, coords)
            m.AddMinEquality(lo, coords)
            obj.append(hi - lo)
    m.Minimize(sum(obj))
    if hint is not None:
        for c in range(n):
            m.AddHint(row[c], hint["row"][c])
            m.AddHint(site[c], hint["site"][c])
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_workers = workers
    t = time.perf_counter()
    status = solver.Solve(m)
    secs = time.perf_counter() - t
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return {"status": solver.StatusName(status), "seconds": secs}
    return {"status": solver.StatusName(status), "seconds": secs,
            "hpwl": solver.ObjectiveValue() / 1000.0, "bound": solver.BestObjectiveBound() / 1000.0,
            "placement": {"row": [solver.Value(v) for v in row], "site": [solver.Value(v) for v in site]}}
