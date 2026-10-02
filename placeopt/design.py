"""Load a synthesized netlist, size the die, place I/O terminals, and build
the C++ Netlist object."""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

from . import _core

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


@dataclass
class Design:
    name: str
    cells: list          # [{"name", "type"}]
    ports: list          # [{"name", "dir"}]
    nl: "_core.Netlist"
    cell_lib: dict
    utilization: float

    @property
    def die(self) -> tuple[float, float]:
        return self.nl.sites * self.nl.site_w, self.nl.rows * self.nl.row_h


def load(name: str, utilization: float = 0.6, aspect: float = 1.0) -> Design:
    lib = json.loads((DATA / "cells.json").read_text())
    cells_lib = lib["cells"]
    src = json.loads((DATA / "netlists" / f"{name}.json").read_text())
    site_w, row_h = lib["site_w"], lib["row_h"]
    widths = [round(cells_lib[c["type"]]["w"] / site_w) for c in src["cells"]]
    area = sum(widths) * site_w * row_h / utilization
    rows = max(1, round(math.sqrt(area / aspect) / row_h))
    sites = max(max(widths), math.ceil(area / (rows * row_h) / site_w))
    nl = _core.Netlist()
    nl.rows, nl.sites, nl.site_w, nl.row_h = rows, sites, site_w, row_h
    nl.width = widths
    W, H = sites * site_w, rows * row_h
    # I/O terminals: inputs spread along the left edge, outputs along the right edge.
    ins = [i for i, p in enumerate(src["ports"]) if p["dir"] == "input"]
    outs = [i for i, p in enumerate(src["ports"]) if p["dir"] != "input"]
    tx, ty = [0.0] * len(src["ports"]), [0.0] * len(src["ports"])
    for group, x in ((ins, 0.0), (outs, W)):
        for k, i in enumerate(group):
            tx[i], ty[i] = x, H * (k + 0.5) / len(group)
    nl.term_x, nl.term_y = tx, ty
    nets = []
    for net in src["nets"]:
        # Clock nets are routed by clock-tree synthesis, not placement-driven; skip them.
        if any(pin == "CLK" for owner, pin in net if owner >= 0):
            continue
        pins = []
        for owner, pin in net:
            if owner >= 0:
                p = cells_lib[src["cells"][owner]["type"]]["pins"][pin]
                pins.append(_core.Pin(owner, p["x"], p["y"]))
            else:
                pins.append(_core.Pin(owner, 0.0, 0.0))
        nets.append(pins)
    nl.nets = nets
    nl.build()
    return Design(name, src["cells"], src["ports"], nl, cells_lib, utilization)


def placement_to_dict(d: Design, p) -> dict:
    return {"row": list(p.row), "site": list(p.site)}
