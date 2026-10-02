"""Synthesize bench/designs.v with Yosys onto the SKY130 HD subset and export
placement netlists (data/netlists/<top>.json)."""
import json, subprocess, sys
from pathlib import Path

TOPS = ["c17", "adder8", "crc16", "counter_cmp", "adder32", "mult6", "alu8", "mult8", "mult12", "alu32", "mult16"]
LIB = "data/sky130_hd_subset.lib"


def synth(top: str) -> dict:
    out = Path(f"data/yosys/{top}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    script = (f"read_verilog bench/designs.v; synth -flatten -top {top}; dfflibmap -liberty {LIB}; "
              f"abc -liberty {LIB}; opt_clean -purge; hilomap -hicell sky130_fd_sc_hd__buf_1 X -locell sky130_fd_sc_hd__buf_1 X; "
              f"stat -liberty {LIB}; write_json {out}")
    log = subprocess.run(["yosys", "-q", "-p", script], capture_output=True, text=True)
    if log.returncode:
        sys.exit(log.stderr)
    return json.loads(out.read_text())


def to_netlist(top: str, yj: dict, cells: dict) -> dict:
    mod = yj["modules"][top]
    insts, nets = [], {}
    for name, c in mod["cells"].items():
        if c["type"] not in cells:
            sys.exit(f"{top}: unmapped cell type {c['type']}")
        idx = len(insts)
        insts.append({"name": name, "type": c["type"]})
        for pin, bits in c["connections"].items():
            for b in bits:
                if isinstance(b, int):
                    nets.setdefault(b, []).append([idx, pin])
    ports = []
    for pname, p in mod["ports"].items():
        for i, b in enumerate(p["bits"]):
            if isinstance(b, int):
                pid = len(ports)
                ports.append({"name": f"{pname}[{i}]" if len(p["bits"]) > 1 else pname, "dir": p["direction"]})
                nets.setdefault(b, []).append([-1 - pid, ""])   # negative index = I/O terminal
    net_list = [pins for pins in nets.values() if len(pins) >= 2]
    return {"top": top, "cells": insts, "ports": ports, "nets": net_list}


def main():
    cells = json.loads(Path("data/cells.json").read_text())["cells"]
    Path("data/netlists").mkdir(parents=True, exist_ok=True)
    for top in TOPS:
        nl = to_netlist(top, synth(top), cells)
        Path(f"data/netlists/{top}.json").write_text(json.dumps(nl))
        area = sum(cells[c["type"]]["w"] * cells[c["type"]]["h"] for c in nl["cells"])
        print(f"{top:>12}: {len(nl['cells']):5d} cells, {len(nl['nets']):5d} nets, {len(nl['ports']):3d} I/O, "
              f"{area:8.1f} um2 cell area")


if __name__ == "__main__":
    main()
