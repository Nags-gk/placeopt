"""Build the SKY130 HD cell subset used for synthesis and placement.

Reads the official sky130_fd_sc_hd LEF files for cell size and signal-pin
locations, and writes:
  data/cells.json             - width/height (um) and pin centers per cell
  data/sky130_hd_subset.lib   - area + boolean function liberty for Yosys/ABC
The liberty has no timing tables: mapping is area-driven only.
"""
import json, re, sys
from pathlib import Path

FUNCS = {  # output pin, boolean function (Liberty syntax)
    "inv_1": ("Y", "!A"), "buf_1": ("X", "A"),
    "nand2_1": ("Y", "!(A&B)"), "nor2_1": ("Y", "!(A|B)"),
    "and2_1": ("X", "(A&B)"), "or2_1": ("X", "(A|B)"),
    "xor2_1": ("X", "(A^B)"), "xnor2_1": ("Y", "!(A^B)"),
    "nand3_1": ("Y", "!(A&B&C)"), "nor3_1": ("Y", "!(A|B|C)"),
    "a21oi_1": ("Y", "!((A1&A2)|B1)"), "o21ai_0": ("Y", "!((A1|A2)&B1)"),
    "a21o_1": ("X", "((A1&A2)|B1)"), "o21a_1": ("X", "((A1|A2)&B1)"),
    "a22oi_1": ("Y", "!((A1&A2)|(B1&B2))"), "maj3_1": ("X", "((A&B)|(A&C)|(B&C))"),
    "mux2_1": ("X", "((A0&!S)|(A1&S))"),
}
FLOPS = {"dfxtp_1": ("Q", "D", "CLK")}


def parse_lef(path: Path) -> dict:
    t = path.read_text()
    w, h = map(float, re.search(r"SIZE\s+([\d.]+)\s+BY\s+([\d.]+)", t).groups())
    pins = {}
    for m in re.finditer(r"PIN (\w+)(.*?)END \1", t, re.S):
        name, body = m.groups()
        if "USE SIGNAL" not in body and "USE CLOCK" not in body:
            continue
        rects = [tuple(map(float, r)) for r in re.findall(r"RECT ([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+)", body.split("PORT")[1])]
        area = sum((x2 - x1) * (y2 - y1) for x1, y1, x2, y2 in rects)
        cx = sum((x1 + x2) / 2 * (x2 - x1) * (y2 - y1) for x1, y1, x2, y2 in rects) / area
        cy = sum((y1 + y2) / 2 * (x2 - x1) * (y2 - y1) for x1, y1, x2, y2 in rects) / area
        pins[name] = {"x": round(cx, 4), "y": round(cy, 4),
                      "dir": "output" if "DIRECTION OUTPUT" in body else "input"}
    return {"w": w, "h": h, "pins": pins}


def main(src: str) -> None:
    src = Path(src)
    cells, lib = {}, ["library (sky130_hd_subset) {", '  delay_model : table_lookup ;',
                      '  time_unit : "1ns" ; voltage_unit : "1V" ; current_unit : "1mA" ;',
                      '  capacitive_load_unit (1,pf) ; area_unit : "1um2" ;']
    for c in list(FUNCS) + list(FLOPS):
        base = c.rsplit("_", 1)[0]
        lef = parse_lef(src / f"cells/{base}/sky130_fd_sc_hd__{c}.lef")
        name = f"sky130_fd_sc_hd__{c}"
        cells[name] = lef
        area = round(lef["w"] * lef["h"], 4)
        lib.append(f"  cell ({name}) {{\n    area : {area} ;")
        if c in FUNCS:
            out, fn = FUNCS[c]
            ins = sorted(set(re.findall(r"[A-Z]\w*", fn)))
            assert set(ins) | {out} == set(lef["pins"]), (c, ins, out, list(lef["pins"]))
            for p in ins:
                lib.append(f"    pin ({p}) {{ direction : input ; capacitance : 0.002 ; }}")
            lib.append(f'    pin ({out}) {{ direction : output ; function : "{fn}" ; }}')
        else:
            q, d, clk = FLOPS[c]
            assert {q, d, clk} == set(lef["pins"]), (c, list(lef["pins"]))
            lib.append(f'    ff (IQ, IQN) {{ next_state : "{d}" ; clocked_on : "{clk}" ; }}')
            lib.append(f"    pin ({d}) {{ direction : input ; capacitance : 0.002 ; }}")
            lib.append(f"    pin ({clk}) {{ direction : input ; clock : true ; capacitance : 0.002 ; }}")
            lib.append(f'    pin ({q}) {{ direction : output ; function : "IQ" ; }}')
        lib.append("  }")
    lib.append("}")
    Path("data").mkdir(exist_ok=True)
    Path("data/cells.json").write_text(json.dumps({"site_w": 0.46, "row_h": 2.72, "cells": cells}, indent=1))
    Path("data/sky130_hd_subset.lib").write_text("\n".join(lib) + "\n")
    print(f"{len(cells)} cells written")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".sc_src")
