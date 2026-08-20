"""
Generate system/solid/topoSetDict that carves the six fuel compacts out of the
solid region into a cellZone named "fuel". The fvOptions heat source is applied to
that cellZone. Geometry comes from fuel_zones.json (written by make_unit_cell.py).

Run:  python make_toposet.py
"""

from __future__ import annotations

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)


def main():
    with open(os.path.join(HERE, "fuel_zones.json")) as fh:
        fz = json.load(fh)
    L = fz["length"]
    rf = fz["r_fuel"]
    centers = fz["centers"]

    lines = [
        "/*--------------------------------*- C++ -*----------------------------------*/",
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    class       dictionary;",
        '    location    "system/solid";',
        "    object      topoSetDict;",
        "}",
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",
        "",
        "// Six fuel compacts -> cellSet 'fuel' -> cellZone 'fuel' (heat-source zone).",
        "actions",
        "(",
    ]

    for i, (cx, cy) in enumerate(centers):
        action = "new" if i == 0 else "add"
        lines += [
            "    {",
            "        name    fuel;",
            "        type    cellSet;",
            f"        action  {action};",
            "        source  cylinderToCell;",
            "        sourceInfo",
            "        {",
            f"            p1      ({cx:.6f} {cy:.6f} 0);",
            f"            p2      ({cx:.6f} {cy:.6f} {L:.6f});",
            f"            radius  {rf:.6f};",
            "        }",
            "    }",
        ]

    # convert the cellSet into a cellZone
    lines += [
        "    {",
        "        name    fuel;",
        "        type    cellZoneSet;",
        "        action  new;",
        "        source  setToCellZone;",
        "        sourceInfo",
        "        {",
        "            set fuel;",
        "        }",
        "    }",
        ");",
        "",
        "// ************************************************************************* //",
        "",
    ]

    out = os.path.join(CASE, "system", "solid", "topoSetDict")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        fh.write("\n".join(lines))
    print(f"wrote {out}")
    print(f"  {len(centers)} fuel cylinders, radius {rf} m, length {L} m")


if __name__ == "__main__":
    main()
