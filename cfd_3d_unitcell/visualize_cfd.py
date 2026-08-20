"""
Visualize the 3D CHT results without ParaView — reads the raw OpenFOAM field files
(cell centres + temperature) and produces two figures with matplotlib:

  1. Mid-plane cross-section temperature map (fuel hot spots, cold coolant channel)
  2. Axial temperature profiles: coolant bulk T and peak fuel T vs. distance along
     the channel (shows the coolant heating up and the fuel running hotter downstream)

This runs on Windows Python (matplotlib), reading results copied back from the WSL
run (default: ./run_hex, written by run.sh). It parses the same ASCII fields as the
energy-balance script, so it needs no OpenFOAM tools and no ParaView.

Run:  python visualize_cfd.py [caseDir] [time]
Outputs: figures/cfd_cross_section.png, figures/cfd_axial_profile.png
"""

from __future__ import annotations

import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(HERE, "figures")
LENGTH = 0.80          # channel length [m]
R_COOL = 0.008         # coolant channel radius [m] (blank hole in the solid)
TRISO_C = 1600.0       # TRISO limit [C]


def _list(path):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)",
                  txt, re.S)
    if m:
        return np.array([float(x) for x in m.group(1).split()])
    m = re.search(r"internalField\s+uniform\s+([-\d.eE+]+)", txt)
    return None if not m else float(m.group(1))


def _cellzone(path, name):
    txt = open(path).read()
    sub = txt[txt.find("\n" + name + "\n"):]
    m = re.search(r"cellLabels\s+List<label>\s*\d+\s*\((.*?)\)", sub, re.S)
    return np.array([int(x) for x in m.group(1).split()])


def _region(case, time, region):
    d = os.path.join(case, time, region)
    x = _list(os.path.join(d, "Cx"))
    y = _list(os.path.join(d, "Cy"))
    z = _list(os.path.join(d, "Cz"))
    T = _list(os.path.join(d, "T")) - 273.15
    return x, y, z, T


def cross_section(case, time):
    # Solid-only view near the outlet (hottest section). Showing the solid alone
    # avoids the cold coolant channel swamping the colour scale, so the radial
    # gradient (cool graphite next to the channel -> hotter fuel compacts outward)
    # is visible. The channel appears as a blank hole (no solid cells there).
    sx, sy, sz, sT = _region(case, time, "solid")
    zc = 0.72 * LENGTH
    band = np.abs(sz - zc) < (LENGTH / 44.0)
    x, y, T = sx[band], sy[band], sT[band]
    tri = mtri.Triangulation(x, y)
    # mask triangles that bridge the coolant channel (centroid inside the bore) or
    # have a very long edge (also spanning the hole)
    xt, yt = x[tri.triangles], y[tri.triangles]
    cr = np.hypot(xt.mean(axis=1), yt.mean(axis=1))
    emax = np.max([np.hypot(xt[:, i] - xt[:, j], yt[:, i] - yt[:, j])
                   for i, j in ((0, 1), (1, 2), (2, 0))], axis=0)
    tri.set_mask((cr < R_COOL * 1.15) | (emax > 3.5 * R_COOL))
    fig, ax = plt.subplots(figsize=(7, 6.2))
    cf = ax.tricontourf(tri, T, levels=40, cmap="inferno")
    ax.tricontour(tri, T, levels=12, colors="k", linewidths=0.2, alpha=0.35)
    ax.set_aspect("equal")
    ax.set_title(f"CFD solid temperature — cross-section near outlet "
                 f"(z = {zc:.2f} m)\n"
                 f"local peak {T.max():.0f} C   (coolant channel = blank centre)",
                 fontsize=10)
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    cb = fig.colorbar(cf, ax=ax, shrink=0.85); cb.set_label("Temperature [C]")
    fig.tight_layout()
    p = os.path.join(FIG_DIR, "cfd_cross_section.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    return p, T.max()


def axial_profile(case, time):
    fx, fy, fz, fT = _region(case, time, "fluid")
    sx, sy, sz, sT = _region(case, time, "solid")
    fuel = _cellzone(os.path.join(case, "constant", "solid", "polyMesh", "cellZones"),
                     "fuel")
    nb = 40
    edges = np.linspace(0, LENGTH, nb + 1)
    zc = 0.5 * (edges[:-1] + edges[1:])
    cool = np.full(nb, np.nan)
    fuelmax = np.full(nb, np.nan)
    fuel_z, fuel_T = sz[fuel], sT[fuel]
    for i in range(nb):
        fm = (fz >= edges[i]) & (fz < edges[i + 1])
        if fm.any():
            cool[i] = fT[fm].mean()
        km = (fuel_z >= edges[i]) & (fuel_z < edges[i + 1])
        if km.any():
            fuelmax[i] = fuel_T[km].max()

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(zc, fuelmax, "o-", color="tab:red", label="peak fuel T")
    ax.plot(zc, cool, "s-", color="tab:blue", label="coolant bulk T (mean)")
    ax.axhline(TRISO_C, ls="--", color="k", lw=1, label="TRISO limit 1600 C")
    ax.set_xlabel("axial position z [m]  (flow direction ->)")
    ax.set_ylabel("Temperature [C]")
    ax.set_title("Axial temperature profile — coolant heats up, fuel peaks downstream")
    ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="center left")
    fig.tight_layout()
    p = os.path.join(FIG_DIR, "cfd_axial_profile.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    return p, np.nanmax(fuelmax)


def main():
    # default: read the in-place run (this case folder). run.sh writes results here.
    case = sys.argv[1] if len(sys.argv) > 1 else HERE
    time = sys.argv[2] if len(sys.argv) > 2 else "4000"
    os.makedirs(FIG_DIR, exist_ok=True)
    p1, peak = cross_section(case, time)
    p2, fpeak = axial_profile(case, time)
    print(f"peak fuel temperature: {peak:.1f} C")
    print(f"figures written:\n  {p1}\n  {p2}")


if __name__ == "__main__":
    main()
