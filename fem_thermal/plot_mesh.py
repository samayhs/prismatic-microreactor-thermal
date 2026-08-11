"""
Render the prismatic-block mesh for inspection.

Colors triangles by material (graphite vs fuel compact) and highlights the two
boundary-condition edge sets (coolant walls, outer surface). Writes
figures/mesh.png.

Run:  python plot_mesh.py
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from materials import MAT_FUEL
from mesh import build_block_mesh

FIG_DIR = os.path.join(os.path.dirname(__file__), "figures")


def _edge_segments(nodes, edges):
    return [[nodes[a], nodes[b]] for a, b in edges]


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    m = build_block_mesh()
    x, y = m.nodes[:, 0], m.nodes[:, 1]

    fig, ax = plt.subplots(figsize=(7.5, 7))

    # filled triangles colored by material, with mesh lines
    facecolors = np.where(m.tri_mat == MAT_FUEL, 1.0, 0.0)
    tpc = ax.tripcolor(x, y, m.tris, facecolors=facecolors,
                       cmap="Oranges", vmin=0, vmax=1.6,
                       edgecolors="k", linewidth=0.15)

    # boundary edge sets
    ax.add_collection(LineCollection(_edge_segments(m.nodes, m.edges_cool),
                                     colors="tab:blue", linewidths=2.0))
    ax.add_collection(LineCollection(_edge_segments(m.nodes, m.edges_out),
                                     colors="tab:red", linewidths=2.0))

    # legend proxies
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor="#fdd0a2", edgecolor="k", label="graphite matrix"),
        Patch(facecolor="#f16913", edgecolor="k", label="fuel compact (q''')"),
        Line2D([0], [0], color="tab:blue", lw=2, label="coolant-wall BC (Robin)"),
        Line2D([0], [0], color="tab:red", lw=2, label="outer-surface BC (conv+rad)"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.9)

    ax.set_aspect("equal")
    ax.set_title(f"Prismatic block mesh — {m.n_nodes} nodes, "
                 f"{m.tris.shape[0]} P1 triangles")
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    fig.tight_layout()

    path = os.path.join(FIG_DIR, "mesh.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Mesh figure written: {path}")
    print(f"  {m.n_nodes} nodes, {m.tris.shape[0]} triangles")
    print(f"  coolant edges: {m.edges_cool.shape[0]}, "
          f"outer edges: {m.edges_out.shape[0]}")


if __name__ == "__main__":
    main()
