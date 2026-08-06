"""
Driver: run the NORMAL and PASSIVE scenarios, emit figures, print the
peak-temperature / TRISO-margin summary table (requirements FR-10, FR-11, SR-2/3).

Run:  python main.py
Outputs: figures/temperature_fields.png, figures/convergence.png
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

from fem import solve_thermal
from materials import NORMAL, PASSIVE, TRISO_LIMIT_K, MAT_FUEL
from mesh import build_block_mesh

FIG_DIR = os.path.join(os.path.dirname(__file__), "figures")


def _triangulation(mesh):
    return mtri.Triangulation(mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.tris)


def _peak_fuel_xy(mesh, T):
    fuel_nodes = np.unique(mesh.tris[mesh.tri_mat == MAT_FUEL])
    j = fuel_nodes[np.argmax(T[fuel_nodes])]
    return mesh.nodes[j]


def plot_temperature_fields(mesh, results):
    tri = _triangulation(mesh)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    for ax, res in zip(axes, results):
        TC = res.T - 273.15
        cf = ax.tricontourf(tri, TC, levels=40, cmap="inferno")
        ax.tricontour(tri, TC, levels=10, colors="k", linewidths=0.25, alpha=0.4)
        px, py = _peak_fuel_xy(mesh, res.T)
        ax.plot(px, py, "co", ms=7, mec="k", mew=0.8)
        ax.annotate(f"peak fuel\n{res.peak_fuel_T-273.15:.0f} C",
                    (px, py), textcoords="offset points", xytext=(8, 8),
                    color="cyan", fontsize=9, weight="bold")
        ax.set_aspect("equal")
        ax.set_title(f"{res.scenario_name}\n"
                     f"peak fuel {res.peak_fuel_T-273.15:.0f} C   |   "
                     f"TRISO margin {res.triso_margin:.0f} K",
                     fontsize=10)
        ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
        cb = fig.colorbar(cf, ax=ax, shrink=0.85)
        cb.set_label("Temperature [C]")
    fig.suptitle("Prismatic gas-cooled block - FEM temperature field "
                 "(nonlinear conduction + radiation)", fontsize=12, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    path = os.path.join(FIG_DIR, "temperature_fields.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_convergence(results):
    fig, ax = plt.subplots(figsize=(7, 5))
    for res in results:
        h = np.array(res.residual_history)
        ax.semilogy(range(len(h)), h / h[0], "o-", label=res.scenario_name)
    ax.set_xlabel("Newton iteration")
    ax.set_ylabel("relative residual  |R| / |R0|")
    ax.set_title("Newton-Raphson convergence (analytic tangent -> quadratic)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(FIG_DIR, "convergence.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def print_summary(results):
    print("\n" + "=" * 72)
    print(f"  {'Scenario':<42}{'Peak fuel':>11}{'TRISO margin':>15}")
    print("  " + "-" * 68)
    for res in results:
        short = res.scenario_name.split(" (")[0]
        print(f"  {short:<42}{res.peak_fuel_T-273.15:>8.0f} C"
              f"{res.triso_margin:>12.0f} K")
    print("  " + "-" * 68)
    print(f"  TRISO integrity limit: {TRISO_LIMIT_K-273.15:.0f} C "
          f"({TRISO_LIMIT_K:.0f} K)")
    print("=" * 72)


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    mesh = build_block_mesh()
    print(f"Mesh: {mesh.n_nodes} nodes, {mesh.tris.shape[0]} triangles\n")

    results = []
    for scen in (NORMAL, PASSIVE):
        print(f"Solving: {scen.name}")
        res = solve_thermal(mesh, scen, verbose=True)
        results.append(res)
        print()

    p1 = plot_temperature_fields(mesh, results)
    p2 = plot_convergence(results)
    print_summary(results)
    print(f"\nFigures written:\n  {p1}\n  {p2}")


if __name__ == "__main__":
    main()
