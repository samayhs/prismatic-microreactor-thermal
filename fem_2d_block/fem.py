"""
From-scratch 2D finite element solver for nonlinear steady heat conduction.

Implements (see docs/02_architecture.md sec. 3):

  strong form   : div(k(T) grad T) + q''' = 0
  coolant walls : -k dT/dn = h_c (T - T_c)                 (Robin)
  outer surface : -k dT/dn = h_inf (T - T_inf)
                            + eps*sigma (T^4 - T_inf^4)     (convection + radiation)

Discretization: linear (P1) triangles. The problem is nonlinear through k(T) and
through the T^4 radiation term, and is solved by Newton-Raphson with a fully
ANALYTIC tangent (Jacobian) -- including the dk/dT conduction term and the
4*eps*sigma*T^3 radiation term. Correctness of that tangent is what gives
quadratic Newton convergence, and verify.py checks it directly.

The assembly is written out explicitly (element loops, closed-form P1 matrices,
2-point Gauss on boundary edges) rather than delegated to a PDE package -- that is
the point of the exercise (requirement SR-5).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from materials import (
    CONDUCTIVITY,
    MAT_FUEL,
    STEFAN_BOLTZMANN,
    TRISO_LIMIT_K,
    Scenario,
)

# 2-point Gauss-Legendre on the reference edge xi in [-1, 1]
_GP = np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
_GW = np.array([1.0, 1.0])


# ---------------------------------------------------------------------------
# Precomputed element geometry (built once; reused every Newton iteration)
# ---------------------------------------------------------------------------
@dataclass
class ElementGeometry:
    area: np.ndarray      # (Ne,) positive element areas [m^2]
    b: np.ndarray         # (Ne,3) shape-function x-gradient numerators
    c: np.ndarray         # (Ne,3) shape-function y-gradient numerators
    S: np.ndarray         # (Ne,3,3) geometric conduction matrix (k=1): (b_m b_n + c_m c_n)/(4A)


def precompute_geometry(nodes: np.ndarray, tris: np.ndarray) -> ElementGeometry:
    xi = nodes[tris[:, 0]]
    xj = nodes[tris[:, 1]]
    xk = nodes[tris[:, 2]]
    # signed 2A
    two_a = ((xj[:, 0] - xi[:, 0]) * (xk[:, 1] - xi[:, 1])
             - (xk[:, 0] - xi[:, 0]) * (xj[:, 1] - xi[:, 1]))
    area = 0.5 * np.abs(two_a)
    # shape-function gradient numerators (cyclic)
    b = np.column_stack([xj[:, 1] - xk[:, 1],
                         xk[:, 1] - xi[:, 1],
                         xi[:, 1] - xj[:, 1]])
    c = np.column_stack([xk[:, 0] - xj[:, 0],
                         xi[:, 0] - xk[:, 0],
                         xj[:, 0] - xi[:, 0]])
    # S[e,m,n] = (b_m b_n + c_m c_n) / (4 A)
    S = (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) \
        / (4.0 * area[:, None, None])
    return ElementGeometry(area=area, b=b, c=c, S=S)


# ---------------------------------------------------------------------------
# Residual + tangent assembly
# ---------------------------------------------------------------------------
def assemble_system(nodes, tris, tri_mat, edges_cool, edges_out,
                    geom: ElementGeometry, T: np.ndarray, scen: Scenario):
    """Return (R, J): residual vector R(T) and sparse analytic tangent J = dR/dT.

    Convention: R = internal_conduction + boundary_flux - volumetric_source = 0.
    """
    nn = nodes.shape[0]
    R = np.zeros(nn)
    rows, cols, vals = [], [], []

    def add_block(dofs, ke):
        for a in range(len(dofs)):
            for b_ in range(len(dofs)):
                rows.append(dofs[a]); cols.append(dofs[b_]); vals.append(ke[a, b_])

    # --- volume: conduction (nonlinear via k(T)) + generation ---------------
    Te = T[tris].mean(axis=1)                      # element-mean temperature
    for mat, kfun in CONDUCTIVITY.items():
        pass  # conductivity evaluated per element below

    for e in range(tris.shape[0]):
        dofs = tris[e]
        kfun = CONDUCTIVITY[int(tri_mat[e])]
        k_e, dk_e = kfun(Te[e])
        S_e = geom.S[e]
        Te_vec = T[dofs]

        # residual: internal conduction  Ke(T) @ T_e
        Ke = k_e * S_e
        R[dofs] += Ke @ Te_vec

        # tangent: Ke  +  (dk/dT * dTe/dT_n) * (S @ T_e)_m
        #   dTe/dT_n = 1/3 (element-mean of 3 nodes)
        SdotT = S_e @ Te_vec                       # (3,)
        Jt = Ke + (dk_e / 3.0) * np.outer(SdotT, np.ones(3))
        add_block(dofs, Jt)

        # residual: volumetric heat generation (fuel only)  -> -f
        if int(tri_mat[e]) == MAT_FUEL and scen.q_fuel != 0.0:
            fe = scen.q_fuel * geom.area[e] / 3.0
            R[dofs] -= fe                          # each node gets A/3 * q'''

    # --- boundary edges: Robin (coolant), convection + radiation (outer) ----
    def edge_length(edge):
        p0, p1 = nodes[edge[0]], nodes[edge[1]]
        return np.hypot(p1[0] - p0[0], p1[1] - p0[1])

    # consistent edge mass matrix  (L/6)[[2,1],[1,2]]
    def robin_edge(edge, h, Tref):
        if h == 0.0:
            return
        L = edge_length(edge)
        d = [edge[0], edge[1]]
        Me = (L / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
        Tedge = T[d]
        # residual: h*Me@T - h*Tref*(L/2)*[1,1]
        R[d] += h * (Me @ Tedge) - h * Tref * (L / 2.0)
        # tangent: h*Me  (linear)
        add_block(d, h * Me)

    def radiation_edge(edge, eps, Tinf):
        L = edge_length(edge)
        d = [edge[0], edge[1]]
        T0, T1 = T[d]
        for xi, w in zip(_GP, _GW):
            N0 = 0.5 * (1.0 - xi)
            N1 = 0.5 * (1.0 + xi)
            N = np.array([N0, N1])
            Tg = N0 * T0 + N1 * T1
            jac = L / 2.0                          # ds = (L/2) dxi
            coeff = eps * STEFAN_BOLTZMANN
            # residual: eps*sigma*(Tg^4 - Tinf^4) * N * w * jac
            R[d] += coeff * (Tg**4 - Tinf**4) * N * w * jac
            # tangent: eps*sigma*4*Tg^3 * outer(N,N) * w * jac
            add_block(d, coeff * 4.0 * Tg**3 * np.outer(N, N) * w * jac)

    for edge in edges_cool:
        robin_edge(edge, scen.h_coolant, scen.T_coolant)

    for edge in edges_out:
        robin_edge(edge, scen.h_outer, scen.T_inf)
        if scen.radiation:
            radiation_edge(edge, scen.emissivity, scen.T_inf)

    J = sp.coo_matrix((vals, (rows, cols)), shape=(nn, nn)).tocsr()
    return R, J


# ---------------------------------------------------------------------------
# Newton-Raphson driver
# ---------------------------------------------------------------------------
@dataclass
class SolveResult:
    T: np.ndarray
    converged: bool
    residual_history: list = field(default_factory=list)
    peak_fuel_T: float = 0.0
    peak_graphite_T: float = 0.0
    energy_closure: float = 0.0
    scenario_name: str = ""

    @property
    def triso_margin(self) -> float:
        """Margin of peak fuel temperature below the TRISO limit [K]."""
        return TRISO_LIMIT_K - self.peak_fuel_T


def solve_thermal(mesh, scen: Scenario, T_init: float = 700.0,
                  tol: float = 1e-8, max_iter: int = 30,
                  verbose: bool = True) -> SolveResult:
    nodes, tris, tri_mat = mesh.nodes, mesh.tris, mesh.tri_mat
    edges_cool, edges_out = mesh.edges_cool, mesh.edges_out
    geom = precompute_geometry(nodes, tris)

    T = np.full(nodes.shape[0], float(T_init))
    history = []
    r0 = None
    converged = False

    for it in range(max_iter):
        R, J = assemble_system(nodes, tris, tri_mat, edges_cool, edges_out,
                               geom, T, scen)
        rnorm = np.linalg.norm(R)
        if r0 is None:
            r0 = rnorm if rnorm > 0 else 1.0
        rel = rnorm / r0
        history.append(rnorm)
        if verbose:
            print(f"  Newton it {it:2d} : |R| = {rnorm:.3e}   rel = {rel:.3e}")
        if rel < tol or rnorm < 1e-9:
            converged = True
            break
        dT = spla.spsolve(J, -R)
        # light damping only if the step is very large (radiation robustness)
        step = np.max(np.abs(dT))
        relax = 1.0 if step < 300.0 else 300.0 / step
        T = T + relax * dT

    peak_fuel, peak_graph = _peak_temperatures(tris, tri_mat, T)
    closure = energy_balance(mesh, geom, T, scen)

    return SolveResult(
        T=T, converged=converged, residual_history=history,
        peak_fuel_T=peak_fuel, peak_graphite_T=peak_graph,
        energy_closure=closure, scenario_name=scen.name,
    )


def _peak_temperatures(tris, tri_mat, T):
    fuel_nodes = np.unique(tris[tri_mat == MAT_FUEL])
    graph_nodes = np.unique(tris[tri_mat != MAT_FUEL])
    peak_fuel = float(T[fuel_nodes].max()) if fuel_nodes.size else float("nan")
    peak_graph = float(T[graph_nodes].max()) if graph_nodes.size else float("nan")
    return peak_fuel, peak_graph


# ---------------------------------------------------------------------------
# Global energy balance  (requirement FR-12 / V&V test V5)
# ---------------------------------------------------------------------------
def energy_balance(mesh, geom: ElementGeometry, T: np.ndarray, scen: Scenario):
    """Return relative closure error |Q_gen - Q_out| / Q_gen (per unit depth)."""
    tris, tri_mat = mesh.tris, mesh.tri_mat
    nodes = mesh.nodes

    # generated power [W/m depth]
    fuel_mask = tri_mat == MAT_FUEL
    Q_gen = scen.q_fuel * geom.area[fuel_mask].sum()

    def edge_len(edge):
        p0, p1 = nodes[edge[0]], nodes[edge[1]]
        return np.hypot(p1[0] - p0[0], p1[1] - p0[1])

    Q_out = 0.0
    for edge in mesh.edges_cool:
        L = edge_len(edge)
        Tm = 0.5 * (T[edge[0]] + T[edge[1]])
        Q_out += scen.h_coolant * L * (Tm - scen.T_coolant)
    for edge in mesh.edges_out:
        L = edge_len(edge)
        Tm = 0.5 * (T[edge[0]] + T[edge[1]])
        Q_out += scen.h_outer * L * (Tm - scen.T_inf)
        if scen.radiation:
            for xi, w in zip(_GP, _GW):
                N0, N1 = 0.5 * (1 - xi), 0.5 * (1 + xi)
                Tg = N0 * T[edge[0]] + N1 * T[edge[1]]
                Q_out += scen.emissivity * STEFAN_BOLTZMANN * \
                    (Tg**4 - scen.T_inf**4) * w * (L / 2.0)

    if Q_gen == 0.0:
        return abs(Q_out)
    return abs(Q_gen - Q_out) / abs(Q_gen)


if __name__ == "__main__":
    from mesh import build_block_mesh
    from materials import NORMAL, PASSIVE

    m = build_block_mesh()
    for scen in (NORMAL, PASSIVE):
        print(f"\n=== {scen.name} ===")
        res = solve_thermal(m, scen)
        print(f"  converged      : {res.converged}")
        print(f"  peak fuel T    : {res.peak_fuel_T - 273.15:8.1f} C")
        print(f"  peak graphite T: {res.peak_graphite_T - 273.15:8.1f} C")
        print(f"  TRISO margin   : {res.triso_margin:8.1f} K")
        print(f"  energy closure : {res.energy_closure*100:8.3f} %")
