"""
Verification suite (docs/03_verification_validation.md).

Runs code-verification tests with quantitative pass/fail metrics:

  V0  Patch test        -- linear field reproduced to machine precision
  V2  MMS order         -- observed order of accuracy p ~ 2 for P1 elements   (NR-1)
  V4a Jacobian check    -- analytic tangent vs finite-difference tangent       (NR-4)
  V4b Newton quadratic  -- residual drops super-linearly near the solution     (NR-4)
  V5  Energy balance     -- generation = rejection closure                      (NR-2)

Each test prints  [PASS]/[FAIL]  the measured metric, and its acceptance criterion.

V0/V2 use a self-contained Poisson assembler on a structured unit-square mesh (an
INDEPENDENT implementation path from fem.py's nonlinear assembler), which is good
verification practice. V4/V5 exercise the actual production assembler in fem.py.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from fem import (assemble_system, precompute_geometry, solve_thermal)


# ===========================================================================
# Structured unit-square mesh + constant-k Poisson assembler (for V0, V2)
# ===========================================================================
def unit_square_mesh(n: int):
    """n x n grid of nodes on [0,1]^2, split into 2 triangles per cell."""
    xs = np.linspace(0.0, 1.0, n)
    X, Y = np.meshgrid(xs, xs)
    nodes = np.column_stack([X.ravel(), Y.ravel()])

    def idx(i, j):
        return j * n + i

    tris = []
    for j in range(n - 1):
        for i in range(n - 1):
            a, b, c, d = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            tris.append([a, b, c])
            tris.append([a, c, d])
    tris = np.array(tris, dtype=int)

    tol = 1e-12
    on_boundary = ((np.abs(nodes[:, 0]) < tol) | (np.abs(nodes[:, 0] - 1) < tol)
                   | (np.abs(nodes[:, 1]) < tol) | (np.abs(nodes[:, 1] - 1) < tol))
    return nodes, tris, np.where(on_boundary)[0]


def assemble_poisson(nodes, tris, source_fn):
    """Assemble K, F for -div(grad T) = f with unit conductivity (P1)."""
    geom = precompute_geometry(nodes, tris)
    nn = nodes.shape[0]
    rows, cols, vals = [], [], []
    F = np.zeros(nn)
    for e in range(tris.shape[0]):
        d = tris[e]
        Ke = geom.S[e]                       # constant k = 1
        for a in range(3):
            for b in range(3):
                rows.append(d[a]); cols.append(d[b]); vals.append(Ke[a, b])
        centroid = nodes[d].mean(axis=0)     # 1-point load quadrature
        fe = source_fn(centroid[0], centroid[1]) * geom.area[e] / 3.0
        F[d] += fe
    K = sp.coo_matrix((vals, (rows, cols)), shape=(nn, nn)).tocsr()
    return K, F


def apply_dirichlet(K, F, bnodes, gfun, nodes):
    K = K.tolil()
    for nd in bnodes:
        K.rows[nd] = [nd]
        K.data[nd] = [1.0]
        F[nd] = gfun(nodes[nd, 0], nodes[nd, 1])
    return K.tocsr(), F


def l2_error(nodes, tris, Th, Tstar_fn):
    geom = precompute_geometry(nodes, tris)
    err2 = 0.0
    for e in range(tris.shape[0]):
        d = tris[e]
        ex = Th[d] - np.array([Tstar_fn(*nodes[n]) for n in d])
        err2 += geom.area[e] / 3.0 * np.sum(ex**2)   # lumped element integral
    return np.sqrt(err2)


# ===========================================================================
# V0 -- Patch test (linear exactness)
# ===========================================================================
def test_patch():
    Tstar = lambda x, y: 3.0 + 2.0 * x - 1.5 * y      # linear -> P1-exact
    source = lambda x, y: 0.0                          # laplacian(linear)=0
    nodes, tris, bnodes = unit_square_mesh(9)
    K, F = assemble_poisson(nodes, tris, source)
    K, F = apply_dirichlet(K, F, bnodes, Tstar, nodes)
    Th = spla.spsolve(K, F)
    err = l2_error(nodes, tris, Th, Tstar)
    ok = err < 1e-10
    print(f"[{'PASS' if ok else 'FAIL'}] V0 patch test        : "
          f"L2 err = {err:.2e}   (accept < 1e-10)")
    return ok


# ===========================================================================
# V2 -- MMS order of accuracy
# ===========================================================================
def test_mms_order():
    k = np.pi
    Tstar = lambda x, y: np.sin(k * x) * np.sin(k * y)
    source = lambda x, y: 2.0 * k * k * np.sin(k * x) * np.sin(k * y)  # -lap T*

    hs, errs = [], []
    for n in (9, 17, 33, 65):
        nodes, tris, bnodes = unit_square_mesh(n)
        K, F = assemble_poisson(nodes, tris, source)
        K, F = apply_dirichlet(K, F, bnodes, Tstar, nodes)
        Th = spla.spsolve(K, F)
        errs.append(l2_error(nodes, tris, Th, Tstar))
        hs.append(1.0 / (n - 1))
    hs, errs = np.array(hs), np.array(errs)
    p = np.polyfit(np.log(hs), np.log(errs), 1)[0]
    ok = 1.8 <= p <= 2.2
    print(f"[{'PASS' if ok else 'FAIL'}] V2 MMS order          : "
          f"observed p = {p:.3f}   (accept 1.8-2.2)")
    for h, e in zip(hs, errs):
        print(f"         h = {h:.4f}   L2 err = {e:.3e}")
    return ok


# ===========================================================================
# V4a -- analytic Jacobian vs finite difference (production assembler)
# ===========================================================================
def test_jacobian(mesh, scen):
    geom = precompute_geometry(mesh.nodes, mesh.tris)
    rng = np.random.default_rng(0)
    T = 700.0 + 100.0 * rng.standard_normal(mesh.n_nodes)
    args = (mesh.nodes, mesh.tris, mesh.tri_mat, mesh.edges_cool,
            mesh.edges_out, geom)
    R0, J = assemble_system(*args, T, scen)
    J = J.tocsc()

    eps = 1e-4
    cols = rng.choice(mesh.n_nodes, size=min(25, mesh.n_nodes), replace=False)
    num, den = 0.0, 0.0
    for j in cols:
        Tp = T.copy(); Tp[j] += eps
        Rp, _ = assemble_system(*args, Tp, scen)
        fd = (Rp - R0) / eps
        num += np.sum((fd - J[:, j].toarray().ravel())**2)
        den += np.sum(fd**2)
    rel = np.sqrt(num / den)
    ok = rel < 1e-6
    print(f"[{'PASS' if ok else 'FAIL'}] V4a Jacobian (FD)     : "
          f"rel diff = {rel:.2e}   (accept < 1e-6)")
    return ok


# ===========================================================================
# V4b -- Newton quadratic convergence
# ===========================================================================
def test_newton_quadratic(mesh, scen):
    res = solve_thermal(mesh, scen, verbose=False)
    h = np.array(res.residual_history)
    # look at the last few non-tiny residuals; each should be << previous
    ratios = h[1:] / h[:-1]
    fast = np.sum(ratios < 0.1) >= 2          # at least two steps cut by >10x
    ok = res.converged and fast
    print(f"[{'PASS' if ok else 'FAIL'}] V4b Newton quadratic  : "
          f"converged in {len(h)-1} its, min ratio = {ratios.min():.1e}")
    return ok


# ===========================================================================
# V5 -- energy balance closure
# ===========================================================================
def test_energy(mesh, scen):
    res = solve_thermal(mesh, scen, verbose=False)
    ok = res.energy_closure <= 0.01
    print(f"[{'PASS' if ok else 'FAIL'}] V5 energy balance     : "
          f"closure = {res.energy_closure*100:.3f} %   (accept <= 1.0 %)")
    return ok


def main():
    from mesh import build_block_mesh
    from materials import NORMAL

    print("=" * 62)
    print("  VERIFICATION SUITE  (docs/03_verification_validation.md)")
    print("=" * 62)
    results = []
    results.append(test_patch())
    results.append(test_mms_order())

    print("-" * 62)
    print("  Production assembler checks (block mesh, NORMAL scenario)")
    print("-" * 62)
    mesh = build_block_mesh(lc=0.012)          # coarse mesh keeps FD check fast
    results.append(test_jacobian(mesh, NORMAL))
    results.append(test_newton_quadratic(mesh, NORMAL))
    results.append(test_energy(mesh, NORMAL))

    print("=" * 62)
    npass = sum(results)
    print(f"  RESULT: {npass}/{len(results)} tests passed")
    print("=" * 62)
    return 0 if npass == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
