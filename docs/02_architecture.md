# Architecture — Prismatic Microreactor Thermal Model

**Document ID:** ARC-001  **Rev:** A
**Traces:** implements SRD-001 (requirements in parentheses, e.g. *(FR-1)*).

---

## 1. Overview

Two independent solvers attack the same physics so they can cross-check each other
(*VR-5*):

```
                    ┌─────────────────────────────────────────┐
                    │   Prismatic gas-cooled block (geometry)   │
                    └───────────────┬───────────────┬───────────┘
                                    │               │
                 Part 1: FEM (this repo)     Part 2: CFD (OpenFOAM)
              finite ELEMENT, conduction     finite VOLUME, conjugate
              + Robin + nonlinear radiation  heat transfer + fvDOM
                                    │               │
                                    └──────┬────────┘
                                           ▼
                              Peak fuel T, TRISO margin,
                              cross-verified & discussed
```

The FEM path is the primary deliverable (*SR-5*: written from first principles).
The CFD path resolves the coolant field that the FEM represents through a wall
coefficient, and provides an independent-method comparison.

## 2. Part 1 — FEM solver module design

Small, single-responsibility modules; the numerical core has no meshing or
plotting dependencies so it stays testable.

```
materials.py ──► mesh.py ──► fem.py ──► main.py ──► figures/
                                 ▲
                             verify.py  (exercises fem.py against benchmarks)
```

| Module | Responsibility | Key requirements |
|--------|----------------|------------------|
| `materials.py` | k(T) and dk/dT for graphite & fuel; scenario definitions (q''', h, T_coolant, T_∞, radiation on/off); TRISO limit; Stefan–Boltzmann. | FR-2, FR-3, FR-5, SR-3 |
| `mesh.py` | Build hex block + interleaved fuel/coolant channels in gmsh; return numpy `nodes, tris, tri_mat, edges_cool, edges_out`. | SR-1, FR-7 |
| `fem.py` | Assemble residual R(T) and tangent J(T); apply Robin & radiation boundary integrals; Newton–Raphson driver; energy-balance check. | FR-1,4,5,6,8,9,12; NR-3,4 |
| `verify.py` | Analytical-slab, MMS, and radiation-fin benchmarks; mesh-convergence / order-of-accuracy harness. | VR-1,2,3,4; NR-1,2 |
| `main.py` | Run NORMAL and PASSIVE scenarios; produce contour + convergence figures; print peak T / margin table. | FR-10,11; SR-2,3 |

### 2.1 Data contract (mesh → FEM)

```
nodes      (Nn,2) float   coordinates [m]
tris       (Ne,3) int     P1 triangle connectivity, 0-based
tri_mat    (Ne,)  int     MAT_GRAPHITE=0 | MAT_FUEL=1
edges_cool (Nc,2) int     boundary edges on coolant walls   → Robin BC
edges_out  (No,2) int     boundary edges on outer surface   → convection + radiation
```

## 3. Numerical formulation (the core of *FR-1, FR-6, FR-8*)

### 3.1 Strong form
Steady conduction with generation, k temperature-dependent:

  ∇·(k(T) ∇T) + q''' = 0   in Ω

Boundary conditions:
- Coolant walls Γ_c:  −k ∂T/∂n = h_c (T − T_c)                         *(FR-4)*
- Outer surface Γ_o:  −k ∂T/∂n = h_∞ (T − T_∞) + εσ (T⁴ − T_∞⁴)       *(FR-5)*

### 3.2 Weak form
Multiply by test function v, integrate by parts:

  ∫_Ω k ∇v·∇T dΩ + ∫_Γc h_c v T dΓ + ∫_Γo h_∞ v T dΓ + ∫_Γo εσ v T⁴ dΓ
     = ∫_Ω v q''' dΩ + ∫_Γc h_c v T_c dΓ + ∫_Γo h_∞ v T_∞ dΓ + ∫_Γo εσ v T_∞⁴ dΓ

### 3.3 Discrete residual and tangent (Newton)
With T ≈ Σ T_j N_j on P1 triangles, define the nonlinear residual R(T) = 0 as
(internal + boundary − source). Newton solves J ΔT = −R with the **analytic**
tangent

  J = ∂R/∂T = K_cond(T) + K_cond′(T)·T   (conduction + dk/dT contribution)
             + H_c                        (Robin, linear)
             + H_∞                        (outer convection, linear)
             + 4εσ ∫_Γo N_i N_j T³ dΓ     (radiation tangent)              *(FR-8)*

The `dk/dT` term in the conduction tangent and the `4εσT³` radiation term are the
two pieces that make Newton converge **quadratically** — `verify.py` checks that
convergence rate as a direct test that the Jacobian is correct *(NR-4)*.

### 3.4 Element technology
- **P1 linear triangles**: constant gradient per element → conduction matrix has a
  closed form via the element area and shape-function gradients.
- **Conductivity per element** evaluated at the element-mean temperature (a common,
  well-behaved choice for nonlinear conduction on P1 meshes).
- **Boundary integrals** on 2-node edges via 2-point Gauss for the radiation T⁴
  term; exact for the linear Robin term.

## 4. Part 2 — OpenFOAM CFD case design

| Element | Choice | Requirement |
|---------|--------|-------------|
| Solver | `chtMultiRegionFoam` (steady) — conjugate heat transfer across solid/fluid | CFD demo, VR-5 |
| Regions | `solid` (graphite+fuel with q'''), `fluid` (helium channel) | SR-1 |
| Radiation | `fvDOM` (finite-volume discrete ordinates) on the solid outer patch | radiation bullet |
| Turbulence/flow | laminar or k-ω SST depending on channel Re; helium as ideal gas | CFD |
| Cross-check | compare solid peak T and wall heat flux against FEM sub-problem | VR-5 |

The **FEM (finite element) vs OpenFOAM (finite volume)** contrast is intentional: it
lets the discretization method itself be discussed, and gives two independent
numerical routes to the same margin number.

## 5. Extension points (documented, not implemented this rev)

- Transient term ρc_p ∂T/∂t → implicit (backward-Euler) time integration for
  startup / LOFC transients.
- Enclosure radiation with view factors (replace VF=1 assumption).
- Temperature-dependent q''' coupling to a simple point-kinetics decay-heat curve
  (ANS-5.1) for a time-accurate passive-cooldown study.
- 3D extrusion with axial coolant enthalpy rise.

## 6. Traceability summary

| Requirement | Realized by |
|---|---|
| SR-1 | `mesh.py` geometry |
| SR-2, SR-3 | `materials.py` scenarios + `main.py` reporting |
| SR-5 | `fem.py` explicit assembly |
| FR-1,7 | `fem.py` conduction assembly, P1 triangles |
| FR-2,3 | `materials.py` k(T), q''' |
| FR-4,5,6 | `fem.py` boundary integrals |
| FR-8,9 | `fem.py` Newton + tangent + history |
| FR-10,11,12 | `main.py` outputs, energy balance |
| NR-1..4 | `verify.py` |
| VR-1..6 | `verify.py` + `openfoam/` + discussion in results |
