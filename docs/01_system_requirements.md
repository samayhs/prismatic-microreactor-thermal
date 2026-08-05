# System Requirements — Prismatic Microreactor Thermal Model

**Document ID:** SRD-001  **Rev:** A  **Status:** Baseline

This document defines what the thermal model must represent, compute, and prove.
Requirements are grouped as **System (SR)**, **Functional (FR)**, **Numerical
(NR)**, and **Verification & Validation (VR)**. Each ID is referenced by the
architecture (`02_architecture.md`) and the V&V plan (`03_verification_validation.md`)
so intent → implementation → evidence is traceable.

Requirement language: **shall** = mandatory, **should** = goal.

---

## 1. Purpose and scope

The model shall predict the steady-state temperature field in the cross-section of
a prismatic, helium-cooled, TRISO-fueled reactor block, for the purpose of
evaluating **peak fuel temperature and its margin to the TRISO integrity limit**
under (a) normal forced-cooling operation and (b) loss of forced cooling with
passive radiative heat rejection.

Out of scope for this revision: neutronics/burnup, 3D axial variation, transient
time integration, fuel-particle-resolved (sub-compact) modeling, and mechanical
stress. These are noted as extensions in the architecture.

## 2. System requirements (SR)

| ID | Requirement |
|----|-------------|
| **SR-1** | The model shall represent a hexagonal graphite block with interleaved TRISO fuel compacts and helium coolant channels (prismatic HTGR arrangement). |
| **SR-2** | The model shall evaluate at minimum two load cases: **normal operation** (full power, forced convection) and **loss of forced cooling** (decay heat, passive rejection). |
| **SR-3** | The model shall report peak fuel temperature and the margin to the TRISO limit of **1600 °C (1873 K)** for each load case. |
| **SR-4** | The model shall be self-contained and reproducible from source with the stated open-source toolchain (Python + gmsh; OpenFOAM for the CFD cross-check). |
| **SR-5** | The FEM solver shall be implemented from first principles (assembly written explicitly, not delegated to a black-box PDE package) to demonstrate method knowledge. |

## 3. Functional requirements (FR)

### 3.1 Governing physics

| ID | Requirement |
|----|-------------|
| **FR-1** | The model shall solve steady-state heat conduction: ∇·(k ∇T) + q''' = 0. |
| **FR-2** | Thermal conductivity shall be **temperature-dependent**, k = k(T), distinct for graphite and fuel-compact materials, making the problem nonlinear. |
| **FR-3** | Volumetric heat generation q''' shall be applied only in fuel-compact subdomains, and shall be scalable (full power vs decay-heat fraction). |
| **FR-4** | Coolant channel walls shall impose a **Robin (convective) boundary condition**: −k ∂T/∂n = h(T − T_coolant). |
| **FR-5** | The outer block surface shall support a combined boundary condition: convection h_∞(T − T_∞) **plus** radiation εσ(T⁴ − T_∞⁴), with radiation individually switchable per scenario. |
| **FR-6** | The radiation boundary term shall be treated **nonlinearly** (no linearization of T⁴ about a fixed reference), consistent with FR-2's nonlinear solve. |

### 3.2 Discretization and solution

| ID | Requirement |
|----|-------------|
| **FR-7** | The domain shall be discretized with conforming linear (P1) triangular elements, with fuel/graphite material interfaces node-conforming. |
| **FR-8** | The nonlinear system shall be solved by **Newton–Raphson** using an **analytic tangent (Jacobian)** that includes dk/dT and the 4εσT³ radiation term. |
| **FR-9** | The solver shall report an iteration history (residual norm per iteration) and shall declare convergence against an explicit tolerance (see NR-3). |

### 3.3 Outputs

| ID | Requirement |
|----|-------------|
| **FR-10** | The model shall output the nodal temperature field, peak fuel temperature, peak graphite temperature, and margin to the TRISO limit. |
| **FR-11** | The model shall produce publication-quality figures: temperature contour on the mesh and Newton convergence history, per scenario. |
| **FR-12** | The model shall perform a **global energy balance check** (heat generated = heat rejected through coolant + outer boundary) and report the closure error. |

## 4. Numerical requirements (NR)

| ID | Requirement | Target |
|----|-------------|--------|
| **NR-1** | Spatial discretization error shall be controlled; the L2 temperature error for P1 elements shall exhibit observed order of accuracy p ≈ 2 under mesh refinement. | p ∈ [1.8, 2.2] |
| **NR-2** | The energy-balance closure error (FR-12) shall be small relative to total generated power. | ≤ 1.0 % |
| **NR-3** | Newton iteration shall converge to a relative residual below tolerance within a bounded iteration count. | ‖R‖/‖R₀‖ ≤ 1e-8 in ≤ 15 iterations |
| **NR-4** | Newton convergence should be asymptotically quadratic near the solution when the analytic tangent is correct (a direct check on FR-8). | residual roughly squares per step |

## 5. Verification & validation requirements (VR)

| ID | Requirement |
|----|-------------|
| **VR-1** | **Code verification:** the solver shall reproduce an analytical conduction solution (1D slab with uniform generation, fixed-conductivity) to within discretization tolerance. |
| **VR-2** | **Code verification (MMS):** the solver shall pass a Method of Manufactured Solutions test confirming the expected order of accuracy (NR-1). |
| **VR-3** | **Radiation verification:** the outer-boundary radiation term shall be verified against a 1D conduction–radiation problem with a known solution. |
| **VR-4** | **Solution verification:** peak fuel temperature shall be reported with a mesh-convergence study (≥3 mesh levels) and an extrapolated (Richardson) estimate with a GCI-style uncertainty band. |
| **VR-5** | **Cross-code check:** the FEM result for a representative sub-problem shall be compared against the independent OpenFOAM (finite-volume) solution; agreement shall be reported and discrepancies explained. |
| **VR-6** | **Sanity/validation context:** predicted temperatures shall fall within physically expected ranges for gas-cooled reactor fuel (order 10²–10³ °C) and be discussed against published HTGR fuel-temperature expectations. |

## 6. Acceptance criteria (roll-up)

The study is considered acceptance-complete when:

1. All FR items are implemented and exercised by `main.py`.
2. VR-1..VR-3 pass at the tolerances in NR-1..NR-2.
3. VR-4 produces a peak-fuel-temperature estimate with a stated numerical
   uncertainty band.
4. VR-5 reports FEM-vs-CFD agreement with quantified difference.
5. Both scenarios (SR-2) report peak fuel temperature and TRISO margin (SR-3).

## 7. Assumptions and limitations

- 2D cross-section; axial coolant heat-up and end effects are not resolved.
- Material correlations are representative of nuclear-grade materials, not a
  qualified vendor dataset; absolute temperatures are therefore indicative and the
  emphasis is on **method correctness and margin behavior**, not licensing-grade
  numbers.
- Coolant is modeled through a wall heat-transfer coefficient in the FEM (the CFD
  case resolves the coolant field directly, which is the point of the cross-check).
- Radiation on the outer surface is modeled as gray-body exchange with a large
  isothermal enclosure (view factor = 1); enclosure/view-factor radiation is an
  identified extension.
