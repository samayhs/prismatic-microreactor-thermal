# Verification & Validation Plan and Metrics

**Document ID:** VV-001  **Rev:** A
**Traces:** SRD-001 (VR-*, NR-*), ARC-001.

This plan follows the standard nuclear/thermal analysis distinction
(cf. ASME V&V 10/20):

- **Verification** — *are we solving the equations right?* (math/code correctness,
  independent of reality). Includes **code verification** (exact/manufactured
  solutions, order of accuracy) and **solution verification** (discretization
  uncertainty on the actual problem).
- **Validation** — *are we solving the right equations?* (agreement with physical
  reality / trusted references).

Each test below states the **method**, the **metric**, and the **quantitative
acceptance criterion** so a reviewer can see pass/fail objectively.

---

## 1. Code verification

### V1 — Analytical slab conduction with uniform generation  *(VR-1)*
- **Problem:** 1D slab, fixed k, uniform q''', symmetric, convective ends. Exact
  parabolic profile T(x) = T_wall + q'''/(2k)(L² − x²)-type solution.
- **Metric:** max nodal error `max|T_h − T_exact|` and L2 error.
- **Acceptance:** L2 relative error ≤ 1e-3 on a moderately refined mesh; error → 0
  under refinement.

### V2 — Method of Manufactured Solutions (MMS), order of accuracy  *(VR-2, NR-1)*
- **Method:** choose a smooth manufactured T*(x,y) (e.g. trigonometric), substitute
  into the PDE to derive the required source term, impose T* on boundaries, solve,
  measure error vs mesh size h across ≥4 levels.
- **Metric:** observed order of accuracy p from a log(error) vs log(h) fit.
- **Acceptance:** **p ∈ [1.8, 2.2]** for P1 elements (theory: 2). This is the single
  strongest evidence the assembly is correct.

### V3 — Conduction–radiation boundary verification  *(VR-3)*
- **Problem:** 1D bar, one end fixed, the other radiating εσ(T⁴−T∞⁴), with a
  reference solution (analytical for the linearized limit; high-resolution 1D
  numerical reference for the full nonlinear term).
- **Metric:** end-temperature and boundary-flux error vs reference.
- **Acceptance:** relative error ≤ 1e-3 vs the fine 1D reference.

### V4 — Newton tangent correctness  *(NR-4)*
- **Method (a):** finite-difference check of the analytic Jacobian: `‖J_analytic −
  J_FD‖ / ‖J_FD‖` at a random state.
- **Method (b):** confirm asymptotic quadratic convergence: residual ratio
  r_{n+1}/r_n² approximately constant near the solution.
- **Acceptance:** (a) relative Jacobian error ≤ 1e-6; (b) residual drops
  super-linearly (each late iteration roughly squares the previous residual).

### V5 — Global energy balance  *(FR-12, NR-2)*
- **Method:** integrate generated power ∫q'''dΩ and compare to total rejected power
  through coolant + outer boundaries computed from the converged field.
- **Metric:** closure error `|Q_gen − Q_out| / Q_gen`.
- **Acceptance:** **≤ 1.0 %.**

## 2. Solution verification (discretization uncertainty)

### S1 — Mesh convergence of peak fuel temperature  *(VR-4)*
- **Method:** solve the actual NORMAL scenario on ≥3 systematically refined meshes;
  apply Richardson extrapolation to estimate the mesh-independent peak temperature
  and a **Grid Convergence Index (GCI)**-style uncertainty band.
- **Metric:** extrapolated T_peak, observed order p, GCI (%).
- **Acceptance:** monotonic convergence with p > 1; reported T_peak carries a
  numerical uncertainty band (target GCI ≤ a few %).

## 3. Validation / cross-code

### C1 — FEM vs OpenFOAM (element vs volume)  *(VR-5)*
- **Method:** define a common sub-problem (single coolant channel + surrounding
  solid with generation). Solve in the FEM (wall-h representation) and in OpenFOAM
  `chtMultiRegionFoam` (resolved coolant). Compare solid peak temperature and
  channel-wall heat flux.
- **Metric:** relative difference in peak solid T and in integrated wall flux.
- **Acceptance:** agreement within combined numerical + modeling uncertainty
  (target ≤ ~5–10 %); any gap attributed to the wall-h vs resolved-flow modeling
  difference and discussed, not hidden.

### C2 — Physical plausibility / literature context  *(VR-6)*
- **Method:** confirm predicted fuel/graphite temperatures lie in the expected
  gas-cooled-reactor range and that the NORMAL case sits safely below, while the
  PASSIVE case demonstrates the intended margin behavior for decay-heat rejection.
- **Metric:** qualitative + magnitude comparison to published HTGR fuel-temperature
  expectations (order 600–1250 °C normal; passive peak below 1600 °C TRISO limit).
- **Acceptance:** results are physically reasonable and the margin narrative holds.

## 4. Metrics summary table

| Test | Metric | Acceptance | Requirement |
|------|--------|-----------|-------------|
| V1 Slab | L2 rel. error | ≤ 1e-3 | VR-1 |
| V2 MMS | observed order p | 1.8 – 2.2 | VR-2, NR-1 |
| V3 Radiation | rel. error vs 1D ref | ≤ 1e-3 | VR-3 |
| V4a Jacobian | ‖J−J_FD‖/‖J_FD‖ | ≤ 1e-6 | NR-4 |
| V4b Newton | quadratic residual drop | qualitative pass | NR-4 |
| V5 Energy | closure error | ≤ 1.0 % | FR-12, NR-2 |
| S1 Mesh | GCI on T_peak | reported, small | VR-4 |
| C1 FEM–CFD | rel. diff peak T & flux | ≤ ~5–10 % | VR-5 |
| C2 Physical | magnitude vs literature | reasonable | VR-6 |

## 5. Reporting

`verify.py` shall print a pass/fail line per test with the measured metric next to
its criterion, so the V&V status is reproducible with a single command. `main.py`
shall print the final peak-temperature / TRISO-margin table with the S1 numerical
uncertainty band attached.

## 6. Why this matters for a safety-relevant thermal model

Peak fuel temperature is a **licensing-relevant** quantity: the entire TRISO safety
case rests on staying below ~1600 °C, including under loss-of-forced-cooling. A
prediction is only as trustworthy as its verification. This plan makes the
numerical error **bounded and stated** (V2, S1), the physics terms **individually
checked** (V1, V3, V5), and the result **independently corroborated** by a second
discretization method (C1) — which is exactly the evidence chain a reviewer of a
thermal safety analysis looks for.
