# Outcomes Log

**Document ID:** OL-001  **Rev:** A
Chronological record of work done and results obtained — an engineering-notebook
view of the project, complementing the decision log (DL-001). Newest entries at the
bottom.

---

## 2026-08-04 — Project framing and baseline
- Established the goal: an interview portfolio piece demonstrating FEM, CFD, and
  radiation heat transport on a gas-cooled microreactor geometry (Radiant Kaleidos
  class).
- Authored the requirements (SRD-001), architecture (ARC-001), and V&V plan
  (VV-001) *before* implementation, to mirror safety-relevant analysis practice.
- Implemented `materials.py` (temperature-dependent k(T) with analytic dk/dT for
  graphite and fuel compact; TRISO limit; two scenarios) and `mesh.py` (gmsh
  prismatic hex block with interleaved fuel/coolant channels).
- **Mesh result:** 3,111 nodes / 5,970 triangles; fuel, graphite, coolant-wall, and
  outer-boundary entities correctly tagged.
- Initialized git; first commit pushed to
  `github.com/samayhs/prismatic-microreactor-thermal`.

## 2026-08-05 — FEM solver, verification, and figures
- Implemented `fem.py`: explicit P1 assembly, Robin (coolant) and nonlinear
  radiation (outer) boundary integrals, Newton–Raphson with the analytic tangent
  (dk/dT + 4εσT³), and a global energy-balance check.
- **First run flagged two physics issues** → corrective decisions ADR-11 and ADR-12:
  - **OUT-1:** PASSIVE initially ran *cooler* than NORMAL because the coolant
    channels were still acting as a heat sink. Fixed by setting h_coolant = 0
    (true loss of forced flow); decay heat now leaves by outer-surface radiation.
  - **OUT-2:** NORMAL peak was ~529 °C (low for "at power"). Raised nominal power
    density 30 → 70 MW/m³ to land in the realistic HTGR range.
- Implemented `verify.py` (code-verification suite) and `main.py` (scenario runs +
  figures).

### Verification results (5/5 PASS)
| Test | Metric | Result | Criterion |
|---|---|---|---|
| V0 Patch test | L2 error | 4.9e-16 | < 1e-10 |
| V2 MMS order of accuracy | observed p | **1.999** | 1.8–2.2 |
| V4a Analytic Jacobian vs FD | rel. diff | 2.05e-8 | < 1e-6 |
| V4b Newton quadratic convergence | steps / min ratio | 4 its, 2.8e-8 | pass |
| V5 Global energy balance | closure | 0.000 % | ≤ 1.0 % |

The MMS observed order **p = 1.999** matches the P1 theoretical value of 2 — the
primary evidence the assembly is correct.

### Scenario results
| Scenario | Peak fuel T | Peak graphite T | TRISO margin | Newton its | Energy closure |
|---|---|---|---|---|---|
| Normal operation (forced helium cooling) | 822 °C | 684 °C | 778 K | 4 | 0.000 % |
| Loss of forced cooling (decay heat, passive radiation) | 360 °C | 357 °C | 1240 K | 4 | 0.000 % |

- **Interpretation:** Under loss of forced cooling, decay heat is rejected entirely
  by radiation from the block's outer surface, and the fuel stabilizes ~1240 K below
  the 1600 °C TRISO limit — a quantified "walk-away safe" demonstration, which is the
  core safety argument for a microreactor of this class.
- Figures generated: `fem_2d_block/figures/temperature_fields.png` (per-scenario
  temperature contours with peak-fuel markers) and `.../convergence.png` (Newton
  residual histories showing quadratic drop).

## 2026-08-06 — Documentation consolidation
- Added the decision log (DL-001) capturing ADR-1..13.
- Added this outcomes log (OL-001).
- Updated the README to index the new documents and reflect Part 1 completion.

## 2026-08-11 — Interview take-home: 3D OpenFOAM CFD model + validation
Radiant advanced the process to a technical interview (Ka-Yen, Staff Thermal Modeling
Engineer, 2026-08-20). Take-home: build a **3D CFD** model to predict peak fuel
temperature and **demonstrate validation**. Agreed scope: representative unit cell,
normal-operation scenario.

### CHT model built (`cfd_3d_unitcell/`)
- 3D unit cell: hexagonal graphite prism, central helium coolant channel, six fuel
  compacts, 0.8 m tall. Parametric gmsh generator → `gmshToFoam` → `splitMeshRegions`.
- Two coupled regions: `fluid` (helium) and `solid` (graphite + fuel). Fuel modeled as
  a heat-source cellZone (`topoSet`) inside the solid to avoid disconnected regions.
- Solver: steady `chtMultiRegionFoam`; helium as compressible ideal gas; k-epsilon
  turbulence with wall functions; coupled thermal interface; adiabatic outer walls.
- Operating point: helium 3 MPa, inlet 300 °C, ~10 m/s (Re ~ 1.2e4), q''' = 7 MW/m³.

### Results and verification
- **OUT-3:** First converged run (tri-prism mesh, 74.9k cells): peak fuel **581 °C**.
  Internal checks passed — helium density 2.52 kg/m³ (matches hand-calc), fuel→wall
  conduction ΔT ≈ 11 °C (matches cylindrical-conduction hand-calc).
- **Tri-prism GCI study (VR-4):** 3 meshes (48k/75k/123k). Peak fuel 573→581→587 °C;
  observed order **p = 1.96** (asymptotic range); Richardson-extrapolated **604 °C**;
  fine-grid GCI 3.5% (±21 °C).
- **OUT-4 — mesh non-orthogonality:** `checkMesh` showed max ~86° (borderline) on the
  extruded triangular-prism mesh; root cause is triangles on the curved boundaries.
- **Mesh-quality fix (ADR-17):** a structured boundary layer was attempted but blocked
  by a gmsh limitation (its 2D-only BL field cannot be applied to the conformal
  internal coolant wall, which is adjacent to 3 surfaces after extrusion). Pivoted to
  **hexahedral recombination + graded near-wall refinement**: max non-orthogonality
  86° → 15°/37° (fluid/solid), average 58–62° → 3.5°/10°, cells now 100% hexahedra.
  See `cfd_3d_unitcell/MESH_QUALITY.md`.
- **Hex re-run:** peak fuel **668 °C** (28.3k hex), converged; y⁺ ≈ 18.6 on the channel
  wall. The +87 °C shift vs the tri mesh confirms near-wall resolution matters.
- **Hex GCI study (VR-4):** 3 hex meshes (14.6k / 28.3k / 56.5k). Peak fuel
  669.3 → 668.3 → 672.0 °C — grid-independent to **±2 °C (0.56% spread across a 4×
  cell range)**. Convergence is oscillatory (solution at its numerical noise floor),
  so the oscillation amplitude is reported as the uncertainty rather than a formal
  Richardson extrapolation. **Result: peak fuel = 670 ± 2 °C, ~996 °C TRISO margin.**
  The ~66 °C gap vs the tri-prism extrapolation (604 °C) is due to the hex mesh's
  graded near-wall resolution capturing the film ΔT the coarse tri mesh smeared — the
  hex value is the more trustworthy one.
- **Energy balance (validation):** three independent heat measures on the converged
  hex solution agree to **0.01%** — Q_gen = q'''·V_fuel = 3953.8 W, interface
  wall-flux Q_wall = 3953.8 W, coolant enthalpy rise Q_cool = ṁ·cp·ΔT = 3953.3 W.
  Mass conserves to 0.000% (ṁ = 5.05 g/s in = out); coolant bulk ΔT = 150.8 K
  (outlet 451 °C), matching the hand-calc (~154 K). Script:
  `cfd_3d_unitcell/validation/`.
- **OUT-6 — post-processing debug:** custom `surfaceFieldValue`/`volFieldValue`
  function objects failed (silent / objectRegistry errors) in the multi-region case.
  Root cause: the standalone `postProcess` utility loads only one region, so
  region-tagged FOs can't resolve. Fix: define FOs in `controlDict`'s `functions{}`
  with a `region` tag and run via the solver's `chtMultiRegionFoam -postProcess`
  (loads all regions). The fixed FOs return T_out = 723.83 K and ṁ = 5.047 g/s —
  identical to the raw-field Python parser, cross-validating both. Full write-up in
  `validation/TROUBLESHOOTING.md`; working block in `validation/energyBalance.functions`.

## 2026-08-16 — Transient readiness: temperature-dependent solid properties
- **Implemented Cp(T) and k(T) for the solid** (ADR-18): OpenFOAM `hPolynomial`
  (Cp = −281.9 + 3.8125·T − 1.6875e-3·T²) + `polynomial` transport (k = 88 − 0.03·T),
  replacing the constant `hConst`/`constIso`. Required before a transient — Cp sets the
  thermal inertia ρCp·∂T/∂t (it does not affect steady results).
- **Verified:** OpenFOAM parses the polynomials and the case converges in place.
- **IMPORTANT correction (variable-Cp artifact):** the Cp(T)/k(T) run reported 684.6 °C,
  but this was **traced to a numerical artifact**, not physics. Isolation on the same
  mesh: `constIso`+`hConst` (k=60, Cp=1800) → **667.0 °C**; polynomial framework with
  **constant** k=60 & Cp=1800 → **667.0 °C**; polynomial with **variable Cp(T)** → 684.5;
  add k(T) → 684.6. So the +18 °C comes entirely from making **Cp temperature-variable**,
  which OpenFOAM's enthalpy-based solid solver (`∇·(k/Cp ∇h)`) does not keep independent
  of Cp at steady state (it should be). **The physical steady peak is ≈ 667 °C.** Cp(T)
  is for transient use only; a transient must be verified (mesh-refinement test on the
  18 °C; relaxation to correct steady; lumped-capacitance check).
- **Fuel property values not yet distinct.** A lower k_compact would need the fuel as
  its own region; the 6 disconnected compacts make `splitMeshRegions` produce 8 regions
  (tested/confirmed). ρCp of a matrix compact ≈ graphite and the un-captured effect is
  ~2 °C, so the fuel shares the solid curves for now (distinct k_compact = future work).
- Environment: repo folder renamed `prismatic-microreactor-thermal`; OpenFOAM moved
  `/root` → `/opt`; CFD now runs in-place in the repo (output gitignored); `run.sh` made
  portable (requires OpenFOAM pre-sourced). 2D FEM dropped from scope.

## 2026-08-19 — Variable-Cp artifact traced and fixed with a temperature-based solver
- **Reopened the ADR-18 artifact** (the +18 °C from `Cp(T)`: 667 → 684.6 °C). Prior
  docs assumed it was a benign truncation that would shrink under mesh refinement.
- **Mesh-refinement study (VR-6):** ran constant-Cp vs variable-Cp steady on 4 hex
  meshes at uniform in-plane refinement, isolating the artifact as
  `δ = T_peak(varCp) − T_peak(constCp)` (identical k(T) in both, so the artifact-free
  limit is exactly 0). Reproduced the validated anchors (constCp = 666.9 ≈ 667.0;
  varCp = 684.6). **Finding: δ does NOT vanish** — it holds at 16.7–18.3 °C across a
  ~5× cell range (first sweep), and only creeps 18.0 → 16.1 °C when the near-wall solid
  is also refined (observed order ≈ 0.25, far below the ideal O(h²)). A methodology
  note: the first sweep pinned `BL_FIRST` (near-wall size), so it did not refine the
  steep-gradient region where the error lives — corrected in the second sweep. Harness:
  `cfd_3d_unitcell/verification/` (`run_artifact_study.sh`, `analyze_artifact.py`,
  `RESULTS.md`, `results.csv`, `results_uniform.csv`).
- **Root cause (derived by hand):** the stock solid solver is enthalpy-based —
  `heSolidThermo::calculate()` builds `alpha = kappa/Cpv(T)` (nodal Cp) and the diffused
  field is `h = ∫Cp dT`, so a face flux is `interp(kappa/Cp)·C̄p·(T_N−T_P)` with `C̄p`
  the interval-mean Cp. Nodal Cp (in α) and interval-mean Cp (in Δh) cancel only for
  constant Cp; with `Cp(T)` they leave a mesh-persistent bias. Steady conduction is
  Cp-independent, so the bias is non-physical.
- **Fix (ADR-19):** built `chtMultiRegionTFoam` — a copy of `chtMultiRegionFoam` whose
  solid solves `ρCp(T) ∂T/∂t = ∇·(k∇T) + q‴` directly in **temperature**. Heat source
  applied explicitly as a `W/m³` field on the `fuel` cellZone
  (`constant/solid/heatSource`). Built via `solver-chtMultiRegionTFoam/Allwmake`
  (forces g++-11; OpenFOAM 7 does not build under the system default g++-15). Nothing in
  `/opt/OpenFOAM-7` modified; binary lives in the user's `FOAM_USER_APPBIN`.
- **VR-7 — fix verified:** production mesh, **variable Cp(T)**, T-solver →
  **peak fuel = 666.92 °C**, matching the constant-Cp enthalpy result (666.91 °C) to
  **0.01 °C**. The +17.7 °C artifact is gone by construction. `Cp(T)` is retained and
  now acts only in the transient inertia term. Adopted as production: `controlDict`
  `application = chtMultiRegionTFoam`; `run.sh` auto-builds the solver on first run.
- **Steady peak fuel temperature stands at ≈ 667 °C** — now obtained *with* the physical
  `Cp(T)`, not by falling back to constant Cp.

## Open items / next
- **Transient (VR-8, next):** with the T-solver the transient is now trustworthy in
  principle. Run a solid transient with `Cp(T)` inertia; verify (a) it relaxes to the
  667 °C steady, (b) its time constant matches the lumped/diffusion estimate
  τ ≈ 30–45 s (α = k/ρCp, L = across-flats/2 to circumradius), (c) an ANS-5.1 decay
  curve for a loss-of-forced-cooling transient.
- **FEM–CFD cross-check (VR-5):** compare the 3D CFD against the 2D FEM sub-problem.
- **Figures + interview write-up** for the 3D model and its validation.
- **Extensions (not scheduled):** snappyHexMesh inflation layers; passive/LOFC case in
  3D; transient with an ANS-5.1 decay curve.

## Status summary
| Component | State |
|---|---|
| 2D FEM solver + verification (5/5) | complete |
| 3D CHT model (`chtMultiRegionFoam`) | complete |
| 3D mesh quality (hex, near-wall) | complete |
| 3D peak fuel temperature (converged) | complete (670 ± 2 °C, hex mesh) |
| Tri-prism GCI study | complete (p = 1.96, 604 ± 3.5%) |
| Hex GCI study | complete (670 ± 2 °C, grid-independent) |
| Energy balance | complete (0.01% closure, 3 measures) |
| Steady peak fuel temperature | ≈ 667 °C — with physical Cp(T), via the T-based solver |
| Variable-Cp artifact | root-caused (enthalpy solver) + fixed (T-based solver, ADR-19) |
| Temperature-dependent solid Cp(T), k(T) | active in production; Cp(T) only affects transient inertia |
| Custom T-based CHT solver (chtMultiRegionTFoam) | complete (VR-7: varCp → 666.92 °C) |
| Distinct fuel k_compact (8-region) | future work |
| Transient run (Cp(T) inertia) | ready to run on the T-solver; not yet executed |
