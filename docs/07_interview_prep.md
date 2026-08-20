# Interview Cheat Sheet — Prismatic Microreactor Thermal Model

**Document ID:** IP-001
One page to walk into the room with. Covers **both** projects — the 2D FEM full-block
model and the 3D OpenFOAM CFD unit cell — end to end: geometry, why the geometry
changed, governing equations, assumptions, meshing, mesh-independence, solvers, results,
the variable-Cp solver bug I found and fixed, and future transient work. Interview Q&A
(with honest redirects) is at the end.

---

## 0. 30-second pitch

> "I built two thermal models of a helium-cooled prismatic microreactor fuel block — the
> Radiant/Kaleidos class. Both predict **peak fuel temperature**, the key safety number
> (TRISO fuel must stay under ~1600 °C). The first is a from-scratch **2D finite-element**
> solver in Python for the full block cross-section — nonlinear conduction plus radiation,
> Newton–Raphson with an analytic tangent — verified to a measured order of accuracy of
> 1.999. The second is a **3D conjugate-heat-transfer CFD** model in OpenFOAM of a unit
> cell, with a full validation suite: energy balance closing to 0.01% and grid
> independence. Along the way I found and fixed a real solver bug — OpenFOAM's
> enthalpy-based solid solver biases the steady temperature by ~18 °C when heat capacity
> is temperature-dependent — by writing a custom temperature-based solver."

---

## 1. The two projects at a glance

| | **Project 1 — FEM** (`fem_2d_block/`) | **Project 2 — CFD** (`cfd_3d_unitcell/`) |
|---|---|---|
| Method | 2D finite element, from scratch (Python) | 3D finite volume, OpenFOAM 7 (CHT) |
| Domain | **Full block cross-section** | **Unit cell**, extruded 0.8 m |
| Geometry | hex circumradius **0.18 m**, 1+6 channels, 6 compacts | hex circumradius **0.03 m**, 1 channel, 6 compacts |
| Physics | nonlinear conduction + surface radiation | conjugate conduction + turbulent forced convection |
| Coolant | modeled (Robin BC, heat-transfer coeff.) | **resolved** (helium flow, k-ε RANS) |
| Peak fuel | 822 °C (normal), 360 °C (loss-of-cooling) | **≈ 667 °C** (normal) |
| Status | complete, verified (5/5, MMS p=1.999) | complete, validated (energy balance 0.01%) |

**They are separate models, not one pipeline.** Same subject, different geometry, tools,
and discretization. A number-for-number FEM↔CFD comparison was **dropped on purpose** —
different geometries make it not apples-to-apples.

---

## 2. Geometry — and why it changed from FEM to CFD

**FEM (full block):** the real prismatic block cross-section — a large graphite hexagon
(circumradius 0.18 m) with a **central coolant channel + a ring of 6 coolant channels**
and **6 fuel compacts** interleaved. This is the right domain for the FEM question:
radial conduction across the *whole* block out to the radiating outer surface, which is
what governs the loss-of-cooling margin.

**CFD (unit cell):** a much smaller tile — **one** hexagon (circumradius 0.03 m) with
**one central coolant channel** and its **6 surrounding fuel compacts**, extruded 0.8 m
along the flow.

**Why the change:**
1. **Cost & validatability.** A full 3D block with resolved turbulent helium in every
   channel is huge and hard to validate in the available time. The block is a periodic
   lattice, so **unit-cell symmetry** lets one channel's tile stand in for all of them at
   a fraction of the cost.
2. **The CFD question is different.** CFD's job is to *resolve* the coolant — the
   near-wall film ΔT, the turbulent convection, the axial heat-up — which the FEM only
   approximated with a wall heat-transfer coefficient. That resolution is what you buy
   with CFD, and a unit cell is enough to get it.
3. **Adiabatic outer walls become exact.** Unit-cell symmetry means neighboring tiles are
   identical, so no net heat crosses the outer hex faces → adiabatic BC is *exact*, which
   makes the steady energy balance close exactly (a clean validation lever).

**Consequence / honest caveat:** the unit cell's fuel:coolant area ratio is ~6:1 vs ~2:1
in a real block (under-cooled); graphite:fuel (~2) is realistic. So the CFD absolute peak
is representative, not licensing-grade. A 6-fuel/3-coolant rebuild would fix the ratio.

---

## 3. Governing equations

### FEM (2D, solid only)
Steady nonlinear heat conduction with a volumetric source:

  ∇·(k(T) ∇T) + q‴ = 0     in the block cross-section

Boundary conditions:
- **Coolant walls (Robin/convective):**  −k ∂T/∂n = h (T − T_coolant)
- **Outer surface (radiation, kept fully nonlinear):**  −k ∂T/∂n = εσ (T⁴ − T∞⁴)
  (plus weak convection in the accident case)

Nonlinearity: `k(T)` (graphite conductivity falls with T) **and** the `T⁴` radiation term.
Solved by Newton–Raphson on the residual `R(T)=0` with the **analytic tangent**
`J = ∂R/∂T`, which includes the `dk/dT` conduction contribution and the `4εσT³`
radiation contribution.

### CFD (3D, conjugate heat transfer)
**Fluid — helium, steady compressible RANS:**
- Continuity:  ∇·(ρU) = 0
- Momentum:  ∇·(ρU⊗U) = −∇p + ∇·τ_eff   (τ_eff includes the turbulent stress)
- Energy (enthalpy):  ∇·(ρU h) = ∇·(α_eff ∇h)
- Turbulence: standard **k–ε** transport equations + high-Re **wall functions**
- Equation of state: ideal gas, ρ = pM/(RT)

**Solid — graphite + fuel:**
  ρ Cp(T) ∂T/∂t = ∇·(k(T) ∇T) + q‴    →  steady:  ∇·(k(T) ∇T) + q‴ = 0

**Fluid–solid interface (conjugate coupling):** temperature continuity `T_f = T_s` and
heat-flux continuity `k_f ∂T/∂n|_f = k_s ∂T/∂n|_s` (OpenFOAM
`turbulentTemperatureCoupledBaffleMixed`).

**Key point for the bug story (§9):** the steady solid equation has **no Cp** — Cp only
multiplies the time derivative. So a *steady* peak temperature must be Cp-independent.

---

## 4. Model assumptions & setup

**FEM:**
- 2D cross-section (ADR-1): captures the governing radial conduction; no axial coolant
  heat-up. Absolute peak is indicative.
- P1 (linear) triangular elements; fuel/graphite interfaces node-conforming.
- Coolant is *modeled*, not resolved: one Robin heat-transfer coefficient per channel.
- Radiation to a large isothermal enclosure, view factor = 1 (gray body).
- Two bounding scenarios: **NORMAL** (full power, forced helium, q‴ = 70 MW/m³) and
  **PASSIVE** (loss of forced cooling: decay heat ≈ 3%, `h_coolant = 0`, reject by
  outer-surface radiation only). Setting `h_coolant = 0` is the physically correct
  loss-of-flow BC (an early version wrongly left the channels as a heat sink).

**CFD:**
- Unit cell + **adiabatic outer hex walls and axial ends** (symmetry) → exact energy
  balance.
- **Two mesh regions only** (fluid, solid); the 6 disconnected fuel compacts are a
  **heat-source cellZone** inside the solid (a separate fuel region would fragment into
  one region per compact → 8 regions). Volumetric source via an explicit `W/m³` field
  (T-solver) / `fvOptions` (stock).
- Operating point: helium **3 MPa**, inlet **300 °C (573 K)**, **10 m/s** (Re ≈ 1.2×10⁴,
  turbulent), q‴ = **7 MW/m³** (tuned to a realistic peak, not derived from a rated
  power).
- Solid properties **temperature-dependent**: k(T) = 88 − 0.03 T (≈ 71→60 W/m·K) and
  Cp(T) quadratic (≈ 1350→1810 J/kg·K). Representative graphite fits.
- Normal operation only; steady only (transient is future work).

---

## 5. Meshing strategy

**FEM:** gmsh (OpenCASCADE) builds the hex-with-holes geometry by boolean operations and
triangulates it; the mesh is exported as plain numpy arrays so the FEM core has **no
mesher dependency**. Material and boundary tags assigned by centroid classification.
Production mesh ≈ **3,111 nodes / 5,970 triangles**.

**CFD:** the hard part, and a real decision (ADR-17).
- Build the 2D cross-section with gmsh's OCC kernel and **fragment** it so fluid / fuel /
  graphite areas are **conforming** (share faces), then **extrude** to 3D.
- **Recombine triangles → quads → hexes.** Triangular prisms on the curved channel/fuel
  boundaries gave max non-orthogonality ~86° (borderline — degrades the conduction term).
  Hexes are far more orthogonal.
- **Graded near-wall refinement** via a gmsh **distance field**: cells are fine at the
  coolant wall (first-cell target **y⁺ ≈ 30**, valid for high-Re wall functions) and
  coarsen into the core. A *structured* boundary-layer field couldn't be used — gmsh's is
  2D-only and the coolant wall is a conformal internal interface adjacent to 3 surfaces
  after extrusion; snappyHexMesh `addLayers` is the route to true inflation layers.
- Result: **100% hexahedra**, non-orthogonality **86° → 15°/37°** (fluid/solid), average
  58–62° → 3.5°/10°. Production **28,292 cells** (fluid 4,312 + solid 23,980), y⁺ avg 18.6.

---

## 6. Mesh-independence criteria & results

**Criteria used:** (1) peak fuel temperature invariant under refinement to within the
solver's numerical noise; (2) a formal grid-convergence (Richardson/observed-order) where
convergence is monotone; (3) an independent **energy-balance closure** check on the final
mesh (see §8).

**CFD hex grid-convergence study** — 3 meshes at ~2× cell steps:

| mesh | cells | peak fuel |
|---|---|---|
| coarse | 14.6k | 669.3 °C |
| medium | 28.3k | 668.3 °C |
| fine | 56.5k | 672.0 °C |

→ **grid-independent to ±2 °C (0.56 % over a 4× cell range).** Convergence is oscillatory
(solution at its numerical noise floor), so the oscillation amplitude is the reported
uncertainty rather than a formal extrapolation: **670 ± 2 °C.**

An earlier **triangular-prism** study (48k/75k/123k) was monotone, observed order
**p = 1.96**, Richardson-extrapolated 604 °C — but the hex value (~670) is the trustworthy
one: the graded hex mesh resolves the near-wall film ΔT that the coarse tri mesh smeared
(a ~66 °C difference, all near-wall resolution).

**FEM verification** (the FEM's analog of grid independence) — 5/5 tests:
patch test L2 = 4.9×10⁻¹⁶; **MMS observed order p = 1.999** (target 2); analytic Jacobian
vs finite-difference 2×10⁻⁸; Newton quadratic convergence; global energy balance 0.000%.

---

## 7. Solvers — what was solved, and how

**FEM (Python, from scratch):**
- *What:* the nonlinear steady conduction+radiation system for nodal temperatures.
- *How:* hand-coded P1 element matrices, Robin and radiation boundary integrals
  (2-point Gauss), assembled global residual and **analytic Jacobian**, solved by
  **Newton–Raphson**. Converges in ~4 iterations (residual drops ~12 orders — quadratic).
  k(T) evaluated at the element-mean temperature; the FD Jacobian is retained only as a
  *verification* tool, not the solver.

**CFD (OpenFOAM):**
- *What:* the coupled fluid (RANS momentum/continuity/energy + k-ε) and solid conduction
  fields, with the conjugate interface.
- *How:* steady `chtMultiRegionFoam`-class solver, `steadyState` ddt, SIMPLE-type
  under-relaxation; GAMG for pressure, PBiCGStab for U/energy/k/ε; `limited corrected
  0.33` on the Laplacian/snGrad for the residual non-orthogonality. ~4000 iterations to
  residuals ~1×10⁻⁶.
- **Custom solver (`chtMultiRegionTFoam`)** — see §9 — solves the **solid in temperature**
  instead of enthalpy; the fluid side is unchanged.

---

## 8. Results

**FEM (full block):**
| Scenario | Peak fuel | Peak graphite | TRISO margin |
|---|---|---|---|
| Normal (forced helium) | **822 °C** | 684 °C | 778 K |
| Loss of forced cooling (decay heat, radiation) | **360 °C** | 357 °C | 1240 K |

The loss-of-cooling case is *cooler* — power drops ~30× to decay heat, and a small load
radiating off a large surface equilibrates low. That's the passive-safety story: nothing
has to work and it still doesn't overheat.

**CFD (unit cell):**
- **Steady peak fuel ≈ 667 °C** (grid-independent 670 ± 2 °C), margin ~930 °C to the
  1600 °C TRISO limit. Obtained *with* physical Cp(T)/k(T) via the T-solver (§9).
- **Energy balance closes to 0.01%:** Q_gen = q‴·V_fuel = 3953.8 W, interface wall flux
  Q_wall = 3953.8 W, coolant enthalpy rise Q_cool = ṁ·cp·ΔT = 3953.3 W. Mass conserves to
  0.000% (ṁ = 5.05 g/s). Coolant bulk ΔT = 150.8 K (300 → 451 °C), matching the hand-calc.
- Hand-check agreements: helium density 2.52 kg/m³ (ideal gas), fuel→wall conduction
  ΔT ≈ 11 °C (cylindrical-conduction hand-calc).

---

## 9. The variable-Cp artifact — root cause and the fix (strong story)

**Symptom.** Turning on temperature-dependent solid `Cp(T)` moved the *steady* peak from
**667 °C → 684.6 °C (+18 °C).** But steady conduction `∇·(k∇T)+q‴=0` has **no Cp** — so
that shift is non-physical.

**Root cause (derived by hand from `heSolidThermo::calculate()`).** The stock solver
solves the solid in **enthalpy**: `∇·(α∇h)` with `α = kappa/Cp` and `h = ∫Cp dT`. On a
face the discrete flux is

  F = interp(kappa/Cp) · (h_N − h_P) = interp(kappa/Cp) · **C̄p** · (T_N − T_P)

where **C̄p is the interval-mean** Cp, but the `α` factor uses **nodal** Cp. Those two Cp
averages cancel only when Cp is constant. With `Cp(T)` they don't → a spurious effective
conductivity → the bias.

**Verified it doesn't wash out.** A 4-mesh study (9k→57k solid cells, near-wall included)
tracked `δ = T(varCp) − T(constCp)` with identical k(T). δ held at **16–18 °C** and did
**not** shrink like an O(h²) truncation (observed order ≈ 0.25). So it can't be
mesh-refined away — earlier notes that assumed it was benign were wrong.

**Fix (ADR-19).** Wrote a custom solver, **`chtMultiRegionTFoam`**, that solves the solid
directly in **temperature**: `ρCp(T) ∂T/∂t = ∇·(k∇T) + q‴`. Cp leaves the steady
operator entirely, and now acts *only* in the transient inertia term — where it
physically belongs. The fission heat is added as an explicit `W/m³` source field on the
fuel cellZone (the ρ-weighted `fvOptions` source is dimensioned for the enthalpy
equation). Built into user space with `Allwmake` (forces g++-11 — OpenFOAM 7 won't compile
under the machine's default g++-15); **nothing in the OpenFOAM install is modified**.

**Result:** with **variable Cp(T)**, the T-solver gives **666.92 °C** — matching the
constant-Cp answer to **0.01 °C**. Artifact eliminated by construction, and I keep the
physical Cp(T).

*Why this is a good story:* trusted a physical invariant over solver output, isolated the
term with a controlled study, read the source to pin the mechanism, and fixed it at the
formulation level rather than papering over it.

---

## 10. Future / transient work

- **Transient (the reason Cp(T) exists).** The T-solver makes a `Cp(T)` transient
  trustworthy. Next: run a solid transient and verify (a) it relaxes to the 667 °C steady,
  (b) its time constant matches the lumped/diffusion estimate **τ ≈ 30–45 s**
  (α = k/ρCp, L = across-flats/2 to circumradius), (c) drive it with an **ANS-5.1 decay
  curve** for a loss-of-forced-cooling transient in 3D (the CFD analog of the FEM PASSIVE
  case).
- **Geometry realism:** 6-fuel / 3-coolant rebuild to fix the fuel:coolant ratio.
- **Distinct fuel material:** lower `k_compact` needs the fuel as its own region → the
  8-region split (~2 °C effect).
- **Near-wall fidelity:** snappyHexMesh inflation layers for a true structured
  boundary layer and a y⁺ sensitivity sweep.

---

## 11. Interview Q&A — tiers

**Tier 1 (nail these): physics, story, results**
- *What did you model and why?* Peak fuel temperature and its margin to the 1600 °C TRISO
  limit, two ways — a 2D FEM full block and a 3D CFD unit cell.
- *Why does 1600 °C matter?* TRISO layers retain fission products up to ~1600 °C; peak
  fuel T is the fundamental safety number.
- *Why is the accident case cooler (FEM)?* Power drops ~30× to decay heat; a small load
  radiating off a big surface equilibrates low — passive safety.
- *How do you know it's right?* FEM: MMS order 1.999, energy balance 0%. CFD: energy
  balance closes 0.01% three independent ways, grid-independent 670 ± 2 °C.

**Tier 2 (handle with this sheet): method choices**
- *Why FEM for one, CFD for the other?* FEM to demonstrate the method and get the block-
  wide conduction+radiation margin cheaply; CFD to *resolve* the coolant (film ΔT,
  turbulence) that FEM only approximated with a wall coefficient.
- *Why a unit cell for CFD?* Symmetry + cost + exact adiabatic BC (see §2).
- *Why hex recombination?* Cut non-orthogonality 86° → 15/37° on curved boundaries.
- *FEM vs FVM?* Both discretize the PDE on a mesh; finite volume enforces conservation
  cell-by-cell (why it's standard for flow), FEM uses weighted residuals. I used each
  where it fits.

**Tier 3 (honest redirect — don't bluff):** derivations of the weak form or the Newton
tangent from scratch. *"I understand conceptually what it does and I verified it —
analytic Jacobian vs finite difference at 2×10⁻⁸, quadratic Newton convergence. I leaned
on references/AI for the detailed derivation. What I own is the physics setup, the
verification, and — for the CFD — reading the solver source to find and fix the
variable-Cp bug."* The variable-Cp story (§9) is the strongest thing to steer toward: it
shows real depth.

---

## 12. Ask them / night-before

**Ask them:**
- "For Kaleidos, how much margin do you design to on peak TRISO temperature under the
  bounding loss-of-cooling case?"
- "Is decay-heat rejection conduction-dominated to the vessel, or does radiation carry
  most of it — and how does that change during transport?"

**Night-before checklist:**
- [ ] Say the 30-second pitch out loud 3×.
- [ ] Be able to state each governing equation and its BCs (§3) from memory.
- [ ] Walk the variable-Cp story (§9) in ~90 seconds — symptom, invariant, root cause, fix.
- [ ] Narrate each figure in one sentence.
- [ ] Have `python verify.py && python main.py` (FEM) ready to run live.
