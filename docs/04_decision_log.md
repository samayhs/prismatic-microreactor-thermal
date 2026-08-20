# Decision Log

**Document ID:** DL-001  **Rev:** A
Architecture Decision Records (ADRs) for the thermal-model project. Each entry
captures a decision, why it was made, the alternatives weighed, and the
consequences. Newest decisions are appended at the bottom.

Status legend: **Accepted** · **Superseded** · **Proposed**

---

## ADR-1 — 2D cross-section (not 3D) for this revision
**Status:** Accepted
**Context:** The thermal question of interest (peak fuel temperature and margin to
the TRISO limit) is dominated by radial conduction from fuel compacts to coolant
channels and to the outer surface. Full 3D adds axial coolant enthalpy rise and end
effects at large cost.
**Decision:** Model a representative 2D cross-section of the prismatic block.
**Alternatives:** Full 3D block; 1D annular unit cell.
**Consequences:** Fast to build, mesh, and verify; captures the governing radial
physics. Axial coolant heat-up is not resolved and is listed as an extension
(ARC-001 §5). Absolute peak temperature is indicative, not licensing-grade.

## ADR-2 — Write the FEM from first principles (no PDE library)
**Status:** Accepted
**Context:** The deliverable must *demonstrate* FEM method knowledge for an
interview, not just produce a number. Libraries (FEniCS, scikit-fem) would hide the
assembly.
**Decision:** Hand-code element geometry, assembly, boundary integrals, and the
Newton solver in numpy/scipy (requirement SR-5).
**Alternatives:** FEniCS/dolfinx; scikit-fem; MOOSE.
**Consequences:** Every term is defensible line-by-line, including the analytic
tangent. Cost: no built-in higher-order elements or adaptivity. Mitigated by
independent verification (MMS) proving correctness.

## ADR-3 — Linear (P1) triangular elements
**Status:** Accepted
**Context:** Need a simple, robust element with a closed-form conduction matrix for
a hand-written solver.
**Decision:** P1 constant-gradient triangles; fuel/graphite interfaces
node-conforming.
**Alternatives:** P2 (quadratic) triangles; Q1 quads.
**Consequences:** Closed-form element matrices; O(h²) accuracy confirmed by MMS
(p = 1.999). Peak-temperature capture at channel walls needs adequate refinement,
addressed by the mesh-convergence study (VR-4).

## ADR-4 — Treat radiation fully nonlinearly (no T⁴ linearization)
**Status:** Accepted
**Context:** The passive (loss-of-forced-cooling) case rejects decay heat mainly by
radiation; the σε(T⁴−T∞⁴) term drives the result. Linearizing T⁴ about a reference
would bias the very quantity we care about.
**Decision:** Keep the T⁴ term exact and resolve it in the Newton iteration with a
2-point Gauss boundary integral.
**Alternatives:** Linearize about a fixed T_ref (radiation "heat-transfer
coefficient"); Picard/successive-substitution.
**Consequences:** Physically faithful across a wide temperature range; requires an
analytic radiation tangent (4εσT³) — see ADR-5.

## ADR-5 — Analytic Newton tangent (not finite-difference Jacobian)
**Status:** Accepted
**Context:** Nonlinearity comes from k(T) and the T⁴ radiation BC. A correct
tangent gives quadratic convergence; a numerical Jacobian is slow and noisy.
**Decision:** Derive and assemble the analytic Jacobian, including the dk/dT
conduction term and the 4εσT³ radiation term.
**Alternatives:** Finite-difference Jacobian; Picard iteration.
**Consequences:** Quadratic convergence (~4 iterations, residual ↓ ~12 orders).
Correctness is independently checked against a finite-difference Jacobian (V4a, rel.
diff 2e-8) — the FD Jacobian is kept as a *verification* tool, not the solver.

## ADR-6 — Evaluate k(T) at the element-mean temperature
**Status:** Accepted
**Context:** P1 temperature is linear per element; conductivity k(T) must be
reduced to a per-element value for the closed-form matrix.
**Decision:** Use k(T_e) with T_e the mean of the element's nodal temperatures;
propagate dT_e/dT_n = 1/3 into the tangent.
**Alternatives:** Gauss-point evaluation of k; nodal-averaged k.
**Consequences:** Simple, consistent, preserves O(h²) (MMS passes). The tangent
term is exact for this choice, keeping Newton quadratic.

## ADR-7 — Coolant modeled by a wall heat-transfer coefficient (Robin BC)
**Status:** Accepted
**Context:** In the FEM the helium is represented, not resolved. A Robin BC
−k∂T/∂n = h(T−T_c) captures convective removal with one parameter.
**Decision:** Impose Robin BCs on channel walls with h from a Dittus-Boelter-order
estimate; resolve the actual coolant field in the separate OpenFOAM case instead.
**Alternatives:** Conjugate heat transfer in the FEM (couple a fluid region).
**Consequences:** Keeps the FEM a pure conduction/BC problem (clean to verify). The
wall-h vs resolved-flow modeling difference is exactly what the FEM–CFD cross-check
(VR-5) quantifies.

## ADR-8 — Gray-body radiation to a large isothermal enclosure (view factor = 1)
**Status:** Accepted
**Context:** The outer surface radiates to the vessel/containment. Full
enclosure/view-factor radiation is expensive and not needed to demonstrate the
margin behavior.
**Decision:** Model outer-surface radiation as εσ(T⁴−T∞⁴) to a large isothermal
sink (VF = 1).
**Alternatives:** Enclosure radiation with computed view factors; participating
media.
**Consequences:** Bounds the passive rejection with a transparent assumption.
View-factor radiation is an identified extension (ARC-001 §5).

## ADR-9 — gmsh (OpenCASCADE) for meshing
**Status:** Accepted
**Context:** The geometry (hex block with interleaved circular fuel/coolant
channels) needs boolean construction and quality unstructured triangulation.
**Decision:** Build geometry with gmsh's OCC kernel; export plain numpy arrays so
the FEM core has no mesher dependency.
**Alternatives:** Hand-rolled structured mesh; Triangle; pygmsh.
**Consequences:** Clean separation (mesh.py → arrays → fem.py). Professional mesh;
material/boundary tagging via centroid classification.

## ADR-10 — Two bounding scenarios: normal operation and loss of forced cooling
**Status:** Accepted
**Context:** The margin story needs a design-basis normal case and a safety-relevant
off-normal case.
**Decision:** NORMAL (full power + forced helium convection) and PASSIVE (decay
heat + passive rejection).
**Alternatives:** A single nominal case; a full transient sweep.
**Consequences:** Directly answers "margin at power" and "margin if cooling is
lost." Time-accurate transients (ANS-5.1 decay curve, thermal capacitance) are an
identified extension.

## ADR-11 — PASSIVE case models true loss of flow: coolant sink removed
**Status:** Accepted (supersedes the initial PASSIVE parameterization)
**Context:** The first PASSIVE setup left the coolant channels acting as a 600 K
heat sink (h_coolant = 25). That produced fuel temperatures *below* the normal case
and *below* the coolant inlet — physically misleading for a "loss of cooling" event,
because stagnant channels with no flow should not remove heat to a fixed-temperature
sink.
**Decision:** Set h_coolant = 0 in PASSIVE so decay heat leaves only through the
outer surface (weak convection + radiation) to the vessel — the physically correct
loss-of-forced-cooling boundary condition.
**Alternatives:** Keep a small natural-convection h into stagnant helium (but the
gas itself heats up, so a fixed-T sink is unjustified).
**Consequences:** The result becomes the correct "radiation-dominated passive
rejection" case (peak fuel 360 °C, margin 1240 K). This is the change recorded in
the outcomes log (OUT-1, 2026-08-05).

## ADR-12 — Nominal power density tuned to a realistic HTGR fuel-temperature range
**Status:** Accepted (supersedes the initial Q_NOMINAL)
**Context:** The initial Q_NOMINAL = 30 MW/m³ gave a normal-operation peak of
~529 °C — physically valid but low for a reactor "at power" (HTGR fuel typically
600–1250 °C).
**Decision:** Raise Q_NOMINAL to 70 MW/m³ so the normal-operation peak (~822 °C)
sits in the expected range; decay-heat fraction (3%) is applied on top for PASSIVE.
**Alternatives:** Leave at 30 MW/m³; back out q''' from a stated core power and
block count.
**Consequences:** More credible normal-case number and a meatier margin narrative.
Absolute values remain indicative given representative (non-vendor) material data.

## ADR-13 — Independent OpenFOAM (finite-volume) case for cross-verification
**Status:** Accepted
**Context:** The role calls for CFD; a second, independent discretization also
strengthens the result (VR-5).
**Decision:** Build an OpenFOAM `chtMultiRegionFoam` case (helium channel + solid,
`fvDOM` radiation) and compare peak solid temperature and wall flux against the FEM.
**Alternatives:** Commercial CFD (STAR-CCM+/Fluent); no cross-check.
**Consequences:** Demonstrates finite-element vs finite-volume fluency and provides
an independent corroboration path. Runs in WSL2/Ubuntu. (Superseded in scope by
ADR-14: the CFD became a full 3D model per the interview take-home.)

## ADR-14 — Represent the 3D CFD as a prismatic unit cell
**Status:** Accepted
**Context:** The interview take-home asks for a 3D CFD model predicting peak fuel
temperature. Modeling a whole block in 3D is expensive and hard to validate in the
available time.
**Decision:** Model a representative **unit cell** — one hexagonal graphite prism with
a central helium coolant channel and six surrounding fuel compacts, extruded 0.8 m.
**Alternatives:** Full extruded hex block (heavy); single fuel-pin cell (too reduced).
**Consequences:** Standard prismatic-HTGR thermal-hydraulics practice; tractable mesh
and run times; clean to validate. Unit-cell symmetry justifies adiabatic outer walls,
which makes the steady energy balance exact (all heat leaves via the coolant).

## ADR-15 — Fuel as a heat-source cellZone inside a single solid region
**Status:** Accepted
**Context:** The six fuel compacts are geometrically disconnected. Making each a
separate mesh region caused `splitMeshRegions` to spawn one region per compact
(fuel + region1..5), which is unworkable.
**Decision:** Use **two** mesh regions only — `fluid` and `solid` — and apply the
volumetric fission heat via `fvOptions` to a `fuel` cellZone carved inside the solid
by `topoSet`.
**Alternatives:** Separate fuel regions (fragmented mesh); concentric fuel annulus
(physically unfaithful).
**Consequences:** Conformal, connected solid region (verified: 1 region); standard way
to embed a heat source; peak fuel temperature = max solid temperature (fuel is the
only source).

## ADR-16 — Steady chtMultiRegionFoam with k-epsilon wall functions
**Status:** Accepted
**Context:** Only the steady peak temperature is of interest. Transient-to-steady is
infeasible (fluid CFL forces ~1e-4 s steps while the graphite equilibrates over
minutes → millions of steps). Channel Re ~ 1.2e4 is turbulent.
**Decision:** Solve steady (`steadyState` ddt, SIMPLE-type relaxation) with k-epsilon
+ wall functions; helium as a compressible ideal gas.
**Alternatives:** Transient run to steady (too slow); LTS/localEuler (more setup);
laminar (wrong — would grossly over-predict the film ΔT); low-Re turbulence (finer
near-wall, but y⁺ floor issues with wall functions).
**Consequences:** Converges in ~500–1000 iterations. Requires a wall-function mesh
(first cell at y⁺ ~ 30), which sets the near-wall meshing strategy (ADR-17) and how
the grid-refinement study is run (refine core/axial, hold the wall-normal first cell).

## ADR-17 — Hex recombination + graded near-wall refinement (no structured BL)
**Status:** Accepted (supersedes the intended structured boundary layer)
**Context:** The extruded triangular-prism mesh had max non-orthogonality ~86°
(borderline), which degrades the diffusion (conduction) term. A structured
boundary-layer mesh was the intended fix.
**Decision:** gmsh's structured BoundaryLayer field is **2D-only and cannot be applied
to the conformal internal coolant wall** (adjacent to 3 surfaces after extrusion), so
instead: **recombine** the mesh to hexahedra (`Mesh.RecombineAll` + `recombine=True`)
and add **distance-based graded near-wall refinement** targeting y⁺ ~ 30.
**Alternatives:** snappyHexMesh `addLayers` (true inflation layers, larger effort);
accept the tri-prism mesh with correctors.
**Consequences:** Max non-orthogonality 86° → 15°/37° (fluid/solid), average 58–62° →
3.5°/10°, 100% hexahedra (see `MESH_QUALITY.md`). Peak fuel prediction shifted 581 →
668 °C, showing the near-wall resolution matters. y⁺ ≈ 18.6 (buffer/log range);
snappyHexMesh inflation layers remain the route to a true structured near-wall mesh.

## ADR-18 — Temperature-dependent solid properties (Cp(T), k(T)); fuel stays a cellZone
**Status:** Accepted
**Context:** Transient analysis needs varying thermal capacity — the graphite thermal
inertia ρCp(T) governs the transient response and graphite Cp varies ~2–3× over the
temperature range (constant Cp would mis-predict the transient). Steady results are
unaffected (Cp does not enter steady conduction).
**Decision:** Give the solid region **temperature-dependent** properties via OpenFOAM
`hPolynomial` (Cp(T)) + `polynomial` transport (kappa(T)), replacing the constant
`hConst`/`constIso`. Keep the fuel as a heat-source cellZone inside the single solid
region (it shares the solid's Cp(T)/k(T)).
**Alternatives:** Constant properties (blocks credible transients). A distinct, lower
`k_compact` for the fuel — which requires promoting the fuel to its own solid region.
**Consequences:** Solid Cp and k now vary with T (representative graphite fits;
Butland-Maddison / MHTGR-350 for licensing-grade). The steady peak temperature is
essentially unchanged (k(T) ≈ 60–71 over the operating band vs the previous constant
60). Verified: OpenFOAM parses the polynomials and the case still converges.
**Why the fuel is not its own region:** the 6 fuel compacts are geometrically
disconnected, so `splitMeshRegions` produces one region per compact — 8 regions total
(tested and confirmed). ρCp of a graphite-matrix compact is close to graphite, so the
shared curves are a fair approximation for thermal inertia; the only distinct-fuel
effect not captured is the lower k_compact, worth ~2 °C. The 8-region split (or a
single-compact / connected-fuel geometry) remains the route to distinct fuel material.

## ADR-19 — Custom temperature-based solid solver (chtMultiRegionTFoam)
**Status:** Accepted (resolves the ADR-18 variable-Cp artifact)
**Context:** With ADR-18's `Cp(T)`, the steady peak fuel temperature rose from the
constant-property **667 °C** to **684.6 °C** — physically impossible, since steady
conduction `∇·(k∇T)+q‴=0` contains **no Cp**. Root cause, traced by hand to
`heSolidThermo::calculate()`: the stock solid solver is **enthalpy-based**, solving
`∇·(α∇h)` with `α = kappa/Cp` and `h = ∫Cp dT`. On a face the flux is
`interp(kappa/Cp)·(h_N−h_P) = interp(kappa/Cp)·C̄p·(T_N−T_P)`, where `C̄p` is the
*interval-mean* Cp but the `α` factor uses *nodal* Cp. The two Cp averages cancel only
when Cp is constant; with `Cp(T)` they leave a bias. A four-mesh study (9k→57k solid
cells, near-wall included) showed the bias is **not** a vanishing O(h²) truncation —
it holds at ~17–18 °C and does not clear at any affordable resolution. So `Cp(T)` was
unsafe for both steady and the transient it was meant to enable.
**Decision:** Build a custom solver — a copy of `chtMultiRegionFoam` whose solid energy
equation is solved directly in **temperature**:
`ρCp(T) ∂T/∂t = ∇·(k(T)∇T) + q‴` (`solveSolid.H`). The fluid side is untouched.
The volumetric fission heat is applied as an explicit `W/m³` source field on the `fuel`
cellZone (`constant/solid/heatSource`), because the ρ-weighted `fvOptions` source is
dimensioned for the enthalpy equation. Built into the user bin via
`solver-chtMultiRegionTFoam/Allwmake` (which forces g++-11: OpenFOAM 7 will not compile
under the machine's default g++-15).
**Alternatives:** (a) constant band-representative Cp ≈ 1626 J/kg·K in the stock solver
— steady-correct by construction, transient inertia good to ~15 %, zero code, but drops
true Cp(T); (b) mesh-refine the artifact away — disproven, it does not vanish; (c)
`laplacianFoam` — temperature-based but constant D_T, no k(T)/source/coupling, unusable.
**Consequences:** With variable `Cp(T)`, the T-solver returns **666.92 °C** on the
production mesh — matching the constant-Cp enthalpy result (666.91 °C) to 0.01 °C. The
artifact is eliminated **by construction**: Cp no longer enters the steady operator, and
`Cp(T)` acts only in the transient inertia term where it belongs. Cost: a maintained
custom solver (rebuild via `Allwmake` on any OpenFOAM change) and the g++-11 build
requirement. The stock `h`-based files (`constant/solid/fvOptions`, the `h` entries in
`fvSolution`) are retained so the enthalpy solver still runs for comparison. Verified in
`cfd_3d_unitcell/verification/` (`RESULTS.md`, `results_uniform.csv`).
