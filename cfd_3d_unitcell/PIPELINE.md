# CFD Pipeline — 3D Unit-Cell CHT Model (end to end)

Full start-to-finish account of the `cfd_3d_unitcell` model: geometry → mesh →
physics/case setup → solve → post-processing → results. Companion to `MODEL.md`
(design decisions + rationale) and `MESH_QUALITY.md` (mesh metrics). File paths are
relative to `cfd_3d_unitcell/`.

**One-line summary:** a steady conjugate-heat-transfer model of one hexagonal graphite
tile with a central helium coolant channel and six fuel compacts, extruded 0.8 m, that
predicts a **peak fuel temperature ≈ 667 °C** (~930 °C margin to the 1600 °C TRISO
limit), validated by an energy balance that closes to **0.01 %** and by grid
independence.

---

## 0. How to run it

```bash
# 1) mesh  (Windows / Python + gmsh)
cd geometry && python make_unit_cell.py && python make_toposet.py
# 2) solve (WSL, OpenFOAM 7)   -- run.sh auto-builds the custom solver on first run
cd .. && bash run.sh
# 3) visualize (Windows / Python) or open *.foam in ParaView
python visualize_cfd.py . <latestTime>
```
Requirements: OpenFOAM 7, a **space-free** path, and the custom solver
`solver-chtMultiRegionTFoam/` (built once via its `Allwmake`; see §4). The `.msh` and all
run output are gitignored — tracked content is source only (`constant/*`, `system/*`,
`0.orig/*`, scripts, docs).

---

## 1. Geometry

A **unit cell** of a prismatic HTGR block — one coolant channel's worth of material,
justified by the lattice's periodicity (see §3, adiabatic walls).

| Parameter | Value | Notes |
|---|---|---|
| Hexagon circumradius | 30 mm | across-flats ≈ 52 mm; pointy-top (vertices at 0,60,…°) |
| Coolant channel | 1 central, r = 8 mm | Ø16 mm |
| Fuel compacts | 6, r = 6.2 mm | Ø12.4 mm, on a ring at r = 18.8 mm, every 60° |
| Fuel–coolant pitch | 18.8 mm | GA-class |
| Axial length | 800 mm | ≈ GA block height (793 mm) |

**Cross-section area fractions:** fuel 31 %, coolant 8.6 %, graphite 60 %.
*Caveat:* fuel:coolant is **6:1** vs ~2:1 in a real block — the cell is **under-cooled**
(coolant 8.6 % vs ~19 %). Graphite:fuel (~2) is realistic. The absolute peak is therefore
representative-but-conservative, not licensing-grade (a full-ratio rebuild is future work).

**Construction** (`geometry/make_unit_cell.py`, gmsh OpenCASCADE kernel):
1. Build the 2D cross-section: a hexagon surface, one coolant **disk**, six fuel disks.
2. `occ.fragment(...)` — boolean imprint that makes the fluid/fuel/graphite surfaces
   **conforming** (shared edges), so the fluid↔solid interface is a matching mesh.
3. Export `.step` / `.geo_unrolled` for viewing.

---

## 2. Meshing

**Strategy: build in 2D, extrude to 3D hexes** (`make_unit_cell.py`).

1. **Size field / near-wall grading.** A gmsh **Distance** field from the coolant-channel
   wall feeds a **Threshold** field: cell size ramps from `BL_FIRST = 1.3 mm` at the wall
   (first-cell centre near **y⁺ ≈ 30**, valid for high-Re k-ε wall functions) to
   `lc = 3 mm` in the core over a `BL_THICK = 6 mm` band. `MeshSizeMin = BL_FIRST`,
   `MeshSizeMax = lc`; size-from-curvature/points/boundary are disabled so only the field
   controls size. (A structured boundary-layer field can't be used — gmsh's is 2D-only and
   the coolant wall is a conformal internal interface adjacent to 3 surfaces after
   extrusion; snappyHexMesh `addLayers` is the route to true inflation layers.)
2. **Recombine → quads.** `RecombineAll` + blossom recombination → quad-dominant 2D mesh.
3. **Extrude → hexes.** `occ.extrude(..., numElements=[40], recombine=True)` sweeps the
   quads into hexahedra, 40 axial layers.
4. **Classify volumes** by centroid/size: fuel = centroid on a compact centre; of the two
   axis-centred volumes, the smaller is coolant, the larger graphite.
5. **Classify boundary patches** by position: coolant `z=0`→`inlet`, `z=L`→`outlet`,
   lateral hex faces→`outerWall`, solid `z=0/z=L`→`solidEnds`. Fluid/solid interfaces are
   left **untagged** — `splitMeshRegions` creates the coupled patches automatically.
6. **Physical groups:** `fluid` (coolant) and `solid` (graphite + fuel together). Export
   `fuel_zones.json` (centres/radius/length) for the heat-source `topoSet`.

**Mesh quality** (production; see `MESH_QUALITY.md`): 100 % hexahedra; non-orthogonality
max **15° (fluid) / 37° (solid)**, avg 3.5° / 10°; skewness ≤ 0.66; channel-wall y⁺ min
7.3 / **avg 18.6** / max 28. Recombination cut non-orthogonality from ~86° (the earlier
triangular-prism mesh) to 15°/37°.

**Cell count:** production ≈ **28,292 cells** (fluid 4,312 + solid 23,980).
*Note:* the current `--lc 0.003 --nz 40` defaults regenerate ~20.6 k cells; the documented
28.3 k came from a slightly different parameterization — reconcile if you need the headline
number to be exactly reproducible.

**Heat-source zone** (`geometry/make_toposet.py` → `system/solid/topoSetDict`): six
`cylinderToCell` entries (r = 6.2 mm, z 0→0.8 m at the fuel centres) carve the `fuel`
cellSet/cellZone inside the solid.

---

## 3. Physics & case setup

**Two coupled regions:** `fluid` (helium) and `solid` (graphite + fuel cellZone).

### Governing equations
- **Fluid — steady compressible RANS (helium ideal gas):** continuity `∇·(ρU)=0`;
  momentum `∇·(ρU⊗U) = −∇p + ∇·τ_eff`; energy `∇·(ρU h) = ∇·(α_eff ∇h)`; standard
  **k-ε** transport + high-Re wall functions; `ρ = pM/RT`.
- **Solid — temperature form (see §4):** `ρ Cp(T) ∂T/∂t = ∇·(k(T) ∇T) + q‴`;
  steady ⇒ `∇·(k∇T) + q‴ = 0` (Cp-independent).
- **Interface:** conjugate coupling — `T` and heat flux continuity.

### Fluid properties & BCs (`constant/fluid/*`, `0.orig/fluid/*`)
- Helium: 3 MPa, Cp 5193 J/kg·K, μ 3.4e-5 Pa·s, Pr 0.68, molWeight 4.0026.
- Inlet: fixed velocity **(0 0 10) m/s** (Re ≈ 1.2×10⁴), fixed **T = 300 °C (573 K)**;
  k = 0.375, ε = 34.
- Outlet: fixed pressure 3 MPa.
- Channel wall: coupled (`turbulentTemperatureCoupledBaffleMixed`), no-slip;
  nut/alphat/k/ε wall functions.

### Solid properties & BCs (`constant/solid/*`, `0.orig/solid/*`)
- ρ = 1750 kg/m³; **k(T) = 88 − 0.03 T** W/m·K (`polynomial`); **Cp(T) = −281.9 + 3.8125 T
  − 1.6875e-3 T²** J/kg·K (`hPolynomial`). Fuel shares the graphite curves.
- **Volumetric fission heat q‴ = 7 MW/m³** in the `fuel` cellZone. For the T-based solver
  this is an explicit W/m³ source field driven by `constant/solid/heatSource` (`qVol`,
  `zone`); the stock enthalpy solver uses `constant/solid/fvOptions` (`scalarSemiImplicit`
  on `h`) instead.
- **Outer hexagon wall + axial ends: adiabatic** (`zeroGradient`). Justified by unit-cell
  symmetry — neighbouring tiles are identical, so no net heat crosses the shared faces;
  the hex faces are imaginary cuts inside continuous graphite (no radiating surface). This
  makes the **only** heat exit the coolant wall → the steady energy balance is exact.
- Fuel–coolant wall: coupled thermal interface.

**Acceptance criterion:** peak fuel temperature vs the **1600 °C (1873 K)** TRISO limit.

---

## 4. Solver

**Custom `chtMultiRegionTFoam`** (`solver-chtMultiRegionTFoam/`) — a copy of the stock
`chtMultiRegionFoam` whose **solid energy equation is solved in temperature, not
enthalpy** (`solid/solveSolid.H`):

```
ρ Cp(T) ∂T/∂t = ∇·(k(T) ∇T) + q‴        (steadyState ddt ⇒ the ddt term is zero)
```

*Why:* the stock solver solves the solid in enthalpy, `∇·((k/Cp)∇h)`, which inserts Cp
as a **nodal** value in `alpha = k/Cp` but as an **interval-mean** in `Δh = ∫Cp dT`; with
variable `Cp(T)` these don't cancel, biasing the steady peak by a **mesh-persistent
~+18 °C** (684.6 vs the true 667 °C). Steady conduction is Cp-independent, so that bias is
non-physical. Solving in T removes Cp from the steady operator entirely; `Cp(T)` then acts
only in the transient inertia term. Root cause, verification, and the fix are in
`verification/` (`RESULTS.md`) and `docs/04` ADR-19 / `docs/05` (2026-08-19).

**Build:** `bash solver-chtMultiRegionTFoam/Allwmake` — forces **g++-11** (OpenFOAM 7 does
not compile under the system's g++-15); builds into the user's `FOAM_USER_APPBIN`
(nothing in `/opt` is modified). `run.sh` auto-builds it on first run.

**Numerics** (`system/*`): `steadyState` ddt; bounded Gauss upwind convection;
`limited corrected 0.33` Laplacian/snGrad (for residual non-orthogonality); GAMG (p_rgh) +
PBiCGStab (U, h/T, k, ε); relaxation p_rgh 0.7 / U 0.5 / (h,T) 0.7 fluid, 1.0 solid /
k,ε 0.5. **Convergence:** ~4000 pseudo-iterations to residuals ~1×10⁻⁶.

**Pipeline** (`run.sh`): `gmshToFoam` → `splitMeshRegions -cellZones -overwrite` (→ fluid,
solid regions + coupled interface patches) → `topoSet -region solid` (carve fuel zone) →
`chtMultiRegionTFoam`.

---

## 5. Post-processing

- **Peak fuel temperature:** the max solid temperature (the fuel is the only source).
  `run.sh` greps the solver log's `Min/max T` lines (`awk`, solid stream min > 574 K).
- **Energy balance (validation)** — `validation/`, run via
  `chtMultiRegionFoam -postProcess` (the standalone `postProcess` loads only one region,
  so region-tagged function objects can't resolve — see `validation/TROUBLESHOOTING.md`).
  Three independent measures, defined in `controlDict`'s `functions{}` /
  `validation/energyBalance.functions`:
  - `Q_gen  = q‴·V_fuel` (fvOptions/heat-source volume),
  - `Q_wall = ∮ wallHeatFlux` at the coupled interface,
  - `Q_cool = ṁ·cp·ΔT` (inlet/outlet `surfaceFieldValue` FOs).
  `validation/energy_balance.py` cross-checks against a raw-field parser.
- **Fields for visualization:** `postProcess -func writeCellCentres/writeCellVolumes`;
  `touch fluid.foam solid.foam` as ParaView entry points.
- **Plots:** `visualize_cfd.py` (cross-section + axial profile → `figures/`).

---

## 6. Results

**Peak fuel temperature ≈ 667 °C** (steady, grid-independent 670 ± 2 °C), obtained **with
the physical `Cp(T)`/`k(T)`** via the T-based solver — production-mesh value **666.92 °C**.
Margin to the 1600 °C TRISO limit: **~930 °C**.

**Grid convergence (solution verification):**

| Study | Meshes | Peak fuel | Finding |
|---|---|---|---|
| Hexahedral | 14.6k / 28.3k / 56.5k | 669.3 / 668.3 / 672.0 °C | grid-independent, 0.56 % spread; oscillatory (noise floor) → **670 ± 2 °C** |
| Triangular (older) | 48k / 75k / 123k | 573 / 581 / 587 °C | monotone, observed order **p = 1.96**, extrapolated 604 °C |

The ~66 °C hex-vs-tri gap is near-wall resolution (the hex mesh resolves the film ΔT the
coarse tri mesh smeared); the hex value is trusted.

**Energy balance (validation) — closure 0.01 %:**

| Quantity | Value |
|---|---|
| Q_gen = q‴·V_fuel | 3953.8 W |
| Q_wall (interface flux) | 3953.8 W |
| Q_cool = ṁ·cp·ΔT | 3953.3 W |
| Mass conservation (in vs out) | 0.000 % |

**Supporting hand-checks:** helium density 2.52 kg/m³ (ideal-gas); ṁ = 5.047 g/s;
coolant bulk ΔT = 150.8 K (300 → 451 °C, matches Q/(ṁ·cp)); fuel-zone volume 564.8 cm³;
fuel→wall conduction ΔT ≈ 11 °C (cylindrical-conduction hand-calc).

**Variable-Cp artifact (resolved):** the enthalpy solver reported 684.6 °C with `Cp(T)`;
a 4-mesh study showed the +18 °C is a **mesh-persistent, non-physical** solver bias (does
not vanish under refinement), root-caused to the nodal-vs-interval-mean Cp mismatch and
fixed by the T-based solver (§4). The T-solver returns 666.92 °C with `Cp(T)` — matching
the constant-Cp result to 0.01 °C.

---

## 7. Known limitations (see `MODEL.md §8`)

Fuel:coolant ratio 6:1 (under-cooled); fuel shares graphite k (no distinct `k_compact`,
~2 °C); normal operation only; steady only; single 0.8 m block (not a full column);
power density tuned, properties representative; not benchmark-anchored (MHTGR-350 is the
natural target); y⁺ ≈ 18.6 (buffer/log range, not swept). Next step: a `Cp(T)` transient
on the T-solver (VR-8) — relax-to-steady + time-constant check (τ ≈ 30–45 s).

---

## 8. Validation — demonstrating the solution is trustworthy

The deliverable is two things: **predict the peak fuel temperature** (≈ 667 °C, ~930 °C
margin to the 1600 °C TRISO limit) **and demonstrate the prediction can be trusted.** A
converged CFD run always prints *a* number; the work below is what separates that from a
*defensible* number. Two ideas are kept distinct:

- **Verification** — "are we solving the equations right?" (grid convergence, code/solver
  correctness). Internal, math-facing.
- **Validation** — "are we solving the right equations / does it obey the physics?"
  (conservation, physical cross-checks, and ultimately experiment/benchmark).

The model is strongly **verified** and **internally validated** (conservation + physical
hand-checks); it is **not yet benchmark-validated** against experimental data (§8.5). Four
independent legs support the number, each of which would fail loudly if the solution were
wrong.

### 8.1 Energy conservation — the primary internal validation

**Why it's a strong test here.** Every external boundary is adiabatic (`outerWall`,
`solidEnds`) except the coolant wall. So at steady state, by first-law conservation, the
fission heat has exactly one way out, and three quantities that are computed by **entirely
different paths** must be equal:

```
Q_gen  (volumetric source in the fuel)                    [what goes in]
   =  Q_wall (conductive flux across the fluid–solid interface)   [what crosses the wall]
   =  Q_cool (enthalpy rise of the helium, ṁ·cp·ΔT)              [what comes out]
```

These share no intermediate result — `Q_gen` is a cell-volume integral of the source,
`Q_wall` is a surface integral of `k∇T` at the interface, `Q_cool` is an inlet/outlet mass-
and-enthalpy flux. Agreement is therefore a genuine, non-circular check on the whole
coupled solution (mesh, coupling BC, turbulence, thermo).

**Result — closure 0.01 %:**

| Measure | How it is computed | Value |
|---|---|---|
| `Q_gen`  | `q‴·V_fuel` (fvOptions/heat-source cell-volume integral) | 3953.8 W |
| `Q_wall` | `∮ wallHeatFlux` over the coupled interface (`wallHeatFlux` FO) | 3953.8 W |
| `Q_cool` | `ṁ·cp·ΔT` (inlet/outlet `surfaceFieldValue` FOs) | 3953.3 W |
| Mass in vs out | inlet vs outlet mass flux | **0.000 %** |

Spread ≈ 0.5 W on ~3954 W → **0.01 %**. Mass conserving to 0.000 % separately confirms the
compressible flow field is consistent (a mass leak would corrupt `Q_cool`).

**Mechanics & robustness.** The function objects are defined in `controlDict`'s
`functions{}` (see `validation/energyBalance.functions`) and run with the **all-region**
`chtMultiRegionFoam -postProcess` — the standalone `postProcess` loads only one region, so
region-tagged FOs silently fail to resolve (root-caused in `validation/TROUBLESHOOTING.md`).
As a second, independent implementation, `validation/energy_balance.py` re-derives the same
numbers from a **raw-field parser** and matches the FO output (T_out = 723.83 K,
ṁ = 5.047 g/s) — so the post-processing itself is cross-validated, not trusted blindly.

**What a failure would look like:** a bug in the interface coupling would break
`Q_wall = Q_cool`; an under-resolved fuel source or wrong `V_fuel` would break `Q_gen`; a
non-conservative scheme would break all three. None do.

### 8.2 Grid independence — solution verification

The peak must be a property of the *physics*, not of the mesh. A hexahedral
grid-convergence study over a **4× cell-count range** (near-identical near-wall treatment):

| mesh | cells | peak fuel |
|---|---|---|
| coarse | 14.6k | 669.3 °C |
| medium | 28.3k | 668.3 °C |
| fine | 56.5k | 672.0 °C |

The values do not drift monotonically — they **oscillate within ±2 °C** about ~670. That is
the signature of a solution sitting at its **numerical noise floor**: further refinement
buys noise, not a trend, so a formal Richardson extrapolation is not meaningful and the
**oscillation amplitude is reported as the uncertainty → 670 ± 2 °C (0.56 % spread)**.

An earlier **triangular-prism** study (48k/75k/123k → 573/581/587 °C) *was* in the
asymptotic range — observed order **p = 1.96** (≈ the theoretical 2), Richardson-
extrapolated 604 °C — which verifies the discretisation order, but under-predicts because
coarse triangles on the curved wall smear the near-wall film ΔT. The graded hex mesh
resolves that film (the ~66 °C hex-vs-tri difference), so the hex value is the trustworthy
one. Together: the *method* converges at the expected order (tri study), and the
*production answer* is mesh-independent (hex study).

### 8.3 Analytic hand-checks — physical cross-validation

Independent closed-form estimates confirm the CFD isn't internally consistent-but-wrong.
Each uses only inputs and first principles, not the CFD result:

| Quantity | Hand calculation | Hand value | CFD |
|---|---|---|---|
| Helium density | `ρ = pM/RT = 3e6·4.0026e-3 / (8.314·573)` | 2.52 kg/m³ | 2.52 |
| Mass flow | `ṁ = ρ·U·A_channel = 2.52·10·π(0.008)²` | 5.07 g/s | 5.047 |
| Coolant bulk ΔT | `ΔT = Q/(ṁ·cp) = 3953 / (5.047e-3·5193)` | 150.8 K | 150.8 K (out 451 °C) |
| Fuel-zone volume | `6·π(0.0062)²·0.8` (geometric) | 564.8 cm³ | 564.8 |
| Fuel→wall conduction ΔT | cylindrical conduction, per-compact power through the graphite web | ~11 °C | ~11 °C |

The coolant-ΔT check is the tightest: it ties the *outlet temperature* the CFD computes to
the *heat it generates* through a one-line energy balance — and they agree exactly, closing
the loop between §8.1 and the temperature field.

### 8.4 Solver / code verification — the variable-Cp artifact

The strongest verification step was *catching a wrong-but-plausible number.* With
temperature-dependent `Cp(T)`, the enthalpy solver reported **684.6 °C** — which looks like
a reasonable +18 °C physical effect. It was tested against a **physical invariant**: steady
conduction `∇·(k∇T)+q‴=0` contains no Cp, so a converged steady peak **must** be
Cp-independent. A controlled 4-mesh study (`verification/`, isolating
`δ = T(varCp) − T(constCp)` with identical `k(T)`) showed δ ≈ **17–18 °C and
mesh-persistent** — it does *not* vanish under refinement, so it is not truncation but a
genuine formulation bias. It was root-caused by hand in the solver source
(`heSolidThermo::calculate()`: nodal Cp in `alpha = kappa/Cp` vs interval-mean Cp in
`Δh = ∫Cp dT`, which cancel only for constant Cp) and fixed by the temperature-based solver
(§4), which returns **666.92 °C** with `Cp(T)` — matching the constant-Cp result to
**0.01 °C**. This closes the "are we solving the equations right?" question: the reported
number is physics, not a solver artifact. Full write-up: `verification/RESULTS.md`,
`docs/04` ADR-19, `docs/05` (2026-08-19).

### 8.5 What is *not* yet validated (honest scope)

- **No experimental / benchmark comparison.** §8.1–8.4 are self-consistency and code
  verification, not agreement with measured data. The natural external validation is the
  OECD/NEA **MHTGR-350** benchmark (real block geometry, He ~6.4 MPa, ~259 → 687 °C); the
  model is not yet anchored to it.
- **Absolute number is representative, not licensing-grade** — the power density is tuned
  and properties are representative fits, and the unit cell's fuel:coolant ratio is 6:1
  (under-cooled) vs ~2:1 real. So the *method and margin behaviour* are the defensible
  claims; the exact 667 °C is indicative.
- **Steady, normal operation only.** No transient and no loss-of-forced-cooling in 3D.

**Bottom line.** The peak-temperature prediction is **conservation-consistent** (energy
closes to 0.01 % three independent ways), **mesh-converged** (670 ± 2 °C over 4× cells),
**physically cross-checked** (five analytic hand-calcs), and **solver-verified** (a
plausible-looking bias was caught, root-caused, and fixed). What remains for a
licensing-grade claim is external benchmark validation — a stated next step, not a hidden
gap.
