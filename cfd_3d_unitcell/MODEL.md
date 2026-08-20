# 3D CFD Model — Design Decisions & Key Metrics

Single-page reference for the OpenFOAM conjugate-heat-transfer (CHT) model of a
prismatic gas-cooled microreactor fuel-block unit cell. Consolidates the geometry,
modeling choices (with rationale), numerics, results to date, and known limitations.

Related: `docs/04_decision_log.md` (ADRs), `docs/05_outcomes_log.md` (dated results),
`MESH_QUALITY.md`, `validation/TROUBLESHOOTING.md`.

---

## 1. What the model is

A steady-state 3D conjugate heat-transfer simulation of a **unit cell** of a prismatic
HTGR fuel block: a hexagonal graphite prism with a central helium coolant channel and
surrounding fuel compacts, extruded along the flow axis. It predicts the **peak fuel
temperature** and its margin to the 1600 °C TRISO limit under normal operation.

Solver: **`chtMultiRegionFoam`** (OpenFOAM 7), steady.

---

## 2. Geometry

| Parameter | Value | Notes |
|---|---|---|
| Hexagon circumradius | 30 mm | across-flats ≈ 52 mm |
| Axial length (extrusion) | **800 mm** | ≈ GA prismatic block height (793 mm) |
| Coolant channel | 1 central, radius 8 mm | Ø16 mm |
| Fuel compacts | 6, radius 6.2 mm | Ø12.4 mm, on a ring at r = 18.8 mm |
| Fuel–coolant pitch | 18.8 mm | GA-class |

**Cross-section composition (area fractions):**

| | Model | Real GA block | Assessment |
|---|---|---|---|
| Fuel | 31% | ~23% | slightly high |
| Coolant | 8.6% | ~19% | **too low** |
| Graphite | 60% | ~58% | good |
| Graphite : fuel | ~2.0 | ~2.5 | roughly OK |
| Fuel : coolant | 6 : 1 | ~2 : 1 | **under-cooled** |

Built with gmsh (OpenCASCADE kernel) in `geometry/make_unit_cell.py`; the 2D
cross-section is fragmented (conforming fluid/fuel/graphite) and extruded to 3D.

---

## 3. Design decisions (with rationale)

| # | Decision | Rationale |
|---|---|---|
| Unit cell, not full block | one central coolant channel + associated fuel, extruded | the block is a periodic lattice; symmetry lets one tile stand in for all channels at a fraction of the cost |
| Adiabatic outer walls | zero heat flux on the hexagon faces + solid ends | unit-cell symmetry: neighbors are identical, so no net heat crosses the mid-plane. Makes the energy balance exact |
| Two mesh regions (fluid + solid) | fuel is a heat-source cellZone *inside* the solid, not its own region | the 6 fuel compacts are geometrically disconnected; separate regions would fragment the mesh (splitMeshRegions makes one region per disconnected piece) |
| Heat source via `fvOptions` | `scalarSemiImplicitSource` on the `fuel` cellZone (topoSet) | standard way to embed a volumetric source |
| Steady solver | custom `chtMultiRegionTFoam`, `steadyState` ddt | only the steady peak is of interest; transient-to-steady is infeasible (fluid CFL vs slow solid). Custom solver solves the **solid in temperature** (not enthalpy) to avoid the variable-Cp bias (§6, §8, ADR-19); fluid side unchanged |
| k-ε turbulence + wall functions | RANS, high-Re wall functions | channel Re ≈ 1.2e4 is turbulent; wall functions are robust and standard for internal forced convection |
| Temperature-dependent solid properties | `hPolynomial` Cp(T) + `polynomial` k(T), shared by graphite and fuel | Cp(T) is required for a credible transient (thermal inertia); a distinct lower `k_compact` for the fuel would need an 8-region split (§8) |
| Hex mesh + graded near-wall | recombine triangles → hexes; distance-field size grading at the wall | hexes cut non-orthogonality from ~86° to ≤37°; gmsh's structured boundary-layer field can't be applied to the conformal internal wall |

---

## 4. Physics & boundary conditions

**Fluid — helium (compressible ideal gas):**
- Pressure 3 MPa · Cp 5193 J/kg·K · μ 3.4e-5 Pa·s · Pr 0.68 · molWeight 4.0026
- Inlet: fixed velocity **(0 0 10) m/s** (Re ≈ 1.2e4), fixed **T = 300 °C (573 K)**
- Outlet: fixed pressure 3 MPa
- Channel wall: coupled (`turbulentTemperatureCoupledBaffleMixed`), no-slip
- Turbulence: k-ε, k = 0.375, ε = 34 at inlet; nut/alphat/k/ε wall functions

**Solid — graphite matrix + fuel (temperature-dependent):**
- ρ 1750 kg/m³
- **k(T)** = 88 − 0.03·T W/m·K (≈ 71 → 60 over the operating band) — `polynomial`
- **Cp(T)** = −281.9 + 3.8125·T − 1.6875e-3·T² J/kg·K (≈ 1350 → 1810) — `hPolynomial`
- Fuel heat source: **q‴ = 7 MW/m³** in the `fuel` cellZone (shares the solid's k(T)/Cp(T))
- Outer hexagon wall + axial ends: adiabatic (`zeroGradient`) — unit-cell symmetry
- Fuel–coolant wall: coupled thermal interface

**Acceptance criterion:** peak fuel temperature vs the **TRISO limit 1600 °C (1873 K)**.

---

## 5. Mesh & numerics

- **Type:** 100% hexahedra (recombined + extruded), graded near the coolant wall.
- **Production size:** **28,292 cells** (fluid 4,312 + solid 23,980).
- **Quality:** non-orthogonality max **15° (fluid) / 37° (solid)**, avg 3.5° / 10°;
  skewness ≤ 0.66; y⁺ on the channel wall min 7.3, **avg 18.6**, max 28.
- **Schemes:** `steadyState` ddt; bounded Gauss upwind convection; `limited corrected
  0.33` Laplacian/snGrad (for non-orthogonality); GAMG (p_rgh) + PBiCGStab (U,h,k,ε).
- **Relaxation:** p_rgh 0.7, U 0.5, h 0.7 (fluid) / 1.0 (solid), k/ε 0.5.
- **Convergence:** ~4000 iterations; residuals ~1e-6.

---

## 6. Key metrics (results to date)

**Peak fuel temperature ≈ 667 °C** (steady; grid-independent 670 ± 2 °C on the finer
study mesh). Obtained **with the physical temperature-dependent `Cp(T)`/`k(T)`**, via the
custom temperature-based solver (below). Margin to the 1600 °C TRISO limit: **~930 °C**.

> **The old 684.6 °C was a solver artifact — now fixed, not just avoided.** The stock
> `chtMultiRegionFoam` solves the solid in **enthalpy** (`∇·(kappa/Cp ∇h)`), which
> mis-inserts a temperature-dependent Cp: at a face, Cp appears as a *nodal* value in
> `alpha = kappa/Cp` but as an *interval-mean* in `h_N−h_P = ∫Cp dT`; the two cancel
> only for constant Cp. With `Cp(T)` this leaves a **mesh-persistent ~+18 °C bias**
> (verified across a ~5× cell range — it does **not** vanish under refinement, contrary
> to an earlier assumption). Since steady conduction `∇·(k∇T)+q‴=0` is Cp-independent,
> the bias is non-physical. **Fix (ADR-19):** the custom `chtMultiRegionTFoam` solves the
> solid directly in temperature — `ρCp(T) ∂T/∂t = ∇·(k∇T) + q‴` — so Cp is absent from
> the steady operator. With variable `Cp(T)` it returns **666.92 °C** (matching the
> constant-Cp result to 0.01 °C). `Cp(T)` is retained and now affects only the transient
> inertia term. See `verification/` and §8.

**Grid-convergence (solution verification):**

| Study | Meshes | Peak fuel T | Finding |
|---|---|---|---|
| Hexahedral | 14.6k / 28.3k / 56.5k | 669.3 / 668.3 / 672.0 °C | grid-independent, 0.56% spread; oscillatory (noise-floor) → **670 ± 2 °C** |
| Triangular (older mesh) | 48k / 75k / 123k | 573 / 581 / 587 °C | monotonic, observed order **p = 1.96**, extrapolated 604 °C |

The ~66 °C hex-vs-tri difference is the near-wall resolution (the hex mesh resolves the
film ΔT the coarse tri mesh smeared); the hex value is the trustworthy one.

**Energy balance (validation) — closure 0.01%:**

| Quantity | Source | Value |
|---|---|---|
| Q_gen = q‴·V_fuel | fvOptions volume | 3953.8 W |
| Q_wall (interface flux) | `wallHeatFlux` FO | 3953.8 W |
| Q_coolant = ṁ·cp·ΔT | surfaceFieldValue FOs | 3953.3 W |
| Mass conservation | inlet vs outlet | 0.000% |

**Supporting numbers / hand-check agreement:**

| Quantity | Value | Cross-check |
|---|---|---|
| Helium density | 2.52 kg/m³ | matches ideal-gas hand-calc |
| Mass flow ṁ | 5.047 g/s | fixed-inlet ρUA |
| Coolant bulk ΔT | 150.8 K (300 → 451 °C) | matches Q/(ṁ·cp) |
| Fuel-zone volume | 564.8 cm³ (7744 cells) | geometric |
| Fuel → channel-wall conduction ΔT | ~11 °C | matches cylindrical-conduction hand-calc |

---

## 7. Reproduce

```bash
# 1) mesh (Windows/Python + gmsh)
cd geometry && python make_unit_cell.py && python make_toposet.py
# 2) solve (WSL, OpenFOAM 7 loaded)
cd .. && bash run.sh
# 3) visualize (Python) or open *.foam in ParaView
python visualize_cfd.py
```
See `../../RUNNING.md` for prerequisites and the full workflow.

---

## 8. Known limitations & assumptions

Honest list of where the model is representative-but-not-exact:

1. **Fuel-to-coolant ratio (6:1) is too high** — a real block is ~2:1, so the cell is
   under-cooled (coolant 8.6% vs ~19%). The graphite:fuel ratio (~2) is reasonable.
   A 6-fuel / 3-coolant rebuild would fix the ratio (adds multiple fluid regions).
2. **Solid properties are temperature-dependent** — Cp(T) and k(T) polynomials
   (representative fits; Butland-Maddison / MHTGR-350 for licensing-grade). The **fuel
   shares the graphite curves** (no distinct `k_compact` ≈ 15–30 W/m·K). Giving the
   compacts their own material needs the 8-region split (the 6 compacts are
   disconnected); the un-captured effect is ~2 °C, and ρCp of a matrix compact ≈
   graphite anyway. Valid Cp(T) range ~300–1100 K; refit for wider transients.
   **Note:** `Cp(T)` is only physically meaningful with the temperature-based solver
   (`chtMultiRegionTFoam`); under the stock enthalpy solver it produces a non-physical
   ~+18 °C steady bias (§6, ADR-19). This is now resolved, not merely worked around.
3. **Normal operation only** — the loss-of-forced-cooling (passive/radiation) case is
   not modeled in 3D.
4. **Steady-state only** — no transient (real accident analysis is transient).
5. **Single block (0.8 m)**, not a full core column (~8 m); coolant heats over one
   block's worth, not the full-height ~430 K rise.
6. **Power density tuned** (q‴ = 7 MW/m³) to a reasonable peak temperature, not derived
   from a specified reactor power. Material properties are representative, not a
   qualified vendor dataset.
7. **Not benchmark-anchored.** The OECD/NEA **MHTGR-350** benchmark (GA block geometry,
   helium ~6.4 MPa, ~259 → 687 °C) is the natural target if literature-matching is
   required.
8. **y⁺ ≈ 18.6** sits in the buffer/log range; blended wall functions handle it, but
   near-wall (y⁺) sensitivity was not swept in the grid study.
