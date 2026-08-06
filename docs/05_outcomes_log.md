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
- Figures generated: `fem_thermal/figures/temperature_fields.png` (per-scenario
  temperature contours with peak-fuel markers) and `.../convergence.png` (Newton
  residual histories showing quadratic drop).

## 2026-08-06 — Documentation consolidation
- Added the decision log (DL-001) capturing ADR-1..13.
- Added this outcomes log (OL-001).
- Updated the README to index the new documents and reflect Part 1 completion.

## Open items / next
- **Part 2 (OpenFOAM):** build the `chtMultiRegionFoam` + `fvDOM` helium-channel
  case in WSL2 and execute the FEM–CFD cross-check (VR-5). *(in progress)*
- **Solution verification:** run the ≥3-mesh convergence / GCI study on peak fuel
  temperature (VR-4) and attach the numerical uncertainty band to the reported
  value.
- **Extensions (not scheduled):** transient time integration with an ANS-5.1 decay
  curve; enclosure/view-factor radiation; 3D with axial coolant enthalpy rise.

## Status summary
| Component | State |
|---|---|
| Requirements / architecture / V&V plan | complete |
| Materials + mesh | complete |
| FEM solver core | complete |
| Verification suite | complete (5/5) |
| Scenario runs + figures | complete |
| Decision & outcomes logs | complete |
| Solution-verification (GCI) study | pending |
| OpenFOAM CFD cross-check | in progress |
