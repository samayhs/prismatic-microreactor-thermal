# CLAUDE.md — Orientation (read this first)

This repository holds **two separate, independent thermal-modeling projects**. They
share a subject area (a prismatic gas-cooled microreactor fuel block) but use
**different geometries, different tools, and different discretizations**, and are
developed independently. Do not conflate them.

```
prismatic-microreactor-thermal/
├── fem_2d_block/       ← Project 1: 2D FEM, full block, Python
├── cfd_3d_unitcell/    ← Project 2: 3D CFD, unit cell, OpenFOAM 7
├── docs/               ← shared history (requirements, V&V, decision & outcomes logs)
├── README.md           ← two-project overview
└── RUNNING.md          ← how to run (both)
```

---

## Project 1 — `fem_2d_block/` (2D FEM)

- **What:** a from-scratch 2D finite-element solver — nonlinear heat conduction +
  radiation, Newton–Raphson with an analytic tangent — for a **full prismatic block
  cross-section**. Pure Python + gmsh; **no OpenFOAM**.
- **Geometry:** the large multi-channel block (hexagon, circumradius **0.18 m**,
  **1 central + 6 ring coolant channels, 6 fuel compacts**).
- **Status:** complete & verified — MMS order of accuracy p = 1.999, energy balance
  closed, 5/5 verification tests. Two scenarios: normal operation (peak fuel ≈ 822 °C)
  and loss-of-forced-cooling (peak fuel ≈ 360 °C, radiation-cooled).
- **Run:** `cd fem_2d_block && python verify.py && python main.py && python plot_mesh.py`
- **Key files:** `materials.py`, `mesh.py`, `fem.py`, `verify.py`, `main.py`, `plot_mesh.py`.

## Project 2 — `cfd_3d_unitcell/` (3D CFD)

- **What:** a steady **conjugate heat-transfer** model (`chtMultiRegionFoam`, OpenFOAM 7)
  of a **unit cell** of the block, extruded 0.8 m along the coolant flow. Predicts peak
  fuel temperature with a validation suite.
- **Geometry:** a small unit cell (hexagon, circumradius **0.03 m**, **1 coolant channel
  + 6 fuel compacts**) — **different from the FEM block. This geometric difference is the
  main reason the two are separate projects.**
- **Status:** steady model complete & validated — energy balance closes to 0.01%,
  grid-independent. **Steady peak fuel temperature ≈ 667 °C** (with physical `Cp(T)`/
  `k(T)`, via the custom temperature-based solver — see open findings below).
- **Authoritative reference:** `cfd_3d_unitcell/MODEL.md` (design decisions + metrics).
- **Run:** see `RUNNING.md`. Needs OpenFOAM 7, a **space-free path**; runs in place and
  all output is gitignored. Regenerate the mesh first (`python geometry/make_unit_cell.py`).

---

## Current state & open findings — READ before trusting CFD numbers

- **Steady peak ≈ 667 °C — now obtained *with* physical `Cp(T)`, via a custom solver.**
  The stock `chtMultiRegionFoam` solves the solid in **enthalpy** (`∇·(kappa/Cp ∇h)`),
  which mis-inserts a temperature-dependent Cp (nodal Cp in `alpha` vs interval-mean Cp
  in `Δh`; they cancel only for constant Cp). With `Cp(T)` this gave **684.6 °C** — a
  **+18 °C non-physical bias** that a 4-mesh study showed does **NOT** vanish under
  refinement (it holds at ~17 °C; an earlier doc wrongly assumed it was benign
  truncation). Steady conduction is Cp-independent, so it's spurious. **Fixed (ADR-19):**
  the custom **`chtMultiRegionTFoam`** solves the solid directly in **temperature**
  (`ρCp(T) ∂T/∂t = ∇·(k∇T) + q‴`) — Cp leaves the steady operator; with variable `Cp(T)`
  it returns **666.92 °C** (matches constant-Cp to 0.01 °C). This is now the production
  solver (`system/controlDict`); `run.sh` auto-builds it. Source in
  `cfd_3d_unitcell/solver-chtMultiRegionTFoam/` (build: `Allwmake` — forces **g++-11**,
  since OF7 won't compile under the system's g++-15); nothing in `/opt` is modified.
  Investigation + fix verified in `cfd_3d_unitcell/verification/`. `Cp(T)` now affects
  only the transient inertia. See `docs/05_outcomes_log.md` (2026-08-19) and `MODEL.md`.
- **Geometry ratio (CFD):** fuel:coolant is 6:1; a real block is ~2:1 (under-cooled,
  coolant 8.6% vs ~19%). Graphite:fuel ~2 is fine. Fixing needs a 6-fuel/3-coolant
  rebuild → multiple fluid regions. See `MODEL.md`.
- **Distinct fuel material** (lower `k_compact`) needs the fuel as its own region; the
  6 disconnected compacts make that **8 regions** (tested/confirmed). Deferred; ~2 °C.
- **FEM ↔ CFD cross-check was dropped on purpose** — the two use different geometries,
  so a number-for-number comparison is not apples-to-apples.

## Shared docs (`docs/`)
Cross-project history: `01` requirements, `02` architecture, `03` V&V plan, `04`
decision log (ADR 1–13 = FEM/early, 14–18 = CFD), `05` outcomes log (dated, both),
`06` study guide (FEM concepts), `07` interview prep (legacy framing — this is now a
technical model the owner is building and owns, not interview prep). *These predate the
two-project split and still mix both; splitting them per-project is a pending cleanup.*

## Environment & gotchas
- **OpenFOAM 7** at `/opt/OpenFOAM-7`; load with the `of7` alias. Runs as the normal
  user (not root).
- **OpenFOAM rejects paths with spaces** — keep the repo path space-free (this folder,
  `prismatic-microreactor-thermal`, is fine).
- **CFD runs in place** in `cfd_3d_unitcell/`; run output is **gitignored** (tracked
  source = `constant/*`, `system/*`, `0.orig/*`, scripts, docs). The `.msh` is
  gitignored — regenerate it before the first run.
- **Commit messages:** the owner provides them **verbatim**, no `Co-Authored-By` trailer.

## Keeping this durable (not a one-time fix)
**This file is the first thing to read; keep it true.** Whenever the model, results, or
structure change, update this file *and* `docs/05_outcomes_log.md`. Prefer verifying a
result (energy balance, grid convergence) over trusting solver output.

## Pending cleanup (recommended next steps)
1. Split `docs/` per project (the 04/05 logs currently mix FEM and CFD).
2. Move the FEM-specific study guide (`docs/06`) into `fem_2d_block/`, and give each
   project a short project-level README.
3. Confirm no stale "684.6 °C" remains in the docs (should read ≈ 667 °C + the ADR-19
   T-solver fix, not the old "constant-Cp workaround" framing).
4. CFD transient: the T-solver makes `Cp(T)` transients trustworthy — run one and verify
   it relaxes to the 667 °C steady with τ ≈ 30–45 s (see outcomes log VR-8). Optionally
   build the realistic 6-fuel/3-coolant geometry if a representative ratio is needed.
