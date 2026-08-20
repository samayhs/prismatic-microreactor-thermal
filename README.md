# Prismatic Gas-Cooled Microreactor — Thermal Modeling

This repository contains **two separate, independent thermal-modeling projects** of a
prismatic gas-cooled microreactor fuel block. They study the same kind of reactor but
use **different geometries, different tools, and different discretizations**, and are
developed independently.

> **New here (agent or human)? Read [`CLAUDE.md`](CLAUDE.md) first** for a full
> orientation. To run things, see [`RUNNING.md`](RUNNING.md).

---

## The two projects

### 1. [`fem_2d_block/`](fem_2d_block/) — 2D FEM (Python)

A from-scratch 2D finite-element solver — nonlinear heat conduction + radiation,
Newton–Raphson with an analytic tangent — for a **full prismatic block cross-section**
(hexagon, circumradius 0.18 m, 1 central + 6 ring coolant channels, 6 fuel compacts).
Pure Python + gmsh; **no OpenFOAM needed**.

- Verified: MMS order of accuracy **p = 1.999**, energy balance closed, 5/5 tests.
- Scenarios: normal operation (peak fuel ≈ **822 °C**), loss-of-forced-cooling
  (radiation-cooled, peak fuel ≈ **360 °C**).
- Run: `cd fem_2d_block && python verify.py && python main.py`

### 2. [`cfd_3d_unitcell/`](cfd_3d_unitcell/) — 3D CFD (OpenFOAM 7)

A steady **conjugate heat-transfer** model (`chtMultiRegionFoam`) of a **unit cell** of
the block — one coolant channel + 6 fuel compacts in a small hexagon (circumradius
0.03 m), extruded 0.8 m along the coolant flow. Predicts the peak fuel temperature with
a validation suite.

- Validated: energy balance closes to **0.01%**, grid-independent.
- **Steady peak fuel temperature ≈ 667 °C** (~930 °C margin to the 1600 °C TRISO limit).
- Authoritative reference: [`cfd_3d_unitcell/MODEL.md`](cfd_3d_unitcell/MODEL.md).

## Why two projects (not one)

Different **geometries** (full multi-channel block vs a single-channel unit cell),
different **discretizations** (finite element vs finite volume), different **tools**
(Python vs OpenFOAM), different **scope**. A number-for-number FEM↔CFD comparison would
not be apples-to-apples, which is exactly why they are kept separate.

## Documentation

- [`CLAUDE.md`](CLAUDE.md) — orientation, current state, open findings, gotchas.
- [`RUNNING.md`](RUNNING.md) — prerequisites and run steps for both projects.
- [`docs/`](docs/) — shared history: requirements, architecture, V&V plan, decision log
  (ADR 1–13 FEM, 14–18 CFD), outcomes log, study guide.
- Project specifics live with each project (`cfd_3d_unitcell/MODEL.md`,
  `MESH_QUALITY.md`, `validation/`).
