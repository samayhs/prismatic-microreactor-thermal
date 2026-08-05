# Prismatic Gas-Cooled Microreactor — Thermal Modeling Study

A self-contained thermal-modeling portfolio project demonstrating **finite element
method (FEM)**, **computational fluid dynamics (CFD)**, and **radiation heat
transport** on a geometry representative of a helium-cooled, TRISO-fueled
microreactor (the class of reactor Radiant is building).

The project is deliberately structured like a small analysis campaign rather than a
script: it opens with **requirements**, an **architecture**, and a **verification &
validation (V&V) plan**, and only then implements and runs the models. That
ordering mirrors how safety-relevant nuclear thermal analysis is actually produced
and reviewed.

---

## What this demonstrates (mapped to the role)

| Job requirement | Where it is demonstrated |
|---|---|
| Finite element modeling principles & methods | Part 1 — a from-scratch 2D FEM conduction solver (weak form, P1 assembly, Newton–Raphson with analytic tangent) |
| Computational fluid dynamics | Part 2 — an OpenFOAM conjugate heat-transfer case for a helium coolant channel |
| Radiation heat transport modeling | Part 1 — nonlinear σε(T⁴−T∞⁴) boundary condition for passive decay-heat rejection; Part 2 — `fvDOM` discrete-ordinates model |
| Nuclear reactor design knowledge | Prismatic block geometry, TRISO 1600 °C limit, normal-operation vs loss-of-forced-cooling scenarios, decay-heat rejection |

## Repository layout

```
Project For Radiant/
├── README.md                     ← you are here
├── docs/
│   ├── 01_system_requirements.md ← what the models must do + acceptance criteria
│   ├── 02_architecture.md        ← module design, data flow, numerical methods
│   └── 03_verification_validation.md ← V&V plan: benchmarks, metrics, tolerances
├── fem_thermal/                  ← Part 1: Python FEM solver
│   ├── materials.py              ← k(T), scenarios, TRISO limit
│   ├── mesh.py                   ← gmsh prismatic-block mesh → numpy arrays
│   ├── fem.py                    ← (to build) assembly + nonlinear solve
│   ├── verify.py                 ← (to build) code-verification benchmarks
│   └── main.py                   ← (to build) run scenarios, emit figures
└── openfoam/                     ← Part 2: CFD case (to build)
```

## Status

- [x] Requirements, architecture, V&V plan
- [x] Materials model and mesh generator
- [ ] FEM solver core (assembly + Newton)
- [ ] Verification benchmarks
- [ ] Scenario runs + figures
- [ ] OpenFOAM CFD case

## Environment

Python 3.13, numpy, scipy, matplotlib, gmsh (all present). OpenFOAM runs in WSL2
(Ubuntu). No other dependencies.

## How to read this repo

Start with `docs/01_system_requirements.md`. Each requirement carries an ID that is
referenced by the architecture and by the V&V acceptance criteria, so the model,
the reasons it exists, and the evidence it works are traceable end to end.
