# Prismatic Gas-Cooled Microreactor — Thermal Modeling Study

## What this project is

This project **models the temperature inside the fuel block of a helium-cooled,
TRISO-fueled nuclear microreactor** — the class of reactor Radiant is building — to
answer one safety-critical question:

> **Does the fuel stay below its 1600 °C integrity limit — both at full power, and
> if all active (forced) cooling is lost?**

Concretely, it takes a 2D cross-section of a hexagonal graphite fuel block with
interleaved fuel compacts and helium coolant channels, and computes the steady-state
temperature field for two cases:

1. **Normal operation** — full fission power, helium actively cooling the channels.
2. **Loss of forced cooling** — the blower stops; only decay heat is produced and the
   block must shed it *passively*, by thermal **radiation** from its outer surface to
   the vessel.

The physics is solved **two independent ways** — a from-scratch finite-element (FEM)
conduction solver written in Python, and an OpenFOAM finite-volume (CFD) case — so the
two methods cross-check each other.

### What it found

| Case | Peak fuel temperature | Margin to the 1600 °C limit |
|---|---|---|
| Normal operation | **822 °C** | 778 K |
| Loss of forced cooling (passive radiation only) | **360 °C** | 1240 K |

Even with all forced cooling lost, radiation alone holds the fuel ~1240 K below the
limit — a quantified demonstration of the **"walk-away safe"** behavior that defines
this reactor class.

### Why it's built like a formal analysis (not just a script)

It opens with **requirements**, an **architecture**, and a **verification &
validation (V&V) plan**, and only then implements and runs the models — mirroring how
safety-relevant nuclear thermal analysis is actually produced and reviewed. The FEM
solver is verified against exact/manufactured solutions (order of accuracy **p =
1.999**, energy balance closed to **0.000 %**), so the numbers above come with
evidence they're correct, not just plausible.

Along the way it demonstrates **finite element method (FEM)**, **computational fluid
dynamics (CFD)**, and **radiation heat transport** modeling — see the mapping below.

> **Want to run it yourself?** See **[RUNNING.md](RUNNING.md)** for prerequisites,
> installation (Python + OpenFOAM 7), and step-by-step commands. The 2D FEM part runs
> in pure Python with no OpenFOAM required.

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
│   ├── 03_verification_validation.md ← V&V plan: benchmarks, metrics, tolerances
│   ├── 04_decision_log.md        ← ADRs: key technical decisions + rationale
│   ├── 05_outcomes_log.md        ← chronological results / engineering notebook
│   ├── 06_study_guide.md         ← plain-language walkthrough (understand it)
│   └── 07_interview_prep.md      ← likely interview Q&A + honest answers
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
- [x] FEM solver core (assembly + Newton with analytic tangent)
- [x] Verification benchmarks — **5/5 pass** (MMS order p = 1.999)
- [x] Scenario runs + figures
- [ ] OpenFOAM CFD case

## Results snapshot (Part 1 — FEM)

| Scenario | Peak fuel T | TRISO margin |
|---|---|---|
| Normal operation (forced helium cooling) | 822 °C | 778 K |
| Loss of forced cooling (decay heat, passive radiation) | 360 °C | 1240 K |

TRISO integrity limit: 1600 °C. The loss-of-forced-cooling case rejects decay heat
**entirely by radiation** from the block's outer surface and still holds the fuel
~1240 K below the limit — a quantified "walk-away safe" demonstration.

**Verification** (`python fem_thermal/verify.py`):

| Test | Metric | Result | Criterion |
|---|---|---|---|
| Patch test | L2 error | 4.9e-16 | < 1e-10 |
| MMS order of accuracy | observed p | **1.999** | 1.8–2.2 |
| Analytic Jacobian vs FD | rel. diff | 2.1e-8 | < 1e-6 |
| Newton quadratic convergence | — | pass | — |
| Global energy balance | closure | 0.000 % | ≤ 1.0 % |

Figures: `fem_thermal/figures/temperature_fields.png`, `.../convergence.png`.

Reproduce everything:

```bash
python fem_thermal/verify.py   # V&V suite
python fem_thermal/main.py     # scenarios + figures
```

## Environment

Python 3.13, numpy, scipy, matplotlib, gmsh (all present). OpenFOAM runs in WSL2
(Ubuntu). No other dependencies.

## How to read this repo

Start with `docs/01_system_requirements.md`. Each requirement carries an ID that is
referenced by the architecture and by the V&V acceptance criteria, so the model,
the reasons it exists, and the evidence it works are traceable end to end. Then:

- `docs/04_decision_log.md` — *why* the model is built the way it is (ADR-1..13),
  including the two mid-build physics corrections (ADR-11, ADR-12).
- `docs/05_outcomes_log.md` — *what happened*: a dated record of results, the V&V
  metrics achieved, and the peak-temperature / margin numbers.
