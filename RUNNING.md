# Getting Started — Run It Yourself

This repo has two independent parts:

- **Part 1 — 2D FEM** (`fem_thermal/`): a from-scratch finite-element conduction+radiation
  solver in pure Python. Runs anywhere Python does — **no OpenFOAM needed.**
- **Part 2 — 3D CFD** (`openfoam/unit_cell_cht/`): a conjugate heat-transfer model in
  **OpenFOAM 7**. Needs a Linux environment (or WSL on Windows).

You can run Part 1 on its own if you just want the Python solver. Part 2 needs OpenFOAM.

---

## Prerequisites

### For Part 1 (and all plotting) — Python 3
```bash
pip install numpy scipy matplotlib gmsh
```

### For Part 2 — OpenFOAM 7
The case dictionaries target **OpenFOAM 7** (from openfoam.org). Other major versions
(8+, or the ESI `.com` builds) use different keywords and may need edits.

- **Linux (Ubuntu/Debian):**
  ```bash
  sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
  sudo add-apt-repository http://dl.openfoam.org/ubuntu
  sudo apt update && sudo apt install openfoam7
  ```
  Installs to `/opt/openfoam7`. Load it with `source /opt/openfoam7/etc/bashrc`.
- **Windows:** install **WSL2** with Ubuntu, then follow the Linux steps inside WSL.
- **macOS:** use the OpenFOAM 7 Docker image, or a Linux VM.

### Optional — ParaView
For interactive 3D viewing of the CFD results (the Python visualizer covers the key
figures without it).

### Two hard requirements for Part 2
1. **No spaces in the repo path.** OpenFOAM rejects paths containing spaces. Clone into
   something like `~/projects/prismatic-microreactor-thermal`, not `~/My Projects/...`.
2. **OpenFOAM must be loaded** in your shell (`source .../etc/bashrc`) before running.

---

## Part 1 — 2D FEM (Python, ~30 s)

From the repo root:
```bash
cd fem_thermal
python verify.py     # verification suite -> 5/5 PASS (MMS order of accuracy p = 1.999)
python main.py       # both scenarios -> figures/temperature_fields.png, convergence.png
python plot_mesh.py  # figures/mesh.png
```
Outputs land in `fem_thermal/figures/`. See `docs/06_study_guide.md` for what it does.

---

## Part 2 — 3D CFD (OpenFOAM)

### Step 1 — generate the mesh (Python + gmsh)
The `.msh` is not committed (it's regenerated output), so build it first:
```bash
cd openfoam/unit_cell_cht/geometry
python make_unit_cell.py     # -> unit_cell.msh (+ unit_cell.step for gmsh viewing)
python make_toposet.py       # -> system/solid/topoSetDict (fuel heat-source zone)
```

### Step 2 — run the solver
In a shell with OpenFOAM 7 loaded:
```bash
source /opt/openfoam7/etc/bashrc            # adjust to your install path
cd openfoam/unit_cell_cht
bash run.sh
```
`run.sh` runs `gmshToFoam → splitMeshRegions → topoSet → chtMultiRegionFoam` **in place**,
prints the **peak fuel temperature**, and leaves all output in this folder (gitignored).
Takes a few minutes.

### Step 3 — visualize (Python, no ParaView needed)
```bash
cd openfoam/unit_cell_cht
python visualize_cfd.py      # -> figures/cfd_axial_profile.png, cfd_cross_section.png
```
Or open `solid.foam` / `fluid.foam` in **ParaView**, click *Apply*, and color by `T`.

### Step 4 (optional) — energy balance & mesh quality
```bash
# energy balance (pure OpenFOAM function objects):
cp system/controlDict system/controlDict.bak
cat validation/energyBalance.functions >> system/controlDict
chtMultiRegionFoam -postProcess -latestTime 2>&1 \
  | grep -iE "with volume|integ\(fluid_to_solid\)|weightedAverage\(outlet\)|sum\(inlet\)"
cp system/controlDict.bak system/controlDict

# mesh quality:
./check_mesh_quality.sh
```
Interpretation is in `validation/TROUBLESHOOTING.md` and `MESH_QUALITY.md`.

---

## Expected results

| Output | Value |
|---|---|
| 2D FEM — normal operation peak fuel T | ~822 °C |
| 2D FEM — loss-of-cooling peak fuel T | ~360 °C |
| 3D CFD — peak fuel T (grid-independent) | ~670 °C |
| Energy-balance closure | ~0.01 % |

Numbers may vary slightly with mesh resolution. See `docs/05_outcomes_log.md` for the
full record and the grid-convergence study.

## Repo layout

```
fem_thermal/            2D FEM solver (Python)
openfoam/unit_cell_cht/ 3D CFD case (OpenFOAM 7)
  geometry/             mesh generator (gmsh)
  constant/ system/ 0.orig/   case source (tracked)
  validation/           energy balance, troubleshooting
docs/                   requirements, architecture, V&V, decision & outcomes logs
```

## Troubleshooting

- **`run.sh` says the path has spaces** — move the repo to a space-free path.
- **`run.sh` says OpenFOAM isn't loaded** — `source .../etc/bashrc` first.
- **Function objects fail in the multi-region case** — see
  `openfoam/unit_cell_cht/validation/TROUBLESHOOTING.md` (use `chtMultiRegionFoam
  -postProcess`, not the standalone `postProcess`).
