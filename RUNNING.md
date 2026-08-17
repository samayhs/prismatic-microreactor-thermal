# Running & Visualizing Everything Yourself

Two parts: the **2D FEM** (pure Python, runs on Windows) and the **3D CFD** (OpenFOAM
in WSL, visualized on Windows). Everything is reproducible from the commands below, and
the CFD now runs **in place inside the repo** — generated output stays in the case
folder and is gitignored, so the repo is self-contained (tracked source + ignored
artifacts).

## Prerequisites

- **Windows Python 3** with numpy, scipy, matplotlib, gmsh — runs the FEM, the mesh
  generator, and all plotting.
- **ParaView** (Windows, `C:\Program Files\ParaView…`) — interactive 3D CFD viewer.
- **WSL / Ubuntu** with **OpenFOAM 7** at `/opt/OpenFOAM-7`, loaded with the `of7`
  alias (added to `~/.bashrc`): typing `of7` runs `source /opt/OpenFOAM-7/etc/bashrc`.
- The repo path must be **space-free** (OpenFOAM rejects spaces). This folder is
  `prismatic-microreactor-thermal`, which satisfies that.

---

## Part 1 — 2D FEM solver (Windows, ~30 s)

```bash
cd fem_thermal
python verify.py     # verification suite -> 5/5 PASS (MMS order p = 1.999)
python main.py       # both scenarios -> figures/temperature_fields.png, convergence.png
python plot_mesh.py  # figures/mesh.png
```

Interactive mesh in the gmsh GUI:
```bash
cd fem_thermal
python -c "from mesh import build_block_mesh; build_block_mesh(show=True)"
```

---

## Part 2 — 3D CFD (conjugate heat transfer)

### 2a. (Re)generate the mesh — Windows

The `.msh` is gitignored, so generate it before the first run (and whenever you change
the geometry). Also writes `.step` / `.geo_unrolled` for viewing in gmsh:

```bash
cd openfoam/unit_cell_cht/geometry
python make_unit_cell.py     # hex mesh, near-wall grading -> unit_cell.msh, unit_cell.step
python make_toposet.py       # rebuild the fuel-zone topoSetDict
```

### 2b. Run the solver in place — WSL

```bash
of7                                    # load OpenFOAM
cd .../openfoam/unit_cell_cht          # (use your path; no spaces)
./run.sh
```

`run.sh` runs `gmshToFoam → splitMeshRegions → topoSet → chtMultiRegionFoam` **in this
folder**, prints the **peak fuel temperature**, and leaves all output here (gitignored).
A few minutes on the `/mnt/c` mount.

### 2c. Visualize — Windows (no ParaView needed)

```bash
cd openfoam/unit_cell_cht
python visualize_cfd.py            # reads the in-place run -> figures/cfd_*.png
```

- **cfd_axial_profile.png** — coolant bulk T and peak fuel T along the channel.
- **cfd_cross_section.png** — solid temperature near the outlet (channel = blank hole).

### 2d. Interactive 3D — ParaView (Windows)

1. Open ParaView → **File ▸ Open** → `openfoam/unit_cell_cht/solid.foam` (and/or
   `fluid.foam`).
2. Click **Apply**, set coloring to **T**.
3. Use **Slice** / **Clip** to cut the block; **Glyph** on the fluid `U` for the flow.

### 2e. Energy balance — WSL (pure OpenFOAM)

From the case folder after a run:
```bash
cp system/controlDict system/controlDict.bak
cat validation/energyBalance.functions >> system/controlDict
chtMultiRegionFoam -postProcess -latestTime 2>&1 | grep -iE "with volume|integ\(fluid_to_solid\)|weightedAverage\(outlet\)|sum\(inlet\)"
cp system/controlDict.bak system/controlDict
```
Gives V_fuel, Q_wall, T_out, and mass flow; do the arithmetic per
`validation/TROUBLESHOOTING.md`. (`validation/energy_balance.py .` is a Python
cross-check.)

---

## Mesh quality & convergence

```bash
of7
cd openfoam/unit_cell_cht
./check_mesh_quality.sh        # non-orthogonality, skewness, cell types (see MESH_QUALITY.md)
```

For the grid-convergence study, the refined meshes are `geometry/mesh_H1.msh …
mesh_H3.msh`; run each through the pipeline and compare peak fuel temperatures
(see `docs/05_outcomes_log.md`).

## Troubleshooting

Post-processing function objects misbehaving in the multi-region case? See
`openfoam/unit_cell_cht/validation/TROUBLESHOOTING.md` — the short version: region
function objects must run via `chtMultiRegionFoam -postProcess`, not the standalone
`postProcess` utility.
