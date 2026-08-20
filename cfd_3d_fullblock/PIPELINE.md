# CFD Pipeline — 3D Full-Block CHT Model (in progress; blocked at region split)

Start-to-finish account of the `cfd_3d_fullblock` model: why it exists, geometry → mesh →
case scaffolding, and the **modeling blocker** that stops it short of a runnable CFD case,
plus the three-way **fork** for getting past it. File paths are relative to
`cfd_3d_fullblock/`.

> **Status: geometry complete & validated; CFD not yet runnable.** The 360 mm block with
> the correct fuel:coolant ratio is built and checked, but `splitMeshRegions` turns the
> ~109 disconnected coolant channels into ~109 separate fluid regions (§4), which is not
> practical to configure or solve as-is. No physics has been run and there are no results
> yet — the model is paused at a modeling-approach decision (§5).

---

## 0. Why this model exists

The `cfd_3d_unitcell` model is verified and validated, but its geometry is **under-cooled**:
one coolant channel serves six fuel compacts (**fuel:coolant 6:1**, coolant area 8.6 %),
whereas a real GA/MHTGR block is ~**2:1** (coolant ~19 %). That makes the unit cell's
absolute peak conservative rather than representative, and it blocks any comparison to a
published benchmark.

**Goal of this model:** keep the hexagon and hole shapes the same, scale the block to a real
**360 mm across-flats** GA/MHTGR-350 size, and populate it with a full triangular lattice of
holes at the correct **~2:1** ratio — so the absolute peak becomes representative and the
model can be **benchmarked** (the OECD/NEA MHTGR-350 benchmark: real block geometry, He
~6.4 MPa, ~259 → 687 °C).

---

## 1. Geometry

A **full prismatic block** (not a unit cell), built to real GA/MHTGR-350 proportions.

| Parameter | Value | Notes |
|---|---|---|
| Across-flats | **360 mm** | circumradius 207.85 mm; hex area ≈ 1123 cm² |
| Lattice | triangular, **18.8 mm** pitch | same GA-class pitch as the unit cell |
| Fuel compacts | r = 6.2 mm (Ø12.4 mm) | same as unit cell |
| Coolant channels | r = 8 mm (Ø16 mm) | same as unit cell |
| Edge margin | ≥ 4 mm graphite web to the boundary | holes filtered to fit fully inside |
| Axial length | 800 mm | same as unit cell |

**Coolant placement — the key idea.** Holes sit on a triangular lattice; **coolant channels
occupy the √3×√3 sublattice** (the sites where lattice index `(i−j) mod 3 == 0`). That
superlattice contains exactly **1/3** of the sites, so the block is automatically **1/3
coolant, 2/3 fuel → fuel:coolant = 2:1** — no hand-placement, and it is symmetric and
non-overlapping (nearest coolant–coolant = √3·pitch = 32.6 mm; every coolant's six nearest
neighbours are fuel, web 4.6 mm).

**Validated result** (`geometry/make_full_block.py`):

| | This block | Real GA / MHTGR-350 |
|---|---|---|
| Coolant channels | **109** | ~108 |
| Fuel compacts | **222** | ~210 |
| Fuel : coolant | **2.04 : 1** | ~2 : 1 |
| Fuel area | 22.6 % | ~23 % |
| Coolant area | 19.3 % | ~19 % |
| Graphite area | 58.1 % | ~58 % |

The area fractions land almost exactly on a real block — this *is* the MHTGR-350 cross-section.

**Construction** is the same gmsh-OCC pipeline as the unit cell: hexagon surface + 331
disks (109 coolant + 222 fuel) → `occ.fragment(...)` (conforming) → extrude.
`fragment` on 331 disks is fast (~0.1 s).

---

## 2. Meshing

Identical **2D-build-then-extrude → hexes** strategy as the unit cell (§2 of
`../cfd_3d_unitcell/PIPELINE.md`), scaled up:
- distance/threshold near-wall grading from **all** coolant walls (`BL_FIRST/lc` ratios
  matched to the unit cell, y⁺ ~ 30 target),
- recombine → quads → extrude to hexes.

**Volume classification** needed one change from the unit cell: the graphite matrix's
centroid sits at the block centre, which **coincides with the central coolant hole**, so it
can't be found by centroid. Fix: pull out graphite as the single **largest-volume** region
first, then classify the remaining small disk-volumes by centroid against the known hole
centres. Patches classified by position exactly as the unit cell (`inlet`/`outlet` = all
channel ends, `outerWall` = the 6 hex faces, `solidEnds` = solid ends); interfaces untagged.

**Validated (coarse test, `--lc 0.006 --nz 8`):** fluid 109, fuel 222, graphite 1;
inlet/outlet 109 each, outerWall 6, ~98.8 k cells. A production `--lc 0.003 --nz 40` mesh
would be ≈ 1 M+ cells (~48× the unit cell's cross-section). *Only the coarse mesh was
generated* — the region blocker (§4) made a fine mesh moot for now.

---

## 3. Case scaffolding

The OpenFOAM case is **patch-name-driven**, so the unit cell's tracked case files carry over
**unchanged**: same region names (`fluid`, `solid`), same patch names, same BCs. `setup_case.sh`
copies `0.orig/`, `constant/`, `system/` from `../cfd_3d_unitcell`, then runs
`gmshToFoam → splitMeshRegions → topoSet`. The **custom `chtMultiRegionTFoam` solver is
reused as-is** (it is geometry-independent; the unit-cell build in
`../cfd_3d_unitcell/solver-chtMultiRegionTFoam/` is untouched). Multiple channels need no BC
changes in principle — the `inlet`/`outlet` patches just contain many faces instead of one,
all at the same fixed inlet condition (identical parallel channels).

---

## 4. The blocker — `splitMeshRegions` explodes the fluid into ~109 regions

All 109 coolant channels are tagged as **one** physical group `fluid`, i.e. one cellZone
containing 109 **geometrically disconnected** blobs. But `splitMeshRegions -cellZones` splits
a disconnected cellZone into **one region per connected piece**:

```
regions created:  fluid, region1, region2, … region108, solid
                  (109 fluid regions, ~256–272 cells each = one channel apiece)  + solid
```

This is the same behaviour documented in ADR-15 (why the six disconnected fuel compacts
can't be their own region) — here it hits the coolant side at full-block scale.

**Why it blocks a runnable case.** A `chtMultiRegion` case needs, **per region**: an entry
in `constant/regionProperties`, a `constant/<region>/` dir (thermo, sources), a `0/<region>/`
field set (T, U, p, k, ε, …), a `system/<region>/` (schemes, solution), and coupled interface
patches. ~110 regions means ~110× that setup and a ~110-region coupled solve — impractical to
configure and very heavy to run. There is no `splitMeshRegions` flag to keep disconnected
same-zone cells as a single region.

**What is *not* the problem:** the geometry, the 2:1 ratio, the meshing, the classification,
the patch naming, and the solver are all correct and verified. The wall is specific to
*resolving* many disconnected channels with this toolchain.

---

## 5. The fork — three ways past the blocker

The coolant has to be handled one of three ways. The choice is a genuine modeling decision
(resolution vs cost vs how faithful to the full block), so the model is paused here.

### Option A — Correct-ratio **unit cell** (1 coolant + 2 fuel)
Build the natural tile of the coolant sublattice — one coolant channel and its two
associated fuel compacts (six neighbours × 1/3 share). **One fluid region.** Directly fixes
the under-cooled ratio, reuses the existing unit-cell case unchanged, cheapest and fastest to
run.
*Trade-off:* it is still a unit cell, so it does not model the full block or capture any
block-scale radial gradient — but for "make the peak representative by fixing the ratio," it
is the minimal correct answer.

### Option B — **Solid-only + convective wall BC** on the full 360 mm block
Keep the full geometry, but **don't resolve the coolant flow.** Mesh the coolant channels as
voids and apply a convective (Robin) BC on the 109 channel walls: `−k∂T/∂n = h(T − T_bulk(z))`,
with `h` from a Dittus-Boelter-class correlation and `T_bulk(z)` from a per-channel axial
energy balance. **Zero fluid regions** — solve only the solid conduction. This is standard
full-block prismatic thermal practice and is cheap and runnable; it is the 3D analog of the
FEM's coolant treatment.
*Trade-off:* the coolant is a boundary condition, not resolved CFD — you lose the turbulent
near-wall film the CHT resolves, and `h`/`T_bulk(z)` become modeling inputs. Requires
re-meshing the channels as voids and a script to write the per-channel wall BCs.

### Option C — **Symmetry sector** with resolved CHT
Model a 30–60° sector of the block with symmetry-plane BCs on the radial cut faces (the same
adiabatic-by-symmetry argument as the outer hex faces). A 30° sector carries ~9–18 coolant
channels → ~9–18 fluid regions. **Keeps resolved turbulent coolant** and stays
benchmark-representative at a fraction of the full-block cost.
*Trade-off:* still a multi-region build/solve (~10–18 regions), and the sector cut must be
placed through graphite (not through holes); more setup than A or B.

**Recommendation.** For the fastest *representative* peak, **A** (fixes the exact complaint
that motivated this model). For a *full-block* answer that is actually runnable, **B**
(uses the geometry already built here). **C** only if resolved turbulent coolant across the
real block is specifically required.

---

## 6. What is done vs. not done

| Stage | Status |
|---|---|
| Geometry generator (`geometry/make_full_block.py`) | ✅ done, parametric (`--lc --nz --af`) |
| Geometry validation (ratio, area fractions, counts) | ✅ done — 2.04:1, matches MHTGR-350 |
| Coarse mesh + gmshToFoam + splitMeshRegions | ✅ done — **revealed the region blocker** |
| Case scaffolding copy + solver reuse | ✅ in place (patch-name-driven) |
| Production mesh (`--lc 0.003`) | ⛔ not generated (moot until the fork is resolved) |
| Multi-region case setup | ⛔ blocked (§4) |
| Physics run / peak prediction | ⛔ not run — **no results yet** |
| Validation (energy balance, grid independence, …) | ⛔ not started |

Next action is a **decision on §5**, after which the chosen path goes through the same
solve → post-process → validate pipeline documented for the unit cell
(`../cfd_3d_unitcell/PIPELINE.md`).
