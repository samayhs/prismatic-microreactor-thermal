# Mesh Quality Metrics — Unit-Cell CHT Model

**Document:** mesh-quality record for the 3D conjugate-heat-transfer mesh.
Generated with `checkMesh` (OpenFOAM 7). Reproduce with `./check_mesh_quality.sh`
from a set-up case directory.

---

## Why mesh quality matters here

The model is a **conjugate heat-transfer** problem, so the accuracy of the diffusion
(Laplacian) term — heat conduction and the near-wall temperature gradient — depends
directly on mesh quality. The two metrics that matter most:

- **Non-orthogonality** — the angle between the cell-center-to-cell-center vector and
  the shared face normal. High values force large explicit corrections on the
  diffusion term and can reduce accuracy/stability. Desirable: **max < 70°**,
  **average < 30°**.
- **Skewness** — offset between the face center and the cell-center line's crossing
  point. Desirable: **< 4**.
- **Connectivity** — each region must be a single connected mesh, or the conduction
  path is broken. Confirmed: **1 region** per part.

## The improvement: triangular prisms → hexahedra

The original mesh extruded a 2D **triangular** surface into **prisms**. Triangles
wrapping the curved channel/fuel boundaries produced high non-orthogonality
(max ~86°). Recombining the 2D mesh into **quads** and extruding to **hexahedra**
(`Mesh.RecombineAll` + `recombine=True`) dramatically improved orthogonality, and a
distance-based **graded near-wall refinement** resolves the coolant-wall boundary
layer and controls y⁺.

## Metrics — before vs after

| Metric | Tri-prism (before) | **Hex (after)** | Target |
|---|:---:|:---:|:---:|
| **Fluid region** | | | |
| Cell type | 100% prism | **100% hex** | — |
| Non-orthogonality — max | 86.0° | **15.2°** | < 70° |
| Non-orthogonality — avg | 61.8° | **3.5°** | < 30° |
| Max skewness | 0.25 | **0.52** | < 4 |
| Max aspect ratio | 42.7 | **36.6** | (reported) |
| Connected regions | 1 | **1** | 1 |
| **Solid region** | | | |
| Cell type | 100% prism | **100% hex** | — |
| Non-orthogonality — max | 86.0° | **36.9°** | < 70° |
| Non-orthogonality — avg | 57.6° | **10.1°** | < 30° |
| Max skewness | 0.26 | **0.66** | < 4 |
| Max aspect ratio | 47.4 | **27.6** | (reported) |
| Connected regions | 1 | **1** | 1 |

Production hex mesh: **28,292 cells** (fluid 4,312 + solid 23,980), lc = 2.5 mm,
44 axial layers.

## Near-wall resolution (y⁺)

`yPlus` on the coolant-channel wall (fluid side), converged solution:

| | min | avg | max |
|---|:---:|:---:|:---:|
| y⁺ | 7.3 | **18.6** | 28.0 |

y⁺ sits in the buffer-to-log-layer range. OpenFOAM's blended wall functions
(`nutkWallFunction`, `epsilonWallFunction`) handle this; centering y⁺ nearer 30 by
increasing the first-layer thickness is an identified refinement.

## Effect on the result

The higher-quality, near-wall-resolved hex mesh changes the predicted peak fuel
temperature — evidence the mesh matters:

| Mesh | Peak fuel T |
|---|:---:|
| Tri-prism, medium (74.9k) | 581 °C |
| Tri-prism, Richardson-extrapolated | 604 °C |
| **Hex, near-wall-resolved (28.3k)** | **668 °C** |

The hex mesh resolves the steep near-wall thermal gradient that the coarser tri mesh
smeared, giving a higher and more trustworthy peak. A grid-convergence study on the
hex family is the next step to bracket the mesh-independent value.

## Reproduce

```bash
# from a converted + split case directory
source /root/OpenFOAM-7/etc/bashrc
./check_mesh_quality.sh            # writes mesh_quality_metrics.txt
# y+:
chtMultiRegionFoam -postProcess -region fluid -func yPlus -latestTime
```

## Honest limitations

- **Aspect ratio** (28–37) is moderate, from the long axial cells — acceptable for
  an axial-flow channel but worth noting.
- **y⁺ ≈ 19** is slightly below the ideal log-layer range (see above).
- A fully structured inflation layer (thin, wall-normal cells) was not used:
  gmsh's structured boundary-layer field cannot be applied to the conformal internal
  coolant wall (it is adjacent to 3 surfaces after extrusion). The hex recombination
  plus graded refinement achieves most of the benefit; `snappyHexMesh addLayers`
  would be the route to true inflation layers.
