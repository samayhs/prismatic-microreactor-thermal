"""
3D full-block mesh for the OpenFOAM CHT model — MHTGR-350-class prismatic block.

Same idea as the unit cell (`cfd_3d_unitcell`), scaled up to a REAL block so the
fuel:coolant ratio is representative instead of the unit cell's under-cooled 6:1:

  * hexagonal graphite prism, **360 mm across flats** (GA/MHTGR-350 block size),
  * a triangular lattice of holes at **18.8 mm pitch** (same hole sizes as the cell),
  * coolant channels on the sqrt(3)xsqrt(3) sublattice -> exactly 1/3 of sites are
    coolant, 2/3 fuel  => fuel:coolant = 2:1 (area fractions ~ fuel 23% / coolant 19%
    / graphite 58%, matching a real block),
  * extruded 0.8 m along the flow axis.

Everything else matches the unit-cell case so the SAME OpenFOAM case files
(0.orig, constant, system) work unchanged — same region names (fluid, solid) and
patch names (inlet, outlet, outerWall, solidEnds); the fluid<->solid interface is
left untagged for splitMeshRegions to couple. All coolant channels share ONE
physical group "fluid" (one region, disconnected) so the case stays two-region.

Output: full_block.msh (+ fuel_zones.json for topoSet).

Run:  python make_full_block.py [--lc 0.003] [--nz 40] [--af 0.360]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time

import gmsh

# ---- geometry [m] ----------------------------------------------------------
AF_DEFAULT = 0.360        # across-flats (GA/MHTGR-350 block)
PITCH = 0.0188            # triangular lattice pitch (fuel-coolant), GA-class
R_FUEL = 0.0062           # fuel compact radius (Ø12.4 mm)
R_COOL = 0.008            # coolant channel radius (Ø16 mm)
MARGIN = 0.004            # min graphite web from a hole edge to the hex boundary [m]
LENGTH = 0.80             # axial length

# ---- mesh controls (match the unit cell's ratios) --------------------------
LC = 0.0030
N_AXIAL = 40
BL_FIRST_RATIO = 0.0013 / 0.0030    # near-wall first size / lc
BL_THICK_RATIO = 0.0060 / 0.0030    # grading band / lc


def hex_apothem(R):
    return R * math.cos(math.pi / 6.0)


def inside_hex(x, y, R, shrink):
    """Pointy-top hexagon (vertices at 0,60,...): edge normals at 30,90,...°.
    Returns True if (x,y) is inside the hexagon shrunk inward by `shrink`."""
    a = hex_apothem(R) - shrink
    for k in range(6):
        th = math.radians(30.0 + 60.0 * k)
        if x * math.cos(th) + y * math.sin(th) > a:
            return False
    return True


def lattice_sites(R, pitch):
    """Triangular lattice sites (i, j, x, y) covering the hexagon bounding box."""
    a1 = (pitch, 0.0)
    a2 = (pitch * 0.5, pitch * math.sqrt(3.0) / 2.0)
    n = int(2.0 * R / pitch) + 3
    out = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            x = i * a1[0] + j * a2[0]
            y = i * a1[1] + j * a2[1]
            out.append((i, j, x, y))
    return out


def make_full_block(path, lc=LC, n_axial=N_AXIAL, af=AF_DEFAULT, verbose=False):
    t0 = time.time()
    R_HEX = af / (2.0 * math.cos(math.pi / 6.0))
    bl_first = BL_FIRST_RATIO * lc
    bl_thick = BL_THICK_RATIO * lc

    # --- decide hole positions/types on the lattice -------------------------
    cool_centers, fuel_centers = [], []
    for (i, j, x, y) in lattice_sites(R_HEX, PITCH):
        is_cool = ((i - j) % 3 == 0)          # sqrt(3) sublattice -> 1/3 coolant
        r = R_COOL if is_cool else R_FUEL
        if not inside_hex(x, y, R_HEX, r + MARGIN):
            continue
        (cool_centers if is_cool else fuel_centers).append((x, y))

    gmsh.initialize()
    if not verbose:
        gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("full_block")
    occ = gmsh.model.occ

    # --- 2D cross-section ---------------------------------------------------
    hex_pts = [occ.addPoint(R_HEX * math.cos(math.radians(60 * k)),
                            R_HEX * math.sin(math.radians(60 * k)), 0.0)
               for k in range(6)]
    hex_lines = [occ.addLine(hex_pts[k], hex_pts[(k + 1) % 6]) for k in range(6)]
    hex_surf = occ.addPlaneSurface([occ.addCurveLoop(hex_lines)])

    cool_disks = [occ.addDisk(cx, cy, 0.0, R_COOL, R_COOL) for (cx, cy) in cool_centers]
    fuel_disks = [occ.addDisk(cx, cy, 0.0, R_FUEL, R_FUEL) for (cx, cy) in fuel_centers]

    print(f"  holes: coolant={len(cool_disks)} fuel={len(fuel_disks)} "
          f"ratio fuel:cool={len(fuel_disks)/max(1,len(cool_disks)):.2f}")
    print(f"  fragmenting {len(cool_disks)+len(fuel_disks)} disks into the block ...")
    occ.fragment([(2, hex_surf)],
                 [(2, d) for d in cool_disks] + [(2, d) for d in fuel_disks])
    occ.synchronize()
    print(f"  fragment done ({time.time()-t0:.1f}s)")

    base_surfs = [t for (d, t) in gmsh.model.getEntities(2)]

    # --- extrude to 3D (recombine -> hexes) --------------------------------
    occ.extrude([(2, t) for t in base_surfs], 0.0, 0.0, LENGTH,
                numElements=[n_axial], recombine=True)
    occ.synchronize()
    print(f"  extruded ({time.time()-t0:.1f}s)")

    # --- classify volumes ---------------------------------------------------
    # The graphite matrix is the single LARGEST volume, and its centroid sits at
    # the block centre — which coincides with the central coolant hole — so it
    # must be pulled out by SIZE first, then the small disk-volumes classified by
    # centroid against the known hole centres.
    def near(cx, cy, centers, r):
        return any(math.hypot(cx - px, cy - py) < r for (px, py) in centers)

    vols = [(occ.getMass(3, t), t, occ.getCenterOfMass(3, t))
            for (d, t) in gmsh.model.getEntities(3)]
    vols.sort(key=lambda v: v[0], reverse=True)

    fluid_vols, fuel_vols = [], []
    graph_vols = [vols[0][1]]                       # largest = graphite matrix
    for (m, t, (cx, cy, cz)) in vols[1:]:
        if near(cx, cy, cool_centers, R_COOL * 0.9):
            fluid_vols.append(t)
        elif near(cx, cy, fuel_centers, R_FUEL * 0.9):
            fuel_vols.append(t)
        else:
            graph_vols.append(t)                    # fallback (unexpected)

    # --- classify exterior boundary surfaces --------------------------------
    all_vols = fluid_vols + fuel_vols + graph_vols
    bnd = gmsh.model.getBoundary([(3, t) for t in all_vols], combined=True,
                                 oriented=False)
    fluid_set = set(fluid_vols)
    tol = 1e-6
    inlet, outlet, outer, ends = [], [], [], []
    for (d, t) in bnd:
        t = abs(t)
        cx, cy, cz = occ.getCenterOfMass(2, t)
        up, _ = gmsh.model.getAdjacencies(2, t)
        is_fluid = any(v in fluid_set for v in up)
        if abs(cz) < tol:
            (inlet if is_fluid else ends).append(t)
        elif abs(cz - LENGTH) < tol:
            (outlet if is_fluid else ends).append(t)
        else:
            if not is_fluid:            # lateral graphite exterior = hex faces
                outer.append(t)
            # lateral fluid faces are internal (channel walls) -> untagged

    # --- physical groups ----------------------------------------------------
    def phys(dim, tags, name):
        g = gmsh.model.addPhysicalGroup(dim, tags)
        gmsh.model.setPhysicalName(dim, g, name)

    phys(3, fluid_vols, "fluid")                 # ONE fluid group (disconnected)
    phys(3, graph_vols + fuel_vols, "solid")     # graphite + fuel = one solid region
    phys(2, inlet, "inlet")
    phys(2, outlet, "outlet")
    phys(2, outer, "outerWall")
    phys(2, ends, "solidEnds")

    # --- near-wall graded refinement on ALL coolant walls -------------------
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
    gmsh.option.setNumber("Mesh.MeshSizeMin", bl_first)

    wall_surfs = []
    for (d, t) in gmsh.model.getBoundary([(3, v) for v in fluid_vols],
                                         combined=True, oriented=False):
        _, _, cz = occ.getCenterOfMass(2, abs(t))
        if abs(cz - LENGTH / 2.0) < LENGTH * 0.25:
            wall_surfs.append(abs(t))

    field = gmsh.model.mesh.field
    dist = field.add("Distance")
    field.setNumbers(dist, "SurfacesList", wall_surfs)
    field.setNumber(dist, "Sampling", 100)
    thr = field.add("Threshold")
    field.setNumber(thr, "InField", dist)
    field.setNumber(thr, "SizeMin", bl_first)
    field.setNumber(thr, "SizeMax", lc)
    field.setNumber(thr, "DistMin", 0.0)
    field.setNumber(thr, "DistMax", bl_thick)
    field.setAsBackgroundMesh(thr)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)

    print(f"  meshing ({time.time()-t0:.1f}s) ...")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(3)
    gmsh.write(path)

    # fuel-zone export for topoSet
    with open(os.path.join(os.path.dirname(path), "fuel_zones.json"), "w") as fh:
        json.dump({"length": LENGTH, "r_fuel": R_FUEL, "centers": fuel_centers}, fh,
                  indent=2)

    ntag, _, _ = gmsh.model.mesh.getNodes()
    etypes, etags, _ = gmsh.model.mesh.getElements(3)
    ncells = sum(len(t) for t in etags)
    print(f"wrote {path}")
    print(f"  across-flats {af*1e3:.0f} mm | circumradius {R_HEX*1e3:.1f} mm")
    print(f"  volumes: fluid={len(fluid_vols)} fuel={len(fuel_vols)} "
          f"graphite={len(graph_vols)}")
    print(f"  patches: inlet={len(inlet)} outlet={len(outlet)} "
          f"outerWall={len(outer)} solidEnds={len(ends)}")
    print(f"  nodes={len(ntag)}  3D cells={ncells}  ({time.time()-t0:.1f}s total)")
    gmsh.finalize()


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description="Generate the 3D full-block mesh.")
    ap.add_argument("--lc", type=float, default=LC)
    ap.add_argument("--nz", type=int, default=N_AXIAL)
    ap.add_argument("--af", type=float, default=AF_DEFAULT)
    ap.add_argument("--out", default=os.path.join(here, "full_block.msh"))
    args = ap.parse_args()
    make_full_block(args.out, lc=args.lc, n_axial=args.nz, af=args.af)
