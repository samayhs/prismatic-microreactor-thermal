"""
3D unit-cell mesh for the OpenFOAM conjugate-heat-transfer (CHT) model.

Geometry (prismatic-HTGR unit cell, extruded axially):
  * a hexagonal graphite prism (the moderator associated with one coolant channel),
  * a central helium COOLANT channel (fluid region),
  * a ring of TRISO fuel compacts (solid region with a volumetric heat source).

The 2D cross-section is built with gmsh's OpenCASCADE kernel and fragmented so the
fluid / fuel / graphite areas are conforming (share faces). It is then extruded to a
3D prism mesh. Three physical VOLUMES (fluid, graphite, fuel) become OpenFOAM
cellZones via gmshToFoam; splitMeshRegions then separates them into coupled regions.

External boundary patches:
  inlet      - coolant channel inlet  (z = 0)
  outlet     - coolant channel outlet (z = L)
  outerWall  - hexagon lateral surface (unit-cell symmetry -> adiabatic)
  solidEnds  - top/bottom faces of graphite + fuel (adiabatic)
The fluid<->graphite and fuel<->graphite interfaces are internal faces and are left
untagged; splitMeshRegions creates the coupled interface patches automatically.

Output: unit_cell.msh  (MSH 2.2 ASCII, for gmshToFoam)

Run:  python make_unit_cell.py
"""

from __future__ import annotations

import json
import math
import os

import gmsh

# ---- geometry parameters [m] ----------------------------------------------
R_HEX = 0.030        # hexagon circumradius (center-to-vertex)
R_COOL = 0.008       # coolant channel radius (~15.9 mm dia, GA-class)
R_FUEL = 0.0062      # fuel compact radius (~12.45 mm dia)
R_RING = 0.0188      # fuel-compact ring radius (fuel-coolant pitch)
N_FUEL = 6           # fuel compacts around the channel
LENGTH = 0.80        # axial length (active height) [m]

# ---- mesh controls --------------------------------------------------------
LC = 0.0030          # target in-plane element size [m]
LC_COOL = 0.0015     # finer size on the coolant-channel wall
N_AXIAL = 40         # extrusion layers along z

# ---- boundary-layer controls (fluid side of the coolant-channel wall) ------
# Anisotropic prism/hex layers that resolve the near-wall thermal & momentum
# boundary layer. First-layer thickness is chosen so the first cell centre sits
# at y+ ~ 30 (valid range for the k-epsilon high-Re wall functions); layers grow
# geometrically into the channel core.
BL_FIRST = 0.0013    # first-layer thickness [m] (~y+ 30 for this flow)
BL_RATIO = 1.25      # geometric growth ratio
BL_THICK = 0.0060    # total boundary-layer thickness [m]

# Recombine triangles->quads (2D) and extrude to hexes (3D). Hex/quad cells are
# much more orthogonal than triangular prisms on the curved boundaries, which is
# the dominant source of mesh non-orthogonality here.
RECOMBINE = True


def make_unit_cell(path: str, lc: float = LC, n_axial: int = N_AXIAL,
                   verbose: bool = False, boundary_layer: bool = True):
    gmsh.initialize()
    if not verbose:
        gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("unit_cell")
    occ = gmsh.model.occ

    # --- 2D cross-section ---------------------------------------------------
    hex_pts = []
    for i in range(6):
        a = math.radians(60 * i)
        hex_pts.append(occ.addPoint(R_HEX * math.cos(a), R_HEX * math.sin(a), 0.0))
    hex_lines = [occ.addLine(hex_pts[i], hex_pts[(i + 1) % 6]) for i in range(6)]
    hex_loop = occ.addCurveLoop(hex_lines)
    hex_surf = occ.addPlaneSurface([hex_loop])

    fluid_disk = occ.addDisk(0.0, 0.0, 0.0, R_COOL, R_COOL)

    fuel_centers = []
    fuel_disks = []
    for i in range(N_FUEL):
        a = math.radians(360.0 / N_FUEL * i)
        cx, cy = R_RING * math.cos(a), R_RING * math.sin(a)
        fuel_centers.append((cx, cy))
        fuel_disks.append(occ.addDisk(cx, cy, 0.0, R_FUEL, R_FUEL))

    # conforming split of hex into {fluid, fuel*, graphite}
    occ.fragment([(2, hex_surf)],
                 [(2, fluid_disk)] + [(2, d) for d in fuel_disks])
    occ.synchronize()

    base_surfs = [t for (d, t) in gmsh.model.getEntities(2)]

    # Identify the coolant-channel wall curves (boundary of the fluid surface) so
    # a boundary-layer field can be attached before meshing. The fluid surface is
    # the smallest-area surface centred on the axis.
    coolant_curves = []
    fluid_2d = None
    fmin_area = 1e30
    for (d, t) in gmsh.model.getEntities(2):
        cx, cy, _ = occ.getCenterOfMass(2, t)
        if math.hypot(cx, cy) < R_COOL * 0.5:
            area = occ.getMass(2, t)
            if area < fmin_area:
                fmin_area, fluid_2d = area, t
    if fluid_2d is not None:
        coolant_curves = [abs(t) for (d, t) in
                          gmsh.model.getBoundary([(2, fluid_2d)], oriented=False)]

    # --- extrude to 3D ------------------------------------------------------
    # recombine=True extrudes the (recombined) 2D quads into HEXES, which are far
    # more orthogonal than triangular prisms on the curved channel/fuel boundaries.
    occ.extrude([(2, t) for t in base_surfs], 0.0, 0.0, LENGTH,
                numElements=[n_axial], recombine=RECOMBINE)
    occ.synchronize()

    # --- classify volumes ---------------------------------------------------
    # Fuel volumes: centroid sits on a known fuel-compact center.
    # The remaining two volumes (fluid channel, graphite matrix) BOTH have their
    # centroid on the axis by symmetry, so they cannot be told apart by position --
    # distinguish them by size instead: the coolant channel is far smaller than the
    # graphite matrix.
    fluid_vols, fuel_vols, graph_vols = [], [], []
    non_fuel = []  # (mass, tag)
    for (d, t) in gmsh.model.getEntities(3):
        cx, cy, cz = occ.getCenterOfMass(3, t)
        if any(math.hypot(cx - fx, cy - fy) < R_FUEL for (fx, fy) in fuel_centers):
            fuel_vols.append(t)
        else:
            non_fuel.append((occ.getMass(3, t), t))
    non_fuel.sort()                         # ascending by volume
    fluid_vols = [non_fuel[0][1]]           # smallest = coolant channel
    graph_vols = [t for _, t in non_fuel[1:]]  # remainder = graphite matrix

    # --- classify exterior boundary surfaces -------------------------------
    all_vols = fluid_vols + fuel_vols + graph_vols
    bnd = gmsh.model.getBoundary([(3, t) for t in all_vols],
                                 combined=True, oriented=False)
    tol = 1e-6
    inlet, outlet, outer, ends = [], [], [], []
    fluid_set = set(fluid_vols)
    for (d, t) in bnd:
        t = abs(t)
        cx, cy, cz = occ.getCenterOfMass(2, t)
        up, _ = gmsh.model.getAdjacencies(2, t)
        is_fluid = any(v in fluid_set for v in up)
        if abs(cz) < tol:                      # z = 0 face
            (inlet if is_fluid else ends).append(t)
        elif abs(cz - LENGTH) < tol:           # z = L face
            (outlet if is_fluid else ends).append(t)
        else:                                  # lateral (graphite exterior)
            outer.append(t)

    # --- physical groups ----------------------------------------------------
    def phys(dim, tags, name):
        g = gmsh.model.addPhysicalGroup(dim, tags)
        gmsh.model.setPhysicalName(dim, g, name)

    # Two mesh regions only: fluid (coolant) and solid (graphite + fuel together).
    # The fuel compacts are geometrically disconnected, so making them their own
    # region would fragment the mesh; instead they live inside the single connected
    # SOLID region and get their heat source through a topoSet cellZone (see
    # system/solid/topoSetDict). This is the standard way to embed a volumetric heat
    # source in a CHT solid.
    phys(3, fluid_vols, "fluid")
    phys(3, graph_vols + fuel_vols, "solid")
    phys(2, inlet, "inlet")
    phys(2, outlet, "outlet")
    phys(2, outer, "outerWall")
    phys(2, ends, "solidEnds")

    # --- graded near-wall refinement on the coolant-channel wall -----------
    # gmsh's structured BoundaryLayer field cannot be used here: the coolant wall
    # is a conformal internal interface (shared by the fluid and solid regions), so
    # after extrusion its curves are adjacent to 3 surfaces and gmsh's 2D-only
    # boundary-layer generator rejects it. Instead we grade the in-plane cell size
    # from BL_FIRST at the wall (first cell centre near y+ ~ 30, valid for the
    # k-epsilon wall functions) up to lc in the channel core, using a distance field
    # from the wall surface. This resolves the near-wall thermal/momentum gradient
    # and smooths the size transition that drove the high non-orthogonality.
    # A fully structured inflation layer would require snappyHexMesh addLayers.
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
    gmsh.option.setNumber("Mesh.MeshSizeMin", BL_FIRST)

    # the fluid lateral surface (channel wall): a fluid-volume boundary face whose
    # centroid is off the z = 0 / z = L end planes
    wall_surfs = []
    for (d, t) in gmsh.model.getBoundary([(3, v) for v in fluid_vols],
                                         combined=True, oriented=False):
        _, _, cz = occ.getCenterOfMass(2, abs(t))
        if abs(cz - LENGTH / 2.0) < LENGTH * 0.25:
            wall_surfs.append(abs(t))

    field = gmsh.model.mesh.field
    if boundary_layer and wall_surfs:
        dist = field.add("Distance")
        field.setNumbers(dist, "SurfacesList", wall_surfs)
        field.setNumber(dist, "Sampling", 200)
        thr = field.add("Threshold")
        field.setNumber(thr, "InField", dist)
        field.setNumber(thr, "SizeMin", BL_FIRST)   # fine at the wall
        field.setNumber(thr, "SizeMax", lc)         # coarse in the core / far field
        field.setNumber(thr, "DistMin", 0.0)
        field.setNumber(thr, "DistMax", BL_THICK)
        field.setAsBackgroundMesh(thr)
    else:
        fid = field.add("Ball")
        field.setNumber(fid, "Radius", R_COOL * 1.6)
        field.setNumber(fid, "Thickness", R_COOL * 0.8)
        field.setNumber(fid, "VIn", LC_COOL)
        field.setNumber(fid, "VOut", lc)
        field.setAsBackgroundMesh(fid)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    if RECOMBINE:
        # quad-dominant 2D meshing so the extrusion yields hexes
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)  # blossom

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(3)
    gmsh.write(path)

    # export fuel-zone geometry so topoSetDict can carve the heat-source cellZone
    fuel_zones = {
        "length": LENGTH,
        "r_fuel": R_FUEL,
        "centers": fuel_centers,
    }
    with open(os.path.join(os.path.dirname(path), "fuel_zones.json"), "w") as fh:
        json.dump(fuel_zones, fh, indent=2)

    # report
    ntag, ncoord, _ = gmsh.model.mesh.getNodes()
    n_nodes = len(ntag)
    _, tet_tags, _ = _count_elements()
    print(f"wrote {path}")
    print(f"  volumes : fluid={len(fluid_vols)} graphite={len(graph_vols)} "
          f"fuel={len(fuel_vols)}")
    print(f"  patches : inlet={len(inlet)} outlet={len(outlet)} "
          f"outerWall={len(outer)} solidEnds={len(ends)}")
    print(f"  nodes   : {n_nodes}")
    print(f"  3D cells: {tet_tags}")
    gmsh.finalize()


def _count_elements():
    # helper: total 3D elements across types
    etypes, etags, _ = gmsh.model.mesh.getElements(3)
    total = sum(len(t) for t in etags)
    return etypes, total, None


if __name__ == "__main__":
    import argparse

    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description="Generate the 3D unit-cell mesh.")
    ap.add_argument("--lc", type=float, default=LC,
                    help="target in-plane element size [m]")
    ap.add_argument("--nz", type=int, default=N_AXIAL,
                    help="number of axial extrusion layers")
    ap.add_argument("--out", default=os.path.join(here, "unit_cell.msh"),
                    help="output .msh path")
    args = ap.parse_args()
    make_unit_cell(args.out, lc=args.lc, n_axial=args.nz, verbose=False)
