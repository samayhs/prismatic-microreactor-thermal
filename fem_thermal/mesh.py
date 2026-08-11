"""
Mesh generation for a prismatic gas-cooled reactor fuel block (2D cross-section).

Geometry: a hexagonal nuclear-graphite block containing, on a single ring,
interleaved TRISO fuel compacts (heat-generating solid subdomains) and helium
coolant channels (voids with a convective wall BC), plus a central coolant
channel. This is the classic prismatic-HTGR arrangement, scaled to a compact
microreactor block.

We use gmsh (OpenCASCADE kernel) to build the geometry with boolean operations
and mesh it with linear (P1) triangles, then hand back plain numpy arrays so the
FEM assembly in fem.py has no gmsh dependency:

    nodes      (Nn, 2)  float   node coordinates [m]
    tris       (Ne, 3)  int     triangle connectivity (0-based node indices)
    tri_mat    (Ne,)    int     material tag per triangle (MAT_GRAPHITE / MAT_FUEL)
    edges_cool (Nc, 2)  int     boundary edges on coolant-channel walls
    edges_out  (No, 2)  int     boundary edges on the outer block surface
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import gmsh
import numpy as np

from materials import MAT_FUEL, MAT_GRAPHITE


@dataclass
class BlockMesh:
    nodes: np.ndarray
    tris: np.ndarray
    tri_mat: np.ndarray
    edges_cool: np.ndarray
    edges_out: np.ndarray

    @property
    def n_nodes(self) -> int:
        return self.nodes.shape[0]


def build_block_mesh(
    circumradius: float = 0.18,   # hexagon center-to-vertex [m]
    ring_radius: float = 0.095,   # radius of the interleaved channel ring [m]
    r_fuel: float = 0.014,        # fuel compact radius [m]
    r_cool: float = 0.011,        # coolant channel radius [m]
    n_pairs: int = 6,             # fuel/coolant pairs around the ring
    lc: float = 0.006,            # target mesh size [m]
    verbose: bool = False,
    show: bool = False,           # open the interactive gmsh viewer before returning
) -> BlockMesh:
    gmsh.initialize()
    if not verbose:
        gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("prismatic_block")
    occ = gmsh.model.occ

    # --- hexagonal outer boundary (flat-top orientation) -------------------
    hex_pts = []
    for i in range(6):
        ang = math.radians(60 * i)
        hex_pts.append(occ.addPoint(circumradius * math.cos(ang),
                                    circumradius * math.sin(ang), 0.0, lc))
    hex_lines = [occ.addLine(hex_pts[i], hex_pts[(i + 1) % 6]) for i in range(6)]
    hex_loop = occ.addCurveLoop(hex_lines)
    hex_surf = occ.addPlaneSurface([hex_loop])

    # --- channel centers ---------------------------------------------------
    fuel_centers = []
    cool_centers = [(0.0, 0.0)]  # central coolant channel
    for i in range(n_pairs):
        a_fuel = math.radians(360.0 / n_pairs * i)
        a_cool = a_fuel + math.radians(360.0 / n_pairs / 2.0)  # interleaved
        fuel_centers.append((ring_radius * math.cos(a_fuel),
                             ring_radius * math.sin(a_fuel)))
        cool_centers.append((ring_radius * math.cos(a_cool),
                             ring_radius * math.sin(a_cool)))

    fuel_disks = [(2, occ.addDisk(cx, cy, 0.0, r_fuel, r_fuel))
                  for (cx, cy) in fuel_centers]
    cool_disks = [(2, occ.addDisk(cx, cy, 0.0, r_cool, r_cool))
                  for (cx, cy) in cool_centers]

    # Cut the coolant channels out (they become voids), then fragment the
    # remaining graphite with the fuel disks so the fuel compacts become their
    # own conforming subdomains.
    cut_res, _ = occ.cut([(2, hex_surf)], cool_disks,
                         removeObject=True, removeTool=True)
    frag_res, _ = occ.fragment(cut_res, fuel_disks)
    occ.synchronize()

    # --- classify surfaces into fuel vs graphite by centroid ---------------
    fuel_surf_tags, graph_surf_tags = [], []
    for dim, tag in gmsh.model.getEntities(2):
        com = occ.getCenterOfMass(dim, tag)
        matched_fuel = any(math.hypot(com[0] - cx, com[1] - cy) < 0.5 * r_fuel
                           for (cx, cy) in fuel_centers)
        (fuel_surf_tags if matched_fuel else graph_surf_tags).append(tag)

    # --- classify boundary curves into coolant vs outer --------------------
    cool_curve_tags, out_curve_tags = [], []
    boundary = gmsh.model.getBoundary([(2, t) for t in fuel_surf_tags + graph_surf_tags],
                                      combined=True, oriented=False, recursive=False)
    for dim, tag in boundary:
        com = occ.getCenterOfMass(1, abs(tag))
        matched_cool = any(math.hypot(com[0] - cx, com[1] - cy) < 0.5 * r_cool
                           for (cx, cy) in cool_centers)
        (cool_curve_tags if matched_cool else out_curve_tags).append(abs(tag))

    # --- physical groups ---------------------------------------------------
    pg_graph = gmsh.model.addPhysicalGroup(2, graph_surf_tags)
    pg_fuel = gmsh.model.addPhysicalGroup(2, fuel_surf_tags)
    pg_cool = gmsh.model.addPhysicalGroup(1, cool_curve_tags)
    pg_out = gmsh.model.addPhysicalGroup(1, out_curve_tags)

    # --- mesh --------------------------------------------------------------
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
    gmsh.option.setNumber("Mesh.MeshSizeMin", lc * 0.3)
    gmsh.model.mesh.generate(2)

    # --- extract nodes -----------------------------------------------------
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_coords = np.array(node_coords, dtype=float).reshape(-1, 3)
    tag2idx = {int(t): i for i, t in enumerate(node_tags)}
    nodes = node_coords[:, :2].copy()

    # --- extract triangles per material -----------------------------------
    tris_list, mat_list = [], []
    for pg, mat in ((pg_graph, MAT_GRAPHITE), (pg_fuel, MAT_FUEL)):
        for surf in gmsh.model.getEntitiesForPhysicalGroup(2, pg):
            etypes, _, enodes = gmsh.model.mesh.getElements(2, surf)
            for et, en in zip(etypes, enodes):
                if et != 2:  # 3-node triangle
                    continue
                conn = np.array(en, dtype=int).reshape(-1, 3)
                conn = np.vectorize(tag2idx.get)(conn)
                tris_list.append(conn)
                mat_list.append(np.full(conn.shape[0], mat, dtype=int))
    tris = np.vstack(tris_list)
    tri_mat = np.concatenate(mat_list)

    # --- extract boundary edges per physical curve ------------------------
    def collect_edges(pg):
        out = []
        for cv in gmsh.model.getEntitiesForPhysicalGroup(1, pg):
            etypes, _, enodes = gmsh.model.mesh.getElements(1, cv)
            for et, en in zip(etypes, enodes):
                if et != 1:  # 2-node line
                    continue
                conn = np.array(en, dtype=int).reshape(-1, 2)
                out.append(np.vectorize(tag2idx.get)(conn))
        return np.vstack(out) if out else np.zeros((0, 2), dtype=int)

    edges_cool = collect_edges(pg_cool)
    edges_out = collect_edges(pg_out)

    if show:
        gmsh.fltk.run()          # interactive viewer: zoom/pan/inspect elements

    gmsh.finalize()

    return BlockMesh(nodes=nodes, tris=tris, tri_mat=tri_mat,
                     edges_cool=edges_cool, edges_out=edges_out)


if __name__ == "__main__":
    m = build_block_mesh(verbose=False)
    print(f"nodes           : {m.n_nodes}")
    print(f"triangles       : {m.tris.shape[0]}")
    print(f"  fuel elements : {(m.tri_mat == MAT_FUEL).sum()}")
    print(f"  graph elements: {(m.tri_mat == MAT_GRAPHITE).sum()}")
    print(f"coolant edges   : {m.edges_cool.shape[0]}")
    print(f"outer edges     : {m.edges_out.shape[0]}")
