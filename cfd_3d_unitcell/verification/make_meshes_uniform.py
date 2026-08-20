"""
Generate UNIFORMLY-refined unit-cell meshes for the corrected artifact study.

The first study refined only the bulk size `lc` while pinning the near-wall size
`BL_FIRST` (fixed for y+). That leaves the steep-gradient near-coolant-wall solid
UNREFINED -- exactly where the variable-Cp truncation is generated -- so the
artifact delta could not shrink. Here we scale BL_FIRST and BL_THICK with lc
(same ratios as production: BL_FIRST/lc = 0.4333, BL_THICK/lc = 2.0) so the WHOLE
domain, near-wall included, refines uniformly. Fluid y+ then changes with mesh,
but that cancels in delta = varCp - constCp (same fluid at each mesh).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.join(HERE, "..", "geometry")
sys.path.insert(0, GEO)
import make_unit_cell as mkuc

BL_FIRST_RATIO = 0.0013 / 0.0030   # = 0.4333, production ratio
BL_THICK_RATIO = 0.0060 / 0.0030   # = 2.0

# (name, lc, nz) -- same lc/nz as the first study for direct comparison
MESHES = [("ucoarse",    0.00378, 32),
          ("uprod",      0.00300, 40),
          ("ufine",      0.00238, 50),
          ("uextrafine", 0.00189, 63)]

outdir = os.path.join(HERE, "work")
os.makedirs(outdir, exist_ok=True)
for name, lc, nz in MESHES:
    mkuc.BL_FIRST = BL_FIRST_RATIO * lc
    mkuc.BL_THICK = BL_THICK_RATIO * lc
    out = os.path.join(outdir, name + ".msh")
    print(f">>> {name}: lc={lc} nz={nz} BL_FIRST={mkuc.BL_FIRST:.5f} "
          f"BL_THICK={mkuc.BL_THICK:.5f}")
    mkuc.make_unit_cell(out, lc=lc, n_axial=nz, boundary_layer=True)
