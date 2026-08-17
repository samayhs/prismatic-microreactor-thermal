#!/bin/bash
###############################################################################
# Run the 3D CHT unit-cell case IN PLACE, inside the repo.
#
# All generated output (mesh, regions, time steps, logs, postProcessing) stays in
# THIS folder and is gitignored — the repo is self-contained: tracked source +
# ignored run artifacts.
#
# Requirements:
#   * a space-free repo path (OpenFOAM's fileName class rejects spaces)
#   * OpenFOAM 7 at /opt/OpenFOAM-7
#   * geometry/unit_cell.msh present (regenerate on Windows:
#       python geometry/make_unit_cell.py)
#
# Usage (from a WSL shell):
#     cd .../openfoam/unit_cell_cht && ./run.sh
###############################################################################
cd "$(dirname "$0")" || exit 1
CASE="$(pwd)"

case "$CASE" in
  *" "*) echo "ERROR: OpenFOAM cannot run under a path containing spaces:";
         echo "  $CASE"; echo "Rename the folder to remove spaces."; exit 1;;
esac

source /opt/OpenFOAM-7/etc/bashrc

if [ ! -f geometry/unit_cell.msh ]; then
    echo "ERROR: geometry/unit_cell.msh not found."
    echo "Generate it first (on Windows):  python geometry/make_unit_cell.py"
    exit 1
fi

# Clean previous RUN ARTIFACTS only. NOTE: constant/fluid and constant/solid also
# hold SOURCE (thermophysicalProperties, fvOptions) — remove only their polyMesh,
# never the whole directory.
rm -rf constant/polyMesh constant/fluid/polyMesh constant/solid/polyMesh
rm -rf 0 [1-9]* processor* postProcessing cellToRegion log.* *.foam

echo ">> gmshToFoam";        gmshToFoam geometry/unit_cell.msh   > log.gmshToFoam 2>&1
echo ">> splitMeshRegions";  splitMeshRegions -cellZones -overwrite > log.splitMeshRegions 2>&1
echo ">> topoSet (fuel)";    topoSet -region solid              > log.topoSet 2>&1
# splitMeshRegions may have created an empty 0/; remove it before seeding fields,
# otherwise cp nests 0.orig as 0/0.orig and the solver can't find 0/fluid/p.
rm -rf 0
cp -r 0.orig 0

echo ">> solving (chtMultiRegionFoam) ... a few minutes"
chtMultiRegionFoam > log.chtMultiRegionFoam 2>&1
TIME=$(foamListTimes | tail -1)

# fields for the Python visualizer, and ParaView entry points
postProcess -region fluid -func writeCellCentres -latestTime > /dev/null 2>&1
postProcess -region solid -func writeCellCentres -latestTime > /dev/null 2>&1
postProcess -region solid -func writeCellVolumes -latestTime > /dev/null 2>&1
touch fluid.foam solid.foam

PEAK=$(awk '/Min.max T/{s=$0; sub(/.*Min.max T:/,"",s); split(s,a," "); if(a[1]+0>574) m=a[2]} END{printf "%.1f", m-273.15}' log.chtMultiRegionFoam)
echo ""
echo "=================================================================="
echo "  DONE.  peak fuel temperature = ${PEAK} C   (latest time ${TIME})"
echo "  Output is in this folder (gitignored)."
echo ""
echo "  Visualize (Windows):   python visualize_cfd.py . ${TIME}"
echo "                         or open fluid.foam / solid.foam in ParaView"
echo "  Energy balance:        see validation/ and validation/TROUBLESHOOTING.md"
echo "=================================================================="
