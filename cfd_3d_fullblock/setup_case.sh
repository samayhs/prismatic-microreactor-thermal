#!/bin/bash
# Set up the full-block CHT case: copy the patch-name-driven scaffolding from the
# unit cell, convert the mesh, split regions, report region cell counts. Stops
# BEFORE the solver. Edit MSH below to point at the mesh to use.
source /opt/OpenFOAM-7/etc/bashrc
if [ -z "$WM_PROJECT_DIR" ]; then
    echo "ERROR: OpenFOAM 7 did not load (WM_PROJECT_DIR empty)."; exit 1
fi
echo "OpenFOAM $WM_PROJECT_VERSION loaded."

U=/mnt/c/Users/samayhs/Documents/PythonProjects/prismatic-microreactor-thermal/cfd_3d_unitcell
cd /mnt/c/Users/samayhs/Documents/PythonProjects/prismatic-microreactor-thermal/cfd_3d_fullblock || exit 1
MSH=geometry/full_block_coarse.msh

cp -r $U/0.orig  0.orig
cp -r $U/constant constant
cp -r $U/system   system
cp $U/geometry/make_toposet.py geometry/ 2>/dev/null
rm -rf constant/polyMesh constant/fluid/polyMesh constant/solid/polyMesh
rm -rf 0 [1-9]* processor* postProcessing cellToRegion

echo ">> gmshToFoam $MSH"
gmshToFoam $MSH > log.gmshToFoam 2>&1
if [ $? -ne 0 ]; then echo "gmshToFoam FAILED:"; tail -15 log.gmshToFoam; exit 1; fi

echo ">> splitMeshRegions -cellZones"
splitMeshRegions -cellZones -overwrite > log.split 2>&1
if [ $? -ne 0 ]; then echo "splitMeshRegions FAILED:"; tail -15 log.split; exit 1; fi

echo "=== regions created ==="
for d in constant/*/polyMesh; do
    r=$(echo "$d" | sed 's#/polyMesh##; s#constant/##')
    n=$(grep -a -o 'nCells:[0-9]*' "$d/owner" 2>/dev/null | head -1 | cut -d: -f2)
    echo "   region '$r' : $n cells"
done
echo ">> DONE (no solver run)"
