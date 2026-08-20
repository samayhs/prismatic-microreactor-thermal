#!/bin/bash
###############################################################################
# Verify the temperature-based solid solver removes the variable-Cp artifact.
#
# One steady run: production mesh + VARIABLE Cp(T) + chtMultiRegionTFoam.
# The enthalpy solver gives 684.6 C here (the +18 C artifact). If the T-based
# solver gives ~667-672 C with the SAME variable Cp(T), the artifact is
# eliminated by construction (steady conduction is Cp-independent).
#
# Runs in place in cfd_3d_unitcell/; restores production state (varCp thermo,
# h-keyed fvOptions) at the end. Run output is gitignored.
###############################################################################
source /opt/OpenFOAM-7/etc/bashrc
cd /mnt/c/Users/samayhs/Documents/PythonProjects/prismatic-microreactor-thermal/cfd_3d_unitcell
VDIR=verification

# mesh (production resolution)
rm -rf constant/polyMesh constant/fluid/polyMesh constant/solid/polyMesh
rm -rf 0 [1-9]* processor* postProcessing cellToRegion log.* *.foam
gmshToFoam $VDIR/work/prod.msh              > $VDIR/log.tcheck.gmshToFoam 2>&1
splitMeshRegions -cellZones -overwrite      > $VDIR/log.tcheck.split      2>&1
topoSet -region solid                       > $VDIR/log.tcheck.topoSet    2>&1

# variable Cp(T) thermo. The T solver reads the volumetric source from
# constant/solid/heatSource (tracked); fvOptions is left as-is and ignored by it.
cp $VDIR/thermo_varCp   constant/solid/thermophysicalProperties

rm -rf 0 && cp -r 0.orig 0
echo ">> running chtMultiRegionTFoam (variable Cp, T-based) ..."
chtMultiRegionTFoam > $VDIR/log.tcheck.solve 2>&1
echo "solver exit=$?"

PEAK=$(awk '/Min\/max T:/{s=$0; sub(/.*Min\/max T:/,"",s); split(s,a," "); if(a[1]+0>700) p=a[2]} END{if(p>0) printf "%.2f", p-273.15; else print "NA"}' $VDIR/log.tcheck.solve)
echo "=================================================================="
echo "  T-SOLVER peak fuel (variable Cp) = $PEAK C"
echo "  enthalpy solver on same setup    = 684.63 C  (artifact)"
echo "  correct Cp-independent steady    ~ 667-672 C"
echo "=================================================================="
echo "--- last solver lines ---"; tail -6 $VDIR/log.tcheck.solve

# restore production source (enthalpy-era files); adoption edits happen separately
cp $VDIR/thermo_varCp constant/solid/thermophysicalProperties
