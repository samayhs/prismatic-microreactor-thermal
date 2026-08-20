#!/bin/bash
###############################################################################
# Variable-Cp artifact — mesh-refinement verification (CHECK #1)
#
# For each of 3 hex meshes (coarse/prod/fine, uniform lc ratio ~1.26) run the
# steady solid CHT twice: constant Cp and variable Cp(T). The artifact at each
# resolution is
#       delta(mesh) = T_peak(varCp) - T_peak(constCp)
# whose artifact-free limit is exactly 0 (k(T) is identical in both variants).
# If delta shrinks ~O(h^2) toward 0, the +18 C is a benign spatial
# discretization error of OpenFOAM's enthalpy-based solid solver, not physics.
#
# Runs IN PLACE in the parent case dir (cfd_3d_unitcell/). All run artifacts are
# gitignored. The production thermophysicalProperties (variable Cp) is restored
# at the end. Meshes must already exist in verification/work/{coarse,prod,fine}.msh
# (regenerate on Windows:  python geometry/make_unit_cell.py --lc .. --nz .. --out ..).
#
# Usage (WSL, from cfd_3d_unitcell/):  bash verification/run_artifact_study.sh
###############################################################################
# NOTE: no `set -u` -- OpenFOAM's etc/bashrc references unbound vars and would abort.
cd "$(dirname "$0")/.." || exit 1          # -> cfd_3d_unitcell/
CASE="$(pwd)"
VDIR="$CASE/verification"
WORK="$VDIR/work"
RESULTS="${RESULTS_CSV:-$VDIR/results.csv}"

case "$CASE" in *" "*) echo "ERROR: space in path: $CASE"; exit 1;; esac
if [ -z "${WM_PROJECT_DIR:-}" ]; then
    source /opt/OpenFOAM-7/etc/bashrc 2>/dev/null || {
        echo "ERROR: OpenFOAM 7 not loadable"; exit 1; }
fi

# lc used to generate each mesh (representative cell size h for the O(h^2) fit)
declare -A LC=( [coarse]=0.00378 [prod]=0.00300 [fine]=0.00238 [extrafine]=0.00189
                [ucoarse]=0.00378 [uprod]=0.00300 [ufine]=0.00238 [uextrafine]=0.00189 )

# Meshes to run (default: the 3-mesh study). Pass names as args to run a subset,
# e.g. `run_artifact_study.sh extrafine` appends one mesh to an existing results.csv.
MESHLIST="${*:-coarse prod fine}"
if [ ! -f "$RESULTS" ] || [ -z "$*" ]; then
    echo "mesh,thermo,lc,solid_cells,peak_fuel_C" > "$RESULTS"
fi

peak_from_log() {   # last solid Min/max T (min > 700 K => fuel-bearing solid) -> peak C
    awk '/Min\/max T:/{s=$0; sub(/.*Min\/max T:/,"",s); split(s,a," ");
         if(a[1]+0>700) p=a[2]} END{if(p>0) printf "%.2f", p-273.15; else print "NA"}' "$1"
}

for mesh in $MESHLIST; do
    MSH="$WORK/$mesh.msh"
    [ -f "$MSH" ] || { echo "MISSING mesh $MSH"; exit 1; }
    echo "======================================================================"
    echo ">>> MESH: $mesh  (lc=${LC[$mesh]})"
    # --- mesh once per resolution -------------------------------------------
    rm -rf constant/polyMesh constant/fluid/polyMesh constant/solid/polyMesh
    rm -rf 0 [1-9]* processor* postProcessing cellToRegion log.* *.foam
    gmshToFoam "$MSH"                       > "$VDIR/log.$mesh.gmshToFoam"   2>&1
    splitMeshRegions -cellZones -overwrite  > "$VDIR/log.$mesh.split"        2>&1
    topoSet -region solid                   > "$VDIR/log.$mesh.topoSet"      2>&1
    # solid cell count from the polyMesh/owner header note ("... nCells:N ...")
    SCELLS=$(grep -a -o 'nCells:[0-9]*' constant/solid/polyMesh/owner 2>/dev/null \
             | head -1 | cut -d: -f2)
    [ -z "$SCELLS" ] && SCELLS=NA

    for thermo in constCp varCp; do
        echo "    -- thermo: $thermo"
        cp "$VDIR/thermo_$thermo" constant/solid/thermophysicalProperties
        rm -rf 0 [1-9]* processor* postProcessing cellToRegion
        cp -r 0.orig 0
        LOG="$VDIR/log.$mesh.$thermo"
        chtMultiRegionFoam > "$LOG" 2>&1
        PEAK=$(peak_from_log "$LOG")
        echo "$mesh,$thermo,${LC[$mesh]},${SCELLS:-NA},$PEAK" >> "$RESULTS"
        echo "       peak_fuel = $PEAK C   (solid cells ~ ${SCELLS:-?})"
    done
done

# restore production state: variable Cp
cp "$VDIR/thermo_varCp" constant/solid/thermophysicalProperties

echo "======================================================================"
echo "DONE. results -> $RESULTS"
column -t -s, "$RESULTS"
