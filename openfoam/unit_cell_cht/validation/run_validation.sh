#!/bin/bash
#
# Generate the fields the energy-balance check needs, then run it.
# Run from a converged case directory:
#     source /root/OpenFOAM-7/etc/bashrc
#     .../validation/run_validation.sh <caseDir> <time>
#
CASE="${1:-.}"
TIME="${2:-4000}"
cd "$CASE" || exit 1

# cell volumes (for V_fuel) and interface heat flux (Q_wall)
postProcess -region solid -func writeCellVolumes -latestTime > /dev/null 2>&1
echo "=== interface heat flux (Q_wall) ==="
chtMultiRegionFoam -postProcess -region fluid -func wallHeatFlux -latestTime 2>&1 \
    | grep -iE "integ\(fluid_to_solid\)"

echo "=== energy balance ==="
python3 "$(dirname "$0")/energy_balance.py" "$CASE" "$TIME"
