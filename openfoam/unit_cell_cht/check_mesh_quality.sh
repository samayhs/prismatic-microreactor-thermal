#!/bin/bash
#
# Extract mesh-quality metrics for both regions and write a metrics report.
# Run from a set-up case directory (mesh converted + split into regions):
#     source /root/OpenFOAM-7/etc/bashrc
#     ./check_mesh_quality.sh [outputFile]
#
# Reports the checkMesh metrics that matter for a conjugate-heat-transfer solve:
# cell type, non-orthogonality (max/avg), skewness, aspect ratio, connectivity.

cd "${0%/*}" || exit 1
OUT="${1:-mesh_quality_metrics.txt}"

{
    echo "Mesh quality metrics  ($(date))"
    echo "Case: $(pwd)"
    echo "==========================================================="
    for r in fluid solid; do
        echo ""
        echo "----- region: $r -----"
        checkMesh -region "$r" -latestTime 2>/dev/null | grep -iE \
          "^\s+cells:|hexahedra|prisms|tet wedges|pyramids|polyhedra|Number of regions|non-orthogonality Max|Max skewness|Max aspect ratio|Mesh OK|failed"
    done
    echo ""
    echo "Guidance: non-orthogonality < 70 deg (max) is desirable; average"
    echo "< 30 deg is good. Skewness < 4. A single connected region per part"
    echo "confirms the heat-conduction path is intact."
} | tee "$OUT"

echo ""
echo "written: $OUT"
