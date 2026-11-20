#!/bin/bash
# Build the reference drivers against geometry-central and write their data.
# Run from the repository root: pixi run -e ref ref-data
# Each tests/data/<area>/reference.cpp is compiled (with whatever
# geometry-central sources it needs, listed in SOURCES below) and run with
# its output going to tests/data/<area>/reference.txt. Each driver gets
# the geometry-central path as its argument (for the test meshes there).
set -e
GC=${GEOMETRY_CENTRAL_PATH:-../../c++/geometry-central}
if [ ! -d "$GC/include/geometrycentral" ]; then
    echo "geometry-central not found at $GC (set GEOMETRY_CENTRAL_PATH)" >&2
    exit 1
fi
echo "geometry-central: $GC at $(git -C "$GC" rev-parse --short HEAD)"
CXX=${CXX:-$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++}
# -ffp-contract=off: no fused multiply-adds, so the reference rounds every
# operation as the C++ is written
FLAGS="-O2 -std=c++17 -ffp-contract=off -I$GC/include -I$CONDA_PREFIX/include/eigen3"
mkdir -p build/ref

declare -A SOURCES=(
    [utilities]="$GC/src/utilities/elementary_geometry.cpp"
    [surface]="$GC/src/surface/surface_mesh.cpp $GC/src/surface/manifold_surface_mesh.cpp $GC/src/utilities/disjoint_sets.cpp"
)
for area in "${!SOURCES[@]}"; do
    $CXX $FLAGS "tests/data/$area/reference.cpp" ${SOURCES[$area]} -o "build/ref/$area"
    "build/ref/$area" "$GC" > "tests/data/$area/reference.txt"
    echo "wrote tests/data/$area/reference.txt"
done
