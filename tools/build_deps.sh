#!/bin/bash
# Precompile MojoDEC, its linear algebra package molina, and its vendored
# NuMojo into build/, with this workspace's Mojo (a .mojoc only loads in
# the Mojo version that built it).
# Run from the repository root: pixi run deps
# MOJODEC_PATH points at the MojoDEC checkout (default: ../MojoDEC).
set -e
MOJODEC_PATH=${MOJODEC_PATH:-../MojoDEC}
if [ ! -d "$MOJODEC_PATH/mojodec" ]; then
    echo "error: MojoDEC not found at $MOJODEC_PATH (set MOJODEC_PATH)" >&2
    exit 1
fi
mkdir -p build
mojo precompile "$MOJODEC_PATH/third_party/numojo" -o build/numojo.mojoc
mojo precompile "$MOJODEC_PATH/molina" -o build/molina.mojoc
mojo precompile -I build "$MOJODEC_PATH/mojodec" -o build/mojodec.mojoc
