#!/bin/sh
# Reproduction driver for consolidation/ACF_ONE_SURFACE.md
cd "$(dirname "$0")"
mkdir -p out
for f in 00_mirror 01_cusp_move 02_sources 03_companion 06_units_and_level12 07_orbits 08_census_check 05_padic; do
  echo "=== $f ==="
  timeout 560 gp -q $f.gp </dev/null > out/$f.log 2>&1
  echo "  -> out/$f.log"
done
python3 04_cuspmove_proof.py > out/04_cuspmove_proof.log 2>&1
echo "  -> out/04_cuspmove_proof.log"
