#!/bin/bash
# 02_periods.sh -- re-verify the full period matrices (fold constants at all three
# finite singular points) of the six mixed-exponent hosts with the independent
# ODE-continuation code of lattice/four_term_deep/11_foldmix.py (the same Frobenius
# + Taylor-continuation machinery as lattice/k3_period).  Two precisions are run and
# the agreeing digit count is reported.
D=/home/ubuntu/code/math-modular-sources/lattice/four_term_deep
O=/home/ubuntu/code/math-modular-sources/lattice/four_term_cdt/out
CLS="-1 2 0 1 1 0 0"
for spec in "R1 8 16 8 48 0 -128" "R2 2 14 8 28 4 8" "R3 -2 6 4 -32 -8 32" \
            "R4 8 16 8 68 8 32" "R5 1 17 10 32 8 16" "R6 -1 13 8 -13 -1 -1" \
            "R7 4 8 4 32 8 64"; do
  set -- $spec
  lab=$1; shift
  echo "=== $lab  r=$1 (a,c,d,f,C)=($2,$3,$4,$5,$6) ==="
  python3 $D/11_foldmix.py row $CLS "$@"
done
