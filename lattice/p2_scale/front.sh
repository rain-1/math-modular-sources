#!/bin/bash
# lattice/p2_scale/front.sh -- with the interleaved split, report the largest n
# such that every n' <= n has all four k done (the contiguous front).
cd "$(dirname "$0")"
cat data/scan*_c*.csv | awk -F, '$2 ~ /^[0-9]+$/ {c[$2]++} END {n=3; while (c[n+1]==4) n++; print "contiguous front: n =", n}'
