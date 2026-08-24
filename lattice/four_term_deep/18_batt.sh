#!/bin/bash
# 18_batt.sh FILE  -- run the Catalan battery on a "label value" file
cd /home/ubuntu/code/math-modular-sources/lattice/four_term_deep
echo "XIFILE=\"$1\"; read(\"/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/07_catalan.gp\")" | gp -q -s 2000000000
