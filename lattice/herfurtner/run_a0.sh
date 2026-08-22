#!/bin/bash
cd "$(dirname "$0")"
> out/a0.txt
while read M J1 J2; do
  [ -z "$M" ] && continue
  CMAX=$(( 12000 / (M*M) )); [ $CMAX -gt 4000 ] && CMAX=4000
  ./02_hscan $M $J1 $J2 0 0 400 $CMAX 30 >> out/a0.txt 2>/dev/null
done < classes.txt
echo "A=0 hits: $(wc -l < out/a0.txt)"
