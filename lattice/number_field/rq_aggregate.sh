#!/bin/bash
# usage: rq_aggregate.sh <shape:R3|R2> <verfile...>
cd /home/ubuntu/code/math-modular-sources/lattice/number_field
S=$1; shift
cat "$@" | grep "^ROW|" | awk -F'|' -v s="$S" '
  { if((s=="R3"&&$2==3)||(s=="R2"&&$2==2)) print }' | sort -u > /tmp/agg_$S.raw
# dedupe by (D,a0,a1,b0,b1,c0,c1)
awk -F'|' '!seen[$3"|"$4"|"$5"|"$6"|"$7"|"$8"|"$9]++' /tmp/agg_$S.raw > /tmp/agg_$S.uniq
echo "shape $S: total distinct integral classes across passes: $(wc -l < /tmp/agg_$S.uniq)"
echo "  non-degenerate (A^2!=C / A^2!=4C, C!=0, not A=B=0): $(awk -F'|' '$12==0 && !($4==0&&$5==0&&$6==0&&$7==0)' /tmp/agg_$S.uniq | wc -l)"
echo "  with Mode-I margin > 0 : $(awk -F'|' '$12==0 && !($4==0&&$5==0&&$6==0&&$7==0) && $19>0' /tmp/agg_$S.uniq | wc -l)"
echo "  with Mode-II margin > 0: $(awk -F'|' '$12==0 && !($4==0&&$5==0&&$6==0&&$7==0) && $21>0' /tmp/agg_$S.uniq | wc -l)"
