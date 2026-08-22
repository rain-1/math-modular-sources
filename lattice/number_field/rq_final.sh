#!/bin/bash
cd /home/ubuntu/code/math-modular-sources/lattice/number_field
for S in 3 2; do
  cat rq_out_R${S}_*.txt 2>/dev/null | grep "^R$S " | sort -u > rq_all_R$S.txt
  ./rq_run_verify.sh rq_all_R$S.txt rq_ver_all_R$S.txt
  echo "R$S raw rows: $(wc -l < rq_all_R$S.txt)"
  tail -1 rq_ver_all_R$S.txt
  echo "  false positives: $(grep -c FALSEPOS rq_ver_all_R$S.txt)"
  echo "  non-degenerate: $(grep '^ROW|' rq_ver_all_R$S.txt | awk -F'|' '$12==0'|wc -l)"
  echo "  Mode-I  > 0   : $(grep '^ROW|' rq_ver_all_R$S.txt | awk -F'|' '$12==0 && $19>0'|wc -l)"
  echo "  Mode-II > 0   : $(grep '^ROW|' rq_ver_all_R$S.txt | awk -F'|' '$12==0 && $21>0'|wc -l)"
done
