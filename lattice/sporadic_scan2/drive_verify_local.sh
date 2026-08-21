#!/bin/bash
cd /home/ubuntu/code/math-modular-sources/lattice/sporadic_scan2
H="$1"; OUT="$2"; NS="${3:-11}"
rm -f ${OUT}_*.jsonl
mkdir -p logs
for i in $(seq 0 $((NS-1))); do
  nohup python3 verify.py --hits $H --out ${OUT}_$i.jsonl --nv 210 --nshard $NS --shard $i \
      > logs/ver_${OUT}_$i.log 2>&1 &
done
