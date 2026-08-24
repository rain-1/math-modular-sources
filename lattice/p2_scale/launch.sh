#!/bin/bash
# lattice/p2_scale/launch.sh -- run the exact scan on the 12 cost-balanced
# chunks in parallel.  Each worker builds its own concatenated gp source
# (never read(), per the project convention) and appends CSV lines as it goes,
# so partial results survive an interruption.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
LAT="$HERE/.."
TMP="${TMPDIR:-/tmp}/p2scale"
mkdir -p "$TMP" "$HERE/log"
KL="${KLIST:-[22.35, 22.4, 23.0, 23.9]}"
TAG="${TAG:-scan}"
while read -r J LO HI; do
  P="$TMP/params_$J.gp"
  cat > "$P" <<EOF
NL = $LO; NH = $HI; NSTEP = ${NSTEP:-1};
KLIST = $KL;
OUTF = "$HERE/data/${TAG}_c$J.csv";
ROWF = "$HERE/data/rows_c$J.txt";
EOF
  rm -f "$HERE/data/${TAG}_c$J.csv"
  cat "$LAT/positivity/rows_pos.gp" "$LAT/p2_structure/p2core.gp" \
      "$HERE/scorel.gp" "$HERE/srun.gp" "$P" "$HERE/run_scan.gp" > "$TMP/run_$J.gp"
  nohup gp -q "$TMP/run_$J.gp" > "$HERE/log/${TAG}_worker_$J.log" 2>&1 &
  echo "worker $J : n = $LO .. $HI   pid $!"
done < "$HERE/data/chunks.txt"
wait
echo "all workers finished"
