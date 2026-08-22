#!/bin/sh
# run_all.sh -- regenerate every log in out/ for consolidation/CUSP_MOVE_PROGRAM.md
set -e
cd "$(dirname "$0")"
mkdir -p out
timeout  600 gp -q 00_selftest.gp          </dev/null > out/00_selftest.log        2>&1
timeout 1800 python3 01_general_move.py                > out/01_general_move.log     2>&1
timeout  900 gp -q 02_orbits.gp            </dev/null > out/02_orbits.log          2>&1
timeout 3000 gp -q 03_census.gp            </dev/null > out/03_census.log          2>&1
timeout 7200 gp -q 04_periods.gp           </dev/null > out/04_periods.log         2>&1
timeout  300 gp -q 05_placements.gp        </dev/null > out/05_placements.log      2>&1
timeout 7200 gp -q 06_padic.gp             </dev/null > out/06_padic.log           2>&1
timeout  300 gp -q 07_measures.gp          </dev/null > out/07_measures.log        2>&1
timeout  900 gp -q 08_quadratic.gp         </dev/null > out/08_quadratic.log       2>&1
timeout  900 gp -q 09_conj_xi.gp           </dev/null > out/09_conj_xi.log         2>&1
timeout  600 gp -q 10_domb_ident.gp        </dev/null > out/10_domb_ident.log      2>&1
timeout  600 gp -q 11_companion_general.gp </dev/null > out/11_companion_general.log 2>&1
timeout  900 gp -q 12_domb_pairs.gp        </dev/null > out/12_domb_pairs.log       2>&1
echo "done; logs in out/"
