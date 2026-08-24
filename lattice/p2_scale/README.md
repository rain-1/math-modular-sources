# `lattice/p2_scale/` — the P2$'$ programme at scale ($4\le n\le10\,000$)

Companion code for `consolidation/P2_SCALE.md`.  Extends
`lattice/p2_structure/` ($n\le120$) and `lattice/p2_holonomic/` ($n\le200$).
Everything decisive is exact integer arithmetic; Catalan's constant enters
only through a **dyadic bracket** on the metric ratio, and every instance is
computed twice (once at each end of the bracket) and certified by agreement.

Run pattern (concatenate, never `read()`):

    cat lattice/positivity/rows_pos.gp lattice/p2_structure/p2core.gp \
        lattice/p2_scale/scorel.gp lattice/p2_scale/srun.gp \
        PARAMS lattice/p2_scale/run_scan.gp > run.gp && gp -q run.gp

| file | prepend | what |
|---|---|---|
| `build_rows.gp` | `rows_pos.gp`, `p2core.gp`, `p2_holonomic/rowrec.gp` | exact rows $X_n,Y_n,V_n,U_n$ to $n=N$ from the order-3 P-recursion; re-verifies the operator on the $n\le120$ cache and the rows against the cached $n\le200$ |
| `scorel.gp` | `rows_pos.gp`, `p2core.gp` | `hermd` (closed-form Hermite data + three exactness checks), `gredx` (exact integer Gauss reduction), `conex`/`eightx` (exact cone scan), `balx` (two-pass streamed balance index) |
| `srun.gp` | + `scorel.gp` | `anal2`: one certified CSV line per $(k,n)$ |
| `run_scan.gp` | + `srun.gp` | driver; needs `NL, NH, KLIST, OUTF, ROWF` |
| `verify.gp` | `rows_pos.gp`, `p2core.gp` | the three exact 2-adic identities that validate the extended rows independently of the fitted recurrence |
| `split_rows.py` | — | per-worker row split: contiguous equal-cost blocks (fastest to a complete range) or interleaved, `-i` (the completed range is contiguous at every moment) |
| `launch.sh` | — | runs the 12 chunks in parallel; `TAG=` names the output, `NSTEP=` for an interleaved split |
| `front.sh` | — | the contiguous front of an interleaved run |
| `checkwin.gp` | + `scorel.gp` | brute-force control: re-evaluates the balance objective at every rung, with no window, and compares |
| `stats.py`, `stats2.py`, `stats3.py` | — | all statistics of `P2_SCALE.md` |

| data file | content |
|---|---|
| `data/rows_scale.txt` | exact $n,X_n,Y_n,V_n,U_n$ (2.5 GB at $n=10\,000$; gitignored, regenerable in 474 s) |
| `data/rows_c*.txt`, `data/chunks.txt` | the per-worker split |
| `data/scan_c*.csv`, `scan2_c*.csv`, `scan3_c*.csv` | the scan, 39 988 lines, 41 columns, header in each file |
| `data/verify*.csv` | the 2-adic row validation |
| `data/summary*.json` | everything the three `stats` scripts compute |

**Why it is fast.**  Three replacements for the $n\le200$ pipeline:
`matkerint` (13.5 s at $n=1000$) by the closed forms of `P2_HOLONOMIC.md` §1.2
plus one modular inverse ($<10^{-3}$ s); the convergent ladder by its own
three-term recurrence $e_i=a_ie_{i-1}+e_{i-2}$ streamed twice, with exact
squares only in a proof-carrying window of about ten indices around the
balance; and the weighted Gauss reduction by a reduction started from the
ladder pair $(w_{i},w_{i-1})$, which Theorem 1 of `P2_STRUCTURE.md` predicts is
already reduced — the loop is still run, and its verdict recorded in the
`redst` column, so the theorem is re-verified rather than assumed.
