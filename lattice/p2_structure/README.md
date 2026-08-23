# lattice/p2_structure — structure of the Catalan congruence lattices

Scripts and data for `consolidation/P2_STRUCTURE.md`.

Run pattern (PARI/GP; concatenate, do not `read()`, because the builders are brace blocks):

    cat lattice/positivity/rows_pos.gp lattice/p2_structure/p2core.gp \
        lattice/p2_structure/p2run.gp  lattice/p2_structure/run_main.gp > run.gp && gp -q run.gp

`rows_pos.gp` (from `lattice/positivity/`) supplies the exact rows `zudrow`, `nestrow`, the
moment family `mom(m,j)` and `klat2`; `p2core.gp` raises `parisizemax` to 8 GB.  Drivers set
`\p 3000` — the weighted Gauss reduction needs `G` at >= 3000 digits.

| file | prepend | entry points |
|---|---|---|
| `build_rows.gp` | `rows_pos.gp` | `buildrows(NLO,NHI,OUTF)` — caches exact `n X Y V U` |
| `p2core.gp` | `rows_pos.gp` | `rdrows`, `kfull`, `kcong`, `gred`, `coneiv`, `conescan`, `eight` |
| `p2run.gp` | + `p2core.gp` | `anal`, `analS`, `cfclass`, `cfhist` |
| `kernel.gp` | + `p2core.gp` | `kerchk`, `bcrit` |
| `verify.gp` | + `p2core.gp` | `chk(NHI)`, `qpchk(NHI)` — the exact identities of §2.1, §5.1 |
| `run_main.gp` | + `p2run.gp` | the 27-column structural table -> `data/struct_n120.csv` |
| `run_hist.gp` | + `p2run.gp` | partial-quotient histograms -> `data/cfhist_n120.csv` |
| `run_cong.gp` | + `p2run.gp` | single-congruence lattice -> `data/structS_n120.csv` |
| `run_kernel.gp` | + `kernel.gp` | kernel identity and `b_crit` -> `data/kernel_n120.csv` |
| `run_surrogate.gp` | + `p2run.gp`, and a line `EEXP=320;` | the whole table for `bestappr(G,10^EEXP)` |
| `stats.py` | — | verification tallies, P2' regressions, Gauss-Kuzmin, Haar-random Monte Carlo |

`anal` **appends** to its dump file, so delete `data/vectors_n120.txt` before re-running
`run_main.gp` (the CSV on stdout is regenerated from scratch either way).

Rebuilding the row cache is the only expensive step (~1.8 h CPU to n=120, 6.8 GB peak at
n=120; split the range over several jobs).  `data/rows_all.txt` holds it for 4 <= n <= 120;
everything else runs in under a minute.

Conventions: oriented coordinates `(a,b) = (sign(XG-Y) c_Z, sign(VG-U) c_N)`, so the positive
cone is the closed first quadrant and the metric is `diag(lz,ln)`.  No builtin names are
shadowed (no `psi`, `M`, `Phi`, `S`, `cmp`; matrices `MAT`/`BB`/`HH`, moduli `MOD`).
