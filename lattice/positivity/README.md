# lattice/positivity — the positivity programme

Scripts and data for `consolidation/POSITIVITY_PROGRAM.md`.

Run pattern (PARI/GP; concatenate, do not `read()`, because the row builders are brace blocks):

    cat lattice/positivity/rows_pos.gp [deps...] X.gp > run.gp && gp -q run.gp

`rows_pos.gp` sets `default(parisizemax, 12000000000)`; put a `default(parisize, ...)` in the
driver to size the initial stack (600 MB is enough up to m = 60 / n = 60; n = 80 peaks at
2.3 GB).

| file | prepend | entry points |
|---|---|---|
| `rows_pos.gp` | — | `zud(M)`, `mom(m,j)`, `zudrow(n)`, `nestrow(n)`, `klat2` |
| `signs.gp` | — | `census(name,ord,a,b,c,xi,N)`, `tailcheck(...)` |
| `hankel.gp` | `signs.gp` | `hank(name,ord,a,b,c,xi,dmax)` |
| `intreps.gp` | `rows_pos.gp` | `momnum(m,j)`, `beuk2(n)`, `beuk3(n)`, `chi3ker()` |
| `cone80.gp` | `rows_pos.gp` | `sweepn(n, klist, G)`, `conemin`, `redu` |
| `grid.gp` | `rows_pos.gp` | `gridrow(m)` |
| `pairs.gp` + `pairs_run.gp` | `rows_pos.gp`, `cone80.gp` | `pairstat(...)`, `runm(m)` |
| `triples.gp` | + `pairs.gp` | `runtrip(m, R)` |
| `control.gp` | + `pairs.gp` | `ctrl_pairs(m, Gstar)` |

Python post-processing: `fit.py`, `fit2.py`, `pnt.py` (cone data), `gan.py` (grid),
`pan.py`, `pan2.py` (pairs).  They read the CSVs from the session scratchpad by default;
point them at `data/` to reproduce.

`data/` holds the three tables the document quotes: `cone_n80.csv` (231 rows),
`grid_m60.csv` (1890), `pairs_m44.csv` (13357).  Column meanings are in
`POSITIVITY_PROGRAM.md` §7.

**Precision.**  `sweepn` needs `G` at >= 3000 digits (production used `\p 4000`); `pairs.gp`
and `triples.gp` run at `\p 700`.  `redu`/`redu2` size their LLL rounding to the precision
actually present in their arguments, so the files can be concatenated in any order without a
`roundr` crash — but a low-precision `G` degrades the answers silently, so set the precision
in the driver.
