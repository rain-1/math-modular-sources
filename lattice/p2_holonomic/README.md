# lattice/p2_holonomic — is the Catalan congruence lattice holonomic?

Scripts and data for `consolidation/P2_HOLONOMIC.md`.  Sequel to `lattice/p2_structure/`.

Run pattern (PARI/GP; concatenate, do not `read()` — the builders are brace blocks):

    cat lattice/positivity/rows_pos.gp lattice/p2_structure/p2core.gp \
        lattice/p2_holonomic/hcore.gp lattice/p2_holonomic/run_h2b.gp > run.gp && gp -q run.gp

`rows_pos.gp` supplies the exact rows and `klat2`; `p2core.gp` supplies `rdrows`, `kfull`,
`gred`, `conescan`, `eight`; `hcore.gp` adds the Hermite data and the exact balance-index
interval; `adelic.gp` adds the 2-adic period.  Drivers that need `G` set `\p 3000`.

| file | prepend | entry points |
|---|---|---|
| `hcore.gp` | rows_pos, p2core | `hdat`, `hpred`, `hchk12`, `ladder`, `balance` |
| `rowrec.gp` | — | `fitrec2`, `extend` (order-3 recurrence for the Nesterenko row) |
| `adelic.gp` | + hcore | `xi2at` (xi2 = zeta_2(2) from the Zudilin row), `vform` |
| `build_rows200.gp` | + rowrec | rebuilds `data/rows_n200.txt` (0.2 s; reproduces `p2_structure/data/rows_all.txt` exactly on n <= 120) |
| `verify.gp` | + hcore | the three exact identities of Proposition 1 (0/591 failures) |
| `run_h1.gp` / `run_h1b.gp` | + hcore | Hermite closed forms, n <= 120 / n <= 200 |
| `run_h2.gp` / `run_h2b.gp` | + hcore | balance index, exact stability interval, d*(n) |
| `run_h3.gp` | + hcore | truncation experiment: which digits of G matter |
| `run_h4.gp` | + adelic | adelic qualities of the cone minimisers |
| `run_h5.gp` | + adelic | target (non-)convergence; `algdep` on xi2 |
| `run_h7.gp`, `run_h8.gp` | + hcore | slide the 2-adic exponent E by one bit |
| `parity.py` | — | all parity statistics (needs numpy, scipy) |
| `control_slide.py` | — | fair-coin control for the E-slide |

**Caveat on `data/h4_adelic.csv`.**  It was produced with xi2 at 4935 2-adic bits (m = 620), so
its `v2NestForm` and `v2Nesttail` columns saturate for n >~ 175.  The reliable valuations are in
`data/h6_val.csv` (xi2 at 7011 bits, m = 880); all statements in the write-up use those.  The
`v2qxip`, `logq`, `archrate`, `d2`, `dinf`, `dad` columns of `h4_adelic.csv` are unaffected
(they are governed by the smaller Zudilin valuation).

Conventions: oriented coordinates `(u,v) = (sign(XG-Y) c_Z, sign(VG-U) c_N)`; no builtin names
shadowed (no `psi`, `M`, `Phi`, `S`, `cmp`; moduli are `MOD`, matrices `BB`/`HH`/`B0`).
