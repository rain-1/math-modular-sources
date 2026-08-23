# `emn_projector` — the EMN adelic-projector search

Document: `consolidation/EMN_PROJECTOR.md`.  Names: the note's $S,H$ are `ang`, `Hfun`/`H`
in code (no PARI builtins are shadowed).  Everything is exact rational arithmetic; the two
slow sweeps write artifacts so they need not be rerun.

| script | what it does | artifact |
|---|---|---|
| `emn_core.py` | master ODE `z(2-z)F'' + (1-z)F' = Psi`, series solver, basis `Psi_{a,b}`, Chebyshev, denominator/valuation probes | — |
| `emn_span.py` | exact residue matrix of the basis, rational nullspaces, threshold sweep for `theta_min` | — |
| `emn_family.py A` | the optimal-geometry family (all poles with `theta<1/2` cancelled) | `family24.out` |
| `emn_cand.py` | series of an arbitrary combination; the `F3=-3H`, `F9=9H` and full-trace checks | — |
| `emn_table.py N` | candidate table: singular set, `t2`, `sigma_F`, `sigma_B`, slopes, entry | `table.json` |
| `emn_slopes.py N` | p-adic slope census for `a<=12`, both shift classes | `slopes.json`, `slopes.out` |
| `emn_meas.py N` | denominators/slopes of `H`, `K`, `Htil`, `B`, `A` | — |
| `emn_fold.py` | 30-digit verification of the fold at `z=2` (`B + 2G log(2-z)` converges) | — |
| `quick_type.py`, `quick_shape.py` | `n^2 C(2n,n) h_n` integral; minimal CDT shape of `B` | — |
| `emn_tau.py`, `emn_ledger.py`, `emn_opt.py` | tau(b;e) and the entry ledger (reuses `catalan_mu4/mu4_tau.py`, `cdt_finder/cdt_bound.py`) | — |
| `emn_bc.py` | Bost-Charles numerators for `phi = lambda/2` (reuses `catalan_mu4/mu4_bc.py`) | — |

Typical runtimes: `emn_span.py 24` ~ 15 min, `emn_slopes.py 1000` ~ 1 min, everything else
seconds to minutes.
