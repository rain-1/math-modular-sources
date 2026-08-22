# lattice/catalan_two_classes — two fold-regular chi_{-4} classes on one host

Write-up: `consolidation/CATALAN_TWO_CLASSES.md`.

| file | what it does | log |
|---|---|---|
| `00_setup.gp` | hosts (level-8 `t8,F8`; level-16 `x16,F16` as eta quotients), the inner/outer weight-3 chi_{-4} Eisenstein coefficient vectors, `mkPhi`/`Peval`, companion builder `Bof` | — |
| `01_check.gp` | sanity: reproduces `t(q)=q-4q^2+12q^3-...`, `A_n=1,4,20,112,676,...`, `b_n=0,1,7,404/9,...` | — |
| `02_classes.gp`,`02run.gp` | 11 classes x 2 hosts: P(0..3), predicted period, `b_n/a_n`, `|b_n|^{1/n}`, `|c_n|^{1/n}`, `v_2(b_n)/n` | `02run.log` |
| `03_level16.gp`,`03run.gp` | level-16 host in depth: basis companions for `E,V2E,V4E` and `T,V2T,V4T`, the fold-regular subspace, `Phi_0`, denominators | `03run.log` |
| `04_symmetrise.gp`,`04run.gp` | exact fold-regularity functional `8c2+c4=0`; normaliser descent `sigma(x)=-x/(4x+1)`, `y=4x^2/(4x+1)`; 2-adic slopes and radii of `Sym^+` in `y` | `04run.log` |
| `05_margin.py` | adelic margins on the level-16 host (reuses `lattice/adelic_holonomy/adelic_bound.py`) | stdout |
| `06_level32src.gp`,`06run.gp` | level-32 sources (`V8 E`) on the level-16 host: not fold-regular | `06run.log` |
| `07_ceiling.py` | thrice-punctured hyperbolic ceiling for `{-1/2,-1,infty}` at 0: `pi^2/Gamma(3/4)^4 = 4.376879` | stdout |
| `08_canon.gp`,`08run.gp` | canonical-source identities `F8*Dq(t8)=(1-8V2)E`, `F16*Dq(x16)=E`, `t=x(1+2x)/(1+4x)^2` | `08run.log` |

Run as `gp -q -s 4G 0Nrun.gp` from this directory (the `NTERM` in `0Nrun.gp` sets
the q-expansion length; 170 terms takes a few minutes).  Note: PARI's `read()`
needs multi-line statements wrapped in `{ }`, and `default(parisizemax,...)`
inside a `read` file aborts the rest of the file — set the stack with `gp -s`.
