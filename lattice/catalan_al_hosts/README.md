# lattice/catalan_al_hosts

Scoring of the weight-3 chi_{-4} Eisenstein class (Catalan's constant) on every
genus-0 host Gamma_0(N)+W with 4 | N.  Document: `consolidation/CATALAN_AL_HOSTS.md`.

| script | what it does | log |
|---|---|---|
| `01_hosts.gp` | genus-0 Atkin-Lehner quotients with 4\|N, N<=120 (validated against the Fricke genus-0 list) | — |
| `04_geom.gp` | cusps of the quotient, order-2 elliptic points, #Sigma | — |
| `05_rational.gp` | rationality of the AL eigenspaces on M_1, M_3 (the "Q must be a perfect square" criterion) | — |
| `lib2.gp` | Ligozat orders, modular-unit exponent solver, hauptmodul from D log t | — |
| `06_rows.gp` | hauptmoduls per cusp placement, solve for F with F*Dt in M_3, A(t), integrality | `06run.log` |
| `07_sing.gp` | fit q2 A'' + q1 A' + q0 A = 0, read off the singular set, lambda_1, lambda_2 | `07run.log` |
| `08_full.gp` | growth \|a_n\|^(1/n), Eisenstein companions, b_n/a_n, k | `08run.log` |
| `09_foldreg.gp` | fold-regular subspace, xi, lindep against (1,zeta(2),G), decay | `09run.log`, `09run2.log` |
| `11_po4.gp` | LLL search for a class decaying below lambda_1 on the (3,1) placement of Gamma_0(12) | — |
| `10_score.py` | conformal ceilings, entry, margins (uses `lattice/cdt_finder/cdt_bound.py`) | — |

`lib.gp` is an earlier shared header, superseded by the self-contained scripts.
Note: `default(parisizemax,...)` must NOT appear inside a file that is `read()`,
it aborts the rest of the parse; set it in the top-level script instead.
