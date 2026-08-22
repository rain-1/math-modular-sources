# lattice/holonomy_lindep

Scripts for `consolidation/HOLONOMY_LINDEP.md`.

**Task A — linear independence via the CDT holonomy bound for multi-period rows.**

| file | what it does |
|---|---|
| `bound_ml.py` | the whole Task-A table: `tau`, entry, margin, deficit, headroom, ranking. Needs `../cdt_finder/conformal.py` and `../cdt_noncongruence/delta_table.json`. `python3 bound_ml.py` |
| `delta_ext.py`, `delta_ext.json` | extends `Delta(r) = BC(phi_r) - log|phi_r'(0)|` for the Kodaira family to `r in [0.80,0.85]` (the quadrature degenerates from `r=0.86`) |
| `row_aesz207.gp` | AESZ 207: char. poly, `k=4` sharp for `B,C,D`, Apery limits to 1354–1361 digits (`N=6000`), period rank over Q and Q(sqrt17), fold-regularity (`|h_n| = c*53248^n/n^4`) |
| `row_beta4_24.gp` | `beta(4)@24 = Sym^3(Zagier E)`: `k=4` sharp for all five companions, limits by asymptotic fit (55 digits; Cauchy gives 2–4), period rank 3 with two explicit Q-relations, and the **failure** of fold-regularity (triple dominant root, codimension-3 log tower) |
| `x15a.gp`, `x15b.gp` | `X_1(5) Sym^2`: minimal recurrence (order 4, degree 3), `chi = (x^2-11x-1)^2`, `k=3` sharp, and the **failure** of fold-regularity (double dominant root) — the correction to `CDT_FINDER.md` §4/§6 |
| `survey_multiperiod.gp`, `s12b.gp`, `s4_denexp.gp`, `badprime_all.gp`, `relverify.gp`, `step12.csv`, `ranked.json` | the 297-operator MUM sweep: z-degree, root multiplicity, LCM-normalisation, period rank, score |
| `topgeom.gp` | singular sets and first `A_n` of the top MUM rows |

**Task B — irrationality measures.**

| file | what it does |
|---|---|
| `measure.py` | `kappa_*` from ICM eq. `withintegrations`; calibrated by reproducing CDT's printed `24781` for `L(2,chi_-3)`; then Beukers' row, Apery's `zeta(2)` row, and `zeta(3)` as control. `python3 measure.py` |

Conventions: absolute paths; `default(parisizemax, ...)` in every `.gp`;
PARI builtin names avoided. Logs sit beside the scripts.
