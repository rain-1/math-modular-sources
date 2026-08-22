# `lattice/cdt_finder` — scoring modular Apéry hosts against the CDT holonomy bound

Report: `consolidation/CDT_FINDER.md`.

| file | what it does |
|---|---|
| `cdt_bound.py` | exact implementation of CDT (arXiv:2408.15403) Def. 6.0.1 ($I_u^v(w)$), $\tau^\flat$ (6.0.4), $\tau^\sharp$ (6.0.5), and the bound (7.0.1). Running it reproduces CDT's $191/49$, $I_2^{14}(2)=21.075$, $27/80$, $16603/3920$, $|\varphi'(0)|=161.081$, $m\le13.9938<14$. |
| `conformal.py` | Poincaré (conformal) radius of punctured planes; exact for $\mathbb C\setminus\{a,b\}$ via the modular $\lambda$; unit test reproduces the archive's level-10 $\tau_0$ to 22 digits and $r_u=1.30651071$. |
| `indep_check2.py` | rebuilds CDT's 14 functions exactly from their ODE (Prop. 11.1.4) and their explicit pure series, and checks $\mathbb Q(y)$-linear independence by rank mod $2^{61}-1$ (no relation with $\deg P_i\le5$). |
| `x15_sym2.py` | the $X_1(5)$ $\mathrm{Sym}^2$ (Beukers 1987 Thm 4) host: $A_n=[t^n]F^2$, minimal recurrence, characteristic roots, absence of a free integration, sharp $k$. |
| `x15_margin.py` | conformal sizes at both archimedean places and the entry/margin table for that host, in both number-field normalisations. |
| `hosts.py`, `final_table.py` | the ranked census table over the paper census + `lattice/sporadic_scan2/table.json`. |

All scripts are pure Python (`fractions`, `sympy`, `mpmath`) and run in seconds to a
few minutes; none needs PARI or `snake`.
