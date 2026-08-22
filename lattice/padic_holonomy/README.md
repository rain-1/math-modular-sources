# `lattice/padic_holonomy` — CDT's arithmetic holonomy bound applied to $p$-adic $L$-values

Report: `consolidation/PADIC_HOLONOMY_CENSUS.md`.

Runs the recipe of Calegari–Dimitrov–Tang, ICM survey §6.2 (their proof that
$\zeta_2(5)\in\mathbb Q_2$ is irrational) on every $p$-adic near miss of
`consolidation/PADIC_IRRATIONALITY_CENSUS.md`.

## Exact arithmetic (PARI/GP)

| file | what it does |
|---|---|
| `cdt_ab.gp` | library: hauptmoduln of $X_0(p),X_0(p^2),X_1(4)$, Eisenstein series, reversion of $x(q)$, the coefficients $a_n,b_n$ of $H=E^*_{2k}(E'_{-2k}+\eta)=\sum(a_n+\eta b_n)x^n$ |
| `cdt_audit.gp`, `cdt_task1.gp`, `cdt_task23.gp` | integrality of $b_n$ and sharpness of $[1..n]^{2k+1}a_n\in\mathbb Z$, for $X_0(p)$ and the two character rows |
| `cdt_task4.gp`, `cdt_task4b.gp` | measured $p$-adic slope ($=\log_pR_p$) and archimedean growth |
| `cdt_branch.gp`, `cdt_x09cusp.gp` | exact location of the nearest archimedean singularity |
| `cdt_task5.gp`, `cdt_task5b.gp` | rank / $\mathbb Q(x)$-independence of $\{1,H,\dots,H^{(2k)}\}$ mod $2^{61}-1$ |
| `cdt_zetacheck.gp`, `cdt_Lcheck.gp` | identification of $\eta$ against Kubota–Leopoldt / Bernoulli |
| `cdt_sharp.gp` | which primes obstruct the smaller LCM exponent |
| `data_X0_p_k*.txt` | the exact $a_n,b_n$ for $n\le60$ |

## Complex analysis (python + numpy/scipy/mpmath)

| file | what it does |
|---|---|
| `haupt.py` | evaluate the hauptmodul anywhere in $|q|<1$: $\Gamma_0(N)$-reduction by Gauss lattice reduction; branch-free $\log|x|$ via $\log|\Delta|$ and $\mathrm{SL}_2(\mathbb Z)$ reduction |
| `family.py` | composable structured templates SCALE / OFF / TANG / BITE (CDT's lunes, L2chi App. A.1.1) |
| `outer.py` | **the exact search space**: every admissible $\psi$ is $z\cdot(\text{outer function of }u)$ for a boundary log-modulus $u\le0$, with $\log|\psi'(0)|=\bar u$.  Rearrangement and Bost–Charles quadratures |
| `bcint.py` | earlier structured-family versions of the same integrals |
| `targets.py` | the eleven cells: $d$, $m$, $\log R_p$, $\tau(\mathbf b)$, budget |
| `optimise.py`, `driver.py`, `scan1.py` | structured-family scans → `res_*.json`, `drv_*.log` |
| `freeopt.py` | free-template optimisation → `free_*.json`, `free_*.log` |
| `certify.py` | conversion to the explicit $\psi(z)=z\,e^{P(z)}$ form + admissibility check → `cert_*.json` |
| `precise5.py` | high-precision two-parameter optimum for $\zeta_5(3)$ → `precise5.log` |
| `calib.py` | reproduces CDT's two $\zeta_2(5)$ contours → `calib.log` |
| `assemble.py` | final table, $\kappa$, contour profiles → `FINAL.json`, `assemble.log` |

All python scripts insert their own directory on `sys.path`, so they can be run
from anywhere:

```
python3 /home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/calib.py
python3 /home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/assemble.py
```

Headline: for $\zeta_5(3)$ the bound comes out at $3.9727<4=m$ (rearrangement
numerator) resp. $3.9504<4$ (Bost–Charles), i.e. a contradiction; see the report
§9 for the ledger of what is certified and what is assumed.
