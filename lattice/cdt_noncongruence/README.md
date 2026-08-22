# `lattice/cdt_noncongruence` — the CDT bound on second-order (non-congruence) Apéry rows

Companion to `consolidation/CDT_NONCONGRUENCE.md`.  Calibrated on `../cdt_finder/`
(CDT's own $L(2,\chi_{-3})$ numbers) and `../adelic_holonomy/`.

| script | what it does |
|---|---|
| `hosts_nc.py` | the eleven second-order rows: $\lambda_{1,2}$, $\delta=\lambda_1\lambda_2$, $\lambda_2^{\rm norm}$, $|x_2|$, score |
| `verify_rows.gp` | PARI: $a_n\in\mathbf Z$, $k=2$ sharp, characteristic roots, $\xi$ to 30 digits, $n\le160$ |
| `geom_denom.py` | Lemma G: rescaling $t\mapsto t/\lambda$ $\equiv$ a uniform $p$-adic slope $-v_p(\lambda)$; margins agree to $5\cdot10^{-15}$ |
| `arch_k.py` | architectures (K)/(D): admissible inventories and their $\tau^\flat,\tau^\sharp$ |
| `bc_multivalent.py` | $\Delta(r)=\mathrm{BC}(\varphi_r)-\log|\varphi_r'(0)|$ for $\varphi_r=x_2\lambda(rz)$ |
| `delta_table.py` | fine $\Delta(r)$ grid ($r=0.20\ldots0.80$) $\to$ `delta_table.json`, plus a convergence check |
| `optimum.py` | optimise $(m,\text{inventory},r)$; the thresholds $-0.71963$, $-0.83852$, $-1.10315$ |
| `entries.py` | entry conditions in `CDT_FINDER.md` presentation (entryC/entryR) |
| `table_nc.py` | margins, deficits, the final ranking; calibration $+0.0053$ on CDT's host |
| `ingredients.py` | value of $c$ fold-regular classes, of a $p$-adic slope, of a Galois trace |
| `index910.py` | index-9/10 non-congruence groups; the window $-0.8385<\text{score}\le0$ |

Run order: `delta_table.py` (slow, ~40 min) then anything else.  `*.log` are the recorded runs.
