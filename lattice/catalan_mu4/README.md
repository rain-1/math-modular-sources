# lattice/catalan_mu4 — Catalan from the dilogarithm module at fourth roots of unity

Companion scripts for `consolidation/CATALAN_MU4.md`.  Python 3 + mpmath + numpy.

| script | what it does |
|---|---|
| `calib.py` | reproduces CDT's tau(b;e) for their Theorem A, their `logsmain`, and their §7 P^1-{0,1,oo} remark |
| `mu4_series.py` | exact rational series of the mu_4 inventory; measured denominator types; fold-regularity check for H - G*A |
| `mu4_geom.py` | conformal ceilings (mu_4 template |phi'(0)| = 16^{1/d}), Schwarz-Pick bounds, fold preimage radii |
| `mu4_bc.py` | Bost-Charles numerator BC(phi) for the three hosts on concentric contours (sanity: BC(rho z) = log rho) |
| `mu4_tau.py` | tau^flat / tau^sharp for the measured and relaxed descriptions; like-for-like level-8 comparison |
| `mu4_indep.py` | Q(x)-independence rank check mod 2^61-1 |
| `mu4_table.py` | the entry / margin ledger |

Run order: `calib.py`, `mu4_series.py`, `mu4_geom.py`, `mu4_bc.py`, `mu4_indep.py`, `mu4_table.py`.
