# lattice/adelic_holonomy — the adelic arithmetic holonomy bound

Companion to `consolidation/ADELIC_HOLONOMY.md`.  Builds on `lattice/cdt_finder/`
(imported directly: `cdt_bound.tau_flat`, `tau_sharp`, `I_uvw`).

| file | what it does | log |
|---|---|---|
| `adelic_bound.py` | the bound itself: `gamma_p(m,slopes,p)` and `adelic(...)`.  Self-test = CDT calibration + the scale-invariance (product-formula) check | `calibration.log` |
| `final_table.py` | the level-8 Catalan and X_1(5) Sym^2 margin tables of §4–5 | `final_table.log` |
| `pure_2adic.py` | v_2 of the y-expansion of Sym^{+/-} Li_j(4x) on the s=1/4 host: verifies the pure module's 2-adic slope 2 | `pure_2adic.log` |
| `rowE_sources.gp` | rebuilds t(q), F(q) on the level-8 host from Phi = F D_q t, then the companions B_Psi for every chi_{-4} class; checks B_{Phi_E} = the census b_n, measures the archimedean limits and v_2 | (run with N=32; ~2 min) |
| `hosts_padic.gp` | 2-adic slope of row E (-> 5) and the p-adic slopes of Zagier D (c=-1, all zero) | `hosts_padic.log` |

Run: `python3 lattice/adelic_holonomy/final_table.py` (from the repo root),
`gp -q lattice/adelic_holonomy/hosts_padic.gp`.
