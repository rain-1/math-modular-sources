# lattice/ — p-adic alignment experiments (PARI/GP)

Run any file with `gp -q <file>` (each ends with `\q`). All arithmetic exact.

- `zagier_padic.gp` — six Zagier second-order rows: $v_p(a_n)$ bounds, single-row $p$-adic slopes, cross-row determinant valuations.
- `zagier_limits.gp`, `zagier_F.gp` — archimedean limits identified by `lindep`.
- `zagier_CF.gp` — 3-adic coincidence of B, C, F with factors ½, ½, ⅝.
- `third_order.gp` — six third-order sporadics: limits, slopes.
- `zeta3_align.gp` — Domb vs (12,4,16) 2-adic alignment (slope 4n, factor 4/3); Apéry row does not align.
- `level12*.gp`, `fitrec.gp` — the book's level-12 purified $L(2,\chi_{-3})$ row: integral coefficients $1,-7,15,-57,\dots$; Picard–Fuchs order 2 degree 9; two cusps on $|t|=1/3$, no fold (not a decaying partner).

Findings are written up in `../consolidation/THEORY_NOTES_03_lattices.md`.
