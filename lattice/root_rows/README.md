# lattice/root_rows — the general w-th root construction

Companion code for `consolidation/ROOT_ROWS.md` and `paper/sections/02_root_rows_draft.tex`.

| script | what it establishes |
|---|---|
| `01_lambda_w.py` | Theorem R1: `v_p(binom(1/w,k))`, minimality of `lambda_w = w*rad(w)`, direct series test |
| `02_descent.py`  | Theorem R2: `Sym^w(th^2+p th+r)` for w=2..5 and the identity for `L_w(g^w)` |
| `lib.gp`         | `fitrec`, `minrec`, `rootrow` (minimal lambda + a_n), `charroots`, `emin` |
| `03_sporadic.gp` | the nine w=2 sporadic rows re-verified to n=420 |
| `04_build.gp`    | 400 exact `A_n` for every w>=3 system -> `rows_*.txt` |
| `05_analyse.gp`, `06_run.gp` | lambda / integrality / first recurrence attempt |
| `07_scan.gp`, `08_widescan.gp` | mod-(2^61-1) (order,degree) kernel scans |
| `09_zeta7_s.gp`  | zeta(7) level 24 in the s-coordinate (w=6): order-4 recurrence, k=2 |
| `10_oriented.gp` | beta(4) -> Zagier E; L(4,chi_-3) -> 2F1(1/3,1/3;1;-27 t0) |
| `11_companion.gp`| companion from the order-2 differential operator; `d_n^k b_n` tests |
| `12_zeta5_zero.gp` | the zero of F that breaks hypothesis (H3) for zeta(5) level 16 |
| `13_zeta7_convergence.gp` | the n^{-1/2} convergence that breaks (H4) for zeta(7) |

`rows_*.txt` hold 401 exact integers `A_0..A_400` per system (regenerate with `04_build.gp`).
Run everything from the repository root, e.g. `gp -q -s 8G lattice/root_rows/03_sporadic.gp`.
Note: `read()`ing a file that calls `default(parisize,...)` silently discards the
definitions — pass the stack size on the command line instead.
