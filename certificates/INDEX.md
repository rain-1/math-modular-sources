# Certificate intake index

Log of files River has dropped and what was done with them. One row per artefact.
See certificates/WANTED.md for what's still outstanding.

| Date | File(s) received | Filed to | Strand | Ledger ID / claim | Ran? | Result |
|---|---|---|---|---|---|---|
| 2026-08-21 | `padic-apery-limits.tex` | `certificates/padic/` | padic | p-adic Apéry limits programme: tower limits $\Lambda_a=\lim_s\chi(p)^sp^{ws}B(ap^s)/A(ap^s)$ at **unramified** $p\ge5$; twisted descent law $(\mathrm{LB}_w^\chi)$; exclusion theorem; identification via $h_p(a)$ + Kazandzidis limits | latexmk | PASS, 31 pp. **Directly relevant — see analysis note below.** |
| 2026-08-21 | `main(20260813-132000).tex` | `certificates/padic/` | padic | The fifteen sporadic Apéry-like pairs (canonical table) | latexmk | PASS, 21 pp |
| 2026-08-21 | `main(7).tex` | `certificates/padic/` | padic | A two-level digit law for Apéry's pair | latexmk | PASS, 5 pp |
| 2026-08-21 | `main(10).tex` | `certificates/eisenstein/` | eisenstein | Modular anchors for Apéry-like companions | latexmk | PASS, 9 pp |
| 2026-08-21 | `main(20260813-131951).tex` | `certificates/eisenstein/` | eisenstein | Eisenstein sources of the sporadic companions (E01/E02) | latexmk | PASS, 12 pp |
| 2026-08-21 | `main(8).tex` | `certificates/harmonic/` | harmonic | Harmonic Jets of Hypergeometric Sums (H05) | latexmk | PASS, 20 pp |
| 2026-08-21 | `main(4).tex` | `certificates/harmonic/` | harmonic | Harmonic companion rows for binomial powers | latexmk | PASS, 16 pp |
| 2026-08-21 | `CATALAN_LITERAL_WINDOW_PAPER.tex` | `certificates/catalan/` | catalan | literal-window construction-quality bound 0.581983… | latexmk | PASS, 12 pp |
| 2026-08-21 | `CATALAN_LITERAL_WINDOW_PROOF.tex` | `certificates/catalan/` | catalan | long-form proof companion to the above | latexmk | PASS, 15 pp |
| 2026-08-21 | `frobenius.tex` | `certificates/misc/` | misc | Frobenius constants / middle-root obstruction (*may belong in padic — cross-ref*) | latexmk | PASS, 26 pp |
| 2026-08-21 | `main(3).tex` | `certificates/misc/` | misc | Why Apéry is alone | latexmk | PASS, 14 pp |
| 2026-08-21 | `main(5).tex` | `certificates/misc/` | misc | The companion problem, narrative/provenance survey | latexmk | PASS, 8 pp |
| 2026-08-21 | `main(6).tex` | `certificates/misc/` | misc | A manufactured Apéry apparatus for a cusp-form L-value | latexmk | PASS, 13 pp |
| 2026-08-21 | `main(9).tex` | `certificates/misc/` | misc | Foundations of the Λ-algebra | latexmk | PASS, 21 pp |

Originals left untouched in `certificates/documents/`. Nothing overwritten, nothing deleted.

## Intake note — 2026-08-21

**All fourteen items listed under WANTED P0 §1 are now received** (the header there says
"15 canonical TeX sources" but enumerates 14; either the count is off by one or one source
is unlisted — worth confirming with River).

All fourteen compile standalone under `latexmk -pdf` with no missing inputs, so they are
self-contained certificates, not fragments.

### Cross-checks against work done this session (Fable's `consolidation/ZETA3_TWO_LATTICE.md`)

`padic-apery-limits.tex` Table `tab:families` (the canonical fifteen) **independently
confirms three things derived from scratch this session**:

* family $(\varepsilon)$ is $A(n)=\sum_k\binom nk^2\binom{2k}n^2$, Apéry limit $7\zeta(3)/32$,
  weight $w=3$ — i.e. exactly the "T $=(12,4,16)$" row, matching the binomial form,
  the limit and the denominator exponent $k=3$ found here.
* family $(\alpha)$ is Domb, limit $7\zeta(3)/24$, $w=3$ — matching.
* $(\alpha)$, $(\gamma)$ (Apéry, $\zeta(3)/6$) and $(\varepsilon)$ are the **only** $\zeta(3)$
  families among the canonical fifteen. This upgrades the census conclusion
  ("Domb+T is the only viable pair") from "only pair my search found" to
  "only pair that can exist among the canonical families".

### Caution flagged for the theory notes

`Finding \ref{find:cross}` ("no cross-family relations", all pairs of the fifteen,
$p=7$, 12–21 digits, nothing found) **bounds** the realisation-rigidity conjecture of
`ZETA3_TWO_LATTICE.md` §7. It does not refute it — the objects differ (see the analysis
appended to §7.2 of that file) — but it does kill the strong form. Recorded there.
