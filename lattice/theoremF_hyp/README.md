# lattice/theoremF_hyp — hypotheses (b) and (d) of Theorem F

Report: `consolidation/THEOREM_F_HYPOTHESES.md`.

| script | what it does | log |
|---|---|---|
| `growth.gp`, `growth_run.gp` | profiles $v_p(a_n)$, $n\le3000$, for all census rows at their slope primes | `growth.log` |
| `vlaw.gp`, `vlaw_run.gp` | digit laws, Lucas for the normalised sequence, supercongruence slack, sharp $\lambda_p s_p(n)$ bounds | `vlaw.log` |
| `zudilin_v2.gp` | proves $v_2(a_n)=-4n+2s_2(n)$ for Zudilin's Catalan row, $n\le1200$ | — |
| `descent.gp`, `descent_run.gp`, `descent_twins_run.gp` | Ligozat divisors, $\deg t$ on $X_0(N)$, Atkin–Lehner scalars $t\vert W_Q/t$ and $F\vert_wW_Q/F$ | `descent.log` |
| `aleigen.gp`, `aleigen_run.gp` | oldform Atkin–Lehner eigenvalues of $\Phi$ (weight $w+2$) and $\Theta$ (weight $-w$); confirms Lemma 1 | `descent.log` |
| `exactlaws.gp` | pushes the two exact digit laws ($\mathbf B$ at 3, $\mathbf E$ at 2) to $n\le10000$ | — |
| `zeta5_check*.gp` | the level-12 $\zeta(5)$ cautionary example: (a) holds, no 2-adic slope | `zeta5_check*.log` |
| `twins.gp`, `twins_run.gp` | the level-lowering $q\mapsto-q$ twins of rows B (level 9), $\eta$ (level 5), $\delta$ (level 6) | `twins.log` |
