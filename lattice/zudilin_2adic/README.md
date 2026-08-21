# `lattice/zudilin_2adic/` — scripts for `consolidation/ZUDILIN_2ADIC.md`

All computations are exact over $\Q$ (PARI/GP 2.15).  Run from the repo root.

| script | what it does |
|---|---|
| `caso2fit.gp`  | symbolic proof-checks of **Lemma A** (mixed Casoratian, $X^2\Omega_n=-(X^2/n^2-2X/n+2)$), **Lemma D** (skew Casoratian $W_n=X^2(X-1)^2\Omega'_n$, closed form, and $W_n(x_n)=q(n)/(64n^2(n+1)^2)$), and $W_n(1/2)$; verified as polynomial identities in $X$ for $1\le n\le26$ |
| `caso2.gp`     | prints $W_n(X)$ for $n\le16$ (raw data behind the fit) |
| `qvals.gp`     | the closed forms $q_m(0)=1$, $q_m(-1)=2m+1$, $q_m(-2)=3m^2+3m+1$, $q_m'(0)=-H_m$, $q_m'(-1)=-(2m+1)H_m$ used to pin the low coefficients of $W_n$ |
| `verify_bridge.gp` | exact verification of $Q_m=q_m(\frac12-m)$ and of the **bridge identity** $P_m=\frac{(-1)^m}{8}p_m(\frac12-m)+G_mQ_m$ for $1\le m\le160$ |
| `bridge_padic.gp`  | 2-adic check: $v_2\bigl(\hat P_m/Q_m-\zeta_2(2)\bigr)=v_2\bigl(P_m/Q_m-\zeta_2(2)\bigr)$, $D_m=0$ |
| `final_check.gp`   | exact valuations $v_2(Q_m)=-4m+2s_2(m)$, $v_2(r_m(x_m))=4m+2-2s_2(m)$, $v_2(\zeta_2(2)-P_m/Q_m)=8m-1-4s_2(m)$ for $m\le60$ at precision $2^{700}$ |
| `tseries.gp`   | resolves the sign conventions: the Laurent series $\sum(2^{n+1}-2)B_n(-1/x)^{n+1}$ (PARI $B_1=-\frac12$) equals $-\Theta$, and $\Theta_2(1/2)=+8\zeta_2(2)$, to $2^{200}$ |
| `doc_check.gp` | re-verifies **every** closed form asserted in `ZUDILIN_2ADIC.md` (Lemmas A and D, the $q$-value formulas, the derived three-term relation, Zudilin's Casoratian) for $n\le24$ — run this after any edit to the write-up |
