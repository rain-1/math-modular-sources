# Theory notes 03 — lattice match-ups: the $p$-adic slope law and extension-class rigidity

*Claude (Fable), 2026-08-21. Computations in `lattice/*.gp` (PARI, exact rationals). Tags as before.*

## 1. Where the $p$-adic gain comes from (exact)

Take a second-order row in Zagier's normalisation
$$(n+1)^2u_{n+1}=(an^2+an+b)\,u_n-c\,n^2u_{n-1},\qquad a_0=1,\ a_1=b,\qquad b_0=0,\ b_1=1 .$$
The Casoratian is exact:
$$a_nb_{n-1}-a_{n-1}b_n=-\frac{c^{\,n-1}}{n^2},$$
hence
$$v_p\!\Big(\frac{b_n}{a_n}-\frac{b_{n-1}}{a_{n-1}}\Big)=(n-1)\,v_p(c)-2v_p(n)-v_p(a_n)-v_p(a_{n-1}). \tag{1}$$
So if $v_p(a_n)=O(\log n)$ (true for all six Zagier rows: max $v_2(a_n)\le 22$, $v_3\le13$ for $n\le260$), the ratio $b_n/a_n$ is $p$-adically Cauchy with **slope $\sigma_p=v_p(c)$** and converges to some $\xi_p\in\mathbb Q_p$. [proved, trivial given (1)]

Since $c=\lambda_1\lambda_2$ is the product of the characteristic roots,
$$\boxed{\ \log|\lambda_1\lambda_2|=\sum_p\sigma_p\log p\ } \tag{2}$$
— the product formula. **Archimedean badness is exactly the sum of the $p$-adic slopes.** Apéry's rows have $c=\pm1$: nothing to harvest, perfect archimedean decay. Catalan's E row has $c=32$: $\lambda_2=4>1$, and a 2-adic slope of 5 to compensate. The two-lattice method is the device that cashes in (2).

Heuristic **budget**: if all $p$-adic slope could be harvested, a row's effective score would be $\log\lambda_1-k$ instead of $\log(1/\lambda_2)-k$:

| row | period | $c$ | $\lambda_1$ | $\lambda_2$ | $\sigma_2$ | $\sigma_3$ | score $=\log(1/\lambda_2)-2$ | budget $=\log\lambda_1-2$ |
|---|---|---|---|---|---|---|---|---|
| A | $\zeta(2)/4$ | $-8$ | 8 | $-1$ | 3 | 0 | $-2$ | $+0.079$ |
| B | $L(2,\chi_{-3})/2$ (no arch. limit) | 27 | $\sqrt{27}$ | $\sqrt{27}$ | 0 | 3 | $-3.65$ | $-0.35$ |
| C | $L(2,\chi_{-3})/2$ | 9 | 9 | 1 | 0 | 2 | $-2$ | $+0.197$ |
| D | $\zeta(2)/5$ | $-1$ | 11.09 | $-0.090$ | 0 | 0 | $+0.41$ | $+0.41$ |
| E | $G/2$ | 32 | 8 | 4 | 5 | 0 | $-3.39$ | $+0.079$ |
| F | $5L(2,\chi_{-3})/8$ | 72 | 9 | 8 | 3 | 2 | $-4.08$ | $+0.197$ |

[verified] slopes measured numerically to $n=260$ match $v_p(c)$ to within $O(\log n)/n$; limits identified by `lindep` (A, C, D, E to 90 digits, F to 30). The Catalan construction reaching quality $1-\varepsilon$ is consistent with E's budget being (barely) positive. Whether "budget" is the right invariant is [conj]; the $m=2$ quality formula $Q$ would decide.

## 2. Extension-class rigidity of the $p$-adic limit (the main finding)

Cross-row test: $v_p(r\,a^X_nb^Y_n-r'\,a^Y_nb^X_n)$ with the archimedean factors $r,r'$ built in.

| pair | common period | $p$ | $v_p$ at $n=100,200,300$ | slope |
|---|---|---|---|---|
| B,C | $\tfrac12L(2,\chi_{-3})$ | 3 | 199, 398, 600 | $2n=\min(3,2)\,n$ |
| C,F ($4\xi^F=5\xi^C$) | $L(2,\chi_{-3})$ | 3 | 198, 401, 601 | $2n$ |
| B,F ($4\xi^F=5\xi^C$) | $L(2,\chi_{-3})$ | 3 | 199, 398, 600 | $2n$ |
| A,E | $\zeta(2)$ vs $G$ | 2 | 4, 4 | none |
| A,F; E,F | different | 2 | 4, 4 | none |

So:

> **Empirical law [verified, 3 positive + 3 negative instances].** If two integral Apéry rows have archimedean limits $r\Lambda$, $r'\Lambda$ for the same period $\Lambda$ and both have positive $p$-adic slope, then their $p$-adic limits satisfy the *same* rational relation, $\xi_p=r\Lambda_p,\ \xi'_p=r'\Lambda_p$, and the cross-determinant valuation grows like $\min(\sigma_p,\sigma'_p)\,n$. Rows with different periods do not align, even at a prime where both have slope.

Two striking corollaries:
- **B has no archimedean Apéry limit** (complex-conjugate roots; $b_n/a_n$ oscillates) **yet its 3-adic limit is exactly C's.** The $p$-adic limit is an invariant of the extension class that survives when the archimedean limit is invisible.
- The alignment prime is not "a prime where slopes exist"; it is a prime where slopes exist *and* the extension class is shared. The conjecture in Notes 02 §3 should be restated: **$p$-adic alignment $\Leftrightarrow$ same extension class + positive slopes at $p$; slopes at $p$ occur iff $p\mid\lambda_1\lambda_2$, i.e. iff the singular values of the Hauptmodul are $p$-adically non-units — which for modular rows happens at primes of the level.**

Conceptually this is what one expects if $\xi_p$ is a $p$-adic $L$-value (Calegari's overconvergence mechanism; the 5:8 paper's $G_2=\zeta_2(2)$): the $p$-adic limit is the $p$-adic realisation of the same Beilinson–Kings extension class. A proof for the Zagier rows would be a clean theorem: *the $p$-adic Apéry limits of B, C, F are $\tfrac12,\tfrac12,\tfrac58$ times one $3$-adic period $L_3(2,\chi_{-3})$-type value.* [open; the method is the level-9/level-6 analogue of Calegari's $p=3$ argument]

## 3. Consequences for engineering

1. **$L(2,\chi_{-3})$ is the calibration case for $m=3$.** Three aligned rows at $p=3$ with slopes (3,2,2) and known irrationality (CDT). But none decays archimedeanly: C has $\lambda_2=1$ exactly (polynomial-size linear form), F has $\lambda_2=8$, B oscillates. The level-12 fully purified row from the book (coefficients $1,-7,15,-57,159,\dots$; `lattice/level12g.gp`) has an order-2, degree-9 Picard–Fuchs operator with **two cusps on the circle $|t|=1/3$ and no fold**, so it is not the decaying partner either. Conclusion: for $L(2,\chi_{-3})$ the $p$-adic side is rich and the archimedean side is poor; a two-lattice attack needs a hypergeometric/cellular row with $\lambda_2<1$, which is not known — consistent with the fact that CDT's holonomy proof was the first. [negative result, useful]
2. **Catalan is special** because Zudilin's hypergeometric row supplies decay ($\lambda_2=\varphi^{-5}$) while the modular E row supplies 2-adic slope. **Cross-world pairs are the template.** Systematic search: for each Dirichlet period with a conductor prime $p$, list (i) hypergeometric rows (Beukers/Rivoal Padé, Rhin–Viola integrals) with $\lambda_2<1$ and (ii) modular rows with $p$-adic slope; test alignment by the cross-determinant.
3. **Odd zeta values — corrected.** I first wrote "Fricke-symmetric rows have $c=\pm1$, so odd zeta has nothing to align". That is true of the *Apéry-perfect* row only. The third-order census (`lattice/third_order.gp`, `zeta3_align.gp`) shows $\zeta(3)$ carried by three rows:

| row | limit | $c$ | $\lambda_1$ | $\lambda_2$ | $\sigma_2$ | budget $\log\lambda_1-3$ |
|---|---|---|---|---|---|---|
| Apéry (17,5,1) | $\zeta(3)/6$ | 1 | 33.97 | 0.029 | 0 | $+0.53$ |
| Domb (10,4,64) | $7\zeta(3)/24$ | 64 | 16 | 4 | 6 | $-0.23$ |
| T $=(12,4,16)$ | $7\zeta(3)/32$ | 16 | 23.3 | 0.686 | 4 | $+0.15$ |

   and **Domb and T share their 2-adic limit with the predicted factor**: $v_2(3a^T_nb^D_n-4a^D_nb^T_n)=4n+O(1)$ ($=\min(6,4)\,n$) checked to $n=400$; Apéry's row does not align with either (slope 0). Other weight-2 rows: $(9,3,-27)\to L(3,\chi_{-3})/3$ with $\sigma_3=3$, $\lambda_2=-1.39$; AZ $(7,3,81)$ and $(11,5,125)$ have complex roots (no archimedean limit) with $\sigma_3=4$, $\sigma_5=3$.

   **So the Domb + T pair is a Catalan-type configuration for $\zeta(3)$**: one row with growing linear form and large 2-adic slope (Domb $\leftrightarrow$ E), one with decaying linear form and smaller slope (T $\leftrightarrow$ Zudilin's row), sampling ratio to balance slopes $6:4$. Since $\zeta(3)\notin\mathbb Q$ is known, running the two-lattice machine on it is a *calibration with a known answer* — and if it crosses 1 it is an independent, second-kind proof of Apéry's theorem, which would validate the method far more convincingly than any conditional Catalan number. The common 2-adic limit should be $\tfrac{7}{24}\zeta_2(3)$-type (Calegari proved $\zeta_2(3)\notin\mathbb Q$ by overconvergence — same mechanism). [open, **highest priority next**]
4. What remains true for $\zeta(5)$: the archive's $\zeta(5)$ rows (level 12, 25/144 and 11/144) should be checked for $c$ and slopes; the "$2+3$ fusion" failure is then explained by different-prime slopes (no common alignment prime), not by absence of slopes.

## 4. Next computations
- Tabulate $c$, slopes, limits for the six third-order sporadics and the three Cooper rows; look for shared periods with slope.
- Write the $m=2$ quality formula $Q(\lambda_1,\lambda_2,k;\lambda'_1,\lambda'_2,k';\sigma_p,\sigma'_p)$ from the 5:8 paper's lattice argument and check it reproduces $0.9025$.
- Prove the $3$-adic coincidence of B, C, F (Calegari-style), as the first instance of the rigidity law.

## 5. Correction (later on 2026-08-21): denominators are a second resource

Section 1 assumed $v_p(a_n)\ge0$. If $a_n$ carries $p$-power denominators at rate $\kappa_p$ (i.e. $v_p(a_n)\approx-\kappa_p n$), the two $-v_p(a_n)$ terms in (1) each contribute $+\kappa_p n$, so
$$\sigma_p = v_p(c) + 2\kappa_p,\qquad \text{budget}=\log\Lambda-k+\sum_p\kappa_p\log p .$$
Zudilin's Catalan row: $c=1$, $v_2(Q_m)=-4m+2s_2(m)$, so $\kappa_2=4$, slope $8$ (measured), budget $+1.18$ — the actual mechanism behind the Catalan $1-\varepsilon$ result, invisible to the product formula (2). Consequences: (i) the "$\log\Lambda$ ceiling" holds for integral rows only; (ii) the design criterion is *a partner with $\kappa_p>0$*, which is what hypergeometric (factorial-normalised) rows supply and modular rows never do — this is the precise content of "cross-world pairs"; (iii) for $\zeta(3)$ the live move is a Rhin–Viola / $\binom{2k}{k}$-type row with $\kappa_2>0$ against Domb as engine. Full derivation: `ZETA3_TWO_LATTICE.md` §14.

## 6. Audit note (2026-08-21, late): what the two-lattice method can and cannot prove

An apparent "$|q_nG-p_n|\to0$" from a box-constrained lattice on the exact Zudilin×Nesterenko rows was audited independently (`CATALAN_AUDIT.md`). Rows are exact; the arithmetic is right; the conclusion is void: replacing $G$ by a rational $G^*$ gives identical numerics because the lattice never uses $G$. Past the knife-edge ($\sigma>E_Z+E_N$, i.e. $F<0$) Minkowski's box equals the covolume, the argument degenerates to Dirichlet, and the short vector may be the rationality-kernel vector. Hence the method proves worthiness $\to1$ (Lean: `catalan_worthiness_one_sub_eps_nesterenko`, all 2-adic inputs proved, axioms standard) and **irrationality iff nonvanishing of the selected linear form** — which is the only open statement, and the one to attack.
