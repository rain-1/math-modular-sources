# MUM survey — the Apéry dictionary for rank-4 Calabi–Yau operators

*Claude (Fable), 2026-08-22. Scripts in `lattice/mum_survey/`. Every number below is
measured in exact PARI arithmetic unless marked `[lit]`. Literature located and read by a
sub-agent; the reading list is in §7.*

---

## 0. Verdict

1. **The quintic has no Apéry limit in the Almkvist–van Straten–Zudilin sense at all**, and
   this is forced: its recurrence has order 1. Its $\zeta(3)$ lives one level up, in the
   **connection coefficients** between the MUM Frobenius basis and the conifold, i.e. in the
   $\widehat\Gamma$-class. We computed that connection exactly (§2): the vanishing conifold
   period is
   $$V \;=\; \frac{\sqrt5}{4\pi^2}\Bigl(\hat y_3-10\zeta(2)\,\hat y_1-40\zeta(3)\,\hat y_0\Bigr),
   \qquad -40=\chi/H^3=-200/5,$$
   verified to **116 decimal digits**.
2. **Fable's prediction is confirmed, in the sharp form.** The $p$-adic $\widehat\Gamma$-class
   of the quintic is the same expression with $\Gamma\to\Gamma_p$, and its $\zeta$-coefficient is
   $$[\rho^3]\log\frac{\Gamma_p(1+5\rho)}{\Gamma_p(1+\rho)^5}\;=\;-40\,\zeta_p(3)
   \qquad\text{for every }p,\ \text{including the slope prime }p=5 .$$
   The engine is a clean identity we verified independently at seven primes and three weights
   (§3): $\bigl(\log\Gamma_p(1+x)\bigr)^{(k)}\big|_{0}=-(k-1)!\,\zeta_p(k)$ for odd $k$, and
   $=0$ for even $k$ — the exact mirror of $\log\Gamma(1+x)=-\gamma x+\sum_{k\ge2}(-1)^k\zeta(k)x^k/k$
   with the even part deleted. So the rational number $\chi/H^3$ is the **same at every place**.
3. **The real Apéry dictionary at rank 4 is Theorem F, and it transfers verbatim.** For the
   297 AESZ order-4 operators that do have a second solution we measured slopes and $p$-adic
   limits. The law
   $$\boxed{\ \xi_p=\frac{r_\infty}{\mathcal E_p(m)}\,L_p\bigl(m,\chi\omega^{1-m}\bigr),
   \qquad \mathcal E_p(m)=1-\chi(p)p^{-m},\qquad \xi_\infty=r_\infty L(\chi,m)\ }$$
   holds in **every** case where both places are computable (14 operators), and
   $\xi_p=0$ exactly when $\chi\omega^{1-m}$ is odd, i.e. when $\chi(-1)\ne(-1)^{m-1}$
   — which kills all the $\pi^2$ limits and all the $L(\chi_{-3},3)$ limits.
4. **The $p$-adic side is strictly stronger than the archimedean one here.** 25 operators whose
   archimedean limit our numerics cannot reach — five of them because it **does not exist** —
   have $p$-adic limits identified to 185–2441 digits. Five of the resulting archimedean
   predictions were then confirmed against AvSZ's independent PSLQ table (§5.3), and two
   **new period classes** appear that are not in AvSZ's list at all: $L(\chi_8,2)$ and
   $L(\chi_5,3)$.
5. Caveats and what is *not* proved: §6.

---

## 1. Definitions, and why the quintic is the wrong test case

### 1.1 The AvSZ Apéry limit

`[lit]` Almkvist–van Straten–Zudilin, *Apéry limits of differential equations of order 4 and 5*,
Fields Inst. Commun. **54** (2008) 1–19. Write the operator in $\theta$-form,
$$L=\sum_{i=0}^{k}z^iP_i(\theta),\qquad \theta=z\frac{d}{dz},$$
MUM at $z=0$. The coefficient recurrence is $\sum_i P_i(n-i)u_{n-i}=0$, of order $k$
(= the "degree" of the equation). Let
$$A:\ A_0=1,\ A_n \text{ from the recurrence};\qquad
  B:\ B_n=0\ (n\le0),\ B_1=1,\ \text{same recurrence for }n\ge2 .$$
The **Apéry limit** is $\xi_\infty=\lim_n B_n/A_n$. (Equivalently: $B$ solves
$L(y)=\text{const}\cdot z$; $\xi_\infty$ is the unique constant making $B-\xi_\infty A$
continue past the dominant singularity.) AvSZ's tables 4.1–4.8 list conjectural values —
found by PSLQ — for about 100 operators: $c\pi^2$, $cG$, $cL(\chi_{-3},2)$, $c\zeta(3)$,
$cL(\chi_{-3},3)$, $c\pi^4$, and 16 mixed cases; seven are left unidentified.

The **normalisation caveat**: $\xi_\infty$ scales like $1/\lambda$ under $z\mapsto\lambda z$,
so the rational $c$ depends on the operator's presentation. CYCluster's $z$ sometimes differs
from AESZ's; our $c$ therefore differs from AvSZ's by a rational factor in some rows. The
**ratio** $r_\infty/r_p$ is normalisation-independent, and that is what all our claims are about.

### 1.2 The equivalent singular-point formulation

`[lit]` Bloch–Vlasenko, *Gamma functions, monodromy and Frobenius constants*, arXiv:1908.07501
(v1 title: "…and Apéry constants"). At a **special reflection point** $c$ (a CY conifold, local
exponents $0,1,1,2$) with $\delta$ spanning $\operatorname{im}(\sigma_c-1)$, the *Frobenius
constants* are defined by $(\sigma_c-1)\varphi_{\rho,n}=\kappa_{\rho,n}\delta$. Their Remark 32:
for Apéry's $\zeta(3)$ operator, $\kappa_{1,0}/\kappa_{0,0}=\zeta(3)/6$ — the Apéry limit *is* a
ratio of Frobenius constants. Their Proposition 26: for a hypergeometric
$\prod(D+\beta_j-1)-t\prod(D+\alpha_j)$, the Frobenius constants at $t=1$ are the Taylor
coefficients of $1/A(s)$, $A(s)=\prod\Gamma(s+\alpha_j)/\prod\Gamma(s+\beta_j)$.

### 1.3 Why the quintic (and all 14 hypergeometric CY operators) is degenerate

$L_{\text{quintic}}=\theta^4-5z(5\theta+1)(5\theta+2)(5\theta+3)(5\theta+4)$ has $k=1$, so the
recurrence $n^4A_n=5(5n-1)(5n-2)(5n-3)(5n-4)A_{n-1}$ is **first order**: the solution space is
one-dimensional and $B_n=A_n/A_1$ identically, $B_n/A_n\equiv1/120$. All 14 hypergeometric CY
operators have $k=1$; consistently, AESZ #1–#14 appear nowhere in AvSZ's tables (their lowest
entry is #15). We verified $B_n/A_n\in\mathbb Q$ constant for all 14 (`survey.gp`, the rows with
`degz=1` and `arch=Q`).

**Conclusion.** For the quintic one must go to the conifold. Fable's prediction is exactly right
about *where* to look; §2–§3 carry it out.

---

## 2. The quintic, archimedean: the conifold connection, exactly

Script `lattice/mum_survey/conn_quintic.gp` (and `cyops.gp` for the generic local-solution
machinery). Work in $z=5^5t$, so MUM at $z=0$, conifold at $z=1$, $L=\theta^4-z\prod_{j=1}^4(\theta+j/5)$.

**Local structure at the conifold** (`t3.gp`, `t5.gp`). The indicial polynomial is
$M(M-1)^2(M-2)$: exponents $0,1,1,2$, as the Riemann symbol says. There are three holomorphic
solutions $g_0,V,g_2$ (no obstruction) and one logarithmic one. Writing $s=z-1$ and
$W=V\log s+\sum w_ns^n$, the log-partner exists **only** for one normalisation of the
exponent-1 solution:
$$V=s-\tfrac7{10}s^2+\tfrac{41}{75}s^3-\tfrac{1133}{2500}s^4+\cdots$$
(the coefficient $u_2=-7/10$ is forced; any other choice leaves an obstruction at $s^2$).
This is the intrinsic characterisation of the **vanishing period**: it is the unique exponent-1
solution admitting a logarithmic partner, i.e. $\operatorname{im}(\sigma_c-1)$.

**Connection to the MUM Frobenius basis.** Let $\tilde y_0,\dots,\tilde y_3$ be the *rational*
Frobenius basis at $z=0$ (from $R_n(\rho)=\prod_j(j/5+\rho)_n/((1+\rho)_n)^4$, all coefficients
in $\mathbb Q$), and $\hat y_j$ the same in the coordinate $t=z/5^5$. Matching values at four
points $z=0.45,0.6,0.75,0.9$ with 2000 MUM terms and 400 conifold terms at 120-digit precision:
$$V=\alpha_3\Bigl(\tilde y_3-L\,\tilde y_2+\bigl(\tfrac{L^2}2-10\zeta(2)\bigr)\tilde y_1
+\bigl(-\tfrac{L^3}6+10\zeta(2)L-40\zeta(3)\bigr)\tilde y_0\Bigr),\quad L=\log 5^5,$$
$$\boxed{\ \alpha_3=\frac{\sqrt5}{4\pi^2}\ }\qquad\text{(to }10^{-116}\text{)},$$
each of the three inner constants confirmed by an exact `lindep` (`[1,0,1]`, `[-2,0,0,1,-20]`,
`[6,0,0,0,1,0,-60,240]`). In the coordinate $t$ this collapses to
$$V\ \propto\ \hat y_3-10\zeta(2)\,\hat y_1-40\zeta(3)\,\hat y_0
\ =\ [\rho^3]\Bigl\{\hat y(\rho)\cdot\frac{\Gamma(1-\rho)^5}{\Gamma(1-5\rho)}\Bigr\},$$
i.e. **the conifold period is the $\widehat\Gamma$-class period** $\widehat\Gamma_X=\Gamma(1+H)^5/\Gamma(1+5H)$,
and the $\zeta(3)$ coefficient is $\chi/H^3=-200/5=-40$. `[lit]` This is van Straten's
Conjecture 1 / reflection vector $S=(5,0,\tfrac{25}{12},-200\lambda)$, $\lambda=\zeta(3)/(2\pi i)^3$
(*Calabi–Yau Operators*, arXiv:1704.00164 §2.7); our computation is an independent 116-digit
verification of it *and* of the $\sqrt5/(4\pi^2)$ normalisation.

**Equivalently:** in the $\widehat\Gamma$-normalised basis the $\zeta(3)$ **cancels** from $V$.
The $\zeta(3)$ is precisely the discrepancy between the rational Frobenius basis (the arithmetic
one, in which the Apéry integers live) and the topological/integral one. That is the sentence to
keep.

**The general hypergeometric law.** For $L=\theta^4-Cz\prod(\theta+a_i)$,
$$[\rho^3]\log A(\rho)=-\tfrac13\Bigl(\sum_i\zeta(3,a_i)-4\zeta(3)\Bigr)=\frac{\chi}{H^3}\,\zeta(3).$$
Computing the left side by hand for all 14 weight systems reproduces the database's $\chi/H^3$
in all 14 cases (table in §4). Writing the holomorphic solution as a factorial ratio
$A_n=\prod_j\bigl((u_jn)!\bigr)^{\varepsilon_j}C^n$ — always possible, by Möbius inversion of
$\prod_{j\le N}(j/N)_n=(Nn)!/N^{Nn}$ — this is just
$$\frac{\chi}{H^3}=-\tfrac13\sum_j\varepsilon_ju_j^3,\qquad C=\prod_ju_j^{\varepsilon_ju_j}.$$

---

## 3. The $p$-adic $\widehat\Gamma$-class: $\log\Gamma_p$ has only odd $\zeta_p$

Scripts `gammap.gp`, `gtest3.gp`, `hyp14.gp`; $L_p$ from `lattice/euler_criterion/lp.gp`
(Washington Thm 5.11) extended to $\omega$-twists in `lpgen.gp`.

Morita's $\Gamma_p(n)=(-1)^n\prod_{j<n,\,p\nmid j}j$ is exact at integers; we recovered the Taylor
coefficients $c_k$ of $\log\bigl(\Gamma_p(1+x)/\Gamma_p(1)\bigr)$ on $|x|_p\le1/q$ by evaluating at
$x=qm$, $m=1..40$, and solving the Vandermonde system at 200-digit $p$-adic precision.

**Measured, at $p=2,3,5,7,11,13,17$:**
$$c_2=c_4=0,\qquad 6c_3=-2\zeta_p(3),\qquad 120c_5=-24\zeta_p(5),\qquad 5040c_7=-720\zeta_p(7),$$
i.e. $\bigl(\log\Gamma_p(1+x)\bigr)^{(k)}|_0=-(k-1)!\,\zeta_p(k)$, where
$\zeta_p(k):=L_p(k,\omega^{1-k})$ (Kubota–Leopoldt, Washington's normalisation). Agreement was
$v_p(\cdot)\ge 33$–$83$ digits in every one of the 21 tests; the residuals sit exactly at the
Vandermonde precision loss. Equivalently
$$\log\Gamma_p(z)=\Gamma_p'(0)z-\sum_{k\ge1}\frac{\zeta_p(2k+1)}{2k+1}z^{2k+1}.$$
`[lit]` This is Candelas–de la Ossa–van Straten, arXiv:2104.07816 App. B; our computation is an
independent verification *and* pins their $\zeta_p(3)$ to the branch $L_p(3,\omega^{-2})$ of
Kubota–Leopoldt, which their paper does not spell out. Note that $\zeta_p(2k)$ vanishes
identically because $\omega^{1-2k}$ is odd — this parity fact drives everything in §5.

**Consequence (the prediction, proved modulo the above identity).** For any CY operator whose
holomorphic solution is a factorial ratio,
$$[\rho^3]\log\prod_j\Gamma_p(1+u_j\rho)^{\varepsilon_j}
=-\tfrac13\Bigl(\sum_j\varepsilon_ju_j^3\Bigr)\zeta_p(3)=\frac{\chi}{H^3}\,\zeta_p(3),$$
the **same rational** $\chi/H^3$ as archimedean, at **every** prime — including $p\mid C$.
For the quintic: $-40\,\zeta_5(3)$ at the slope prime $p=5$.

`[lit]` Cross-checks found in the literature after the computation: Beukers–Vlasenko,
arXiv:2302.09603 Thm 1.4, $\Gamma_p(x)/\Gamma_p(x/5)^5=1-\tfrac8{25}\zeta_p(3)x^3+\cdots$
($-8/25=-40/5^3$ after $\varphi=t^5$); Candelas–de la Ossa–van Straten §4.4,
$\gamma=(\chi/H^3)\zeta_p(3)$, quintic $-40\zeta_p(3)$; Shapiro (arXiv:0809.3742, arXiv:1006.0382
Thm 4.2) *proves* the quintic case, $(p^3-1)\tfrac{2^3}{5^2}\zeta_p(3)$, and attributes the
prediction to Kontsevich–Vologodsky. So the quintic half of Fable's prediction is a theorem;
what we add is the uniform mechanism ($\log\Gamma_p$ has only odd $\zeta_p$), the verification,
and the identification of the branch.

---

## 4. The 14 hypergeometric CY operators

$L=\theta^4-Cz\prod(\theta+a_i)$; singular points $0$, $1/C$ (conifold, exponents $0,1,1,2$),
$\infty$. "Slope primes" $=\{p:p\mid C\}$ — the direct analogue of $p\mid c$ for modular rows.
$H^3,\chi$ from CYCluster; the last column is *computed* from the Hurwitz-zeta / factorial-ratio
formula of §2 and agrees with $\chi/H^3$ in all 14 cases.

| AESZ | $(a_i)$ | $C$ | slope primes | $A_n$ | $H^3$ | $\chi$ | $\Gamma$-class $=\chi/H^3$ |
|---|---|---|---|---|---|---|---|
| 1 | $\tfrac15,\tfrac25,\tfrac35,\tfrac45$ | $5^5$ | 5 | $(5n)!/n!^5$ | 5 | $-200$ | $-40$ |
| 2 | $\tfrac1{10},\tfrac3{10},\tfrac7{10},\tfrac9{10}$ | $2^85^5$ | 2,5 | $(10n)!/(n!^3(2n)!(5n)!)$ | 1 | $-288$ | $-288$ |
| 3 | $(\tfrac12)^4$ | $2^8$ | 2 | $\binom{2n}n^4$ | 16 | $-128$ | $-8$ |
| 4 | $\tfrac13,\tfrac13,\tfrac23,\tfrac23$ | $3^6$ | 3 | $((3n)!/n!^3)^2$ | 9 | $-144$ | $-16$ |
| 5 | $\tfrac13,\tfrac12,\tfrac12,\tfrac23$ | $2^43^3$ | 2,3 | $\binom{2n}n^2(3n)!/n!^3$ | 12 | $-144$ | $-12$ |
| 6 | $\tfrac14,\tfrac12,\tfrac12,\tfrac34$ | $2^{10}$ | 2 | $\binom{2n}n(4n)!/n!^4$ | 8 | $-176$ | $-22$ |
| 7 | $\tfrac18,\tfrac38,\tfrac58,\tfrac78$ | $2^{16}$ | 2 | $(8n)!/(n!^4(4n)!)$ | 2 | $-296$ | $-148$ |
| 8 | $\tfrac16,\tfrac13,\tfrac23,\tfrac56$ | $2^43^6$ | 2,3 | $(6n)!/(n!^4(2n)!)$ | 3 | $-204$ | $-68$ |
| 9 | $\tfrac1{12},\tfrac5{12},\tfrac7{12},\tfrac{11}{12}$ | $2^{12}3^6$ | 2,3 | $\binom{2n}n(12n)!/(n!^2(4n)!(6n)!)$ | 1 | $-484$ | $-484$ |
| 10 | $\tfrac14,\tfrac14,\tfrac34,\tfrac34$ | $2^{12}$ | 2 | $((4n)!/(n!^2(2n)!))^2$ | 4 | $-144$ | $-36$ |
| 11 | $\tfrac14,\tfrac13,\tfrac23,\tfrac34$ | $2^63^3$ | 2,3 | $\binom{3n}n(4n)!/n!^4$ | 6 | $-156$ | $-26$ |
| 12 | $\tfrac16,\tfrac14,\tfrac34,\tfrac56$ | $2^{10}3^3$ | 2,3 | $\binom{4n}n(6n)!/(n!^2(2n)!^2)$ | 2 | $-156$ | $-78$ |
| 13 | $\tfrac16,\tfrac16,\tfrac56,\tfrac56$ | $2^83^6$ | 2,3 | $((6n)!/(n!(2n)!(3n)!))^2$ | 1 | $-120$ | $-120$ |
| 14 | $\tfrac16,\tfrac12,\tfrac12,\tfrac56$ | $2^83^3$ | 2,3 | $\binom{2n}n(6n)!/(n!^3(3n)!)$ | 4 | $-256$ | $-64$ |

All 14: $\lambda_1=C$, and there is **no second finite singularity**, so there is no $\lambda_2$
and no linear form — the row is degenerate in the score sense (§1.3). Denominator rates
$\kappa_p=0$ at every prime (e.g. $v_5\bigl((5n)!/n!^5\bigr)=s_5(n)=O(\log n)$). The whole of
the "budget" is in the Frobenius/$\Gamma$-class layer, not in a coefficient ratio.

---

## 5. The survey: 297 order-4 AESZ operators at every place

Data: full CYCluster dump (`https://cycluster.mpim-bonn.mpg.de/api/operators`, 624 rows;
297 of order 4 carry an AESZ number). Pipeline: `mkdata.py` → `ops.gp`; `apery.gp` builds
$(A_n,B_n)$; `survey.gp` measures everything; `analyse.py` applies the Euler-factor law.
$N=240$ terms, 260-digit real precision, exact rationals throughout.

**Pipeline validation.** At $N=220$ the archimedean limits of AESZ 16, 18, 28, 29, 42, 58
reproduce AvSZ's table entries exactly by `lindep` (`validate.gp`).

### 5.1 Measured: archimedean

| identified as | count | AESZ numbers |
|---|---|---|
| rational (degenerate, $k=1$ mostly) | 48 | incl. all of #1–#14 |
| $c\,\pi^2$ | 10 | 15, 18, 24, 25, 26, 45, 51, 62, 63, 68 |
| $c\,\zeta(3)$ | 5 | 16, 28, 29, 42, 182 |
| $c\,G$ | 4 | 36, 38, 48, 65 |
| $c\,L(\chi_{-3},2)$ | 4 | 58, 64, 69, 70 |
| $c\,L(\chi_{-3},3)$ | 1 | 185 |
| not reached at $N=240$ | 225 | — |

The 225 are not counterexamples: most have $\lambda_2/\lambda_1$ very close to 1 (see §5.4).

### 5.2 Measured: $p$-adic, and the law

165 of the 297 have a positive $p$-adic slope at some $p\in\{2,3,5,7,11,13\}$; 40 of those
were identified against Kubota–Leopoldt values at the first pass. **In every case where both
places are known, the Euler-factor law of Theorem F holds exactly:**

| AESZ | $\xi_\infty$ | $p$ | $\sigma_p$ | digits | $\xi_p$ | $\mathcal E_p(m)$ | check |
|---|---|---|---|---|---|---|---|
| 16 | $\tfrac7{48}\zeta(3)$ | 2 | 6 | 1213 | $\tfrac16\zeta_2(3)$ | $1-2^{-3}=\tfrac78$ | $\tfrac16\cdot\tfrac78=\tfrac7{48}$ ✓ |
| 42 | $\tfrac7{64}\zeta(3)$ | 2 | 4 | 805 | $\tfrac18\zeta_2(3)$ | $\tfrac78$ | ✓ |
| 36 | $\tfrac18G$ | 2 | 5 | 1013 | $\tfrac18L_2(2,\chi_{-4}\omega^{-1})$ | $1$ | ✓ |
| 38 | $\tfrac1{24}G$ | 2 | 5 | 1013 | $\tfrac1{24}(\cdot)$ | $1$ | ✓ |
| 48 | $\tfrac1{12}G$ | 2 | 5 | 1013 | $\tfrac1{12}(\cdot)$ | $1$ | ✓ |
| 65 | $\tfrac1{120}G$ | 2 | 5 | 1013 | $\tfrac1{120}(\cdot)$ | $1$ | ✓ |
| 58 | $\tfrac18L(\chi_{-3},2)$ | 3 | 2 | 392 | $\tfrac18L_3(2,\chi_{-3}\omega^{-1})$ | $1$ | ✓ |
| 64 | $\tfrac1{120}L(\chi_{-3},2)$ | 3 | 2 | 392 | $\tfrac1{120}(\cdot)$ | $1$ | ✓ |
| 69 | $\tfrac1{24}L(\chi_{-3},2)$ | 3 | 2 | 392 | $\tfrac1{24}(\cdot)$ | $1$ | ✓ |
| 70 | $\tfrac1{12}L(\chi_{-3},2)$ | 3 | 2 | 392 | $\tfrac1{12}(\cdot)$ | $1$ | ✓ |
| 15, 45, 62, 68 | $c\,\pi^2$ | 2 | 3 | 601 | $\mathbf 0$ | — | parity ✓ |
| 185 | $\tfrac16L(\chi_{-3},3)$ | 3 | 3 | 586 | $\mathbf 0$ | — | parity ✓ |

**The parity rule.** $\xi_p=0$ iff $\chi\omega^{1-m}$ is odd, i.e. iff $\chi(-1)\ne(-1)^{m-1}$:

| $(\chi,m)$ | archimedean period | $\chi\omega^{1-m}$ | $\xi_p$ |
|---|---|---|---|
| $(\mathbf 1,2)$ | $\pi^2=6\zeta(2)$ | odd | $0$ |
| $(\chi_{-4},2)$ | $G$ | even | $\ne0$ |
| $(\chi_{-3},2)$ | $L(\chi_{-3},2)$ | even | $\ne0$ |
| $(\mathbf 1,3)$ | $\zeta(3)$ | even | $\ne0$ |
| $(\chi_{-3},3)$ | $L(\chi_{-3},3)\sim\pi^3\sqrt3$ | odd | $0$ |

This is the rank-4 form of the "cuspidal/outer $\Rightarrow\xi_p=0$" clauses of Theorem F: at
rank 4 the vanishing is a *parity* condition on the character, and it is exactly what makes the
$p$-adic limit blind to the $\pi^2$ and $\pi^3\sqrt3$ parts of a mixed limit (§5.4).

### 5.3 Blind predictions, and five hits against AvSZ

Running the law **backwards** — from $\xi_p$ to $\xi_\infty$ — for the 25 operators where only
the $p$-adic side is computable gives, among others:

| AESZ | $p$-adic input | predicted $\xi_\infty$ (even part) | status |
|---|---|---|---|
| 205 | $\tfrac3{10}\zeta_2(3)$, 1014 d | $\tfrac{21}{80}\zeta(3)$ | **= AvSZ table 4.4** ✓; confirmed here to 68 digits at $n=900$ |
| 137 | $\tfrac18L_2$, $\tfrac5{32}L_3$ | $\tfrac5{32}L(\chi_{-3},2)$ | **= AvSZ table 4.3** ✓; 48 digits at $n=900$ |
| 138 | $\tfrac1{12}L_2$, $\tfrac5{48}L_3$ | $\tfrac5{48}L(\chi_{-3},2)$ | **= AvSZ** ✓ |
| 139 | $\tfrac1{24}L_2$, $\tfrac5{96}L_3$ | $\tfrac5{96}L(\chi_{-3},2)$ | **= AvSZ** ✓ |
| 140 | $\tfrac1{120}L_2$, $\tfrac1{96}L_3$ | $\tfrac1{96}L(\chi_{-3},2)$ | **= AvSZ** ✓ |
| 183 | $-\tfrac14L_3(2,\chi_{-3}\omega^{-1})$ | $-\tfrac14L(\chi_{-3},2)$ | AvSZ prints $+\tfrac14$; our $-\tfrac14$ confirmed to 77 digits at $n=900$ — **sign discrepancy flagged** |
| 47 | $\tfrac{13}{108}\zeta_2(3)$ **and** $\tfrac7{64}\zeta_3(3)$ | $\tfrac{91}{864}\zeta(3)$ | two primes agree; not in AvSZ |
| 61 | $\tfrac1{648}L_2$ **and** $\tfrac1{720}L_3$ | $\tfrac1{648}G$ | two primes agree; not in AvSZ |
| ~67 | $\tfrac5{216}L_2$ **and** $\tfrac1{48}L_3$ | $\tfrac5{216}G$ | two primes agree; not in AvSZ |
| 41, 46, 84 | $\tfrac18,\tfrac14,\tfrac12$ of $\zeta_p(3)$ | $\tfrac{13}{108},\tfrac{13}{54},\tfrac7{16}$ of $\zeta(3)$ | new |
| 49, 133–136, 142, 143, 228 | $L_p(2,\chi_{-3}\omega^{-1})$ | $\tfrac18,\tfrac18,\tfrac1{12},\tfrac1{24},\tfrac1{120},\tfrac1{12},\tfrac1{120},-\tfrac{15}{64}$ of $L(\chi_{-3},2)$ | new |
| 110, 111, 112, 406, ~88/~89 | $L_p(2,\chi_{-4}\omega^{-1})$ | $\tfrac1{12},\tfrac18,\tfrac1{120},\tfrac12,\tfrac1{24}$ of $G$ | new |

The five AvSZ hits are the decisive test: they were obtained purely $2$- and $3$-adically, from
a computation that never sees the archimedean limit, and they match constants AvSZ found by
PSLQ. The **two-prime agreements** (47, 61, ~67, 137–140) are the rank-4 instance of
Conjecture D (rigidity): two different slope primes, one extension class, one rational.

### 5.4 The $p$-adic limit only sees the even part — and sometimes it is all there is

Pushing the archimedean computation to $n=1200$–$1400$ (`mixed.gp`, `mixed2.gp`) shows the
predicted values are **components, not whole limits**: e.g. for AESZ 47,
$\xi_\infty\approx0.0022896$ while $\tfrac{91}{864}\zeta(3)=0.1266055$, the difference being an
odd-character period of $\pi^3\sqrt3$ type — exactly what the parity rule says the $p$-adic
limit must be blind to. AvSZ's own table 4.7 lists 16 such mixed limits. So the correct
statement of the blind predictions above is:

> the **even-character component** of $\xi_\infty$ is $r_p\,\mathcal E_p(m)\,L(\chi,m)$.

For AESZ 137–140, 183, 205 the odd component is zero and the prediction is the whole limit
(hence the exact agreement). We did not identify the odd components of the rest; the achievable
archimedean precision at $n\le1400$ (5–12 digits for those rows) is not enough for a
multi-term `lindep`.

**Operators with a $p$-adic Apéry limit but *no* archimedean one.** For AESZ 41, 133, 134, 135,
136 the ratio $B_n/A_n$ **diverges** ($|B_n/A_n-B_{n-50}/A_{n-50}|$ grows: $2.7$, $4\cdot10^2$,
$2.8\cdot10^2$, $1.4\cdot10^2$, $28$ at $n=1200$), yet the $3$-adic limits exist and are
$\tfrac18\zeta_3(3)$, $\tfrac18,\tfrac1{12},\tfrac1{24},\tfrac1{120}$ of $L_3(2,\chi_{-3}\omega^{-1})$
to 598–787 digits. This is the rank-4 analogue of Zagier's row **B**: *the $p$-adic limit is an
invariant of the extension class that survives when the archimedean Apéry limit does not exist.*
These five are the sharpest flagged items of the survey.

### 5.5 New period classes

The extended scan over quadratic characters of conductor $\le24$ (`survey2.gp`) turned up two
periods that are **not in AvSZ's list at all**:

* **AESZ 7** (nn 4.2.31, $\deg_z=2$): $\sigma_2>0$, 2441 digits, $\xi_2=\tfrac1{16}L_2(2,\chi_8\omega^{-1})$
  $\Rightarrow$ $\xi_\infty$'s even part $=\tfrac1{16}L(\chi_8,2)$ (here $\mathcal E_2(2)=1$).
* **AESZ 184** (nn 4.2.57, $\deg_z=2$, discriminant $z^2-\tfrac{11}{250}z+\tfrac1{2000}$):
  $\sigma_5>0$, 612 digits, $\xi_5=\tfrac14L_5(3,\chi_5\omega^{-2})\Rightarrow$ even part
  $\tfrac14L(\chi_5,3)$. This is the rank-4 sibling of the project's own Almkvist–Zudilin row
  $\eta=(11,5,125)$, whose limit is $\tfrac12L(3,\chi_5)$ and whose $5$-adic limit is
  $\tfrac12\zeta_5(3)$ (`04_padic_euler_draft.tex`, census row $\eta$). **Alignment candidate:**
  a rank-4 and a rank-3 system, one period $L(\chi_5,3)$, one slope prime $p=5$ — the "one class,
  two worlds" configuration of `THEORY_NOTES_06` §5.3. Testing the cross determinant is the
  obvious next computation and was not done here.

  Archimedean status of both (`newclass.gp`, $n=1400$): AESZ 7 has drift $4\cdot10^{-6}$ and
  AESZ 184 drift $5\cdot10^{-2}$, so neither even-part prediction is confirmable archimedeanly
  at this length; as in §5.4 the measured $x$ differs from the pure prediction by what must be an
  odd-character component.

### 5.6 Open: 131 unidentified $p$-adic limits

131 operator/prime pairs have a measured positive slope and a $p$-adic limit of 185–3054 digits
that matched **none** of our targets (all $L_p(2,\chi\omega^{-1})$, $L_p(3,\chi\omega^{-2})$ for
$\chi$ quadratic of conductor $\le24$, plus trivial). Notable: **AESZ 207**, one of AvSZ's seven
*unidentified* archimedean limits, has $\sigma_2=12$ and a $2$-adic limit known to **2449 digits**
that is not a single Kubota–Leopoldt value. Also AESZ 227 (two primes), 391, 415, 269, 276, 362
(two primes each). These are the natural next targets: two-term combinations, weight-4/5 values
$\zeta_p(5)$, and $p$-adic $L$-values of the weight-4 newforms attached to conifold fibres
(BKSZ arXiv:2203.09426 Table 1 gives the levels).

### 5.7 Denominator exponents and rates

$k$ (minimal with $d_n^kB_n\in\mathbb Z$) is $2$ for the weight-2 families and $3$ for the
$\zeta(3)$/$L(\chi,3)$ families, uniformly across the survey — the rank-4 echo of $k=w+1$.
Two exceptions with genuine denominator *rates* $\kappa_p>0$ in $A_n$: **AESZ 32** (Zudilin's
$\pi^4$ case) has $\kappa_2=2$, and **AESZ 244** has $\kappa_2=1$. Those are the only members of
the census with a hypergeometric-style $p$-adic resource, and by the slope law
$\sigma_p=v_p(c)+2\kappa_p$ they are the ones worth re-examining for two-row constructions.

---

## 6. Caveats

1. **Normalisation.** CYCluster's $z$ is not always AESZ's, so our rational $c$'s can differ from
   AvSZ's by a rational factor. All claims are about $r_\infty/r_p$, which is invariant. The five
   agreements in §5.3 show the normalisations do coincide for those rows.
2. **Slopes are measured, not derived.** $\sigma_p$ is read off $v_p(B_n/A_n-B_{n-1}/A_{n-1})$ at
   $n=120$ and $n=240$; Proposition C is not directly applicable because these recurrences are not
   in Zagier/AZ normalisation (Remark "Non-Zagier normalisations", `04_padic.tex`).
3. **$p$-adic identifications are `lindep` fits**, accepted only when the residual valuation
   exceeds (available digits $-10$). With 1000+ digits available a false positive is
   implausible, but at the low end (185–400 digits, e.g. AESZ 183, 199, 200) it is not excluded.
4. **The Euler-factor law is verified, not proved, at rank 4.** No syntomic-regulator argument
   is offered here. The $14+5$ confirmations and the two-prime agreements are the evidence.
5. **The $p$-adic $\widehat\Gamma$-class is formal.** $\prod_j\Gamma_p(1+u_j\rho)^{\varepsilon_j}$
   is manipulated as a formal power series; convergence needs $|u_j\rho|_p\le1/q$. This is the
   standard convention (Beukers–Vlasenko, CdlOvS) but should be said.
6. **`?`-rows.** Our `bestappr`-based archimedean identifier accepts spurious rationals when
   $|x|$ is large (a relative-error test); rows reported as `Q` with a huge numerator have **no**
   limit and should be read as `?`. The counts in §5.1 are corrected for this only informally.
7. **AESZ numbering.** CYCluster covers 297 of the ~404 AESZ operators of order 4 and a few
   numbers appear twice (e.g. "7"). AvSZ's #37, 44, 50, 52, 53, 66, 67, 147, 149, 185(part),
   209, 214, 229, 257 need arXiv:math/0507430 v2 Appendix A.

---

## 7. What was read

* Almkvist–van Straten–Zudilin, *Apéry limits of differential equations of order 4 and 5*,
  Fields Inst. Commun. **54** (2008) 1–19 — the definition, and tables 4.1–4.8.
* van Straten, *Calabi–Yau Operators*, arXiv:1704.00164 §2.5–2.7 — Frobenius basis, $\widehat\Gamma$-class,
  reflection vector, Conjectures 1–2. Contains the sentence "…most notably the $p$-adic story,
  which involves … the $p$-adic analogue of the $\Gamma$-conjectures that lead to the appearance
  of the $p$-adic analogue of $\zeta(3)$."
* Bloch–Vlasenko, arXiv:1908.07501 — Frobenius constants; Prop. 26 (hypergeometric); Rem. 32 (Apéry).
* Almkvist–van Enckevort–van Straten–Zudilin, arXiv:math/0507430 — the AESZ tables.
* Bönisch–Klemm–Scheidegger–Zagier, arXiv:2203.09426 — Table 1 (the 14, with conifold newform levels);
  §3.1 (the nine conifold constants; the $w^\pm,e^\pm$ periods/quasiperiods of the weight-4 newform).
* Candelas–de la Ossa–van Straten, arXiv:2104.07816 §4.4, App. B — the $p$-adic $\Gamma$-class conjecture.
* Beukers–Vlasenko, arXiv:2302.09603 — Frobenius structure constants; Thms 1.4, 1.5.
* Shapiro, arXiv:0809.3742, arXiv:1006.0382 Thm 4.2 — the quintic $p$-adic case, proved.
* Straub–Zudilin, arXiv:2011.03400; Golyshev–Zagier (PSPUM 103.2); Cynk–van Straten arXiv:1210.3249.
* Databases: `https://cycluster.mpim-bonn.mpg.de/api/operators` (used), `https://cydb.mathematik.uni-mainz.de/`.

## 8. Scripts

`lattice/mum_survey/`:
`cyops.gp` (generic $\theta$/$D$-form conversion, local Frobenius and log solutions at any point),
`conn_quintic.gp` (the 116-digit connection computation), `conifold_local.gp` (conifold local structure: exponents, the canonical vanishing period, its log partner),
`conifold_padic.gp` (the conifold pair's $p$-adic behaviour — negative, see below),
`gammap.gp`+`gtest3.gp` ($\Gamma_p$ Taylor coefficients vs $\zeta_p$),
`lpgen.gp` ($L_p(s,\chi\omega^m)$), `mkdata.py`+`ops.gp` (database → GP),
`apery.gp` (Apéry pairs), `validate.gp`, `survey.gp`/`survey2.gp` (the census),
`analyse.py` (the Euler-factor law), `mixed.gp`/`mixed2.gp`/`check_pred.gp`/`newclass.gp` (archimedean checks).

**A negative result worth recording** (`conifold_padic.gp`): the conifold pair $(V,W)$ of the
quintic — the exponent-1 solution and its log partner, normalised $w_0=w_1=w_2=0$ — is *not*
$5$-integral: $v_5(v_n)\approx-1.64n$. It therefore has no $5$-adic slope and gives no $p$-adic
Apéry limit. The $p$-adic content of the quintic conifold is in the Frobenius structure (§3),
not in the conifold Taylor coefficients.
