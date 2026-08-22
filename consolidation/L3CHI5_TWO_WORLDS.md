# $L(3,\chi_5)$ at $p=5$: AESZ 184 $\times$ $\eta$, and why the pair is not a cross-world pair

*Claude (Fable), 2026-08-22.  Scripts: `lattice/followups/a1_structure.gp`,
`a2_align.gp`, `a3_padic.gp`, `a5_foldraw.gp`, `a6_design.gp`.
Tags: **[proved]** exact identity; **[verified]** exact rational / $p$-adic computation over a
stated range; **[measured]** floating-point rate.  Follows up `MUM_SURVEY.md` §5.5,
`EULER_CRITERION.md` §4.1, `THEOREM_B_EXACT.md` §5, `THEORY_NOTES_03_lattices.md` §5,
`ONE_CLASS_TWO_WORLDS.md` §5–6, `paper/sections/05_two_row.tex` §5.2–5.3.*

---

## 0. Verdict

1. **AESZ 184 is not a new period class.  It is the row $\eta=(11,5,125)$ multiplied by
   $\binom{2n}{n}$** — i.e. the Hadamard product of $\eta$'s Picard–Fuchs operator with
   $(1-4z)^{-1/2}$:
   $$\boxed{\;A^{184}_n=\binom{2n}{n}A^{\eta}_n,\qquad
     B^{184}_n=\tfrac12\binom{2n}{n}B^{\eta}_n,\qquad\text{hence}\quad
     \frac{B^{184}_n}{A^{184}_n}=\frac12\,\frac{B^{\eta}_n}{A^{\eta}_n}\ \text{for every }n.\;}$$
   **[proved]** (§1.2), **[verified]** exactly for $0\le n\le1200$.  The "$\tfrac14L(\chi_5,3)$
   new class at rank 4" of `MUM_SURVEY.md` §5.5 is therefore $\eta$ in rank-4 clothing, not an
   independent realisation.  This should be corrected in the survey.

2. **The alignment at $p=5$ is exact, not asymptotic.**  $\xi_5^{184}=\tfrac14\zeta_5(3)$ and
   $\xi_5^{\eta}=\tfrac12\zeta_5(3)$ (here $\zeta_5(3):=L_5(3,\mathbf 1)$, see the
   notational warning in §2), each to $5^{2993}$ at $N=1000$ (full available precision),
   and the cross-determinant $D_{n,m}=B^{\eta}_nA^{184}_m-2A^{\eta}_nB^{184}_m$ has
   $v_5(D_{n,m})=3\min(n,m)+O(1)$ on the whole grid, with $D_{n,n}\equiv0$ **identically**. §2–§3.

3. **New: the archimedean complex fold of AESZ 184.**  AESZ 184 has *no* Apéry limit — its two
   characteristic roots are $44\pm8i$, of equal modulus $20\sqrt5$ — but its complex fold
   constant is exactly half $\eta$'s:
   $$\boxed{\;\xi^{184}=\tfrac14L(\chi_5,3)+i\,\tfrac{\pi}{20}L(\chi_5,2)
     =\tfrac12\,\xi^{\eta}.\;}$$
   **[proved]** from (1) plus `THEOREM_B_EXACT.md` §5, and **[measured]** independently to
   $2.4\cdot10^{-6}$ by a *recurrence-only* estimator that never uses modularity (§4).  So the
   even-part prediction $\tfrac14L(\chi_5,3)$ of the $p$-adic survey is confirmed, and it is the
   **real part of a genuine complex constant**, exactly as for $\eta$.  The same estimator
   reproduces `THEOREM_B_EXACT`'s modularly computed $\xi^{\eta}$ to $3.6\cdot10^{-6}$ — a first
   independent check of that value by a completely different method.

4. **The pair is not a cross-world pair and $\delta=0$.**  Both rows are complex-fold rows:
   $|\lambda_2|=|\lambda_1|$, no decay, and the master formula gives $H=F>0$, hence
   $\delta=1-F/H=0$ at every sampling ratio $r$ (minimum $F=2.1931$ at $r=1$).  §5.
   The failure is **not** the one `ONE_CLASS_TWO_WORLDS.md` §6.1 diagnoses: the corrected
   criterion $w+\kappa_p>0$ *is* satisfied, with the **largest net $p$-adic resource in the
   project**, $(w+\kappa_5)\log5=3\log5=4.828$ per index (against Zudilin's Catalan
   $\kappa_2\log2=2.773$ and the $\chi_{-3}$ hypergeometric family's $\le2.5$,
   `ONE_CLASS_TWO_WORLDS.md` §5).  What is missing is
   the archimedean half: **no row carrying $L(\chi_5,3)$ decays, and none has even a
   polynomial-size linear form** — there is no "row $\mathbf C$ analogue" on this line. §5.2.

5. **Sharp design deficit.**  Against the engine $\eta$ ($k=3$, $\sigma_5=3$,
   $\log\rho_2^{\rm eng}=\tfrac32\log5$), a hypothetical $\chi_5$-decayer with denominator rates
   $(\kappa_5,\nu,k_d)$ and Casoratian rate $w$ crosses $\delta=1$ iff, at the balanced ratio
   $r=\sigma_5^{\rm dec}/3$ and assuming $k_d\le\sigma^{\rm dec}_5$,
   $$\boxed{\;\log\tfrac1\lambda\;>\;\bigl(1-\tfrac12\log5\bigr)w\;+\;2\,\kappa_5\;+\;\nu
     \;=\;0.19528\,w+2\kappa_5+\nu.\;}$$
   (The coefficient $2$ is exact: $2(1-\tfrac12\log5)+\log5=2$.)  §5.3.

---

## 1. AESZ 184 reconstructed, and its Hadamard origin

### 1.1 The row

CYCluster `4.2.57`, AESZ 184, $\deg_z=2$ (`lattice/mum_survey/ops.gp` line 99):
$$L_{184}=\theta^4-2z(2\theta+1)^2\bigl(11\theta^2+11\theta+5\bigr)
  +500\,z^2(2\theta+1)(\theta+1)^2(2\theta+3),$$
recurrence (order two, from $\sum_iP_i(n-i)u_{n-i}=0$)
$$n^4u_n=2(2n-1)^2\bigl(11n^2-11n+5\bigr)u_{n-1}-500(2n-3)(n-1)^2(2n-1)u_{n-2}. \tag{1.1}$$
Apéry pair in the AvSZ normalisation ($A_0=1$; $B_0=0$, $B_1=1$):
$$A_n:\;1,\,10,\,210,\,5500,\,159250,\,4852260,\,151466700,\,4755808200,\dots$$
$$B_n:\;0,\,1,\,\tfrac{243}{8},\,\tfrac{103775}{108},\,\tfrac{54787775}{1728},
  \tfrac{103126471}{96},\tfrac{1181533815}{32},\tfrac{7001276746035}{5488},\dots$$
Discriminant $z^2-\tfrac{11}{250}z+\tfrac1{2000}$, i.e. the two finite singularities are
$z_\pm=(11\pm2i)/500$, **complex conjugate**, $|z_\pm|=1/(20\sqrt5)$.

### 1.2 Theorem (Hadamard twist).  **[proved]**

Let $\eta=(11,5,125)$ be the Almkvist–Zudilin row
$(n+1)^3u_{n+1}=(2n+1)(11n^2+11n+5)u_n-125n^3u_{n-1}$ with $A_0=1,A_1=5$, $B_0=0,B_1=1$.  Then
$$A^{184}_n=\binom{2n}{n}A^{\eta}_n,\qquad B^{184}_n=\tfrac12\binom{2n}{n}B^{\eta}_n .$$

*Proof.*  Put $w_n=\binom{2n}{n}u_n$.  Using
$\binom{2n}{n}/\binom{2n-2}{n-1}=2(2n-1)/n$ and
$\binom{2n}{n}/\binom{2n-4}{n-2}=4(2n-1)(2n-3)/\bigl(n(n-1)\bigr)$, multiplying $\eta$'s
recurrence at index $n$ by $\binom{2n}{n}$ turns it into exactly (1.1).  So
$u\mapsto\binom{2n}{n}u$ is an isomorphism of the two solution spaces.  It sends $A^\eta$
(value $1$ at $n=0$) to $A^{184}$, and $B^\eta$ (values $0,1$) to the solution with values
$0,2$, i.e. to $2B^{184}$. $\square$

Equivalently $L_{184}$ is the **Hadamard product** of $\eta$'s Picard–Fuchs operator with the
order-one operator whose solution is $(1-4z)^{-1/2}=\sum\binom{2n}{n}z^n$.  Consequences,
all immediate:

* $\dfrac{B^{184}_n}{A^{184}_n}=\dfrac12\dfrac{B^{\eta}_n}{A^{\eta}_n}$ **for every $n$**
  — at every place at once. **[verified]** $n\le1200$, `a2_align.gp`.
* characteristic roots $4\cdot(11\pm2i)=44\pm8i$; $|\lambda_{1,2}|=20\sqrt5=44.72136$;
  $c=\lambda_1\lambda_2=2000=2^4\cdot5^3$;
* singularities $z_\pm^{184}=\tfrac14 z^{\eta}_\pm=(11\pm2i)/500$;
* $A^{184}_n\in\mathbb Z$ (as $A^\eta_n\in\mathbb Z$), so $\kappa_p=0$ at every $p$.

### 1.3 Measured invariants (`a1_structure.gp`)

| invariant | $\eta=(11,5,125)$ | AESZ 184 |
|---|---|---|
| char. roots | $11\pm2i$ | $44\pm8i$ |
| $\lambda_1=\lambda_2$ | $\sqrt{125}=11.18034$ | $20\sqrt5=44.72136$ |
| $c=\lambda_1\lambda_2$ | $125=5^3$ | $2000=2^4\cdot5^3$ |
| archimedean Apéry limit | **none** (oscillates) | **none** (oscillates, exactly $\tfrac12\times$) |
| denominator exponent $k$ (sharp) | $3$ | $3$ |
| $\kappa_5$ | $0$ | $0$ |
| $\max_{n\le400}v_5(A_n)$ | $10$ | $14$ |
| $w=\lim\frac1nv_5(\mathrm{Cas}_n)$ | $3$ | $3$ |
| $\sigma_5=w+2\kappa_5$ | $3$ | $3$ |
| $\sigma_2$ | $0$ | $0$ |

$k=3$ is sharp in both rows: $d_n^2B_n\notin\mathbb Z$, $d_n^3B_n\in\mathbb Z$
**[verified]** $n\le200$.  $\kappa_5=0$: the $A_n$ are integers and $v_5(A_n)=O(\log n)$
(max $14$ resp. $10$ over $n\le400$), so **the rows are integral and carry no $\kappa$
resource**.

**A product-formula anomaly worth recording.**  For AESZ 184
$\log|\lambda_1\lambda_2|=\log2000=7.6009$ but $\sum_p\sigma_p\log p=3\log5=4.8283$: the
$2$-part of $c$ is **not harvestable**.  The reason is the non-Zagier normalisation: from
(1.1) the Casoratian satisfies $\mathrm{Cas}_n=\bigl(P_2(n-2)/n^4\bigr)\mathrm{Cas}_{n-1}$ with
$P_2(n-2)=500(2n-3)(n-1)^2(2n-1)$, so
$$v_2(\mathrm{Cas}_n)=2n+0+2n-4n=0,\qquad v_5(\mathrm{Cas}_n)=3n+\tfrac n2+\tfrac n2-n=3n,$$
the $-4v_p(n!)$ from the leading $n^4$ eating exactly the $4\log2$.  Measured
$v_2(\mathrm{Cas}_n)\in\{-4,\dots,8\}$ bounded and oscillating, $v_5(\mathrm{Cas}_n)/n\to3$
($2.9875$ at $n=400$).  **The "$\Lambda$-ceiling" $\sum_p\sigma_p\log p=\log|\lambda_1\lambda_2|$
of `THEORY_NOTES_03` §1 holds only in Zagier normalisation**; for a rank-4 row it can fail
downward, and here it fails by exactly $4\log2$.

---

## 2. The $5$-adic limits, with the character normalisation nailed down

$\chi_5$ = Legendre symbol mod $5$ (even).  At $p=5$, $\omega^2=\chi_5$, so
$\chi_5\omega^{-2}=\mathbf1$ and $L_5(3,\chi_5\omega^{-2})=L_5(3,\mathbf1)$ — the number
`EULER_CRITERION.md` writes $\zeta_5(3)$.  (Beware: `MUM_SURVEY.md` §3 uses
$\zeta_p(k)=L_p(k,\omega^{1-k})$, which at $p=5$, $k=3$ is $L_5(3,\chi_5)$, a *different*
number; we verified $v_5\bigl(L_5(3,\mathbf1)-L_5(3,\omega^{-2})\bigr)=-1$, i.e. they differ.
The survey's hit label "`L_p(3,chi5om^-2)`" means $L_5(3,\mathbf1)$.  **This notational clash
should be fixed in the two documents.**)

`a3_padic.gp`, $N=1000$, Kubota–Leopoldt to $5^{3050}$:

| row | Cauchy $v_5(\xi_N-\xi_{N-1})$ | $\xi_5$ | agreement |
|---|---|---|---|
| AESZ 184 | $2974$ | $\tfrac14L_5(3,\mathbf1)$ | $v_5=2993$ **[verified]** |
| $\eta$ | $2974$ | $\tfrac12L_5(3,\mathbf1)$ | $v_5=2993$ **[verified]** |

`lindep` returns exactly $1/4$ and $1/2$.  $\mathcal E_5(3)=1-\chi_5(5)5^{-3}=1$, so
$r_5=r_\infty$: $\tfrac14$ and $\tfrac12$ at both places, ratio $2$ — Theorem F and
Conjecture D, trivially, by §1.2.

---

## 3. The cross-determinant

With $\xi^\eta_5=2\xi^{184}_5$ the alignment test is
$D_{n,m}=B^{\eta}_nA^{184}_m-2A^{\eta}_nB^{184}_m$.  `a2_align.gp`, exact:

| $n\backslash m$ | 100 | 200 | 300 | 400 | 600 | 800 | 1000 |
|---|---|---|---|---|---|---|---|
| **100** | $\infty$ | 301 | 299 | 301 | 305 | 300 | 301 |
| **200** | 301 | $\infty$ | 598 | 600 | 604 | 599 | 600 |
| **300** | 300 | 599 | $\infty$ | 901 | 905 | 900 | 901 |
| **400** | 301 | 600 | 900 | $\infty$ | 1204 | 1199 | 1200 |
| **600** | 304 | 603 | 903 | 1203 | $\infty$ | 1796 | 1797 |
| **800** | 301 | 600 | 900 | 1200 | 1798 | $\infty$ | 2400 |
| **1000** | 301 | 600 | 900 | 1200 | 1798 | 2399 | $\infty$ |

($\infty$ = $D_{n,n}=0$ exactly.)  So
$v_5(D_{n,m})=3\min(n,m)+O(1)=\min(\sigma_5^\eta n,\sigma_5^{184}m)+O(1)$, the balanced
sampling being $1:1$ — the rigidity law of `THEORY_NOTES_03` §2 in its degenerate, exact form.
This is the strongest possible alignment and the weakest possible evidence: the two rows are the
same row.

---

## 4. The archimedean side: two complex folds, one constant

Both rows have complex-conjugate characteristic roots, so $B_n/A_n$ oscillates and there is **no
Apéry limit at either row** (`a1_structure.gp`: $B_n/A_n$ at $n=50,\dots,400$ runs
$-0.032,-61.2,0.457,0.310,0.246,0.197,0.141,0.036$ for AESZ 184, exactly half the corresponding
$\eta$ values).

`THEOREM_B_EXACT.md` §5 computes $\eta$'s fold constant modularly (fold $\tau_*=(1+2i)/10$ on
the Fricke circle of level $20$) and gets
$\xi^\eta=\tfrac12L(\chi_5,3)+i\tfrac{\pi}{10}L(\chi_5,2)$ to $10^{-94}$.

**A recurrence-only estimator.**  If the characteristic roots are $\rho_\pm$ with
$|\rho_+|=|\rho_-|$ and $A_n=\alpha s^+_n+\bar\alpha s^-_n$, $B_n=\beta s^+_n+\bar\beta s^-_n$
with $s^\pm_n\sim\rho_\pm^n n^\theta$, then
$$\xi_n:=\frac{B_{n+1}-\rho_-B_n}{A_{n+1}-\rho_-A_n}\;\longrightarrow\;\frac\beta\alpha
  \qquad\text{with error }O(1/n)\ \text{times an oscillating factor }(\rho_-/\rho_+)^n .$$
Richardson extrapolation therefore fails (we checked: it diverges), but a long window average
kills the oscillation.  `a5_foldraw.gp`, average over $n=1000..2400$ at $40$ digits:

| row | window average | target | difference |
|---|---|---|---|
| $\eta$, $\rho_-=11-2i$ | $0.4274114945-0.2218593616\,i$ | $\overline{\xi^\eta}=\tfrac12L(\chi_5,3)-i\tfrac{\pi}{10}L(\chi_5,2)$ | $3.6\cdot10^{-6}$ |
| AESZ 184, $\rho_-=44-8i$ | $0.2137056047-0.1109290989\,i$ | $\tfrac14L(\chi_5,3)-i\tfrac{\pi}{20}L(\chi_5,2)$ | $2.4\cdot10^{-6}$ |

(The estimator lands on the conjugate fold; both folds occur, exchanged by $W_{20}$.)  This is
an **independent, non-modular confirmation of `THEOREM_B_EXACT`'s $\xi^\eta$**, and it fixes

$$\xi^{184}=\tfrac14L(\chi_5,3)+i\,\tfrac{\pi}{20}L(\chi_5,2),$$

which also follows **exactly** from §1.2 (the map $u\mapsto\binom{2n}{n}u$ multiplies both
$s^\pm$ by the same real sequence $\binom{2n}{n}$, so it fixes $\beta/\alpha$ up to the factor
$\tfrac12$ from the $B$-normalisation).  Hence `MUM_SURVEY.md` §5.5's blind $p$-adic prediction
"even part $=\tfrac14L(\chi_5,3)$" is **exactly right**, and its unexplained "drift
$5\cdot10^{-2}$ at $n=1400$" is explained: the sequence does not converge at all.

---

## 5. The design rule: what the pair delivers and what it does not

### 5.1 $\delta=0$, at every ratio

Master formula (`05_two_row.tex` eq. (master)), decayer at $\gamma n$, engine at $\alpha n$,
$r=\alpha/\gamma$, $\eta_{\rm cost}=\kappa_5\log5=0$, $S=k\max(\alpha,\gamma)=3\max(r,1)$,
$G=\min(\sigma^{\rm eng}r,\sigma^{\rm dec})\log5=3\min(r,1)\log5$.  With **no decay**,
$\lambda=\Lambda_{\rm dec}$, so $\gamma\log\tfrac1\lambda=-\log\Lambda_{\rm dec}$ and
$H=F+\gamma\log(\Lambda_{\rm dec}/\lambda)=F$:
$$\delta=1-\frac FH=0\qquad\text{identically}.$$
Numerically (`a6_design.gp`), decayer $=$ AESZ 184 / engine $=\eta$:
$F=2.797,\,2.495,\,\mathbf{2.193},\,2.870,\,3.547,\dots$ at $r=0.5,0.75,1,1.25,1.5$; the reverse
assignment gives the same minimum $F=2.193$ at $r=1$.  The cost–resource deficit $2F=4.386$ per
index.

### 5.2 The budget, and what the Hadamard twist costs

With $\lambda_2=\lambda_1$ the order-$m$ score of `05_two_row.tex` eq. (orderm) reads
$\text{budget}=\log\frac1{|\lambda_2|}-k+\sum_p\sigma_p\log p$:

| row | $\log\frac1{|\lambda_2|}$ | $k$ | $\sum_p\sigma_p\log p$ | budget |
|---|---|---|---|---|
| $\eta$ | $-2.4142$ | 3 | $3\log5=4.8283$ | $\mathbf{-0.5858}$ |
| AESZ 184 | $-3.8005$ | 3 | $3\log5=4.8283$ | $\mathbf{-1.9721}$ |

($-0.5858$ reproduces the $-0.586$ recorded for AZ $(11,5,125)$ in `SPORADIC_SCAN2.md` (sporadic table, budget column), an
independent check of the convention.)  The Hadamard twist costs exactly $\log4=1.3863$: it
multiplies $\Lambda$ by $4$ and adds **no** $p$-adic slope, because $\binom{2n}{n}$ is a
$5$-adic unit up to $O(\log n)$.  **Going to rank $4$ buys nothing here.**  As an *engine*
AESZ 184 is also strictly worse than $\eta$ ($\log\rho_2^{\rm eng}=3.80$ vs $2.41$).

### 5.3 The precise deficit, and the shopping list

The corrected criterion of `ONE_CLASS_TWO_WORLDS.md` §6.1, $w+\kappa_p>0$, is **satisfied and
then some**: $w=3$, $\kappa_5=0$, net $5$-adic resource $(w+\kappa_5)\log5=3\log5=4.8283$ per
index — larger than any other row in the project (Zudilin's Catalan row: $\kappa_2\log2=2.773$;
the $\chi_{-3}$ hypergeometric family: at best $(3-\kappa_3)\log3\approx2.5$).  The obstruction is
entirely archimedean:

> **For the class $L(\chi_5,3)$ at $p=5$ every known realisation is a complex fold.**
> $\eta$ and AESZ 184 have $|\lambda_2|=|\lambda_1|$; there is no decayer, and — unlike the
> $L(2,\chi_{-3})$ situation, where row $\mathbf C$ has $\lambda_2=1$ exactly — there is not
> even a row with a polynomial-size linear form to serve as a cheap engine.

Taking $\eta$ as engine ($k_{\rm eng}=3$, $\sigma_5^{\rm eng}=3$,
$\log\rho_2^{\rm eng}=\tfrac32\log5=2.4142$) and a hypothetical decayer with linear-form rate
$\lambda<1$, denominator rates $\kappa_5$ (5-power), $\nu$ (prime-to-5 in $Q$), $k_d$ (in $P$)
and Casoratian rate $w$, the master formula at the balanced ratio $r=\sigma^{\rm dec}_5/3$ gives
$F<0$ iff
$$\log\tfrac1\lambda\;>\;\max(\sigma^{\rm dec}_5,k_d)+\nu+\kappa_5\log5
   -\tfrac12\log5\cdot\sigma^{\rm dec}_5,$$
i.e., when $k_d\le\sigma^{\rm dec}_5$ and $\sigma^{\rm dec}_5=w+2\kappa_5$,
$$\log\tfrac1\lambda\;>\;\bigl(1-\tfrac12\log5\bigr)\,w+2\,\kappa_5+\nu
  \;=\;0.195281\,w+2\kappa_5+\nu .$$
The coefficient of $\kappa_5$ is exactly $2$ because $2(1-\tfrac12\log5)+\log5=2$.
**Reading**: a $\chi_5$ decayer with $\kappa_5=0$ needs only $\log\frac1\lambda>0.1953\,w$ —
extremely cheap; the cost of $5$-power denominators is exactly $2$ per unit of $\kappa_5$, which
is *worse* than in the $\chi_{-3}$ case, so the target to build is an **integral or
nearly-integral** decayer.

**Shopping list (the concrete next construction).**  The $\chi_{-3}$, weight-2 decayer of
`ONE_CLASS_TWO_WORLDS.md` §2 came from residue sums of rational functions with *double* poles at
$\mathbb Z+\tfrac13,\mathbb Z+\tfrac23$, antisymmetric under $t\mapsto-t-b-1$.  The object
needed here is its $(\chi_5,\text{weight }3)$ analogue: **triple** poles at
$\mathbb Z+\tfrac j5$, $j=1,2,3,4$, with the antisymmetry and the two extra linear conditions
that kill $\zeta(3)$, $\pi^2\log$, $\pi/\sqrt5$-type terms and leave a two-term
$Q\,L(\chi_5,3)-P$.  This is Rivoal–Zudilin territory (the same enlargement §9 of that note
sketches for $\beta(4)$) and is **not attempted here**.

---

## 6. Corrections this note forces on earlier documents

1. `MUM_SURVEY.md` §0.4 and §5.5: "$L(\chi_5,3)$ is a **new period class** not in AvSZ's list"
   is wrong as stated — AESZ 184 carries the class of the project's own row $\eta$, by an exact
   Hadamard identity.  The claim survives only in the weak form "$L(\chi_5,3)$ does not appear in
   AvSZ's tables", which remains true because AESZ 184 has no Apéry limit.
   (`MUM_SURVEY.md` §5.5's own suggestion "testing the cross determinant is the obvious next
   computation" is answered: it is identically zero.)
2. The $\zeta_p$ notation clash of §2 above between `MUM_SURVEY.md` and `EULER_CRITERION.md`.
3. `THEORY_NOTES_03_lattices.md` §1 eq. (2): the product formula
   $\sum_p\sigma_p\log p=\log|\lambda_1\lambda_2|$ needs the hypothesis
   $v_p(\mathrm{Cas}_n)=n\,v_p(c)+O(\log n)$, i.e. Zagier normalisation.  AESZ 184 is an explicit
   counterexample without it, by exactly $4\log2$ (§1.3).
4. `MUM_SURVEY.md` §5.5's "archimedean drift $5\cdot10^{-2}$ at $n=1400$" for AESZ 184: the row
   has no archimedean limit at all; the correct statement is the complex fold value of §4.

## 7. Scripts

`lattice/followups/`: `a1_structure.gp` (Hadamard identity, roots, $k$, $\kappa$, Casoratian
rates), `a2_align.gp` (exact identities, cross-determinant grid), `a3_padic.gp` ($5$-adic limits
and the character normalisation), `a4_fold.gp` (Richardson attempt — *fails*, recorded as a
negative), `a5_foldraw.gp` (window-averaged fold estimator), `a6_design.gp` (master formula and
threshold).
