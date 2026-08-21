# Good-prime towers: is the tower limit pure-Frobenius data?

*Exploratory. Scripts: `lattice/good_prime_towers/` (`towers.gp`, `lll.gp`, drivers
`drv_*.gp`, raw output `cert_full.txt`, `census_raw.txt`, `disc.out`, `*.log`).
Answers the "good-prime dichotomy" item left open in
`consolidation/THEORY_NOTES_04_realisations.md` §4.4.*

---

## 0. Verdict, first

The hypothesis under test was: *at a good prime $p\nmid N$ the tower limit
$\Lambda_a$ is a period of the **pure** local system, read off from the unit-root
eigenvector of the Frobenius matrix of the rank-2 (resp. rank-3) Picard–Fuchs
module; in particular $a_{ap^{s+1}}/a_{ap^s}\to$ the Dwork unit root, and the
$\Gamma_p$-values appear because Frobenius on CM fibres is a Gauss sum
(Gross–Koblitz).*

**Split verdict.**

| claim | verdict |
|---|---|
| (H1) The tower limit is Frobenius data of the pure system, not a $p$-adic $L$-value | **CONFIRMED**, and sharpened to an exact eigenvalue statement (§2) |
| (H2) $u_a=\lim_s a_{ap^{s+1}}/a_{ap^s}$ is the Dwork unit root of the fibre | **REFUTED.** $u_a=1$ *identically*, at every row, prime and tower base tested (§3) |
| (H3) $\Lambda_a$ is a rational function of $u_a$, $p$ and an off-diagonal entry | **VACUOUS as stated** ($u_a=1$); the correct statement is (§2)/(§6): $\Lambda$ is the *ratio of the two Frobenius eigenvectors*, and the eigenvectors — not the eigenvalues — carry all the arithmetic |
| (H4) $\Lambda_a$ relates to $a_p$ of an elliptic/CM fibre at a Teichmüller point | **REFUTED.** The relevant fibre is the *degenerate* fibre over the cusp $t=0$; no $a_p$ enters. At $p=11$ the motive attached to the Apéry $\zeta(3)$ row is **non-ordinary** ($a_{11}=-44$, $v_{11}=1$), so no unit root exists in $\mathbf Q_{11}$ — yet the tower is perfectly ordinary with eigenvalues $(1,11^{-3})$ (§4) |
| (H5) $\Gamma_p$ appears via Gross–Koblitz / Gauss sums | **REFUTED** (§5): Gross–Koblitz concerns $\Gamma_p$ at rationals of small denominator, which is *algebraic* and is separately excluded by Paper C's algebraicity tests; the $\Gamma_p$ here is at $p$-power arguments (Morita/Kazandzidis) |
| (H6) $\Lambda_a$ is determined by the residue disc, i.e. depends only on $a\bmod p$ | **REFUTED** (§4) |
| (H7) The dichotomy "pure Frobenius at good primes, extension regulator at slope primes" | **CONFIRMED**, and upgraded to a **trichotomy** with an exact criterion (§6) |

So: Fable's *headline* is right and now has an exact, cheaply-checkable form.
Fable's *mechanism* is wrong in every particular — there is no unit root, no
$a_p$, no Gauss sum, no residue disc. The tower sees Frobenius **at the cusp
$t=0$**, where the fibre is the degenerate (Tate/nodal) one, and the eigenvalues
are the trivial pair $(1,\chi(p)p^{-w})$ forced by that degeneration. All the
arithmetic sits in the eigen*vectors*, which are exactly the $\Gamma_p$-tower
assemblies of Paper C / the $\Lambda$-algebra paper.

---

## 1. Set-up and method

Rows in Zagier/Malik–Straub normalisation, $A(0)=1,A(1)=b$, $B(0)=0,B(1)=1$:

* (R2) $(n+1)^2u_{n+1}=(an^2+an+b)u_n-cn^2u_{n-1}$,
* (R3) $(n+1)^3u_{n+1}=(2n+1)(an^2+an+b)u_n-cn^3u_{n-1}$;

$w$ is the minimal exponent with $d_n^wB(n)\in\mathbf Z$ ($w=2$ for the six (R2)
rows, $w=3$ for the (R3) rows) and $\chi$ the quadratic character of
Table (families) in `padic-apery-limits.tex`.

$$\Lambda_a=\lim_{s\to\infty}\chi(p)^sp^{ws}\frac{B(ap^s)}{A(ap^s)}.$$

**Computational device (new, and the reason the runs are cheap).** Put
$\widehat u_n=(n!)^{w}u_n$. Then the recurrences become *division-free integer*
recurrences

$$\widehat u_{n+1}=(an^2+an+b)\widehat u_n-c\,n^4\widehat u_{n-1}\quad(\mathrm{R2}),
\qquad
\widehat u_{n+1}=(2n+1)(an^2+an+b)\widehat u_n-c\,n^6\widehat u_{n-1}\quad(\mathrm{R3}),$$

and $B(n)/A(n)=\widehat B_n/\widehat A_n$ exactly, so the tower ratios need no
factorials at all. Everything is done in $\mathbf Z/p^{K_0}$ with
$K_0=w\,v_p(N!)+K+40$; there is **no precision loss anywhere** (no division is
ever performed) and each step costs one big$\times$small multiply. Cross-level
ratios ($A(ap^{s+1})/A(ap^s)$ and the same for $B$) do need the factorial
correction, and the $p$-unit part of $n!$ is carried alongside at one extra
big$\times$small multiply per step. This is an $O(N\cdot K_0)$ pass; $N=10^6$
runs in minutes.

---

## 2. The headline: the two tower Frobenius eigenvalues, exactly

Write $\rho^A_s=A(ap^{s+1})/A(ap^s)$ and $\rho^B_s=B(ap^{s+1})/B(ap^s)$. These are
the two eigenvalue-ratios of the operator $(U_pX)(n)=X(pn)$ along the tower.

**Measurement (`cert(...)` in `towers.gp`, output `cert_full.txt`).**
For each cell the script measures $v_p(\rho^A_s-1)$ and, for the $B$-row, the
three competing hypotheses $v_p(\rho^B_s-1)$, $v_p(\chi(p)p^{w}\rho^B_s-1)$,
$v_p(\rho^B_s)$, across all available levels $s$.

**Result.** In **every** cell:

$$\boxed{\ \lim_s \rho^A_s\ =\ 1\ \text{ exactly}\ }$$

and for the $B$-row exactly one of three alternatives holds, decided by whether
$p\mid c$ and by the placement of the character in the Eisenstein source:

| case | $\lim_s\rho^B_s$ | $\xi_p=\lim_{n\to\infty}b_n/a_n$ |
|---|---|---|
| $p\nmid c$ (**good prime**) | $\chi(p)\,p^{-w}$ — *ordinary, distinct eigenvalues* | does not exist; tower limits $\Lambda_a$ exist |
| $p\mid c$ (**slope prime**), inner placement $\varphi=\mathbf 1$ | $1$ — *unipotent, repeated eigenvalue* | $=-Q(w+1)\kappa_p$, a Kubota–Leopoldt value (Thm F) |
| $p\mid c$, outer placement $\varphi\neq\mathbf 1$ or cuspidal source | $0$ — *degenerate* | $=0$ |

The convergence is linear in $s$ with slope $w$ per level, e.g. for row $C$ at
$p=5$: $v_5(\chi(5)5^{2}\rho^B_s-1)=2,4,6,8,10$ for $s=0,\dots,4$; for the Apéry
$\zeta(3)$ row at $p=5$: $3,6,9,12,15$.

The unit in the second eigenvalue is **exactly the tabulated quadratic
character**, verified digit-by-digit: row $C$ ($\chi=\chi_{-3}$) gives
$\lim\rho^B\cdot p^{2}=-1$ at $p=2,5,11$ and $+1$ at $p=7,13$; row $E$
($\chi=\chi_{-4}$) gives $-1$ at $p=3,7$ and $+1$ at $p=5$; row $\eta$
($\chi=\chi_5$) gives $-1$ at $p=2,3,7,13$ and $+1$ at $p=11$.

**Scope of the census.** 12 rows ($A,B,C,D,E,F,\gamma,\alpha,\delta,\varepsilon,
\zeta,\eta$) $\times$ 6 primes ($2,3,5,7,11,13$) $\times$ 2 tower bases = **144
cells, zero exceptions**, at 6–25 certified $p$-adic digits per cell.

Note this census includes $p=2,3$, where the descent law of Paper A is stated
only for $p\ge5$: the eigenvalue dichotomy holds there too, empirically.

**Consequences that follow formally from the eigenvalue statement.**
Writing $\widetilde A_a=\lim_sA(ap^s)$ and $\widetilde B_a=\lim_s\chi(p)^sp^{ws}B(ap^s)$,

$$\widetilde A_{pa}=\widetilde A_a,\qquad
  \widetilde B_{pa}=\chi(p)p^{-w}\widetilde B_a,\qquad
  \Lambda_a=\widetilde B_a/\widetilde A_a ,$$

i.e. $\widetilde A$ and $\widetilde B$ are $U_p$-eigenvectors with eigenvalues
$1$ and $\chi(p)p^{-w}$; the scaling law $\Lambda_{pa}=\chi(p)p^{-w}\Lambda_a$
(Paper C Prop. 2.5(iii)) is exactly the statement that $\Lambda$ is the *ratio of
two Frobenius eigenvectors*. This is the correct form of the hypothesis: the
tower limit is the coordinate of the canonical Frobenius splitting, and the
splitting exists **because the two eigenvalues are distinct**, which is precisely
$p\nmid c$.

---

## 3. $u_a=1$: the "Dwork unit root" carries no information

The specific prediction $a_{ap^{s+1}}/a_{ap^s}\to$ (unit root of Frobenius on the
fibre) is false. Measured $v_p(\rho^A_s-1)$:

| row | $p$ | $a$ | $v_p(\rho^A_s-1)$, $s=0,1,2,\dots$ |
|---|---|---|---|
| $\gamma$ (Apéry $\zeta(3)$) | 5 | 1 | 2, 5, 8, 11, 14 |
| $\gamma$ | 5 | 2 | 3, 6, 9, 12, 15 |
| $C$ | 5 | 1 | 2, 4, 6, 8, 10 |
| $D$ | 5 | 1 | 4, 7, 10, 13, 16 |
| $D$ | 7 | 1 | 3, 6, 9, 12 |

The valuation grows linearly without bound, so $u_a=1$ on the nose. This is not
an accident of these rows: it is the Beukers–Coster supercongruence
$A(mp^{r})\equiv A(mp^{r-1})\pmod{p^{wr}}$, which is exactly the statement that
the unit-root eigenvalue *along the tower* is $1$. A genuine Dwork unit root
$\lambda_0(t)=F(t)/F(t^\sigma)$ at a Teichmüller point $t_0\ne0$ is a unit
congruent to $a_p$ mod $p$ and is generally $\not\equiv1$; the tower forces
$t\to0$, and at the cusp the fibre is the Tate/nodal curve, whose Frobenius unit
root **is** $1$.

The second eigenvalue $\chi(p)p^{-w}$ is likewise the cusp value: $p^{-w}$ is the
Tate twist (the $d_n^w$ denominators of the second solution), and $\chi(p)$ is the
quadratic twist of the modular parametrisation. Both are "trivial" Frobenius
data; neither sees a fibre.

---

## 4. Two clean falsifications

**(a) A non-ordinary prime where the tower is ordinary.** Beukers attaches to the
Apéry $\zeta(3)$ numbers the newform $\eta(2z)^4\eta(4z)^4\in S_4(\Gamma_0(8))$;
its Hecke eigenvalues are $a_5=-2$, $a_7=24$, $a_{11}=-44$, $a_{13}=22$. At
$p=11$, $v_{11}(a_{11})=1$: the motive is **non-ordinary**, $X^2-a_{11}X+11^3$ has
both roots of positive slope, and there is *no unit root in $\mathbf Q_{11}$*.
Nevertheless $\Lambda_1,\Lambda_2,\Lambda_3$ exist at $p=11$ and the tower
Frobenius pair is the ordinary $(1,11^{-3})$. A tower limit governed by the
unit-root eigenvector of the motive's Frobenius could not exist at $p=11$. The
same happens for the $\chi_{-3}$-rows against the CM form
$\eta(2z)^3\eta(6z)^3\in S_3(\Gamma_0(12),\chi_{-3})$, whose $a_p$ vanishes at
the inert primes $p=5,11$.

**(b) $\Lambda_a$ is not a function of the residue disc.** If $\Lambda_a$ were
read off at the Teichmüller point of the $a$-th residue disc it would depend only
on $a\bmod p$. It does not — even the valuation differs. At $p=5$
(`disc.out`):

| row | $v_5(\Lambda_1)$ | $v_5(\Lambda_6)$ | $v_5(\Lambda_{11})$ | $v_5(\Lambda_2)$ | $v_5(\Lambda_7)$ | $v_5(\Lambda_{12})$ |
|---|---|---|---|---|---|---|
| $\gamma$ | $-1$ | $-4$ | $-2$ | $0$ | $-4$ | $-3$ |
| $C$ | $0$ | $-2$ | $-3$ | $-1$ | $-2$ | $-2$ |
| $D$ | $0$ | $-2$ | $0$ | $2$ | $-2$ | $0$ |

The $a$-dependence is instead the one Paper C proves: $\Lambda_a$ is built from
$h_p(a)=\sum_{x\in\mathbf Z[1/p],\,0<x\le a}x^{-w}$ and Kazandzidis limits
$c_p(ap^i,j)$, which are genuinely functions of the integer $a$, not of $a\bmod p$.

---

## 5. On Gross–Koblitz

The $\Gamma_p$ that occurs is Morita's $\Gamma_p$ at **$p$-power arguments**,
through the Kazandzidis limits $c_p(A,B)=\lim_t\binom{Ap^t}{Bp^t}$ and the tower
Betas $\beta_p(B,D)=\gamma_p(B+D)/(\gamma_p(B)\gamma_p(D))$, whose logarithms are
infinite $\mathbf Q$-combinations of *all* the odd $\zeta_p(m)$. Gross–Koblitz
evaluates $\Gamma_p$ at rationals of small denominator as Gauss sums; those
values are **algebraic**, and both papers already exclude algebraicity of degree
$\le3$ for $\Lambda_a$ at $p=5,7,11,13$. So the two occurrences of $\Gamma_p$ are
disjoint, and the CM/Gauss-sum explanation of the $\Gamma_p$-values is not
available. (Independently: §2 shows no CM fibre is involved.)

---

## 6. The single dichotomy statement

Let $(\Gamma,t,F)$ be a modular Apéry row of weight $w$ with Eisenstein source
$\Phi=P(V)E^{\psi,\varphi}_{w+2}$, and let $\Phi_p$ be Frobenius on the
Picard–Fuchs crystal in the formal neighbourhood of the cusp $t=0$, read through
$U_p$ on the coefficient sequences. Then, in the normalisation above:

> **$U_p$ has eigenvalue $1$ on the $A$-row always. On the $B$-row it has
> eigenvalue $\chi(p)p^{-w}$ iff $p\nmid c$; eigenvalue $1$ iff $p\mid c$ and
> $\varphi=\mathbf 1$; eigenvalue $0$ iff $p\mid c$ and the source is
> outer-placed or cuspidal.**
>
> * **Distinct eigenvalues** ($p\nmid c$, good prime): the extension splits
>   canonically over $\mathbf Q_p$, $b_n/a_n$ has no limit, and the invariant is
>   the splitting coordinate $\Lambda_a=\widetilde B_a/\widetilde A_a$ — pure
>   Frobenius-eigenvector data, living in the $\Gamma_p$-tower ($\Lambda$-)algebra,
>   with no relation to any $p$-adic $L$-value and no cross-row rigidity.
> * **Repeated eigenvalue $1$** ($p\mid c$, inner): the crystal is unipotent, the
>   splitting is obstructed, and the obstruction *is* the invariant: $b_n/a_n$
>   converges over all $n$ to $\xi_p=-Q(w+1)\kappa_p$, the $p$-adic (syntomic)
>   regulator of the extension class, a Kubota–Leopoldt value. Conjecture D
>   (rigidity) holds here and only here.
> * **Eigenvalue $0$** ($p\mid c$, outer/cuspidal): the extension is split by the
>   vanishing of the constant term and $\xi_p=0$.

*"$p$-adic Apéry data = pure Frobenius at good primes, extension regulator at
slope primes"* is therefore right, with the sharpening that the "pure Frobenius"
in question is the **cusp** Frobenius $(1,\chi(p)p^{-w})$, which is the same for
every row of the same $(w,\chi)$; the rows differ only in the eigenvectors.
That is also the correct explanation of Paper C's Finding `cross` (no relations
between rows sharing a period at good primes): the eigen*values* — the only
things a period ratio of the *eigenvalue* kind could see — are identical across
those rows, so any relation would have to come from the eigenvectors, and those
are unrelated infinite $\Gamma_p$-assemblies.

It also explains, in one line, the negative half of `EULER_CRITERION.md`: when
the Euler factor fails to divide $P$, the geometric series
$\sum_k\psi(p)^kp^{-k(w+1)}V_p^k$ diverges — but *its rate of divergence is
exactly the second Frobenius eigenvalue*, which is why the tower increments
decrease **linearly with slope $w$** ($-2,-4,-6,\dots$ for $C$ at $p=2$;
$-3,-6,-9,\dots$ for Domb at $p=3$) rather than erratically. Dividing that
divergence out is precisely the $\chi(p)^sp^{ws}$ normalisation, and what is left
is $\Lambda_a$.

---

## 7. Honest limitations

* The requested precision $p^{40}$ is **not reachable** and was not reached; see
  Paper C Remark (precision budget). The rate law caps new digits at
  $\min(3,r(\chi,w))$ per tower level, so $40$ digits needs $s\approx14$, i.e.
  $n\approx p^{14}$. The runs here reach $s=8$ at $p=5$, $s=7$ at $p=7$, $s=5$–$6$
  at $p=11,13$; the eigenvalue census (§2), which is the load-bearing result, does
  not need precision — it needs only that the measured valuations increase, and
  they increase by $w$ per level in every cell.
* The eigenvalue statements are **verified**, not proved. For $p\ge5$ the good-prime
  case is a restatement of the master descent law $(\mathrm{LB}^\chi_w)$ of Paper A
  (itself a Finding for all fifteen rows, proved for four), and the $A$-row case
  is Beukers–Coster. The slope-prime cases $(1,1)$ and $(1,0)$ are, as far as this
  file goes, new empirical statements; they are consistent with Theorem F but are
  not derived from it here.
* $p=2,3$ are included in the census but lie outside the range of every theorem
  quoted.
* Two cells of the first census pass ($A$ at $p=2$, $\zeta$ at $p=3$) initially
  reported garbage because $B(ap^s)$ underflows to $0$ in the working precision;
  they are exactly the eigenvalue-$0$ cells, and were re-run separately.
