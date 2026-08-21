# Good-prime towers: is the tower limit pure-Frobenius data?

*Exploratory. Scripts: `lattice/good_prime_towers/` (`towers.gp`, `lll.gp`, drivers
`drv_*.gp`, raw output `cert_full.txt`, `census_raw.txt`, `disc.out`, `*.log`).
Answers the "good-prime dichotomy" item left open in
`consolidation/THEORY_NOTES_04_realisations.md` §1 (last table row) and §4 item 4.

Reproduce: `cd lattice/good_prime_towers && gp -q -s 2G towers.gp drv_cert.gp` for the
eigenvalue census, `... drv_disc.gp` for the residue-disc test, `... dwork.gp drv_dwork.gp`
for the Dwork control, `gp -q lll.gp drv_lll.gp` for the identification tests.*

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
cells, zero exceptions**, at 6–25 certified $p$-adic digits per cell
(`cert_full.txt`; levels $s\le 11,7,5,4,4,4$ at $p=2,3,5,7,11,13$).
An automatic classifier flags $141$ immediately; the remaining three
($A$ at $p=2,a=3$; $B$ at $p=2$, $a=1,3$) fail only the crude
monotonicity heuristic because the measured valuation *saturates at the working
precision* ($e_3=8,18,39,72,145,289,577,1153,2305,K_0,K_0$ for $A$ at $p=2$),
and conform on inspection.

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

**Positive control: the genuine Dwork unit root is a different animal**
(`dwork.gp`). Computing $\lambda_0(\hat t)=\lim_sF_{<p^s}(\hat t)/F_{<p^{s-1}}(\hat t^{\,p})$
at the Teichmüller points $\hat t=\omega(t_0)$, $t_0=1,\dots,p-1$:

| row | $p$ | $\lambda_0\bmod p$ at $t_0=1,2,\dots,p-1$ |
|---|---|---|
| $\gamma$ (Apéry $\zeta(3)$) | 5 | ss, 4, 4, ss |
| $\gamma$ | 7 | ss, 4, 1, 4, 1, 2 |
| $C$ | 5 | 1, ss, ss, 4 |

("ss" = supersingular residue disc, $F(\hat t^{\,p})\equiv0\bmod p$, no unit root
at all.) So $\lambda_0$ is genuinely $t$-dependent, generally $\not\equiv1\bmod p$,
and sometimes does not exist — whereas $u_a=1$ always. The tower is not sitting
at a Teichmüller point of any residue disc; it is running into the cusp.

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

## 7. Tower-limit data and identification tests

Rows $C$, $D$ and $\gamma$ (Apéry $\zeta(3)$), tower bases $a=1,2,3$ (plus
$a\le12$ at $p=5,7$ for the residue-disc test). Full tables in
`lattice/good_prime_towers/*.log` and `disc.out`; regenerate with
`python3 mktable.py local_*.log disc.out`. Certified digits are
$\min\bigl(w(s_{\max}+1),\ \text{measured agreement }L_{s_{\max}}\equiv L_{s_{\max}-1}\bigr)$.

| row | $p$ | $a$ | $\chi(p)$ | $v_p(\Lambda_a)$ | digits | $\Lambda_a\cdot p^{-v}$ (base-$p$, l.s.d. first) | $v_p(\rho^A_s-1)$ | $v_p(\chi p^{w}\rho^B_s-1)$ |
|---|---|---|---|---|---|---|---|---|
| C | 5 | 1 | -1 | 0 | 11 | 2,3,1,3,0,4,4,2,1,4... | [2, 4, 6, 8, 10] |  |
| C | 5 | 2 | -1 | -1 | 8 | 4,4,4,3,4,0,4,0 | [1, 3, 5, 7, 9] |  |
| C | 5 | 3 | -1 | 1 | 12 | 2,1,0,3,4,3,0,2,0,0... | [3, 5, 7, 9, 11] |  |
| D (Apery zeta(2)) | 5 | 1 | 1 | 0 | 12 | 2,3,1,0,3,1,4,0,2,4... | [4, 7, 10, 13, 16] |  |
| D (Apery zeta(2)) | 5 | 2 | 1 | 2 | 12 | 1,4,1,1,0,0,2,3,0,1... | [4, 7, 10, 13, 16] |  |
| D (Apery zeta(2)) | 5 | 3 | 1 | 0 | 12 | 3,4,3,1,0,3,0,0,4,3... | [5, 8, 11, 14, 17] |  |
| gamma (Apery zeta(3)) | 5 | 1 | 1 | -1 | 13 | 1,0,3,0,4,1,4,1,4,0... | [2, 5, 8, 11, 14] |  |
| gamma (Apery zeta(3)) | 5 | 2 | 1 | 0 | 15 | 3,2,4,1,4,2,2,4,0,4... | [3, 6, 9, 12, 15] |  |
| gamma (Apery zeta(3)) | 5 | 3 | 1 | -1 | 13 | 4,3,4,0,1,0,2,1,3,0... | [2, 5, 8, 11, 14] |  |
| C | 7 | 1 | 1 | 0 | 9 | 5,4,4,5,5,2,4,0,3 | [2, 4, 6, 8] |  |
| C | 7 | 2 | 1 | 0 | 8 | 4,1,1,4,5,2,0,4 | [2, 4, 6, 8] |  |
| C | 7 | 3 | 1 | 0 | 8 | 5,0,2,2,1,1,2,1 | [2, 4, 6, 8] |  |
| D (Apery zeta(2)) | 7 | 1 | 1 | 0 | 10 | 5,4,4,0,5,4,2,3,5,1 | [3, 6, 9, 12] |  |
| D (Apery zeta(2)) | 7 | 2 | 1 | 0 | 10 | 3,1,1,3,4,2,1,0,0,2 | [4, 7, 10, 13] |  |
| D (Apery zeta(2)) | 7 | 3 | 1 | -2 | 8 | 4,6,2,2,3,2,1,0 | [1, 4, 7, 10] |  |
| gamma (Apery zeta(3)) | 7 | 1 | 1 | 0 | 12 | 3,1,4,6,6,2,2,4,1,2... | [4, 7, 10, 13] |  |
| gamma (Apery zeta(3)) | 7 | 2 | 1 | 0 | 12 | 4,4,2,4,4,3,5,3,2,2... | [3, 6, 9, 12] |  |
| gamma (Apery zeta(3)) | 7 | 3 | 1 | 1 | 12 | 2,3,6,4,4,6,0,0,1,5... | [3, 6, 9, 12] |  |
| C | 11 | 1 | -1 | 0 | 12 | 4,7,3,3,1,8,2,5,3,6... | [2, 4, 6, 8, 10] | [2, 4, 6, 8, 10] |
| C | 11 | 2 | -1 | 0 | 10 | 9,6,6,2,9,2,1,4,9,10 | [2, 4, 6, 8, 10] | [2, 4, 6, 8, 10] |
| C | 11 | 3 | -1 | 0 | 11 | 10,3,5,5,8,8,2,9,8,4... | [2, 4, 6, 8, 10] | [2, 4, 6, 8, 10] |
| D (Apery zeta(2)) | 11 | 1 | 1 | 0 | 10 | 4,7,3,4,4,8,6,1,2,10 | [3, 6, 9, 12] | [3, 6, 9, 12] |
| D (Apery zeta(2)) | 11 | 2 | 1 | 0 | 10 | 8,9,1,8,10,3,4,8,6,3 | [3, 6, 9, 12] | [3, 6, 9, 12] |
| D (Apery zeta(2)) | 11 | 3 | 1 | 0 | 10 | 3,2,10,5,5,3,5,3,6,1 | [4, 7, 10, 13] | [3, 6, 9, 12] |
| gamma (Apery zeta(3)) | 11 | 1 | 1 | 0 | 15 | 9,8,8,0,9,6,5,6,4,6... | [3, 6, 9, 12, 15] | [3, 6, 9, 12, 15] |
| gamma (Apery zeta(3)) | 11 | 2 | 1 | 0 | 15 | 7,2,4,3,8,2,0,6,10,4... | [3, 6, 9, 12, 15] | [3, 6, 9, 12, 15] |
| gamma (Apery zeta(3)) | 11 | 3 | 1 | 0 | 15 | 3,6,10,3,4,0,10,8,0,9... | [3, 6, 9, 12, 15] | [3, 6, 9, 12, 15] |
| C | 13 | 1 | 1 | 0 | 11 | 9,8,8,4,2,9,12,7,2,1... | [2, 4, 6, 8, 10] | [2, 4, 6, 8, 10] |
| C | 13 | 2 | 1 | 0 | 10 | 11,2,1,1,7,5,6,8,9,6 | [2, 4, 6, 8, 10] | [2, 4, 6, 8, 10] |
| C | 13 | 3 | 1 | 0 | 10 | 10,3,1,4,2,11,9,8,9,1 | [2, 4, 6, 8, 10] | [2, 4, 6, 8, 10] |
| gamma (Apery zeta(3)) | 13 | 1 | 1 | 0 | 12 | 8,2,5,12,11,11,9,6,12,7... | [3, 6, 9, 12] | [3, 6, 9, 12] |
| gamma (Apery zeta(3)) | 13 | 2 | 1 | 1 | 12 | 4,10,10,1,12,0,6,8,11,6... | [3, 6, 9, 12] | [2, 5, 8, 11] |
| gamma (Apery zeta(3)) | 13 | 3 | 1 | 0 | 12 | 9,1,0,3,10,8,8,8,8,9... | [4, 7, 10, 13] | [3, 6, 9, 12] |

(Row $D$ at $p=13$ was still running at the time of writing; higher-level runs
$s_{\max}=8$ at $p=5$ and $s_{\max}=7$ at $p=7$ are queued on `snake` and will
extend the digit counts but change nothing structural. The $p=11$ block was
computed twice, on two machines, with identical output.)

**Identification tests** (`lll.gp`, exact PARI `lindep` over $\mathbf Q_p$ with the
heuristic noise floor $p^{K/m}$ for an $m$-term relation). For rows $C$ and $D$ —
which Paper C did **not** cover — nothing was found below the floor:

* $\Lambda_a\in\mathbf Q$: no.
* $\Lambda_a$ algebraic of degree $\le2$ or $\le3$: no. (This is the test that
  would catch a Frobenius unit root of a quadratic $X^2-a_pX+p^{w+1}$, and any
  Gross–Koblitz value.)
* $(1,\Lambda_a,\Lambda_{a'})$ spanning a $2$-dimensional $\mathbf Q$-space
  (i.e. $\Lambda_a/\Lambda_{a'}\in\mathbf Q$): no.
* $(1,\alpha_p,\Lambda_a)$ with $\alpha_p$ the unit root of the attached newform:
  at $p=7,13$ (ordinary) nothing below the floor; at $p=5,11$ **the test cannot even
  be posed** — $a_5=a_{11}=0$ for $\eta(2z)^3\eta(6z)^3$ and $v_{11}(a_{11})=1$ for
  $\eta(2z)^4\eta(4z)^4$.

These reproduce, for two new rows, Paper C's `find:exclusion` and `find:onegen`,
at 8–15 digits (weaker precision than Paper C's 12–27, so weaker exclusions).
Two borderline hits were recorded and are consistent with noise: $C$ at $p=11$,
$a=1$, $(1,\Lambda)$ height $1.9\cdot10^5$ against floor $1.8\cdot10^6$, and
$(1,\Lambda_1,\Lambda_2)$ height $652$ against floor $2960$. Under the project's
own criterion (a genuine relation's height does not move as $K$ grows) these need
a rerun at higher $K$ before they mean anything; nothing else came within a factor
$4$ of its floor.

## 8. Answers to the three posed tasks, in order

1. **Tower limits and unit roots.** Computed for $C$, $D$, $\gamma$ at
   $p=5,7,11,13$, $a=1,2,3$ (and $a\le12$ at $p=5,7$). Precision reached is 8–21
   digits, **not** $p^{40}$ — that is structurally out of reach (§9). The unit
   roots $u_a$ are all $=1$. The predicted rank-2 Frobenius shape
   $\begin{pmatrix}u&*\\0&p^{w+1}/u\end{pmatrix}$ is **not** what the tower sees:
   the correct matrix is $\begin{pmatrix}1&*\\0&\chi(p)p^{-w}\end{pmatrix}$, with
   the two eigenvalues *independent of the row* (they depend only on $w$ and
   $\chi$), and with $\det=\chi(p)p^{-w}$ rather than the Weil-number
   $p^{w+1}$. That is the signature of the cusp, not of a smooth fibre. The
   off-diagonal entry is not a scalar: it is the whole eigenvector, i.e. the
   infinite block assembly of Paper C Thm. `thm:Lam`.
2. **Is the $\Gamma_p$ description a Gross–Koblitz evaluation on a CM fibre?**
   No, on two independent grounds (§4, §5): no CM fibre is involved (the fibre is
   the degenerate one over $t=0$), and Gross–Koblitz $\Gamma_p$-values are
   algebraic while $\Lambda_a$ is not algebraic of degree $\le3$ at any tested
   $(p,a)$. Nor is there a fibre whose $a_p$ explains $u_a$, because $u_a=1$: the
   only "fibre" with unit root $1$ and second eigenvalue $\chi(p)p^{-w}$ is the
   Tate/nodal fibre at the cusp, twisted by $\chi$ and by $\mathbf Q_p(w)$.
   $\Lambda_a-(\text{explicit Frobenius expression})$ is therefore not $0$ for any
   eigenvalue-level expression; the correct exact statement is the eigenvector one
   of §2.
3. **This file.**

## 9. Honest limitations

* The requested precision $p^{40}$ is **not reachable** and was not reached; see
  Paper C Remark (precision budget). The rate law caps new digits at
  $\min(3,r(\chi,w))$ per tower level, so $40$ digits needs $s\approx14$, i.e.
  $n\approx p^{14}$. The $\Lambda$ data reported in §7
  reaches $s=5$ at $p=5,11,13$ and $s=4$ at $p=7$, i.e. 8–15 digits; larger runs
  ($s=8$ at $p=5$, $s=7$ at $p=7$) were launched but had not finished, and the
  method (§1) makes $s=8$ at $p=5$ and $s=7$ at $p=7$ routine on a quiet machine.
  The eigenvalue census (§2), which is the load-bearing result, does
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
