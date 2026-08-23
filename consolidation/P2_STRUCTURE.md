# The structure of the Catalan congruence lattices $\mathcal K_n$, and what P2$'$ really says

*Fable, 2026-08-23.  Scripts: `lattice/p2_structure/`, data in `lattice/p2_structure/data/`.
Tags: **[proved]** = mathematical proof given here or cited; **[verified]** = exact
integer/rational computation over a stated range; **[measured]** = numerical fit or statistic;
**[open]**.  All lattice arithmetic is exact in $\mathbb Z^2$; only the two metric weights are
real ($G$ at $3\,000$ digits, `\p 3000`).*

**No claim of irrationality is made anywhere below.**  Every column that depends on $G$ only
through the ratio $(X_nG-Y_n)/(V_nG-U_n)$ is reproduced verbatim by a rational surrogate up to
an explicit horizon, quantified exactly in §5.3.

---

## 0. Verdict

1. **$\mathcal K_n$ is a continued-fraction lattice, and its first minimum is a convergent.**
   In oriented coordinates $\mathcal K_n=\mathbb Z(h_{11},0)+\mathbb Z(h_{12},h_{22})$, and the
   first minimum $v_1(n)$ in the weighted metric is, for **every** one of the 351 instances
   $4\le n\le120$, $k\in\{22.4,23.0,23.9\}$, a *convergent* $p_i/q_i$ of $h_{12}/h_{11}$ — never
   an intermediate fraction.  This is a theorem (relative minima of a monotone norm are best
   approximations of the second kind), and the balance index sits at
   $i(n)/L(n)=0.5195\pm0.0165$ of the continued fraction, against the predicted
   $\tfrac12+\tfrac{E_2-E_1}{2\kappa}=0.524$.  **[proved + verified]**
2. **The cone condition is a parity condition on that index, and nothing else.**
   The errors $h_{12}q_i-p_ih_{11}$ alternate in sign with $i$, so
   $$\boxed{\ v_1(n)\in\pm\mathcal P\iff i(n)\ \text{is even}\ }$$
   — verified with *no exception* in all 351 instances (joint counts:
   $(\text{even},\text{in-cone})=186$, $(\text{odd},\text{outside})=165$, the two mixed cells
   empty).  The
   $\approx50\%$ in-cone frequency of `06_threshold.tex` and `POSITIVITY_PROGRAM.md` is
   therefore not a coincidence to be explained but the parity of a continued-fraction index,
   and the $460\,365$ partial quotients of the $h_{12}/h_{11}$ obey Gauss–Kuzmin to three
   decimals ($\chi^2=5.09$ on $8$ df).  **[proved + measured]**
3. **An exact dichotomy replaces the ratio $\rho$.**  With $(b_1,b_2)$ Gauss-reduced in the
   weighted metric and $\rho=(\text{cone-min})/\lambda_1$:
   $$\pm b_1\in\mathcal P\ \Longrightarrow\ 1\le\rho\le\sqrt2;\qquad
     \pm b_1\notin\mathcal P\ \Longrightarrow\ \tfrac{\sqrt3}{2}\tfrac{\lambda_2}{\lambda_1}\le\rho
       \le\sqrt2\Bigl(1+\tfrac{\lambda_2}{\lambda_1}\Bigr).$$
   **[proved]**, no violations in 351 instances.  The median $\rho=\sqrt2$ of
   `POSITIVITY_PROGRAM.md` §3.3 is thus *exactly* the statement that the parity is even about
   half the time; $\sqrt2$ is the boundary of the dichotomy, not a fitted constant.
4. **P2$'$ is equivalent to a statement about one partial quotient.**  Modulo the (proved)
   dichotomy, $$\text{P2}'\iff\text{for infinitely many }n:\ i(n)\ \text{is even, or }
   a_{i(n)+1}(n)=e^{o(n)} .$$  ($\iff$ is proved for the first, dichotomy, form; the
   partial-quotient form uses the measured $\lambda_2/\lambda_1\asymp a_{i+1}$.)  Since $i(n)$ is
   even for $53\%$ of $n\le120$ and the longest observed run of odd $n$ is $5$, P2$'$ holds in the
   computed range along a set of density $\approx\tfrac12$, not just sporadically.
   **[proved + measured]**
5. **The lattices are statistically indistinguishable from Haar-random lattices.**  Against
   $2\times10^5$ Haar-random unimodular planar lattices with a uniformly random cone: in-cone
   $0.4993$ vs $0.5299$; median $\rho$ $1.4112$ vs $1.4023$; $\Pr[\rho>20]$ $0.0270$ vs
   $0.0399$; median $\lambda_2/\lambda_1$ $1.9336$ vs $2.0221$; Kolmogorov–Smirnov distance
   $0.0593$ against a $5\%$ critical value $0.0725$ — **not rejected**.  **[measured]**
6. **Two corrections to the record.**
   (a) $51$ of the $231$ cone minima of `lattice/positivity/data/cone_n80.csv` are **not
   minimal**; the exact scan here returns strictly smaller values in every one of those cases
   (e.g. $k=22.4$, $n=12$: $-0.6793\to-0.7024$, $\rho\ 3.24\to2.46$).  All corrections lower
   $\rho$, i.e. they strengthen P2.  **[verified]**
   (b) The index deficit is identified exactly: $M_n/[\mathbb Z^2:\mathcal K_n]=\gcd(X_n,Y_n,V_n,U_n)$
   for all 351 instances, and $\log\gcd\in[9.0,29.3]$ over $4\le n\le120$ with no linear growth
   ($v_2$ of it grows like $\approx4.3\log n$).  So the gcd-loss constant $c$ of
   `POSITIVITY_PROGRAM.md` §3.4 is $c=0$: $\kappa_n=\sigma(k)-O(\log n/n)$ and the achievable
   rate is $F(k)$, not $F(k)-c/2$.  This settles ranked next step #1 of that document.
   **[verified $n\le120$]**
7. **The kernel-vector proposition of `06_threshold.tex` needs amending, and the honest version
   is stronger.**  The kernel direction is *always* in $\mathcal K_n$ (no multiple needed) and
   its two scaled coordinates are exactly equal, so its length is
   $\sqrt2\,b\,(M_n/[\mathbb Z^2:\mathcal K_n])\operatorname{covol}$; it becomes the first
   minimum only when $|F|>c/2$, which at the measured deficits **fails at $k=22.4$**.  But
   P2$'$ fails for a rational $G=a/b$ at every $k$ anyway, by the trivial bound
   cone-min $\ge1/b$, at the rate $|F_{\rm eff}|$.  **[proved]**
8. **The rational surrogate cannot be pushed into the failure regime — by a fixed factor.**
   The horizon at which $G^*=\mathrm{bestappr}(G,10^E)$ stops reproducing the $G$-table is
   measured to be $n_f=13,25,50,99$ for $E=40,80,160,320$ (identical for all three $k$),
   i.e. $n_f=2E\log10/(\xi-E_1)$ with $\xi-E_1=14.43$; the $n$ at which its rationality would
   become visible is $n^*=\log b/|F_{\rm eff}|$.  Their ratio is $b$-independent:
   $$\frac{n^*}{n_f}=\frac{\xi-E_1}{2|F_{\rm eff}|}=11.0\ (k=23.9),\quad 21.1\ (23.0),\quad
     53.8\ (22.4).$$
   The $n=98$ of Proposition (Control) is exactly one below the horizon $n_f=99$ of the
   surrogate used there.  **[measured, exact horizon]**

---

## 1. Setup

Rows as in `CATALAN_POSITIVITY.md` §1 and `POSITIVITY_PROGRAM.md` §3.1:
$$X_n=2^{e_{3n}}D_{6n}^2Q_{3n},\quad Y_n=2^{e_{3n}}D_{6n}^2P_{3n},\quad
V_n=4^{7n+1}D_{6n}^2B_n,\quad U_n=4^{7n}D_{6n}^2C_n,$$
$S_n=D_{6n}^2$, $T_n=2^{\lfloor kn\rfloor}$, $M_n=S_nT_n$, and the honest congruence lattice
$$\mathcal K_n=\{c\in\mathbb Z^2:\ M_n\mid c_ZX_n+c_NV_n,\ M_n\mid c_ZY_n+c_NU_n\},$$
computed by `matkerint`/`mathnf` — this is the lattice of `06_threshold.tex` §5.1 ("the
combination $(q,p)=(c_ZX+c_NV,\ c_ZY+c_NU)/M_n$ is integral") and of `POSITIVITY_PROGRAM.md` §3.
The *single*-congruence lattice $\{c:c_ZY_n+c_NU_n\equiv0\bmod S_n\}$, which the same paper quotes
in Theorem (Absorption)(4), is a different object; it is treated in §7, where every structural
statement below is seen to hold for it verbatim.  Weights
$\lambda_Z=|X_nG-Y_n|/M_n$, $\lambda_N=(V_nG-U_n)/M_n$; signs $s_Z=\operatorname{sign}(X_nG-Y_n)=(-1)^{3n}$,
$s_N=+1$.

**Oriented coordinates.**  Put $(u,v)=(s_Zc_Z,\ s_Nc_N)$ (the letters $a,b$ are reserved for a
hypothetical $G=a/b$).  Then the positive cone $\mathcal P$ is the closed first quadrant of
$\mathbb Z^2$, the metric is $\mathrm{diag}(\lambda_Z,\lambda_N)$, and the linear form is
$\ell(u,v)=u\lambda_Z+v\lambda_N=|q_nG-p_n|$ on $\mathcal P$.  All of
§§2–4 are statements about *a sublattice of $\mathbb Z^2$ and the coordinate quadrant*; the
number theory has been fully absorbed into $(h_{11},h_{12},h_{22};\lambda_Z,\lambda_N)$.

Constants: $E_1=13.0995887908$, $E_2=14.3931452672$, $k_*=22.3512905953$, $\sigma(k)=12+k\log2$,
$F(k)=\tfrac{\log2}{2}(k_*-k)$, so $F(22.4)=-0.01688$, $F(23.0)=-0.22483$, $F(23.9)=-0.53674$.
Measured row rates (at $n=60$): $\xi:=\tfrac1n\log X_n=27.521$, $\tfrac1n\log V_n=29.246$,
$\tfrac1n\log|X_nG-Y_n|=13.096$, $\tfrac1n\log(V_nG-U_n)=14.297$.  We write
$F_{\rm eff}(n)=\tfrac1{2n}\log\operatorname{covol}$ for the measured half-log-covolume.

## 2. Exact description of $\mathcal K_n$

### 2.1 Hermite form, index, contents

$\mathcal K_n=\mathbb Z(h_{11},0)+\mathbb Z(h_{12},h_{22})$ with $0\le h_{12}<h_{11}$,
$[\mathbb Z^2:\mathcal K_n]=h_{11}h_{22}$.  Verified for all 351 instances (`verify.gp`):

| statement | status |
|---|---|
| $S_n\mid X_n$ and $S_n\mid V_n$ | 0 failures |
| $[\mathbb Z^2:\mathcal K_n]\mid M_n$ | 0 failures |
| $M_n/[\mathbb Z^2:\mathcal K_n]=\gcd(X_n,Y_n,V_n,U_n)$ | **0 failures** |
| content $\gcd(h_{11},h_{12},h_{22})=1$ | 0 failures |
| directional content $m_Z:=\gcd(h_{11},h_{12})>1$ | 0 failures; mean $\tfrac1n\log m_Z=2.3911$, **independent of $k$** |

Three readings.

* **The gcd loss is the row content.**  $M_n/[\mathbb Z^2:\mathcal K_n]=\gcd(X_n,Y_n,V_n,U_n)=2^{v}\cdot(\text{small odd})$
  with $v\in[13,32]$ and odd part $\le1309$ over the whole range; $\log$ of it lies in
  $[9.01,29.27]$ with mean $20.8$ and **no linear trend** ($v\approx7.3+4.3\log n$).  Hence
  $\tfrac1n\log(M_n/[\mathbb Z^2:\mathcal K_n])=O(\log n/n)\to0$: the constant $c$ left undecided in
  `POSITIVITY_PROGRAM.md` §3.4 (there fitted equally well by $0.196+11.7/n$ and
  $-0.196+7.55\log n/n$) is **$c=0$**, the second fit being the right one.  **[verified $n\le120$]**
* **$m_Z$ is the Absorption(4) directional mass.**  Every $c\in\mathcal K_n$ has $m_Z\mid c_Z$
  with $\tfrac1n\log m_Z\approx2.39$; the congruences force this on the whole lattice, which is
  exactly clause (4) of Theorem (Absorption) in `06_threshold.tex` (`CATALAN_DIRECTIONAL.md` §3):
  no $\ell$-adic device can exploit the divisibility because the kernel vector has it too.
  That $m_Z$ is the same for $k=22.4,23.0,23.9$ shows it is a property of the rows mod $S_n$,
  not of the $2$-adic factor $T_n$.
* Consequently the continued fraction that governs everything is that of the *reduced* fraction
  $(h_{12}/m_Z)/(h_{11}/m_Z)$, of denominator $e^{(27.19-2.39)n}$; its length is measured at
  $L(n)=20.95\,n$ against the Lévy prediction $0.8427\times24.80\,n=20.90\,n$.  **[measured]**

### 2.2 Gauss reduction, $v_1$ and $v_2$

`gred` performs exact Lagrange reduction in the metric $\mathrm{diag}(\lambda_Z,\lambda_N)$; in
dimension $2$ the reduced basis realises both successive minima, so $\lambda_1=|b_1|$,
$\lambda_2=|b_2|$ exactly.  Rounding margins are recorded (`mrg` column, min $3\cdot10^{-6}$ over
all instances, working precision $3\,000$ digits), so the reduction is exact.

Sample ($k=23.9$; full table `data/struct_n120.csv`, 351 rows):

| $n$ | $\tfrac1n\log[\mathbb Z^2{:}\mathcal K_n]$ | $F_{\rm eff}$ | $\tfrac1n\log\lambda_1$ | $\tfrac1n\log\lambda_2$ | $\lambda_2/\lambda_1$ | $v_1\in\pm\mathcal P$ | $\tfrac1n\log$cone-min | $\rho$ | $i(n)$ | $L(n)$ | $a_{i+1}$ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 27.6117 | $-1.0672$ | $-1.0669$ | $-1.0663$ | 1.01 | yes | $-1.0561$ | 1.24 | 224 | 417 | 3 |
| 40 | 27.9203 | $-0.8949$ | $-0.8947$ | $-0.8935$ | 1.05 | yes | $-0.8885$ | 1.28 | 438 | 845 | 2 |
| 60 | 28.1311 | $-0.7787$ | $-0.7864$ | $-0.7710$ | 2.52 | no | $-0.7662$ | 3.36 | 691 | 1343 | 4 |
| 80 | 28.3597 | $-0.6620$ | $-0.6701$ | $-0.6537$ | 3.71 | no | $-0.6514$ | 4.49 | 915 | 1752 | 8 |
| 100 | 28.2572 | $-0.7075$ | $-0.7081$ | $-0.7064$ | 1.19 | no | $-0.7046$ | 1.42 | 1107 | 2177 | 2 |
| 120 | 28.3614 | $-0.6542$ | $-0.6565$ | $-0.6519$ | 1.74 | yes | $-0.6537$ | 1.41 | 1366 | 2602 | 3 |

and at $k=22.4$ the two extreme cases of the dichotomy sit next to each other:

| $n$ | $\lambda_2/\lambda_1$ | $v_1\in\pm\mathcal P$ | $\rho$ |
|---|---|---|---|
| 110 | 118.57 | **yes** | 1.35 |
| 120 | 49.25 | **no** | 54.89 |

### 2.3 $v_1(n)$ is a continued-fraction convergent

> **Theorem 1 [proved].**  Let $\|\cdot\|$ be a norm on $\mathbb R^2$ strictly monotone in
> $(|x_1|,|x_2|)$ — in particular the $\mathrm{diag}(\lambda_Z,\lambda_N)$-weighted $\ell^2$ norm.
> A shortest nonzero $x\in\mathcal K_n$ is a *relative minimum*: no
> $y\in\mathcal K_n\setminus\{0,\pm x\}$ has $|y_1|\le|x_1|$ and $|y_2|\le|x_2|$ (else
> $\|y\|\le\|x\|$ with $y\ne\pm x$, contradicting minimality after replacing $x$ by $y$).
> Write $\theta=h_{12}/h_{11}$ and $x=(h_{11}(j\theta-p),\ jh_{22})$, $j=x_2/h_{22}\in\mathbb Z$,
> $p=-i\in\mathbb Z$; the relative-minimum property says $|j'\theta-p'|>|j\theta-p|$ for all
> $0<j'<j$ and all $p'$, i.e. $p/j$ is a best approximation of the second kind to $\theta$, hence
> (Lagrange) a convergent.

**Verified: $v_1(n)$ is a convergent in all $351$ instances** ($k$-independent statement, checked
per $k$), and the cone minimiser is likewise always a convergent — never an intermediate
fraction, although the theory allows those.  The balance index obeys
$$\frac{i(n)}{L(n)}=0.5195\pm0.0165\qquad\text{(predicted }\tfrac12+\tfrac{E_2-E_1}{2\kappa}
=0.5+\tfrac{1.294}{2\cdot27.19}=0.5238),$$
i.e. **the first minimum is the convergent at the midpoint of the continued fraction**, displaced
by exactly the $\log$-anisotropy of the metric.  **[measured]**

### 2.4 Is there a recursion $v_1(n-1)\to v_1(n)$?  No.

* Sizes at $n=120$: $\tfrac1n\log|c_Z(v_1)|=14.26,\ 14.46,\ 14.80$ and
  $\tfrac1n\log|c_N(v_1)|=13.02,\ 13.24,\ 13.56$ at $k=22.4,23.0,23.9$; the ratio
  $\tfrac1n\log(c_N/c_Z)=-1.24$ against the prediction $E_1-E_2=-1.2936$.  The coefficients grow
  like $e^{14.3n}$–$e^{14.8n}$: the short vectors do **not** come from a fixed finite set of
  directions.  **[measured]**
* $\gcd(c_Z(v_1(n)),c_Z(v_1(n+1)))$ is huge (median $\sim10^{50}$), but dividing by the forced
  directional contents $m_Z(n),m_Z(n+1)$ leaves a residual gcd with median $1$–$3$ and value $1$
  in about half the pairs.  **All the apparent correlation between consecutive $n$ is the
  Absorption divisibility $m_Z$, and nothing else survives it.**  **[verified $n\le120$]**
* The in-cone indicator $\{i(n)\bmod2\}_n$ passes a Wald–Wolfowitz runs test against a fair coin
  at all three $k$: $z=-0.26,\,+0.42,\,+1.03$.  **[measured]**
* The $460\,365$ partial quotients of $h_{12}/h_{11}$ over the 351 instances:

  | $a$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | $\ge9$ |
  |---|---|---|---|---|---|---|---|---|---|
  | observed | .4146 | .1696 | .0937 | .0592 | .0403 | .0299 | .0229 | .0178 | .1520 |
  | Gauss–Kuzmin | .4150 | .1699 | .0931 | .0589 | .0406 | .0297 | .0227 | .0179 | .1520 |

  $\chi^2=5.09$ on $8$ degrees of freedom ($p\approx0.75$).  **[measured]**

So: no recursion, no finite set of directions, no arithmetic memory between consecutive $n$
beyond the forced content.  What replaces a recursion is §2.3: *the same universal object (a
continued fraction) at every $n$, sampled at its midpoint.*

## 3. The cone, exactly

Write $\mathcal Q$ for the closed first quadrant (the image of $\mathcal P$), $\ell(x)=x_1+x_2$
in scaled coordinates, $\lambda_1^{\mathcal Q}=\min\{|x|:x\in L\cap\mathcal Q\setminus0\}$,
$\rho=\text{cone-min}/\lambda_1$, $\rho_2=\lambda_1^{\mathcal Q}/\lambda_1$.  On $\mathcal Q$ one has
$|x|\le\ell(x)\le\sqrt2\,|x|$, hence $\rho_2\le\rho\le\sqrt2\rho_2$ **[proved]** — the whole
difference between the linear form and the norm is a bounded factor.

> **Lemma A [proved].**  If $(b_1,b_2)$ is Gauss-reduced then at least one of the eight vectors
> $\pm b_1,\pm b_2,\pm(b_1+b_2),\pm(b_1-b_2)$ lies in $\mathcal Q$.  Hence
> $\lambda_1^{\mathcal Q}\le\lambda_1+\lambda_2$.
>
> *Proof.*  Normalise $b_1=\lambda_1(1,0)$, $b_2=\lambda_2(\cos\theta,\sin\theta)$, $\sin\theta>0$,
> $t=\lambda_2/\lambda_1\ge1$, $|t\cos\theta|\le\tfrac12$ (reduction), so $\theta\in[60^\circ,120^\circ]$.
> The four lines through $b_1,b_1+b_2,b_2,b_1-b_2$ have directions $0<\alpha<\theta<\beta<\pi$ with
> $\alpha=\arctan\frac{t\sin\theta}{1+t\cos\theta}$, $\pi-\beta=\arctan\frac{t\sin\theta}{1-t\cos\theta}$.
> Since $1\pm t\cos\theta\ge\tfrac12>0$, both $\alpha$ and $\pi-\beta$ are $<90^\circ$; and
> $\theta-\alpha<90^\circ$, $\beta-\theta<90^\circ$ — trivially when $\theta\le90^\circ$, and when
> $\theta>90^\circ$ because $\alpha\ge\arctan(t\sin\theta)\ge\arctan\sin120^\circ=40.9^\circ>\theta-90^\circ$
> (symmetrically for $\beta$).  All four gaps
> being $<90^\circ$, the eight rays cut the plane into eight open sectors of angle $<90^\circ$, so
> a closed $90^\circ$ sector cannot avoid all of them. $\square$

> **Proposition B (exact cone formula) [proved].**  Every $x\in L\cap\mathcal Q$ with
> $\ell(x)\le\ell_0:=\min\{\ell(\cdot)\ \text{over the eight vectors of Lemma A that lie in }\mathcal Q\}$
> satisfies $|x|\le\ell(x)\le\ell_0\le2\sqrt2\lambda_2$ and therefore, writing $x=ib_1+jb_2$,
> $$|j|=\frac{|\det(b_1,x)|}{\operatorname{covol}}\le\frac{|x|}{\lambda_2\sin\theta}\le
> \frac{2\sqrt2}{\sin60^\circ}<3.27 .$$
> Hence the cone minimum and $\lambda_1^{\mathcal Q}$ are obtained by scanning $|j|\le3$ and, for
> each $j$, the two endpoints of the integer interval
> $I_j=\{i:\ ib_1+jb_2\in\mathcal Q\}$ (a linear form on an interval is minimised at an endpoint;
> for the norm one adds the two roundings of the interior minimiser).  This is a **complete**
> $O(1)$ algorithm, in contrast to the $\pm R$ sweep of `lattice/positivity/cone80.gp`.

> **Proposition C (dichotomy) [proved; 0 violations in 351 instances].**
> $$\pm b_1\in\mathcal Q:\qquad \rho_2=1,\qquad 1\le\rho\le\sqrt2;$$
> $$\pm b_1\notin\mathcal Q:\qquad \tfrac{\sqrt3}{2}\tfrac{\lambda_2}{\lambda_1}\le\rho_2\le
> 1+\tfrac{\lambda_2}{\lambda_1},\qquad
> \tfrac{\sqrt3}{2}\tfrac{\lambda_2}{\lambda_1}\le\rho\le\sqrt2\bigl(1+\tfrac{\lambda_2}{\lambda_1}\bigr).$$
> *Proof.*  First case: $\pm b_1$ is a cone vector, so $\lambda_1^{\mathcal Q}=\lambda_1$ and
> $\text{cone-min}\le\ell(\pm b_1)\le\sqrt2|b_1|$; and $\ell\ge|\cdot|\ge\lambda_1$ on $\mathcal Q$.
> Second case: every $x\in L\cap\mathcal Q$ has $j\ne0$, so
> $|x|\ge|\det(b_1,x)|/\lambda_1\ge\operatorname{covol}/\lambda_1=\lambda_2\sin\theta\ge\tfrac{\sqrt3}2\lambda_2$;
> the upper bound is Lemma A. $\square$

> **Corollary D [proved].**  $\rho\le\sqrt2$ forces $\pm b_1\in\mathcal Q$ **or**
> $\lambda_2/\lambda_1\le2\sqrt2/\sqrt3=1.633$.  Measured: 18 of the 351 instances have
> $\rho\le\sqrt2$ with $v_1$ outside the cone, and all 18 have $\lambda_2/\lambda_1\le1.306$.

**Why the median of $\rho$ is $\sqrt2$.**  Proposition C says $\rho\le\sqrt2$ essentially *iff*
the parity is even, an event of probability $\tfrac12$; so the median of $\rho$ sits at the
boundary $\sqrt2$ of the two branches.  Measured $\Pr[\rho\le\sqrt2]=0.581$ (Catalan, 351) and
$0.5424$ (Haar-random, $2\times10^5$); median $\rho=1.4023$ vs $1.4112$.  The "$\sqrt2$ to four
figures" of `POSITIVITY_PROGRAM.md` §3.3 and §4.3(a) is this, not a numerical coincidence.

## 4. P2$'$, structurally

### 4.1 Equivalent forms

With $\rho_2\le\rho\le\sqrt2\rho_2$ (§3) and Proposition C:

| form | statement | status |
|---|---|---|
| P2$'$ | $\liminf_n\frac1n\log\rho_n=0$ | the conjecture |
| P2$'_{\ \rm norm}$ | $\liminf_n\frac1n\log\bigl(\lambda_1^{\mathcal Q}/\lambda_1\bigr)=0$ | **equivalent** [proved] |
| P2$'_{\ \rm dich}$ | i.o. $n$: $v_1(n)\in\pm\mathcal P$, **or** $\frac1n\log\frac{\lambda_2}{\lambda_1}\to0$ along it | **equivalent** [proved] |
| P2$'_{\ \rm par}$ | i.o. $n$: $i(n)$ even, **or** $a_{i(n)+1}(n)=e^{o(n)}$ | equivalent given $\lambda_2/\lambda_1\asymp a_{i+1}$ [measured: $\lambda_2/\lambda_1\in[0.019,1.38]\cdot a_{i+1}$, correlation of logs $0.79$–$0.82$] |
| P2$'_{\ \rm short}$ | i.o. $n$: $\lambda_1(\mathcal K_n)\ge e^{-o(n)}\sqrt{\operatorname{covol}}$ | **sufficient**, not necessary [proved] |

The last line is worth isolating: $\lambda_2/\lambda_1=\operatorname{covol}/\lambda_1^2$, so
"sub-exponentially skew" $=$ "no exponentially short vector".  P2$'$ does **not** require it —
an exponentially short $v_1$ inside the cone is harmless, and this is not vacuous: at $k=22.4$,
$n=110$ the lattice has $\lambda_2/\lambda_1=118.6$ and still $\rho=1.35$, because the parity is
even.

### 4.2 The data to $n=120$ (`data/struct_n120.csv`, 351 instances, all `ok`)

| $k$ | in-cone | median $\rho$ | $\max\rho$ | $\Pr[\rho\le\sqrt2]$ | $\frac1n\log\rho$ fit | $\log\rho$ fit | longest run of bad $n$ |
|---|---|---|---|---|---|---|---|
| 22.4 | 0.521 | 1.409 | 54.9 | 0.573 | $(-0.00044\pm0.00009)n+0.048$ | $(+0.0034\pm0.0023)n+0.548$ | 5 |
| 23.0 | 0.556 | 1.385 | 424.7 | 0.598 | $(-0.00062\pm0.00013)n+0.060$ | $(-0.0024\pm0.0026)n+0.806$ | 4 |
| 23.9 | 0.513 | 1.410 | 42.3 | 0.573 | $(-0.00036\pm0.00006)n+0.040$ | $(+0.0013\pm0.0022)n+0.623$ | 4 |

Two remarks on the regressions.

* The model used in `POSITIVITY_PROGRAM.md` §3.3, "$\tfrac1n\log\rho$ linear in $n$", is
  mis-specified: if $\log\rho=O(1)$ then $\tfrac1n\log\rho\sim c/n$ and a straight-line fit
  *must* return a small negative slope, which is what happens here now that the range is longer
  and the cone minima are correct ($-0.0004$ to $-0.0006$, nominally $4\sigma$).  The
  well-specified model is $\log\rho$ against $n$; its slope is $0\pm0.003$ at all three $k$,
  i.e. **any exponential growth of $\rho$ is below $0.008$ per unit $n$ at $2\sigma$** — against
  $|F(23.9)|=0.537$ (factor $67$) and $|F(22.4)|=0.0169$ (factor $2$).  As before, the mechanism
  is best supported exactly where $F$ is most negative.  **[measured]**
* P2$'$ only asks for infinitely many $n$.  The relevant statistic is the run structure: the
  longest run of consecutive $n\le120$ with $\rho>\sqrt2$ is $5,4,4$ at the three $k$, so **in
  every window of six consecutive $n$ there is an $n$ with $\rho\le\sqrt2$** — a far stronger
  statement than the average behaviour, and exactly what P2$'$ needs.  **[verified $n\le120$]**

### 4.3 The random-lattice model

`stats.py` samples $2\times10^5$ Haar-random unimodular planar lattices (shape drawn from
$\tfrac3\pi y^{-2}dx\,dy$ on the modular fundamental domain, orientation uniform) and runs the
*same* exact cone algorithm:

| statistic | Haar ($2\times10^5$) | Catalan (351) |
|---|---|---|
| $\Pr[v_1\in\pm\mathcal Q]$ | 0.4993 | 0.5299 |
| median $\rho$ | 1.4112 | 1.4023 |
| $\Pr[\rho\le\sqrt2]$ | 0.5424 | 0.5812 |
| $\Pr[\rho>20]$ | 0.0270 | 0.0399 |
| median $\lambda_2/\lambda_1$ | 1.9336 | 2.0221 |
| $\Pr[\lambda_2/\lambda_1>20]$ | 0.0482 | 0.0598 |

Kolmogorov–Smirnov distance between the two $\rho$-distributions: $0.0593$, against the $5\%$
critical value $1.358/\sqrt{351}=0.0725$.  **The Catalan lattices are not distinguishable from
random lattices by any of these statistics**, and the Gauss–Kuzmin fit of §2.4 says the same
thing on the continued-fraction side.  Under that model P2$'$ holds with room to spare
($\Pr[\rho>e^{\varepsilon n}]$ decays like $e^{-\varepsilon n}$, Borel–Cantelli); the Haar tail
$\Pr[\rho>R]\approx0.54/R$ predicts $3.2$ instances with $\rho>20$ per $k$ (observed $5,5,4$ out
of $117$) and a maximum of order $63$ (observed $54.9,\ 424.7,\ 42.3$).  **[measured]**

## 5. The kernel vector

### 5.1 Exact form

If $G=a/b$ then $c^{\rm ker}=(aV_n-bU_n,\ -(aX_n-bY_n))$ annihilates the form identically.
Verified exactly (`kerchk`, surrogates $b\le10^6$, $6\le n\le24$, both $k$; `data/kernel_n120.csv`):

* $c^{\rm ker}\!\cdot\!(X_nG-Y_n,\ V_nG-U_n)=0$ **exactly** (rational arithmetic) — 0 failures;
* $c^{\rm ker}\in\mathcal K_n$ **without any multiple**: indeed
  $(q,p)=-\bigl(h_n/T_n\bigr)(b,a)$ with $h_n/T_n\in\mathbb Z$ the odd part of the mixed minor —
  0 failures.  (`06_threshold.tex` states the kernel output as $M_n(b,a)$; the correct
  normalisation is the odd part of $h_n$, which is much smaller.)
* its two scaled coordinates are **equal**: $|c^{\rm ker}_Z|\lambda_Z=|c^{\rm ker}_N|\lambda_N=b\,M_n\lambda_Z\lambda_N$
  — 0 failures.  So the kernel direction is the *anti-diagonal* of the oriented quadrant: it is
  the single worst direction for the cone, at $45^\circ$ outside it.
* Length: $\ \|c^{\rm ker}\|=\sqrt2\,b\,M_n\lambda_Z\lambda_N=\sqrt2\,b\,\dfrac{M_n}{[\mathbb Z^2:\mathcal K_n]}\operatorname{covol}
  =\sqrt2\,b\,\gcd(X_n,Y_n,V_n,U_n)\operatorname{covol}$, by §2.1.  **[proved + verified]**
* The *primitive* kernel vector is shorter by $g/t$ where $g=\gcd(aV-bU,aX-bY)$ and $t$ is the
  smallest multiplier returning it to $\mathcal K_n$; measured $t\in\{1,2,7,17,31\}$ and
  $\tfrac1n\log(g/t)\in[0.177,0.597]$, so the true kernel generator is *shorter* than the estimate
  above — the effect strengthens the conclusion below and is $(a,b)$-dependent, so it is not used
  in the rigorous statements.

### 5.2 When does it become the first minimum?  The corrected proposition

Set $b_{\rm crit}(n)=\lambda_1/\bigl(\sqrt2\,\gcd(X,Y,V,U)\operatorname{covol}\bigr)$; the kernel
vector beats $\lambda_1$ iff $b<b_{\rm crit}$, i.e.
$$\frac1n\log b_{\rm crit}=\frac1n\log\lambda_1-2F_{\rm eff}-\frac1n\log\gcd(X,Y,V,U)
\ \longrightarrow\ |F(k)|-\tfrac{c}{2},\qquad c=\lim\tfrac1n\log\gcd(X,Y,V,U)$$
$$\text{(because }\operatorname{covol}=e^{2F(k)n}/\gcd(X,Y,V,U),\ \text{ i.e. }F_{\rm eff}=F-c/2) .$$
Measured (`bcrit`, $n\le120$):

| $k$ | $\tfrac1n\log b_{\rm crit}$ at $n=20,60,100,120$ | $\max_{n\le120}\log b_{\rm crit}$ |
|---|---|---|
| 22.4 | $-0.193,\ -0.165,\ -0.005,\ -0.068$ | $3.11$ (at $n=37$) |
| 23.0 | $+0.009,\ +0.050,\ +0.200,\ +0.152$ | $20.16$ (at $n=91$) |
| 23.9 | $+0.357,\ +0.356,\ +0.516,\ +0.466$ | $55.91$ (at $n=120$) |

Two consequences.

* **The claim of `06_threshold.tex` ("scaled length $\asymp be^{2Fn}$, first minimum for
  $n\gtrsim\log b/|F|$") is right in shape and needs the deficit correction**: the length is
  $b\gcd(X,Y,V,U)\operatorname{covol}$, and the crossover is at $n\approx\log b/(|F|-c/2)$.
  With the $n\le120$ deficits still at $\tfrac1n\log\gcd\approx0.18$, the criterion **fails at
  $k=22.4$** in the computed range ($b_{\rm crit}<1$: no rational at all is excluded there),
  holds comfortably at $k=23.9$.  Since §2.1 shows $c=0$ asymptotically, the paper's statement is
  correct in the limit at every $k>k_*$, but it is *not* the mechanism by which P2$'$ fails at
  finite $n$.  **[verified + measured]**
* **The mechanism that does work at every $k$ is trivial and stronger.**  If $G=a/b$ then every
  cone vector has $q_nG-p_n>0$ (Zudilin/Nesterenko positivity, a theorem about the number $G$)
  and $q_nG-p_n\in\tfrac1b\mathbb Z$, so cone-min $\ge1/b$ while $\lambda_1\le(2/\sqrt3)^{1/2}
  e^{F_{\rm eff}n}$.  Hence $$\rho_n\ \ge\ \tfrac1b\,\bigl(\tfrac{\sqrt3}{2}\bigr)^{1/2}e^{|F_{\rm eff}|n}
  \quad\text{for all }n,$$ and P2$'$ fails **at every $k>k_*$, at the rate $|F(k)|$**, for every
  rational $G$.  **[proved]**
* Read backwards, the same inequality is a rigorous (if weak) by-product of the computation: were
  $G=a/b$ with $b\le b_{\rm crit}(n)$, $\mathcal K_n$ would contain a vector shorter than the
  measured $\lambda_1$.  Hence $G$ is not a rational of denominator $\le10^{24.3}$
  ($k=23.9$, $n=120$).  The continued fraction of $G$ gives incomparably more; the point is only
  that this is *all* the lattice sees.  **[verified]**

### 5.3 The horizon: how far a rational surrogate reproduces the table

Everything in §§2–4 depends on $G$ only through the ratio $r_n=(X_nG-Y_n)/(V_nG-U_n)$.  For
$G^*$ with $\delta=|G-G^*|$, the ratio changes at index $n$ once $|X_n|\delta\gtrsim|X_nG-Y_n|$,
i.e. at $n_f=\log(1/\delta)/(\xi-E_1)$ with $\xi-E_1=27.52-13.10=14.43$.  Running the entire
analysis with $G^*=\mathrm{bestappr}(G,10^E)$ and comparing $v_1(n)$ vector by vector
(`run_surrogate.gp`):

| $E$ | $\log(1/\delta)$ | predicted $n_f$ | **measured** first $n$ with $v_1(G^*)\ne v_1(G)$ |
|---|---|---|---|
| 40 | 184.2 | 12.8 | **13** |
| 80 | 368.4 | 25.5 | **25** |
| 160 | 736.8 | 51.0 | **50** |
| 320 | 1472.8 | 102.1 | **99** |

identical for all three $k$ (the horizon is a property of $r_n$, not of $T_n$).  Proposition
(Control) of `06_threshold.tex` uses $E=320$ at $n=98$: **one index below the horizon.**
**[verified]**

Combining with §5.2, for a surrogate of denominator $b$,
$$\frac{n^*}{n_f}=\frac{\log b/|F_{\rm eff}|}{2\log b/(\xi-E_1)}=\frac{\xi-E_1}{2|F_{\rm eff}|}
=11.0\ (k{=}23.9),\quad21.1\ (23.0),\quad53.8\ (22.4),$$
independent of $b$: **the regime in which a rational surrogate still reproduces the $G$-table and
the regime in which its rationality is visible are separated by a fixed factor**, so no choice of
surrogate can exhibit the failure.  Directly: past $n_f$ a fixed rational $G^*$ degenerates the
*other* way — its covolume grows ($F_{\rm eff}\to+13.6$ at $E=10$), and $\rho$ stays $O(1)$ for
its own reasons; moreover, since $r_n\to X_n/V_n$ for *any* fixed rational, all far surrogates
produce literally the same table (checked: $E=10$ and $E=20$ agree from $n=8$ on).  A rational
number cannot be made to behave like the hypothesis "$G$ rational".  **[verified]**

## 6. Back to the rows

The cone-minimising vector, written in the original coordinates
$c=(s_Zu,\,s_Nv)$, gives
$$q_n=\frac{c_ZX_n+c_NV_n}{M_n},\qquad p_n=\frac{c_ZY_n+c_NU_n}{M_n},\qquad
q_nG-p_n=u\lambda_Z+v\lambda_N>0,$$
with $q_n,p_n\in\mathbb Z$ and the reported value re-derived from the returned integers in all
351 instances (`ok` column).  Sizes: $\tfrac1n\log q_n\approx14.3$–$14.8$, $\tfrac1n\log|q_nG-p_n|
=F_{\rm eff}+O(1/n)\approx-0.13$ to $-0.65$; so $\delta=-\log|q_nG-p_n|/\log q_n\approx0.01$–$0.045$.

* **$v_1$ is not a Hermite–Padé convergent of $G$**: with $\log q_n\approx14.3n$ and
  $|q_nG-p_n|\approx e^{-0.65n}$, $p_n/q_n$ is a far worse approximation to $G$ than a continued-
  fraction convergent of the same height by a factor $e^{13.6n}$.  It is a convergent — but of
  the *congruence datum* $h_{12}/h_{11}$ carried by the rows modulo $M_n$ (for the
  single-congruence lattice of §7, literally $-U_nY_n^{-1}\bmod S_n$), not of $G$.
  **[measured + proved (Thm 1)]**
* **Sign pattern.**  In the oriented coordinates the cone is $\{u\ge0,v\ge0\}$; in the original
  ones it is $\{(-1)^nc_Z\ge0,\ c_N\ge0\}$.  Normalising $v_1$ to have $v>0$, its first
  coordinate is the convergent error $h_{11}(j\theta-p)$, whose sign alternates: the pattern is
  $(+,+)$ for even $i(n)$ and $(-,+)$ for odd.  There is **no structural reason for
  $v_1$ to be in the cone** — the question is the parity of a continued-fraction index, and
  §2.4 shows that index behaves like a random one.  This is the precise answer to "is $v_1$
  always in the cone, or its negative?": it is, for exactly the even half of $n$.
* **Absorption is visible in the coefficients.**  Every $c\in\mathcal K_n$ has $m_Z\mid c_Z$ with
  $\tfrac1n\log m_Z\approx2.39$ (§2.1) — the Zudilin coefficient of $v_1$ carries $\approx2.4n$
  nats of forced divisibility, the same for all $k$, and the kernel vector carries it too.

## 7. The single-congruence lattice of Theorem (Absorption)(4)

For $\mathcal K^S_n=\{c:\ c_ZY_n+c_NU_n\equiv0\bmod S_n\}$ (`analS`, $4\le n\le120$,
`data/structS_n120.csv`): $\gcd(Y_n,U_n,S_n)>1$ at every $n$, index $e^{11.85n}$ at $n=120$
(against $S_n=e^{12n}$), $\lambda_1=e^{+7.66n}$ — the forms are enormous, so this lattice is
useless for approximation.  Its *structure* is identical:

* $v_1$ is a continued-fraction convergent of $h_{12}/h_{11}$: **117/117**;
* $v_1\in\pm\mathcal P\iff i(n)$ even: **117/117**;
* in-cone frequency $0.479$, median $\rho=1.4122$, $\Pr[\rho\le\sqrt2]=0.521$, median
  $\lambda_2/\lambda_1=2.086$.

So the $\sqrt2$/parity phenomenon is a property of *congruence lattices with a coordinate cone*,
not of the particular modulus — as `POSITIVITY_PROGRAM.md` §4.3(a) already suspected from the
13 332 pair lattices.

## 8. The cleanest equivalent of P2$'$, and whether it is attackable

> **P2$'$ (final form).**  For infinitely many $n$, the convergent index $i(n)$ of
> $h_{12}(n)/h_{11}(n)$ at which the weighted balance
> $\lambda_Z|h_{12}q_i-p_ih_{11}|\asymp\lambda_Nh_{22}q_i$ occurs is **even**; or, failing that,
> the next partial quotient $a_{i(n)+1}$ is $e^{o(n)}$.

Equivalently (Prop. C): for infinitely many $n$ the shortest vector of $\mathcal K_n$ lies in the
closed coordinate quadrant, or the lattice is sub-exponentially skew.  Both disjuncts are
statements about one explicit continued fraction, and the first is a **parity**.

Is it attackable?  Three observations, in decreasing order of comfort.

1. **The lattice side is now completely explicit and needs no analysis of $G$**: $h_{11},h_{12},h_{22}$
   are computable integers determined by the rows modulo $M_n$, and Gauss reduction in dimension
   two is a continued fraction.  Everything that is *not* about $G$ has been reduced to: at which
   index does the balance occur, and what is the parity there.
2. **$G$ enters at exactly one place, and it is the place that decides.**  The index $i(n)$ is
   determined by the balance, hence by $\log(\lambda_N/\lambda_Z)=\log\frac{V_nG-U_n}{|X_nG-Y_n|}$;
   moving $G$ by $\delta$ moves this quantity as soon as $\delta>e^{-(\xi-E_1)n}$ (§5.3), i.e.
   changes the index and therefore the parity.  A proof of P2$'$ must know $\log(\lambda_N/\lambda_Z)$
   to within one continued-fraction step at index $\approx11n$ — that is, must know $G$ to
   $\approx14.4n$ nats and must control which of two consecutive integers a real number falls
   between.  This is precisely the paper's "must distinguish $G$ from $G^*$ at precision growing
   with $n$", now with the constant: **$(\xi-E_1)=14.43$ nats of $G$ per unit $n$**, and with the
   quantitative reason why no finite computation helps (the factor $11$–$54$ of §5.3).
3. **It is true for almost every target, and false for every rational one.**  Under the
   random-lattice model, which the data supports at every statistic tested (§4.3, §2.4), P2$'$
   holds with density $\tfrac12$; and by §5.2 it fails at rate $|F(k)|$ for every rational $G$.
   So P2$'$ separates $\mathbb Q$ from a full-measure set — it is an irrationality statement,
   confirming `06_threshold.tex`, and no argument that is insensitive to $e^{-14.4n}$-perturbations
   of $G$ can prove it.

**Honest assessment.**  The reduction achieved here is genuine but is a reduction of *form*, not
of difficulty: P2$'$ is now a parity statement about one explicit continued fraction instead of a
statement about lattice minima, and the parity is equidistributed in every test.  There is no
known technique for forcing the parity of a continued-fraction index at a prescribed depth for an
explicit family; and by point 2 any such technique would have to consume $14.43n$ nats of
information about $G$, which is exactly the analytic input the paper says is unavoidable.  I do
not believe this is attackable as stated.  What it *does* make possible, and what is new, is a
sharp finite-$n$ statement: for each $n$, "$i(n)$ is even" is decidable in seconds, and the set of
$n$ for which it holds has density $\tfrac12$ with no run longer than $5$ — so a proof of P2$'$
needs only a **lower bound on the density of even indices**, not any control of individual $n$.
That is the only shape of statement I can see that is both sufficient and not obviously as hard as
irrationality itself. **[open]**

## 9. Corrections to the record

1. `lattice/positivity/data/cone_n80.csv`: $51$ of $231$ cone minima are not minimal (the $\pm R$
   sweep of `cone80.gp` scans the wrong index — $|i|$ can be as large as $\lambda_2/\lambda_1$ —
   and it never tests the lattice point adjacent to an interval endpoint that is the zero vector).
   The corrected values are in `data/struct_n120.csv`; all differences lower $\rho$.  The largest
   corrections: $k=22.4,n=4$ ($\rho\ 2.18\to1.27$), $k=22.4,n=23$ ($2.34\to1.62$),
   $k=23.0,n=22$ ($2.39\to1.53$), $k=22.4,n=51$ ($1.83\to1.14$).
2. `POSITIVITY_PROGRAM.md` §3.4: the gcd loss is $\gcd(X_n,Y_n,V_n,U_n)$, of size $e^{O(\log n)}$,
   so $c=0$; the $\log n/n$ fit is the correct one and $\kappa_n\to\sigma(k)$.
3. `POSITIVITY_PROGRAM.md` §3.3: the regression model $\tfrac1n\log\rho=\alpha n+\beta$ is
   mis-specified; use $\log\rho=\alpha n+\beta$ (slopes $0\pm0.003$).
4. `06_threshold.tex`, Proposition after Lemma P2: the kernel output is
   $(q,p)=-(h_n/T_n)(b,a)$, not $M_n(b,a)$; the scaled length is
   $\sqrt2\,b\gcd(X,Y,V,U)\operatorname{covol}$; and the reason P2 fails for rational $G$ should be
   the bound cone-min $\ge1/b$ (valid at every $k$) rather than the kernel becoming the first
   minimum (which additionally needs $|F|>c/2$, false in the computed range at $k=22.4$).
5. `06_threshold.tex` uses two different lattices under one name: §5.1 defines $\mathcal K_n$ by
   integrality of the combination after division by $M_n$ (two congruences), while Theorem
   (Absorption)(4) argues from "the defining congruence $c_ZY+c_NU\equiv0\pmod{S_n}$" (one
   congruence).  §7 shows the structural statements coincide, but the numbers do not
   (index $e^{11.85n}$ against $e^{27.3n}$, $\lambda_1=e^{+7.66n}$ against $e^{-0.65n}$).

## 10. Status of every claim

**Proved.**  Theorem 1 ($v_1$ is a convergent); the parity criterion; Lemma A; Proposition B
(complete $|j|\le3$ cone algorithm); Proposition C (dichotomy) and Corollary D; the equivalences
P2$'\iff$P2$'_{\rm norm}\iff$P2$'_{\rm dich}$; the kernel identity and its exact length; the
failure of P2$'$ for rational $G$ at rate $|F_{\rm eff}|$ via cone-min $\ge1/b$.

**Verified (exact, $4\le n\le120$, $k\in\{22.4,23.0,23.9\}$, 351 instances).**  All four identities
of §2.1; $v_1$ and the cone minimiser are convergents (351/351), never intermediate fractions; the
parity criterion (351/351, no exceptions); no violation of Prop. C; integrality of every returned
$(q_n,p_n)$ and agreement of $|q_nG-p_n|$ with the reported cone minimum; the surrogate horizons
$n_f=13,25,50,99$; the exclusion $G\ne a/b$ for $b\le10^{24.3}$.

**Measured.**  All regressions of §4.2; the Gauss–Kuzmin fit; the balance-index ratio
$0.5195\pm0.0165$; the Haar comparison and the KS test; $\lambda_2/\lambda_1$ versus $a_{i+1}$;
$\tfrac1n\log m_Z\approx2.39$; $\xi-E_1=14.43$.

**Open.**  P2$'$ itself; the density of even indices $i(n)$ (measured $0.53\pm0.03$, conjecturally
$\tfrac12$); whether $\log\gcd(X_n,Y_n,V_n,U_n)=O(\log n)$ provably; a proof that
$\lambda_2/\lambda_1\le Ca_{i(n)+1}$.

## 11. Scripts and data

`lattice/p2_structure/` — PARI/GP unless stated.  Run pattern (concatenate, do not `read()`):

    cat lattice/positivity/rows_pos.gp lattice/p2_structure/p2core.gp \
        lattice/p2_structure/p2run.gp lattice/p2_structure/run_main.gp > run.gp && gp -q run.gp

| file | prepend | what |
|---|---|---|
| `build_rows.gp` | `rows_pos.gp` | caches the exact integer rows $X,Y,V,U$ for a range of $n$ |
| `p2core.gp` | `rows_pos.gp` | `rdrows`, `kfull`, `kcong`, weighted Gauss reduction `gred`, the complete cone scan `coneiv`/`conescan`, `eight` (Lemma A) |
| `p2run.gp` | + `p2core.gp` | `anal` (the full per-$(k,n)$ line), `analS` (single-congruence lattice), `cfclass`, `cfhist` |
| `kernel.gp` | + `p2core.gp` | `kerchk` (exact kernel identity for a surrogate), `bcrit` |
| `verify.gp` | + `p2core.gp` | `chk` (the five exact identities of §2.1), `qpchk` (the kernel output identity of §5.1) |
| `run_main.gp`, `run_hist.gp`, `run_cong.gp`, `run_kernel.gp`, `run_surrogate.gp` | drivers | |
| `stats.py` | — | §§1–5 statistics, regressions, Gauss–Kuzmin, Haar-random Monte Carlo |
| `horizon.py` | — | §5.3: first $n$ at which a surrogate's $v_1$ differs from the true one |

| data file | rows | content |
|---|---|---|
| `data/rows_all.txt` | 117 | exact $n,X_n,Y_n,V_n,U_n$ for $4\le n\le120$ |
| `data/struct_n120.csv` | 351 | the full structural table (27 columns; header in the file) |
| `data/vectors_n120.txt` | 351 | exact $v_1,v_2$ and the cone minimiser in oriented coordinates |
| `data/cfhist_n120.csv` | 351 | partial-quotient histograms |
| `data/structS_n120.csv` | 117 | the single-congruence lattice |
| `data/kernel_n120.csv` | 545 | kernel identity checks and $b_{\rm crit}$ |
| `data/g0_n120.csv` | 117 | $\gcd(X_n,Y_n,V_n,U_n)$, its $2$-valuation and odd part |
| `data/horizon.csv` | 4 | the surrogate horizons of §5.3 (the $1.8$ MB surrogate dumps themselves are regenerable with `run_surrogate.gp` and are not kept) |

**Cost.**  Building the rows to $n=120$ (the Nesterenko partial-fraction solve is $O(n^{4.1})$)
is $\approx1.8$ h of CPU, done once and cached; every table above is then produced in under a
minute.  Peak memory $6.8$ GB at $n=120$.  Precision: $G$ at `\p 3000` throughout; the smallest
Gauss-reduction rounding margin met anywhere is $3\cdot10^{-6}$, so the reductions are exact.
