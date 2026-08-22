# Non-congruence Beukers rows: theory, a finite classification, and two scans

*Claude (Opus 5), 2026-08-22. Scripts and data: `lattice/noncongruence_scan/`.
Executes the programme of `THEORY_NOTES_06_deeper_object.md` §4(b): "a scan of
genus-zero **non-congruence** groups with integral-after-scaling Hauptmoduls, scored by
$\log(1/\lambda_2)-k-\log\lambda$, is a well-posed search no one has done; the known
positive instance is Beukers'." Builds on `SQRT_APERY.md` (the $w=2$ theory and the
attribution to Beukers 1987 Thm 3), `ROOT_ROWS.md` (Theorems R1–R4), `SPORADIC_SCAN2.md`
(scan method).*

---

## 0. Verdict first

Three theorems, one finite classification, two independent exhaustive scans, and a
negative answer with a sharp reason.

* **Theorem N1 (the score splits).** For a genus-zero group $\Gamma$ commensurable with
  $\mathrm{PSL}_2(\mathbf Z)$, a Hauptmodul-type parameter $t=q+O(q^2)$ at a MUM cusp, and a
  **weight-two eta quotient** $F$ with $F(0)=1$, the row
  $a_n=\lambda^n[t^n]\sqrt F$ is a second-order Apéry row with $k=2$, and
  $$\boxed{\ \operatorname{score}=\log|t_2|-2-\log\lambda\ }$$
  where $t_2$ is the **second-smallest** of the values $t(P)$, $P$ running over the
  cusps and elliptic points of $\Gamma$. The archimedean factor $|t_2|$ depends only on
  $(\Gamma,t)$ — **not on $F$** — and $\lambda\in\{1,2,4\}$ depends only on the $2$-adic
  structure of $F$. Non-congruence $\Leftrightarrow\lambda>1$ (Calegari–Dimitrov–Tang).
  Hence a positive score needs
  $$|t_2|>e^2\lambda=\begin{cases}7.389 & \lambda=1\ \text{(congruence)}\\
  14.778 & \lambda=2\\ 29.556 & \lambda=4.\end{cases}$$
  Beukers' row attains $|t_2|=(\sqrt2+1)^4=33.9706$ with $\lambda=4$: margin $1.149$.

* **Theorem N2 (four special points; finiteness).** A **three-term** such row forces
  $\Gamma$ to have exactly four special points (cusps $+$ elliptic points). By
  Gauss–Bonnet its covolume index in $\mathrm{PSL}_2(\mathbf Z)$ is
  $\mu=6\bigl(2-\sum_i 1/e_i\bigr)\le 12$. **There are therefore only finitely many
  hosts.** Inside $\mathrm{PSL}_2(\mathbf Z)$ there are exactly **28** of them
  (up to conjugacy), of which **16 are congruence and 12 are non-congruence**; the
  non-congruence ones have index $7$, $9$, $10$ and signatures $(0;2,3;2)$, $(0;2;3)$,
  $(0;3;3)$ — four groups each (§3, `06_groups.py`). The index-$12$ four-cusp case is
  exactly Beauville's six families, all congruence, recovered as a check.

* **Theorem N3 (unequal exponents kill the score).** If the two finite singular points
  $t_1,t_2$ carry **different** exponent differences $1/e_1\ne1/e_2$ then they are
  rational, hence (with the integral normalisation $P(0)=1$, $P,Q,R\in\mathbf Z[t]$)
  $\lambda_1,\lambda_2\in\mathbf Z$, so $|\lambda_2|\ge1$ and
  $\operatorname{score}\le-k<0$. Combined with N2 this leaves, for a non-congruence
  **host group**, only the "both finite points are cusps" case — the Zagier class —
  and that class is scanned exhaustively in §4.

* **Scan A (recurrence classes, `01_class_scan.c`).** Complete integer search of the
  five classes $e_1=e_2=e\in\{2,3,4,6,\infty\}$ of second-order MUM rows, with the
  Casoratian (non-vanishing) filter and the measured $k$. Result:
  **Apéry's/Beukers' row $(\alpha,\gamma,\delta,\zeta)=(136,10,16,4)$ is the unique
  positive score in every box searched**, the runner-up being its own rescaling at
  $-0.554$ and then the level-8 $T$ row at $-1.624$. Beyond it, every positive score in
  the scan is either a **Casoratian degeneracy** ($Q(n_0)=0$, so $W_n\equiv0$ and the
  companion is a rational multiple of the row — the criterion is vacuous) or a **classical
  Legendre/Padé $\log$ row** (e.g. $a_n=\binom{2n}nP_n(75)$, $\exp(4\xi)=\tfrac{38}{37}$),
  or **Apéry's $\zeta(2)$ row** $(11,3,-1,0)$ at $+0.406$. Exactly **fourteen** rows have
  positive score *and* $k\ge2$; twelve are the $\log$ family, one is Apéry's $\zeta(2)$,
  one is Beukers'.

* **Scan B (modular, `03b_exact.gp`).** All eta-quotient pairs $(t,F)$ on $\Gamma_0(N)$,
  $N\le60$: $t$ of weight $0$ with $\operatorname{ord}_\infty t=1$ and divisor degree
  $\le4$, $F$ holomorphic of weight $2$ with $\operatorname{ord}_\infty F=0$ —
  $3\,036\,854$ pairs. Result: $739$ pairs admit a recurrence of order $\le4$ and degree
  $\le3$ ($485$ with $\lambda=1$, $127$ with $\lambda=2$, $127$ with $\lambda=4$), falling
  into $104$ distinct invariant classes. **Exactly two have $|\lambda_2|<1$**: Beukers'
  ($+0.1392$, $\lambda=4$, non-congruence) and the level-$8$ $T$ row ($-1.6235$,
  $\lambda=1$, congruence). Everything else has $|\lambda_2|\ge1$ — no decay at all.

* **Is it already in Beukers 1987?** Yes. `papers/beukers1987.txt` Theorem 3 **is** the
  unique winner, and his §3 ($\Gamma_1(5)$) is the same construction at the other perfect
  level. Nothing in the scan improves on it, and no new irrationality claim is made here.

---

## 1. The construction, and why the score splits

**Setup.** $\Gamma$ a genus-zero Fuchsian group commensurable with
$\mathrm{PSL}_2(\mathbf Z)$ with a distinguished cusp $\infty$; $t$ a Hauptmodul of
$\Gamma$ normalised $t=q_h+O(q_h^2)$ in the local parameter $q_h=e^{2\pi i\tau/h}$ at
$\infty$, with integral $q$-coefficients; $F$ a weight-two form on $\Gamma$ (possibly with
a character) with $F=1+O(q)$ and **divisor supported on the cusps**. Put
$$g:=\sqrt F,\qquad a_n:=\lambda^n[t^n]g .$$

**(i) $g$ is a weight-one form on an index-$\le2$ subgroup.** $F$ is nowhere zero on
$\mathfrak H$ (its divisor is cuspidal), so $\sqrt F$ is single-valued on $\mathfrak H$;
its automorphy factor is $\pm$ the weight-one one, i.e. $g\in M_1(\Gamma',\chi)$ for
$\Gamma'=\ker\chi$, $[\Gamma:\Gamma']\le2$. Taking $F$ an **eta quotient** guarantees the
divisor hypothesis automatically; this is exactly Beukers' "$E$ has only zeros and poles
in the cusps".

**(ii) The row is second order.** For a genus-zero $\Gamma$ and a Hauptmodul $t$, a
weight-$w$ form written as a function of $t$ satisfies $\operatorname{Sym}^w$ of the
weight-one equation. Hence $F(t)$ satisfies $\operatorname{Sym}^2L_1$ and, by the descent
Theorem R2 of `ROOT_ROWS.md`, $g=F^{1/2}$ satisfies the **second-order** $L_1$. So the row
$a_n$ is Apéry-like of order two — no hypothesis beyond genus zero is needed.

**(iii) $\lambda\in\{1,2,4\}$ and it is the only arithmetic input.** $t\in q\mathbf Z[[q]]$
monic gives $q\in t\mathbf Z[[t]]$ (Lagrange), so $F(t)\in1+t\mathbf Z[[t]]$; Theorem R1
($w=2$) gives $4^n[t^n]\sqrt{F}\in\mathbf Z$, and R1$'$ sharpens it to
$$\lambda=\max\bigl(1,2^{\,2-e}\bigr),\qquad e:=\min_{n\ge1}v_2\bigl([t^n]F\bigr).$$

**(iv) $k=2$.** Theorem R3: with $\Psi=g^3u$ of weight three, $B=g\cdot\theta_q^{-2}\Psi$
and $d_n^2b_n\in\mathbf Z$; the Apéry limit is $L(\Psi,2)$. (Measured, not assumed, in both
scans below.)

**(v) The singular points of $L_1$ are the special points of $\Gamma$ — and nothing else.**
This is the key structural point, and it is what makes the search finite.

> **Lemma 1.1.** The finite singular points of $L_1$ in the $t$-line are exactly the values
> $t(P)$ for $P$ a cusp or elliptic point of $\Gamma$, $t(P)\notin\{0,\infty\}$.
> They do **not** depend on $F$.

*Proof.* At a point $P\in\mathfrak H$ that is not elliptic, $t$ is a local isomorphism and
$g$ is holomorphic and non-vanishing, so the two solutions $g$, $g\int W$ of $L_1$ are
holomorphic and independent: $t(P)$ is an ordinary point. At an elliptic point of order
$e$, $t-t(P)$ vanishes to order $e$ in the uniformiser, producing exponent difference
$1/e$. At a cusp $c$ with $t(c)$ finite, the second solution carries $\log q_c$, i.e.
$\log(t-t(c))$: exponents $(0,0)$. Zeros or poles of $F$ at a cusp change the **exponents**
there, not the **location**. $\square$

Consequently $|t_1|\le|t_2|\le\cdots$ are determined by $(\Gamma,t)$ alone,
$\lambda_1=\lambda/|t_1|$, $\lambda_2=\lambda/|t_2|$, and

> **Theorem N1.** $\operatorname{score}=\log(1/\lambda_2)-k=\log|t_2|-2-\log\lambda$.

Two corollaries worth stating.

* **The brief's family (a) is empty.** $f=\sqrt{F\cdot U}$ with $U$ a modular unit of
  weight $0$, $U(0)=1$: $FU$ is again a weight-two form with cuspidal divisor, so by
  Lemma 1.1 the singular set — hence $|t_2|$, hence the archimedean half of the score — is
  **unchanged**. Multiplying by a unit can only change $\lambda$, and only through
  $e=\min v_2[t^n](FU)$. So (a) is a search over $\lambda$, already covered by scanning all
  weight-two eta quotients (the set of which is closed under multiplication by units).
* **Non-eta $F$ is strictly worse.** If $F$ has a zero in $\mathfrak H$ then $\sqrt F$
  branches there, adding a singular point inside the disc of convergence and shrinking
  $|t_2|$. (This is exactly hypothesis (H3) of Theorem R4, and exactly how the
  $\zeta(5)$ level-16 root row dies in `ROOT_ROWS.md` §4.)

## 2. What CDT does and does not give

The unbounded-denominators theorem (Calegari–Dimitrov–Tang, arXiv:2109.09040) is used here
in exactly one place and in exactly one direction.

* **What it gives.** A modular form on a finite-index subgroup with algebraic Fourier
  coefficients and **bounded** denominators is a form on a congruence subgroup. So
  $\lambda=1\Rightarrow$ congruence, i.e. **$\lambda>1$ is a certificate of
  non-congruence**. Their proof is entirely archimedean: an arithmetic holonomicity bound
  $$\dim_{\mathbf Q(p(x))}V(U,x(t),\mathbf Z)\ \le\ e\cdot
  \frac{\int_{\mathbf T}\log^+|p\circ\varphi|\,d\mu}{\log|\varphi'(0)|},\qquad|\varphi'(0)|>1,$$
  with the conformal radius made explicit by
  $|F_N'(0)|=16^{1/N}\Gamma(1+\tfrac1{2N})^2\Gamma(1-\tfrac1N)/
  \bigl(\Gamma(1-\tfrac1{2N})^2\Gamma(1+\tfrac1N)\bigr)
  =16^{1/N}(1+\zeta(3)/2N^3+\cdots)$, against a doubling of the excess
  $[R_{2N}:M_{2N}]$ at each new prime; the clash is
  $2^{(1-\varepsilon)X/\log X}$ against $O(X)$.
* **What it does not give.** *No $\lambda$.* CDT is qualitative: bounded vs unbounded.
  There is no level bound, no set of primes, no rate anywhere in the paper. The one
  explicit modulus in it, $M=c^2N^4M_1^2$, bounds the **inverse** of the uniformising
  coordinate and is discarded in the $\kappa\to0$ limit. CDT themselves flag the gap
  (§7.5.5), and note that finiteness of the deficient set comes from **Eisenstein's
  theorem** on algebraic power series, not from their method.
* **Where our $\lambda$ actually comes from.** The elementary binomial estimate
  $(1+u)^{1/N}\in\mathbf Z[[u/N^2]]$ — which CDT quote inside their Siegel lemma as "a
  simple denominator estimate". At $N=2$ this is $\lambda\mid4$, i.e. Theorem R1, and it
  is sharp ($(1-4x)^{1/2}=1-2\sum \mathrm C_{n-1}x^n$). **So Beukers' $4^n$ is binomial,
  not CDT.** CDT's role is only the converse implication $\lambda=1\Rightarrow$ congruence,
  which is what lets us call the winner a non-congruence row at all.
* **Sequel note.** The companion paper (arXiv:2408.15403) contains the inequality that our
  score is a shadow of: with $\tau=b\sigma$ the growth rate of $\operatorname{lcm}$
  denominators, the entry condition is $|\varphi'(0)|>e^{\tau}$ and the bound is
  $m\le e\int\log^+|\varphi|/(\log|\varphi'(0)|-\tau)$. Our
  $\log(1/\lambda_2)>k$ is literally "conformal size beats $e^{\text{denominator rate}}$"
  for a single function; CDT's version applies to the whole extension at once. That is the
  route by which a *negative*-score row could still be attacked, and it is the one
  `THEORY_NOTES_06` §4(c) recommends.

## 3. Four special points: the host groups are a finite list

A **three-term** recurrence with quadratic coefficients means the operator $L_1$ has exactly
four singular points on $\mathbf P^1_t$: the MUM point $t=0$, two finite points $t_1,t_2$,
and $t=\infty$. By Lemma 1.1 these are values of $t$ at special points of $\Gamma$, and $t$
takes each value once (it is a Hauptmodul), so $\Gamma$ has exactly **four** special points:
$r$ elliptic of orders $e_1,\dots,e_r$ and $s$ cusps, $r+s=4$, $s\ge1$ (the MUM point is a
cusp). Gauss–Bonnet at genus zero gives the covolume index in $\mathrm{PSL}_2(\mathbf Z)$
$$\mu=6\Bigl(2g-2+s+\sum_i(1-1/e_i)\Bigr)=6\Bigl(2-\sum_{i=1}^r\frac1{e_i}\Bigr)\ \le\ 12 ,$$
with $\mu=12$ exactly in the four-cusp case. So **the hosts form a finite list**.
The possible $(e_2$-count, $e_3$-count, cusps$)$ and indices inside $\mathrm{PSL}_2(\mathbf Z)$
(where elliptic orders are $2$ and $3$ only) are

| $(e_2,e_3,\text{cusps})$ | $(3,0,1)$ | $(2,1,1)$ | $(1,2,1)$ | $(0,3,1)$ | $(2,0,2)$ | $(1,1,2)$ | $(0,2,2)$ | $(1,0,3)$ | $(0,1,3)$ | $(0,0,4)$ |
|---|---|---|---|---|---|---|---|---|---|---|
| index $\mu$ | 3 | 4 | 5 | 6 | 6 | **7** | 8 | **9** | **10** | 12 |

`06_groups.py` enumerates all of them exactly (subgroups of index $n$ of
$\mathrm{PSL}_2(\mathbf Z)=\langle S,L\mid S^2=L^3=1\rangle$ $\leftrightarrow$ transitive
pairs $(s,t)\in S_n$ with $s^2=t^3=1$, canonicalised by BFS relabelling; cusps $=$ cycles of
$st$; congruence tested by Wohlfahrt — $N=\operatorname{lcm}$ of cusp widths, then does the
coset action factor through $\mathrm{PSL}_2(\mathbf Z/N)$, decided by a graph-subgroup order
comparison in `sympy`).

> **Result.** Up to conjugacy there are exactly **28** genus-zero subgroups of
> $\mathrm{PSL}_2(\mathbf Z)$ with four special points: **16 congruence**, **12
> non-congruence**. The non-congruence ones are
>
> | index | signature | cusp widths |
> |---|---|---|
> | 7 | $(0;2,3;2)$ | $[3,4]$, $[2,5]$, $[1,6]$, $[1,6]$ |
> | 9 | $(0;2;3)$ | $[1,3,5]$, $[1,2,6]$, $[1,1,7]$, $[2,3,4]$ |
> | 10 | $(0;3;3)$ | $[1,4,5]$, $[1,2,7]$, $[1,1,8]$, $[2,3,5]$ |
>
> **Check:** the index-$12$ four-cusp row of the table returns exactly six groups, widths
> $[1,1,5,5]$, $[1,2,3,6]$, $[1,1,2,8]$, $[1,1,1,9]$, $[2,2,4,4]$, $[3,3,3,3]$, **all
> congruence** — Beauville's six families, i.e. Zagier's six sporadic weight-one rows,
> recovered from scratch. So *no* non-congruence group has four cusps: every
> non-congruence host has an elliptic point.

**Theorem N3 (unequal exponents kill the score).** The exponent difference of $L_1$ at
$t_i$ is $1/e_i$ ($e_i=\infty$ at a cusp). Write $L_1=P\theta^2+Q\theta+R$, $P=1-\alpha
t+\delta t^2$, $P(0)=1$. Then $Q(t_i)=(1-1/e_i)\,t_iP'(t_i)$, two linear conditions on $Q$.
If $t_1,t_2$ are Galois-conjugate irrationals, applying the nontrivial automorphism carries
the condition at $t_1$ to the *same* condition at $t_2$, forcing $e_1=e_2$. Hence
$$e_1\ne e_2\ \Longrightarrow\ t_1,t_2\in\mathbf Q\ \Longrightarrow\
\lambda_1,\lambda_2\in\mathbf Q,$$
and if moreover $P,Q,R\in\mathbf Z[t]$ (the leading recurrence coefficient is exactly
$(n+1)^2$ — true for every row in the project's census) then $\lambda_{1,2}$ are roots of a
monic integer quadratic, hence **integers**; a positive score needs $|\lambda_2|<e^{-2}$,
so $\lambda_2=0$ and the recurrence degenerates. $\square$

**Consequence for the brief's questions (b), (d).** Read the table against N3:

* **index 7**, $(0;2,3;2\text{ cusps})$. The MUM point is a cusp; the remaining three
  special points are $\{$cusp, $e{=}2$, $e{=}3\}$, and $t$ sends one of them to $\infty$.
  Every choice leaves a **mixed** pair $\{2,3\}$, $\{2,\text{cusp}\}$ or
  $\{3,\text{cusp}\}$ at $t_1,t_2$. By N3 the score is negative. **All four index-7
  non-congruence groups are excluded outright.**
* **index 9** $(0;2;3)$ and **index 10** $(0;3;3)$. The only non-mixed option is to send the
  elliptic point to $t=\infty$, leaving two **cusps** at $t_1,t_2$: exponent differences
  $0,0$, i.e. the **Zagier class** $\beta=\alpha$, $\varepsilon=0$. That class is scanned
  exhaustively in §4 ($e=\infty$).
* So: **a non-congruence *host group* contributes nothing beyond what the $e=\infty$
  recurrence scan already covers.** The interesting non-congruence-ness is not in the group
  of $t$ but in the multiplier system of $g=\sqrt F$ — which is exactly Beukers' situation
  (his $t$ lives on the *congruence* group $\Gamma_0(6)+6$; his $\sqrt E$ does not).
* **Quadratic pullbacks are always fatal.** If $t=t'^2$ then each singular value $t_i$ has
  two preimages $\pm\sqrt{t_i}$ of equal modulus, so
  $\lambda_1'=\lambda_2'=\lambda_1^{1/2}>1$ and the score is
  $-\tfrac12\log\lambda_1-k<0$. (This is the brief's "$t'=\sqrt t$, NOT what we want",
  now with a one-line proof.)

## 4. Scan A: the recurrence classes, ansatz-free

The construction of §1 is *not* assumed here. By §3 any three-term second-order MUM row with
integer coefficients and leading term $(n+1)^2$ lies in
$$(n+1)^2a_{n+1}=(\alpha n^2+\beta n+\gamma)a_n-(\delta n^2+\varepsilon n+\zeta)a_{n-1},
\qquad a_0=1,\ a_1=\gamma,$$
with $\lambda_{1,2}$ the roots of $x^2-\alpha x+\delta$; by N3 a positive score forces
$e_1=e_2=e$, and then
$$\beta=\frac{e-1}{e}\,\alpha,\qquad \varepsilon=-\frac{2}{e}\,\delta,\qquad
\zeta\ \text{free}\quad(e=\infty:\ \beta=\alpha,\ \varepsilon=0).$$
$e=2$ is exactly the R3 normalisation $Q=\tfrac12\theta P$ of `ROOT_ROWS.md` §3 — the
free-integration class, i.e. the class of the square-root rows; $\zeta=\delta/4$ recovers the
three-parameter "sqrt class" $(A,B,C)$ of `SQRT_APERY.md` Thm 2. $e=\infty$ is the class of
Zagier's six (which further impose $\zeta=0$; here $\zeta$ is free).

`01_class_scan.c` enumerates integer $(\alpha,\gamma,\delta,\zeta)$ and tests integrality
division-free via $A_n=a_n(n!)^2$, $A_{n+1}=P(n)A_n-n^2Q(n)A_{n-1}$, prime by prime modulo
$p^K$ with $K>2v_p(n!)$ — no big integers, $\approx3\times10^7$ candidates/s/core.
`04_score.py` then applies the two filters that matter:

* **Casoratian.** $W_n=a_nb_{n+1}-a_{n+1}b_n$ satisfies $W_n=Q(n)W_{n-1}/(n+1)^2$, so if
  $Q(n)=0$ for some integer $n\ge1$ then $W_m\equiv0$ for $m\ge n$: the companion is a
  rational multiple of the row and the criterion is vacuous. **This filter alone removes
  $98\%$ of the apparent winners** — every "score $>3$" family found by the raw scan is of
  this degenerate type (e.g. $\gamma=-\alpha/2$, $\zeta=-2\delta$, where
  $P(n)=\tfrac\alpha2(2n-1)(n+1)$ and $Q(n)=\delta(n-2)(n+1)$).
* **Measured $k$.** $k=\min\{k: d_n^kb_n\in\mathbf Z\}$ is *measured*, not assumed.

**Boxes searched** ($N=26$ terms of exact integrality per prime, then re-verified to $n=40$
and, for the finalists, to $n=240$):

| class | $\alpha$ | $|\delta|$ | $|\gamma|$ | $|\zeta|$ |
|---|---|---|---|---|
| $e=2$ (R3 / square-root class) | $\le3000$ | $\le150$ | $\le150$ | $\le300$ |
| $e=\infty,3,4,6$ | $\le3000$ | $\le150$ | $\le150$ | $\le300$ |

plus the wider boxes of §4.3.

### 4.1 Every positive-score row with $k\ge2$

There are exactly **fourteen**, and they are all classical:

| score | class | $(\alpha,\gamma,\delta,\zeta)$ | $\lambda_1$ | $\lambda_2$ | $k$ | period $\xi$ |
|---|---|---|---|---|---|---|
| $+0.52533$ | $\infty$ | $(200,-150,16,-36)$ | $199.9200$ | $0.0800320$ | 2 | $-\tfrac1{12}\log\tfrac{13}{12}$ |
| $+0.52533$ | $\infty$ | $(200,-150,16,-100)$ | $199.9200$ | $0.0800320$ | 2 | $\mathbf Q$-comb. of the same log |
| $+0.44187$ | $\infty$ | $(184,-138,16,-36)$ | $183.9130$ | $0.0869977$ | 2 | $-\tfrac1{12}\log\tfrac{12}{11}$ |
| $+0.44187$ | $\infty$ | $(184,-138,16,-100)$ | $183.9130$ | $0.0869977$ | 2 | same log |
| $\mathbf{+0.40606}$ | $\infty$ | $(11,3,-1,0)$ | $\tfrac{11+5\sqrt5}2$ | $\tfrac{11-5\sqrt5}2$ | 2 | $\boldsymbol{\zeta(2)/5}$ — **Apéry** |
| $+0.35081$ | $\infty$ | $(168,-126,16,-36/-100)$ | $167.9047$ | $0.0952921$ | 2 | $-\tfrac1{12}\log\tfrac{11}{10}$ |
| $+0.25060$ | $\infty$ | $(152,-114,16,-36/-100)$ | $151.8947$ | $0.1053362$ | 2 | $-\tfrac1{12}\log\tfrac{10}{9}$ |
| $\mathbf{+0.13920}$ | $\mathbf 2$ | $\mathbf{(136,10,16,4)}$ | $4(17{+}12\sqrt2)$ | $4(17{+}12\sqrt2)^{-1}$ | 2 | **Beukers 1987 Thm 3** |
| $+0.13920$ | $\infty$ | $(136,-102,16,-36/-100)$ | $135.8823$ | $0.1177490$ | 2 | $-\tfrac1{12}\log\tfrac9{8}$ |
| $+0.01379$ | $\infty$ | $(120,-90,16,-36/-100)$ | $119.8665$ | $0.1334818$ | 2 | $-\tfrac1{12}\log\tfrac87$ |

All twelve $e=\infty$ entries with $\gamma=-\tfrac34\alpha$, $\delta=16$ are the classical
**Legendre/Padé** family: $\alpha=8x$ with $x$ odd, $\lambda_{1,2}=4(x\mp\sqrt{x^2-1})$,
$\xi=-\tfrac1{12}\log\frac{x+1}{x-1}$ up to a rational shift — identified by
$\exp(12\xi)\in\mathbf Q$ at $130$ digits (`08_period.py`). The two genuinely modular
entries are Apéry's $\zeta(2)$ row and **Beukers' Theorem 3**, and the latter is the only
one in the $e=2$ (square-root, non-congruence) class. Integrality of $a_n$ and $k=2$ were
re-verified exactly to $n=240$ for all fourteen.

### 4.2 Positive score with $k\le1$

$3310$ further rows have positive score with $k=1$: $213$ in class $e=\infty$, $877$ in
$e=3$, $1822$ in $e=4$, $398$ in $e=6$. They cannot come from the §1 construction —
Theorem R3 forces $k=2$ there — and the ones we identified are again Legendre rows:
e.g. $(600,150,16,-4)$ is $a_n=\binom{2n}nP_n(75)$ (the recurrence factors as
$(n+1)a_{n+1}=75(2n+1)\cdot2a_n-\ldots$, i.e. Legendre's recurrence at $x=75$) with
$\exp(4\xi)=\tfrac{38}{37}$. Three representatives from $e=3,4,6$ resisted a
$\log$/algebraic battery ($\exp(k\xi)$ rational for $k\le36$; `algdep` to degree $6$); they
are Padé-type limits close to $1/\gamma$ and **we make no claim about them** — they are
outside the question this note asks.

## 5. The runner-up: a second non-congruence Beukers row, at level 5

The best negative-score row in the $e=2$ class that is **not** a rescaling of a known one is
$$(\alpha,\gamma,\delta,\zeta)=(88,6,-64,-12):\qquad
(n+1)^2a_{n+1}=2(44n^2+22n+3)\,a_n+4(4n-1)(4n-3)\,a_{n-1},\quad a_0=1,\ a_1=6 .$$
It is **not** one of the nine $\operatorname{Sym}^1$ square roots of `SQRT_APERY.md` §7
($\zeta\ne\delta/4$), and it is primitive ($(44,3,-16,-3)$ is not integral). We identify it
(`09_level5.gp`):

> **The level-$5$ Fricke row.** Let $u=(\eta_5/\eta_1)^6$ (a Hauptmodul of $\Gamma_0(5)$,
> $u\circ W_5=1/(125u)$), and
> $$t=\frac{u}{1+22u+125u^2}=q-16q^2+122q^3-568q^4+1555q^5-\cdots$$
> the resulting Hauptmodul of the Fricke group $\Gamma_0(5)+5$; let
> $F=E_{2,5}=\tfrac14\bigl(5E_2(5\tau)-E_2(\tau)\bigr)=1+6q+18q^2+24q^3+\cdots$
> be the weight-two Eisenstein series of $\Gamma_0(5)$ ($F|_2W_5=-F$, divisor $=$ the cusp
> $0$). Then $g=\sqrt F$ is weight one on an index-two subgroup and
> $$a_n=2^n[t^n]\sqrt F=1,\,6,\,210,\,10500,\,615510,\,39474036,\,2682201396,\dots$$
> \[verified integral for $n\le202$; $\lambda=2$ is **minimal** — $[t^n]\sqrt F\notin\mathbf Z$,
> $v_2(\mathrm{den})=0,1,1,3,3,4,4,7,7,8,8,10,\dots$\]. The minimal recurrence is the one
> above, with characteristic polynomial $-L^2+88L+64$,
> $$\lambda_1=44+20\sqrt5=88.72135955\ldots,\qquad \lambda_2=44-20\sqrt5=-0.72135955\ldots,$$
> $k=2$ **sharp** ($d_nb_n\notin\mathbf Z$, $d_n^2b_n\in\mathbf Z$ for $n\le200$),
> Casoratian $W_n=1,-3,\tfrac{140}3,-1155,36036,\dots\ne0$, and
> $$\xi=\lim\frac{b_n}{a_n}=0.16430670106434215862993215002712765277112781369140994089554832534\ldots$$
> \[$260$ digits computed; the Casoratian tail at $n=1200$ is $O(10^{-2425})$.\]

$\lambda=2>1$, so by Calegari–Dimitrov–Tang the weight-one $g$ and the weight-three source
$\Psi=g^3t/2$ live on a **non-congruence** group: this is a genuine second instance of
Beukers' Theorem-3 construction, at the other perfect level. Its host $\Gamma_0(5)+5$ has
covolume index $3$ and signature $(0;2,2,2;1\ \text{cusp})$ — one cusp (the MUM point) and
three elliptic points of order two, one of which $t$ sends to $\infty$ (whence the exponents
$(\tfrac14,\tfrac34)$ there, i.e. $\zeta/\delta=3/16$).

**But the score is negative:** $|t_2|=\lambda/|\lambda_2|=2.7726$, far below
$2e^2=14.778$, so
$$\operatorname{score}=\log\tfrac1{0.7213595}-2=-1.6734 .$$
No irrationality follows. $\xi$ resisted identification: two-term `lindep` at $120$ digits
against $1,\pi,\pi^2,\pi^3,\zeta(3),\log2,\log5,\log\varphi,\sqrt5$, and against $L(f,s)$,
$s=1,2,3$, for **every** weight-three newform of every level $M\le120$ with $5\mid M$ —
zero hits (`10_ident5.gp`). Exactly as for Beukers' $\xi$, this is what one expects of the
critical value of a non-congruence weight-three form.

**Two remarks.** (i) This row is invisible to any eta-quotient scan: neither $t$ nor $F$ is
an eta quotient (at level $5$ the only weight-two eta quotient with $F(0)=1$ is
$\eta_1^5/\eta_5$, which pairs with $u$ to give the Almkvist–Zudilin $(11,5,125)$ square
root, characteristic roots $-44\pm8i$, no real limit). It is the ansatz-free Scan A that
finds it. (ii) $\lambda_1=44+20\sqrt5=8\varphi^5\cdot\ldots$, $\lambda_1\lambda_2=-64$: the
Fricke pairing here is $t_1t_2=-\tfrac14$, not $t_1t_2=1$ as at level $6$. That factor of
$4$ is precisely why level $6$ wins and level $5$ loses.

## 6. The winner, and its full Apéry-ingredient status

The unique positive-score non-congruence row in every scan is

> $t=(\eta_1\eta_6/\eta_2\eta_3)^{12}$ on $\Gamma_0(6)+6$,
> $F=\eta_2^7\eta_3^7/(\eta_1^5\eta_6^5)$ (Apéry's weight-two form, $F|_2W_6=-F$,
> divisor $=(1/2)+(1/3)$), $g=\sqrt F$,
> $a_n=4^n[t^n]\sqrt F=1,10,534,40900,3672550,\dots$,
> $\lambda_1=4(17+12\sqrt2)$, $\lambda_2=4(17-12\sqrt2)$, $k=2$, $\lambda=4$,
> $\operatorname{score}=+0.13920$.

This is **Beukers, Astérisque 147–148 (1987), Theorem 3** — already recorded at the head of
`SQRT_APERY.md`. Its four Apéry ingredients are *proved* there, not merely measured
(`SQRT_APERY.md` Theorems 1–4): $a_n\in\mathbf Z$ (Theorem 1 with $\lambda=4$ sharp),
$d_n^2b_n\in\mathbf Z$ with $k=2$ sharp (Theorem 3, via $\Psi=f^3t/4$ and
$b_n=4^n[t^n](f\,\theta_q^{-2}\Psi)$), $r_n=a_n\xi-b_n>0$ for every $n$ (Casoratian
$W_n=\mathrm C_n^2$), and $\limsup r_n^{1/n}\le\lambda_2$ (Poincaré). Hence
$\xi=L(\Psi,2)=0.100187449229339406\ldots\notin\mathbf Q$, $\mu(\xi)\le50.654$. The scans
of this note add nothing to that statement and take nothing away; what they add is that
**it is the only one**.

Two archimedean remarks that put the margin in context.

* The criterion is $\lambda e^{k}<1/|t_1\lambda^{-1}|$, i.e. $4e^2<(\sqrt2+1)^4$:
  $29.5562<33.9706$, a margin of $1.1494$. Beukers writes the condition as
  $4e<(\sqrt2+1)^4$ (his denominator claim for the companion is one $\operatorname{lcm}$
  factor); with the $k=2$ that this project verifies sharply ($d_5b_5\notin\mathbf Z$) the
  condition is $4e^2<(\sqrt2+1)^4$ and it still holds — the theorem is safe either way, but
  the measure $\mu\le50.654$ uses $k=2$.
* Level $6$ wins because $t\circ W_6=t$ **with $t_1t_2=1$**: the Fricke pairing is exact, so
  $|t_2|=1/|t_1|=\lambda_1/\lambda=33.97$ is as large as the fastest-growing row allows.
  At level $5$ (§5) the pairing is $t_1t_2=-\tfrac14$ and the same $\lambda_1$-for-$|t_2|$
  exchange loses a factor $4$; that is the whole difference between the two levels.

## 7. Scan B: the modular cross-check, $N\le60$

`02_eta_enum.py` enumerates, level by level, via the (invertible) Ligozat matrix
$A_{c,d}=\frac{N}{24\gcd(c^2,N)}\frac{\gcd(c,d)^2}{d}$ acting on **cusp-order vectors**
rather than exponent vectors:

* parameters $t=\prod_{d\mid N}\eta(d\tau)^{r_d}$ of weight $0$ with
  $\operatorname{ord}_\infty t=1$ and divisor degree $\le4$ (so $t$ is a Hauptmodul of
  $\Gamma_0(N)$ or of an Atkin–Lehner quotient), $r\in\mathbf Z^{\tau(N)}$;
* forms $F=\prod\eta(d\tau)^{s_d}$ of weight $2$, **holomorphic**, with
  $\operatorname{ord}_\infty F=0$, i.e. $\sum_c\varphi(\gcd(c,N/c))\operatorname{ord}_cF=\mu/6$.

$N\le60$: **3 036 854 pairs**. `03b_exact.gp` computes, exactly over $\mathbf Q$,
$q(t)=\mathrm{serreverse}(t)$, $F(t)$, $g=\sqrt{F(t)}$, the minimal
$\lambda=2^{\lceil v_2/n\rceil}$, the row $a_n=\lambda^n[t^n]g$ to $n=34$, and the minimal
recurrence (order $\le4$, degree $\le3$) by a mod-$(2^{31}-1)$ gate followed by an exact
$\mathbf Q$ kernel; $\lambda_{1,2}$ from the leading symbol; $k$ measured for every row with
$|\lambda_2|<1.05$.

**Result.**

| | |
|---|---|
| pairs scanned | $3\,036\,854$ |
| pairs with a recurrence of order $\le4$, degree $\le3$ | $739$ |
| of those, $\lambda=1$ / $\lambda=2$ / $\lambda=4$ | $485$ / $127$ / $127$ |
| distinct invariant classes $(\lambda_1,\lambda_2,\lambda,\text{ord},\deg)$ | $104$ |
| levels represented | $4,5,6,8,9,10,12,14,15,16,18,20,24,27,28,30,32,36,40,48$ |
| **classes with $|\lambda_2|<1$** | **2** |

| score | $k$ | $\lambda$ | $N$ | $\lambda_1$ | $\lambda_2$ | $t$ | $F$ | $F$ exponents even? |
|---|---|---|---|---|---|---|---|---|
| $\mathbf{+0.1392}$ | 2 | **4** | 6 | $135.8823$ | $0.1177490$ | $(\eta_1\eta_6/\eta_2\eta_3)^{12}$ | $\eta_2^7\eta_3^7/(\eta_1^5\eta_6^5)$ | **no** $\Rightarrow$ non-congruence |
| $-1.6235$ | 2 | 1 | 8 | $12{+}8\sqrt2$ | $12{-}8\sqrt2$ | $(\eta_1\eta_8/\eta_2\eta_4)^{8}$ | $\eta_2^{16}/(\eta_1^8\eta_4^6)\cdot\ldots$ | yes $\Rightarrow$ congruence |

Everything else has $|\lambda_2|\ge1$: the linear form does not decay, whatever $k$ is. The
next classes by $|\lambda_2|$ are $\lambda_2=1$ exactly (levels $6$, $12$), then $1.0635$,
$1.0718$, $1.6180$ — i.e. the gap between "decays" and "does not" is populated by nothing.

**Two documented blind spots of Scan B** (which is why Scan A is the load-bearing one):
1. $F$ is required to be an eta quotient. Weight-two forms with cuspidal divisor need not be
   eta quotients — the level-$5$ row of §5 has $F=E_{2,5}$, which is not one.
2. $t$ is required to be an eta quotient. Fricke Hauptmoduls of $\Gamma_0(p)+p$ are rational
   functions of one ($t=u/(1+au+p^3u^2)$ at $p=5$), as is Cooper's $s_7$ parameter.
   Both blind spots are the same two already recorded in `SQRT_APERY.md` §7 correction 3.

## 8. Answers to the six sub-questions of the brief

| | question | answer |
|---|---|---|
| (a) | $f=\sqrt{F\cdot U}$, $U$ a modular unit | **Empty as a new family.** $FU$ is again weight two with cuspidal divisor, so by Lemma 1.1 the singular set — hence $|t_2|$ — is unchanged; only $\lambda$ can move, and the set of weight-two eta quotients is closed under multiplication by units, so Scan B already sweeps it. |
| (b) | non-congruence Hauptmoduls $t'$, not roots of $t$ | **Finitely many hosts, all excluded.** §3: only $12$ non-congruence genus-zero groups have four special points (index $7,9,10$). Index $7$ dies by Theorem N3 (mixed exponents $\Rightarrow$ rational singular points $\Rightarrow|\lambda_2|\ge1$); index $9,10$ reduce to the $e=\infty$ recurrence class, scanned exhaustively in §4 with no positive-score $k\ge2$ row except Apéry's $\zeta(2)$ and the classical $\log$ family. |
| (c) | all weight-two eta quotients $F$ on genus-zero $\Gamma_0(N)$, $N\le60$, with odd eta exponents | **Done (Scan B).** $3.0\times10^6$ pairs, $739$ rows, $104$ invariant classes, exactly $2$ with $|\lambda_2|<1$, exactly $1$ with positive score: Beukers'. $254/739$ of the rows are non-congruence ($\lambda\in\{2,4\}$); none of them decays except Beukers'. |
| (d) | index-two subgroups from $\sqrt{\text{unit}}$ | Same as (a)+(b): adjoining $\sqrt U$ gives a degree-two cover; the induced parameter is an algebraic function of $t$, and if it is the quadratic pullback $t=t'^2$ the dominant singularity acquires a partner of equal modulus and $\lambda_1'=\lambda_2'=\lambda_1^{1/2}>1$ (§3, last bullet). For a general degree-two cover the host is one of the finitely many groups of §3. |
| — | CDT's $\lambda$ | **CDT does not determine $\lambda$** (§2). $\lambda$ is the binomial/Eisenstein constant, here $\lambda\mid4$ (Theorem R1). CDT supplies only $\lambda=1\Rightarrow$ congruence. |
| — | is the winner already in Beukers 1987? | **Yes**, Theorem 3, and his §3 ($\Gamma_1(5)$) is the neighbouring level. Our §5 exhibits the $\Gamma_0(5)+5$ member of the same construction; it fails by $1.67$ nats. |

## 9. What this changes

1. **A finiteness theorem where there was a search.** `THEORY_NOTES_06` §4(b) proposed "a
   scan of genus-zero non-congruence groups … no one has done". Theorem N2 turns it into a
   *finite* problem: four special points $\Rightarrow$ covolume index $\le12$ $\Rightarrow$
   28 groups inside $\mathrm{PSL}_2(\mathbf Z)$, 12 of them non-congruence, all enumerated.
   Theorem N3 then kills every one of them that has unequal exponent differences.
2. **The score splits archimedean $\times$ arithmetic** (Theorem N1): $|t_2|$ is a property of
   $(\Gamma,t)$ alone and $\lambda$ of $F$ alone. That is why "level $6$ is perfect" is a
   statement about the Fricke pairing $t_1t_2=1$ and nothing else, and it explains
   quantitatively why level $5$ (pairing $t_1t_2=-\tfrac14$) loses by exactly $\log4$.
3. **A second explicit non-congruence Beukers row** (§5), on $\Gamma_0(5)+5$, with
   $\lambda=2$, $k=2$ sharp, and an unidentified period
   $0.16430670106434215863\ldots$ — the natural companion to
   `SQRT_APERY.md`'s $0.100187449229\ldots$, and equally resistant to identification.
   It is a new entry for the census (`ROW_LEDGER.md`).
4. **The Casoratian filter should be standard.** In the raw recurrence scan, *every* row
   scoring above $+1$ was a $Q(n_0)=0$ degeneracy with $W_n\equiv0$; without that filter one
   "finds" thousands of spurious winners. Likewise the measured-$k$ filter: the surviving
   positive scores are dominated by the classical Legendre/Padé $\log$ rows, which have
   $k=1$ and are not Beukers rows at all.
5. **Conjecture `conj:barrier` is untouched but its $w=2$ boundary is now sharp.** Among
   *all* integral second-order MUM rows in the boxes searched, exactly two have positive
   score and a modular (non-$\log$) period: Apéry's $\zeta(2)$ row ($+0.406$) and Beukers'
   Theorem 3 ($+0.139$). Both are congruence-hosted; only the second is non-congruence in
   the multiplier.
6. **No new irrationality theorem.** No positive-score row was found that is not already
   Apéry's, Beukers', or classical.

## 10. Status table

| statement | class |
|---|---|
| Lemma 1.1: singular points of $L_1$ $=$ special points of $\Gamma$, independent of $F$ | **proved** |
| Theorem N1: $\operatorname{score}=\log|t_2|-2-\log\lambda$ | **proved** (given R1, R3) |
| Theorem N2: four special points $\Rightarrow$ covolume index $\le12$; 28 groups, 12 non-congruence | **proved** (Gauss–Bonnet) + **computed** exhaustively (`06_groups.py`) |
| Beauville's six recovered as the index-12 four-cusp case, all congruence | **computed** (Wohlfahrt test) |
| Theorem N3: $e_1\ne e_2\Rightarrow t_1,t_2\in\mathbf Q\Rightarrow$ score $<0$ | **proved** *modulo* the normalisation $P,Q,R\in\mathbf Z[t]$, $P(0)=1$ (true in every census row) |
| quadratic pullback $t=t'^2$ always has $\lambda_1'=\lambda_2'>1$ | **proved** |
| Scan A: Apéry $\zeta(2)$, Beukers Thm 3, and the Legendre $\log$ family are the only positive-score rows with $k\ge2$ in the boxes | **computed**; integrality and $k$ re-verified exactly to $n=240$ |
| the level-5 $\Gamma_0(5)+5$ row: $\lambda=2$ minimal, $k=2$ sharp, $\operatorname{score}=-1.6734$ | **computed** ($a_n\in\mathbf Z$ to $n=202$, $k$ to $n=200$) |
| its period $0.164306701064\ldots$ | **open**; excluded against $L(f,s)$ for all weight-3 newforms of level $M\le120$, $5\mid M$, $s=1,2,3$, and a 9-constant battery |
| Scan B: only $2$ of $3.0\times10^6$ eta-quotient pairs give $|\lambda_2|<1$ | **computed** |
| Beukers' row is the unique positive-score non-congruence row | **computed within the stated boxes**; **proved** to be the only possible *host geometry* by N2+N3 modulo the box on $(\gamma,\zeta)$ |
| a positive-score row with $c>1$ (leading coefficient $c(n+1)^2$, $c\ne1$) | **not excluded**; see §11 |

## 11. Scope, honestly

* Scan A fixes the leading recurrence coefficient to be exactly $(n+1)^2$ (equivalently
  $P,Q,R\in\mathbf Z[t]$ with $P(0)=1$). A row whose primitive integral operator has
  $P(0)=c>1$ would have $\lambda_{1,2}$ with denominator dividing $c$ and is not covered;
  no such row occurs anywhere in the project's census, but it is not excluded a priori.
* Scan A's boxes are finite. The accessory parameters $(\gamma,\zeta)$ are unbounded in
  principle; the boxes reach $|\gamma|\le150$, $|\zeta|\le300$ at $\alpha\le3000$ and
  $|\gamma|\le1000$, $|\zeta|\le2000$ in the wider runs of §4.3.
* The group enumeration of §3 is complete for $\Gamma\subseteq\mathrm{PSL}_2(\mathbf Z)$. A
  genus-zero group *commensurable with but not inside* $\mathrm{PSL}_2(\mathbf Z)$ (Fricke
  and Atkin–Lehner quotients — including Beukers' own host $\Gamma_0(6)+6$ and the
  $\Gamma_0(5)+5$ of §5) is not in that list; such groups normalise congruence subgroups, so
  the *host* is congruence and the non-congruence-ness sits in the multiplier of $\sqrt F$.
  Scan A covers them anyway, being ansatz-free.
* Rows with $k\le1$ and positive score are *not* claims of new irrationality here. They are
  classical Padé/Legendre families; the three representatives we could not identify are
  recorded and left open.

## 12. Reproduction

```
lattice/noncongruence_scan/01_class_scan.c   # Scan A: integer recurrence classes
        #   ./01_class_scan CLASS ALMIN ALMAX DEMAX GAMAX ZEMAX NTEST
        #   CLASS = e (2,3,4,6), 0 = e=infinity, -1 = free exponent rho (beta free)
lattice/noncongruence_scan/04_score.py       # Casoratian filter, measured k, score, ranking
lattice/noncongruence_scan/11_summary.py     # merge + rank all classes
lattice/noncongruence_scan/06_groups.py      # the 28 four-special-point genus-0 subgroups
        #   of PSL_2(Z), with the Wohlfahrt congruence test (sympy)
lattice/noncongruence_scan/02_eta_enum.py    # Scan B: Ligozat enumeration of (t,F), N<=60
lattice/noncongruence_scan/mkgp.py mkjobs.py runjob.sh   # sharding
lattice/noncongruence_scan/03_modscan.gp     # eta products, fitmin, mod-p gated fitminp
lattice/noncongruence_scan/03b_exact.gp      # per-level scan: lambda, recurrence, k, score
lattice/noncongruence_scan/07_modtable.py    # Scan B aggregation and ranking
lattice/noncongruence_scan/08_period.py      # Apery limit by the exact Casoratian series
lattice/noncongruence_scan/09_level5.gp      # the Gamma_0(5)+5 row, verified
lattice/noncongruence_scan/10_ident5.gp      # its period to 260 digits + identification battery
lattice/noncongruence_scan/run_r3.sh run_e.sh run_free.sh run_more.sh snake_wide.sh
```
Raw output under `lattice/noncongruence_scan/out/`; the enumeration in `eta_pairs.json`
and `lev/L*.gp`.
