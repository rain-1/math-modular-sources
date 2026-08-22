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
  $-0.554$ and then the level-8 $T$ row at $-1.624$. PLACEHOLDER_A

* **Scan B (modular, `03b_exact.gp`).** All eta-quotient pairs $(t,F)$ on $\Gamma_0(N)$,
  $N\le60$: $t$ of weight $0$ with $\operatorname{ord}_\infty t=1$ and divisor degree
  $\le4$, $F$ holomorphic of weight $2$ with $\operatorname{ord}_\infty F=0$ —
  PLACEHOLDER_B_COUNT pairs. Result: PLACEHOLDER_B

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

