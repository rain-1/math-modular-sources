# Roots and free integrations: the general $w$-th root theorem and its census

*Claude (Opus 5), 2026-08-22. Scripts: `lattice/root_rows/`. Generalises
`SQRT_APERY.md` / `SQRT_APERY_FORMAL.md` (the $w=2$ case) to arbitrary $w$, and runs the
resulting construction over every $\operatorname{Sym}^w$ row in the project.*

---

## 0. Verdict first

Three theorems and one census.

* **Theorem R1** (integrality). For every $G\in1+x\mathbf Z[[x]]$ and every $w\ge2$,
  $$\boxed{\;\lambda_w^{\,n}\,[x^n]\,G^{1/w}\in\mathbf Z,\qquad \lambda_w=w\cdot\operatorname{rad}(w)=\prod_{p\mid w}p^{\,v_p(w)+1}\;}$$
  and $\lambda_w$ is **minimal**: $\lambda$ works for all $G$ iff $\lambda_w\mid\lambda$.
  ($w=2$: $\lambda=4$, recovering `SQRT_APERY_FORMAL` Thm 3.2. $w=3$: $9$; $w=4$: $8$;
  $w=5$: $25$; $w=6$: $36$; $w=8$: $16$.)
  Graded form: $\lambda(G)=\prod_{p\mid w}p^{\max(0,\,v_p(w)+1-e_p)}$ with
  $e_p=\min_{n\ge1}v_p([x^n]G)$, recovering Theorem $1'$ of `SQRT_APERY.md` at $w=2$.

* **Theorem R2** (descent). If $L_1=\tht^2+p\tht+r$ with $p,r\in x\mathbf Q[[x]]$ and
  $L_w=\operatorname{Sym}^wL_1$, then for $g\in1+x\mathbf Q[[x]]$
  $$L_w(g^w)=g^w\Bigl(w\,\tht^{w-1}\psi+\sum_{b<w-1}c_b\,\tht^b\psi+N(\psi)\Bigr),
  \qquad \psi:=\frac{L_1g}{g},$$
  with $c_b\in x\mathbf Q[[x]]$ and $N$ a polynomial all of whose monomials are of degree
  $\ge2$ in $(\psi,\tht\psi,\dots)$. Hence $L_w(g^w)=0\Rightarrow L_1g=0$. At $w=2$ the
  $N$-term is absent and this is exactly $L_{\rm Ap}(g^2)=2g\tht\omega+6(\tht g)\omega$.

* **Theorem R3** (free integration, $w$-independent). The companion of the root row costs
  $d_n^2$, not $d_n^{w+1}$: the second-order theory
  $L_1(g\,h(q))=\tht_q^2h/g^3$, $\Psi=g^3u$, $B=g\cdot\tht_q^{-2}\Psi$ applies verbatim,
  and $\Psi$ has weight $3$, so $\xi_{\rm root}=L(\Psi,2)$.
  **Verified, sharply ($d_nb_n\notin\mathbf Z$, $d_n^2b_n\in\mathbf Z$), on the $w=4$
  $\zeta(5)$ level-16 root row and the $w=6$ $\zeta(7)$ level-24 root row to $n=400$–$600$.**

* **Theorem R4** (score). Under hypotheses (H1)–(H4) of §4,
  $$\boxed{\;\operatorname{score}_{\rm root}=\operatorname{score}_{\rm parent}+(k_{\rm parent}-2)-\log\lambda\;}
  \qquad\text{generically}\quad =\operatorname{score}_{\rm parent}+(w-1)-\log\lambda .$$
  With the worst-case $\lambda=\lambda_w$ the gain $(w-1)-\log(w\operatorname{rad}w)$ is
  $-0.386$ ($w{=}2$), $-0.197$ ($w{=}3$), $\mathbf{+0.921}$ ($w{=}4$),
  $+0.781$ ($w{=}5$), $\mathbf{+1.417}$ ($w{=}6$), $+4.227$ ($w{=}8$).
  **So for $w\ge4$ the root is always a strict improvement** — the higher the weight, the
  bigger the prize.

**Census verdict.** The prize is not collected anywhere in this project. Every
$\operatorname{Sym}^w$ system with $w\ge3$ that the project owns fails one of the two
archimedean hypotheses of Theorem R4:

* $\beta(4)$ level $24$ and $L(4,\chi_{-3})$ level $24$ **are** clean $\operatorname{Sym}^3$
  systems, and their cube roots are clean second-order rows — but only in the *oriented*
  elliptic coordinate. There the $\beta(4)$ root row **is literally Zagier's Catalan row E**
  $(12,4,32)$, score $-3.386$, and the $L(4,\chi_{-3})$ root row is the hypergeometric term
  ${}_2F_1(\tfrac13,\tfrac13;1;-27t_0)$, whose recurrence has order $1$ — no second
  characteristic root, hence no Apéry companion at all.
* $\zeta(5)$ level $16$ and $\zeta(7)$ level $24$: the $w$-th root is integral with a
  surprisingly small $\lambda$ ($2$ and $3$), **is** second order as a differential operator,
  and has $k=2$ sharp — but every characteristic root is **double**, because taking the
  $w$-th root introduces fractional local exponents and clearing denominators squares the
  leading coefficient. A scalar conditional combination then cannot separate the
  two-dimensional dominant direction, and the linear form does **not** decay: measured,
  $b_n/a_n-\xi\asymp n^{-1/2}$ for $\zeta(7)$ (exactly the exponent difference
  $\tfrac7{12}-\tfrac1{12}$ carried by $(1-34s+s^2)^{1/12}$), and only polynomial
  convergence for $\zeta(5)$. For $\zeta(5)$ level $16$ hypothesis (H3) fails as well:
  a *new* dominant root $94.1902\ldots$ appears, coming from a **real zero of $F$** inside
  its disc of convergence.

**Apéry's square root remains the only positive-score root row in the project**
($+0.1392$; `SQRT_APERY.md` Theorem E, re-verified here). No new irrationality theorem.
Two by-products worth keeping: the identification $\beta(4)_{24}=\operatorname{Sym}^3(\text{Zagier E})$,
and the fact that measured $\lambda$ is far below the worst case in every case
($1,1,2,3,8,9$ against $9,9,8,36,8,36$).

---

## 1. Theorem R1: minimal $\lambda_w$

**Lemma 1.1.** Let $p$ be a prime and $a=v_p(w)$. Then for $k\ge1$
$$v_p\binom{1/w}{k}=\begin{cases}\;-ak-\dfrac{k-s_p(k)}{p-1}, & p\mid w,\\[2mm]
\;\ge 0,& p\nmid w,\end{cases}$$
where $s_p(k)$ is the digit sum of $k$ in base $p$.

*Proof.* $\binom{1/w}{k}=\frac1{w^kk!}\prod_{j=0}^{k-1}(1-jw)$. If $p\mid w$ each factor
$1-jw\equiv1\pmod{p^{a}}$ is a $p$-adic unit, so the product has $v_p=0$; and
$v_p(k!)=(k-s_p(k))/(p-1)$ (Legendre). If $p\nmid w$ then $1/w\in\mathbf Z_p$ and
$\binom{\cdot}{k}$ maps $\mathbf Z_p\to\mathbf Z_p$ (Mahler: $\binom{\cdot}{k}$ is
continuous and integer-valued on $\mathbf Z$, which is dense). $\square$
\[Verified exactly for all $2\le w\le12$, all $p\le13$, $k\le60$: `01_lambda_w.py`.\]

**Theorem R1.** Let $w\ge2$ and $\lambda\in\mathbf Z_{\ge1}$. Then
$\lambda^n[x^n]G^{1/w}\in\mathbf Z$ for every $G\in1+x\mathbf Z[[x]]$ and every $n\ge0$
**iff** $\lambda_w\mid\lambda$, where $\lambda_w=w\cdot\operatorname{rad}(w)$.

*Proof.* Write $G=1+S$, $S\in x\mathbf Z[[x]]$, $G^{1/w}=\sum_k\binom{1/w}kS^k$ (formal
binomial series; it is the unique $w$-th root with constant term $1$). Since
$[x^n]S^k=0$ for $k>n$,
$$\lambda^n[x^n]G^{1/w}=\sum_{k=0}^n\lambda^{\,n-k}\bigl(\lambda^k\tbinom{1/w}k\bigr)[x^n]S^k .$$
*Sufficiency.* For $p\nmid w$, Lemma 1.1 gives $v_p\ge0$ already. For $p\mid w$ with
$a=v_p(w)$ and $v_p(\lambda)\ge a+1$,
$$v_p\Bigl(\lambda^k\tbinom{1/w}k\Bigr)\ \ge\ k(a+1)-ak-\frac{k-s_p(k)}{p-1}
=k-\frac{k-s_p(k)}{p-1}\ \ge\ 0,$$
since $\frac{k-s_p(k)}{p-1}\le k-s_p(k)\le k$ for every $p\ge2$ (with equality throughout
only when $p=2$ and $s_2(k)=0$, which never happens for $k\ge1$; for $p=2$ the bound is
exactly $s_2(k)\ge1$). 
*Necessity.* Take $G=1+x$, so $[x^n]G^{1/w}=\binom{1/w}{n}$. If $v_p(\lambda)=v\le a$ then
$v_p\bigl(\lambda^n\binom{1/w}n\bigr)=n(v-a)-\frac{n-s_p(n)}{p-1}\to-\infty$
(take $n=p^j$: the value is $n(v-a)-\frac{p^j-1}{p-1}$). $\square$

Numerically $\lambda_w$: $w=2,3,4,5,6,7,8,9,10,11,12\mapsto
4,9,8,25,36,49,16,27,100,121,72$; sufficiency verified for $k\le120$ and minimality
(each $\lambda_w/p$ fails) for $k\le400$, plus a direct series test on random integral $G$
for $w\in\{2,3,4,5,6,8,9,12\}$ to $n=40$ (`01_lambda_w.py`).

**Theorem R1$'$ (graded).** With $e_p:=\min_{n\ge1}v_p([x^n]G)$ for $p\mid w$,
$$\lambda(G)=\prod_{p\mid w}p^{\max(0,\;v_p(w)+1-e_p)}$$
also works. *Proof.* Every coefficient of $S$ has $v_p\ge e_p$, so
$v_p([x^n]S^k)\ge ke_p$; run the same estimate with $v_p(\lambda)\ge a+1-e_p$ and
$n\ge k$. $\square$ At $w=2$ this is $\max(1,2^{2-e})$, i.e. Theorem $1'$ of
`SQRT_APERY.md`.

**This bound is sufficient, not sharp for a given $G$.** All six higher-weight rows of §5
have $\lambda$ strictly smaller than $\lambda(G)$; e.g. the $\beta(4)$ level-24 row has
$e_3=0$, hence $\lambda(G)=9$, but the true minimal $\lambda$ is $1$.

## 2. Theorem R2: the descent identity for every $w$

Work in $R=\mathbf Q[[x]]$, $\tht=x\,d/dx$. Let $L_1=\tht^2+p\tht+r$ with
$p,r\in xR$ (equivalently: $L_1=P\tht^2+Q\tht+R$ with $P(0)=1$, $Q(0)=R(0)=0$; exponents
$(0,0)$ at $x=0$). Let $L_w:=\operatorname{Sym}^wL_1$, the unique **monic** operator of
order $w+1$ annihilating $y^w$ for every solution $y$ of $L_1$.

Fix $g\in1+xR$ and put
$$u:=\frac{\tht g}{g}\in xR,\qquad \psi:=\frac{L_1g}{g}\in R,\qquad
Q_j:=\frac{\tht^j(g^w)}{g^w}.$$
Then $Q_0=1$, $Q_{j+1}=\tht Q_j+w\,u\,Q_j$, and $\tht u=\psi-u^2-pu-r$.

**Lemma 2.1 (weights).** Grade $\mathbf Q[u,p,r,\psi;\tht]$ by
$\operatorname{wt}(\tht^iu)=\operatorname{wt}(\tht^ip)=1+i$,
$\operatorname{wt}(\tht^ir)=\operatorname{wt}(\tht^i\psi)=2+i$. Then $Q_j$ is homogeneous
of weight $j$, and the coefficients $c_j$ of $L_w=\tht^{w+1}+\sum_{j\le w}c_j\tht^j$ are
homogeneous of weight $w+1-j$.
*Proof.* $\tht$ raises weight by $1$ (the substitution $\tht u=\psi-u^2-pu-r$ is homogeneous
of weight $2$), and $wuQ_j$ has weight $1+j$; induct. The $c_j$ are determined uniquely by
the $\psi=0$ specialisation $\sum_jc_jQ_j\big|_{\psi=0}=0$, a linear system with
homogeneous data. $\square$

**Theorem R2.** $\displaystyle \frac{L_w(g^w)}{g^w}
= w\,\tht^{w-1}\psi+\sum_{b=0}^{w-2}c^{(b)}\,\tht^b\psi+N(\psi,\tht\psi,\dots)$,
where each $c^{(b)}$ is a weighted-homogeneous polynomial in $u,p,r$ and their
$\tht$-derivatives of weight $w-1-b>0$ — hence $c^{(b)}\in xR$ — and every monomial of
$N$ has total degree $\ge2$ in the $\tht^i\psi$.

*Proof.* $\sum_jc_jQ_j$ vanishes identically at $\psi=0$ (definition of $L_w$), so every
monomial contains some $\tht^i\psi$. The highest $\psi$-derivative reachable is
$\tht^{w-1}\psi$, produced only along the chain $Q_2\ni w\psi\to\tht^{w-1}$ inside
$Q_{w+1}$, with coefficient $w$; the terms $c_jQ_j$ with $j\le w$ reach at most
$\tht^{w-2}\psi$. By Lemma 2.1 the coefficient of $\tht^b\psi$ has weight
$(w+1)-(b+2)=w-1-b$, which is $0$ only for $b=w-1$ and positive otherwise; a
positive-weight polynomial in $u,p,r,\dots\in xR$ lies in $xR$. $\square$

**Corollary 2.2 (descent).** If $g\in1+xR$ and $L_w(g^w)=0$ then $L_1g=0$.
*Proof.* $g^w$ is a unit, so the bracket vanishes. $\psi(0)=(L_1g)(0)=p(0)\tht g|_0+r(0)=0$,
so if $\psi\ne0$ its order $m$ satisfies $m\ge1$. Then $\tht^{w-1}\psi$ has order $m$ with
leading coefficient $m^{w-1}\ne0$; every $c^{(b)}\tht^b\psi$ has order $\ge m+1$; every
monomial of $N$ has order $\ge2m\ge m+1$. So $[x^m]$ of the bracket is $w\,m^{w-1}\cdot
[x^m]\psi\ne0$: contradiction. Hence $\psi=0$, i.e. $L_1g=0$. $\square$

At $w=2$ the statement is stronger — $N\equiv0$, and the identity is exactly linear:
$$\frac{L_2(g^2)}{g^2}=2\,\tht\psi+4(p+2u)\psi,\qquad\text{i.e.}\qquad
L_2(g^2)=2g\,\tht\omega+6(\tht g)\,\omega+4pg\,\omega,\quad\omega=L_1g,$$
which for $p=0$ (or, in the un-normalised form $L_1=P\tht^2+Q\tht+R$, whenever
$Q=\tfrac12\tht P$) is `SQRT_APERY_FORMAL` Lemma 2.3.

Explicitly (all verified symbolically, `02_descent.py`; the $\psi$-free part is exactly $0$
in each case):
$$\operatorname{Sym}^2 = \tht^3+3p\,\tht^2+(2p^2+\tht p+4r)\tht+2(2pr+\tht r),$$
$$\operatorname{Sym}^3 = \tht^4+6p\,\tht^3+(11p^2+4\tht p+10r)\tht^2
+(6p^3+7p\tht p+30pr+\tht^2p+10\tht r)\tht+3(6p^2r+5p\tht r+2r\tht p+3r^2+\tht^2r),$$
$$\frac{L_3(g^3)}{g^3}=3\tht^2\psi+15(p+2u)\tht\psi
+6\bigl(3p^2+10pu+\tht p-2r+10u^2\bigr)\psi+21\psi^2,$$
$$\frac{L_4(g^4)}{g^4}=4\tht^3\psi+36(p{+}2u)\tht^2\psi
+4\bigl(26p^2{+}90pu{+}7\tht p{-}14r{+}90u^2\bigr)\tht\psi+8\bigl(\cdots\bigr)\psi
+136\,\psi\bigl(2p\psi+4u\psi+\tht\psi\bigr).$$
$w=5$ is also printed by the script. The weight pattern
$\operatorname{wt}(c^{(b)})=w-1-b$ is visible in every line.

## 3. Theorem R3: the companion is $w$-independent — $k=2$ always

This is pure second-order theory, so it is inherited by the root row of **any**
$\operatorname{Sym}^w$ family.

Let $L_1=P\tht_u^2+Q\tht_u+R$ over $\mathbf Q[[u]]$ with $P(0)=1$, exponents $(0,0)$ at
$u=0$, and $g\in1+u\mathbf Z[[u]]$ the analytic solution. Let $W$ be the Wronskian-type
quantity $W=g^2\,\tht_u\log q$ where $q$ is the nome; $\tht_uW=-(Q/P)W$.

**Lemma 3.1.** For any $h\in\mathbf Q[[q]]$, viewed in $u$ via $q=q(u)$,
$$L_1(g\,h)=P\,W^2\,g^{-3}\,\tht_q^2h .$$
*Proof.* $L_1(gh)=Pg\tht_u^2h+(2P\tht_ug+Qg)\tht_uh$ since $L_1g=0$; substitute
$\tht_uh=\mu\tht_qh$, $\mu=\tht_u\log q=W/g^2$, $\tht_u^2h=\mu^2\tht_q^2h+(\tht_u\mu)\tht_qh$
and use $\tht_u\mu=-\mu(2\tht_ug/g+Q/P)$; the $\tht_qh$ bracket cancels identically. $\square$

**Normalisation.** If $Q=\tfrac12\tht_uP$ then $W=P^{-1/2}$ and $PW^2=1$, so
$$\boxed{\,L_1(g\,h)=g^{-3}\tht_q^2h\,}$$
— the identity of `SQRT_APERY_FORMAL` Lemma 5.3, with **no reference to $w$**.
(This condition is automatic for the operators here: $P=1-136u+16u^2$, $Q=-68u+16u^2$
for Apéry; $P=1-12t+32t^2$, $Q=-12t+64t^2$ for the $\beta(4)$/Catalan root; etc.)

**Theorem R3.** With $\Psi:=g^3u=\sum_{m\ge1}\psi(m)q^m$ and $\Theta:=\tht_q^{-2}\Psi$,
the census companion $B$ ($L_1B=u$, $b_0=0$, $b_1=1$) is $B=g\Theta$, and if
(i) $u(q)\in q\mathbf Z[[q]]$ and $q(u)\in u\mathbf Z[[u]]$ up to the scale $\lambda$, and
(ii) $\lambda^m\psi(m)\in\mathbf Z$, then $d_n^2b_n\in\mathbf Z$ — the proof of
`SQRT_APERY_FORMAL` §6 transposes verbatim. $\Psi=g^3u$ has weight $3$ ($g$ weight $1$,
$u$ weight $0$), so
$$\xi_{\rm root}=\lim_n\frac{b_n}{a_n}=L(\Psi,2)$$
whenever the Fricke/endpoint hypotheses of Theorem B$^*$ hold for $\Psi$.

**Measured, this run.** $k=2$ and $k=2$ *sharp* ($d_nb_n\notin\mathbf Z$) for:
all nine sporadic $w=2$ root rows ($n\le400$; re-verified), the $w=6$ $\zeta(7)$
level-24 root row in the $s$-coordinate ($n\le600$) and in the $z$-coordinate ($n\le400$),
and the $w=4$ $\zeta(5)$ level-16 root row ($n\le400$). **This is the first evidence that
the free integration $k:w+1\to2$ survives at $w=4$ and $w=6$.**

## 4. Theorem R4: the score, and exactly when it applies

Let the parent row be $F=\sum A_nt^n$ with $L_{\rm parent}=\operatorname{Sym}^wL_1$ over
$\mathbf Q(t)$, and let $g=F^{1/w}$, $a_n=\lambda^n[t^n]g$, i.e. the root row lives in
$u=t/\lambda$.

**Hypotheses.**
* **(H1)** $L_{\rm parent}=\operatorname{Sym}^wL_1$ with $L_1$ of order $2$ over
  $\mathbf Q(t)$, exponents $(0,0)$ at $t=0$.
* **(H2)** $\lambda^n[t^n]g\in\mathbf Z$ (Theorem R1/R1$'$, or better).
* **(H3)** $g$ has no zero and no extra branch point in $|t|\le1/|\lambda_2^{\rm parent}|$,
  so the finite singular points of $L_1$ in that disc are exactly those of $L_{\rm parent}$.
* **(H4)** the dominant singular point of $L_1$ is simple with exponents $(0,\rho)$,
  $\rho\notin\mathbf Z_{<0}$, so that the conditional combination $B-\xi g$ continues past it
  and $\limsup|a_n\xi-b_n|^{1/n}=\lambda\,|\lambda_2^{\rm parent}|$.

**Theorem R4.** Under (H1)–(H4): $\operatorname{Sym}^w$ does not move singular points, so
the characteristic roots of the root row are $\lambda\cdot$(those of the parent); in
particular $\lambda_1^{\rm root}=\lambda\lambda_1$, $\lambda_2^{\rm root}=\lambda\lambda_2$.
With $k_{\rm root}=2$ (Theorem R3),
$$\operatorname{score}_{\rm root}=\log\frac1{\lambda|\lambda_2|}-2
=\Bigl(\log\frac1{|\lambda_2|}-k_{\rm parent}\Bigr)+(k_{\rm parent}-2)-\log\lambda
=\operatorname{score}_{\rm parent}+(k_{\rm parent}-2)-\log\lambda,$$
$$\operatorname{budget}_{\rm root}=\operatorname{budget}_{\rm parent}+(k_{\rm parent}-2)+\log\lambda .$$
Generically $k_{\rm parent}=w+1$, so the gain is $(w-1)-\log\lambda$.

**Checks.** Apéry: $0.5255+(3-2)-\log4=+0.1392$ ✓. Zagier E as the cube root of
$\beta(4)_{24}$ in the oriented coordinate: parent
$\operatorname{score}=\log\frac14-4=-5.386$, root $=-5.386+2-0=-3.386$ ✓ (Zagier E's
tabulated score). $\zeta(7)$ level $24$ in $s$: parent
$\log\frac1{17-12\sqrt2}-7=-3.4745$, root $=-3.4745+5-\log9=-0.6717$, and indeed the fitted
root recurrence has $\lambda_2=9(17-12\sqrt2)=0.264935$ ✓ (though (H4) fails, see below).

**Best-possible gain.** With $\lambda=\lambda_w=w\operatorname{rad}(w)$:
| $w$ | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| $\lambda_w$ | 4 | 9 | 8 | 25 | 36 | 49 | 16 | 27 | 72 |
| gain $(w-1)-\log\lambda_w$ | $-0.386$ | $-0.197$ | $\mathbf{+0.921}$ | $+0.781$ | $\mathbf{+1.417}$ | $+2.108$ | $+4.227$ | $+4.704$ | $+6.723$ |

So **the root is a strict improvement for every $w\ge4$**, and the improvement grows with
$w$: a $\zeta(2k+1)$ system of weight $w=2k$ pays $\log(2k\operatorname{rad}2k)$ and gains
$2k-1$ units of denominator. This is the sharpest form of Problem 7.6 ("free integration")
the project has.

**Why (H3) and (H4) matter — both fail in practice.**
* (H3) fails for the $\zeta(5)$ level-16 row. $F$ has a **real zero inside its disc of
  convergence**, at
  $$t_0=0.021233630900370578190741143414062588154\ldots,\qquad
  1/t_0=47.0951013838404614321710604062754995511\ldots,$$
  the largest root of the irreducible quartic $x^4-38x^3-396x^2-1480x-2016$ (all four of its
  roots $47.0951$, $-3.1381$, $-2.9785\pm2.1839i$ are zeros of $F$, and all four appear,
  doubled, among the root row's characteristic roots). $F^{1/4}$ therefore has a branch point
  at $t_0$, well inside the parent's radius of convergence $1/(2+4\sqrt2)=0.1306$, and the
  root row acquires a *new* dominant characteristic root
  $$\lambda_1^{\rm root}=2/t_0=94.19020276768092286\ldots$$
  (root of $x^4-76x^3-1584x^2-11840x-32256$), demoting the parent's dominant root to
  $\lambda_2^{\rm root}=2(2+4\sqrt2)=4+8\sqrt2=15.3137>1$. The remaining characteristic roots
  $2\cdot4=8$, $2\cdot2=4$ and $8\sqrt2-4=7.3137$ are $\lambda$ times the parent's
  $-4,-2,2-4\sqrt2$, as Theorem R4 predicts for the singularities the parent does have.
* (H4) fails whenever the $w$-th root involves an algebraic gauge factor. For
  $\zeta(7)$ level $24$, $F=\sqrt{1-34s+s^2}\,\mathcal A(s)^3$, so
  $g=F^{1/6}=(1-34s+s^2)^{1/12}\sqrt{\mathcal A}$ has exponents $(\tfrac1{12},\tfrac7{12})$
  at each finite singular point instead of $(0,\tfrac12)$. Clearing denominators squares
  the leading coefficient: the fitted recurrence has leading coefficient
  $(1-306s+81s^2)^2$ and **every characteristic root is double**. A scalar $\xi$ cannot
  kill a two-dimensional dominant direction, and indeed
  $$b_n/a_n-\xi\ \asymp\ n^{-1/2}=n^{-(7/12-1/12)}$$
  is measured (successive-difference exponent $\to-1.5$ at $n=600$), while
  $\tfrac1n\log|a_n\xi-b_n|\uparrow\log\lambda_1$. The linear form grows.

## 5. The census

### 5.1 The nine sporadic $w=2$ rows (re-verified from scratch)

`03_sporadic.gp`: $A_n$ built from the order-3 recurrences, $g=\sqrt F$, minimal $\lambda$
measured, $a_n$ integral to $n=420$, minimal recurrence fitted and verified exactly on all
$n\le420$, characteristic roots to 250 digits. **Every entry of `SQRT_APERY.md` §7 and
`SPORADIC_SCAN2.md` §§7–8 is reproduced exactly**, including $\lambda\in\{1,2,4\}$ and the
graded prediction $\lambda=\max(1,2^{2-e_2})$.

| parent $(a,b,c)$ | $e_2$ | $\lambda$ | root recurrence $(n{+}1)^2a_{n+1}=p_1a_n-p_0a_{n-1}$ | $\lambda_1$ | $\lambda_2$ | $k$ | score | budget | $\xi_{\rm root}$ |
|---|---|---|---|---|---|---|---|---|---|
| **Apéry $(17,5,1)$** | 0 | **4** | $(136,10;4(2n{-}1)^2)$ | $4(17{+}12\sqrt2)=135.882$ | $4(17{-}12\sqrt2)=0.117749$ | 2✓ | $\mathbf{+0.1392}$ | $+2.9118$ | $L(\Psi,2)=0.1001874492\ldots$ |
| $T\,(12,4,16)$ | 2 | 1 | $(24,2;4(2n{-}1)^2)$ | $12{+}8\sqrt2=23.3137$ | $12{-}8\sqrt2=0.686292$ | 2✓ | $-1.6235$ | $+1.1490$ | $L(f_8,2)$ |
| Domb $(10,4,64)$ | 2 | 1 | $(20,2;16(2n{-}1)^2)$ | $16$ | $4$ | 2✓ | $-3.3863$ | $+0.7726$ | $L(f_{12},2)$ |
| AZ $(9,3,-27)$ | 0 | 4 | $(72,6;-108(2n{-}1)^2)$ | $77.5692$ | $-5.5692$ | 2✓ | $-3.7173$ | $+2.3512$ | $0.14551448\ldots$ unid. |
| AZ $(11,5,125)$ | 0 | 4 | $(88,10;500(2n{-}1)^2)$ | $44\pm8i$, $|\cdot|=44.721$ | — | 2✓ | n/a | $+1.8005$ | none |
| AZ $(7,3,81)$ | 0 | 4 | $(56,6;324(2n{-}1)^2)$ | $28\pm22.627i$, $|\cdot|=36$ | — | 2✓ | n/a | $+1.5835$ | none |
| Cooper $s_7$ | 2 | 1 | $p_1=26n^2{+}13n{+}2$, $p_0=-3(3n{-}1)(3n{-}2)$ | $27$ | $-1$ | 2✓ | $-2.0000$ | $+1.2958$ | $L(f_7,2)$ |
| Cooper $s_{10}$ | 1 | 2 | $p_1=24n^2{+}12n{+}2$, $p_0=-4(8n{-}3)(8n{-}5)$ | $32$ | $-8$ | 2✓ | $-4.0794$ | $+1.4657$ | $0.31692536\ldots$ unid. |
| Cooper $s_{18}$ | 1 | 2 | $p_1=56n^2{+}28n{+}6$, $p_0=12(8n{-}3)(8n{-}5)$ | $32$ | $24$ | 2✓ | $-5.1781$ | $+1.4657$ | $0.48423755\ldots$ unid. |

### 5.2 The $w\ge3$ systems

`04_build.gp` builds $400$ terms of every $A$-row; `05–08` find the minimal $\lambda$ and the
minimal recurrence $(\text{order},\deg_n)$ by a mod-$(2^{61}{-}1)$ kernel scan, verified
exactly over $\mathbf Q$ on all $n\le400$; `10_oriented.gp` does the oriented coordinates;
`11_companion.gp` builds the companion from the order-$2$ **differential** operator and
tests $d_n^kb_n$.

*(In the "recurrence" column $(\rho,\delta)$ means order $\rho$ in $n$, degree $\delta$ in
$n$; $\delta$ **is** the order of the differential operator, $\rho$ its $t$-degree.)*

| system | $w$ | $k_{\rm par}$ | coordinate | $\lambda$ (measured) | $\lambda$ (Thm R1$'$) | root recurrence | $\lambda_1$ | $\lambda_2$ | $k_{\rm root}$ | score |
|---|---|---|---|---|---|---|---|---|---|---|
| $\beta(4)$ level 24 | 3 | 4 | oriented $t=\eta_1^4\eta_4^2\eta_8^4/\eta_2^{10}$ | **1** | 9 | $(2,2)$ | $8$ | $4$ | 2 | $-3.386$ |
| $\beta(4)$ level 24 | 3 | 4 | rational $z=r/(1{+}r)$ | **1** | 9 | none with $\rho\le26,\delta\le8$ | — | — | — | direct image |
| $L(4,\chi_{-3})$ level 24 | 3 | 4 | oriented $t_0=(\eta_3/\eta_1)^{12}$ | **1** | 1 | $(1,2)$ | $-27$ | none | — | degenerate |
| $L(4,\chi_{-3})$ level 24 | 3 | 4 | rational $x$ (level 24) | **1** | 1 | $\ge(19,6)$ | — | — | — | direct image |
| $\zeta(5)$ level 16 | 4 | 5 | $t=x/(8x^2{+}2x{+}1)$ | **2** | 2 | $(14,\mathbf 2)$ | $94.1902$ (double) | $4{+}8\sqrt2=15.3137$ | **2**✓ | n/a |
| $\zeta(5)$ level 12 $(+)$ | 4 | 5 | Domb $t=w/(1{+}w)^2$ | **8** | 8 | $\ge(15,6)$ | — | — | — | not order 2 |
| $\zeta(5)$ level 12 $(-)$ | 4 | 5 | Domb $t=w/(1{+}w)^2$ | **8** | 8 | $\ge(15,6)$ | — | — | — | not order 2 |
| $\zeta(7)$ level 24 | 6 | 7 | $z=r/(1{+}r)$ | **3** | 36 | $(10,\mathbf 2)$ | $3(4{+}2\sqrt2)=20.485$ (dbl) | $3(2{+}2\sqrt2)=14.485$ | **2**✓ | n/a |
| $\zeta(7)$ level 24 | 6 | 7 | $s=(z/(1{-}z))^2$ | **9** | 36 | $(4,\mathbf 2)$ | $9(17{+}12\sqrt2)=305.735$ (dbl) | $9(17{-}12\sqrt2)=0.264935$ | **2**✓ | $(-0.672)$† |

† the formula's value; **not attained**, because $\lambda_1$ is a double root — see §4 (H4).

**Positive score:** Apéry's square root only. **$\lambda=1$:** $T$, Domb, Cooper $s_7$
(sporadic); $\beta(4)_{24}$ and $L(4,\chi_{-3})_{24}$ in both coordinates.

### 5.3 Two structural identifications (new)

**(a) $\beta(4)$ at level $24$ is the symmetric cube of the Catalan row.**
With $t=\eta_1^4\eta_4^2\eta_8^4/\eta_2^{10}$ and $F=\eta_2^{10}/(\eta_1^4\eta_4^4)=\theta_3(q)^2$
(verified as a $q$-identity to $q^{300}$), the cube root of the $\beta(4)$ Sym$^3$ row is
$$[t^n]F=1,\,4,\,20,\,112,\,676,\,4304,\,28496,\,194240,\,1353508,\dots$$
which is **exactly Zagier's row E $(12,4,32)$** — verified as an equality of integer
sequences for $n\le300$. Its operator is $(1-12t+32t^2)\tht^2+(-12t+64t^2)\tht+(-4t+32t^2)$
(note $Q=\tfrac12\tht P$, so the clean companion identity of §3 applies), $\lambda_1=8$,
$\lambda_2=4$, $k=2$, period $G/2$. The Sym$^3$ parent in $t$ is
$[t^n]F^3=1,12,108,880,6876,52752,\dots\in\mathbf Z$.

**(b) $L(4,\chi_{-3})$ at level $24$ is the symmetric cube of a hypergeometric term.**
With $t_0=(\eta_3/\eta_1)^{12}$ and $U=\eta_1^9/\eta_3^3=b(q)^3$, $b=\eta_1^3/\eta_3$, the
cube root in $t_0$ is
$$[t_0^n]\,b = (-27)^n\Bigl(\frac{(1/3)_n}{n!}\Bigr)^2 = 1,-3,36,-588,11025,\dots,
\qquad b={}_2F_1\bigl(\tfrac13,\tfrac13;1;-27t_0\bigr)$$
(verified as an identity of rationals for $n\le300$). Its recurrence is **first order**
($n^2a_n+(27n^2-36n+12)a_{n-1}=0$): the level-3 curve has only three singular points, so
there is a single characteristic root $-27$, no second solution with a finite Apéry limit,
and no companion. This is a genuine obstruction, not a failure of search.

**(c) The rational coordinates are direct images, and the root sees it.** In the coordinate
$z=r/(1+r)$ the $\beta(4)$ cube root admits *no* recurrence of order $\le26$ with
coefficients of degree $\le8$ (mod-$p$ scan on $401$ exact terms), consistent with the
$\operatorname{rank}f_*=d\cdot\operatorname{rank}$ principle of `book/v8/03_three_layers.tex`
(degree $4$ cover, rank $2$ upstairs $\Rightarrow$ rank $8$, with a $t$-degree well above the
scanned range). By contrast the $\zeta(5)$ level-16 and $\zeta(7)$ level-24 roots *are*
order-$2$ differential objects in their rational coordinates — their covers are the
identity resp. degree $2$.

### 5.4 The Apéry limits of the new root rows, and why they cannot be identified

| root row | $\xi_{\rm root}$ (best available) | convergence | identification |
|---|---|---|---|
| $\zeta(7)_{24}$ in $s$ | $-0.06209067158642\ldots$ ($n=600$) | $b_n/a_n-\xi\asymp n^{-1/2}$ | impossible: $\approx4$ correct digits |
| $\zeta(7)_{24}$ in $z$ | $-0.7936421895\ldots$ ($n=400$) | polynomial | impossible |
| $\zeta(5)_{16}$ | $-0.01069093429\ldots$ ($n=400$) | polynomial | impossible |
| $\beta(4)_{24}$ in $t$ | $G/2$ (Catalan) | geometric, $\lambda_2/\lambda_1=1/2$ | **known** (Zagier E) |
| $L(4,\chi_{-3})_{24}$ in $t_0$ | — | no companion | n/a |

The polynomial convergence is the signature of the doubled dominant root of §4(H4): the
exponent is exactly the exponent difference at the dominant singularity
($\tfrac7{12}-\tfrac1{12}=\tfrac12$ for $\zeta(7)$). Four digits cannot be fed to `lindep`,
so no identification was attempted; the honest statement is that these rows have an Apéry
limit but no usable one. The three unidentified $w=2$ limits
($0.14551448\ldots$, $0.31692536\ldots$, $0.48423755\ldots$) remain open exactly as in
`SQRT_APERY.md` §7 — nothing new was found for them here.

### 5.5 Catalan-type systems: nothing to take a root of

Every Catalan/$\beta(2)$ object in the project — Zagier E (level 8), the level-16 "K leaf",
the level-144 double-Eichler system, Zudilin's non-modular row, the Nesterenko Padé
construction — is a **weight-one** ($w=1$) row: $\beta(2)=L(2,\chi_{-4})$ sits at $r=2$ in
the $\operatorname{Sym}^{r-1}$ ladder of `book/v8/03_three_layers.tex`. There is no $w$-th
root to take. The one $\chi_{-4}$ object with $w\ge2$ is $\beta(4)$ at level 24, and §5.3(a)
shows its cube root is Zagier E itself — so the Catalan row is already the bottom of that
tower.

## 6. What this changes

1. `SQRT_APERY.md` Theorems 1 and $1'$ are the $w=2$ case of Theorems R1, R1$'$; Lemma 2.3
   / Prop 2.4 of `SQRT_APERY_FORMAL` are the $w=2$ case of Theorem R2. The
   formalization blueprint can be stated for general $w$ at no extra cost (the descent
   argument is one order-of-vanishing estimate).
2. Paper Problem 7.6 gains an exact statement: the free integration is worth
   $(k_{\rm parent}-2)-\log\lambda$, i.e. $(w-1)-\log(w\operatorname{rad}w)$ generically,
   **positive for every $w\ge4$**. The $\operatorname{Sym}^w$ ladder of
   `book/v8/03_three_layers.tex` is therefore *not* barren in principle: what kills it in
   practice is that the rational Apéry coordinate is a direct image (so the root is not a
   scalar row) or that the $w$-th root carries an algebraic gauge factor (so the
   characteristic roots double).
3. Conjecture `conj:barrier` (scalar barrier, $w\ge3$) is **strengthened, not threatened**:
   the natural attack "take the $w$-th root and harvest $w-1$ units of $d_n$" now has a
   proved cost ($\log\lambda_w$) and two identified obstructions ((H3), (H4)), and every
   $w\ge3$ system in the project hits one of them.
4. Census additions: $\beta(4)_{24}=\operatorname{Sym}^3(\text{Zagier E})$;
   $L(4,\chi_{-3})_{24}=\operatorname{Sym}^3({}_2F_1(\tfrac13,\tfrac13;1;\cdot))$;
   $k=2$ for the $\zeta(5)$ and $\zeta(7)$ root rows (first $w>2$ instances of the free
   integration).
5. **No irrationality claim is made anywhere in this note.** The only row with all four
   Apéry ingredients is Apéry's own square root, already proved in `SQRT_APERY.md`.

## 7. Reproduction

```
lattice/root_rows/01_lambda_w.py     # Theorem R1: valuations, minimality, series test
lattice/root_rows/02_descent.py      # Theorem R2: Sym^w operators and the identity, w=2..5
lattice/root_rows/lib.gp             # fitrec / minrec / rootrow / charroots
lattice/root_rows/03_sporadic.gp     # the nine w=2 rows, re-verified to n=420
lattice/root_rows/04_build.gp        # 400 exact A_n for every w>=3 system  -> rows_*.txt
lattice/root_rows/05_analyse.gp 06_run.gp   # lambda, integrality, first recurrence attempt
lattice/root_rows/07_scan.gp 08_widescan.gp # mod-p (order,degree) kernel scans
lattice/root_rows/09_zeta7_s.gp      # zeta(7) level 24 in the s-coordinate (w=6)
lattice/root_rows/10_oriented.gp     # beta(4) -> Zagier E ; L(4,chi_-3) -> 2F1(1/3,1/3;1;.)
lattice/root_rows/11_companion.gp    # companion from the order-2 operator; d_n^k b_n tests
```
