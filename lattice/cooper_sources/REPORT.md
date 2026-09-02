# Cooper's three sources: meromorphy, magnetism, and the congruence behind the companion Lucas law

*Working note, 2026-09-02.  All computations in PARI/GP 2.15.4, exact over **Q** unless a
precision is quoted.  Scripts `lib.gp`, `01_series.gp` … `05_magnetic.gp` with logs
`0*.log`; data files `s7_c.txt`, `s7_cp.txt`, … contain `c(m)` and `c'(m)=c(m)/m` for
`m = 1..400`.  Nothing in `/home/ubuntu/code/math-modular-sources` was modified.*

Claims are tagged **[proved]**, **[exact]** (an exact rational/series identity verified by
computation over a stated range), **[num]** (floating point, digits quoted), or
**[conj]**.

---

## 0. Summary of findings

1. **The sources are computed and are magnetic.**  $c(m)$ for $m\le400$ is produced from the
   eta quotients by two independent formulas which agree identically; $m\mid c(m)$ for all
   $m\le400$ in all three rows, and the companion formula
   $b_n=\sum_m c(m)m^{-3}e_{n,m}$ reproduces the recurrence's $b_n$ exactly for $n\le45$.

2. **All three sources are meromorphic on $\mathbf H$, $s_{18}$ included.**  The task's
   working hypothesis that $\Phi_{s_{18}}$ might be holomorphic is **false**.  The polar
   locus is, in every case, a *CM point*:

   | row | $\tau_0$ | min. poly | disc | what it is | $R=e^{2\pi\operatorname{Im}\tau_0}$ |
   |---|---|---|---|---|---|
   | $s_7$    | $(5\pm\sqrt{-3})/14$ | $7\tau^2-5\tau+1$ | $-3$ | the two elliptic points of order $3$ of $X_0(7)$ | $e^{2\pi\sqrt3/14}=2.175682\ldots$ |
   | $s_{10}$ | $(3+i)/10,\ (7+i)/10$ | $10\tau^2-6\tau+1$ | $-4$ | the two elliptic points of order $2$ of $X_0(10)$ | $e^{\pi/5}=1.874456\ldots$ |
   | $s_{18}$ | $(3+i)/6$ | $18\tau^2-18\tau+5$ | $-36$ | the fixed point of the **Atkin–Lehner involution $W_9$** ($X_0(18)$ has *no* elliptic points and the roots are *not* cusp values) | $e^{\pi/3}=2.849653\ldots$ |

   In each case $\Phi$ has a **double pole in $\tau$** with a **non-zero residue**, and the
   coefficient asymptotics are pinned exactly (§2.4).

3. **Task 4 is the headline.**  For every row and *every* prime $p\le43$ — good and
   nominally bad alike —
   $$\boxed{\;c'(p^im)\equiv\psi(p)^i\,c'(m)\pmod{p^{2}},\qquad p\nmid m,\ i=1,2,3\;}$$
   with $\psi=\mathbf 1,\mathbf 1,\chi_{-3}$, and the modulus $p^2$ is **sharp** (except
   $p=2,5$ for $s_7$, where it is $p^3$, and the two exact cells $p=7\mid s_7$,
   $p=5\mid s_{10}$).  Exact multiplicativity of $c'$ **fails** (818–972 failures out of 972
   coprime pairs).  This is precisely hypothesis (a) of `COMPANION_ARITHMETIC.md`
   Theorem 5.1 with $r-1=2$ in place of $r$ — in the *weakened, congruence-only* form that
   the proof actually consumes — so **Remark 5.1(3) is confirmed and the companion Lucas
   law $\beta_{np+m}\equiv\psi(p)\beta_na_m$, $\beta_n=p^{2\lfloor\log_pn\rfloor}b_n$,
   for Cooper's rows is now a consequence of Theorem 5.1's proof.**

4. **Equivalently, in Paşol–Zudilin's language:** $\Phi|U_p\equiv\psi(p)\,p\,\Phi\pmod{p^3}$
   for every $p\le43$ (sharp).  This implies their *strong $p$-magnetic property*
   $p^n\mid m\Rightarrow p^n\mid c(m)$, hence magnetism $m\mid c(m)$.  **The free integration
   and the companion Lucas law are the same congruence.**

5. **Identification.**  $\Phi$ is *not* a holomorphic form of weight 4 at any level; what is
   exactly identified is $\Phi=F^2\cdot u\,\frac{dx}{du}$ with $F\in M_2(\Gamma_0(N))$ the
   explicit $E_2$-combination and $F^2\in M_4(\Gamma_0(N))$ decomposed into Eisenstein +
   newform pieces (§3).  For $s_{18}$ the cuspidal part is $12(\eta(3\tau)^8+4\eta(6\tau)^8)$
   — the CM newform of level $9$ for $\mathbf Q(\sqrt{-3})$, which is presumably where
   $\psi=\chi_{-3}$ comes from.  The Eisenstein/character structure of the source itself
   survives **only** at the critical slot: $L(\Xi,2)=\zeta(2)/7,\ \zeta(2)/5,\
   \tfrac12L(2,\chi_{-3})$ [num, 60 digits], reconfirming $\psi$ independently of §4.

6. **Three corrections to the repository** (§6): $u$ is a degree-one Hauptmodul only for
   $s_7$ (degrees $1,2,4$); `WEIGHT_DROP.md` §4.3's sign $F|_2W_9=-F$ is wrong (it is $+1$)
   and its conclusion "the residue vanishes" is false; `WEIGHT_DROP.md` §3.1's measured
   growth rate $1.178$ for $s_{10}$ is an artefact of sampling $m\equiv0\ (5)$.

---

## 1. Task 1 — the coefficients $c(m)$

Conventions as in the task: $u$ the eta quotient, $F=D\log u=Du/u$,
$x=u/(1+Bu+Cu^2)$, $\Phi=F\,Dx=\sum_{m\ge1}c(m)q^m$, $D=q\,d/dq$.

**[exact]** $F\cdot Dx \;=\; x\sqrt{P(x)}\,F^2 \;=\; \dfrac{u(1-Cu^2)F^2}{(1+Bu+Cu^2)^2}$
to $O(q^{406})$, all three rows (`01_series.log`).  Both were computed from
$u=\prod_d\eta(d\tau)^{r_d}$ as $q$-series.

**[exact]** $c(m)\in\mathbf Z$ and **$m\mid c(m)$** for all $1\le m\le400$, all three rows.
Data: `s7_c.txt`, `s10_c.txt`, `s18_c.txt` (and `*_cp.txt` for $c'(m)=c(m)/m$).

| $m$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| $c^{s_7}$ | 1 | $-14$ | 30 | 100 | $-620$ | 1308 | 7 | $-9016$ | 28521 | $-31320$ |
| $c'^{s_7}$ | 1 | $-7$ | 10 | 25 | $-124$ | 218 | 1 | $-1127$ | 3169 | $-3132$ |
| $c^{s_{10}}$ | 1 | $-6$ | $-24$ | 116 | 5 | $-720$ | 1036 | 2280 | $-8091$ | $-30$ |
| $c'^{s_{10}}$ | 1 | $-3$ | $-8$ | 29 | 1 | $-120$ | 148 | 285 | $-899$ | $-3$ |
| $c^{s_{18}}$ | 1 | $-10$ | 54 | $-236$ | 870 | $-3132$ | 10640 | $-35368$ | 115182 | $-368700$ |
| $c'^{s_{18}}$ | 1 | $-5$ | 18 | $-59$ | 174 | $-522$ | 1520 | $-4421$ | 12798 | $-36870$ |

**[exact]** cross-checks (`01_series.log`): $a_n=[x^n]F$ agrees with the recurrence for
$n\le45$; $e_{n,m}=[x^n](Fq^m)\in\mathbf Z$; and
$b_n=\sum_{m\le n}c(m)m^{-3}e_{n,m}$ equals the recurrence's companion for $n\le45$.

**[exact]** two exact $U_p$-eigen relations at the special primes:
$$c(7m)=7\,c(m)\ (s_7),\qquad c(5m)=5\,c(m)\ (s_{10})\qquad (m\le400/p),$$
i.e. $\Xi|U_p=\Xi$ where $\Xi=\sum c'(m)q^m$.  No such relation for $s_{18}$ at $2$ or $3$.
(This is why $c(m)$ is tiny on multiples of $7$, resp. $5$; see §2.4.)

---

## 2. Task 2 — meromorphy: the exact polar data

### 2.1 Where the poles are

$\Phi=F^2\,\rho(u)$ with $\rho(u)=u\,\frac{dx}{du}=\dfrac{u(1-Cu^2)}{g(u)^2}$,
$g(u)=1+Bu+Cu^2$.  Since $F^2$ is holomorphic on $\mathbf H$, the polar locus of $\Phi$ is
exactly $\{\tau: g(u(\tau))=0\}$, i.e. the two $\Gamma_0(N)$-orbits with $u=u_\pm$. **[proved]**

**[num, 60 digits]** the roots $u_\pm$ and the points attaining them:

* $s_7$: $u_\pm=(-13\pm3\sqrt{-3})/98$ attained at $\tau_0=(5+\sqrt{-3})/14$ and
  $(9+\sqrt{-3})/14=-\overline{\tau_0}+1$: **the two elliptic points of order $3$** of
  $X_0(7)$ (fixed points of $\binom{\ 3\ -1}{\ 7\ -2}$, trace $1$), exchanged by $W_7$.
* $s_{10}$: $u_\pm=(-3\pm4i)/25$ attained at $(3+i)/10$ and $(7+i)/10$: **the two elliptic
  points of order $2$** of $X_0(10)$, exchanged by $W_{10}$.
* $s_{18}$: $u_+=-7+4\sqrt3$ at $\tau_0=(3+i)/6$ and $u_-=-7-4\sqrt3$ at
  $W_{18}\tau_0=-\tfrac1{10}+\tfrac{i}{30}$.  $X_0(18)$ has **no** elliptic points and
  $\tau_0$ is **not** a cusp; it is the fixed point of $W_9$
  ($\binom{\ 9\ -5}{18\ -9}$, det $9$).  **[num]** $u|W_9=u$, $x|W_9=x$, $F|_2W_9=+F$,
  $\Phi|_4W_9=+\Phi$ (48 digits).

All three $\tau_0$ are CM points, of discriminants $-3$, $-4$, $-36$ respectively
(the last is the conductor-$3$ order in $\mathbf Q(i)$).

The subdominant pole: for $s_7,s_{10}$ the second pole is the complex conjugate, at the
*same* height; for $s_{18}$ the $\Gamma_0(18)$-orbit of $W_{18}\tau_0$ has maximal height
$1/30$ **[num search over $|c|\le6$]**, so it contributes $O((e^{\pi/15}/e^{\pi/3})^m)$.

### 2.2 A correction: the degree of $u$

**[exact, Ligozat]** divisors of the three eta quotients on $X_0(N)$
(`02_poles.log`; $[c,\#\text{cusps},\operatorname{ord}]$):

| row | divisor | $\deg u$ |
|---|---|---|
| $s_7$ | $[1,1,-1],[7,1,1]$ | $1$ |
| $s_{10}$ | $[1,1,-1],[2,1,-1],[5,1,1],[10,1,1]$ | $2$ |
| $s_{18}$ | $[1,1,-1],[2,1,1],[3,2,1],[6,2,-1],[9,1,-1],[18,1,1]$ | $4$ |

So **$u$ is a Hauptmodul of $X_0(N)$ only for $s_7$**.  Consequently $x$ has degree $2,4,8$
on $X_0(7),X_0(10),X_0(18)$.  See §6.1 for the (nil) impact on the repository's theorems.

### 2.3 The order of vanishing $\nu$ and the exact principal part

Put $\nu:=\operatorname{ord}_{\tau=\tau_0}\bigl(u-u_0\bigr)$ (order in $\tau-\tau_0$, not in
the local coordinate of the curve).  **[num, 60 digits, 5 step sizes]**
$$\nu = 3\ (s_7),\qquad \nu=4\ (s_{10}),\qquad \nu=4\ (s_{18}).$$
($s_7$: $\deg u=1$ but $\tau_0$ is an order-3 elliptic point, so $3\cdot1$.  $s_{10}$:
$\deg u=2$, order-2 elliptic point, total ramification, $2\cdot 2$.  $s_{18}$: $\deg u=4$,
total ramification over $u_\pm$, $4\cdot1$.)

**[proved]** With $y=\tau-\tau_0$, $\Phi=A_2y^{-2}+A_1y^{-1}+O(1)$ and
$$A_2=\frac{\nu^2}{4\pi^2\,g'(u_0)},\qquad g'(u_0)=B+2Cu_0,\qquad
g'(u_0)^2=B^2-4C=\lambda_1\lambda_2=c\ (\text{the recurrence constant}).$$
*Proof.*  $u-u_0=\alpha y^\nu+\cdots$ gives $F=Du/u=\nu\alpha y^{\nu-1}/(2\pi iu_0)+\cdots$
and $g(u)^2=g'(u_0)^2\alpha^2y^{2\nu}+\cdots$, so
$\Phi=F^2u(1-Cu_0^2)/g^2\to -\nu^2(1-Cu_0^2)/(4\pi^2u_0g'(u_0)^2y^2)$, and
$u_0g'(u_0)=-(1-Cu_0^2)$ because $g(u_0)=0$.  $\square$
Numerically $A_2^{\rm num}/A_2^{\rm formula}=1$ to 13 digits in all three rows.

**[num, 14 digits]** the residue is **non-zero** and equals
$$A_1=\frac{i\,A_2}{\operatorname{Im}\tau_0}=-\frac{2A_2}{\tau_0-\overline{\tau_0}} .$$
This is exactly the statement that the quadratic differential $\Phi\,(d\tau)^2$ is **even**
in the elliptic coordinate $w=(\tau-\tau_0)/(\tau-\overline{\tau_0})$ (equivalently
$\Phi|_4\sigma=+\Phi$ for the elliptic $\sigma$ fixing $\tau_0$), which forces the pole order
to be even and forbids a $w^{-1}$ term — but a $w^{-2}$ term produces a $y^{-1}$ term.
**[proved, given $\Phi|_4\sigma=+\Phi$ and $\operatorname{ord}_w\Phi=-2$]**

Explicit values: $g'(u_0)=3i\sqrt3,\ 8i,\ 8\sqrt3$ for $s_7,s_{10},s_{18}$
(so $A_2=\frac{9}{4\pi^2\cdot3i\sqrt3},\ \frac{16}{4\pi^2\cdot8i},\
\frac{16}{4\pi^2\cdot8\sqrt3}$).

### 2.4 Growth of $c(m)$: exponential, with the constant pinned

**[proved from §2.3 + [num] to 60 digits at $m=396,\dots,399$]**
$$c(m)\;=\;-\nu^2\Bigl(\frac{q_0^{-m}}{g'(u_0)}+\overline{\phantom{x}\cdots\phantom{x}}\Bigr)
\Bigl(m-\frac1{2\pi\operatorname{Im}\tau_0}\Bigr)\;+\;O\bigl(R_2^{\,m}\bigr),
\qquad q_0=e^{2\pi i\tau_0},$$
the conjugate term being present for $s_7,s_{10}$ (two conjugate poles) and absent for
$s_{18}$ ($q_0=-e^{-\pi/3}$ real, $A_2$ real).  The measured ratio
$c(m)/(\text{right side})$ is $1.000\ldots0$ to **60 digits** at $m=396,397,398$ for all
three rows.  So the whole coefficient sequence is governed by one double pole.

Consequences:
* $R:=\limsup|c(m)|^{1/m}=e^{2\pi\operatorname{Im}\tau_0}$: $2.175682\ldots$, $1.874456\ldots$,
  $2.849653\ldots$ — i.e. **polynomial growth is excluded in all three rows; $\Phi$ is
  meromorphic with poles inside $\mathbf H$ in every case, and holomorphic weight-4 identification
  is impossible.**
* Specialised: $s_{18}$: $c(m)\sim-\tfrac{2}{\sqrt3}(-1)^m m\,e^{\pi m/3}$;
  $s_7$: $c(m)\sim 2\sqrt3\,m\,R^m\sin(5\pi m/7)$;
  $s_{10}$: $c(m)\sim 4\,m\,R^m\sin(3\pi m/5)$.
  The last two vanish to leading order exactly when $7\mid m$, resp. $5\mid m$ — which is
  the analytic shadow of the exact relations $c(7m)=7c(m)$, $c(5m)=5c(m)$ of §1.
* $\operatorname{Im}\tau_0<1/\sqrt N$ in all three cases
  ($0.1237<0.3780$, $0.1<0.3162$, $0.1667<0.2357$), so the Fricke geodesic $(0,i\infty)$
  avoids the poles and `WEIGHT_DROP.md` Theorem WD3 is unaffected. **[proved]**

---

## 3. Task 3 — identification

### 3.1 What is *not* true

**[exact]** `mftobasis` on the full spaces $M_4(\Gamma_0(N))$ (dimensions $3,7,13$) returns a
best fit whose defect against $\Phi$ at coefficient $\le150$ is of size
$10^{52},10^{43},10^{70}$: $\Phi\notin M_4(\Gamma_0(N))$ for **all three rows**, $s_{18}$
included.  Likewise no Eisenstein-plus-cusp expression can exist, because $|c(m)|$ grows like
$R^m$ with $R>1$ (§2.4).  Task 3(a) as posed therefore has no solution; the premise
"for $s_{18}$ the roots sit at cusps" is false (§2.1).

### 3.2 What is exactly true

**[exact, $O(q^{406})$]** the identification of the source is
$$\boxed{\ \Phi \;=\; F^2\cdot\frac{u\,(1-Cu^2)}{(1+Bu+Cu^2)^2}\ },$$
a holomorphic weight-4 form times an explicit degree-$(3,4)$ rational function of the
eta quotient $u$, with the double poles of §2 at $u=u_\pm$.  Partial fractions:
$\dfrac{u(1-Cu^2)}{g(u)^2}=\dfrac1u-\dfrac{2B+(B^2+3C)u+2BCu^2+C^2u^3}{g(u)^2}$.

**[exact]** the two factors, identified in standard bases (`03_ident.log`):

$F\in M_2(\Gamma_0(N))$ (defect $0$ over 120 coefficients):
$$F_{s_7}=\tfrac16\bigl(7E_2(7\tau)-E_2(\tau)\bigr),\qquad
F_{s_{10}}=\tfrac1{12}\bigl(10E_2(10\tau)+5E_2(5\tau)-2E_2(2\tau)-E_2(\tau)\bigr),$$
$$F_{s_{18}}=\tfrac14\bigl(18E_2(18\tau)-9E_2(9\tau)-12E_2(6\tau)+6E_2(3\tau)+2E_2(2\tau)-E_2(\tau)\bigr).$$

$F^2\in M_4(\Gamma_0(N))$ (defect $0$ over 199 coefficients, verified against explicit
series in `03_ident.log`):
$$F_{s_7}^2=\tfrac1{50}\bigl(E_4(\tau)+49E_4(7\tau)\bigr)+\tfrac{16}5 f_{7.4.a.a},$$
$$F_{s_{10}}^2=\tfrac1{130}\bigl(E_4(\tau)+4E_4(2\tau)+25E_4(5\tau)+100E_4(10\tau)\bigr)
+\tfrac{28}{13}\bigl(f_{5.4.a.a}(\tau)+4f_{5.4.a.a}(2\tau)\bigr),$$
$$F_{s_{18}}^2=\tfrac15\bigl(E_4(3\tau)+4E_4(6\tau)\bigr)
+12\bigl(\eta(3\tau)^8+4\,\eta(6\tau)^8\bigr).$$
In each case the Eisenstein part is $\bigl(\sum_{d\in S}d^2E_4(d\tau)\bigr)/\sum_{d\in S}d^2$
with $S=\{1,7\},\{1,2,5,10\},\{3,6\}$, each $S$ stable under $d\mapsto N/d$.
$f_{7.4.a.a}=q-q^2-2q^3-7q^4+16q^5+\cdots$ and $f_{5.4.a.a}=q-4q^2+2q^3+8q^4-5q^5+\cdots$
are the level-$7$ and level-$5$ weight-$4$ newforms; $\eta(3\tau)^8$ is the level-$9$
weight-$4$ **CM** newform for $\mathbf Q(\sqrt{-3})$ ($a_p=0$ for $p\equiv2\bmod3$) — the
only visible source of the character $\chi_{-3}$ that governs $s_{18}$.

### 3.3 The character data, and where the Eisenstein shape does live

The $L$-function of $\Phi$ has **no** Euler product and no Dirichlet-polynomial
factorisation: $c'$ is not multiplicative (§4), and $\sum c'(m)m^{-s}$ does not even
converge.  What *is* Eisenstein is the **critical slot** of $\Xi=D^{-1}\Phi$.
Using $\Lambda(\Phi,s)=I(s)-N^{2-s}I(4-s)$ with
$I(s)=\sum_{m\le400}c(m)(2\pi m)^{-s}\Gamma(s,2\pi m/\sqrt N)$ **[num, 60 digits]**:
$$\Lambda(\Phi,2)=0\ (<10^{-79}),\qquad
L(\Xi,2)=4\pi^3\Lambda(\Phi,3)=\frac{\zeta(2)}7,\ \frac{\zeta(2)}5,\
\frac{L(2,\chi_{-3})}2 .$$
This reproduces `WEIGHT_DROP.md` (V3) from scratch and reconfirms
$\psi_{s_7}=\psi_{s_{10}}=\mathbf 1$, $\psi_{s_{18}}=\chi_{-3}$ — *independently of the
congruences of §4*, which give the same answer by a completely different route.

---

## 4. Task 4 — the congruences behind the companion Lucas law

Let $\psi=\mathbf1,\mathbf1,\chi_{-3}$ for $s_7,s_{10},s_{18}$ and $c'(m)=c(m)/m$.

### 4.1 (i)+(ii): the main congruence

**[exact, all $m\le400/p^i$, $i=1,2,3$]** for every prime $p\le43$
$$c'(p^i m)\equiv\psi(p)^i\,c'(m)\pmod{p^{2}},\qquad p\nmid m .$$
Table of $e=v_p\bigl(\gcd\text{ of all defects}\bigr)$ (`04_congr.log`; "inf" = identically
zero; #tests $175,129,80,57,36,30,23,21,17,13,12,10,9,9$ for $p=2,\dots,43$):

| $p$ | 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $s_7$    | **3** | 2 | **3** | inf | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| $s_{10}$ | 2 | 2 | inf | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| $s_{18}$ | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |

So the congruence holds **modulo $p^2$ and no better** in 38 of the 42 cells, including the
nominally bad primes ($p=7$ for $s_7$; $p=2,5$ for $s_{10}$; $p=2,3$ for $s_{18}$), and
better in exactly two cells.  For $s_{18}$ at $p=3$ one has $\psi(3)=0$ and the statement
reads $c'(3^im)\equiv0\pmod 9$.

### 4.2 (iii): exact multiplicativity fails

**[exact]** number of coprime pairs $(a,b)$, $ab\le400$, with $c'(ab)\ne c'(a)c'(b)$:
$864/972$ ($s_7$), $818/972$ ($s_{10}$), $972/972$ ($s_{18}$).  As expected for a
meromorphic source.  *(The successes for $s_7,s_{10}$ are exactly the pairs involving the
special prime, where $c'(p^am)=c'(m)$ holds identically.)*

### 4.3 (iv): $c'(p)$ against $\psi(p)$

**[exact]** $v_p\bigl(c'(p)-\psi(p)\bigr)=2$ for every $p\le43$ and every row, except
$s_7$ at $p=2,5$ (value $3$) and the two "inf" cells.  Samples:
$c'^{s_{10}}(13)=-7266=-43\cdot169+1$; $c'^{s_{18}}(11)=106116=877\cdot121-1$ with
$\chi_{-3}(11)=-1$; $c'^{s_7}(11)=-6896=-57\cdot121+1$.

The next digit is *not* Eisenstein: $\lambda(p):=(c'(p)-\psi(p))/p^2\bmod p$ takes the
values (e.g. $s_{18}$) $1,2,2,3,8,12,2,5,16,14,16,6,15,39$ for $p=2,\dots,43$ — not a
character, so there is **no** $\chi$ with $c'(p)\equiv\psi(p)+\chi(p)p^2\pmod{p^3}$.  The
mod-$p^2$ shape is exactly as much Eisenstein structure as the source carries.

### 4.4 Consequence for Theorem 5.1

`COMPANION_ARITHMETIC.md` Theorem 5.1 assumes (a) $c(p^im)=c(p^i)c(m)$ and
$c(p^i)\equiv\psi(p)^i\pmod{p^r}$.  Its proof, however, only uses these through
$$\Theta_i\equiv\psi(p)^i\,\Theta_0\pmod{p^{r}},\qquad
\Theta_i(q)=\sum_{p\nmid m}c(p^im)m^{-r}q^m .$$
Writing $c(m)=m\,c'(m)$ one has $\Theta=D^{-3}\Phi=\sum_m c'(m)m^{-2}q^m$, so the grouping
by $p$-power reads $\Theta=\sum_i p^{-2i}\Theta_i(q^{p^i})$ with
$\Theta_i=\sum_{p\nmid m}c'(p^im)m^{-2}q^m$, and the hypothesis needed is precisely
$$\Theta_i\equiv\psi(p)^i\Theta_0\pmod{p^{2}}\iff c'(p^im)\equiv\psi(p)^ic'(m)\pmod{p^{2}},$$
verified in §4.1.  **[proved, modulo §4.1 and Malik–Straub]** the proof of Theorem 5.1 then
goes through verbatim with $r$ replaced by $r-1=2$ and yields
$$\beta_n\equiv\psi(p)^j\,b_{n_j}\prod_{i<j}a_{n_i}\pmod p,\qquad
\beta_n=p^{2\lfloor\log_pn\rfloor}b_n,$$
i.e. exactly the census's $k_p=2$ Lucas law for Cooper's three rows, with
$\psi=\mathbf1,\mathbf1,\chi_{-3}$ — and, because §4.1 holds at the *bad* primes too, it
explains the census's observation that $k_p=2$ works "at every prime including the nominally
bad $3$, resp. $2$", and that $u=0$ at $p=3$ for $s_{18}$ ($\psi(3)=0$).
**Remark 5.1(3) is confirmed**, in the corrected form: the hypothesis is a congruence, not
exact multiplicativity.  (Note that "$c'$ is the coefficient system of a weight-$(r-1)$
Eisenstein series" is *false* as an identity — §4.2, §4.3 — but true modulo $p^2$, which is
all the theorem needs.)

---

## 5. Task 5 — why $m\mid c(m)$

### 5.1 The two elementary equivalences **[proved]**

Since $x(q)\in q+q^2\mathbf Z[[q]]$ and hence $q(x)\in x+x^2\mathbf Z[[x]]$,
$$m\mid c(m)\ \forall m
\iff \Xi:=D^{-1}\Phi=\sum_m\tfrac{c(m)}m q^m\in\mathbf Z[[q]]
\iff \int_0^x\!A\,dx=\sum_n\tfrac{a_n}{n+1}x^{n+1}\in\mathbf Z[[x]]
\iff (n+1)\mid a_n .$$
The last is Bogner's theorem; **[exact]** we re-verify $(n+1)\mid a_n$ for $n\le3000$ in all
three rows.  It is *sharp*: $\min_n\bigl(v_p(a_n)-v_p(n+1)\bigr)=0$ for every $p\le23$ and
every row, except $s_7$ at $p=2$ (min $1$) and $s_{18}$ at $p=3$ (min $2$).

### 5.2 The literature: magnetic modular forms

V. Paşol and W. Zudilin, *Magnetic (quasi-)modular forms*, Nagoya Math. J. **248** (2022)
(arXiv:2009.14609), following D. Broadhurst–W. Zudilin and Y. Li–M. Neururer:

* a meromorphic modular form of weight $2k$ is **magnetic** if $n^{k-1}\mid c(n)$ — for
  $k=2$ exactly our $m\mid c(m)$ — and **strongly magnetic** if $\delta^{-1}F\in q\mathbf Z[[q]]$,
  exactly our $\Xi\in\mathbf Z[[q]]$;
* **folklore conjecture**: *no holomorphic* $F\in q\mathbf Z[[q]]$ is magnetic.  **This is
  the conceptual reason Cooper's sources have to be meromorphic**: the free integration
  ($d_n^{\,2}b_n\in\mathbf Z$ rather than $d_n^{\,3}b_n$) forces magnetism, and magnetism
  forbids holomorphy.  The repository's "meromorphic, unidentified" is not an accident of
  the parametrisation; it is compulsory. **[conj, conditional on the folklore conjecture]**
* their level-$1$ examples of weight $4$ are $\Delta/E_4^2$ (double pole at the order-$3$
  elliptic point $\rho$, disc $-3$) and $E_4\Delta/E_6^2$ (double pole at $i$, disc $-4$);
  the proof route is a Shimura–Borcherds lift, and the mechanism extracted is the
  **strong $p$-magnetic property** $F|U_{p^n}\equiv0 \pmod{p^n}$, i.e. $p^n\mid m\Rightarrow
  p^n\mid c(m)$.  They add: "constructing magnetic modular forms — meromorphic ones with
  poles at quadratic irrationalities from the upper half-plane — is a routine on the basis
  of the Shimura–Borcherds lift", and "there is good reason to believe that all such
  magnetic forms originate from suitable Shimura–Borcherds lifts."

**Cooper's three sources are the higher-level members of exactly this family**: weight-4
meromorphic forms on $\Gamma_0(N)$ whose only pole in $\mathbf H$ is a double pole at a CM
point, of discriminant $-3$ ($s_7$, matching $\Delta/E_4^2$), $-4$ ($s_{10}$, matching
$E_4\Delta/E_6^2$) and $-36$ ($s_{18}$, a conductor-$3$ order in $\mathbf Q(i)$).

### 5.3 The mechanism, made explicit here

**[exact, all $m\le400/p$, every prime $p\le43$]**
$$\boxed{\ \Phi\big|U_p\ \equiv\ \psi(p)\,p\,\Phi \pmod{p^{3}}\ }$$
with $e=3$ sharp in 38 of 42 cells, $e=4$ for $s_7$ at $p=2,5$, and identically for
$p=7\,|\,s_7$, $p=5\,|\,s_{10}$ (`05_magnetic.log`).

**[proved]** this is *equivalent* to §4.1 ($c(p^im)-\psi(p)^ip^ic(m)=p^im\,[c'(p^im)-\psi(p)^ic'(m)]$),
and it gives at once — **[exact]**, verified cell by cell for $p\le19$, $n\le3$, with
$e=n+2$ *exactly* in every cell except the identical ones and $s_7$ at $p=2,5$ ($e=n+3$)
and $s_{18}$ at $p=3$ ($e=3n$) —
$$\Phi\big|U_{p^n}\equiv(\psi(p)p)^n\,\Phi\pmod{p^{\,n+2}}\ \Longrightarrow\
p^n\mid m\Rightarrow p^n\mid c(m)\ \Longrightarrow\ m\mid c(m),$$
i.e. **Paşol–Zudilin's strong $p$-magnetic property, with two extra powers of $p$ and a
character**.  **[exact]** the intermediate statement $p^n\mid m\Rightarrow p^n\mid c(m)$ was
checked directly for all $m\le400$.

So there is one congruence, $\Phi|U_p\equiv\psi(p)p\,\Phi\ (p^3)$, and it does both jobs:

* dividing by $m$ it is hypothesis (a) of Theorem 5.1 with $r-1=2$, hence the **companion
  Lucas law** with $k_p=2$ and character $\psi$;
* summing over $n$ it is the strong $p$-magnetic property, hence **Cooper's free
  integration**.

That the free integration and the companion Lucas law have the same single source is, as far
as I can see, new; it also explains why the census found $k_p=2$ (not $3$) and why it found
it uniformly at good *and* bad primes.

### 5.4 What can be proved unconditionally about $(n+1)\mid a_n$

* **[proved]** the case $p\,\|\,(n+1)$: if $p^{\,}\mid(n+1)$ then the base-$p$ digit
  $n_0=p-1$, and the Malik–Straub Lucas congruence gives $a_n\equiv a_{p-1}\prod_{i\ge1}a_{n_i}\pmod p$;
  **[exact]** $v_p(a_{p-1})\ge1$ for every $p\le43$ and every row
  (equal to $1$ except $v_2(a_1^{s_7})=2$, $v_3(a_2^{s_{10}})=2$, $v_3(a_2^{s_{18}})=3$),
  so $p\mid a_n$.  For $s_{10}=\sum_k\binom nk^4$ the input is one line:
  $\binom{p-1}k\equiv(-1)^k$, so $a_{p-1}\equiv\sum_{k=0}^{p-1}1=p\equiv0\pmod p$.
* **[proved]** the whole $2$-part for $s_{10},s_{18}$: the exact digit law
  $v_2(a_n)=s_2(n)$ (`COMPANION_ARITHMETIC.md` §4.3) plus the elementary fact that
  $2^e\,\|\,(n+1)$ forces the last $e$ binary digits of $n$ to be $1$, whence
  $s_2(n)\ge e=v_2(n+1)$.
* the general $p^e$ case would need a Dwork-type congruence modulo $p^e$; not attempted.

---

## 6. Corrections to the repository

1. **`ASYMPTOTIC_CONSTANTS.md` Thm 3.4 / `lattice/asymptotic_constants/REPORT.md` §2.1**
   call $u$ a "normalised degree-one Hauptmodul of $\Gamma_0(N)$".  **[exact]** $\deg u=1,2,4$
   for $s_7,s_{10},s_{18}$ (Ligozat, §2.2), so $u$ is a Hauptmodul only for $s_7$ (and, of
   the other four weight-two rows, the claim was not re-checked here).  *No consequence for
   the constants*: Theorem 3.4 uses only $u|W_N=1/(Cu)$, $F=D\log u$, $r'(u_c)=0$ and
   $Du(\tau_c)\ne0$, all of which are degree-independent.  What does change is the geometry
   of $x$: it has degree $2,4,8$ on $X_0(N)$, so $x$ is a Hauptmodul of $\Gamma_0(N)^+$ only
   for $s_7$.
2. **`WEIGHT_DROP.md` §4.3** asserts $F|_2W_9=-F$ for $s_{18}$ and concludes that the residue
   of $\Phi\,d\tau$ at $\tau_0=(3+i)/6$ vanishes, so that "$\Xi$ is single-valued near
   $\tau_0$".  **[num, 48 digits]** $F|_2W_9=+F$, $\Phi|_4W_9=+\Phi$, and **[num, 14 digits]**
   the residue is $A_1=iA_2/\operatorname{Im}\tau_0\ne0$.  (The sign as stated is in fact
   impossible: $\Phi|_4\sigma=-\Phi$ would force an odd-order pole, while the pole has order
   $2$.)  The correct statement is the *even-order* one of §2.3; the corrected residue is
   visible in the $c(m)$ asymptotics as the shift $m\mapsto m-1/(2\pi\operatorname{Im}\tau_0)$,
   confirmed to 60 digits.  Nothing else in `WEIGHT_DROP.md` depends on the vanishing.
3. **`WEIGHT_DROP.md` §3.1** reports measured $1/|q_0|\approx2.26,\,1.178,\,2.928$ for
   $s_7,s_{10},s_{18}$.  The middle value is an artefact: $c(5m)=5c(m)$ makes
   $|c(m)|^{1/m}\approx R^{1/5}$ on multiples of $5$.  The true rates are
   $2.175682\ldots$, $1.874456\ldots$, $2.849653\ldots$ (§2.4).  All three remain
   $<e^{2\pi/\sqrt N}$, so the section's conclusion is unaffected.
4. **`COMPANION_ARITHMETIC.md` Remark 5.1(3)** predicts that Cooper's $c'$ satisfies (a) with
   $r-1$ in place of $r$, "i.e. the meromorphic source is $D$ of a weight-$(r-1)$
   Eisenstein-type series".  The congruence half is confirmed (§4.1); the "is $D$ of an
   Eisenstein series" half is **false** as an identity — $\Xi$ has exponentially growing
   coefficients and $c'$ is not multiplicative — and should be weakened to
   "$c'$ is Eisenstein-like modulo $p^2$ at every prime".

---

## 7. Open

* Prove $\Phi|U_p\equiv\psi(p)p\,\Phi\pmod{p^3}$.  The Paşol–Zudilin route (Shimura–Borcherds
  lift of a weight-$5/2$ weakly holomorphic form in the Kohnen plus space) is level $1$;
  the level-$N$ analogue with the CM points of disc $-3,-4,-36$ would be the natural
  attack, and would give the companion Lucas law for Cooper's rows unconditionally.
* Why $p^3$ and not $p^4$?  And why $p^4$ exactly at $(s_7,p=2)$ and $(s_7,p=5)$?
* Identify $\lambda(p)=(c'(p)-\psi(p))/p^2\bmod p$ (§4.3).  It is not a character;
  is it the reduction of a Hecke eigenvalue of the cuspidal part of $F^2$ ($f_{7.4.a.a}$,
  $f_{5.4.a.a}$, $\eta(3\tau)^8$)?
* Are Cooper's three the only Apéry-like rows whose source is magnetic?  The criterion
  "$\Phi$ has a double pole at a CM point and is otherwise holomorphic" is checkable across
  the fifteen sporadic rows: for the twelve Eisenstein rows $g(u)$ has no root in $\mathbf H$
  in the relevant sense (their sources are holomorphic), so the answer is presumably yes.
* $\nu=3,4,4$: is the total ramification of $u$ over $u_\pm$ forced by $u|W_N=1/(Cu)$
  together with the fold structure?

---

## 8. Files

| file | contents |
|---|---|
| `lib.gp` | rows, eta quotients, $F$, $x$, $\Phi$, two formulas for $\Phi$, recurrence |
| `01_series.gp/.log` | Task 1: $c(m)$, magnetism, companion cross-check |
| `02_poles.gp/.log` | Task 2: Ligozat divisors, pole locations, $\nu$, principal parts, growth |
| `03_ident.gp/.log` | Task 3: $F$, $F^2$, $\Phi\notin M_4$, $W_9$/$W_{18}$ signs, $L(\Xi,2)$ |
| `04_congr.gp/.log` | Task 4: the mod-$p^2$ congruences, $c'(p)$, multiplicativity |
| `05_magnetic.gp/.log` | Task 5: $\Phi|U_p$ and $\Phi|U_{p^n}$, strong $p$-magnetic, $(n+1)\mid A_n$, $v_p(A_{p-1})$ |
| `s{7,10,18}_c.txt`, `s{7,10,18}_cp.txt` | $c(m)$ and $c'(m)=c(m)/m$, $m=1..400$ |
| `pz.txt` | text of arXiv:2009.14609 (Paşol–Zudilin), for reference |

Run with `gp -q <file>`; total runtime under 5 s.
