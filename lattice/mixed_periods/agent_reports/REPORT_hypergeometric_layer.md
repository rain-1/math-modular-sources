# Denominator layers and fold periods on the host $H_A=(1-4mx)^{-1/2}$

All computations were done in this directory (see the **Files** table at the end).
PARI/GP 2.15.4 for exact rational power series; python3 with mpmath 1.3 and `fractions` for
numerics and for an independent re-derivation. Every claim is tagged
`[verified exact, range]` or `[numerical, digits]`.

**Conventions.** $m\ge 1$ integer, $D'=4m-1$, $a=\sqrt{D'}$, $\theta_m=2\arctan(1/a)$,
$\delta=1/(4m)$.
$$H_A(x)=\frac{1}{\sqrt{1-4mx}}=\sum_{n\ge0}\binom{2n}{n}m^nx^n,\qquad
H[k](x)=H_A(x)\int_0^x H_A(t)k(t)\,dt .$$
$L(n):=\operatorname{lcm}(1,\dots,\lfloor n\rfloor)$.  A series $\sum a_nx^n$ has *type* $T$
if $T(n)\,a_n\in\mathbf Z$ for all $n$; the **excess** over $T$ is
$\max_n \operatorname{den}(a_n)/\gcd(\operatorname{den}(a_n),T(n))$; excess $1$ = exact type.
$H_B=H[(1-t)^{-1}]$, $H_C=H[\log(1-t)/t]$, $H_D=H[\log(1-t)/(1-t)]$.

---

## Summary

1. **[verified exact, $n\le 400$, $m\in\{1,2,3,4,5,6,7,11,17,41\}$]**
   $H_A$ is integral; $H_B$ has exact type $[n]$; $H_D$ has exact type $[n][n/2]$;
   $H_C$ has exact type $[n]^2$ and its excess over $[n][n/2]$ is **unbounded**,
   $e_n=\exp\!\big(n/3+o(n)\big)$ — in fact $e_n=\big(\prod_{2n/3<p\le n}p\big)\cdot(\text{small cofactor})$.

2. **[verified exact, $n\le 200$/$400$, $m=1,2$]** Over the kernel family
   $k=t^i(1-t)^{-j}\log(1-t)^{\ell}H_A(t)^{2s}$ the denominator type of $H[k]$ is
   **completely independent of the outer factor** ($H_A$, $H_A^3$, $H_A/(1-x)$ all give the
   same type, in all 176 cases tested), and follows an exact ladder
   $$1\ \subset\ [n]\ \subset\ [n][n/2]\ \subset\ [n][n/2][n/3]:$$
   for $j\ge1$, $H[t^i(1-t)^{-j}\log^\ell(1-t)]$ has **exact** type $[n][n/2]\cdots[n/(\ell+1)]$.
   The $1/t$ kernel is the exception: $H_C$-type kernels land on exact type $[n]^2$, which is
   *not* in the ladder ($H_C$ is not of type $[n][n/2][n/3]$ either).
   The stated expectation "$\log^2$ kernels fail" is **refuted for $j=0$**: $H[\log^2(1-t)]$
   has exact type $[n][n/2]$.

3. **[numerical, $\ge110$ digits]** $c_B=\theta_m/a$ and
   $c_D=(\theta_m\log(D'/m)-2\operatorname{Cl}_2(\pi-\theta_m))/a$ confirmed for $m=1,2,3,5,11$.
   **The whole $[n][n/2]$ layer has period module exactly $\mathbf Q+\mathbf Q c_B+\mathbf Q c_D$**
   (3-dimensional): *every* period of a $\le[n][n/2]$ kernel is a rational combination of
   $1,c_B,c_D$ — including the new one, for which
   $$c[\log^2(1-t)]=\frac{4-2D'c_B+D'c_D}{m}\qquad[\text{verified }10^{-150},\ m=1,2,3,5,11,17].$$
   **No period of the $[n][n/2]$ layer lies outside $\mathbf Q+\mathbf Q c_B+\mathbf Q c_D$.**
   By contrast the next layer does escape: $c[\log^2(1-t)/(1-t)^j]$, $j\ge1$, admits no relation
   with $\{1,c_B,c_D\}$ up to `maxcoeff` $10^{10}$ at $10^{-120}$.

4. **[numerical, 900-digit arithmetic]** Fold-regularity confirmed for $m=1,2,3,5$:
   $|a_n(H_B-c_BH_A)|^{1/n}$ and $|a_n(H_D-c_DH_A)|^{1/n}\to1$ (0.978–0.987 at $n=400$),
   while the unsubtracted coefficients have $|a_n|^{1/n}\to 4m$. The subtracted coefficients
   decay like $C/n$ ($H_B$ case, exactly: $n|a_n|$ is constant to 4 digits) and like
   $C\log n/n$ ($H_D$ case) — radius exactly 1, with a $\log(1-x)$ resp. $\log^2(1-x)$
   singularity at $x=1$. 900-digit arithmetic is *required*: the cancellation at $n=400$,
   $m=3$ is about 430 decimal digits, so the requested 80 digits would return pure noise.

5. **[numerical, 56–120 digits; also derived in closed form]** Cassaigne–Maillot confirmed
   (17 triples, worst case 56.0 digits) and Smyth's $m(1+x+y)=\frac{3\sqrt3}{4\pi}L(2,\chi_{-3})$
   to 56.5 digits. The integer-triple search $a,b,c\le60$ is **null** for $m=2,3,5$ (31 079
   triples, 0 hits) — but the right triple is algebraic, and for **every** $m$
   $$Q_m:=2\mathrm{Cl}_2(\pi-\theta_m)-\theta_m\log\tfrac{D'}{m}
   =2\pi\,m\Big(1+x+\sqrt{\tfrac{D'}{m}}y\Big)-\pi\log\tfrac{D'}{m}
   =\pi\Big[m\big(m(1+x)^2-D'y^2\big)-\log D'\Big],$$
   verified to $10^{-121}$ for $m=1,2,3,5,7,11,25$ and re-checked independently to 60 digits.
   Since $c_D=-Q_m/a$, the $[n][n/2]$ fold-period module is
   $\mathbf Q+\mathbf Q\frac{\theta_m}{\sqrt{D'}}+\mathbf Q\frac{Q_m}{\sqrt{D'}}$ with the
   last generator a Mahler measure of an **integer** polynomial. For $m=1$ this recovers CDT's
   $4\pi m(1+x+y)-3Q_1-\pi\log3=0$.

6. **[verified exact]** Independent `fractions` re-derivation reproduces the PARI series.
   $H_C$'s coefficients are $-1,-\tfrac{13}{4},-\tfrac{197}{18},-\tfrac{5525}{144},\dots$;
   already $\operatorname{den}(a_2)=4\nmid L(2)L(1)=2$. **The claim in
   arXiv:2408.15403 (Calegari–Dimitrov–Tang), Lemma 2.11.13, that $H_C$ has denominator type
   $[1,\dots,n][1,\dots,n/2]$ is false**; the correct exact type is $[1,\dots,n]^2$.
   The claims for $H_A$ ($1$), $H_B$ ($[n]$) and $H_D$ ($[n][n/2]$) are all correct.
   $\int_0^{1/4}\frac{\log(1-t)}{t\sqrt{1-4t}}dt=-\pi^2/18$ confirmed to 80 digits.
   Consequences for their Theorem 2.11.17 are discussed in §6.3; note their own
   **Question 2.8.3** asks exactly whether Theorem 2.8.4's conclusion survives under the
   hypothesis $a_n[1,\dots,n]^2\in\mathbf Z$ with a fourth puncture $\delta\in[-1/2,1/2]$ —
   which is precisely the type $H_C$ does have, and $\delta=1/4$ is in range.

7. §7 collects the consequence: $H_A,H_B,H_D$ satisfy the $[n][n/2]$ denominator hypothesis for
   **every** $m$ tested, so the layer, its 3-dimensional period module and the fold-regularity
   mechanism all exist verbatim for the whole family $D'=4m-1$, with the third period generator
   equal to $\frac{\pi}{\sqrt{D'}}[\,m(m(1+x)^2-D'y^2)-\log D'\,]$.

---

## 1. Denominator types of $H_A,H_B,H_C,H_D$ (Task 1)

Script `task1.gp`; exact rational power series to $O(x^{402})$, all $n\le 400$.

| $m$ | $D'$ | $H_A$ exc/$1$ | $H_B$ exc/$[n]$ | $H_D$ exc/$[n][n/2]$ | $H_C$ exc/$[n]^2$ |
|---|---|---|---|---|---|
| 1 | 3 | 1 | 1 | 1 | 1 |
| 2 | 7 | 1 | 1 | 1 | 1 |
| 3 | 11 | 1 | 1 | 1 | 1 |
| 4 | 15 | 1 | 1 | 1 | 1 |
| 5 | 19 | 1 | 1 | 1 | 1 |
| 6 | 23 | 1 | 1 | 1 | 1 |
| 7 | 27 | 1 | 1 | 1 | 1 |
| 11 | 43 | 1 | 1 | 1 | 1 |
| 17 | 67 | 1 | 1 | 1 | 1 |
| 41 | 163 | 1 | 1 | 1 | 1 |

**[verified exact, $1\le n\le400$, all ten $m$]**  Types are *exact* (excess 1) and *minimal*:

* $H_A\in\mathbf Z[[x]]$ (excess over $1$ equals $1$).
* $H_B$: excess over $1$ is huge (e.g. $1.97\cdot10^{102}$ at $m=1$) — so $H_B$ is genuinely
  non-integral — and excess over $[n]$ is $1$. Also excess over $[n][n/2]$ is 1 (weaker type).
* $H_D$: excess over $[n]$ is huge ($4.36\cdot10^{49}$ at $m=1$), excess over $[n][n/2]$ is $1$.
* $H_C$: excess over $[n][n/2]$ is huge, excess over $[n]^2$ is $1$.

### $H_C$: unboundedness of the excess over $[n][n/2]$ (m = 1)

$e_n := \operatorname{den}(a_n)/\gcd(\operatorname{den}(a_n),L(n)L(n/2))$, $n=1,\dots,40$:

```
1, 2, 3, 6, 5, 5, 35, 14, 21, 21, 11, 33, 429, 143, 143, 286, 221, 221, 4199, 323,
323, 323, 7429, 7429, 37145, 2185, 6555, 1311, 667, 10005, 310155, 206770, 20677,
20677, 4495, 13485, 498945, 33263, 99789, 99789
```
$n=41,\dots,80$:
```
1363783, 1363783, 58642669, 2022161, 2022161, 2022161, 3065857, 3065857, 21460999,
21460999, 21460999, 21460999, 162490421, 162490421, 162490421, 30741431, 30741431,
30741431, 1813744429, 259106347, 15805487167, 385499687, 2698497809, 5396995618,
62755763, 62755763, 600662303, 600662303, 600662303, 4204636121, 6351684353,
6351684353, 463672957769, 66238993967, 66238993967, 66238993967, 66238993967,
66238993967, 5232880523393, 98733594781
```

**Structure [verified exact, $n\le400$].** $e_n$ is always squarefree-dominated by the primes
just below $n$: writing $P_n:=\prod_{2n/3<p\le n}p$, one has $P_n \mid e_n$ for every $n\le 400$
tested, with a small cofactor $e_n/P_n\in\{1,2,3,5,7,10,13,15,133,357,\dots\}$ coming from prime
powers. Examples: $e_{48}=37\cdot41\cdot43\cdot47$, $e_{43}=29\cdot31\cdot37\cdot41\cdot43$,
$e_{256}=2\cdot3\cdot\prod_{173\le p\le251}p$.

**Growth [numerical, from exact data].** $\log(e_n)/n$:

| $n$ | 20 | 40 | 80 | 120 | 200 | 280 | 320 | 400 |
|---|---|---|---|---|---|---|---|---|
| $\log e_n/n$ | 0.2889 | 0.2878 | 0.3164 | 0.3160 | 0.3717 | 0.3346 | 0.3315 | 0.3310 |

i.e. $e_n=\exp(n/3+o(n))$, consistent with $e_n\approx P_n$ and the prime number theorem
($\log P_n\sim n-\tfrac{2n}{3}=\tfrac n3$). **The excess is unbounded**, and it involves
infinitely many primes, so it cannot be absorbed into any factor $A^{n+1}$ with $A\in\mathbf N$.

### $H_D$ at $m=1$: $\operatorname{den}(a_n)$ against $L(n)L(n/2)$

| $n$ | $\operatorname{den}(a_n)$ | $L(n)L(\lfloor n/2\rfloor)$ | ratio |
|---|---|---|---|
| 1 | 1 | 1 | 1 |
| 2 | 2 | 2 | 1 |
| 3 | 6 | 6 | 1 |
| 4 | 24 | 24 | 1 |
| 5 | 30 | 120 | 4 |
| 6 | 120 | 360 | 3 |
| 7 | 210 | 2520 | 12 |
| 8 | 1120 | 10080 | 9 |
| 9 | 1680 | 30240 | 18 |
| 10 | 6300 | 151200 | 24 |
| 11 | 12600 | 1663200 | 132 |
| 12 | 184800 | 1663200 | 9 |
| 13 | 102960 | 21621600 | 210 |
| 14 | 1261260 | 151351200 | 120 |
| 15 | 1401400 | 151351200 | 108 |
| 16 | 201801600 | 605404800 | 3 |
| 17 | 1715313600 | 10291881600 | 6 |
| 18 | 2333760 | 30875644800 | 13230 |
| 19 | 155195040 | 586637251200 | 3780 |
| 20 | 387987600 | 586637251200 | 1512 |
| 21 | 142942800 | 586637251200 | 4104 |
| 22 | 59750090400 | 6453009763200 | 108 |
| 23 | 22904201320 | 148419224553600 | 6480 |
| 24 | 1099401663360 | 148419224553600 | 135 |
| 25 | 8883336000 | 742096122768000 | 83538 |
| 26 | 25521824328000 | 9647249595984000 | 378 |
| 27 | 1822987452000 | 28941748787952000 | 15876 |
| 28 | 7308522421200 | 28941748787952000 | 3960 |
| 29 | 211947150214800 | 839310714850608000 | 3960 |
| 30 | 706490500716000 | 839310714850608000 | 1188 |

The ratio is always an integer (type $[n][n/2]$ holds) and is $>1$ from $n=5$ on, but it is
*not* systematically large: the type $[n][n/2]$ is exact because for other $n$ the ratio is
small (3 at $n=16$, 1 at $n\le4$) — and, over the tested range $n\le400$, no smaller product of
$L$'s works ($H_D$ is not of type $[n]$: excess $4.36\cdot10^{49}$).

---

## 2. Systematic scan of the $[n][n/2]$ layer (Task 2)

Kernels $k=t^i(1-t)^{-j}\log(1-t)^{\ell}H_A(t)^{2s}$ with
$i\in\{-1,0,1,2\}$ (the value $i=-1$ was **added** to the requested range so that
$H_C$ itself is inside the family; $i=-1$ needs $\ell\ge1$ for $H_Ak$ to be a power series),
$j\in\{0,1,2,3\}$, $\ell\in\{0,1,2\}$, $s\in\{0,1\}$.
Outer factors: $H[k]=H_A\!\int\!H_Ak$, $H'[k]=H_A^3\!\int\!H_Ak$, $H''[k]=\frac{H_A}{1-x}\!\int\!H_Ak$.
Scripts `task2.gp` (N=200, full 3-outer table), `task2b.gp` (N=200, + kernel's own type),
`task2c.gp`/`task2d.gp` (N=400, refinements). $m=1$ and $m=2$.

### 2.1 First structural finding [verified exact, $n\le200$, $m=1,2$]

> **The outer factor is irrelevant.** In all $88\times2\ (i,j,\ell,s)$ kernels $\times\,2$ values of $m$,
> $H[k]$, $H'[k]=H_A^3\!\int\!H_Ak$ and $H''[k]=\frac{H_A}{1-x}\!\int\!H_Ak$ have **identical**
> denominator type. (The only place the three differ at all is the size of a *bounded* excess in
> one row, $(i,j,\ell,s)=(2,0,2,0)$ at $m=1$: excess over $[n]^2$ is $15,5,15$ for
> $H,H',H''$ respectively — same class, different constant.)

Since $H_A^3=H_A\cdot(1-4mx)^{-1}$ and $\frac{H_A}{1-x}$ multiply by an *integral* series, this
is expected; it is recorded because it means the whole scan reduces to the kernel.

### 2.2 The exact ladder [verified exact, $n\le400$ where marked]

Write $\mathcal T_r:=[n][n/2]\cdots[n/r]$ ($\mathcal T_0=1$, $\mathcal T_1=[n]$, $\mathcal T_2=[n][n/2]$,
$\mathcal T_3=[n][n/2][n/3]$).

| kernel $k$ | type of $k$ | **exact** type of $H[k]$ | verified |
|---|---|---|---|
| $1$, $(1-4mt)^{-1}$ | $1$ | $1$ | $n\le200$ |
| $t^i$, $i=1,2$ | $1$ | $[n]$ | $n\le200$ |
| $t^i(1-t)^{-j}$, $j\ge1$ | $1$ | $[n]$ | $n\le200$ |
| $\log(1-t)$ | $[n]$ | $[n]$ | $n\le200$ |
| $t^i\log(1-t)$, $i\ge1$ | $[n]$ | $[n][n/2]$ | $n\le200$ |
| $t^i(1-t)^{-j}\log(1-t)$, $j\ge1$ | $[n]$ | $[n][n/2]$ | $n\le200$ |
| $\log^2(1-t)$ | $[n][n/2]$ | $[n][n/2]$ | $n\le400$ |
| $t^i\log^2(1-t)$, $i=1,2$ | $[n][n/2]$ | $[n][n/2]$ up to bounded factor $3$ ($i=1$), $15$ ($i=2$) | $n\le400$ |
| $t^i(1-t)^{-j}\log^2(1-t)$, $j\ge1$ | $[n][n/2]$ | $[n][n/2][n/3]$ | $n\le400$ |
| $t^{-1}(1-t)^{-j}\log(1-t)$ (= $H_C$ type), $j=0,1,2,3$ | $[n{+}1]$ | $[n]^2$ | $n\le400$ (`task2e.gp`) |
| $t^{-1}(1-t)^{-j}\log^2(1-t)$ | $[n]^2$ | worse than $[n]^2$ (unbounded) | $n\le400$ |

Attaching the factor $H_A(t)^{2}=(1-4mt)^{-1}$ ($s=1$) to a kernel **never changes the ladder
rung**: all 88 $(i,j,\ell)$ rows give the same class for $s=0$ and $s=1$, at both $m=1$ and $m=2$.
It can change a *bounded* excess into $1$ or vice versa in two $m=1$ accidents,
$(i,j,\ell)=(1,0,2)$ (excess $3$ vs $1$) and $(2,2,2)$ (excess $10^{16}$ vs $3$ — at $m=2$ both
are $\approx10^{18}$, so this is not structural).

### 2.3 Answers to the stated expectations

* *"$H[k]$ has type $[n][n/2]$ whenever $k$ has type $[n]$ coefficientwise, e.g.
  $k=\log(1-t)/(1-t)^j$."* — **Confirmed** for $j\ge1$ and for $k=t^i\log(1-t)$, $i\ge1$
  (excess exactly 1, $n\le200$, $m=1,2$). **One exception**: $k=\log(1-t)$ itself
  ($i=j=0$) has type $[n]$ but $H[\log(1-t)]$ has the *better* type $[n]$, not $[n][n/2]$.
  More generally the pure-$\log$ column $(i,j)=(0,0)$ sits one rung *below* the ladder
  ($\ell=0\to1$, $\ell=1\to[n]$, $\ell=2\to[n][n/2]$) — a Catalan-type integrality
  ($\binom{2n}{n}/(n+1)\in\mathbf Z$).
* *"$H[k]$ fails for $k=\log(1-t)/t$."* — **Confirmed**, and sharpened: it does not merely fail
  $[n][n/2]$, it fails $[n][n/2][n/3]$ as well (excess $2.3\cdot10^{55}$ at $n\le400$), and its
  exact type is $[n]^2$.
* *"$\log(1-t)^2$ kernels should fail."* — **Refuted for $j=0$**: $H[\log^2(1-t)]$ has
  **exact type $[n][n/2]$** (excess 1 over $[n][n/2]$, at $n\le400$, $m=1,2$).
  Confirmed for $j\ge1$: those have exact type $[n][n/2][n/3]$, one rung too high.
  (Consistently, $\log^2(1-x)$ itself has exact type $[n][n/2]$ — this is the function CDT
  single out as "the first new function after $\log(1-x)$" in the $\tau=3/2$ layer.)

### 2.4 The layer boundary is sharp

For $k=\log^2(1-t)/(1-t)^j$, $j=1,2,3$ (m = 1 and 2, $n\le400$):
excess over $[n][n/2]$ is $\approx10^{31}$–$10^{34}$ ($m=1$) / $\approx10^{33}$ ($m=2$);
excess over $[n]^2$ is $\approx10^{29}$–$10^{31}$;
excess over $[n]^2[n/2]$ is $1$; **excess over $[n][n/2][n/3]$ is $1$**.
So these functions sit exactly on rung 3 of the ladder.

The $i=-1$ row deserves a remark: $\log(1-t)/(1-t)^j/t$ has coefficients with denominators
requiring exactly $L(n{+}1)$ — the excess over $L(n)$ is exactly $401$ at $N=400$, i.e. the single
prime $n{+}1$ when it is prime. That one-step shift is enough to push $H[k]$ from the ladder rung
$[n][n/2]$ all the way to $[n]^2$ ($j=0,1,2,3$ and $m=1,2$ all give exact type $[n]^2$ at
$n\le400$, with unbounded excess over $[n][n/2]$).

Full machine-generated tables: `task2.out` (216 rows/$m$, three outer factors),
`task2b.out` (88 rows/$m$ with kernel types and excess growth at $n\le50,100,200$),
`task2d.out`.

---

## 3. Fold periods and the $\mathbf Q$-span of the $[n][n/2]$ layer (Task 3)

Scripts `task3.py`, `task3b.py`, `task3c.py` (mpmath, `mp.dps=120`–$150$).

### 3.0 Reduction

With $u=\sqrt{1-4mt}$ ($t=(1-u^2)/(4m)$, $1-t=(D'+u^2)/(4m)$, $H_A=1/u$, $dt=-u\,du/(2m)$):
$$c[k]=\int_0^{\delta}H_A(t)k(t)\,dt=\frac{1}{2m}\int_0^1 k\big(t(u)\big)\,du ,$$
so for $k=t^i(1-t)^{-j}\log^\ell(1-t)$,
$$c(i,j,\ell)=\frac{1}{2m}\int_0^1\Big(\tfrac{1-u^2}{4m}\Big)^{\!i}\Big(\tfrac{D'+u^2}{4m}\Big)^{\!-j}
\Big[\log\tfrac{D'+u^2}{4m}\Big]^{\ell}du .$$

**Note on $s=1$.** The kernels carrying $H_A(t)^2$ have **no finite fold period**:
$H_A\cdot H_A^2 k\sim(1-4mt)^{-3/2}$ is not integrable at $t=\delta$. So the period list below
covers exactly the $s=0$ kernels.

### 3.1 Closed forms verified [numerical, $>110$ digits]

| $m$ | $D'$ | $\lvert c_B-\theta_m/a\rvert$ | $\lvert c_D-(\theta_m\log(D'/m)-2\mathrm{Cl}_2(\pi-\theta_m))/a\rvert$ |
|---|---|---|---|
| 1 | 3 | $9.7\cdot10^{-122}$ | $9.7\cdot10^{-122}$ |
| 2 | 7 | $0$ | $5.5\cdot10^{-122}$ |
| 3 | 11 | $2.4\cdot10^{-122}$ | $2.6\cdot10^{-122}$ |
| 5 | 19 | $0$ | $2.7\cdot10^{-122}$ |
| 11 | 43 | $0$ | $1.4\cdot10^{-122}$ |

(working precision 120 digits; both closed forms confirmed to the full working precision,
well past the requested 50 digits). Values:
$c_B(1)=0.60459978807807261686469275254738524409\ldots$,
$c_D(1)=-0.11708165598778083879284817191973371781\ldots$

### 3.2 The span [numerical, PSLQ at 150 digits, tol $10^{-110}$, maxcoeff $10^{8}$]

$\{1,c_B,c_D\}$ admits **no** integer relation with `maxcoeff` $10^{30}$ at tol $10^{-130}$,
for $m=1,2,3,5,11$: they are a basis of a 3-dimensional space (numerically).

**Result: every fold period of the $\le[n][n/2]$ layer lies in
$\mathbf Q+\mathbf Q\,c_B+\mathbf Q\,c_D$.** No period escapes.
All 27 kernels ($\ell=0$: 12; $\ell=1$: 12; $\ell=2$, $j=0$: 3) were resolved for
$m=1,2,3,5,11$ — full output in `task3b.out`. Sample tables ($m=1$, $D'=3$):

$\ell=0$, $c=q_0+q_1c_B$:

| $(i,j)$ | $(0,0)$ | $(0,1)$ | $(0,2)$ | $(0,3)$ | $(1,0)$ | $(1,1)$ | $(1,2)$ | $(1,3)$ | $(2,0)$ | $(2,1)$ | $(2,2)$ | $(2,3)$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $q_0$ | $1/2$ | $0$ | $1/3$ | $1/2$ | $1/12$ | $-1/2$ | $1/3$ | $1/6$ | $1/60$ | $-7/12$ | $5/6$ | $-1/6$ |
| $q_1$ | $0$ | $1$ | $2/3$ | $2/3$ | $0$ | $1$ | $-1/3$ | $0$ | $0$ | $1$ | $-4/3$ | $1/3$ |

$\ell=1$, $c=q_0+q_1c_B+q_2c_D$:

| $(i,j)$ | $(0,0)$ | $(0,1)$ | $(0,2)$ | $(0,3)$ | $(1,0)$ | $(1,1)$ | $(1,2)$ | $(1,3)$ | $(2,0)$ | $(2,1)$ | $(2,2)$ | $(2,3)$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $q_0$ | $-1$ | $0$ | $1/3$ | $13/36$ | $-17/36$ | $1$ | $1/3$ | $1/36$ | $-497/1800$ | $53/36$ | $-2/3$ | $-11/36$ |
| $q_1$ | $3/2$ | $0$ | $-2/3$ | $-7/9$ | $3/4$ | $-3/2$ | $-2/3$ | $-1/9$ | $9/20$ | $-9/4$ | $5/6$ | $5/9$ |
| $q_2$ | $0$ | $1$ | $2/3$ | $2/3$ | $0$ | $1$ | $-1/3$ | $0$ | $0$ | $1$ | $-4/3$ | $1/3$ |

Note $c(0,1,0)=c_B$ and $c(0,1,1)=c_D$ by definition, and the $\ell=0$ row has $q_2=0$
throughout: **the $\ell=0$ sub-layer spans only $\mathbf Q+\mathbf Q c_B$.**
For general $m$: $c(0,0,0)=\frac1{2m}$, $c(1,0,0)=\frac1{12m^2}$, $c(2,0,0)=\frac1{60m^3}$,
$c(0,2,0)=\frac1{D'}+\frac{2m}{D'}c_B$, and $c(i,j,1)$ has the *same* $c_D$-coefficient as
$c(i,j,0)$ has $c_B$-coefficient (visible in the two tables above; verified for all five $m$).

### 3.3 The new $\log^2$ period is not new

$$\boxed{\;c\big[\log^2(1-t)\big]\;=\;\frac{4-2D'\,c_B+D'\,c_D}{m}\;}$$

[numerical, verified to $<10^{-150}$ absolute for $m=1,2,3,5,11,17$]:

| $m$ | $c[\log^2(1-t)]$ | $\lvert\text{pred}-c\rvert$ |
|---|---|---|
| 1 | $0.0211563035682217824332989689565$ | $3.1\cdot10^{-152}$ |
| 2 | $0.00233190836287103055529303894992$ | $3.3\cdot10^{-151}$ |
| 3 | $0.000664591336825554771210166749876$ | $2.2\cdot10^{-151}$ |
| 5 | $0.000139290467008761745494831902653$ | $3.0\cdot10^{-151}$ |
| 11 | $0.000012770450515714387928528414225$ | $1.4\cdot10^{-151}$ |
| 17 | $0.00000343563994018741604099604094737$ | $6.7\cdot10^{-152}$ |

Equivalently, in terms of $Q_m:=2\mathrm{Cl}_2(\pi-\theta_m)-\theta_m\log(D'/m)$ (so $c_D=-Q_m/a$):
$$c[\log^2(1-t)]=\frac{1}{m}\Big(4-2a\,\theta_m-a\,Q_m\Big),\qquad a=\sqrt{D'} .$$
The periods of $t\log^2(1-t)$ and $t^2\log^2(1-t)$ likewise lie in the span, e.g. at $m=1$
$c(1,0,2)=\frac{124-189c_B+81c_D}{54}$ and $c(2,0,2)=\frac{21107-32535c_B+12150c_D}{13500}$.

### 3.4 Which periods are NOT in $\mathbf Q+\mathbf Q c_B+\mathbf Q c_D$

**Answer: none from the $[n][n/2]$ layer.**  The layer's period module is exactly
$\mathbf Q\oplus\mathbf Q c_B\oplus\mathbf Q c_D$ (dimension 3), so
$\{1,\;\theta_m/a,\;Q_m/a\}$ is a basis.

The escape happens one rung up. For the $[n][n/2][n/3]$-layer kernels
$k=\log^2(1-t)/(1-t)^j$, $j=1,2,3$, PSLQ on $[c,1,c_B,c_D]$ finds **no relation** with
`maxcoeff` $10^{10}$ at tol $10^{-120}$, for $m=1,2,3$ (nine cases). Values at $m=1$:
$c(0,1,2)=0.0270572411888020429117259010361$,
$c(0,2,2)=0.0346806513388123245120244931872$,
$c(0,3,2)=0.0445432180114383572287343927857$.
[numerical, 150 digits — this is evidence of escape, not a proof.]

---

## 4. Fold-regularity (Task 4)

Script `task4.gp`, exact rational coefficients, `realprecision = 900` (needed: the cancellation
is $\sim n\log_{10}(4m)\approx430$ digits at $n=400,m=3$; 80 digits would be far too few).

$|a_n|^{1/n}$ of the indicated series:

**$m=2$** ($4m=8$, $c_B=0.2731678691$, $c_D=-0.0244265737$):

| $n$ | $H_A$ | $H_B$ | $H_B-c_BH_A$ | $H_D$ | $H_D-c_DH_A$ |
|---|---|---|---|---|---|
| 25 | 7.32996 | 6.95919 | **0.81262** | 6.31853 | **0.85739** |
| 50 | 7.60514 | 7.41030 | **0.88924** | 7.06097 | **0.91641** |
| 100 | 7.77318 | 7.67296 | **0.93654** | 7.48992 | **0.95208** |
| 200 | 7.87215 | 7.82123 | **0.96441** | 7.72738 | **0.97299** |
| 300 | 7.90919 | 7.87506 | **0.97482** | 7.81193 | **0.98081** |
| 400 | 7.92895 | 7.90327 | **0.98035** | 7.85571 | **0.98497** |

**$m=3$** ($4m=12$, $c_B=0.1765908357$, $c_D=-0.0102734400$):

| $n$ | $H_A$ | $H_B$ | $H_B-c_BH_A$ | $H_D$ | $H_D-c_DH_A$ |
|---|---|---|---|---|---|
| 25 | 10.99493 | 10.25820 | **0.79809** | 9.15506 | **0.84204** |
| 50 | 11.40771 | 11.01889 | **0.88125** | 10.40957 | **0.90817** |
| 100 | 11.65977 | 11.45934 | **0.93231** | 11.13800 | **0.94779** |
| 200 | 11.80822 | 11.70629 | **0.96224** | 11.54099 | **0.97080** |
| 300 | 11.86379 | 11.79542 | **0.97335** | 11.68412 | **0.97933** |
| 400 | 11.89342 | 11.84198 | **0.97924** | 11.75807 | **0.98386** |

(Also computed for $m=1$ and $m=5$, same picture — see `task4.out`.)

**Conclusion [numerical, 900-digit arithmetic, exact coefficients].** Fold-regularity holds:
the unsubtracted $H_B,H_D$ have $|a_n|^{1/n}\to 4m$ (radius $\delta=1/(4m)$), while the
subtracted series have $|a_n|^{1/n}\to1$. The approach to 1 is slow because the decay is
*polynomial*, not exponential (`task4b.gp`):

| $n$ | $n\lvert a_n(H_B-c_BH_A)\rvert$, $m=2$ | $n\lvert a_n(H_D-c_DH_A)\rvert$, $m=2$ | $n\lvert\cdot\rvert$ $H_B$, $m=3$ | $n\lvert\cdot\rvert$ $H_D$, $m=3$ |
|---|---|---|---|---|
| 50 | 0.141247 | 0.635891 | 0.089930 | 0.404771 |
| 100 | 0.142047 | 0.737049 | 0.090416 | 0.469105 |
| 200 | 0.142450 | 0.837429 | 0.090662 | 0.532955 |
| 300 | 0.142586 | 0.895885 | 0.090744 | 0.570142 |
| 400 | 0.142653 | 0.937274 | 0.090785 | 0.596474 |

$n\lvert a_n\rvert$ is **constant** for $H_B-c_BH_A$ (a $\log(1-x)$ singularity at $x=1$) and
grows like $\log n$ for $H_D-c_DH_A$ (increments $\approx0.100$ per doubling, $=0.145\log2$ —
a $\log^2(1-x)$ singularity). Either way the radius is exactly $1$: the fold at $x=\delta$
has been removed.

---

## 5. Cassaigne–Maillot and Mahler measures attached to $Q_m$ (Task 5)

Scripts `task5a.py`–`task5f.py`, log `task5_results.txt`; independent re-check `task5_check.py`.
Throughout $Q_m:=2\mathrm{Cl}_2(\pi-\theta_m)-\theta_m\log(D'/m)$, so that $c_D=-Q_m/a$.

### 5.1 Verification of the formula [numerical, 56 digits]

Mahler measures were computed by the *exact* Jensen reduction
$$m(a+bx+cy)=\frac{1}{2\pi}\int_0^{2\pi}\log\max\big(|a+be^{is}|,\,c\big)\,ds ,$$
split at $\cos s=(c^2-a^2-b^2)/(2ab)$ so the integrand is smooth on each piece
(mpmath.quad at 55 working digits).

* $L(2,\chi_{-3})=0.7813024128964862968671874296240923563651\ldots$ — direct Dirichlet series,
  Hurwitz-zeta closed form and $\tfrac{4}{3\sqrt3}\mathrm{Cl}_2(\pi/3)$ agree to **60.6–60.8 digits**.
* **Smyth**: $m(1+x+y)=0.3230659472194505140936365107238063940722\ldots$ vs
  $\frac{3\sqrt3}{4\pi}L(2,\chi_{-3})$: **56.5 digits**.
* **Cassaigne–Maillot** tested on 17 triples — triangle cases $(1,1,1)$, $(2,3,4)$, $(3,4,5)$,
  $(5,6,7)$, $(2,2,3)$, $(1,2,2)$, $(7,8,9)$, $(3,5,7)$, $(11,13,17)$, $(1,1,1.5)$,
  $(2.5,3.5,4.5)$, $(1,7,7)$; non-triangle $(1,1,5)$, $(5,1,1)$, $(1,2,10)$, $(2,3,10)$;
  and the degenerate boundary $(1,1,2)$. **Worst-case agreement 56.0 digits**, i.e. limited only
  by the quadrature precision. **The formula is confirmed**, degenerate case included.

### 5.2 Integer-triple search: null [numerical, dps 70, tol $10^{-45}$]

All $31\,079$ triples $1\le a\le b\le c\le60$ with $\gcd=1$, for $m=2,3,5$ ($D'=7,11,19$):

| PSLQ basis | maxcoeff | hits |
|---|---|---|
| $[\pi m(a{+}bx{+}cy),\,Q_m,\,\theta_m,\,\sqrt{D'},\,1]$ | $10^6$ | **0** |
| $[\pi m(a{+}bx{+}cy),\,\theta_m,\,\sqrt{D'},\,1]$ (control, no $Q_m$) | $10^6$ | **0** |
| $[D_{\mathrm{BW}}((b/a)e^{i\gamma}),\,\mathrm{Cl}_2(\pi-\theta_m)]$ | $10^6$ | **0** |
| $[\pi m(a{+}bx{+}cy),\,Q_m,\,\theta_m,\,\sqrt{D'},\,1,\,\pi\log(D'/m)]$ | $10^4$ | **0** |

Non-degeneracy checked first: PSLQ finds no relation among $[Q_m,\theta_m,\sqrt{D'},1]$ nor
among $[\theta_m,\sqrt{D'},1]$ for $m=2,3,5$, so no $M$-free relation could contaminate the scan.
**The null is genuine and its cause is identified**: for an *integer* triple the
Cassaigne–Maillot log-part lies in $\operatorname{span}\{\alpha\log a,\beta\log b,\gamma\log c\}$
and cannot produce $\theta_m\log(D'/m)$ with $D'/m=7/2,\,11/3,\,19/5$.

### 5.3 What does work: an exact identity for every $m$

The right triple is *algebraic*, not integral. Since
$\cos\theta_m=\frac{D'-1}{D'+1}=\frac{2m-1}{2m}$, the isosceles triple
$\big(1,1,\sqrt{D'/m}\big)$ (always a triangle, since $\sqrt{D'/m}<2$) has apex angle
$\gamma$ with $\cos\gamma=1-\tfrac{D'}{2m}=-\tfrac{2m-1}{2m}$, i.e. **$\gamma=\pi-\theta_m$**
exactly. Cassaigne–Maillot with $a=b=1$ then reads
$$\pi\,m\Big(1+x+\sqrt{\tfrac{D'}{m}}\,y\Big)=\tfrac12(\pi-\theta_m)\log\tfrac{D'}{m}+\mathrm{Cl}_2(\pi-\theta_m),$$
hence

$$\boxed{\;Q_m \;=\; 2\pi\,m\Big(1+x+\sqrt{\tfrac{D'}{m}}\,y\Big)\;-\;\pi\log\frac{D'}{m}\;}$$

and, using $2m(1+x+ty)=m\big((1+x)^2-t^2y^2\big)$ and $m(\lambda P)=\log\lambda+m(P)$, the
**integer-polynomial** form

$$\boxed{\;Q_m \;=\; \pi\Big[\,m\big(m(1+x)^2-D'y^2\big)-\log D'\,\Big]\;}$$

**[numerical, 120 digits]** Residuals of the first identity (from `task5_results.txt`):

| $m$ | 1 | 2 | 3 | 5 | 7 | 11 | 25 |
|---|---|---|---|---|---|---|---|
| $D'$ | 3 | 7 | 11 | 19 | 27 | 43 | 99 |
| $Q_m$ | 0.202791376805 | 0.064626639379 | 0.034073145914 | 0.015449534766 | 0.009230422378 | 0.004642311192 | 0.001342738196 |
| residual | $3.9\cdot10^{-121}$ | $-4.8\cdot10^{-121}$ | $0.0$ | $3.9\cdot10^{-121}$ | $-6.8\cdot10^{-121}$ | $1.9\cdot10^{-121}$ | $2.4\cdot10^{-121}$ |

**[independent re-check, `task5_check.py`, 60 digits]** Recomputing the Mahler measure from
scratch (own Jensen-reduced quadrature) gives residuals $\le3.2\cdot10^{-61}$ for
$m=1,2,3,5,7,11$. The Cassaigne–Maillot values were also cross-checked against raw 2-D torus
quadrature (56–57 digits for $m=2,3,5$; a slower direct 2-D quadrature of
$m(m(1+x)^2-D'y^2)$ reproduces the integer form to 5 digits).

The identity is also derivable in closed form (the angle computation above is exact), so this is
**an identity, not a numerical coincidence**.

### 5.4 Consequence

Combining with §3–§4: the fold-period module of the $[n][n/2]$ layer over the host
$(1-4mx)^{-1/2}$ is $\mathbf Q+\mathbf Q\frac{\theta_m}{\sqrt{D'}}+\mathbf Q\frac{Q_m}{\sqrt{D'}}$,
and
$$\frac{Q_m}{\sqrt{D'}}=\frac{\pi}{\sqrt{D'}}\Big[m\big(m(1+x)^2-D'y^2\big)-\log D'\Big].$$
So *"$1$, $\theta_m/\sqrt{D'}$, $Q_m/\sqrt{D'}$ are $\mathbf Q$-linearly independent"* is exactly
a statement about the Mahler measure of the integer polynomial $m(1+x)^2-D'y^2$ together with
$\pi/\sqrt{D'}$ and $\pi\log D'/\sqrt{D'}$.

**$m=1$ control [numerical, 120 digits].** $\theta_1=\pi/3$ exactly,
$\pi m(1+x+y)=\mathrm{Cl}_2(\pi/3)=\frac{3\sqrt3}{4}L(2,\chi_{-3})=1.0149416064096536250212\ldots$,
and PSLQ on $[\pi m(1+x+y),\,Q_1,\,\pi\log3]$ returns $[4,-3,-1]$:
$$4\pi\,m(1+x+y)-3Q_1-\pi\log3=0\qquad(\text{residual }3.9\cdot10^{-121}),$$
i.e. $Q_1=\tfrac43\mathrm{Cl}_2(\pi/3)-\tfrac{\pi}{3}\log3$ and
$\pi m(1+x+y)=\tfrac34Q_1+\tfrac{\pi}{4}\log3$. This is CDT's Smyth-based hit
$m\big((1+x+y)^4/3\big)=4m(1+x+y)-\log3$ — recovered independently, and the general-$m$ identity
of §5.3 specialises to it (for $m=1$ the isosceles triple is $(1,1,\sqrt3)$, and
$\mathrm{Cl}_2(2\pi/3)=\tfrac23\mathrm{Cl}_2(\pi/3)$).
Note PSLQ on the plain 5-vector $[\pi m(1+x+y),Q_1,\theta_1,\sqrt3,1]$ returns **None** at
maxcoeff $10^6$ and $10^4$: $\pi\log3$ must be adjoined. That is why the §5.2 integer scan is
null — it was run without $\pi\log(D'/m)$ in the primary basis (the augmented run at
maxcoeff $10^4$ was also null, but a positive control on the algebraic triple
$(1,1,\sqrt{D'/m})$ recovers $[2,-1,0,0,0,-1]$ in the augmented basis for $m=1,2,3,5$,
confirming the pipeline).

---

## 6. Independent re-derivation and the CDT denominator claim (Task 6)

Script `task6.py` (pure `fractions.Fraction`, no PARI), `task6b.py` (mpmath).

### 6.1 First coefficients ($m=1$) [verified exact]

$$H_C(x)=\frac{1}{\sqrt{1-4x}}\int_0^x\frac{\log(1-t)}{t\sqrt{1-4t}}\,dt
=-x-\frac{13}{4}x^2-\frac{197}{18}x^3-\frac{5525}{144}x^4-\frac{27633}{200}x^5
-\frac{911939}{1800}x^6-\frac{82987349}{44100}x^7-\cdots$$
$$H_D(x)=\frac{1}{\sqrt{1-4x}}\int_0^x\frac{\log(1-t)}{(1-t)\sqrt{1-4t}}\,dt
=-\frac12x^2-\frac{13}{6}x^3-\frac{193}{24}x^4-\frac{881}{30}x^5-\frac{12967}{120}x^6
-\frac{84359}{210}x^7-\cdots$$
These agree exactly with the PARI series of §1–2 (two independent implementations).

Hand check of $a_2(H_C)$: $\frac{\log(1-t)}{t}\cdot\frac1{\sqrt{1-4t}}
=-(1+\tfrac t2+\tfrac{t^2}3)(1+2t+6t^2)=-(1+\tfrac52t+\tfrac{22}3t^2+\cdots)$;
integrating, $-(x+\tfrac54x^2+\cdots)$; multiplying by $1+2x+\cdots$ gives
$a_2=-(\tfrac54+2)=-\tfrac{13}{4}$.

### 6.2 The claim of arXiv:2408.15403

Calegari–Dimitrov–Tang, *The linear independence of $1$, $\zeta(2)$, and $L(2,\chi_{-3})$*,
**Lemma 2.11.13** (§2.11.12 "Some mixed periods") states, for exactly these four functions:

> *"They have the respective denominator types $1$, for $H_A$; $[1,2,\dots,n]$, for $H_B$ and
> $[1,2,\dots,n][1,2,\dots,n/2]$, for $H_C$ and $H_D$."*

Their convention (p. 12): "$[1,\dots,n]$ is used to denote the lowest common multiple of the
first $\lfloor n\rfloor$ positive integers", and the denominator-type condition used in
Theorem 2.5.1 (2.5.3) and in **Theorem 2.8.4** is literally
$[1,\dots,n][1,\dots,n/2]\,a_n\in\mathbf Z$ for all $n$ — no auxiliary $A^{n+1}$ factor.

**Verdict [verified exact].**

* $H_A$: type $1$. **Correct.**
* $H_B$: type $[1,\dots,n]$, exact. **Correct.**
* $H_D$: type $[1,\dots,n][1,\dots,n/2]$, exact. **Correct** ($n\le400$, and for all ten $m$).
* $H_C$: **incorrect.** $H_C$ is *not* of type $[1,\dots,n][1,\dots,n/2]$. It fails at $n=2$
  already: $\operatorname{den}(a_2)=4$ while $[1,2][1,1]=2$. The correct exact type is
  $[1,\dots,n]^2$.

| $n$ | $\operatorname{den}(a_n(H_C))$ | $L(n)L(n/2)$ | divides? | $L(n)^2$ | divides? |
|---|---|---|---|---|---|
| 1 | 1 | 1 | yes | 1 | yes |
| 2 | 4 | 2 | **no** | 4 | yes |
| 3 | 18 | 6 | **no** | 36 | yes |
| 4 | 144 | 24 | **no** | 144 | yes |
| 5 | 200 | 120 | **no** | 3600 | yes |
| 7 | 44100 | 2520 | **no** | 176400 | yes |
| 10 | 3175200 | 151200 | **no** | 6350400 | yes |
| 16 | 19238419200 | 605404800 | **no** | 519437318400 | yes |
| 24 | 14701418922782592 | 148419224553600 | **no** | 28667766899426054400 | yes |

The failure is not an indexing convention: testing $\operatorname{den}(a_n)\mid L(n+s)L((n+s)/2)$
for $s=0,1,2$ fails from $n=2,2,3$ respectively. And it is not absorbable into an $A^{n+1}$:
by §1 the excess is $\prod_{2n/3<p\le n}p\cdot(\text{small})$, which involves infinitely many
primes.

### 6.3 What this does and does not affect

Lemma 2.11.13 feeds **Theorem 2.11.17** ("the four periods $1$, $\pi/\sqrt3$, $\pi^2$,
$3L(2,\chi_{-3})-\frac{\pi}{\sqrt3}\log3$ are $\mathbf Q$-linearly independent"), whose proof
applies part $(*)$ of Theorem 2.8.4 to $f=aH_A+bH_B+cH_C+dH_D$. Since $H_C$ does not satisfy
the denominator hypothesis of Theorem 2.8.4, the argument as written does not cover the
coefficient $c$ of $H_C$.

However, **the corollary survives**: setting $c=0$, the combination $f=aH_A+bH_B+dH_D$ does
satisfy the hypothesis (types $1$, $[n]$, $[n][n/2]$, all verified exact here), and the
overconvergence relation reduces to
$$a+b\frac{\pi}{3\sqrt3}+d\Big(\frac{\pi}{3\sqrt3}\log3-L(2,\chi_{-3})\Big)=0 ,$$
i.e. the $\mathbf Q$-linear independence of $1$, $\pi/\sqrt3$ and
$3L(2,\chi_{-3})-\frac{\pi}{\sqrt3}\log3$ — which is exactly what their Mahler-measure
corollary (2.11.18), $m\big((1+x+y)^4/3\big)\notin\mathbf Q$, needs.
It is the $\pi^2$ entry of Theorem 2.11.17, and only that, which rests on $H_C$.
(Note the monodromy $T(H_C)=H_C+\frac{\pi^2}{9}H_A$ is the only source of $\pi^2$ there.)
This is stated as a fact about the denominators, not as a claim that Theorem 2.11.17 is false:
$H_C$ may still qualify under a different (e.g. $[1,\dots,n]^2$) holonomy bound. Note also that
the paper's main result (Theorem A, linear independence of $1,\zeta(2),L(2,\chi_{-3})$) is proved
in §13 via a Picard–Fuchs construction (§11.1) and does not go through Lemma 2.11.13;
§2.11.12 is explicitly labelled a "showcase application"/"proof-of-concept" in the introduction.
Their own Question 2.8.3 asks precisely whether the conclusion of Theorem 2.8.4 survives under
the weaker hypothesis $a_n[1,\dots,n]^2\in\mathbf Z$ — which is exactly the type $H_C$ has, so
an affirmative answer to Question 2.8.3 would repair Theorem 2.11.17 as stated.

### 6.4 The period $\int_0^{1/4}\frac{\log(1-t)}{t\sqrt{1-4t}}dt$ [numerical, 80 digits]

$$\int_0^{1/4}\frac{\log(1-t)}{t\sqrt{1-4t}}\,dt
=2\int_0^1\frac{\log\frac{3+u^2}{4}}{1-u^2}\,du
=-0.54831135561607547882413838888200839640631663373560\ldots$$
$$-\frac{\pi^2}{18}=-0.54831135561607547882413838888200839640631663373560\ldots$$
Agreement to the full 80-digit working precision (difference $0.0$ at `mp.dps=80`) via the
$u$-substituted integral; the direct $t$-integral agrees to 43 digits (endpoint singularity).
**Confirmed: the period equals $-\pi^2/18$.**  This matches CDT's $T(H_C)=H_C+\frac{\pi^2}{9}H_A$
(the monodromy coefficient is $-2\times$ the fold period).

### 6.5 Cross-validation of our normalisation against CDT (m = 1) [numerical, 120 digits]

The local monodromy at the fold is $T(H[k])=H[k]-2c[k]\,H_A$. CDT's Lemma 2.11.13 gives
$T(H_B)=H_B-2L(1,\chi_{-3})H_A$, $T(H_C)=H_C+\frac{\pi^2}{9}H_A$,
$T(H_D)=H_D-2\big(L(1,\chi_{-3})\log3-L(2,\chi_{-3})\big)H_A$. Our fold periods reproduce all
three exactly:

| | our $c[k]$ at $m=1$ | CDT's coefficient | difference |
|---|---|---|---|
| $c_B$ | $0.604599788078072616864692752547385244094689$ | $L(1,\chi_{-3})=\pi/(3\sqrt3)$ | $0.0$ |
| $c_C$ | $-0.5483113556160754788241383888820083964063$ | $-\pi^2/18$ | $0.0$ |
| $c_D$ | $-0.117081655987780838792848171919733717814186$ | $L(1,\chi_{-3})\log3-L(2,\chi_{-3})$ | $2.4\cdot10^{-121}$ |

($L(2,\chi_{-3})=0.781302412896486296867187429624092356365134$.) So the general-$m$ formulas
$c_B=\theta_m/a$, $c_D=(\theta_m\log(D'/m)-2\mathrm{Cl}_2(\pi-\theta_m))/a$ specialise correctly,
and the series we tested in §6.2 are literally CDT's $H_A,\dots,H_D$.

**Bonus (not identified).** The escapee $c[\log^2(1-t)/(1-t)]$ at $m=1$ equals
$0.02705724118880204291172590103609236932192\ldots$; PSLQ over an 18-element weight-$\le3$,
level-6 basis $\{1,\pi/\sqrt3,\mathrm{Cl}_2(\pi/3)/\sqrt3,\zeta(3),\pi^2\log3,\log^33,
\pi\,\mathrm{Cl}_2(\pi/3),\pi^3,\pi\log^23,\ldots\}$ returned only relations with
5-digit coefficients at 120 working digits, i.e. numerically unresolved. Not identified here.

---

## 7. Synthesis: the same three functions for every $m$

$H_A$, $H_B$, $H_D$ have exact types $1$, $[n]$, $[n][n/2]$ for **every** $m$ tested
($m\in\{1,2,3,4,5,6,7,11,17,41\}$, $n\le400$ — §1), and their fold at $\delta=1/(4m)$ lies in
$(-\infty,1)$, and in $[-1/2,1/2]$ for every $m\ge1$. So the *denominator* hypothesis of CDT's
Theorem 2.8.4 (and of their Question 2.8.3) is met by the whole family, not just $m=1$ —
whereas $H_C$ meets it for **no** $m$ (§1, §6).

The overconvergence / fold-regularity condition for $f=aH_A+bH_B+dH_D$ at the fold is
$$a+b\,c_B+d\,c_D=0,\qquad\text{i.e.}\qquad
a+b\,\frac{\theta_m}{\sqrt{D'}}-d\,\frac{Q_m}{\sqrt{D'}}=0 ,$$
which §4 verifies numerically is exactly when the fold is removed and the radius becomes 1;
§3 shows $\{1,c_B,c_D\}$ is the whole period module of the layer (nothing else to test);
and §5 identifies
$$\frac{Q_m}{\sqrt{D'}}=\frac{\pi}{\sqrt{D'}}\Big[\,m\big(m(1+x)^2-D'y^2\big)-\log D'\,\Big].$$

Whether the remaining (analytic/holonomy) hypotheses of the CDT argument transfer to general $m$
was **not** checked here; only the denominators, the periods and the fold-regularity were.

---

## Files

| file | contents |
|---|---|
| `task1.gp`, `task1.out` | Task 1: types of $H_A,H_B,H_C,H_D$, ten values of $m$, $n\le400$ |
| `task1b.gp`, `task1b.out` | $H_C$ excess growth, $\prod_{2n/3<p\le n}p$ structure |
| `task2.gp`, `task2.out` | Task 2 full table, 3 outer factors $\times$ 88 kernels $\times$ $m=1,2$, $N=200$ |
| `task2b.gp`, `task2b.out` | same + kernel's own type + excess at $n\le50,100,200$ |
| `task2c.gp`/`task2d.gp`, `.out` | $N=400$ refinements; the $[n][n/2][n/3]$ rung |
| `task2e.gp`, `task2e.out` | $N=400$: the $t^{-1}$ ($H_C$-type) kernels, all $j$ |
| `task3.py`, `task3b.py`, `task3c.py`, `.out` | Task 3: closed forms, 27 periods, PSLQ, span |
| `task4.gp`, `task4.out`, `task4b.gp`, `task4b.out` | Task 4: fold-regularity, 900-digit |
| `task5a.py`–`task5f.py`, `task5_results.txt` | Task 5: Cassaigne–Maillot, integer scan, the identity |
| `task5_check.py` | independent 60-digit re-check of the §5.3 identity |
| `task6.py`, `task6b.py`, `task6c.py`, `.out` | Task 6: `fractions` re-derivation, $-\pi^2/18$, CDT cross-check |
| `cdt.pdf`, `cdt.txt` | arXiv:2408.15403 (fetched, for quoting Lemma 2.11.13 verbatim) |
