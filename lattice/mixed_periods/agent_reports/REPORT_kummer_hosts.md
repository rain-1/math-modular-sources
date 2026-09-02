# Kummer-type hosts $H=(1\mp Nx)^{-1/k}$: integrality, fold periods, closed forms

All computations in `/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/kummer`.
Tools: PARI/GP 2.15.4 (exact), python3 + mpmath 1.3.0 + sympy 1.14.0 (numerical, PSLQ).

Claim tags: **[verified exact, range]** = exact rational arithmetic over the stated range;
**[numerical, d digits]** = floating identity confirmed at $d$ significant digits;
**[PSLQ, m/d/c]** = integer relation found (or excluded) with basis size $m$, tolerance $10^{-d}$, coefficient bound $c$.

---

## 0. Conventions, and two sign corrections to the statement

Fix $k\ge 2$ and $N\ge 1$. Two families:

| family | host | fold $\delta$ | $D$ | $Q(u)$ | $\sigma$ | roots $\alpha_j$ of $Q$ |
|---|---|---|---|---|---|---|
| **M** ("minus") | $H=(1-Nx)^{-1/k}$ | $+1/N$ | $N-1$ | $u^k+D$ | $+1$ | $D^{1/k}e^{i\pi(2j+1)/k}$ |
| **P** ("plus", real) | $H=(1+Nx)^{-1/k}$ | $-1/N$ | $N+1$ | $u^k-D$ | $-1$ | $D^{1/k}e^{2\pi i j/k}$ |

so that $Q(u)=u^k+\sigma D$ and $\alpha_j^k=-\sigma D$.

Substituting $u=(1\mp Nt)^{1/k}$ gives, in **both** families,
$$1-t=\frac{|Q(u)|}{N},\qquad t=1-\frac{|Q(u)|}{N},\qquad dt=-\sigma\,\frac{k}{N}\,u^{k-1}du ,$$
and hence, uniformly,
$$\boxed{\;c^{(a)}[\tfrac{1}{1-t}]=k\!\int_0^1\!\frac{u^{a-1}\,du}{Q(u)},\quad
c^{(a)}[\tfrac{\log(1-t)}{1-t}]=k\!\int_0^1\!\frac{u^{a-1}\log\!\big(|Q(u)|/N\big)}{Q(u)}du,\quad
c^{(a)}[\log(1-t)]=\sigma\frac{k}{N}\!\int_0^1\! u^{a-1}\log\frac{|Q(u)|}{N}\,du. \;}$$

**Sign correction 1.** For the real family the prompt's recipe "replace $D+u^k$ by $\widetilde D-u^k$" is right for the *integrand* but the period acquires an **overall minus sign**:
$$c^{(a)}[\tfrac1{1-t}] \;=\; -\,k\int_0^1\frac{u^{a-1}}{\widetilde D-u^k}\,du,\qquad \widetilde D=N+1,$$
and likewise for the other two kernels (for $\log(1-t)$ the sign is the $\sigma$ above). Consequently **every P-family period of the kernels $1/(1-t)$ and $\log(1-t)/(1-t)$ is negative**, whereas the M-family $c^{(a)}[1/(1-t)]$ are positive. *[numerical, 20 digits: direct $t$-quadrature over $[0,\delta]$ vs. the $u$-formulas agree to $\ge 20$ digits on all 15 test hosts $k=2..6$, both families; the residual loss is the algebraic endpoint singularity $(1\mp Nt)^{-(k-a)/k}$ at $t=\delta$.]*

**Sign correction 2** (to the Task 8 host list): for $k=5$ the two cyclotomic-collapse hosts are
$(1-1025x)^{-1/5}$ (M-family, $D=N-1=1024=4^5$) and $(1+7775x)^{-1/5}$ (P-family, $D=N+1=7776=6^5$) —
the signs given in the task were interchanged. All six other hosts in the list are correct as stated.

Throughout, $\rho:=D^{1/k}$, and the three kernels are abbreviated
$$B=\tfrac1{1-t},\qquad \mathcal D=\tfrac{\log(1-t)}{1-t},\qquad L=\log(1-t).$$

---

## 1. Task 1 — integrality and denominator types at $N=k^km$

Denominator type: a series $\sum a_nx^n$ has type $[1..n]$ if $L(n)a_n\in\mathbb Z$ with $L(n)=\operatorname{lcm}(1,\dots,n)$;
type $[1..n][1..\lfloor n/2\rfloor]$ if $L(n)L(\lfloor n/2\rfloor)a_n\in\mathbb Z$. *Excess* $:=\max_n \operatorname{den}(a_n)/\gcd(\operatorname{den}(a_n),T(n))$; excess $1$ $\iff$ the type holds.

**Result [verified exact, $0\le n\le 300$].** For every $k\in\{2,\dots,7\}$, $m\in\{1,2,3\}$, both families:

* $H^j\in\mathbb Z[[x]]$ for all $1\le j\le k-1$ — **no** non-integral coefficient in $n\le 300$;
* $H^{(a)}[B]$ has type $[1..n]$, **excess $=1$**, all sectors $a=1,\dots,k-1$;
* $H^{(a)}[\mathcal D]$ and $H^{(a)}[L]$ have type $[1..n][1..\lfloor n/2\rfloor]$, **excess $=1$**, all sectors.

(252 series; script `task1.gp`, output `task1.out`.)

### 1a. The tests have teeth (controls) *[verified exact, $n\le 120$]*

| control | outcome |
|---|---|
| scaling $N=m$ instead of $k^km$ | $H^j\notin\mathbb Z[[x]]$, first failure at $n=1$ (all $k=3,4,5$) |
| scaling $N=k^km-1$ | $H^j\notin\mathbb Z[[x]]$, first failure at $n=1$ |
| $H^{(a)}[B]$ against type $[1]$ (i.e. integrality) | fails hugely, e.g. $k=3,m=1,a=1$: raw denominator excess $4.8\cdot10^{29}$ at $n=118$ — so $L(n)$ is genuinely needed |
| $H^{(a)}[\mathcal D]$ against type $[1..n]$ alone | fails hugely, e.g. $k=4,m=2,a=1$: excess $2.29\cdot10^{16}$ at $n=118$ — so $L(\lfloor n/2\rfloor)$ is genuinely needed |

### 1b. A sharpening: the $\log(1-t)$ kernel is *almost* of type $[1..n]$

*[verified exact, $n\le 200$, $k\le6$, $m\le3$, both families]* The excess of $H^{(a)}[\log(1-t)]$ against type $[1..n]$ **alone** is either $1$ or exactly $a$, and when it is $a$ the maximum is attained at $n=2a$:

| $a$ | max excess vs $[1..n]$ | attained at |
|---|---|---|
| 1 | 1 | — |
| 2 | 2 | $n=4$ |
| 3 | 3 | $n=6$ |
| 4 | 4 | $n=8$ |
| 5 | 5 | $n=10$ |

Whether excess $=a$ or $=1$ depends on $m$ (e.g. $k=3,a=2$: excess $2$ for $m=1,3$, excess $1$ for $m=2$; $k=6,a=5$: excess $5$ for all $m$). The bad indices $n$ form a sparse set (e.g. $k=3,m=1,a=2$: exactly $n\in\{4,8,16,32,64,128\}$). At $n=2a$ the extra prime power is precisely the one supplied by $L(\lfloor n/2\rfloor)=L(a)$, so **type $[1..n][1..\lfloor n/2\rfloor]$ is sharp but only barely** for this kernel, in contrast with $\mathcal D$ where the second factor is needed massively.

---

## 7. Task 7 — the true integrality condition is $k\cdot\mathrm{rad}(k)\mid N$

**Result [verified exact, $n\le200$].** For $k=2,\dots,7$, $m=1,\dots,4$, both signs, $N=k\,\mathrm{rad}(k)\,m$:
$(1\mp Nx)^{-j/k}\in\mathbb Z[[x]]$ for all $1\le j\le k-1$. Here $k\cdot\mathrm{rad}(k)=4,9,8,25,36,49$ for $k=2,\dots,7$.

**Sharpness [verified exact, $1\le N\le400$, $n\le80$, $k=2,\dots,7$, both signs].** The set of $N\ge1$ for which $(1\mp Nx)^{-j/k}\in\mathbb Z[[x]]$ for all $j=1,\dots,k-1$ is **exactly the set of multiples of $k\cdot\mathrm{rad}(k)$** — verified as an exact set equality on $[1,400]$ in all twelve cases.

Named counterexamples (first non-integral coefficient index, $j=1,\dots,k-1$):

| host | first bad $n$ per $j$ |
|---|---|
| $(1-3x)^{-j/3}$, $(1+3x)^{-j/3}$, $(1-6x)^{-j/3}$ | $[3,3]$ |
| $(1-4x)^{-j/4}$, $(1+4x)^{-j/4}$ | $[2,-1,2]$ ($-1$ = integral; only $j=2$ survives) |
| $(1-2x)^{-j/4}$ | $[1,2,1]$ |
| $(1-5x)^{-j/5}$ | $[5,5,5,5]$ |
| $(1-6x)^{-j/6}$ | $[2,3,2,3,2]$ |
| $(1-12x)^{-j/6}$ | $[3,3,-1,3,3]$ |
| $(1-18x)^{-j/6}$ | $[2,-1,2,-1,2]$ |
| $(1-7x)^{-j/7}$ | $[7,7,7,7,7,7]$ |

Note $N=k^km$ is the special case $m'=k^{k-1}m/\mathrm{rad}(k)$; the enlarged family is genuinely bigger ($k=3$: $N=9m$ vs $27m$; $k=4$: $8m$ vs $256m$).

**Denominator types persist on the enlarged family [verified exact, $n\le200$]:** for $k=3,\dots,6$, $N=k\,\mathrm{rad}(k)\,m$ with $m=1,2,3$, both families, all sectors — excess $1$ for $B$ against $[1..n]$ and for $\mathcal D,L$ against $[1..n][1..\lfloor n/2\rfloor]$.

---

## 2. Task 2 — fold periods, 60 digits

Full table (all $k\in\{3,4,5,6\}$, $m\in\{1,2,3\}$, $N=k^km$, both families, all sectors, kernels $B,\mathcal D,L$; 252 values at 60 significant digits): **`periods60.txt`**.

*[numerical, 80 digits]* Every entry was computed two ways — the $u$-substitution quadrature and the exact closed forms of §3–§4 — and the **maximum relative discrepancy over the whole table is $2.1\times10^{-80}$** (working precision `mp.dps=90`). Independently, direct $t$-quadrature over $[0,\delta]$ reproduces the same values to $\ge20$ digits (§0). Every period is real: $|\Im|<10^{-70}$ throughout.

Sample ($k=3$, $m=1$, $N=27$):

| family | $D$ | $a$ | $c^{(a)}[B]$ | $c^{(a)}[\mathcal D]$ | $c^{(a)}[L]$ |
|---|---|---|---|---|---|
| M | 26 | 1 | $0.1142988941743341896127030\ldots$ | $-0.0032396217092159236520929\ldots$ | $-0.0031365279408123411188579\ldots$ |
| M | 26 | 2 | $0.0568254903139172976626384\ldots$ | $-0.0012879847109140443293953\ldots$ | $-0.0012520695465639033761890\ldots$ |
| P | 28 | 1 | $-0.1081195147990716557157316\ldots$ | $-0.0029452110536450965049828\ldots$ | $-0.0030384904860007066711651\ldots$ |
| P | 28 | 2 | $-0.0543542739556373193342919\ldots$ | $-0.0011849664160306684722689\ldots$ | $-0.0012177594865469411866763\ldots$ |

---

## 3. Task 3 — closed form of $c^{(a)}[1/(1-t)]$

### 3a. General formula

Partial fractions with $Q'(\alpha_j)=k\alpha_j^{k-1}=-\sigma kD/\alpha_j$ give $\dfrac{k}{Q(u)}=\sum_j\dfrac{-\sigma\alpha_j/D}{u-\alpha_j}$. Writing $u^{a-1}=(u-\alpha_j)q_j(u)+\alpha_j^{a-1}$, the polynomial parts contribute $\sum_j\alpha_j^{\,i+1}$ with $1\le i+1\le a-1\le k-2$, and **all these power sums of the roots of $u^k+\sigma D$ vanish**. Hence

$$\boxed{\;c^{(a)}[\tfrac1{1-t}]\;=\;-\frac{\sigma}{D}\sum_{j=0}^{k-1}\alpha_j^{\,a}\,\log\!\Big(1-\frac1{\alpha_j}\Big)\;}\qquad(\text{principal branches}).$$

*[numerical, 111 digits]* Verified against $u$-quadrature on all 15 test hosts ($k=2,\dots,6$, both families, all sectors): max relative error $2.8\times10^{-111}$ at `mp.dps=120`. **No branch corrections are needed** — the principal branch of $\log(1-\alpha_j)-\log(-\alpha_j)$ is correct for every root.

Sanity check $k=2$, $N=4$, $D=3$, $a=1$: $\alpha_\pm=\pm i\sqrt3$, giving $\pi/(3\sqrt3)$. ✓

### 3b. Explicit real forms, $k=3$

Put $\rho=D^{1/3}$.

**M-family** ($D=N-1$), with $A=\log\dfrac{\rho+1}{\sqrt{\rho^2-\rho+1}}$, $T=\sqrt3\,\arctan\dfrac{\sqrt3}{2\rho-1}$:
$$c^{(1)}[\tfrac1{1-t}]=\frac{\rho}{D}\,(A+T),\qquad c^{(2)}[\tfrac1{1-t}]=\frac{\rho^2}{D}\,(T-A).$$

**P-family** ($D=N+1$), with $A'=\log\dfrac{\rho-1}{\sqrt{\rho^2+\rho+1}}$, $T'=\sqrt3\,\arctan\dfrac{\sqrt3}{2\rho+1}$:
$$c^{(1)}[\tfrac1{1-t}]=\frac{\rho}{D}\,(A'-T'),\qquad c^{(2)}[\tfrac1{1-t}]=\frac{\rho^2}{D}\,(A'+T').$$

(The identity $\tfrac\pi3-\arctan\frac{\rho\sqrt3}{\rho+2}=\arctan\frac{\sqrt3}{2\rho+1}$ was used to remove a stray $\pi$.)

### 3c. Explicit real forms, $k=4$

Put $\rho=D^{1/4}$.

**M-family**, with $\Lambda=\tfrac12\log\dfrac{\rho^2+\sqrt2\rho+1}{\rho^2-\sqrt2\rho+1}$ and $\Theta=\arctan\dfrac{\sqrt2\rho}{\rho^2-1}$:
$$c^{(1)}=\frac{\sqrt2\,\rho}{D}(\Lambda+\Theta),\qquad
c^{(2)}=\frac{2\rho^2}{D}\arctan\frac1{\rho^2}=\frac{2}{\sqrt D}\arctan\frac1{\sqrt D},\qquad
c^{(3)}=\frac{\sqrt2\,\rho^3}{D}(\Theta-\Lambda).$$

**P-family**:
$$c^{(1)}=\frac{\rho}{D}\Big(\log\frac{\rho-1}{\rho+1}-2\arctan\frac1\rho\Big),\qquad
c^{(2)}=\frac{\rho^2}{D}\log\frac{\rho^2-1}{\rho^2+1},\qquad
c^{(3)}=\frac{\rho^3}{D}\Big(\log\frac{\rho-1}{\rho+1}+2\arctan\frac1\rho\Big).$$

*[numerical, 79 digits]* All 64 instances ($k=3$: $N=9,18,27,54,81,63,126$; $k=4$: $N=8,16,80,256,512,624$; both families, all sectors) verified against §3a: worst relative error $4.2\times10^{-79}$.

Note the middle sector of $k=4$ M is a pure arctangent, $c^{(2)}=\frac{2}{\sqrt D}\arctan\frac1{\sqrt D}$, i.e. the $k=2$ period at $D$ — the $a=k/2$ sector degenerates.

---

## 4. Task 4 — closed form of $c^{(a)}[\log(1-t)/(1-t)]$

This was derived symbolically (not guessed), then verified numerically.

Since $\prod_l(u-\alpha_l)=Q(u)$ and $\log(1-t)=\log|Q(u)|-\log N$, one has on $u\in[0,1]$
$$\log(1-t)=\sum_{l}\operatorname{Log}(u-\alpha_l)-\log N-i\pi w,\qquad w:=\frac1\pi\sum_l\arg(u-\alpha_l)\in\{0,\pm1\}\ \text{(constant)},$$
with $w=0$ for the M family and $w=1$ for the P family *[numerical, checked at $u=1/2$, all hosts]*. Exactly as in §3a, the polynomial parts die by vanishing power sums, leaving

$$\boxed{\;c^{(a)}[\tfrac{\log(1-t)}{1-t}]\;=\;-\frac{\sigma}{D}\sum_{j,l=0}^{k-1}\alpha_j^{\,a}\,I(j,l)\;-\;\big(\log N+i\pi w\big)\,c^{(a)}[\tfrac1{1-t}]\;}$$

where $I(j,l)=\displaystyle\int_0^1\frac{\operatorname{Log}(u-\alpha_l)}{u-\alpha_j}\,du$ is given by

$$I(j,j)=\tfrac12\Big(\operatorname{Log}^2(1-\alpha_j)-\operatorname{Log}^2(-\alpha_j)\Big),$$
$$I(j,l)=\operatorname{Log}(1-\alpha_l)\operatorname{Log}(1-P_1)+\operatorname{Li}_2(P_1)-\operatorname{Log}(-\alpha_l)\operatorname{Log}(1-P_0)-\operatorname{Li}_2(P_0)\quad(j\ne l),$$
$$P_1=P_1(l,j)=\frac{1-\alpha_l}{\alpha_j-\alpha_l},\qquad P_0=P_0(l,j)=\frac{-\alpha_l}{\alpha_j-\alpha_l}=\frac{\alpha_l}{\alpha_l-\alpha_j}.$$

The antiderivative used is $F(u)=\operatorname{Log}(u-\alpha_l)\operatorname{Log}\frac{\alpha_j-u}{\alpha_j-\alpha_l}+\operatorname{Li}_2\frac{u-\alpha_l}{\alpha_j-\alpha_l}$, with $F'(u)=\operatorname{Log}(u-\alpha_l)/(u-\alpha_j)$.

Likewise, elementarily (no dilogarithms),
$$c^{(a)}[\log(1-t)]=\frac{\sigma k}{N}\cdot\frac1a\left[\sum_l\Big(\operatorname{Log}(1-\alpha_l)-\sum_{i=0}^{a-1}\frac{\alpha_l^{\,i}}{a-i}-\alpha_l^{\,a}\big(\operatorname{Log}(1-\alpha_l)-\operatorname{Log}(-\alpha_l)\big)\Big)-\log N-i\pi w\right].$$

**Verification [numerical, 111 digits].** All three closed forms were checked against $u$-quadrature on 15 hosts spanning $k=2,\dots,6$, both families, all sectors, including all the Task 4 targets ($k=3$, $m=1,2,3$; $k=4$, $m=1,2$): worst relative error $2.8\times10^{-111}$ at `mp.dps=120`. **All principal branches; no $2\pi i$ corrections are required anywhere.**

The dilogarithm arguments are exactly the shapes anticipated: $\alpha_l/(\alpha_l-\alpha_j)$ and $(\alpha_l-1)/(\alpha_l-\alpha_j)$. The coefficient of $\operatorname{Li}_2(P_0(l,j))$ is $+\sigma\alpha_j^a/D$ and of $\operatorname{Li}_2(P_1(l,j))$ is $-\sigma\alpha_j^a/D$ — written out, $\frac{1}{k\alpha_j^{k-1}}=-\frac{\sigma\alpha_j}{kD}$ as required.

---

## 5. Task 5 — the period span of a sector

**Result [numerical/PSLQ, 140 digits, coefficient bound $10^{12}$; $k=3$: $N=9,27$ M, $27,63$ P; $k=4$: $N=8,80,256$ M/P; all sectors].**
All five kernels reduce to $\{1,\,c^{(a)}[B],\,c^{(a)}[\mathcal D]\}$, and the relations are *exact and elementary*, not arithmetic coincidences. PSLQ found them, and they then admit one-line proofs:

1. $\dfrac{t}{1-t}=\dfrac1{1-t}-1$ and $\displaystyle\int_0^\delta H^{k-a}dt=\frac{\sigma k}{Na}$, hence
$$c^{(a)}[\tfrac{t}{1-t}]=c^{(a)}[\tfrac1{1-t}]-\frac{\sigma k}{Na}.$$
2. $\dfrac{\log(1-t)}{1-t}=\log(1-t)+\dfrac{t\log(1-t)}{1-t}$, hence
$$c^{(a)}[\tfrac{t\log(1-t)}{1-t}]=c^{(a)}[\tfrac{\log(1-t)}{1-t}]-c^{(a)}[\log(1-t)].$$
3. Integration by parts in $u$ (the boundary term vanishes because $|Q(1)|=N$ in both families), followed by $u^{a+k-1}/Q=u^{a-1}-\sigma D u^{a-1}/Q$:
$$\boxed{\;c^{(a)}[\log(1-t)]=\frac{kD}{Na}\,c^{(a)}[\tfrac1{1-t}]\;-\;\frac{\sigma k^2}{Na^2}.\;}$$

*[numerical, 140 digits]* Identity 3 reproduces every PSLQ relation found, e.g. $k=4,N=256$ M: $a=1$: $64L=255B-4$; $a=2$: $128L=255B-2$; $a=3$: $576L=765B-4$; and $k=4,N=256$ P: $a=1$: $64L=257B+4$.

**Conclusion.** Exactly as for $k=2$: for every $k$, host and sector,
$$\mathbb Q\text{-span}\big\{c^{(a)}[\kappa]:\kappa\in\{1,\tfrac1{1-t},\tfrac{\log(1-t)}{1-t},\log(1-t),\tfrac t{1-t},\tfrac{t\log(1-t)}{1-t}\}\big\}\;=\;\mathbb Q+\mathbb Q\,c_B^{(a)}+\mathbb Q\,c_{\mathcal D}^{(a)} .$$
In particular $c^{(a)}[\log(1-t)]$ needs **no** $c_{\mathcal D}$ at all: it lies in $\mathbb Q+\mathbb Q\,c_B^{(a)}$.

**Cross-sector [PSLQ, 3 terms / 140 digits / $10^{12}$]: no relation.** For every host tested and every pair $a\ne b$, neither $c^{(b)}[B]\in\mathbb Q+\mathbb Q\,c^{(a)}[B]$ nor $c^{(b)}[\mathcal D]\in\mathbb Q+\mathbb Q\,c^{(a)}[\mathcal D]$ (32 tests, all negative). E.g. for $k=3$ M the two sectors are $\frac{\rho}{D}(A+T)$ and $\frac{\rho^2}{D}(T-A)$, which span the same $\mathbb Q(\rho)$-plane $\langle A,T\rangle$ but are $\mathbb Q$-affinely independent. **Different sectors of the same host are arithmetically independent.**

---

## 8. Task 8 — cyclotomic collapse

### 8a. The complete list of collapse hosts

$D$ is a perfect $k$-th power $n^k$ exactly when $N=n^k+1$ (M) or $N=n^k-1$ (P), subject to $k\,\mathrm{rad}(k)\mid N$. *[verified exact, $N\le3\cdot10^5$]*

| $k$ | condition | hosts (smallest) |
|---|---|---|
| 3 | $n\equiv2\ (3)$ → M, $n\equiv1\ (3)$ → P | $(1-9x)^{-1/3}$ $[8=2^3]$, $(1+63x)^{-1/3}$ $[64=4^3]$, $(1-126x)^{-1/3}$ $[125=5^3]$, $(1+342x)$ $[7^3]$, $(1-513x)$ $[8^3]$, … |
| 4 | **P only** ($n^4\equiv1\ (8)$ for odd $n$; $n^4\equiv-1\ (8)$ impossible) | $(1+80x)^{-1/4}$ $[81=3^4]$, $(1+624x)^{-1/4}$ $[625=5^4]$, $(1+2400x)$ $[7^4]$, $(1+6560x)$ $[9^4]$, … |
| 5 | $n\equiv-1\ (5)$ → M, $n\equiv1\ (5)$ → P | $(1-1025x)^{-1/5}$ $[1024=4^5]$, $(1+7775x)^{-1/5}$ $[7776=6^5]$, $(1-59050x)$ $[9^5]$, $(1+161050x)$ $[11^5]$ |
| 6 | **P only** | $(1+15624x)^{-1/6}$ $[15625=5^6]$, $(1+117648x)$ $[7^6]$ |

Then $\alpha_j=n\zeta_j$ with $\zeta_j\in\mu_{2k}\setminus\mu_k$ (M) or $\zeta_j\in\mu_k$ (P), and the dilogarithm points split into two families:
$$P_0(l,j)=\frac{\zeta_l}{\zeta_l-\zeta_j}=\frac{1}{1-\mu},\ \ \mu=\zeta_j/\zeta_l \quad(\textbf{$n$-independent, purely cyclotomic}),\qquad
P_1(l,j)=\frac{n\zeta_l-1}{n(\zeta_l-\zeta_j)} .$$

For $k=4$ (P) the $P_0$ points are $\{\tfrac{1\pm i}2,\ \tfrac12\}$ with $D(\tfrac{1\pm i}2)=\pm G$ (Catalan);
for $k=3$ they are $\{\tfrac12\pm\tfrac{i}{2\sqrt3}\}$ with $D=\pm\tfrac23\mathrm{Cl}_2(\pi/3)$.

### 8b. **The cyclotomic points cancel identically** — so $G$ does *not* enter through them

$P_0$ depends only on $\mu=\zeta_j/\zeta_l$; the pairs $(l,j)$ realising a given $\mu$ are $j=$ (rotation of) $l$, and the aggregate coefficient is
$$\frac{\sigma}{D}\sum_{l}(n\zeta_l\mu)^a=\frac{\sigma n^a\mu^a}{D}\sum_{l}\zeta_l^{\,a}=0\qquad(1\le a\le k-1),$$
because $\sum_l\zeta_l^a$ is a power sum of a full coset of $\mu_k$.

**[numerical, 320 digits]** Confirmed: for every host in the table and every sector, the aggregate $\operatorname{Li}_2$-coefficient of every $P_0$ point is $0$ to $\le 5\times10^{-321}$ (hosts $k=3$ $N=9,63$; $k=4$ $N=80,624$; $k=5$ $N=1025$; $k=6$ $N=15624$; 19 sector-instances).

> **Answer to the key question: for the $k=4$ hosts $(1+80x)^{-1/4}$ and $(1+624x)^{-1/4}$, Catalan's constant $G$ does NOT appear.** Its coefficient is exactly $0$ in every sector — the $G$-point $\frac{1+i}2$ is a $P_0$ point and cancels. Nor does $G$ reappear through the $P_1$ points: the surviving Bloch–Wigner values are not rational multiples of $G$ *[PSLQ, 2-term and 5-term, 280 digits, bound $10^8$: the only relation found among $\{D(z_i)\}\cup\{G\}$ is the tautology $D(\frac{1-i}2)=-G$]*.

Likewise for $k=3$: the $\mathrm{Cl}_2(\pi/3)$ coefficient coming from $P_0$ vanishes. It can only survive when the **$P_1$ points happen to collapse**, which occurs for exactly one host in our list:

> **$(1-9x)^{-1/3}$, $D=8=2^3$.** Here $D(P_1)=\pm\frac56\mathrm{Cl}_2(\pi/3)$, and the total Bloch part is
> $$\text{Bloch}\big(c^{(1)}[\tfrac{\log(1-t)}{1-t}]\big)=-\tfrac{5}{24}\sqrt3\,\mathrm{Cl}_2(\pi/3),\qquad
> \text{Bloch}\big(c^{(2)}\big)=-\tfrac{5}{12}\sqrt3\,\mathrm{Cl}_2(\pi/3),$$
> i.e. $-\tfrac{5}{24}\cdot\tfrac{3\sqrt3}{4}\sqrt3\,L(2,\chi_{-3})=-\tfrac{15}{32}L(2,\chi_{-3})$ resp. $-\tfrac{15}{16}L(2,\chi_{-3})$.
> *[PSLQ, 2-term, 300 digits, bound $10^{12}$]* **So $\sqrt3\,L(2,\chi_{-3})$ does appear, for this host only.**
> For $(1+63x)^{-1/3}$ and $(1-126x)^{-1/3}$ the $\mathrm{Cl}_2(\pi/3)$ coefficient is $0$ *[PSLQ, 300 digits]*.

### 8c. Exact closed forms

Reducing modulo $\operatorname{Li}_2(z)+\operatorname{Li}_2(1-z)=\frac{\pi^2}6-\log z\log(1-z)$ and $z\mapsto\bar z$ (the point set is closed under both), each period becomes

$$c^{(a)}\big[\tfrac{\log(1-t)}{1-t}\big]\;=\;\sum_{z}\Big(\beta_z\,\Re\operatorname{Li}_2(z)+\gamma_z\,D(z)\Big)\;+\;\mathcal E,$$

with $\mathcal E$ a $\mathbb Q$- (resp. $\mathbb Q(\sqrt3)$-) combination of products of $\pi$, $\log p$ and one extra angle $\theta$. **All coefficients below are exact rationals** ($\beta,\gamma$ read off symbolically; $\mathcal E$ by PSLQ). *[PSLQ, coefficient bound $10^6$. Decisive (relation length $\times$ bound $\ll$ tolerance): $(1-9x)^{-1/3}$ basis 12 at 276 digits; $(1+80x)^{-1/4}$ basis 15 at 270; $(1+624x)^{-1/4}$ basis 21 at 258; $(1+63x)^{-1/3}$ basis 30 at 240. **Marginal**: $(1-126x)^{-1/3}$, basis 42 at 216 digits ($42\times6=252>216$) — the relation there is corroborated by the small denominators ($\le225$) and by the exact proportionality of the $a=1$ and $a=2$ coefficient vectors, but is not certified. Independently, the assembled formulas for $(1-9x)^{-1/3}$ ($a=1,2$) and $(1+80x)^{-1/4}$ ($a=1,2$) were re-evaluated from the printed rationals and agree with §4 to $\ge61$ digits.]*

**Host $(1-9x)^{-1/3}$, $D=8$** (points $z_1=\frac34-\frac{i\sqrt3}4$ and $\frac12$; angle basis $\{\pi\}$; primes $\{2,3\}$):
$$c^{(1)}=-\tfrac34\Re\operatorname{Li}_2(z_1)-\tfrac14\Re\operatorname{Li}_2(\tfrac12)+\tfrac{\sqrt3}4 D(z_1)
-\tfrac12\log^22+\tfrac38\log2\log3-\tfrac1{32}\log^23+\tfrac{7}{96}\pi^2+\tfrac{\sqrt3}8\pi\log2-\tfrac{\sqrt3}{48}\pi\log3,$$
$$c^{(2)}=\tfrac32\Re\operatorname{Li}_2(z_1)+\tfrac12\Re\operatorname{Li}_2(\tfrac12)+\tfrac{\sqrt3}2 D(z_1)
+\log^22-\tfrac34\log2\log3+\tfrac1{16}\log^23-\tfrac{7}{48}\pi^2+\tfrac{\sqrt3}4\pi\log2-\tfrac{\sqrt3}{24}\pi\log3,$$
with $D(z_1)=-\tfrac56\mathrm{Cl}_2(\pi/3)$ (§8b). *[numerical, 62 digits: both verified against §4 to $3\times10^{-61}$.]*

**Host $(1+80x)^{-1/4}$, $D=81=3^4$.** Points $w_1=\frac23-\frac i3$, $w_2=\frac23$, $w_3=\frac23-\frac{2i}3$, $w_4=\frac12-\frac i6$; primes $\{2,3,5\}$; angle $\theta=\arg(2-i)=-\arctan\frac12$.

| $a$ | $\beta_{w_1}$ | $\beta_{w_2}$ | $\beta_{w_3}$ | $\beta_{w_4}$ | $\gamma_{w_1}$ | $\gamma_{w_3}$ | $\gamma_{w_4}$ |
|---|---|---|---|---|---|---|---|
| 1 | $2/27$ | $2/27$ | $2/27$ | $0$ | $2/27$ | $-2/27$ | $-2/27$ |
| 2 | $4/9$ | $0$ | $-4/9$ | $-2/9$ | $0$ | $0$ | $0$ |
| 3 | $2/3$ | $2/3$ | $2/3$ | $0$ | $-2/3$ | $2/3$ | $2/3$ |

$$\mathcal E_1=\tfrac1{54}\log^22-\tfrac5{27}\log2\log3+\tfrac1{18}\log2\log5+\tfrac19\log^23-\tfrac1{27}\log3\log5
+\tfrac1{54}\pi\log2-\tfrac2{27}\pi\log3+\tfrac1{54}\pi\log5+\tfrac2{27}\theta\log2-\tfrac8{27}\theta\log3+\tfrac2{27}\theta\log5-\tfrac1{54}\pi^2-\tfrac1{54}\pi\theta,$$
$$\mathcal E_2=-\tfrac14\log^22+\tfrac59\log2\log3-\tfrac1{18}\log2\log5-\tfrac19\log^23-\tfrac19\log3\log5+\tfrac1{36}\log^25-\tfrac7{432}\pi^2-\tfrac16\pi\theta-\tfrac19\theta^2,$$
$$\mathcal E_3=\tfrac16\log^22-\tfrac53\log2\log3+\tfrac12\log2\log5+\log^23-\tfrac13\log3\log5
-\tfrac16\pi\log2+\tfrac23\pi\log3-\tfrac16\pi\log5-\tfrac23\theta\log2+\tfrac83\theta\log3-\tfrac23\theta\log5-\tfrac16\pi^2-\tfrac16\pi\theta.$$

*[numerical, 62 digits: $a=1$ and $a=2$ reassembled from the printed coefficients agree with §4 to $5\times10^{-63}$ and $1\times10^{-61}$.]*

**Host $(1+624x)^{-1/4}$, $D=625=5^4$**: same shape with $w_1=\frac35-\frac{2i}5$, $w_2=\frac35$, $w_3=\frac35-\frac{3i}5$, $w_4=\frac12-\frac i{10}$, primes $\{2,3,5,13\}$, $\theta=\arg(3-2i)$, and $\beta,\gamma$ obtained from the $D=81$ table by $27\mapsto125$, $9\mapsto25$, $3\mapsto5$. Full data: `t8g.out`.

**Hosts $(1+63x)^{-1/3}$ ($D=64$) and $(1-126x)^{-1/3}$ ($D=125$)**: two surviving points each,
$\beta\in\frac1{16}\mathbb Z$ resp. $\frac1{25}\mathbb Z$ and $\gamma\in\frac{\sqrt3}{16}\mathbb Z$ resp. $\frac{\sqrt3}{25}\mathbb Z$; e.g.
$$c^{(1)}_{D=64}=\tfrac3{16}\Re\operatorname{Li}_2\big(\tfrac58-\tfrac{i\sqrt3}8\big)-\tfrac1{16}\Re\operatorname{Li}_2\big(\tfrac12-\tfrac{i\sqrt3}4\big)
+\tfrac{\sqrt3}{16}D\big(\tfrac58-\tfrac{i\sqrt3}8\big)-\tfrac{\sqrt3}{16}D\big(\tfrac12-\tfrac{i\sqrt3}4\big)+\mathcal E,$$
$\mathcal E=\tfrac14\log^22-\tfrac18\log2\log7-\tfrac1{128}\log^23+\tfrac1{64}\log3\log7+\tfrac1{128}\log^27-\tfrac{19}{1152}\pi^2-\tfrac5{96}\pi\theta-\tfrac1{32}\theta^2-\tfrac{\sqrt3}{16}\pi\log2-\tfrac{3\sqrt3}8\theta\log2+\tfrac{\sqrt3}{192}\pi\log3+\tfrac{\sqrt3}{32}\theta\log3+\tfrac{\sqrt3}{96}\pi\log7+\tfrac{\sqrt3}{16}\theta\log7$,
with $\theta=\arg(\tfrac58-\tfrac{i\sqrt3}8)$. Full data: `t8g.out`.

**Surviving non-Bloch-trivial points (the honest answer to "which points survive"):**

| host | surviving $z$ with $\gamma_z\ne0$ | $D(z)$ | rational multiple of $G$ / $\mathrm{Cl}_2(\pi/3)$? |
|---|---|---|---|
| $(1-9x)^{-1/3}$ | $\frac34-\frac{i\sqrt3}4$ | $-0.845784672008\ldots$ | **yes**, $-\frac56\mathrm{Cl}_2(\pi/3)$ |
| $(1+63x)^{-1/3}$ | $\frac58-\frac{i\sqrt3}8$, $\frac12-\frac{i\sqrt3}4$ | $-0.5473880894$, $-0.8613851257$ | **no** |
| $(1-126x)^{-1/3}$ | $\frac35-\frac{i\sqrt3}5$, $\frac12+\frac{i\sqrt3}{10}$ | $-0.7627027250$, $+0.4497903024$ | **no** |
| $(1+80x)^{-1/4}$ | $\frac23-\frac i3$, $\frac23-\frac{2i}3$, $\frac12-\frac i6$ | $-0.7447789113$, $-0.9797850741$, $-0.4348208089$ | **no** |
| $(1+624x)^{-1/4}$ | $\frac35-\frac{2i}5$, $\frac35-\frac{3i}5$, $\frac12-\frac i{10}$ | $-0.8264930725$, $-0.9659534070$, $-0.2710744469$ | **no** |

**Two structural remarks.**
* For $k=4$, **sector $a=2$ has no Bloch–Wigner content at all**: $\gamma_z=0$ for every $z$, in both hosts. (Consistent with §3c, where $c^{(2)}[B]$ also degenerates.)
* The proposed shape "period $=$ (rational)$\cdot G+\sum$(rational)$D(z_i)+\sum$(rational)(log products)$+$(rational)$\pi^2$" is **false**: there is an irreducible $\Re\operatorname{Li}_2$ part. *[PSLQ, negative, 300 digits, bound $10^{10}$: for $(1-9x)^{-1/3}$ (basis 12), $(1+80x)^{-1/4}$ (basis 15) and $(1+624x)^{-1/4}$ (basis 21), neither $c_{\mathcal D}-\text{Bloch}$ nor the $\Re\operatorname{Li}_2$ part alone lies in the span of $\{\pi^2,\log p\log q,\pi\log p,\theta\log p,\pi\theta,\theta^2\}$ (and its $\sqrt3$-multiples).]* The correct shape is the boxed one of §8c: **Bloch part $+$ $\Re\operatorname{Li}_2$ part $+$ elementary**.

---

## 6. Task 6 — Mahler measures

$\pi\,m(1+x+y)=\mathrm{Cl}_2(\pi/3)=\frac{3\sqrt3}{4}L(2,\chi_{-3})$ *[numerical, 25 digits, control confirmed]*.

**Hit.** For $(1-9x)^{-1/3}$ ($D=8$):
$$\text{Bloch}\big(c^{(1)}[\tfrac{\log(1-t)}{1-t}]\big)=-\tfrac5{24}\sqrt3\;\pi\,m(1+x+y),\qquad
\text{Bloch}\big(c^{(2)}\big)=-\tfrac5{12}\sqrt3\;\pi\,m(1+x+y).$$
*[PSLQ, 2-term, 300 digits, bound $10^{10}$]* — the coefficient is $\sqrt3$-rational, not rational, which is why a plain rational search misses it.

**The general Cassaigne–Maillot picture (positive, and better than a search).** Cassaigne–Maillot: for a genuine triangle with sides $a,b,c$ and opposite angles $\alpha,\beta,\gamma$,
$\pi m(a+bx+cy)=D\!\big(\tfrac ba e^{i\gamma}\big)+\alpha\log a+\beta\log b+\gamma\log c$.

**Every surviving Bloch point of every collapse host is the Cassaigne–Maillot point of the triangle with vertices $1,\ n\zeta_l,\ n\zeta_j$ in $\mathbb C$**, i.e. sides
$$\big(a,b,c\big)=\big(\,|n\zeta_l-n\zeta_j|,\ |n\zeta_l-1|,\ |n\zeta_j-1|\,\big).$$
*[numerical, 60 digits]* Verified for all five hosts $k=3,4$: in each case $D(z_{\text{period}})=\pm D(z_{\text{CM}})$ exactly (sign recorded below), for every $(l,j)$ whose triangle is non-degenerate.

| host | triangle sides | CM point | $\pi m(a+bx+cy)$ |
|---|---|---|---|
| $(1-9x)^{-1/3}$, $n=2$ | $(2\sqrt3,\ 3,\ \sqrt3)$ | $\frac34+\frac{i\sqrt3}4$ | $4.235505913975314\ldots$ |
| $(1+63x)^{-1/3}$, $n=4$ | $(4\sqrt3,\ \sqrt{21},\ 3)$ | $\frac58+\frac{i\sqrt3}8$ | $6.132719357462022\ldots$ |
| | $(4\sqrt3,\ \sqrt{21},\ \sqrt{21})$ | $\frac12+\frac{i\sqrt3}4$ | $6.352232802679864\ldots$ |
| $(1-126x)^{-1/3}$, $n=5$ | $(5\sqrt3,\ 6,\ \sqrt{21})$ | $\frac35+\frac{i\sqrt3}5$ | $6.949409713257435\ldots$ |
| | $(5\sqrt3,\ \sqrt{21},\ \sqrt{21})$ | $\frac12+\frac{i\sqrt3}{10}$ | $6.807184870069179\ldots$ |
| $(1+80x)^{-1/4}$, $n=3$ | $(3\sqrt2,\ \sqrt{10},\ 2)$ | $\frac23+\frac i3$ | $4.705460023765574\ldots$ |
| | $(3\sqrt2,\ 4,\ \sqrt{10})$ | $\frac23+\frac{2i}3$ | $5.223945462488639\ldots$ |
| | $(6,\ \sqrt{10},\ \sqrt{10})$ | $\frac12+\frac i6$ | $5.651658019465864\ldots$ |
| $(1+624x)^{-1/4}$, $n=5$ | $(5\sqrt2,\ \sqrt{26},\ 4)$ | $\frac35+\frac{2i}5$ | $6.379692953946477\ldots$ |
| | $(5\sqrt2,\ 6,\ \sqrt{26})$ | $\frac35+\frac{3i}5$ | $6.692722583426045\ldots$ |
| | $(10,\ \sqrt{26},\ \sqrt{26})$ | $\frac12+\frac i{10}$ | $7.238952502414742\ldots$ |

So **the entire Bloch–Wigner content of these periods is Mahler-measure content**: substituting the table into §8c expresses $\gamma_z D(z)$ as $\pm\gamma_z\big(\pi m(a+bx+cy)-\alpha\log a-\beta\log b-\gamma\log c\big)$, the log terms merging into $\mathcal E$.

**Null results.** A direct PSLQ search over triangles with sides drawn from $\{1,2,3,4,\rho,\rho^2,2\rho,3\rho,\rho\pm1,\rho/2,\sqrt\rho\}$, $\rho=D^{1/3}$, against $\{\pi m(P),\pi\log2,\pi\log3,\pi\log\rho,\pi^2\}$ found **no** hit for any $k=3$ host at 120 digits, bound $10^6$ — the reason being visible above: the correct sides are $|n\zeta_l-n\zeta_j|,|n\zeta_l-1|,|n\zeta_j-1|$ (which involve $\sqrt3\,n$, $\sqrt{n^2+n+1}$, …), not powers of $D^{1/3}$. The control $P=1+x+y$ is a hit only after $\sqrt3$-rational coefficients are allowed.

---

## Files

| file | contents |
|---|---|
| `task1.gp` / `task1.out` | Task 1, $n\le300$ |
| `ctrl.gp` / `ctrl.out` | Task 1 controls and sharpness |
| `sharp.gp`, `sharp2.gp` | §1b, the $\log(1-t)$ kernel |
| `task7.gp`,`task7b.gp` / `.out` | Task 7 |
| `collapse.gp` | §8a host enumeration |
| `kummer.py` | host class, quadratures, closed forms §3a §4 |
| `bloch.py`, `decomp.py`, `ident2.py` | Bloch–Wigner, exact decomposition, PSLQ bases |
| `task2.py` / `periods60.txt` | Task 2, 60-digit tables |
| `task3.py` | Task 3 explicit real forms |
| `task5.py` / `task5.out` | Task 5 |
| `t8e.py`, `t8f.py`, `t8g.py` / `t8g.out` | Task 8 |
| `task6b.py` | Task 6 |
