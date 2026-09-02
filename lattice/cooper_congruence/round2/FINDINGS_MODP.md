# The mod-$p$ structure of Cooper's kernel $\gamma$ / $K$

*Round 2, working note, 2026-09-02.  All computations in PARI/GP 2.15.4, exact over
$\mathbf Z$ (no floating point anywhere except one $\chi^2$ statistic).  Scripts
`20_*.gp` … `29_*.gp` with matching `.log` files in this directory.  Background:
`consolidation/COOPER_CONGRUENCE.md`, `lattice/cooper_congruence/REPORT.md` (round 1),
`lattice/cooper_sources/REPORT.md` §4.3.  Nothing outside
`lattice/cooper_congruence/round2/` was created or modified.*

Objects, as in round 1: $\Phi=\sum c(m)q^m$ (Cooper's rows $s_7,s_{10},s_{18}$, levels
$N=7,10,18$, characters $\psi=\mathbf1,\mathbf1,\chi_{-3}$), $c'(m)=c(m)/m$,
$\beta=c'\star(\mu\psi)$, $\gamma(n)=\beta(n)/n^2\in\mathbf Z$,
$$K(q)=\sum_{n\ge1}\gamma(n)q^n,\qquad \Xi(q)=\sum_m c'(m)q^m,\qquad B(q)=\sum_n\beta(n)q^n
=\textstyle\sum_n n^2\gamma(n)q^n .$$
Throughout, **$P$ denotes the row's distinguished odd prime**
$$P_{s_7}=7,\qquad P_{s_{10}}=5,\qquad P_{s_{18}}=3$$
— for $s_7,s_{10}$ the prime of the two proved Atkin–Lehner cells
($\Phi|U_P=P\,\Phi$, round 1 Thm 2.2 / V14); for $s_{18}$ the prime with $\psi(3)=0$.

Data recomputed here from scratch, exact: $c',\beta,\gamma$ for $n\le 12000$
(`20_extend12000.gp`, files `20c_*.txt`; also `20b_*` at $6000$, `22_*` at $3000$,
`20_*` at $1500$).  **Bonus:** $n^2\mid\beta(n)$ (round-1 conjecture V4 / Conj. 4.1)
now **[verified, $n\le12000$]** for all three rows, up from $n\le1500$.

---

## 0. Verdict table

Tags: **[proved]**, **[verified, range]**, **[refuted, range]**, **[conjectural]**.

| # | statement | status |
|---|---|---|
| **M1** | **Mod-2 law.**  $\gamma(n)\equiv1\pmod 2\iff 2^{1+v_2(N)}\nmid n$ and $P\nmid n$.  Explicitly: $s_7$: $\gcd(n,14)=1$; $s_{10}$: $4\nmid n$ and $5\nmid n$; $s_{18}$: $4\nmid n$ and $3\nmid n$ | **[verified, $n\le12000$]** |
| **M2** | Hence $\gamma\bmod 2$ is a **multiplicative** function of $n$ (although $\gamma$ is not), and $K\bmod 2$ is a **rational** function of $q$ with denominator $1-q^{T}$, $T=14,20,12$ | **[verified, $n\le12000$; $0/6506$ coprime-pair failures]** |
| **M3** | **Mod-3 law for $s_{18}$.**  $(-1)^{n-1}\gamma_{s_{18}}(n)\equiv1\pmod 3$ if $3\nmid n$, $\equiv-1\pmod 3$ if $3\mid n$; equivalently $K_{s_{18}}\equiv(q+2q^2+2q^3)/(1+q)^3\pmod 3$, period $6$ | **[verified, $n\le12000$]** |
| **M4** | These four cells $\bigl((s_7,2),(s_{10},2),(s_{18},2),(s_{18},3)\bigr)$ are the **only** ones: for every other $(\text{row},p)$ with $p\le43$, $K\bmod p$ is **not** a rational function of $q$ of degree $\le250$ and **not** algebraic in the boxes $(\deg_q,\deg_K)=(40,6)$, $(20,20)$, $(60,4)$ | **[refuted, checked against $1500$–$3000$ coefficients]** |
| **M5** | $\gamma(n)\bmod p$ is **not** eventually periodic for any $(\text{row},p)$ outside those four cells | **[refuted, offset $\le40$, period $\le400$, $n\le1500$]** |
| **M6** | The 2-adic refinement stops at finite depth: $\gamma\bmod4$ is periodic (period $40$, $24$) for $s_{10},s_{18}$ but **not** for $s_7$; $\gamma\bmod 8$ is periodic for none.  $\gamma_{s_{18}}\bmod 9$ is not periodic, but **is** periodic (period $18$) on $\{n:3\nmid n\}$ | **[verified/refuted, $n\le12000$, period $\le1000$]** |
| **M7** | $\gamma_{s_{18}}(n)\gamma_{s_{18}}(n+1)\equiv-1\pmod 9$ for every $n\equiv1\pmod3$; sharp (fails mod $27$ at $n=7$) | **[verified, $n\le12000$]** |
| **M8** | **Zeros.**  $\gamma_{s_7}(n)=0\iff7\mid n$; $\gamma_{s_{10}}(n)=0\iff5\mid n$; $\gamma_{s_{18}}(n)\ne0$ always | **[proved]** ($\Leftarrow$, §5.1) + **[verified, $n\le12000$]** ($\Rightarrow$) |
| **M9** | **Extra valuation laws.**  $v_2(\gamma_{s_7}(n))=[2\mid n]$ exactly, for $7\nmid n$; $v_2(\gamma_{s_{10}}(n))\ge v_2(n)-1$ (sharp); $v_3(\gamma_{s_7}(n))\ge v_3(n)-1$ (sharp); $3\nmid\gamma_{s_{18}}(n)$ always; $v_2(\gamma_{s_{18}}(n))$ depends only on $n\bmod24$ | **[verified, $n\le12000$]** |
| **M10** | **$\Xi\bmod2$ is Eisenstein.**  $c'(m)$ is odd $\iff$ ($s_7$) the $14$-free part of $m$ is a perfect square; ($s_{10}$) the $10$-free part of $m$ is a square; ($s_{18}$) $3\nmid m$ and the odd part of $m$ is a square.  In particular $c'\bmod2$ is multiplicative | **[verified, $m\le6000$; $0/14102$ coprime-pair failures]** |
| **N1** | $\gamma(p)\bmod p$ (= round-1's $\lambda(p)=(c'(p)-\psi(p))/p^2$): **no arithmetic formula found**.  Tested against the trace forms of $S_k^{\rm new}(\Gamma_0(M))$ for $k=2,4,6$, $M\le80$; every rational eigenform for $k=2,4,6$, $M\le60$; $\eta(3t)^8$, $\eta(t)^4\eta(7t)^4$, $\eta(t)^{12}\eta(2t)^{12}$, $\tau$; $\chi_{-3},\chi_{-4},\chi_{-7},\chi_{-8},\chi_5$; $B_{p-3}$, $B_{(p-1)/2}$, Fermat quotients $q_p(2),q_p(3),q_p(5),q_p(7)$, the Wilson quotient, $\binom{p-1}{(p-1)/2}$, $h(-p)$, $h(-4p)$, $H_{(p-1)/2}$, $\sum_{j\le(p-1)/2}j^{-2}$, $p\bmod 5,7,9,12$ — all four of $=,\,=-,\,=(\cdot)^{-1},\,=-(\cdot)^{-1}$, over the $42$ primes $11\le p\le199$.  Best score anywhere: $5/42$, i.e. exactly noise | **[refuted for this candidate list]** |
| **N2** | $\gamma(p)\equiv(a_{p-1}-\psi(p))/p\pmod p$ (the natural $x$-side guess from round-1 V9) | **[refuted, $\le2/42$]** |
| **N3** | $\gamma(p)\bmod p$ is equidistributed: $\chi^2=1.17,\ 4.40,\ 0.42$ on $4$ bins over the $426$ primes $11\le p\le3000$ (5% critical value $7.81$) | **[verified]** |
| **N4** | $p\mid\gamma(p)$ happens for $s_7$ at $p=2,5,7,\mathbf{167}$; for $s_{10}$ at $p=5$; for $s_{18}$ never, $p\le3000$.  ($\sum_{p\le3000}1/p=2.34$, so these counts are Poisson-normal; $p=167$ is new but not structural) | **[verified, $p\le3000$]** |
| **L1** | Lucas $\gamma(pn)\equiv\gamma(p)\gamma(n)\pmod p$ | **[refuted]**, massively — e.g. $s_{18}$, $p=11$: $124/136$ fail.  Holds only in the cells where $\gamma(p)\equiv0$, where it degenerates to $\gamma(pn)\equiv0$ |
| **L2** | Lucas $\gamma(pn+r)\equiv\gamma(n)\gamma(r)\pmod p$, $1\le r<p$ | **[refuted]**, massively (e.g. $s_7$, $p=13$: $1213/1373$ fail) |
| **L3** | $\gamma(pn)\equiv\psi(p)\gamma(n)\pmod p$ (and mod $p^2$) | **[refuted]** (e.g. $s_{18}$, $p=3$: $500/500$ fail).  Round-1's $(S)$ does **not** propagate to $\gamma$ |
| **L4** | $\gamma(pn)/\gamma(n)\bmod p$ constant in $n$ | **[refuted]**: the number of distinct ratios is $p$ or $p-1$ for $p\ge5$ (all values occur) |
| **L5** | $K(q)\equiv A(q)K(q^p)+B(q)\pmod p$, $\deg A,\deg B\le60$ | exists **exactly** in the four cells of M1/M3 (where it follows from rationality); **[refuted]** for all other $(\text{row},p)$, $p\le13$, against $1500$ coefficients |
| **L6** | The $p$-kernel $\{n\mapsto\gamma(p^sn+r)\bmod p\}$ is finite (Christol) | finite in the four cells; **[refuted]** elsewhere — at $p\ge5$ all $p+1$ depth-$1$ subsequences are already pairwise distinct, at $p=3$ all $13$ depth-$2$ ones are |
| **E1** | The Lambert exponents $e_n=\beta(n)/n=n\gamma(n)$: $e_n\bmod2$ is periodic ($14,10,6$), $e_n\bmod4$ periodic for $s_{10},s_{18}$ ($20,12$), $e_n\bmod9$ periodic for $s_{18}$ ($18$); nothing beyond, and nothing at $p\ge5$ | **[verified/refuted, $n\le3000$, period $\le600$]** |
| **E2** | The $q$-product $\mathfrak g=\prod_n(1-q^n)^{-e_n}=\exp\bigl(\sum_mc'(m)q^m/m\bigr)$ is integral for $s_7,s_{10}$ | **[verified, $q^{1200}$]** (round 1 had $q^{80}$) |
| **E3** | The coefficients of $\mathfrak g$ are periodic mod $2,3,5$ | **[refuted, period $\le400$, $1200$ coefficients]** |
| **X1** | $\Xi\bmod p$ or $B\bmod p$ algebraic/rational | $B$ behaves exactly like $K$ (it is $D^2K$); $\Xi$ is **[refuted]** at *every* $p\le43$ including $p=2,3$ — consistent with M10, since the indicator of the squares is not automatic |

**Headline.**  There is a genuine, previously unrecorded mod-$p$ structure, and it lives
entirely at $p=2$ (all three rows) and at $p=3$ (the row where $\psi(3)=0$): there
$\gamma$ becomes a *multiplicative* function and $K$ becomes a *rational* function of $q$,
i.e. the master conjecture's kernel degenerates to a genuine Eisenstein series mod $p$.
At every other prime the structure is exactly nil: $K\bmod p$ is not rational, not
algebraic, not automatic, satisfies no Lucas or Dwork law, and $\gamma(p)\bmod p$ is
equidistributed and matches no modular or elementary invariant we could think of.

---

## 1. The mod-2 law  [verified, $n\le12000$]

`25_ratform.gp`, `29_prod.gp`, `21_laws6000.gp`, `21_laws12000.gp`, `22_cprime2.gp`.

**Observation 1.1.**  $\gamma(n)\bmod2$ is purely periodic in $n$, with periods
$14,20,12$ for $s_7,s_{10},s_{18}$; one period is
$$\begin{aligned}
s_7&: [1,0,1,0,1,0,0,0,1,0,1,0,1,0]\\
s_{10}&: [1,1,1,0,0,1,1,0,1,0,1,0,1,1,0,0,1,1,1,0]\\
s_{18}&: [1,1,0,0,1,0,1,0,0,1,1,0]
\end{aligned}$$
and the pattern is exactly the indicator of a **multiplicative** condition:
$$\boxed{\ \gamma(n)\equiv1\pmod2\iff 2^{\,1+v_2(N)}\nmid n\ \text{ and }\ P\nmid n\ }$$
($N=7,10,18$; $2^{1+v_2(N)}=2,4,4$; $P=7,5,3$).  Verified for every $n\le12000$, all
three rows.  Consequently $\gamma\bmod2$ is multiplicative — $0$ failures of
$\gamma(ab)\equiv\gamma(a)\gamma(b)\pmod 2$ over $6506$ coprime pairs with $a\le60$,
$ab\le3000$ (`21_extra.log`) — even though $\gamma$ itself is not multiplicative at all.

**Corollary 1.2 (the rational form).**  With $\mathcal L(t)=t/(1-t)$,
$$K\equiv\mathcal L(q)-\mathcal L(q^{2^{1+v_2(N)}})-\mathcal L(q^{P})+\mathcal L(q^{2^{1+v_2(N)}P})
\pmod2 .$$
Independently found by Hermite–Padé (`25_ratform.log`) as the minimal
$P/Q$, verified to $q^{12000}$ (`26_dbg.gp`, `26_dbg2.gp`):
$$K_{s_7}\equiv\frac{q+q^5+q^7+q^{11}}{1+q^2+q^4+q^6+q^8+q^{10}+q^{12}}
=\frac{q+q^3+q^5+q^9+q^{11}+q^{13}}{1-q^{14}},$$
$$K_{s_{10}}\equiv\frac{q^{15}+q^{14}+q^{13}+q^{11}+q^5+q^3+q^2+q}{1+q^4+q^8+q^{12}+q^{16}},
\qquad
K_{s_{18}}\equiv\frac{q+q^2+q^6+q^7}{1+q^4+q^8}\pmod2 .$$
(Minimal degrees $12,16,8$.)

**Corollary 1.3 ($\Xi\bmod2$ is Eisenstein).**  Mod $2$, $d^2\equiv d$ and
$\psi(k)\in\{0,\pm1\}$, so $c'(m)=\sum_{d\mid m}\psi(m/d)d^2\gamma(d)$ collapses to a
divisor count.  Explicitly **[verified, $m\le6000$, `22_cprime3.log`]**
$$c'_{s_7}(m)\ \text{odd}\iff \tfrac{m}{2^{v_2}7^{v_7}}\ \text{is a square},\qquad
c'_{s_{10}}(m)\ \text{odd}\iff \tfrac{m}{2^{v_2}5^{v_5}}\ \text{is a square},$$
$$c'_{s_{18}}(m)\ \text{odd}\iff 3\nmid m\ \text{and}\ \tfrac{m}{2^{v_2}}\ \text{is a square}.$$
Equivalently $\Xi_{s_7}\equiv\sum_{a,b\ge0}\Theta(q^{2^a7^b})\pmod2$ with
$\Theta(t)=\sum_{\gcd(s,14)=1}t^{s^2}$, and similarly for the other rows.  $c'\bmod 2$ is
therefore multiplicative ($0/14102$ coprime-pair failures).  This is the mod-$2$ shadow of
Conjecture 4.1: the conjecture says $\Xi$ has the *shape* of a weight-three Eisenstein
series with a non-multiplicative kernel; **mod $2$ the kernel becomes multiplicative and
$\Xi$ is honestly Eisenstein**.  It also explains X1: $\Xi\bmod2$ is a theta-type series
(indicator of squares), which is not a $2$-automatic set, hence not algebraic — and indeed
the Hermite–Padé boxes find nothing for $\Xi$ at any $p$.

## 2. The mod-3 law for $s_{18}$  [verified, $n\le12000$]

Only the row with $\psi(3)=0$ has a second periodic prime.  $\gamma_{s_{18}}\bmod3$ has
period $6$, pattern $[1,2,2,2,1,1]$, i.e.
$$K_{s_{18}}\equiv\frac{q+2q^2+2q^3}{1+q^3}=\frac{q+2q^2+2q^3}{(1+q)^3}\pmod3 ,$$
equivalently (using $\operatorname{sign}\gamma_{s_{18}}(n)=(-1)^{n-1}$, itself
**[verified, $n\le12000$]**)
$$|\gamma_{s_{18}}(n)|\equiv\begin{cases}1&3\nmid n\\-1&3\mid n\end{cases}\pmod 3 .$$
In particular $3\nmid\gamma_{s_{18}}(n)$ for every $n$ — the only row/prime with an empty
divisibility set.  One digit deeper: $\gamma_{s_{18}}\bmod9$ is **not** periodic, but its
restriction to $\{3\nmid n\}$ is, with period $18$ and palindromic pattern
$[1,8,\cdot,5,7,\cdot,4,2,\cdot,2,4,\cdot,7,5,\cdot,8,1,\cdot]$, whence
$$\gamma_{s_{18}}(n)\,\gamma_{s_{18}}(n+1)\equiv-1\pmod 9\qquad(n\equiv1\bmod3),$$
**[verified, $n\le12000$]**, sharp (fails mod $27$ already at $n=7$).

## 3. Nothing at any other prime  [refuted]

`23_alg.gp`, `26_deep.gp`, `26_mahler.gp`, `24_lucas.gp`.

*Rationality.*  Kernel of the Hermite–Padé matrix for $\sum_ic_iq^iK^j$, $j\le1$,
$\deg\le120$ at $p=2$ ($1500$ coefficients) and $\deg\le250$ at $3\le p\le43$ ($3000$
coefficients): **empty except** at the four cells above.  When non-empty, the candidate
relation is re-checked against all $1500$/$3000$ coefficients and passes.

*Algebraicity.*  Boxes $(\deg_q,\deg_K)=(40,6)$ for all $p\le43$ ($1500$ coefficients) and
$(20,20)$, $(60,4)$ for $p=3,5,7,11,13$ ($3000$ coefficients): empty except at the four
cells.  (The $(10,30)$ box was abandoned as too slow in interpreted GP; $(20,20)$ already
covers $\deg_K\le20$.)

*Automaticity.*  The $p$-kernel of $(\gamma(n)\bmod p)$: at $p\ge5$ all $p+1$
depth-1 subsequences $n\mapsto\gamma(pn+r)$ are pairwise distinct already on their first
$100$ terms; at $p=3$ (rows $s_7,s_{10}$) all $13$ depth-2 subsequences are distinct.  So
the $p$-kernel is growing like $p^s$ and Christol's criterion fails — independent
confirmation of the algebraicity refutation.

*Mahler/Cartier.*  $K(q)=A(q)K(q^p)+B(q)\bmod p$ with $\deg A,\deg B\le60$: solvable
exactly in the four cells, insoluble for every other $(\text{row},p)$, $p\le13$.

*Lucas/Dwork.*  All refuted, massively, over $n\le1500$ (`24_lucas.log`):

| test | typical failure rate |
|---|---|
| $\gamma(pn)\equiv\gamma(p)\gamma(n)$ | $124/136$ ($s_{18},p=11$), $93/115$ ($s_7,p=13$) |
| $\gamma(pn+r)\equiv\gamma(n)\gamma(r)$ | $1213/1373$ ($s_7,p=13$) |
| $\gamma(pn)\equiv\psi(p)\gamma(n)\pmod p$ | $500/500$ ($s_{18},p=3$), $108/136$ ($s_7,p=11$) |
| same $\bmod p^2$ | $\ge$ the above everywhere |

The only clean successes are the degenerate ones: $\gamma_{s_7}(2n)\equiv0\ (2)$,
$\gamma_{s_7}(5n)\equiv0\ (5)$, $\gamma_{s_7}(7n)\equiv0\ (7)$, $\gamma_{s_{10}}(5n)\equiv0\ (5)$
— all consequences of §1 and §5.1, not of a Lucas law.  Note in particular that round-1's
$(S)$ ($c'(pm)\equiv\psi(p)c'(m)\bmod p^2$) does **not** descend to a congruence for
$\gamma$: L3 is false in every cell.

## 4. $\gamma(p)\bmod p$: no formula  [refuted for the tested list]

`21_gammap.gp`, `28_arith.gp`, `21_extra.gp`.  $\gamma(p)=\lambda(p)$ of
`cooper_sources` §4.3 exactly (both equal $(c'(p)-\psi(p))/p^2$, since
$\beta(p)=c'(p)-\psi(p)$).  Values $\bmod p$ for $p=2,\dots,43$:

| row | $\gamma(p)\bmod p$, $p=2,3,5,7,11,13,17,19,23,29,31,37,41,43$ |
|---|---|
| $s_7$ | $0,1,0,0,9,11,2,5,10,6,26,24,6,20$ |
| $s_{10}$ | $1,2,0,3,5,9,15,5,7,14,26,14,25,41$ |
| $s_{18}$ | $1,2,2,3,8,12,2,5,16,14,16,6,15,39$ |

(the $s_{18}$ line reproduces `cooper_sources` §4.3 exactly).  Full lists to $p=599$ in
`21_gammap.log`.

* $\gamma(p)\equiv0\pmod p$ for $s_7$ at $p=2,5,7$ **and $p=167$**; for $s_{10}$ only at
  $p=5$; for $s_{18}$ never — over all $430$ primes $p\le3000$.  Heuristic expectation
  $\sum_{p\le3000}1/p=2.34$, so $4,1,0$ is Poisson-normal; $p=167$ (new, not in round 1)
  looks like a coincidence, not a fourth exceptional cell — $\gamma_{s_7}(167)$ is not
  $-167$ and no neighbouring cell behaves specially.
* Equidistribution: $\chi^2 = 1.17,\ 4.40,\ 0.42$ on four bins over $426$ primes.
* The candidate sweep of N1 (several hundred candidates: all trace forms of
  $S_{2,4,6}^{\rm new}(\Gamma_0(M))$, $M\le80$; all rational eigenforms, $M\le60$; CM eta
  quotients; Dirichlet characters; Bernoulli/Fermat/Wilson quotients; class numbers;
  harmonic sums; $p\bmod k$) produced a maximum of $5/42$ matches — for a single candidate
  the chance expectation is $\sum1/p\approx1.3$, and with several hundred candidates $5$ is
  noise.  **No formula.**
* The natural $x$-side guess $\gamma(p)\equiv(a_{p-1}-\psi(p))/p$ (the next digit of
  round-1's V9 slot $a_{p-1}\equiv\psi(p)$) scores $1/42$, $0/42$, $0/42$: **refuted**.

Honest reading: $\gamma(p)\bmod p$ carries the *second* $p$-adic digit of a unit-root
period, and there is no reason from round 1's Cartier picture for it to be modular.  The
sweep is consistent with it being a genuinely transcendental $p$-adic invariant.

## 5. Zeros, valuations, and new divisibilities

### 5.1 The zero set is exactly the multiples of $P$  [proved $\Leftarrow$, verified $\Rightarrow$]

**Proposition.**  If $\Phi|U_P=P\,\Phi$ and $\psi=\mathbf1$ then $\beta(n)=0$ for every
$n$ divisible by $P$.

*Proof.*  $\Phi|U_P=P\Phi$ gives $c(Pm)=Pc(m)$, i.e. $c'(Pm)=c'(m)$ for all $m$.  With
$\psi=\mathbf1$, $c'=\mathbf1\star\beta$, so $c'(Pm)-c'(m)=\sum_{d\mid Pm,\,d\nmid m}\beta(d)$.
Writing $m=P^vm_0$, $P\nmid m_0$, the divisors of $Pm$ not dividing $m$ are exactly
$P^{v+1}e$ with $e\mid m_0$; so $\sum_{e\mid m_0}\beta(P^{v+1}e)=0$ for every $m_0$ prime
to $P$ and every $v\ge0$.  Möbius inversion on the divisor lattice of $m_0$ gives
$\beta(P^{v+1}e)=0$ for all such $e,v$, i.e. $\beta(n)=0$ whenever $P\mid n$. $\square$

By round-1 Thm 2.2 this applies to $(s_7,7)$ and $(s_{10},5)$: $\gamma_{s_7}(7n)=0$ and
$\gamma_{s_{10}}(5n)=0$ **unconditionally, for all $n$**.  The converse
— no other zeros — is **[verified, $n\le12000$]**: `27_val.log` finds the zero sets to be
exactly $7\mathbf Z$ ($428$ zeros in $n\le3000$, $\gcd=7$, equal to $7\mathbf Z\cap[1,3000]$),
exactly $5\mathbf Z$, and $\varnothing$.  So the answer to the round-1 question is: the zero
set is **exactly** the multiples of the special prime, no larger.

### 5.2 $\{n:p\mid\gamma(n)\}$

For every other $(\text{row},p)$ this set has $\gcd=1$ and is not a union of arithmetic
progressions in any obvious way; the densities (`27_val.log`, $n\le3000$) are
$\approx1/(p-1)$-ish for $s_{18}$ but noticeably larger for $s_7,s_{10}$ at small $p$
(e.g. $s_7$, $p=3$: density $0.459$ versus the naive $1/3$) — a shadow of M9.

### 5.3 New unbounded divisibilities  [verified, $n\le12000$]

$\min\{v_p(\gamma(n)):v_p(n)=j\}$ for $j=0,\dots,5$ (`27_val.log`) is $0$ for every
$(\text{row},p)$ and every $j$ **except**:

| cell | $j=0,1,2,3,4,5$ | law (sharp) |
|---|---|---|
| $s_7$, $p=2$ | $0,1,1,1,1,1$ | $v_2(\gamma(n))=[\,2\mid n\,]$ exactly, for $7\nmid n$ |
| $s_7$, $p=5$ | $0,1,1,1,1$ | $v_5(\gamma(n))\ge1$ for $5\mid n$ |
| $s_7$, $p=3$ | $0,0,1,2,3,4$ | $\;v_3(\gamma(n))\ge v_3(n)-1$ |
| $s_{10}$, $p=2$ | $0,0,1,2,3,4$ | $\;v_2(\gamma(n))\ge v_2(n)-1$ |
| $s_{18}$, $p=2$ | $0,0,2,1,1,1$ | $v_2(\gamma(n))$ depends only on $n\bmod24$ |
| $s_{18}$, $p=3$ | $0,0,0,0,0,0$ | $3\nmid\gamma(n)$ |

The two entries $v_p(\gamma(n))\ge v_p(n)-1$ — $(s_7,3)$ and $(s_{10},2)$ — are new and
unbounded: they say $\beta(n)$ is divisible by $n^2\cdot p^{v_p(n)-1}$, i.e. **the master
divisibility $n^2\mid\beta(n)$ is not sharp at those two primes**.  Neither prime is one of
the Atkin–Lehner cells, and neither was visible in round 1.

## 6. $e_n$ and the $q$-product

`29_prod.gp`.  $e_n=\beta(n)/n=n\gamma(n)$ (the exponents of Cooper's Lambert product).
Since $e_n\bmod p$ is $n\gamma(n)\bmod p$, its periodicity is that of $\gamma$ twisted by
$n\bmod p$: periods $14,10,6$ mod $2$; $20,12$ mod $4$ for $s_{10},s_{18}$; $6$ and $18$
mod $3,9$ for $s_{18}$.  Nothing at $p\ge5$ (period $\le600$, $n\le3000$).  The mod-$9$
periodicity of $e_n$ together with the non-periodicity of $\gamma\bmod9$ is exactly the
restricted-periodicity statement of M6.

$\mathfrak g=\exp\bigl(\sum_mc'(m)q^m/m\bigr)=\prod_n(1-q^n)^{-e_n}$ is integral to
$q^{1200}$ for $s_7,s_{10}$ (round 1 had checked $q^{80}$):
$$\mathfrak g_{s_7}=1+q-3q^2+0q^3+14q^4-23q^5-17q^6+125q^7-177q^8-166q^9+\cdots$$
Its coefficients are **not** periodic mod $2,3,5$ (period $\le400$ over $1200$
coefficients), so the mod-2 rigidity of $K$ does **not** propagate to the product side.

---

## 7. What this suggests

1. The four periodic cells are exactly the cells attached to the level's 2-part and to the
   distinguished odd prime: $2^{1+v_2(N)}$ and $P$.  A proof of M1 would presumably come
   from reducing Conjecture 4.1 mod 2, where every $d^2$ collapses to $d$ and every
   $\psi(k)$ to $1$ or $0$; the content is then a statement about $\Phi\bmod2$, i.e. about
   the eta quotients themselves.  This looks genuinely attackable and would be the first
   *unconditional* information about $\gamma$ beyond the two Atkin–Lehner cells.
2. M10 identifies the mod-2 reduction of Conjecture 4.1 with a classical theta-type
   identity.  That is a concrete, checkable target: prove that $c'(m)$ is odd $\iff$ the
   $2P$-free part of $m$ is a square.  It implies M1 (invert the $\psi$-convolution mod 2),
   hence $\operatorname{rad}(n)^2\mid\beta(n)$ at $p=2$, i.e. **eq:magnetic at $p=2$ for all
   three rows**.
3. The two new divisibilities $v_3(\gamma_{s_7}(n))\ge v_3(n)-1$ and
   $v_2(\gamma_{s_{10}}(n))\ge v_2(n)-1$ say the master conjecture is not sharp at those
   primes; they should be added to the list of things Conjecture 4.1's eventual proof has
   to explain (or they are a second, independent phenomenon).
4. $\gamma(p)\bmod p$ should probably be regarded as *not* an arithmetic invariant.  The
   Shimura–Borcherds reading of `COOPER_CONGRUENCE.md` §4 predicts
   $a(|D|p^2)=p^3\gamma(p)$; N1 then says the next digit of the half-integral-weight
   coefficient at a square index is not itself modular, which is what
   Ahlgren–Ono/Edixhoven would lead one to expect.

---

## 8. Scripts and logs

| file | content |
|---|---|
| `20_extend.gp` / `22_extend3000.gp` / `20_extend6000.gp` / `20_extend12000.gp` | exact $c',\beta,\gamma$ to $n=1500,3000,6000,12000$; each re-verifies $n^2\mid\beta(n)$.  Data `20_*.txt`, `22_*.txt`, `20b_*.txt`, `20c_*.txt` |
| `21_gammap.gp` | $\gamma(p)\bmod p$ tables, $p\le599$ |
| `21_extra.gp` | $p\mid\gamma(p)$ list, multiplicativity of $\gamma\bmod2$, mod-9 restricted periodicity, $\chi^2$ |
| `21_laws6000.gp`, `21_laws12000.gp` | the thirteen closed-form laws L1–L13, re-verified at $n\le6000$, $n\le12000$ |
| `21_odd.gp` | periodicity on the odd-index subsequence (mod $4,8,16,32$) |
| `22_cprime2.gp`, `22_cprime3.gp` | the $\Xi\bmod2$ Eisenstein/theta laws for $c'(m)$ |
| `23_alg.gp` | Hermite–Padé rationality ($\deg\le120$) and algebraicity $(40,6)$ for $K,\Xi,B$, $p\le43$ |
| `24_lucas.gp` | the Lucas/Dwork/periodicity battery (a)–(d) |
| `25_ratform.gp` | minimal $P/Q$ for $K\bmod p$; periods of $\gamma\bmod p^k$ |
| `25_patterns.gp` | the explicit periodic patterns |
| `25_mod9.gp` | the mod-9 pair law for $s_{18}$ |
| `26_dbg.gp`, `26_dbg2.gp` | independent verification that $QK-P\equiv0$ for all four rational forms, to $q^{12000}$ |
| `26_deep.gp` | rationality $\deg\le250$ and algebraic boxes $(20,20)$, $(60,4)$ |
| `26_mahler.gp` | test (e) $K(q)=A(q)K(q^p)+B(q)$, and the $p$-kernel |
| `27_val.gp` | zero sets, $\{n:p\mid\gamma(n)\}$, $\min v_p(\gamma)$ tables, periodicity of $v_p$ |
| `28_arith.gp` | the $\gamma(p)\bmod p$ candidate sweep |
| `29_prod.gp` | $e_n$, the $q$-product $\mathfrak g$, and the first pass at L1–L11 |

*Caveat recorded in the logs:* the final block of `21_laws6000.log` / `21_laws12000.log`
("restricted to $n$ odd") is void — with an odd shift the condition "$n$ and $n+t$ both
odd" is vacuous, so the loop trivially returns period 1.  It is superseded by
`21_odd.gp` / `21_odd.log`, which tests the odd-index subsequence correctly.
