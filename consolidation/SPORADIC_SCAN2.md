# Scan 2: Apéry-like rows from the full weight-$w$ form space on genus-zero parameters

*Claude (Opus 5), 2026-08-21/22. Scripts and raw data: `lattice/sporadic_scan2/`.
Compute: the box "snake" (8 cores, PARI 2.17.4).*

**STATUS: run in progress — this file is a partial report, updated as results land.**

## 0. What this scan does differently

`SPORADIC_SEARCH.md` (scan 1) enumerated $520\,005$ pairs $(t,F)$ with **both** $t$ and $F$
eta quotients on the $14$ genus-zero levels $N\le 25$. Its own §7 lists the two gaps that
matter, and Problem 7.4 of the paper restates them:

1. **$F$ need not be an eta quotient.** Apéry's and $T$'s weight-2 forms are Eisenstein
   *combinations* (`ZETA3_TWO_LATTICE.md` §9); Zagier's row D needs a weight-1 form with a
   *quartic* nebentypus mod 5. Scan 1 could not see either.
2. **The group need not be $\Gamma_0(N)$ with $N\le25$.**

Scan 2 fixes both, and replaces the "guess $F$, then test" loop by an exact solve.

## 1. Method

### 1.1 The parameter $t$: enumerate the *divisor*, not the exponent vector

$t=\prod_{d\mid N}\eta(d\tau)^{r_d}$ of weight $0$ is determined by its cusp-order vector
$\bigl(\operatorname{ord}_c t\bigr)_{c\mid N}$ through the (invertible) Ligozat matrix
$A_{c,d}=\frac{N}{24\gcd(c^2,N)}\frac{\gcd(c,d)^2}{d}$. So instead of boxing $|r_d|\le24$
(scan 1: $10^5$–$10^6$ vectors per level, almost all of huge degree), scan 2 enumerates the
*divisors*
$$\operatorname{div}(t)=(\infty)+\Delta_+-\Delta_-,\qquad
\deg t=\sum_c\varphi(\gcd(c,N/c))\max(\operatorname{ord}_c t,0)\le 4,$$
and keeps those with $r=A^{-1}(\text{ord})\in\mathbf Z^{\tau(N)}$ (a congruence condition,
solved by meet-in-the-middle on residues). This is **exhaustive** for $\deg t\le4$ and takes
seconds for all $N\le120$.

Why degree $\le4$ is the right box: $\deg t=1$ means $t$ is a Hauptmodul of $\Gamma_0(N)$;
$\deg t=2$ means $t$ is the pull-back of a Hauptmodul of an Atkin–Lehner quotient
$\Gamma_0(N)+W_Q$; $\deg t=4$ covers quotients by a Klein four-group of Atkin–Lehner
involutions. **Both Apéry's and Domb's parameters have $\deg t=2$**, which is exactly why
scan 1's earlier "must be a Hauptmodul of $\Gamma_0(N)$" pass lost them.

Result: **3 979 distinct parameters** over **51 levels** $N\le104$
($N=2,\dots,10,12,\dots,16,18,20,21,22,24,\dots,28,30,32,33,35,36,39,40,42,44,45,46,48,50,52,54,56,60,63,64,66,70,72,75,78,80,81,88,92,96,100,104$),
of which **227 have $\deg t\le2$**. Apéry $(\eta_1\eta_6/\eta_2\eta_3)^{12}$, Domb
$(\eta_2\eta_6/\eta_1\eta_3)^6$ and $T=(\eta_1\eta_8/\eta_2\eta_4)^8$ are all in the list with
$\deg t = 2$.

### 1.2 The form $F$: the full integral lattice, one nebentypus Galois orbit at a time

For each level $M$ and weight $w$ we build, in PARI, the $\mathbf Z$-lattice
$$\Lambda_{M,w,\chi}=\Bigl\{\,f\in\bigoplus_{\sigma}M_w(\Gamma_0(M),\chi^\sigma)\ \Bigm|\
f \text{ has a rational, integral }q\text{-expansion}\,\Bigr\}$$
from the Galois traces $\operatorname{Tr}_{K/\mathbf Q}(\zeta^j f_i)$ of an `mfbasis`, saturated
with `matkerint`. This contains Eisenstein combinations, cusp forms and eta quotients alike —
it is the whole space, which is the point of Problem 7.4. Orbits are kept only when
$\varphi(\operatorname{ord}\chi)\cdot(w\psi(M)/12+2)+10\le 210$, i.e. when the $q$-expansion
precision provably determines the space.

### 1.3 The search over $F$ is a *solve*, not a box

$a_n(F)$ is **linear** in $F$, so for $F=F_0+\sum_j m_jg_j$ the condition
$$\sum_{i\le r}\ \sum_{d\le D} c_{i,d}\,n^d\,a_{n+i}(F)=0$$
is **bilinear** in $(c,m)$. Setting $u_{i,d,j}=c_{i,d}m_j$ ($m_0:=1$) linearises it; scan 2
computes the nullspace of the resulting $\bigl((r{+}1)(D{+}1)(K{+}1)\bigr)$-column system over
$\mathbf F_p$ and extracts its **rank-one points** by a matrix-pencil computation (roots of
$\det(A_1-\lambda A_0)$ by Cantor–Zassenhaus), then reconstructs $m$ rationally and keeps only
$m\in\mathbf Z^K$.

This matters: Apéry's $F$ sits at $m=(5,8)$ in the LLL-reduced basis of
$\Lambda_{6,2,\mathbf 1}$, so **no small box finds it** — a box of radius 2, which is all one
can afford in a scan of this size, misses Apéry, Domb, $T$ and all six Zagier rows. The
bilinear solve finds them with no box at all.

### 1.4 Invariants

For each surviving $(t,M,\chi,w,F)$: the minimal $(r,D)$ recurrence is refitted **exactly over
$\mathbf Q$** and verified for all $n\le200$; $\lambda_i$ are the roots of the leading symbol
$\sum_i c_{i,D}x^i$; $c=\lambda_1\cdots$; the companion is the second solution
$b_0=\cdots=b_{r-2}=0$, $b_{r-1}=1$ of the same recurrence (the census convention), extended
exactly to $n=700$; $k$ is the least exponent with $d_n^kb_n\in\mathbf Z$ for $n\le200$;
$\sigma_p$ is measured from $v_p(a_nb_{n+1}-a_{n+1}b_n)$ at $n=100,300$; the Apéry limit is
$b_n/a_n$ at $n\approx700$, to as many digits as $\log_{10}(\lambda_1/\lambda_2)$ allows
(often several hundred), and is identified by `lindep` against the constant basis **and**
against $L(f,s)$, $s=w,w{+}1,w{+}2$, for every newform $f$ of weight $w+2$ at every level
dividing $M$ or equal to $2M,3M,4M$.

The Eichler-integral companion $B=F\cdot D^{-(w+1)}(F\,Dt)$ of the brief is computed as well
and its own $k$ recorded; for $w=2$ rows it reproduces the census limits exactly
($\zeta(3)/6$, $7\zeta(3)/24$, $7\zeta(3)/32$, …), but for many $w=1$ rows it is **not** a
solution of the fitted recurrence, so the recurrence companion is used as primary.

## 2. Coverage (honest)

**Covered.** 51 levels $N\le104$ carry an eta-quotient parameter of cusp-divisor degree
$\le4$; all $15$ genus-zero $\Gamma_0(N)$ are among them, together with the Atkin–Lehner
quotients they reach by degree-$2$ and degree-$4$ pull-backs
(deg 1: $N=2$–$10,12,13,16,18,20,24,25,32$;
 deg 2: $4,6,8,9,10,12,14,15,16,18,20,21,22,24,26,28,30,32,33,35,36,39,40,46,48,50$;
 deg 3: $4,6,8,10,12,14,15,16,18,20,21,24,27,28,30,32,36,40,45,48,52,64,81$;
 deg 4: 39 levels up to $104$).
Weights $w=1,2,3$; nebentypus: every Galois orbit of Dirichlet characters mod $M$ whose
$\mathbf Q$-dimension fits inside the $q$-precision (so at small level *all* characters,
including the quartic ones mod $5$ that Zagier's row D needs; at large level, orbits with
$\varphi(\operatorname{ord}\chi)\le4$).

**Not covered — the honest list.**

1. **Levels with no low-degree eta-quotient parameter.** At prime level the only weight-$0$
   eta quotients are $(\eta_1/\eta_p)^k$ with $24\mid k(p-1)$, so $p\in\{2,3,5,7,13\}$ only.
   Consequently the genus-zero Fricke groups
   $$N=11,17,19,23,29,31,34,38,41,47,49,51,55,59,62,69,71,84,87,94,95,105,110,119$$
   are **absent** from this scan: their Hauptmoduls are McKay–Thompson series, not eta
   quotients, and scan 2 does not construct them. This is the single largest gap and it is
   inherited from scan 1.
2. **Cusp-divisor degree $>4$**, i.e. quotients by Atkin–Lehner groups of order $\ge8$
   (relevant from $N$ with three prime factors upward).
3. **High-dimensional form spaces.** The bilinear solve uses the first $K$ vectors of an
   LLL-reduced basis of $\{F\in\Lambda_{M,w,\chi}: F(0)=0\}$, with $K+1\le(205-r)/((r{+}1)(D{+}1))$
   (so $K\le20$ for $w=1$, $\le15$ for $w=2$, $\le12$ for $w=3$ at the shortest shape). For
   $M\le36$ this is the whole lattice; above that it is a genuine restriction and the scan is
   **not** exhaustive there.
4. **$F$-levels** are restricted to multiples $M$ of the parameter level $N$ with $M\le4N$
   and $M$ in the covered list.
5. **Recurrence shapes**: length $r+1\le5$ and polynomial degree $D\le w+2$. A row whose
   Picard–Fuchs operator has more than $4$ finite singular points is invisible.
6. Rows over parameters $t$ whose divisor is **not supported on cusps** (poles at elliptic
   points) are outside the eta-quotient family entirely.

## 3. Acceptance test: the known census is recovered

Every row of `paper/tables/census.tex` that lies in the searched family is reproduced, with
its recurrence, its $\lambda_i$, its sharp $k$ and its period recomputed from scratch by
`lindep` (no ledger input):

| row | found at | $a_n$ | recomputed period |
|---|---|---|---|
| Apéry $(17,5,1)$ | $N=6$, $t=(\eta_1\eta_6/\eta_2\eta_3)^{12}$, $w=2$, $M=6$, $m=(5,8)$ | $1,5,73,1445,33001$ | $\zeta(3)/6$ |
| $T=(12,4,16)$ | $N=8$, $t=(\eta_1\eta_8/\eta_2\eta_4)^8$, $w=2$ | $1,4,40,544,8536$ | $7\zeta(3)/32$ |
| Domb $(10,4,64)$ | $N=6$, $t=(\eta_2\eta_6/\eta_1\eta_3)^6$, $w=2$ | $1,-4,28,-256,2716$ | $7\zeta(3)/24$ |
| AZ $(9,3,-27)$ | $N=9$, $w=2$ | $1,3,27,309,4059$ | $L(3,\chi_{-3})/3$ |
| AZ $(11,5,125)$ | $N=5$, $w=2$ | $1,-5,35,-275,2275$ | complex roots, no arch. limit |
| AZ $(7,3,81)$ | $N=6$, $M=12$, $w=2$ | $1,-3,9,-3,-279$ | complex roots |
| Zagier A $(7,2,-8)$ | $N=6$, $w=1$ | $1,2,10,56,346$ | $\zeta(2)/4$ |
| Zagier B $(9,3,27)$ | $N=9$, $w=1$ | $1,-3,9,-21,9$ | complex roots |
| Zagier C $(10,3,9)$ | $N=6$, $w=1$ | $1,3,15,93,639$ | $L(2,\chi_{-3})/2$ |
| Zagier E $(12,4,32)$ | $N=8$, $w=1$ | $1,4,20,112,676$ | $G/2$ |
| Zagier F $(17,6,72)$ | $N=6$, $w=1$ | $1,-6,42,-312,2394$ | $5L(2,\chi_{-3})/8$ |
| $L(f,2)$ row (scan 1 §5) | $N=6$ Domb $t$, $w=1$, $M=12$, $\chi_{-3}$ | $1,-2,12,-104,1078$ | $L(f,2)$, $f=\eta_2^3\eta_6^3$ |

The last line is the point of the exercise: scan 1 found that row only because its $F$ happens
to be an eta quotient *of level 12*, and had to be told to look there. Scan 2 finds it as one
element of the integral lattice $\Lambda_{12,1,\chi_{-3}}$ with no shape hypothesis, and it
also produces the object of §4 which scan 1 could not see at all.

## 4. The main new object: the $\operatorname{Sym}^1$ row on the $T$ curve

$$\boxed{\ (n+1)^2a_{n+1}=(24n^2+12n+2)\,a_n-4(2n-1)^2\,a_{n-1},\qquad a_0=1,\ a_1=2\ }$$
$$a_n=1,\;2,\;18,\;236,\;3634,\;61236,\;1093356,\;20325096,\;389236914,\;7626179396,\dots$$

**Modular source.** The parameter is the one of the $\zeta(3)$ row $T=(12,4,16)$ of
`ZETA3_TWO_LATTICE.md` §9,
$$t=\Bigl(\frac{\eta(\tau)\eta(8\tau)}{\eta(2\tau)\eta(4\tau)}\Bigr)^{8}=q-8q^2+\cdots
\quad(\Gamma_0(8),\ \deg t = 2),$$
and $F$ is the unique integral weight-one form with $F(0)=1$ in
$M_1(\Gamma_0(8),\chi_{-8})$ (which is one-dimensional):
$$F=\sum_{x,y\in\mathbf Z}q^{x^2+2y^2}=1+2q+2q^2+4q^3+2q^4+4q^6+2q^8+6q^9+\cdots$$
— the CM theta series of $\mathbf Z[\sqrt{-2}]$. So $F^2$ is $T$'s weight-two form up to the
$\operatorname{Sym}$-degree: this row is to $T$ exactly what scan 1's $L(f,2)$ row is to Domb.

**Invariants** (all recomputed; $a_n\in\mathbf Z$ verified to $n=900$, recurrence to $n=200$):

| | value |
|---|---|
| char. polynomial | $x^2-24x+16$ |
| $\lambda_1$ | $12+8\sqrt2=23.3137\ldots$ |
| $\lambda_2$ | $12-8\sqrt2=0.68629\ldots\ \mathbf{<1}$ |
| $c=\lambda_1\lambda_2$ | $16$ |
| $k$ (sharp) | $\mathbf 2$ ($d_n^2b_n\in\mathbf Z$ for $n\le200$; $d_n^1b_n\notin\mathbf Z$) |
| $\sigma_2$ (measured) | $0$ (the trailing coefficient is $4(2n-1)^2$, not $cn^2$; the Casoratian is $4^n((2n-1)!!)^2/((n+1)!)^2$, so $v_2=O(\log n)$) |
| score $=\log(1/\lambda_2)-k$ | $-1.624$ |
| budget $=\log\lambda_1-k$ | $\mathbf{+1.149}$ |

**Period.** With the census companion $b_0=0,b_1=1$,
$$\lim_{n\to\infty}\frac{b_n}{a_n}
=0.52874945559576683377595806385986587134115392973305253672438329\ldots
=L(f,2),$$
$$f=\eta(\tau)^2\eta(2\tau)\eta(4\tau)\eta(8\tau)^2=q-2q^2-2q^3+4q^4+4q^6-8q^8-5q^9+14q^{11}-\cdots
\in S_3(\Gamma_0(8),\chi_{-8}),$$
the weight-three CM newform attached to $\mathbf Q(\sqrt{-2})$ — **exactly on the nose**, no
rational factor, agreeing to the $140$ digits at which `lindep` was run (the limit itself is
available to $>1000$ digits).

**Why it is the most interesting row in the scan.**

* Its **budget $+1.149$ is the largest of any row here with a genuine second solution and a
  non-rational period** — larger than the scan-1 $L(f,2)$ row ($+0.773$) and than Apéry's own
  $(+0.526)$.
* Unlike every other cusp-form-period row known in this project it **decays**:
  $\lambda_2=12-8\sqrt2<1$, so $b_n/a_n\to L(f,2)$ geometrically and the *unnormalised* linear
  form $a_nL(f,2)-b_n$ shrinks like $(12-8\sqrt2)^n$. The scan-1 $L(f,2)$ row on the Domb curve
  has $\lambda_2=4>1$ and grows.
* Its score is nevertheless $\log(1/(12-8\sqrt2))-2=0.3763-2=-1.624<0$: $d_n^2$ beats the
  decay, so **nothing is proved and nothing is claimed**. What the row supplies is the
  *decaying* half of a two-lattice pair, for the period $L(f,2)$ of the level-$8$ weight-$3$ CM
  form. It has no $p$-adic slope of its own ($\sigma_2=0$), so as an engine it is worthless;
  as a decayer it wants a partner with the same period and $\sigma_2>0$.

**The structural statement.** Write $\mathcal C$ for a genus-zero parameter carrying an
integral $\operatorname{Sym}^2$ row. Then:

| curve | $\operatorname{Sym}^2$ row | $\operatorname{Sym}^1$ row on the same $t$ | its period | $\lambda_1,\lambda_2$ | $k$ | budget |
|---|---|---|---|---|---|---|
| Domb, $\Gamma_0(6)$ | $(10,4,64)$, $7\zeta(3)/24$ | $F=\eta_1^2\eta_3^2/\eta_2\eta_6$ (CM $\theta$, disc $-3$) | $L(\eta_2^3\eta_6^3,2)$ | $16,\;4$ | 2 | $+0.773$ |
| $T$, $\Gamma_0(8)$ | $(12,4,16)$, $7\zeta(3)/32$ | $F=\theta_{x^2+2y^2}$ (CM $\theta$, disc $-8$) | $L(\eta_1^2\eta_2\eta_4\eta_8^2,2)$ | $12{+}8\sqrt2,\;12{-}8\sqrt2$ | 2 | $\mathbf{+1.149}$ |
| Apéry, $\Gamma_0(6)$ | $(17,5,1)$, $\zeta(3)/6$ | **none** | — | — | — | — |

The third line is a clean negative found by the scan: over Apéry's parameter
$(\eta_1\eta_6/\eta_2\eta_3)^{12}$ the bilinear solve returns **no** weight-one row at any
$F$-level $M\in\{6,12,18,24\}$ and any nebentypus, for any recurrence shape with
$r\le4,\ D\le3$. (Had one existed it would have had $\lambda_1=33.97$, $\lambda_2=1/33.97$ and
$k=2$, hence score $+1.53$; it does not exist.) The $\operatorname{Sym}^1$ descent is available
exactly when $c=\lambda_1\lambda_2$ is a square in the relevant sense — Domb $64$, $T$ $16$ —
and Apéry's $c=1$ row is the Fricke-perfect one, whose weight-one form would have to be
$\sqrt{F}$ with $F=1+5q+13q^2+\cdots$, not a $q$-series over $\mathbf Z$.

## 5. What this does to Problems 7.3 and 7.4

**Problem 7.4 (new sporadic sequences).** The construction asked for is now implemented in the
form the problem specifies — growth filter / full form space, no eta-quotient ansatz for $F$ —
and its output is §6's table. The headline is negative in the same direction as scan 1: **no
new $\operatorname{Sym}^2$ or $\operatorname{Sym}^3$ row in Zagier's or Almkvist–Zudilin's
normalisation appears at any level in range.** Everything new is weight one.

**Problem 7.3 (decaying partners).** The problem asks, among other things, for a decaying row
for $L(f,2)$, $f=\eta(2\tau)^3\eta(6\tau)^3$, at $p=2$ with $\Lambda_{\rm dec}>14.8$. Scan 2
does **not** find one: the only rows in range with that exact period are the scan-1 row and its
reparametrisations, all with $\lambda_2=4$. What it does find is the *analogous configuration
one level over*: the level-$8$ CM form $f_8=\eta_1^2\eta_2\eta_4\eta_8^2$ has a **decaying** row
($\lambda_2=12-8\sqrt2$, $k=2$). The missing half there is the opposite one: a row with period
$L(f_8,2)$ and $\sigma_2>0$ to act as engine. The scan produced no such row either, so the
$(\text{engine},\text{decayer})$ pair is still incomplete — but the two halves now exist in the
census for *different* CM periods, which is a sharper statement of the obstruction than
"no decayer is known".

## 6. The ranked table

(see `lattice/sporadic_scan2/table.md` for the full machine-generated list, and
`table.json` for the raw record of every row including its exact recurrence)
