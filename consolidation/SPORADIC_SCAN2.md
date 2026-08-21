# Scan 2: Apéry-like rows from the full weight-$w$ form space on genus-zero parameters

*Claude (Opus 5), 2026-08-21/22. Scripts and raw data: `lattice/sporadic_scan2/`.
Compute: the box "snake" (8 cores, PARI 2.17.4).*

**Headline.** The scan ($5.5\times10^5$ exact bilinear solves; 3 979 parameters over
51 levels $N\le104$; weights $1,2,3$; the *whole* integral form space, not an eta-quotient
ansatz) reproduces every census row in range and produces **no new $\operatorname{Sym}^2$ or
$\operatorname{Sym}^3$ row**. What it does produce is a family of weight-one rows over
Atkin–Lehner quotients, of which one is genuinely new and good: over the $\zeta(3)$ parameter
$T$ of `ZETA3_TWO_LATTICE.md` §9, the weight-one CM theta series of $\mathbf Z[\sqrt{-2}]$
gives $a_n=1,2,18,236,3634,\dots$ with $\lambda_1=12+8\sqrt2$, **$\lambda_2=12-8\sqrt2<1$**,
$k=2$, budget $+1.149$ and period exactly $L(f,2)$ for the weight-three level-$8$ CM newform
$f=\eta_1^2\eta_2\eta_4\eta_8^2$ (§4).

Following that object upstream gave the main result, **Part II**: *every* sporadic order-three
row (six Almkvist–Zudilin + three Cooper) has an integral $\operatorname{Sym}^1$ **square
root** $a_n=\lambda^n[t^n]\sqrt{\sum A_nt^n}$, satisfying a three-term recurrence with
quadratic coefficients and costing $k=2$ instead of $k=3$ — a free integration in the sense of
Problem 7.6. Apéry's own square root, $a_n=1,10,534,40900,3672550,\dots$, has
$\lambda_1=4(17{+}12\sqrt2)$, $\lambda_2=4(17{-}12\sqrt2)$, $k=2$ and **score $+0.1392>0$**;
its period is *unidentified* and no claim is made about it (§7).

A clean dichotomy falls out. Three of the nine $\operatorname{Sym}^1$ rows have a weight-one
form that is an integral CM theta series of class number one, and their periods are
$L(f,2)$ for the corresponding weight-three CM newform, exactly and with no rational factor:
disc $-7$ (Cooper $s_7$) $\to L(\eta_1^3\eta_7^3,2)$, disc $-8$ ($T$) $\to
L(\eta_1^2\eta_2\eta_4\eta_8^2,2)$, disc $-3$ (Domb) $\to L(\eta_2^3\eta_6^3,2)$. The other
four (Apéry, AZ$(9,3,-27)$, Cooper $s_{10}$, $s_{18}$) resist a $1792$-element `lindep`.

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

*Caveat on "degree".* $\deg t$ above is the divisor degree computed on $X_0(N)$. When
$\prod_d d^{\,r_d}$ is not a rational square the eta quotient carries a quadratic character
$\chi(d)$ and is a function not on $X_0(N)$ but on the index-two subgroup $\ker\chi$; the
number is then a bound on a related quantity rather than a map degree. Nebentypus is
deliberately not excluded (scan 1 records that excluding it kills level $5$ and hence
$(11,5,125)$), so the enumeration is a well-defined finite family that *contains* the
Atkin–Lehner Hauptmoduls and some $\Gamma_0(N)\cap\Gamma_1$-type ones besides; the
recurrence fit, not the degree label, is the actual filter.

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
| Apéry, $\Gamma_0(6)$ | $(17,5,1)$, $\zeta(3)/6$ | $F=\sqrt{F_{\text{Apéry}}}$, integral only after $t\mapsto t/4$ (§7) | unidentified | $4(17{\pm}12\sqrt2)$ | 2 | $+2.912$ |

The third line needs a correction that Part II supplies. Over Apéry's parameter
$(\eta_1\eta_6/\eta_2\eta_3)^{12}$ the bilinear solve returns **no** weight-one row at any
$F$-level $M\in\{6,12,18,24\}$, any nebentypus, and any shape with $r\le4$, $D\le3$ — but
that is an artefact of the scan's normalisation, not a fact about the curve. The
$\operatorname{Sym}^1$ row of Apéry's curve *does* exist; it is integral only after the
substitution $t\mapsto t/4$, so its parameter is $q/4+O(q^2)$ and the scan, which insists on
$t=q+O(q^2)\in\mathbf Z[[q]]$, cannot see it. It is the row of §7, and it is the best object
in this report. The pattern is: $F_{\text{Domb}}$ and $F_T$ are literally squares of
integral weight-one forms (CM theta series of discriminant $-3$ and $-8$), so $\lambda=1$;
$F_{\text{Apéry}}=1+5q+13q^2+\cdots$ is not, its square root has half-integral coefficients,
and the denominators are cleared exactly by $4^n$.

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

85 distinct integral rows survive de-duplication (a row is identified with its
$a_n$ up to $a_n\mapsto(-1)^na_n$; the same sequence typically arises from several
$(t,M,\chi,F)$, and the earliest is listed). Below are the top 26 **by budget**, restricted
to rows with $k\ge1$ — a row with $k=0$ has an integral companion, hence a rational Apéry
limit and a trivialised extension class, and is not interesting however large its
$\lambda_1$ (there are 12 such, at the very top of the raw ranking).
Full list: `lattice/sporadic_scan2/table.md` (all 85 rows, with $\lambda_{\min}$ and the
exact recurrence in `table.json`).

Notation: rec $=(\text{number of terms},\deg p_i)$; $M,\chi,w$ = level, Conrey label and
weight of the form space $F$ was taken from; "--" in the score column means
complex-conjugate leading roots (no archimedean limit).

| # | $N$ | $M,\chi,w$ | $t$ | $\deg t$ | rec | $\lambda_1$ | $\lambda_2$ | $c$ | $k$ | $\sigma_p$ | score | budget | $a_n$ | period |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 2 | $2,1,2$ | $\eta_{2}^{24}/\eta_{1}^{24}$ | 1 | $(4,3)$ | 64.0003 | 64.0003 | -262144 | 2 | $\sigma_{2}{=}6$ | -- | **2.159** | $1,24,-552,19392,-810024$ | -- |
| 8 | 3 | $3,1,2$ | $\eta_{3}^{12}/\eta_{1}^{12}$ | 1 | $(4,3)$ | 27.0001 | 26.9999 | -19683 | 2 | $\sigma_{3}{=}3$ | -5.296 | **1.296** | $1,12,-108,1524,-25884$ | -- |
| 9 | 8 | $8,3,1$ | $\eta_{1}^{8}\eta_{8}^{8}/\eta_{2}^{8}\eta_{4}^{8}$ | 2 | $(3,2)$ | 23.3137 | 0.6863 | 16 | 2 | -- | -1.624 | **1.149** | $1,2,18,236,3634$ | $L(f_{8,3,\chi_{3}},2)$ |
| 10 | 6 | $6,5,1$ | $\eta_{2}^{6}\eta_{6}^{6}/\eta_{1}^{6}\eta_{3}^{6}$ | 2 | $(4,2)$ | 16.0000 | 16.0000 | -1024 | 2 | $\sigma_{2}{=}4$ | -- | **0.773** | $1,6,-36,312,-3210$ | -- |
| 11 | 4 | $16,15,1$ | $\eta_{4}^{8}/\eta_{1}^{8}$ | 1 | $(3,2)$ | 16.0000 | 16.0000 | 256 | 2 | $\sigma_{2}{=}2$ | -4.773 | **0.773** | $1,0,-4,64,-924$ | -- |
| 12 | 6 | $12,5,1$ | $\eta_{2}^{6}\eta_{6}^{6}/\eta_{1}^{6}\eta_{3}^{6}$ | 2 | $(3,2)$ | 16.0000 | 4.0000 | 64 | 2 | $\sigma_{2}{=}2$ | -3.386 | **0.773** | $1,-2,12,-104,1078$ | **scan-1 $L(f,2)$ row**, $-L(f_{12,3,\chi_{5}},2)$ |
| 13 | 12 | $12,5,1$ | $\eta_{1}^{5}\eta_{3}^{1}\eta_{4}^{1}\eta_{12}^{5}/\eta_{2}^{6}\eta_{6}^{6}$ | 2 | $(4,2)$ | 14.9282 | 1.0718 | 16 | 2 | $\sigma_{2}{=}1$ | -2.069 | **0.703** | $1,2,10,80,778$ | unidentified |
| 14 | 12 | $12,5,1$ | $\eta_{2}^{12}\eta_{3}^{4}\eta_{12}^{8}/\eta_{1}^{4}\eta_{4}^{8}\eta_{6}^{12}$ | 2 | $(4,2)$ | 13.9282 | 1.0000 | 1 | 2 | -- | -2.000 | **0.634** | $1,-1,7,-55,511$ | unidentified |
| 15 | 12 | $12,5,1$ | $\eta_{2}^{12}\eta_{3}^{4}\eta_{12}^{8}/\eta_{1}^{4}\eta_{4}^{8}\eta_{6}^{12}$ | 2 | $(5,2)$ | 13.9282 | 1.0000 | 1 | 2 | -- | -2.000 | **0.634** | $1,-2,8,-62,566$ | unidentified |
| 16 | 6 | $6,1,2$ | $\eta_{1}^{12}\eta_{6}^{12}/\eta_{2}^{12}\eta_{3}^{12}$ | 2 | $(3,3)$ | 33.9706 | 0.0294 | 1 | 3 | -- | 0.525 | **0.525** | $1,5,73,1445,33001$ | **Apéry**, $\tfrac{1}{6}\zeta(3)$ |
| 17 | 16 | $16,15,1$ | $\eta_{1}^{4}\eta_{4}^{4}\eta_{16}^{4}/\eta_{2}^{6}\eta_{8}^{6}$ | 2 | $(4,2)$ | 11.6569 | 2.0000 | 8 | 2 | $\sigma_{2}{=}1$ | -2.693 | 0.456 | $1,2,8,48,356$ | unidentified |
| 18 | 16 | $16,15,1$ | $\eta_{1}^{4}\eta_{4}^{4}\eta_{16}^{4}/\eta_{2}^{6}\eta_{8}^{6}$ | 2 | $(5,2)$ | 11.6569 | 2.0000 | 16 | 2 | $\sigma_{2}{=}1$ | -2.693 | 0.456 | $1,0,4,32,260$ | unidentified |
| 19 | 15 | $15,14,1$ | $\eta_{1}^{3}\eta_{15}^{3}/\eta_{3}^{3}\eta_{5}^{3}$ | 2 | $(5,2)$ | 11.0902 | 1.6180 | 1 | 2 | -- | -2.481 | 0.406 | $1,1,5,31,229$ | unidentified |
| 20 | 10 | $40,39,1$ | $\eta_{2}^{4}\eta_{10}^{4}/\eta_{1}^{4}\eta_{5}^{4}$ | 2 | $(4,2)$ | 10.4721 | 4.0000 | -64 | 2 | -- | -3.386 | 0.349 | $1,-2,8,-44,290$ | $L(f_{20,3,\chi_{19}},2) - L(f_{20,3,\chi_{19}}^{[2.1]},2)$ |
| 21 | 10 | $40,39,1$ | $\eta_{2}^{4}\eta_{10}^{4}/\eta_{1}^{4}\eta_{5}^{4}$ | 2 | $(5,2)$ | 10.4721 | 4.0000 | 256 | 2 | -- | -3.386 | 0.349 | $1,0,2,-20,168$ | unidentified |
| 22 | 3 | $9,1,2$ | $\eta_{3}^{12}/\eta_{1}^{12}$ | 1 | $(4,3)$ | 27.0001 | 26.9999 | -19683 | 3 | $\sigma_{3}{=}2$ | -6.296 | 0.296 | $1,3,-54,1038,-21375$ | -- |
| 23 | 3 | $9,1,2$ | $\eta_{3}^{12}/\eta_{1}^{12}$ | 1 | $(3,3)$ | 27.0000 | 27.0000 | 729 | 3 | $\sigma_{3}{=}2$ | -6.296 | 0.296 | $1,-6,81,-1392,26874$ | -- |
| 24 | 16 | $16,15,1$ | $\eta_{1}^{2}\eta_{2}^{3}\eta_{8}^{3}\eta_{16}^{2}/\eta_{4}^{10}$ | 2 | $(5,2)$ | 9.6569 | 2.0000 | -64 | 2 | $\sigma_{2}{=}1$ | -2.693 | 0.268 | $1,2,4,24,148$ | unidentified |
| 25 | 16 | $16,15,1$ | $\eta_{2}^{9}\eta_{8}^{3}\eta_{16}^{2}/\eta_{1}^{2}\eta_{4}^{12}$ | 2 | $(4,2)$ | 9.6569 | 2.0000 | -32 | 2 | $\sigma_{2}{=}2$ | -2.693 | 0.268 | $1,0,4,-16,116$ | unidentified |
| 26 | 6 | $6,5,1$ | $\eta_{2}^{1}\eta_{6}^{5}/\eta_{1}^{5}\eta_{3}^{1}$ | 1 | $(4,2)$ | 9.0000 | 9.0000 | -648 | 2 | $\sigma_{3}{=}2$ | -4.197 | 0.197 | $1,3,-12,66,-414$ | -- |
| 27 | 6 | $24,7,1$ | $\eta_{1}^{4}\eta_{6}^{8}/\eta_{2}^{8}\eta_{3}^{4}$ | 1 | $(5,2)$ | 9.0000 | 9.0000 | 81 | 2 | -- | -4.197 | 0.197 | $1,-4,-12,-56,-316$ | -- |
| 28 | 6 | $6,5,1$ | $\eta_{1}^{4}\eta_{6}^{8}/\eta_{2}^{8}\eta_{3}^{4}$ | 1 | $(4,2)$ | 9.0000 | 9.0000 | 81 | 2 | $\sigma_{3}{=}2$ | -4.197 | 0.197 | $1,-6,-12,-42,-198$ | -- |
| 29 | 6 | $6,5,1$ | $\eta_{2}^{1}\eta_{6}^{5}/\eta_{1}^{5}\eta_{3}^{1}$ | 1 | $(4,2)$ | 9.0000 | 8.0000 | -576 | 2 | $\sigma_{2}{=}3$ | -4.079 | 0.197 | $1,2,-6,24,-102$ | $1/2 + \tfrac{3}{4}\zeta(2)$ |
| 30 | 6 | $24,7,1$ | $\eta_{2}^{1}\eta_{6}^{5}/\eta_{1}^{5}\eta_{3}^{1}$ | 1 | $(4,2)$ | 9.0000 | 8.0000 | -576 | 2 | -- | -4.079 | 0.197 | $1,0,0,-4,60$ | unidentified |
| 31 | 6 | $24,7,1$ | $\eta_{2}^{1}\eta_{6}^{5}/\eta_{1}^{5}\eta_{3}^{1}$ | 1 | $(4,2)$ | 9.0000 | 8.0000 | -576 | 2 | -- | -4.079 | 0.197 | $1,-4,24,-164,1196$ | unidentified |
| 32 | 9 | $9,8,1$ | $\eta_{9}^{3}/\eta_{1}^{3}$ | 1 | $(5,2)$ | 9.0000 | 9.0000 | 2187 | 2 | $\sigma_{3}{=}2$ | -- | 0.197 | $1,6,-18,60,-180$ | -- |
| 33 | 6 | $6,5,1$ | $\eta_{1}^{4}\eta_{6}^{8}/\eta_{2}^{8}\eta_{3}^{4}$ | 1 | $(5,2)$ | 9.0000 | 3.0000 | 81 | 2 | $\sigma_{3}{=}1$ | -3.099 | 0.197 | $1,0,6,48,360$ | unidentified |
| 34 | 6 | $6,5,1$ | $\eta_{1}^{4}\eta_{6}^{8}/\eta_{2}^{8}\eta_{3}^{4}$ | 1 | $(4,2)$ | 9.0000 | 1.0000 | 9 | 2 | -- | -2.000 | 0.197 | $1,2,12,78,546$ | $1/2 - \tfrac{1}{4}\zeta(2)$ |

*Period column.* `ident.gp` ran to completion on all $34$ rows whose limit is known to
$\ge20$ significant digits: `lindep` against $17$ constants and, per row, against
$L(f,s)$ for every newform of weight $w+2$ at every level dividing $M$ or equal to
$2M,3M,4M$ — singly, and in pairs with a constant. "--" means the row has too few reliable
digits (equal or complex-conjugate leading roots); "unidentified" means the search ran and
found nothing. `build_tables.py` re-keys identifications by the canonical $a_n$, so the
table survives a re-run of `dedup.py`.

Reading the table:

* **Rows 4 and 8** ($N=2,3$, $t=(\eta_2/\eta_1)^{24}$, $(\eta_3/\eta_1)^{12}$) have the
  largest budgets but a **double** leading root, so $a_n\sim\lambda_1^n\cdot\text{poly}$ and
  there is no second growth rate to separate — the two-lattice machinery has nothing to work
  with. They are the level-$2$/level-$3$ analogues of scan 1's $\Lambda=64,27$ degenerations.
* **Row 9 is the only row in the whole scan that is simultaneously new, has $k=2$, has
  $|\lambda_2|<1$, and has an identified cusp-form period.** It is §4's object.
* Rows 13–15, 17–19 (levels $12,15,16$) are weight-one rows over Atkin–Lehner quotients with
  budgets $0.4$–$0.7$ whose periods `lindep` does not identify. Row 20 ($N=10$, $M=40$) is
  identified as a **difference of two level-$20$ weight-three newform $L$-values**, which is
  the expected shape when the associated weight-$(w{+}2)$ form is not itself an eigenform —
  further evidence that the family's periods live in $L(S_{w+2},w{+}1)$.
* Row 56 of the full table ($N=14$, the Fricke quotient of $X_0(14)$, $w=1$) has period
  $-\tfrac18-\tfrac9{16}L(f_7,2)$ with $f_7=\eta_1^3\eta_7^3\in S_3(\Gamma_0(7),\chi_{-7})$
  — the discriminant $-7$ member of the CM family of §7.
* Every row with an identified period and $k\ge2$ carries either a Dirichlet $L$-value
  ($\zeta(2)/4$, $L(2,\chi_{-3})/2$, $5L(2,\chi_{-3})/8$, $G/2$, $L(3,\chi_{-3})/3$), a
  $\zeta(3)$ value ($\zeta(3)/6$, $7\zeta(3)/24$, $7\zeta(3)/32$) or a weight-three CM
  cusp-form $L$-value ($L(f_8,2)$, $L(f_{12},2)$, $L(f_7,2)$).

---

# PART II — the $\operatorname{Sym}^1$ census (the main new result)

While ranking the scan output it became clear that the two "new" weight-one rows
(scan 1's $L(f,2)$ row on the Domb curve, and §4's row on the $T$ curve) are
instances of one construction. Running that construction on the whole
Almkvist–Zudilin list gives a complete, checkable table.

## 7. Every sporadic $\operatorname{Sym}^2$ row has an integral $\operatorname{Sym}^1$ square root, and it costs $k=2$ instead of $k=3$

Let $A(t)=\sum_{n\ge0}A_nt^n$ be one of the six order-three sporadic rows
$$(n+1)^3A_{n+1}=(2n+1)(an^2+an+b)A_n-c\,n^3A_{n-1}.$$
Form the power-series square root $\sqrt{A(t)}=\sum c_nt^n$ and rescale,
$$a_n:=\lambda^n c_n .$$

**Observation (verified, $n\le260$; for Apéry to $n\le1200$).** For each of the six rows
there is a $\lambda\in\{1,4\}$ with $a_n\in\mathbf Z$ for all tested $n$, and $a_n$ then
satisfies a **three-term recurrence with quadratic coefficients**
$$\boxed{\ (n+1)^2a_{n+1}=\Bigl(An^2+\tfrac A2 n+B\Bigr)a_n-C(2n-1)^2a_{n-1},\qquad a_0=1,\ a_1=B\ }$$
whose Casoratian is $a_nb_{n+1}-a_{n+1}b_n=\pm\,C^{\,n}\bigl((2n-1)!!\bigr)^2/\bigl((n+1)!\bigr)^2$
(for $C=4$ this is exactly the square of the Catalan number $\mathrm{C}_n$), and whose
companion needs only $d_n^2$:

| $\operatorname{Sym}^2$ parent | period of parent | $\lambda$ | $(A,B,C)$ | $\lambda_1$ | $\lambda_2$ | $c$ | $k$ | $\sigma_p$ | score | budget | $\operatorname{Sym}^1$ period |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Apéry $(17,5,1)$** | $\zeta(3)/6$ | $4$ | $(136,10,4)$ | $4(17{+}12\sqrt2)=135.882$ | $4(17{-}12\sqrt2)=0.11775$ | $16$ | $\mathbf 2$ | none | $\mathbf{+0.1392}$ | $+2.912$ | **unidentified**, $=0.1001874492\ldots$ |
| $T\,(12,4,16)$ | $7\zeta(3)/32$ | $1$ | $(24,2,4)$ | $12{+}8\sqrt2=23.314$ | $12{-}8\sqrt2=0.68629$ | $16$ | $\mathbf 2$ | none | $-1.624$ | $+1.149$ | $L(f_8,2)$, $f_8=\eta_1^2\eta_2\eta_4\eta_8^2$ |
| Domb $(10,4,64)$ | $7\zeta(3)/24$ | $1$ | $(20,2,16)$ | $16$ | $4$ | $64$ | $\mathbf 2$ | $\sigma_2{=}2$ | $-3.386$ | $+0.773$ | $L(f_{12},2)$, $f_{12}=\eta_2^3\eta_6^3$ |
| AZ $(9,3,-27)$ | $L(3,\chi_{-3})/3$ | $4$ | $(72,6,-108)$ | $77.569$ | $-5.569$ | $-432$ | $\mathbf 2$ | $\sigma_3{=}3$ | $-3.717$ | $+2.351$ | unidentified, $=0.1455144820\ldots$ |
| AZ $(11,5,125)$ | — | $4$ | $(88,10,500)$ | complex, $\lvert\lambda\rvert=44.721$ | (same) | $2000$ | $\mathbf 2$ | $\sigma_5{=}3$ | n/a | $+1.801$ | none (complex roots) |
| AZ $(7,3,81)$ | — | $4$ | $(56,6,324)$ | complex, $\lvert\lambda\rvert=36$ | (same) | $1296$ | $\mathbf 2$ | $\sigma_3{=}4$ | n/a | $+1.584$ | none (complex roots) |

Three independent checks that this table is right:

1. **The two rows already in the project's census appear in it, with their known data.**
   The Domb line reproduces scan 1's $L(f,2)$ row exactly (sequence $1,2,12,104,1078,\dots$,
   recurrence $(20n^2+10n+2)$, $16(2n-1)^2$, $k=2$, limit $L(\eta_2^3\eta_6^3,2)$ to
   $420$ digits). The $T$ line is §4's row.
2. **An independent search finds the same list.** A Zagier-style integrality sweep over the
   whole normalisation above — $|A|\le4000$ even, $|B|\le60$, $C$ in a list of $27$ values,
   testing $a_n\in\mathbf Z$ for $n\le45$ exactly (`lattice/sporadic_scan2/family_scan.py`)
   — returns, apart from the scalings $t\mapsto t/\mu$ (which send $(A,B,C)\mapsto(\mu A,\mu B,\mu^2C)$
   and $a_n\mapsto\mu^na_n$) and the degenerate double-root families
   $(8,B,4)$, $(16,B,16)$, $(24,B,36)$, $(32,B,64)$ with $B=2(2j^2+2j+1)\cdot$scale
   (these are $\binom{2n}{n}$ and $(4n+1)\binom{2n}{n}$ and their rescalings),
   **exactly the six triples of the table and nothing else.**
3. Every $a_n$ and every $d_n^2b_n$ was recomputed with exact rational arithmetic.

**Why this matters: it is a free integration.** Problem 7.6 of the paper says
"each free integration lowers $k$ by one, which is worth more than any $p$-adic device".
Passing from $\operatorname{Sym}^2$ to $\operatorname{Sym}^1$ on the *same* curve is exactly
such a step. The bookkeeping is exact: the square root has
$\lambda_1^{\mathrm{Sym}^1}=\lambda\,\lambda_1$, $\lambda_2^{\mathrm{Sym}^1}=\lambda\,\lambda_2$,
so
$$\text{score}_{\mathrm{Sym}^1}=\text{score}_{\mathrm{Sym}^2}+1-\log\lambda,\qquad
\text{budget}_{\mathrm{Sym}^1}=\text{budget}_{\mathrm{Sym}^2}+1+\log\lambda .$$
The square root is a strict improvement in score exactly when $\lambda\le2$.
For Apéry's curve the arithmetic is
$$\text{Sym}^2:\ \log\tfrac1{0.02944}-3=+0.5255,\qquad
\text{Sym}^1:\ \log\tfrac1{0.11775}-2=+0.1392 .$$
Both are positive — the $\operatorname{Sym}^1$ row of Apéry's curve is the **second**
positive-score integral row this project has seen, and the first that is not Apéry's own.

**The caveat, stated plainly, and no claim is made.** The archimedean limit of the Apéry
$\operatorname{Sym}^1$ row is
$$\xi=\lim_{n\to\infty}\frac{b_n}{a_n}
=\sum_{m\ge0}\frac{\mathrm C_m^{\,2}}{a_ma_{m+1}}
=0.10018744922933940616775868213306112162877772050081221171621939873497\ldots$$
($\mathrm C_m$ the Catalan numbers; the series identity is the telescoped Casoratian, hence
exact — and it is what makes the irrationality criterion applicable at all, since it shows
$a_n\xi-b_n\ne0$ for every $n$). `lindep` at $110$–$560$ digits finds **no** relation
between $\xi$ and $1$, $\zeta(2)$, $\zeta(3)$, $\zeta(4)$, $\pi^2$, $\pi^3$, $G$,
$L(s,\chi_D)$ for $D=-3,-4,-7,-8,5,8,12,-11,-24$ and $s=2,3,4$, $\log^32$, $\pi^2\log2$,
$\sqrt2$, $\sqrt3$, $\log(1+\sqrt2)^{1,2,3}$, or against $L(f,s)$ for **any**
newform $f$ of weight $3$ (all $1792$ real parts at the $31$ levels $3\le L\le72$ of
`ident_sym1.gp`, $s=1,2$) or weight $3,4$ at any level dividing $64$ or in
$\{3,4,12,24,48\}$ ($s\le k-1$), singly or in pairs with a constant. **So the number the positive score refers to has no name yet.** Nothing is
claimed about its arithmetic here; the row is reported as data.

## 8. Cooper's three rows have $\operatorname{Sym}^1$ square roots too

The same construction applied to Cooper's $s_7,s_{10},s_{18}$ (which are order-three but not
in Almkvist–Zudilin normalisation) again produces integral three-term rows with quadratic
coefficients — with a *more general* trailing coefficient, a product of two linear forms
rather than a square:

| parent | parent period | $\lambda$ | $\operatorname{Sym}^1$ recurrence $(n+1)^2a_{n+1}=p_1(n)a_n-p_0(n)a_{n-1}$ | $a_n$ | $\lambda_1$ | $\lambda_2$ | $k$ | score | budget | period |
|---|---|---|---|---|---|---|---|---|---|---|
| $s_7$ | $\zeta(2)/7$ | $1$ | $p_1=26n^2{+}13n{+}2$, $p_0=-3(3n{-}1)(3n{-}2)$ | $1,2,22,336,6006,\dots$ | $27$ | $-1$ | $2$ | $-2.000$ | $+1.296$ | $L(f_7,2)$, $f_7=\eta_1^3\eta_7^3$ |
| $s_{10}$ | $\zeta(2)/5$ | $2$ | $p_1=24n^2{+}12n{+}2$, $p_0=-4(8n{-}3)(8n{-}5)$ | $1,2,34,588,12726,\dots$ | $32$ | $-8$ | $2$ | $-4.079$ | $+1.466$ | unidentified $0.31692535921\ldots$ |
| $s_{18}$ | $\tfrac12L(2,\chi_{-3})$ | $2$ | $p_1=56n^2{+}28n{+}6$, $p_0=12(8n{-}3)(8n{-}5)$ | $1,6,90,1716,36774,\dots$ | $32$ | $24$ | $2$ | $-5.178$ | $+1.466$ | unidentified $0.48423755360\ldots$ |

**A CM pattern.** Where the $\operatorname{Sym}^1$ period is identified, the weight-one form
is the theta series of an imaginary quadratic order of class number one and the period is
$L(f,2)$ for the corresponding weight-three CM newform, *exactly on the nose*:

| $\operatorname{Sym}^1$ row | disc | weight-one form | weight-three CM newform | period |
|---|---|---|---|---|
| Cooper $s_7$ | $-7$ | $\theta$, level $7$ | $f_7=\eta_1^3\eta_7^3\in S_3(\Gamma_0(7),\chi_{-7})$ | $L(f_7,2)=0.46721176888\ldots$ |
| $T$ | $-8$ | $\theta_{x^2+2y^2}$, level $8$ | $f_8=\eta_1^2\eta_2\eta_4\eta_8^2\in S_3(\Gamma_0(8),\chi_{-8})$ | $L(f_8,2)=0.52874945559\ldots$ |
| Domb | $-3$ | $\eta_1^2\eta_3^2/\eta_2\eta_6$, level $12$ | $f_{12}=\eta_2^3\eta_6^3\in S_3(\Gamma_0(12),\chi_{-3})$ | $L(f_{12},2)=0.73729299619\ldots$ |

(each verified by `lindep` at $100$ digits against $17$ constants **and** all $1792$ real
parts of $L(f,s)$, $s=1,2$, for every weight-three newform at the $31$ levels
$3,\dots,72$ listed in `ident_sym1.gp` — the relation found is $\pm1$ in every case, no
rational factor.) The four rows whose weight-one form is **not** an integral CM theta series
— Apéry's, AZ$(9,3,-27)$'s, Cooper $s_{10}$'s and $s_{18}$'s — are exactly the four whose
periods the same $1792$-element search does **not** identify. That is a sharp, falsifiable
dichotomy and it is the single most useful thing to attack next.

So **all nine** sporadic order-three rows in the project's census admit an integral
$\operatorname{Sym}^1$ square root. The general shape is

$$(n+1)^2a_{n+1}=\Bigl(An^2+\tfrac A2n+B\Bigr)a_n-p_0(n)\,a_{n-1},\qquad
p_0(n)=C\,(\alpha n-\beta)(\alpha n-\beta'),$$

$a_0=1$, $a_1=B$, with $p_1$ **always** of the special form $A n^2+\tfrac A2 n+B$
(i.e. $p_1(n)=p_1(-n-\tfrac12)$: symmetric about $n=-\tfrac14$, the $\operatorname{Sym}^1$
analogue of Zagier's $p_1(n)=p_1(-n-1)$).

**The general relation between a row and its square root** (immediate, and it explains the
whole table): if the parent has $(\lambda_1,\lambda_2,k)$ and the square root needs the scale
$\lambda$, then the square root has $(\lambda\lambda_1,\lambda\lambda_2,k')$ and
$$\text{score}_{\mathrm{Sym}^1}=\text{score}_{\mathrm{Sym}^2}+(k-k')-\log\lambda .$$
For the six Almkvist–Zudilin rows $k=3\to k'=2$, so the square root **gains one** and pays
$\log\lambda$; it is a strict improvement exactly when $\lambda\le2$. Cooper's rows already
enjoy a free integration ($k=2$), and their square roots stay at $k'=2$, so they only pay
$\log\lambda$ and never improve. Apéry's is the only case where $\lambda=4$ and the parent's
score was positive, so both are positive and the parent is better ($+0.5255$ vs $+0.1392$).

## 9. The scan's own negative results

Over $52\,160$ bilinear solves for $\deg t\le2$ at $q$-precision $205$ (complete),
$493\,670$ for $\deg t\in\{3,4\}$ at the same precision (complete), and a partial
third pass at precision $400$ with a larger sublattice at levels $\le36$ (see the caveat at
the end of this section) — **522 hits, 85 distinct integral rows after de-duplication**
(weights $1{:}59$, $2{:}17$, $3{:}9$; parameter levels
$2,3,4,5,6,7,8,9,10,12,14,15,16,20,32$) — the scan says:

1. **No new $\operatorname{Sym}^2$ or $\operatorname{Sym}^3$ row.** Every order-$3$
   (weight-two) row found is a known one or a Möbius reparametrisation of one; every
   weight-three row found has $\lambda_1=\lambda_2$ (double root) and $k=4$, budget $<0.16$.
   Beauville–Zagier's exhaustion in weight one, and Almkvist–Zudilin's in weight two, are
   not contradicted anywhere in range.
2. **Only three rows in the entire scan have $|\lambda_2|<1$**: Apéry's, $T$'s, and the
   new $\operatorname{Sym}^1$ row on the $T$ curve (§4). Decay is *extremely* rare.
   Exactly one of them has $k=2$: §4's.
3. **The only positive score inside the modular scan is Apéry's own $+0.5255$.** (The
   $+0.1392$ of Part II is *not* inside the modular scan: its parameter is $s/4$ with $s$
   Apéry's Hauptmodul, so $t\ne q+O(q^2)$ and the scan's normalisation excludes it. That is
   a coverage lesson: **requiring $t=q+O(q^2)$ integral loses rescaled rows**, and the
   rescaling is precisely what makes the $\operatorname{Sym}^1$ square roots integral.)
4. **Raising the level does not raise $\lambda_1$**, confirming scan 1. The largest
   $\lambda_1$ of any row with a genuine second solution and a non-rational period is
   Apéry's $33.97$ inside the modular scan; the family search of Part II reaches
   $135.88$, but only by rescaling Apéry's own curve.
5. **The precision-$400$ pass found nothing the precision-$205$ pass had missed.** A third
   pass was run at $q$-precision $400$ with $K\le24$ — a strictly larger sublattice of
   $\Lambda_{M,w,\chi}$ than the main pass could afford — over the levels $\le36$ for which
   `prep600/` was built. It produced $195$ hits and $190$ verified rows, and **every one of
   them is a sequence already in the 85**: zero new distinct rows. That is evidence (not
   proof) that the sublattice truncation of §2(3) is not where the interesting rows hide,
   at least at low level. *Honest caveat*: this pass is **incomplete** — its first attempt
   died on a malformed prep file (a `.part` written twice) and the rerun was stopped
   part-way to free cores for the identification, so it covered roughly the first third of
   the $\deg t\le2$ parameters.
6. **The weight-one rows over Atkin–Lehner quotients form a visible family.** Beyond the
   six Zagier rows the scan produces integral weight-one rows over the degree-two
   parameters at $N=6,8,12,14,15,16$; their $\lambda_1$ are
   $16$, $12{+}8\sqrt2$, $8{+}4\sqrt3$, $7{+}4\sqrt3$, $8$, $\tfrac12(11{+}\sqrt{125})$, $6{+}4\sqrt2$; where the period is identified it is always $L(f,2)$ for a weight-three
   CM newform ($f_8=\eta_1^2\eta_2\eta_4\eta_8^2$ for $\mathbf Q(\sqrt{-2})$;
   $f_{12}=\eta_2^3\eta_6^3$ for $\mathbf Q(\sqrt{-3})$; $f_7=\eta_1^3\eta_7^3$ for
   $\mathbf Q(\sqrt{-7})$, at $N=14$). Several of the level-$12$, $15$, $16$ limits are
   **unidentified** (§6 table).

## 10. What to do next (in order of value)

1. **Identify $\xi=0.10018744922933940616775868213306112162877772050081\ldots$**, the Apéry
   $\operatorname{Sym}^1$ period. This is the sharpest question the scan produces: it is the
   period of a row whose score is $+0.1392>0$ and whose $a_n$ are integers with a clean closed
   form ($a_n=4^n[t^n]\sqrt{\sum A_nt^n}$, $A_n$ Apéry's numbers), and it resists a
   several-hundred-digit `lindep` against everything in the project's basis. Natural next
   places to look: $L$-values of the *weight-three form with quadratic multiplier* attached to
   $\sqrt{F_{\text{Apéry}}}$ (i.e. of level $24$ or $36$ with a quartic nebentypus, which the
   present `mfinit` sweep only partly covers); Eichler integrals of $F_{\text{Apéry}}^{1/2}$;
   or a $\Gamma$-value/Mahler-measure expression.
   The same question is open for the four other unidentified $\operatorname{Sym}^1$ periods
   ($0.14551448\ldots$, $0.46721177\ldots$, $0.31692536\ldots$, $0.48423755\ldots$).
2. **Prove the square-root integrality.** "$\lambda^n[t^n]\sqrt{A(t)}\in\mathbf Z$ with
   $\lambda\in\{1,2,4\}$ for every sporadic order-three row" is, as far as this project's
   ledger goes, a new observation. It is a statement about $\operatorname{Sym}^2$ local
   systems descending to $\operatorname{Sym}^1$ over $\mathbf Z[\tfrac1\lambda]$, and it is
   the mechanism behind Problem 7.6's free integration.
3. **Cover the missing genus-zero Fricke levels.** $N=11,17,19,23,29,31,41,47,59,71$ (and
   $34,38,49,51,55,62,69,84,87,94,95,105,110,119$) have genus-zero Atkin–Lehner quotients whose
   Hauptmoduls are McKay–Thompson series, not eta quotients, so they are invisible to both
   scans. The class-number-one pattern of §4 (disc $-3\to$ Domb curve, $-7\to$ level $14$,
   $-8\to T$ curve) suggests $\mathbf Q(\sqrt{-11})$ at level $11$ and $\mathbf Q(\sqrt{-19})$
   at level $19$ as the next members, with larger $\lambda_1$ and hence larger budget.
4. **Drop the $t=q+O(q^2)$ normalisation from the modular scan.** Allowing $t=q/\lambda+\cdots$
   is exactly what the $\operatorname{Sym}^1$ square roots need; the present scan's largest
   miss (the $+0.1392$ row) is entirely due to this.
5. **A partner for $L(f_8,2)$.** §4's row is a *decayer* ($\lambda_2=12-8\sqrt2<1$, $k=2$) for
   the level-$8$ weight-three CM period, but has $\sigma_p=0$ everywhere, so it cannot be an
   engine. A row with period $L(f_8,2)$ and $\sigma_2>0$ would complete a Catalan-type
   two-lattice pair for a cusp-form period — the first such configuration where the decaying
   half is already in hand.

## 11. Two closing remarks on completeness

* **The $\operatorname{Sym}^1$ census is complete relative to the $\operatorname{Sym}^2$ one.**
  The construction of §7 takes an order-three row and returns an order-two one; since
  Almkvist–Zudilin's list of order-three sporadic rows in their normalisation is believed
  exhaustive, and Cooper's three are the known additions, the table of §7+§8 is the whole
  output of the construction. The independent $(A,B,C)$ sweep of `family_scan.py` confirms
  this inside the normalisation $p_0=C(2n-1)^2$: it returns the six Almkvist–Zudilin square
  roots, their rescalings $t\mapsto t/\mu$, and the degenerate double-root families
  $\binom{2n}{n}$, $(4n{+}1)\binom{2n}{n}$ — nothing else, for $|A|\le4000$, $|B|\le60$ and
  $27$ values of $C$.
* **The converse map is not a bijection.** Squaring one of Zagier's six weight-one rows also
  produces an order-three integral row, but not in Almkvist–Zudilin normalisation (the
  trailing coefficient is not $cn^3$), so those squares are not in the census and the
  $\operatorname{Sym}^1\leftrightarrow\operatorname{Sym}^2$ correspondence is one-way as
  stated here. Classifying the order-three rows with trailing coefficient
  $n\,(\alpha n-\beta)(\alpha n-\beta')$ — Cooper's shape — is a well-posed follow-up that
  this scan did not attempt.
## 12. Reproduction

```
cd lattice/sporadic_scan2
python3 gen_t.py                        # 3979 parameters, all N <= 120     (~40 s)
gp -q w_c{0..7}.gp                      # integral form lattices  -> spaces/  (~30 min, 8 cores)
gp -q p_c{0..7}.gp                      # LLL-reduced cosets      -> prep/
python3 scan2.py --degmax 2 --levcap 4 --nshard 8 --shard $i --out hitsA
python3 verify.py --hits 'hitsA_*.jsonl' --out rowsA_$i.jsonl
python3 dedup.py --rows rowsA_*.jsonl rowsB_*.jsonl rowsC_*.jsonl
gp -q ident.gp                          # periods -> ident.out
python3 build_tables.py                 # -> table.md, top_table.md
python3 sym1.py                         # PART II: the Sym^1 census -> sym1.json
python3 family_scan.py                  # PART II: independent (A,B,C) sweep
```

