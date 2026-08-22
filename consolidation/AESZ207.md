# AESZ 207: the constant, corrected — and what it is not

*Claude (Fable), 2026-08-22.  Scripts: `lattice/followups/b1_probe.gp` … `b21_xi2digits.gp`.
Follows up `MUM_SURVEY.md` §5.6 ("AESZ 207 … a $2$-adic limit known to 2449 digits that is not a
single Kubota–Leopoldt value").  Tags: **[verified]** exact rational / $p$-adic computation over
a stated range; **[measured]** floating-point; **[lit]** located in the literature.*

---

## 0. Verdict

1. **AvSZ's printed value for #207 is wrong, and we give the right one to 1139 digits.**
   `[lit]` Almkvist–van Straten–Zudilin (*Apéry limits …*, Fields Inst. Commun. **54** (2008))
   §4.8 "Unidentified limits with many digits" lists seven, #207 among them (with 17, 34, 147,
   206, 214, 229), printing
   $-0.00050462505145900474057831709244307528529622730007723$
   (PDF located and read by a literature sub-agent:
   `download.uni-mainz.de/mathematik/Algebraische Geometrie/DvS-Publikationen/Apery Limits of
   Differential Equations of Order 4 and 5.pdf`; the operator printed for #207 in AESZ Table A,
   arXiv:math/0507430 v1 and v2, is coefficient-for-coefficient the CYCluster `4.4.38` we use,
   so the $z$-normalisation agrees).  The true limit is
   $$\boxed{\ \xi_\infty(207)=-0.00050455459344136708862542545797071516117215446057914345028\ldots\ }$$
   **[verified]** to $1139$ digits ($N=5000$, `b14_hiprec.gp`; file
   `lattice/followups/b14_xi_arch_1200.txt`).  The two disagree from the 4th significant figure
   ($7\cdot10^{-8}$ absolute).  The reason is visible: AvSZ write "these limits converge very
   fast", but for #207 $|\lambda_2/\lambda_1|=53248/89531.3892=0.5947411$, i.e. **0.2257 digits per
   index** — the slowest in their table.  Their fifty digits are un-converged.  *This alone
   explains why no PSLQ identification was ever found: the target number was wrong.*
   Our pipeline is validated: it reproduces AvSZ's own identified entries
   (AESZ 16, 28, 29, 42, 58, 182, 185 $=\tfrac7{48},\tfrac17,\tfrac1{12},\tfrac7{64}$ of
   $\zeta(3)$, $\tfrac18L(\chi_{-3},2)$, $\tfrac3{11}\zeta(3)$, $\tfrac16L(\chi_{-3},3)$) to
   $10^{-347}$ (`b11_valid.gp`), and AvSZ's #17 and #34 to 50 and 31 digits.

2. **The geometry, located.**  `[lit]` AESZ Table A: #207 *is the reflection at $z=\infty$ of
   AESZ 99 / 4.4.34*, the degree-13 $5\times5$-Pfaffian in $\mathbf P^5$ ($H^3=13$, $c_2H=58$,
   $\chi=-120$, so $\chi/H^3=-120/13$); "a formula for $A_n$ is not known".  #207 itself carries
   **no** topological data in CYCluster or the Mainz CYDB.  Riemann scheme:
   $z=0$ MUM $(0,0,0,0)$; $z=1/53248$, multiplicity **2** in the discriminant, exponents
   $(0,1,3,4)$; the conjugate pair $z_\pm=(349\pm85\sqrt{17})/2^{17}$, exponents $(0,1,1,2)$ —
   **conifolds over the real quadratic field $\mathbf Q(\sqrt{17})$**; $z=\infty$
   $(\tfrac12,\tfrac12,\tfrac12,\tfrac12)$.  We re-derived the whole scheme independently
   (`b7_local.gp`): $q_4(z)=z^4(53248z-1)^2(2^{24}z^2-89344z-1)$.

3. **The two places see different singularities.**  The characteristic roots are
   $$\lambda=-89531.3892\ (=-2^7(349+85\sqrt{17})),\qquad 53248=2^{12}\cdot13\ \text{(double)},
     \qquad 187.3892,$$
   so the archimedean limit is governed by the $\mathbf Q(\sqrt{17})$ conifold
   ($\lambda_1=-2^7(349+85\sqrt{17})$, $v_2=8$): $A_n\sim\lambda_1^n$, the linear form
   $B_n-\xi_\infty A_n$ grows only like $53248^n$, so $B_n/A_n-\xi_\infty\sim0.59474^n$.  The measured
   $2$-adic slope is $\sigma_2=12=v_2(53248)$ **exactly** — the *rational* singularity with the
   non-conifold exponents $(0,1,3,4)$.  **[measured]** increments $v_2(\xi_N-\xi_{N-1})$ per 200
   indices: $11.99,11.99,11.985,11.975,12.035,11.965,12.005$ for $N\le1600$, and
   $v_2(\xi_{1600}-\xi_{1599})=19156$.  $\kappa_2=0$ ($v_2(A_n)\in\{12,\dots,24\}$ for
   $n\le1600$).  No slope at $13$ or $17$ (bounded negative increments).

4. **Not identified — but the exclusions are now sharp.**  With $1139$ archimedean digits and
   $5958$ $2$-adic digits (§2, §3), $\xi_\infty$ and $\xi_2$ are **not**:
   a rational multiple, nor a $2$-term rational combination, of any Dirichlet $L$-value
   $L(\chi,m)$ / $L_2(m,\chi\omega^{j})$ in a catalogue of $\sim\!350$ archimedean and $188$
   $2$-adic targets; not in the $\mathbf Q$- or $\mathbf Q(\sqrt{17})$-span of the
   $\widehat\Gamma$-class basis $\{1,\Lambda,\Lambda^2,\Lambda^3,\zeta(2),\zeta(2)\Lambda,
   \zeta(3)\}$ at the dominant conifold, at either place; not algebraic of degree $\le8$ (arch.)
   or $\le6$ ($2$-adic) at any reachable height; and $\xi_2\notin\mathbf Q(\sqrt{17})$ (which is
   a subfield of $\mathbf Q_2$, since $17\equiv1\bmod 8$).  Details and the full catalogue in §4.

5. **The trichotomy cannot be used here — a negative result.**  The good-prime tower
   eigenvalue $\chi(p)p^{-w}$ of `GOOD_PRIME_TOWERS.md` §6 does **not** exist for AESZ 207:
   the ratios $t_{s+1}/t_s$, $t_s=B_{ap^s}/A_{ap^s}$, do not converge, and their valuations
   are not even constant — at $p=3$ they alternate $-3,-4,-3,-4,\dots$ (mean $-\tfrac72$,
   **non-integral**), at $p=5$ they settle at $-3$ and at $p=7$ at $-4$, with unit parts that
   agree only mod $p$ ($1$, $3$, $5$ resp.) and diverge thereafter.  So **$\chi(3),\chi(5),\chi(7)$
   cannot be read off**, and the constraint the task hoped for is unavailable.  §5.

---

## 1. The operator

CYCluster `4.4.38` = AESZ 207, $\deg_z=4$:
$$L=\theta^4-2^4z\bigl(1072\theta^4-17824\theta^3-10888\theta^2-1976\theta-145\bigr)
 -2^{17}z^2\bigl(51088\theta^4+116368\theta^3-45264\theta^2-14228\theta-1397\bigr)$$
$$\qquad+2^{28}\!\cdot\!13\,z^3\bigl(73104\theta^4+1536\theta^3-488\theta^2+384\theta+97\bigr)
 -2^{44}\!\cdot\!13^2z^4(2\theta+1)^4 .$$
$A_n:\;1,-2320,57601296,-2373661139200,121665506430000400,\dots$;
$B_1=1$, $B_2=-29761$, $B_3=\tfrac{32129592496}{27}$, ….
Denominator exponent $k=4$ **[verified]** $n\le200$ (so the "Hodge depth" is $w+1=4$; note $k=4$
is the modal value in the rank-4 census, so it carries less weight-information than
`MUM_SURVEY.md` §5.7 suggests).

**Discriminant and roots.**  $q_4(z)=z^4(53248z-1)^2(2^{24}z^2-89344z-1)$;
$53248=2^{12}\cdot13$, $89344=2^8\cdot349$, and
$349^2-85^2\cdot17=-1024$, so $z_\pm=(349\pm85\sqrt{17})/2^{17}$ with $z_+z_-=-2^{-24}$.
Characteristic polynomial of the recurrence
$$x^4-17152x^3-6696206336x^2+255108172480512x-47569271064100864=(x-53248)^2(x^2+89344x-2^{24}),$$
roots $-2^7(349+85\sqrt{17})=-89531.389$, $2^{17}/(349+85\sqrt{17})=187.389$, and $53248$ twice.
Hence $\lambda_1=89531.389$, $\lambda_2=53248$, and
$$|\lambda_2/\lambda_1|=0.5947411,\qquad \text{$0.225672$ decimal digits of $\xi_\infty$ per index.}$$

---

## 2. The archimedean constant, to 1139 digits

`b14_hiprec.gp`, $N=5000$, exact rationals, $\mathrm{realprecision}=1250$; last increment
$10^{-1139.2}$.  First 300 digits:

```
-0.000504554593441367088625425457970715161172154460579143450288368704004362680122570507
 26641683575012642830666998981277285243971074365915361974024684407837919073440763147804
 05944481004391073559412495789180024359354034004201901516003082047986788940972950004267
 322872427242641273929807512694076104095696976391767670762468981380422676628674381521
```
(full value: `lattice/followups/b14_xi_arch_1200.txt`).  $1/\xi_\infty=-1981.946082740811\ldots$

---

## 3. The $2$-adic constant

`b19_alg.gp`/`b21_xi2digits.gp`, $N=500$: Cauchy precision $v_2(\xi_N-\xi_{N-1})=5958$, so
$\xi_2$ is known to $5958$ $2$-adic digits (the survey's 2449 was an artefact of its $0.85\sigma N$
cap at $N=240$; at $N=1600$ we have $19156$).
$$v_2(\xi_2)=-4,\qquad 16\,\xi_2\equiv 1275528383008663815257650142033548687 \pmod{2^{120}},$$
$16\xi_2\equiv 360776799631\pmod{2^{40}}$.

---

## 4. What $\xi$ is not

### 4.1 Archimedean (1100–1139 digits, `b3`,`b4`,`b6`,`b8`,`b13`,`b15b`,`b17`,`b18`)

Rejected as a rational multiple **and** as a $2$-term rational combination (heights allowed up to
$10^{9}$–$10^{20}$, residual test at $\ge280$ digits; `b11_valid.gp` plants known relations of the
same shape and recovers them, so these are real negatives):

* $\zeta(3),\zeta(5),\zeta(7)$, $\pi^{2},\pi^{3},\pi^{4},\pi^{5},\pi^{6}$, $\pi^2\zeta(3)$,
  $\zeta(3)^2$, $G$, $G\pi$, $\mathrm{Li}_4(\tfrac12)$, $\mathrm{Li}_3(\tfrac12)$,
  $\log2,\log^22,\log^32,\log^42$, $\zeta(3)\log2$, $\pi^2\log^22$, $\log13$, $\log17$;
* $L(\chi_D,m)$ for $m=2,3,4,5$ and $\sim70$ fundamental discriminants $D$ including all
  $D$ built from $\{2,13,17\}$ ($-4,8,-8,13,-52,104,-104,17,-68,136,-136,221,884,1768,\dots$)
  and everything with $|D|\le56$;
* **all** Dirichlet characters (any order, not just quadratic) of conductor
  $f\in\{3,4,5,7,8,11,12,13,16,17,20,24,26,32,34,39,40,51,52,68,104,136,221,272\}$,
  $m=2,\dots,6$, real *and* imaginary parts — $980$ targets, no single hit (`b15b.log`);
* all of the above times $\sqrt{17}$ or $\sqrt{13}$; and $1$ adjoined to every pair
  (i.e. $\xi=q_0+q_1C$ tested for every $C$ in the catalogue), and all triples over a
  $33$-element core;
* the $\widehat\Gamma$-class span at the dominant conifold,
  $\{1,\Lambda,\Lambda^2,\Lambda^3,\pi^2,\pi^2\Lambda,\zeta(3)\}$ with
  $\Lambda=\log|1/z_-|=11.402344560\ldots$, over $\mathbf Q$ (heights $\le10^{30}$) and over
  $\mathbf Q(\sqrt{17})$ (14 terms, heights $\le10^{25}$), with and without $\pi^4,\zeta(5)$;
  and the same allowing $\xi$ itself an algebraic or $\pi^2$ denominator
  ($\{\xi,\sqrt{17}\xi,\pi^2\xi,\pi^2\sqrt{17}\xi\}$ adjoined);
* $\mathrm{algdep}$ of degree $1..8$: all returned heights $>10^{25}$, i.e. **no algebraic
  relation of small height**;
* the same tests against $\log$ of the algebraic singular values
  $\alpha=349+85\sqrt{17}$, $\beta=349-85\sqrt{17}$, $1/53248$.

### 4.2 $2$-adic (900–5500 digits, `b5`,`b10`,`b17`,`b19`)

* $L_2(m,\chi)$ for $83$ quadratic characters ($D$ as above) and $m=2,3,4,5$ — $188$ nonzero
  targets — singly and in pairs: **no hit** (this extends `MUM_SURVEY.md` §5.6's 34 targets by a
  factor $5.5$, and adds all pairs);
* $\zeta_2(m)=L_2(m,\mathbf 1)$ for $m=2,\dots,7$; $L_2(m,\chi_8)$, $L_2(m,\chi_{13})$,
  $L_2(m,\chi_{17})$ for $m=2,3,4$, singly, in pairs, and with $\sqrt{17}$-coefficients;
* Iwasawa logarithms $\log_2u$ for $u=3,5,7,11,13,17,29,349,\pm$, their pairwise products and
  products with $\zeta_2(m)$;
* the $2$-adic $\widehat\Gamma$-class span
  $\{1,\Lambda_2,\Lambda_2^2,\Lambda_2^3,\zeta_2(2),\zeta_2(2)\Lambda_2,\zeta_2(3)\}$ with
  $\Lambda_2=\log_2$ of the unit part of $\alpha$ (resp. of $\beta$), over $\mathbf Q$ and over
  $\mathbf Q(\sqrt{17})$, and the $11$-term span using both logs;
* $\xi_2\in\mathbf Q$?  no.  $\xi_2\in\mathbf Q(\sqrt{17})$?  no (residual $v_2=5493$ against a
  relation of height $10^{500}$ — i.e. garbage).  $\mathrm{algdep}$ degree $\le6$: heights
  $\ge10^{220}$, no relation.

Useful arithmetic facts recorded along the way: $17\equiv1\pmod 8$, so $\sqrt{17}\in\mathbf Q_2$
and the two conifolds are $2$-adically rational, with
$v_2(349+85\sqrt{17})=1$, $v_2(349-85\sqrt{17})=9$, $v_2(z_-)=-8$;
$v_2(\zeta_2(2))=-1$, $v_2(\zeta_2(3))=-2$, $v_2(\log_2\alpha_{\rm unit})=v_2(\log_2\beta_{\rm unit})=2$.

### 4.3 Weight-4 newforms

`[lit]` BKSZ arXiv:2203.09426 treats only the fourteen hypergeometric families (conifold newform
levels $8,9,16,25,27,32,36,72,108,128,144,200,216,864$); **AESZ 207 appears nowhere**, and its
conifolds are quadratic-irrational, so a rational newform level is not expected in the naive
sense.  We nevertheless scanned $L(f,s)$, $1\le s<k$, over **all** newforms of weight $k\in\{2,4\}$ and
level $N\in\{2,4,8,13,16,17,26,32,34,52,68,104,136,221,272\}$ — **345 nonzero $L$-values**
(all Galois conjugates of every eigenform, `b16b.gp`, log `lattice/followups/b16c.log`,
$105$ digits): as rational multiples of $\xi_\infty$ (heights $\le10^{8}$), paired with
$\{1,\zeta(3),\pi,\pi^2,\pi^3,\pi^4,\sqrt{13},\sqrt{17}\}$, and in all $345\cdot344/2$ pairs
among themselves (heights $\le10^{15}$, residuals at $95$ digits).  **No hit in any of the three
sweeps.**

---

## 5. Frobenius / tower data: why the character cannot be pinned

`b18_ratio.gp`, $N=6000$, exact, $p=3,5,7$, tower bases $a=1,2,3$.  With
$t_s=B_{ap^s}/A_{ap^s}$:

| $p$ | $a$ | $v_p(t_s)$, $s=0,1,2,\dots$ | $v_p(t_{s+1}/t_s)$ | $p^3t_{s+1}/t_s \bmod p$ |
|---|---|---|---|---|
| 3 | 1 | $0,-3,-7,-10,-14,-17,-21,-24$ | $-3,-4,-3,-4,-3,-4,-3$ | $1,\ast,1,\ast,1,\ast,1$ |
| 3 | 3 | $-3,-7,-10,-14,-17,-21,-24$ | $-4,-3,-4,-3,-4,-3$ | — |
| 5 | 1 | $-1,-5,-8,-11,-14,-17$ | $-4,-3,-3,-3,-3$ | $3,3,3,3$ |
| 5 | 2 | $0,-4,-7,-10,-13$ | $-4,-3,-3,-3$ | $3,3,3$ |
| 7 | 1 | $0,-3,-7,-11,-15$ | $-3,-4,-4,-4$ | $5,\ast,\ast,\ast$ |

The trichotomy of `GOOD_PRIME_TOWERS.md` §6 predicts a *constant* $v_p(t_{s+1}/t_s)=-w$ and
$p^wt_{s+1}/t_s\to\chi(p)=\pm1$.  Neither holds:

* at $p=3$ the valuation alternates with mean $-\tfrac72$ — a **half-integral** slope, i.e. the
  cusp Frobenius eigenvalue is not in $\mathbf Q_3$;
* at $p=5$ and $p=7$ the valuation stabilises ($-3$, $-4$) but the unit parts agree only mod $p$
  ($\equiv3$ at $p=5$, $\equiv5$ at $p=7$ where the valuation permits a comparison at all; neither is $\pm1$) and then diverge — no limit;
* the denominators of $B_n$ dominate: $v_p(B_{p^s})$ grows like $-3.5s$, $-3s$, $-4s$ at
  $p=3,5,7$ against the bound $-4s$ from $d_n^4$.

**Diagnosis.**  The trichotomy is a statement about a rank-2 extension of $\mathbf Q(0)$ by an
Eisenstein line with a *single* second eigenvalue $\chi(p)p^{-w}$.  AESZ 207 is a rank-4
$\mathbf Q$-motive of weight 3 whose cusp Frobenius has four eigenvalues; $B$ is a mixture, and
no single eigenvalue governs the tower.  The prediction that "$\sigma_2=12$ plus the tower
eigenvalues at $p=3,5,7$ determine $\chi(p)$" is therefore **not available for this operator**,
and this is a genuine limitation of the tower method at rank 4, not a computational failure.

---

## 6. What is actually known, and the sharpest characterisation

> $\xi_\infty(207)=-0.000504554593441367088625425457970715161\ldots$ (1139 digits) and
> $\xi_2(207)=2^{-4}\cdot(\text{$2$-adic unit})$ (5958 digits) are the Apéry constants at the
> two places of the **second MUM point of the degree-13 $5\times5$-Pfaffian threefold in
> $\mathbf P^5$** (AESZ 99, $H^3=13$, $c_2H=58$, $\chi=-120$).  The archimedean constant is a
> period of the conifold $z_-=(349-85\sqrt{17})/2^{17}$, a conifold defined over
> $\mathbf Q(\sqrt{17})$; the $2$-adic slope $\sigma_2=12=v_2(2^{12}\cdot13)$ is carried by the
> *rational* singular point $z=1/53248$, whose exponents $(0,1,3,4)$ are **not** of conifold type.
> Neither constant is a Dirichlet $L$-value, a Kubota–Leopoldt value, a $\widehat\Gamma$-class
> combination at the conifold, or an algebraic number of small height, over $\mathbf Q$ or over
> $\mathbf Q(\sqrt{17})$.

### 6.1 Why the standard law is not expected to apply

`MUM_SURVEY.md` §5.2's law $\xi_p=r_\infty\mathcal E_p(m)^{-1}L_p(m,\chi\omega^{1-m})$ was
verified on operators whose Apéry limit is a *Dirichlet* (Eisenstein) period.  For AESZ 207 the
relevant degeneration is a conifold over a real quadratic field, i.e. the motive of the vanishing
cycle is not an Eisenstein/Dirichlet object.  Two consistent signals:

1. the archimedean and $2$-adic sides are governed by **different** singular points (§0.3);
2. the good-prime tower has no Eisenstein-type eigenvalue (§5).

### 6.2 The next computations, in priority order

1. **Do the connection problem.**  Following `lattice/mum_survey/conn_quintic.gp`, expand the
   local solutions at $z_-$ (exponents $0,1,1,2$; radius of convergence
   $|z_--z_+|=2\cdot85\sqrt{17}/2^{17}=0.0053477$, and $|z_-|=1.117\cdot10^{-5}$, so a matching
   point at $0.85z_-$ works) and express $\xi_\infty$ *explicitly* as a ratio of Frobenius
   constants.  This converts the identification problem from a blind PSLQ into the identification
   of a handful of individual constants, exactly as the quintic computation did.  This is the
   single highest-value next step and is entirely within reach.
2. **Hilbert modular / $\mathbf Q(\sqrt{17})$ $L$-values.**  If the conifold vanishing cycle is
   a Hilbert modular motive over $\mathbf Q(\sqrt{17})$, the target is $L(g,3)$ for a Hilbert
   newform of parallel weight $(4,4)$ and small level over $\mathbf Q(\sqrt{17})$ — not
   computable in PARI, but computable in Magma/Pari-HMF.  The exclusion of every rational
   newform level in §4.3 is evidence in this direction.
3. **AESZ 99 first.**  #99 has the topological data ($\chi/H^3=-120/13$) and a rational
   $\widehat\Gamma$-class; its own Apéry limit and $p$-adic limits should be measured, and the
   relation between the two MUM points (the "reflection at infinity") used to transport the
   answer.  #99 is `4.4.34` in `ops.gp` and was **not** part of this run.
4. **Tell AvSZ.**  The corrected value of §0.1 should be reported; likewise their #206, #214,
   #229 should be re-checked for the same under-convergence.

---

## 7. Scripts

`lattice/followups/`: `b1_probe.gp` (roots, convergence rate, slopes), `b2_main.gp` (limits,
towers), `b3`,`b4`,`b6`,`b8` (archimedean catalogues), `b5`,`b10` ($2$-adic catalogues),
`b7_local.gp`/`b7b.gp` (Riemann scheme, independently re-derived), `b9_struct.gp` (coefficient
structure), `b11_valid.gp` (**pipeline validation against AvSZ's identified entries — read this
before doubting any negative above**), `b12`/`b16b` (newform $L$-values), `b13`/`b17`/`b18`
($\widehat\Gamma$-class ansatz, ratio ansatz, towers), `b14_hiprec.gp` (the 1139-digit value),
`b15b.gp` (980 Dirichlet characters), `b19`/`b20`/`b21` (algebraicity, slope, digits).
Data: `b14_xi_arch_1200.txt`, `b2_xi_arch.txt`, `b2_xi_2adic.txt`.
