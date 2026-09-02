# Cooper's magnetic congruence, round 3: the vector-valued home of the $s_7$ weight-$5/2$ input, and the integral basis that round 2 was missing

*Working note, 2026-09-02.  PARI/GP 2.15.4.  All scripts and logs in
`lattice/cooper_congruence/round3/`; nothing outside this directory was created or modified.
Background read first: `consolidation/COOPER_CONGRUENCE.md` (v2), `lattice/cooper_congruence/round2/REPORT.md`
§§2.8, 5.4, 5.5, 8 and `round2/FINDINGS_PZ.md`.  Round 2's files were read, never written.*

Claims are tagged **[proved]**, **[verified, range]**, **[num, digits]**, **[refuted, range]**,
**[conjectural]**.

---

## 0. Verdict table

| # | statement | status |
|---|---|---|
| **A** | **THE INPUT IS ANTISYMMETRIC.**  The $s_7$ weight-$5/2$ datum is a *vector-valued* object $c(\beta,d)$, $\beta\in\mathbf Z/14$, and it satisfies $c(-\beta,d)=-c(\beta,d)$ — not $+$.  Proof: $W_7$ carries the Heegner class $\beta$ to $-\beta$, fixes the genus character, and $\widehat f\mid W_7=-\widehat f$ (round 2, R2/R3) | **[proved]** + **[verified, every admissible $d\le120$, both classes, 45 digits]** `09_beta.gp/.log` |
| **B** | **PROVED CONSEQUENCE: $c(d)=0$ for every admissible $d$ with $7\mid d$.**  $7\mid d\iff\beta\in\{0,7\}\iff\beta\equiv-\beta\ (14)$, and A then forces $c=0$.  This proves round 2's P6 observation "$G(d)=0$ exactly when $7\mid d$" (both directions) and **resolves all 25 `DEGEN` entries of `round2/73_cd_s7_2500.txt` — they are $0$** | **[proved]** |
| **C** | **THE REPRESENTATION IS FIXED BY MILGRAM.**  $\sum_{\beta\bmod14}e(\beta^2/28)=\sqrt{14}\,e(1/8)$, so $(\mathbf Z/14,\,x^2/28)$ has signature $\equiv1\pmod 8$.  Weight-$k$ forms are symmetric iff $2k\equiv s\ (4)$, antisymmetric iff $2k\equiv s+2\ (4)$.  With $2k=5$ and A the input transforms under the **dual** Weil representation $\rho_L^{*}$: exponents $|D|/28$ with $D\equiv-\beta^2\ (28)$, $D=D_0d=-3d$ the CM discriminant | **[proved]** `20_diag.gp` |
| **D** | **HENCE $M^{!,-}_{5/2}(\rho^*_L)=J^{!}_{3,7}$**, weakly holomorphic **Jacobi forms of weight $3$, index $7$** (odd Jacobi weight $\Leftrightarrow$ antisymmetric theta components).  This names the space round 2 §8(2)(iii) asked for | **[proved]** |
| **E** | **UNIQUENESS, FOR FREE.**  $J_{3,7}=0$: the holomorphic subspace of $\Delta^{-n}J^{\mathrm{weak}}_{3+12n,7}$ is $0$ for $n=0,1,2$ and the principal-part map is **injective** (rank $=\dim$: $7/7$, $13/13$, $19/19$).  A weakly holomorphic form is determined by its principal part, with **no auxiliary hypothesis** — stronger and cleaner than round 2's P4.  Agrees with Skoruppa–Zagier: $J^{\mathrm{cusp}}_{3,7}\cong S_4(\Gamma_0(7))^{W_7=-1}=0$ | **[verified, $n\le2$]** + **[proved, modulo Skoruppa–Zagier]** `20_diag.gp` |
| **F** | **THE INTEGRAL BASIS — round 2's "single genuinely missing ingredient" — IS CONSTRUCTED.**  $J^{!}_{3,7}=\bigcup_n\Delta^{-n}\varphi_{-1,2}\bigoplus_{b=0}^{5}M_{4+2b+12n}\varphi_{0,1}^{5-b}\varphi_{-2,1}^{b}$.  At pole level $n=1$ ($\dim 13$) the basis $g_{D_0}=q^{D_0/28}(e_\beta-e_{-\beta})+O(q^{1/28})$ exists for $13$ of the $17$ admissible $D_0\ge-72$ and has **integer** coefficients (denominator $1$ over the first $24$ slots in every case) | **[constructed]** + **[verified]** `15_jaclib.gp`, `22_basis.gp/.log` |
| **G** | **Obstructions are real:** $D_0=-72,-65,-60,-57$ are *not* realizable at pole level $n=1$; the principal-part map has rank $13$ into a $17$-dimensional slot space | **[verified]** `22_basis.gp` |
| **H** | **NEGATIVE, and it is not a data problem.**  No element of $J^{!}_{3,7}$ with pole level $n\le3$ has $(D,\beta)=(3m^2,\,5m\bmod14)$ coefficient proportional to $m^{t}\beta_{s_7}(m)$, $t\in\{-3,\dots,3\}$: rank $=\dim+1$, kernel $0$, longest matchable prefix $=\dim$ exactly (generic; no near-miss).  The test uses **only Cooper's $\beta_{s_7}$**, not the CM-trace table, so the failure cannot be blamed on round 2's $c(d)$ computation.  Same at Jacobi weights $\kappa=1,5$ | **[refuted, $n\le3$, $m\le33$]** `18_sqsearch.gp`, `21_final.gp` |
| **I** | The same with the full table: no element of $J^!_{3,7}$ ($n\le1$, $\dim\le13$) matches $c(d)$ at $(3d,\beta_d)$ over the $280$ admissible $d\le940$, with the $3\nmid D$ slots left completely free; longest matchable prefix $=\dim$ | **[refuted]** `23_full.gp` |
| **J** | **Round 2's scalar model was refuted for the wrong reason.**  Its §4.4(a) irrationality argument mixes the level-one normalisation $a(m^2)=m\beta(m)$ with the level-$28$ one $c(m^2)=\beta(m)$ that its own §5.5 asserts.  The correct reason is **A**: the Kohnen plus space of level $28$ is the *symmetric* part.  Independently, round 2's search `38_wh.gp` imposed $a(m^2)=m\beta_{s_7}(m)$, i.e. the level-one normalisation, contradicting §5.5 — doubly misdirected | **[established]** §4 |
| **K** | Direct confirmation: $\sum_d c(d)q^{d}$ (half-weights $1/2,3/2,5/2,7/2$; levels $28,56,112$; $\Delta^{r}$, $r=3..6$; constant term, polar coefficients and `DEGEN` values free) is **not** a scalar modular form.  A control with $h\in S^+_{5/2}(\Gamma_0(28))$ passes the same pipeline | **[refuted]** + **[control passes]** `05`–`08`, `07_ctrl.gp` |
| **L** | The $T,S$ transformation law of the vector-valued series was tested directly (unknown constant $C$, unknown polar vector, $6$ components $\times$ up to $10$ points $\tau$, $60$–$70$ digits) in both representations and four weights: every fit is inconsistent (residuals $O(1)$–$O(10^2)$, $|C|\ne1/\sqrt{14}$) | **[refuted]** `11`,`12`,`14`,`19`; control `13_ctrlvv.gp` |
| **M** | **Explanation of round 2's P5 "trivial middle character".**  $T_{p^2}$ on $\rho^*_L$ sends the component index $\beta\mapsto p\beta$; by **A** the datum then acquires the sign of $p\beta\equiv\pm\beta_{p^2d}$, and that sign is what round 2 absorbed into "$\tau=+1$".  The Kronecker symbol disappears because of the *antisymmetry sign of the component relabelling*, not because "the genus character squares away" | **[established]** §5 |
| **N** | **The remaining candidate home**: the $|D_0|$-rescaled lattice, discriminant form $(\mathbf Z/42,\,x^2/84)$, i.e. **Jacobi index $21$**, where the exponent is $d$ (not $3d$), every admissible $d$ occurs, and the components with $3\nmid r$ should vanish.  Not tested (needs $\approx45$ generators and $\varphi_{0,1}^{19}$) | **[conjectural]** §6 |
| **GAP** | $n^2\mid\beta_{s_7}(n)$: **still open**.  Two of round 2's three inputs are now genuinely in place on the vector-valued side (uniqueness E, integral basis F); what has moved from "done" to "open" is the **identification** — which element of the space the $s_7$ input is | **open** §6 |
| — | $s_{10}$ (level $40$) and $s_{18}$ (level $72$) | **not attempted** |

**Headline.**  Round 2 declared the level-$28$ proof blocked on one missing ingredient, "an integral
basis of $M^!_{5/2}(\rho_L)$ indexed by $(\beta,d)$ together with the tower relation".  That
ingredient is now supplied: the space is $J^{!}_{3,7}$, weakly holomorphic Jacobi forms of weight
$3$ and index $7$; it is spanned explicitly by
$\Delta^{-n}\varphi_{-1,2}\varphi_{0,1}^{5-b}\varphi_{-2,1}^{b}M_{4+2b+12n}$; its holomorphic part
vanishes, so the principal part determines the form (uniqueness for free, replacing round 2's
Atkin–Lehner argument); and the resulting basis $\{g_{D_0}\}$ is integral.  What broke instead is
the *identification*: no element of that space has $\beta_{s_7}(m)$ — or $m^{\pm1,\pm2,\pm3}
\beta_{s_7}(m)$ — as its coefficient at the Heegner slot $(D,\beta)=(3m^2,5m)$.  Since that test
uses only Cooper's $\beta$, the discrepancy is in the *dictionary*, not in round 2's CM-trace
computation.  The natural repair — the $|D_0|$-rescaled lattice, Jacobi index $21$ — is identified
but untested.

---

## 1. What the object is: antisymmetry, and the two things it proves

Round 2 computed $c(d)=\sqrt{-3}\sum_Q\chi_{-3}(Q)\widehat f(\alpha_Q)/\omega_Q$ over the
$\Gamma_0(7)$-Heegner classes of discriminant $-3d$ with $B\equiv\beta\pmod{14}$, choosing **one**
$\beta$ per $d$ (for squares $\beta=5m$, otherwise the smaller of the two).  The first thing to do
is to compute *both* classes.

`09_beta.gp` does this for every admissible $d\le120$ at 45–50 digits (re-using round 2's `lib.gp`,
`heeg.gp`, `maass2.gp` read-only; it must be run from `../round2` because those use relative
`read`s).  The result is uniform and unambiguous:

$$\boxed{\;c(-\beta,d)=-\,c(\beta,d)\;}$$

e.g. $d=1$: $c(5,\cdot)=+1$, $c(9,\cdot)=-1$; $d=4$: $\pm8$; $d=25$: $\pm125$; $d=120$: $\pm12000$;
imaginary parts $<10^{-53}$ throughout.

**Proof.**  $W_7$ acts on $\Gamma_0(7)$-Heegner forms of discriminant $D$ by
$[A,B,C]\mapsto[7C,-B,A/7]$, hence on $\beta$-classes by $\beta\mapsto-\beta$; the genus character
$\chi_{D_0}$ is invariant under it; and $\widehat f\mid W_7=-\widehat f$ (round 2 R2/R3,
**[proved]**).  Summing over a class therefore changes sign. $\square$

**1.1 The $\beta\in\{0,7\}$ components vanish, so $c(d)=0$ whenever $7\mid d$.**
If $7\mid d$ then $\beta^2\equiv-3d\equiv0\pmod 7$, so $7\mid\beta$, i.e. $\beta\in\{0,7\}$ mod
$14$; those are exactly the classes with $\beta\equiv-\beta$, where antisymmetry forces $c=0$.
Conversely $\beta\in\{0,7\}\Rightarrow 7\mid 3d\Rightarrow7\mid d$.  Equivalently: the admissible
$d$ with $7\mid d$ are exactly $d\equiv0,21\pmod{28}$.

This **proves** round 2's P6 bullet ("$G(d)=0$ exactly when $7\mid d$", there verified on 690
entries) and it **resolves the 25 `DEGEN` entries** of `73_cd_s7_2500.txt`
($d=49,196,245,392,441,588,637,784,833,980,1029,1176,1225,1372,1421,1568,1617,1764,1813,1960,2009,2156,2205,2352,2401$,
all of the shape $49e$): the numerical trace there is a spurious purely imaginary residue coming
from the imprimitive classes of discriminant $-3\cdot49e$ (at $d=49$ the run returns
$-223.331\ldots i$); the true component value is $0$.  Round 2's table can be completed with zeros.

**1.2 The plus space is the wrong space, for a reason.**  Kohnen's plus space
$M^{+}_{5/2}(\Gamma_0(28))$ is the *symmetric* part of the vector-valued space: its unique cusp
form $h=-q+3q^4-7q^8+5q^9+5q^{16}-14q^{21}-11q^{25}+7q^{28}+\cdots$ is supported also on the
self-conjugate classes ($a(21)=-14$, $a(28)=7$) that antisymmetry kills.  So no scalar plus form
can be the $s_7$ input.

---

## 2. Which representation, and hence which space

Milgram's formula $\sum_{\beta\in\mathcal D}e(Q(\beta))=\sqrt{|\mathcal D|}\,e(\mathrm{sig}/8)$
applied to $\mathcal D=(\mathbf Z/14,\,x^2/28)$ gives (`20_diag.gp`)
$$\sum_{\beta=0}^{13}e(\beta^2/28)=\sqrt{14}\,e(1/8),\qquad \mathrm{sig}(\mathcal D)\equiv1\pmod 8 .$$
For a Weil representation of a discriminant form of signature $s$, nonzero weight-$k$ forms are
symmetric when $2k\equiv s\pmod4$ and antisymmetric when $2k\equiv s+2\pmod4$.  With $2k=5$:

* $\rho_{\mathcal D}$ ($s\equiv1$): weight-$5/2$ forms are **symmetric**; exponents $\equiv+Q(\beta)$, i.e. $d/28$ with $d\equiv\beta^2\ (28)$ — the Kohnen plus space of level $28$;
* $\rho^{*}_{\mathcal D}$ ($s\equiv-1$): weight-$5/2$ forms are **antisymmetric**; exponents $\equiv-Q(\beta)$, i.e. $|D|/28$ with $D\equiv-\beta^2\ (28)$.

By §1 the $s_7$ input is antisymmetric, hence lives in the second space, with exponents $|D|/28$,
$D=-3d$.  Both descriptions are compatible with the observed support, because $-3\equiv5^2\ (28)$
makes "$d$ a square mod $28$" and "$-3d$ a square mod $28$" the same condition; only the *exponent*
differs, and the signature fixes it to be $3d$.

Antisymmetric weight-$5/2$ forms for $\rho^*_{(\mathbf Z/2m,\,x^2/4m)}$ are exactly the theta
components of **Jacobi forms of weight $3$ and index $m$** (odd Jacobi weight $\Leftrightarrow$
$c(n,-r)=-c(n,r)$).  With $m=N=7$:
$$\boxed{\;M^{!,-}_{5/2}(\rho^*_L)\;=\;J^{!}_{3,7}\;}$$

---

## 3. The space, explicitly: uniqueness, integral basis, obstructions

`15_jaclib.gp` builds the three Eichler–Zagier generators as exact $(q,\zeta)$-expansions,
$$\varphi_{-2,1}=(\zeta-2+\zeta^{-1})\prod_{n\ge1}\frac{(1-q^n\zeta)^2(1-q^n\zeta^{-1})^2}{(1-q^n)^4},\qquad
\varphi_{-1,2}=(\zeta-\zeta^{-1})\prod_{n\ge1}\frac{(1-q^n\zeta^{2})(1-q^n\zeta^{-2})}{(1-q^n)^2},$$
$$\varphi_{0,1}=4\sum_{i=2,3,4}\bigl(\theta_i(\tau,z)/\theta_i(\tau,0)\bigr)^{2},$$
represented as $z^{-R}P(z)$ with $P\in\mathbf Q[z]$ having $q^{1/2}$-series coefficients.
`16_jactest.gp` checks them against the standard expansions
($\varphi_{-2,1}=(\zeta-2+\zeta^{-1})+q(-2\zeta^{2}+8\zeta-12+8\zeta^{-1}-2\zeta^{-2})+\cdots$,
$\varphi_{0,1}=(\zeta+10+\zeta^{-1})+q(10\zeta^{2}-64\zeta+108-\cdots)+\cdots$,
$\varphi_{-1,2}=(\zeta-\zeta^{-1})+q(\zeta^{3}-3\zeta+3\zeta^{-1}-\zeta^{-3})+\cdots$) and verifies
the index law ($c(n,r)$ depends only on $4mn-r^2$ and $r\bmod2m$): $0$ violations, likewise for
every basis element of $J^{\mathrm{weak}}_{3,7}$ (`21_final.gp`).

By Eichler–Zagier the odd-weight weak Jacobi forms are free of rank one on $\varphi_{-1,2}$ over
$M_*[\varphi_{0,1},\varphi_{-2,1}]$, so
$$J^{\mathrm{weak}}_{3+12n,7}=\varphi_{-1,2}\cdot\bigoplus_{b=0}^{5}M_{4+2b+12n}(\mathrm{SL}_2\mathbf Z)\,\varphi_{0,1}^{5-b}\varphi_{-2,1}^{b},
\qquad J^{!}_{3,7}=\bigcup_{n\ge0}\Delta^{-n}J^{\mathrm{weak}}_{3+12n,7},$$
of dimensions $7,13,19,25,31$ for $n=0,1,2,3,4$.

**3.1 Uniqueness [verified, $n\le2$].**  The map "form $\mapsto$ its coefficients at the slots
$D<0$" has rank $=\dim$ for $n=0,1,2$, so its kernel — the holomorphic Jacobi forms — is $0$:
$$J_{3,7}=0 .$$
This is the level-$28$ replacement for $\dim S_4(\mathrm{SL}_2\mathbf Z)=0$, and it is cleaner than
round 2's P4: no Atkin–Lehner bookkeeping and no auxiliary condition "$a(n)=0$ for $7\mid n$",
because the antisymmetry has already removed the $\beta\in\{0,7\}$ components.  It agrees with
Skoruppa–Zagier ($J^{\mathrm{cusp}}_{3,7}\cong S^{\mathrm{new}}_{4}(\Gamma_0(7))^{W_7=(-1)^3}=0$)
and with the absence of odd-weight Jacobi Eisenstein series.

**3.2 The integral basis [constructed].**  `22_basis.gp` solves, for each admissible $D_0<0$, for
the element of $\Delta^{-1}J^{\mathrm{weak}}_{15,7}$ whose principal part is the single slot $D_0$.
Thirteen of the seventeen slots $D_0\ge-72$ are realizable, and each solution has **integer**
coefficients — the denominator over the first $24$ positive slots is $1$ in every case.  The
smallest few, as $c(D)$ at $D=3,12,19,20,24,27,31,40$:

| $D_0$ | $\beta$ | $c(3),c(12),c(19),c(20),c(24),c(27),c(31),c(40)$ |
|---|---|---|
| $-1$ | $1$ | $1,\;12,\;31,\;20,\;32,\;30,\;61,\;152$ |
| $-4$ | $2$ | $-2,\;-13,\;26,\;-78,\;154,\;156,\;-304,\;-338$ |
| $-8$ | $6$ | $2,\;-52,\;214,\;-114,\;-411,\;372,\;1072,\;-4142$ |
| $-9$ | $3$ | $1,\;-52,\;-225,\;324,\;272,\;1054,\;828,\;-4824$ |
| $-16$ | $4$ | $2,\;142,\;-938,\;-1460,\;-1564,\;6324,\;5872,\;48068$ |
| $-25$ | $5$ | $-5,\;260,\;2405,\;5772,\;-20720,\;34410,\;-106555,\;296920$ |
| $-29$ | $1$ | $5,\;896,\;9978,\;5952,\;33600,\;40701,\;194048,\;1727936$ |
| $-32$ | $2$ | $-6,\;-546,\;6942,\;-16992,\;69192,\;125892,\;-434384,\;-1512576$ |
| $-36$ | $6$ | $6,\;-1885,\;27090,\;-16794,\;-105410,\;138508,\;732816,\;-8031510$ |
| $-37$ | $3$ | $4,\;-1664,\;-24831,\;43456,\;68544,\;368373,\;486912,\;-7933632$ |
| $-44$ | $4$ | $4,\;3352,\;-61716,\;-110955,\;-194088,\;1125672,\;1630688,\;31988160$ |
| $-53$ | $5$ | $-15,\;4224,\;101652,\;272448,\;-1510848,\;3407148,\;-15394304,\;92900544$ |
| $-64$ | $6$ | $18,\;-24466,\;778342,\;-534744,\;-4915408,\;8440308,\;62099824,\;-1339959496$ |

Integrality is not an accident of the solve: the generators
$\varphi_{-2,1},\varphi_{0,1},\varphi_{-1,2},E_4,E_6,\Delta$ all have integer coefficients, and
$J^!_{3,7}$ is a module over $\mathbf Z[j]$ (multiplication by $j(\tau)$ shifts $D$ by $-28$ inside
a fixed $\beta$-class), which propagates integrality along each class from a seed — exactly the
mechanism the level-one proof uses on Borcherds–Zagier's $\{g_M\}$, now available at level $28$.

**3.3 Obstructions [verified].**  $D_0=-72,-65,-60,-57$ are *not* realizable at pole level $n=1$:
the principal-part map has rank $13$ into a $17$-dimensional slot space.  Unlike the level-one plus
space, not every principal part occurs at a given pole level, so the tower elements
$g_{3p^{2r}}$ must be checked realizable before the round-2 induction can be run.

---

## 4. Round 2's scalar model: refuted, but not by round 2's argument

Round 2 §4.4 gave two reasons to reject the scalar plus space of level $28$.

* **(a) Irrationality** ("matching the square indices would force $a(d)=\sqrt d\,c(d)$, irrational
  at $d=8$").  This mixes normalisations: the $\sqrt d$ comes from the *level-one* dictionary
  $a(m^2)=m\beta(m)$ (round 2 R10), while round 2 §5.5 states that at level $28$ the normalisation
  is $c(m^2)=\beta(m)$ with **no** extra factor $m$.  Under the latter there is no $\sqrt d$ and no
  irrationality.  The argument as written does not stand.
* **(b) Support law.**  True but not decisive: $-3\equiv5^2\ (28)$, so the observed support
  "$-3d\equiv\square\ (28)$" *is* the level-$28$ Kohnen plus condition "$d\equiv\square\ (28)$".

The actual obstruction is §1: the input is **antisymmetric** and the plus space is the symmetric
part.  For the record we also checked the scalar model directly and negatively
(`05_find.gp`, `06_sweep.gp`, `08_sweep2.gp`): with $F=\sum_dc(d)q^d$ and the polar coefficient
$c(-3)$, the constant term and the `DEGEN` values free, $F\theta^{8-j}\Delta^{r}$ lies in no
$M_{4+12r}(\Gamma_0(L))$ for $j=1,3,5,7$, $r=3,4,5,6$, $L=28,56,112$.  A control (`07_ctrl.gp`)
with the genuine Kohnen form $h$ passes the same pipeline at every $r$, so the negative is not a
bookkeeping artefact.  (Trap worth recording: PARI's `eta(q)` is $\prod(1-q^n)$, **not**
$q^{1/24}\prod$; using it as $\Delta$ silently produces the wrong answer, and did, until the
control caught it.)

Round 2's own search `38_wh.gp` looked for $f$ with $a(m^2)=m\,\beta_{s_7}(m)$ — the level-one
normalisation, contradicting its own §5.5 — inside the scalar plus space.  It was searching the
wrong space for the wrong vector, and its negative should not be quoted as evidence about the plus
space.

---

## 5. What the antisymmetry says about round 2's Hecke fit

On $M^{!}_{5/2}(\rho^*_L)$, for $p\nmid14$,
$$(F\mid T_{p^2})(\beta,D)=c(p\beta,\,p^2D)+p\Bigl(\tfrac{D}{p}\Bigr)c(\beta,D)+p^{3}c(\beta/p,\,D/p^{2}),$$
so the leading term reads the coefficient in the component $p\beta$, not $\beta_{p^2d}$.  These
differ by a sign, because $p\beta\equiv\pm\beta_{p^2d}$ and, by §1,
$c(-\gamma,\cdot)=-c(\gamma,\cdot)$.  Round 2's P5 fitted
$c(p^2d)+\tau\,p\,c(d)+p^3c(d/p^2)=\lambda c(d)$ with $\tau=+1$ and read the outcome as "the genus
character has squared away, $\chi_{-3}(p)^2=1$".  The correct reading: the Kronecker symbol has
been absorbed into the **antisymmetry sign of the component relabelling $\beta\mapsto p\beta$**.
The two readings agree numerically but differ structurally, and only the second survives to a
proof, because it is a statement about $\rho^*_L$ rather than about a scalar space.

---

## 6. Where the proof stands, and the lemma that is missing

| input | level 1 | level 28, after this round | status |
|---|---|---|---|
| principal part determines the form | $\dim S_4(\mathrm{SL}_2\mathbf Z)=0$ | $J_{3,7}=0$ (E) | **done**, cleaner than round 2's P4 |
| integral basis $\{g_M\}$ | Borcherds–Zagier | $\{g_{D_0}\}\subset J^!_{3,7}$, explicit, integral (F) | **done** |
| eigen-identity (T1) + tower (T2) | proved | not yet written for $\rho^*_L$ | open, now formulable (§5) |
| **identification**: which element is the $s_7$ input | trivial | **fails** with the natural dictionary (H, I) | **the new gap** |

The gap has moved.  Round 2 believed the identification settled and the space missing; we have the
space, and the identification fails.  With $(D,\beta)=(3m^2,\,5m\bmod14)$ — the only slot allowed
by $D\equiv-\beta^2\ (28)$ — no element of $J^!_{3,7}$ of pole level $\le3$ has coefficient
$m^{t}\beta_{s_7}(m)$ there for any $t\in\{-3,\dots,3\}$, over $m\le33$, $7\nmid m$.  The test is a
rank computation on an over-determined system, and the longest matchable prefix equals the
dimension exactly, i.e. the system behaves generically — there is no near-miss.

**The most likely repair (N).**  In the twisted (genus-character) Shimura–Borcherds setting the
natural lattice is the $|D_0|$-rescaled one: discriminant form
$(\mathbf Z/2N|D_0|,\,x^2/4N|D_0|)=(\mathbf Z/42,\,x^2/84)$, i.e. **Jacobi index $21$**.  There the
exponent of the trace generating function is $|D_0d|/(4N|D_0|)=d/28$ — the GKZ convention, in which
the generating variable is $d$, not $3d$ — every admissible $d$ occurs, and the components with
$3\nmid r$ should vanish identically (the "$3\mid D$" pattern the index-$7$ model cannot
accommodate: an index-$7$ form has coefficients at all admissible $D$, whereas the twisted traces
populate only $D\in3\mathbf Z$).  Testing this needs
$J^{\mathrm{weak}}_{3+12n,21}=\varphi_{-1,2}\bigoplus_{b=0}^{19}M_{4+2b+12n}\varphi_{0,1}^{19-b}
\varphi_{-2,1}^{b}$ ($\approx45$ generators at $n=0$) and $q$-precision $\sim84n$; `15_jaclib.gp`
supports it verbatim, only the cost is larger.

A second, cheaper possibility to exclude first: that the Shimura-type dictionary between
$\beta_{s_7}$ and the vector-valued coefficients carries a factor which is *not* a power of $m$
(a class number, a $\Gamma_0(7)$-index factor, or an Atkin–Lehner factor of the kind round 2 found
at $2$ for $s_{10}$).  The present test would miss such a factor.

$s_{10}$ (level $40$, $D_0=-4$, index $10$) and $s_{18}$ (level $72$, $D_0=-36$, index $18$) were
not attempted; the Jacobi weight is again $3$, with $J_{3,10}$, $J_{3,18}$ and the $|D_0|$-rescalings
to consider.

---

## 7. Honest ledger

* **Achieved.**  The representation-theoretic home of the $s_7$ input (A, C, D); a proof that
  $c(d)=0$ for $7\mid d$, with the resolution of round 2's `DEGEN` entries (B); unconditional
  uniqueness at level $28$ on the vector-valued side (E); an explicit, verified-integral basis of
  the space round 2 called "the whole gap" (F), with its obstructions (G); corrections to two of
  round 2's arguments (J, M).
* **Not achieved.**  The target theorem $n^2\mid\beta_{s_7}(n)$.  Not even the identification of the
  $s_7$ input inside the space: the natural dictionary is refuted (H, I).  Consequently (T1), (T2)
  and the induction were not run, and $s_{10}$, $s_{18}$ were not started.
* **Regression.**  One thing round 2 recorded as settled — the reading of $c(m^2)=\beta_{s_7}(m)$ as
  a vector-valued Fourier coefficient — is now in doubt.  Round 2's numerical *identity*
  $\beta_{s_7}(m)=\sqrt{-3}\sum_Q\chi_{-3}(Q)\widehat f(\alpha_Q)/\omega_Q$ (R1, 70 digits) is
  untouched; what is in doubt is its reading as "the $(3m^2,5m)$ coefficient of a weight-$5/2$ form
  for $\rho^*_L$".
* **Ranges.**  Antisymmetry: every admissible $d\le120$, 45 digits.  Jacobi searches: pole level
  $n\le3$ (dim $\le25$), $m\le33$ ($30$ equations), $d\le940$ ($280$ equations).  Integral basis:
  $n=1$, first $24$ positive slots.  Scalar searches: $400$–$600$ $q$-coefficients.
* **Reproducibility.**  Every claim is produced by a script here with its `.log`; the controls
  (`07_ctrl.gp`, `13_ctrlvv.gp`, `16_jactest.gp`, `20_diag.gp`) exist precisely so that the
  negatives can be trusted.

---

## 8. Scripts

| file | what it does |
|---|---|
| `01_setup.gp` | PARI at level $28$: $\dim M_{5/2}=8$, $\dim S_{5/2}=4$, $\dim S^+_{5/2}=1$ with generator $-q+3q^4-7q^8+5q^9-14q^{21}+\cdots$; cusps; $\dim M_{29/2}=56$ |
| `03_dims.gp`, `04_dim40.gp` | $\dim M_{53/2}(\Gamma_0(28))=104$ (8.6 min); $M_{77/2}$ overflows 2 GB — half-integral weight is the expensive route.  $\dim M_{40}(\Gamma_0(28))=159$ in 3.7 s — the integral-weight route via $\theta^3$ |
| `05_find.gp`, `06_sweep.gp`, `08_sweep2.gp` | the scalar model: is $\sum c(d)q^d$ a weight-$j/2$ form on $\Gamma_0(L)$?  Negative for $j=1,3,5,7$, $L=28,56,112$, $\Delta$-powers $3..6$ |
| `07_ctrl.gp` | **control**: $h\theta^3\Delta^r\in M_{4+12r}(\Gamma_0(28))$ for $r=0..3$ — validates the multiplier bookkeeping (and caught the `eta` trap) |
| `09_beta.gp` | **the antisymmetry**: both $\beta$-classes at every admissible $d\le120$ (run from `../round2`) |
| `10_twist.gp` | the odd-character twist $\widetilde c(d)=(\beta_d/7)c(d)$ against even characters mod $28$ (superseded) |
| `11_vv.gp`, `12_vv2.gp`, `14_vv3.gp`, `19_vv5.gp` | the numerical $S$-transformation law, $60$–$70$ digits, both representations, weights $1/2,3/2,5/2,7/2$, coefficient powers $d^{-1/2},d^0,d^{1/2},d$, polar vector free |
| `13_ctrlvv.gp` | control of that machinery on the Kohnen form |
| `15_jaclib.gp` | **the weak Jacobi library**: $\varphi_{-2,1},\varphi_{0,1},\varphi_{-1,2},E_4,E_6,\Delta$ as exact $(q,\zeta)$-objects; products, coefficients |
| `16_jactest.gp` | validation of the generators (standard expansions; index consistency, $0$ violations) |
| `17_search.gp` | identification inside $\Delta^{-n}J^{\mathrm{weak}}_{3+12n,7}$ against the full $c(d)$ table with $3\nmid D$ forced to $0$ — negative |
| `18_sqsearch.gp` | identification using **only** $\beta_{s_7}(m)$ at $(3m^2,5m)$, normalisations $m^t$, $t=-3..3$ — negative |
| `20_diag.gp` | Milgram signature; $J_{3,7}=0$ and injectivity of the principal-part map, $n=0,1,2$; antisymmetry of the basis |
| `21_final.gp` | the same search at $q$-precision $140$, Jacobi weights $\kappa=1,3,5$, pole levels $n=0..4$ |
| `22_basis.gp` | **the integral basis** $\{g_{D_0}\}$ of $J^!_{3,7}$ at pole level $1$: realizability and integrality |
| `23_full.gp` | identification against the full table with the $3\nmid D$ slots free; longest matchable prefix |

Run with `gp -q <file>` from this directory, except `09_beta.gp`, which must be run from
`../round2` and writes its output back here.
