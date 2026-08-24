# The deep four-term scan: mixed-exponent classes, the decayer question, and Catalan

*Claude (Opus 5), 2026-08-23/24. Scripts, logs and data: `lattice/four_term_deep/`.
Sequel to `FOUR_TERM_SCAN.md` (Theorems F1–F6, the first census), `K3_ROW_PERIOD.md`
(the period of the one elliptic-$K3$ row), `CATALAN_OBSTRUCTION.md` §3 (why a new
realisation of the Catalan class would matter) and `MULTI_PRIME_LATTICE.md` (what a
decayer with new $\kappa_p$ arithmetic would buy).*

*Tags: **[proved]** = exact argument; **[verified]** = exact finite computation over a
stated range; **[measured]** = numerical at stated precision.*

*Everything below is either proved, or verified by exact / high-precision computation over a
stated range. The scans are $a$-block ordered, so each reported box is a complete scan of the
stated prefix; nothing here depends on a scan that was cut off mid-block.*

---

## 0. The two questions and the verdict

The brief was to look, in the four-term (five-singular-point) world, for either

* **(a)** a row whose Apéry limit involves Catalan's $G=L(2,\chi_{-4})$ — any algebraic
  combination with $\zeta(2)$, $\pi^2$, other $L$-values; or
* **(b)** a genuine five-point row with **real** $|\lambda_2|<1$ (a *decayer*) whose
  period is a classical constant,

either of which would feed the two-row machinery of `MULTI_PRIME_LATTICE.md` with
fresh $\kappa_p$ arithmetic.

**Verdict (a): POSITIVE.** In the **mixed-exponent classes** — the ones the first census
never scanned, with exponent differences $(0,\tfrac12,\tfrac12)$ at the three finite points
— there are genuine five-singular-point integral four-term rows whose Apéry limits are
$$\tfrac14G,\qquad \tfrac12G-\tfrac3{16}\zeta(2),\qquad \tfrac38\zeta(2)-\tfrac12G,$$
and three more with $L(2,\chi_{-3})$ and $\zeta(2)$; all verified to $139$–$141$ digits,
integral to $n=200$, $k=2$, no apparent singularity, and **provably not three-term**. Their
singular points are irrational — conjugate over $\mathbf Q(\sqrt2)$ for the Catalan rows and
over $\mathbf Q(\sqrt3)$ for the $\chi_{-3}$ ones — which is exactly the evasion **(E1)** of
`CATALAN_OBSTRUCTION.md` §3. The best Catalan score is $-2.6931$, $+0.693$ nats over
Zagier $\mathbf E$. See **§6.4**.

**And the whole period matrix is classical, not just the Apéry limit (§6.5).** The fold
constants at *all three* singular points of these rows — real and imaginary parts — are
identified at $125$–$129$ digits: for the $\sqrt2$ family they all lie in
$\mathbf QG+\mathbf Q\zeta(2)$, for the $\sqrt3$ family in
$\mathbf QL(2,\chi_{-3})+\mathbf Q(\sqrt3)\zeta(2)$. One of them,
$(8,4,32,8,64)_{r=4}$, has *no* archimedean Apéry limit (dominant pair $2\pm2\sqrt3\,i$) and
still has three classical fold periods, so this is a property of the local system, not of the
recurrence's convergence.

**Their gauge partners give the lemniscate period.** The same class in the other
half-integral gauge, $(+\tfrac12,0;1,1,1)$, contains the rows with the *same* $(a,d,f,C)$ and a
different accessory $c$, and all three have
$\xi=\Gamma(\tfrac14)^4/(64\pi)$ — verified to $195$–$241$ digits. The two gauges of one
five-point local system over $\mathbf Q(\sqrt2)$ give the $L$-value side and the
algebraic-period side of the Gaussian world. See **§6.6**.

Sweeping both classes to completion over $|r|\le60$, $|a|\le200$, $|c|\le200$, $|d|\le2000$,
$|f|\le400$, $|C|\le600$ gives **exactly nine primitive genuine five-point rows in each**,
in twist-partner bijection; every one with an archimedean limit in the $-\tfrac12$ gauge is
identified, and three in the $+\tfrac12$ gauge share one unidentified constant $\xi^*$ which
is *also* the Apéry limit of a five-term six-point row — the same number to $319$ digits from
two recurrences of different order (§6.7).

**But the arithmetic verdict on them is negative (§6.8).** Their $\sigma_2$ is genuinely new
in *form* — half-integral ($\tfrac72$ for two of them), decoupled from $v_2(g)$, and row 1
even has $\kappa_2=-\tfrac12$, an integral row with $v_2(a_n)=\tfrac n2+O(1)$, none of which
occurs anywhere in the `MULTI_PRIME_LATTICE.md` census — but the *content* is old: the slope
set is $\{2\}$ (or empty), the $2$-adic limits are $\tfrac14\zeta_2(2)$, $\pm\tfrac12\zeta_2(2)$,
$\tfrac38$ and $\tfrac34L_2(2,\chi_{12})$ — the ledger's own two constants — and the engine
yield $\sum_p\sigma_p\log p-\log|\lambda_2|$ is $0.35,0.69,1.73,0.07,2.01,0$ against
$2.079$ for Zagier $\mathbf E$. **Every one is dominated by the engine it would replace: §3
of `CATALAN_OBSTRUCTION.md` has been re-decorated, not evaded.**

By contrast the *disguised* route fails, sharply: the only integral three-term row with a
Catalan period is Zagier $\mathbf E$ ($\xi=\tfrac12G$), and its four-term cusp-move images
*lose* the period — their canonical Apéry limits are $\tfrac12GA(\nu/32)-B(\nu/32)$, interior
values of the Eichler integral, unidentified (Proposition D4, §6.1). And no *equal*-exponent
four-term row, old census or new, has a Catalan period (§6.3).

**Verdict (b): negative, and now for a structural reason.** Theorem D1: $|\lambda_2|<1$
forces an *irreducible* characteristic cubic, so no class with unequal exponents can hold a
decayer at all; Theorem D2: a **real** $\lambda_2$ additionally forces the cubic to be
totally real ($\Delta>0$). Empirically, the extended window (all $27$ equal-$\rho$ classes,
$|c|\le250$, $|f|\le500$, complete in $(d,C)$ for every $a\le180$, about $17\times$ the
$(c,f)$ volume of box W) contains **exactly three** rows with $|\lambda_2|<1$: two disguised
$I_0^*$ Padé rows and the single known genuine one, whose two small roots are a *complex*
pair ($\Delta=-81\,920\,000$). Over all $20$ mixed-class candidates the smallest
$|\lambda_2|$ is **exactly $1$**, never below — Corollary D1.1 in action. See §3, §4.2.

**The five-term side quest produced one more.** In the **six-singular-point** world, the
all-semistable class contains
$$(n+1)^2u_{n+1}=(n^2{+}n)u_n+(14n^2{+}2)u_{n-1}+(20n^2{-}20n{+}8)u_{n-2}+8(n{-}1)^2u_{n-3},$$
$u_n=1,0,4,8,40,144,616,2544,\dots$, integral, primitive, six distinct singular points, all
semistable, $k=2$, score $-2.6931$, whose Apéry limit is
$$\xi=\frac{\sqrt2}3L(g,2)=\tfrac12\xi_{K3}(1),$$
the critical value of the weight-$3$ CM newform **32.3.d.a** — i.e. the period of the
elliptic-$K3$ four-term row, which in the four-term world exists only as a *fold constant*
because that row has no archimedean limit. **[verified to 110 digits]** See §7.3.

---

## 1. What is new here

1. **Theorem D1 (a structural obstruction to decayers).** A four-term row with
   $|\lambda_2|<1$ *strictly* has an **irreducible** characteristic cubic, hence
   (Corollary F2.1) all three finite exponent differences are equal and the three
   finite fibres are Galois-conjugate of one Kodaira type. **[proved]**, §3.
2. **Theorem D2.** If in addition $\lambda_2$ is **real**, the cubic is irreducible and
   **totally real** ($\Delta>0$). In the whole first census there is no such row: the
   only row with $|\lambda_2|<1$ has $\Delta=-81\,920\,000<0$, i.e. its two small roots
   are a complex pair. **[proved + verified]**, §3.
3. **The mixed-exponent classes**, absent from the first census: $\rho_1\ne\rho_2=\rho_3$
   with one or two of the exponent differences $\tfrac12$ and the others $0$. Their
   normalisation forms are derived from the §2 dictionary (Theorem D3, §2) and verified
   against Theorem F1. Twenty such classes are scanned. They are the quadratic twists,
   ramified at two of the three finite points, of the equal-$\rho$ classes; they carry
   integral rows the equal-$\rho$ scan cannot see. **[proved + verified]**
   **They are where the Catalan rows live** — §6.4, the main result.
4. **Corollary D1.1.** No mixed-exponent class contains a decayer at all: a rational
   characteristic root is forced, and D1 forbids it. So the mixed classes bear only on
   question (a). **[proved]**
5. **The complex-fold periods of the first census** (the nine rows with no archimedean
   Apéry limit, including the elliptic $K3$) are computed for the first time — 27 fold
   constants to $\ge154$ digits — and identified. **[verified]**, §5.
6. **A five-term (six-singular-point) probe**: Theorem D5 (the operator dictionary one term
   further), a scanner validated against exact brute force, and
   **an integral six-point row whose Apéry limit is the elliptic-$K3$ $L$-value
   $\tfrac{\sqrt2}3L(g,2)$**, the first weight-$3$ cusp-form critical value to occur as an
   ordinary Apéry limit anywhere in the project, §7.3.
7. **Six genuine five-point rows with Catalan / $L$-value Apéry limits**, plus their three
   gauge partners with $\xi=\Gamma(\tfrac14)^4/(64\pi)$, all in the mixed class
   $(\mp\tfrac12,0;1,\cdot,\cdot)$: $\tfrac14G$, $\tfrac12G-\tfrac3{16}\zeta(2)$,
   $\tfrac38\zeta(2)-\tfrac12G$, and three with $L(2,\chi_{-3})$; verified to $139$–$141$
   digits and shown to satisfy **no** three-term recurrence. Their singular points are
   irrational — the evasion (E1) of `CATALAN_OBSTRUCTION.md` §3 realised explicitly. §6.4.
8. **Proposition D4**: the generic cusp move destroys the Apéry limit. Applied to Zagier
   $\mathbf E$ this is the precise reason why the disguised four-term rows built on the
   Catalan local system do not have Catalan periods. §6.1.

---

## 2. Theorem D3: the mixed-exponent normalisation forms

Keep the notation of `FOUR_TERM_SCAN.md` §1: $\chi(\lambda)=\lambda^3-a\lambda^2+d\lambda-g$
has roots $\lambda_i=1/t_i$, and $T=\mathcal S-t\mathcal R'$ has
$\rho_i=-T(t_i)/(t_i\mathcal R'(t_i))$.

**Lemma D3.0.** *Let $W(\lambda)=\lambda^3+\beta_1\lambda^2+\beta_2\lambda+\beta_3$ be the
reversal of $T$, i.e. $T(t)=1+\beta_1t+\beta_2t^2+\beta_3t^3$. Then*
$$\boxed{\ W(\lambda)=\chi(\lambda)\Bigl[1+\sum_{i}\frac{\rho_i\lambda_i}{\lambda-\lambda_i}\Bigr]
   =\chi(\lambda)+\sum_i\rho_i\lambda_i\prod_{j\ne i}(\lambda-\lambda_j)\ }$$
*and $b=-\beta_1$, $d+e=\beta_2$, $2g+h=-\beta_3$.*

*Proof.* $t_iR'(t_i)$ where $\mathcal R=\prod(1-\lambda_it)$ gives
$t_i\mathcal R'(t_i)=-(\lambda_i-\lambda_j)(\lambda_i-\lambda_k)/\lambda_i^2$, so the
condition $\rho_i=-T(t_i)/(t_i\mathcal R'(t_i))$ reads $W(\lambda_i)=\rho_i\lambda_i\chi'(\lambda_i)$.
$W$ is the unique **monic** cubic with those three values; Lagrange interpolation with
$L_i(\lambda)=\chi(\lambda)/((\lambda-\lambda_i)\chi'(\lambda_i))$ gives the display, and
the leading coefficient is $1$ because the correction term has degree $\le2$. $\square$

Setting all $\rho_i=\rho$ recovers $W=\rho\lambda\chi'+(1-3\rho)\chi$, i.e. Theorem F2.

**Theorem D3 (mixed classes).** *Let $\chi=(\lambda-r)(\lambda^2-\sigma\lambda+\pi)$ with
$r,\sigma,\pi\in\mathbf Z$, and put $\rho_r$ on the rational root and $\rho_p$ on both
roots of the quadratic factor (Corollary F2.1 forces $\rho$ to be constant on Galois
orbits, so this is the general unequal-exponent shape with a quadratic orbit). Then*
$$a=r+\sigma,\qquad d=\sigma r+\pi,\qquad g=r\pi,$$
$$\boxed{\ b=(1-\rho_r)r+(1-\rho_p)\sigma,\qquad
   e=-\rho_p(2\pi+r\sigma)-\rho_r r\sigma,\qquad
   h=-(1+2\rho_p+\rho_r)\,g\ }$$
*and, writing $R(n)=C(Mn-j_1)(Mn-j_2)$ as in §2.2 of `FOUR_TERM_SCAN.md`,*
$$g=CM^2,\qquad j=Cj_1j_2,\qquad \frac{j_1+j_2}{M}=s_1+s_2=1+2\rho_p+\rho_r .$$
*Fuchs is automatic: $\sum\rho_i=2\rho_p+\rho_r=s_1+s_2-1$.*

*Proof.* Expand Lemma D3.0 with $\rho_1=\rho_r$ at $\lambda_1=r$ and $\rho_2=\rho_3=\rho_p$:
$$W=\chi+\rho_p(\lambda-r)(\sigma\lambda-2\pi)+\rho_r r(\lambda^2-\sigma\lambda+\pi),$$
using $\lambda_2\prod_{j\ne2}(\lambda-\lambda_j)+\lambda_3\prod_{j\ne3}(\lambda-\lambda_j)
=(\lambda-r)(\sigma\lambda-2\pi)$, and read off the coefficients. $\square$

`lattice/four_term_deep/01_mixcheck.py` verifies Theorem D3 numerically against the
Theorem F1 dictionary for all $15$ shapes tried (six random parameter points each):
the $\rho_i$ at the three roots and $\delta_\infty$ come out exactly as prescribed,
$0$ failures. **[verified]**

### 2.1 The mixed classes scanned

A class is now $(\rho_p,\rho_r;M,j_1,j_2)$. Twenty are Kodaira-admissible, non-Casoratian
(no positive integer root of $R$) and genuinely mixed:

| $(\rho_p,\rho_r)$ | fibre shape at $t_1,t_2,t_3$ | $(M,j_1,j_2)$, $\delta_\infty$ |
|---|---|---|
| $(\tfrac12,0)$ | $I_n$, $III$, $III$ | $(4,3,5)\ \tfrac12$; $(6,5,7)\ \tfrac13$; $(3,2,4)\ \tfrac23$; $(2,1,3)\ 1$ |
| $(-\tfrac12,0)$ | $I_n$, $III$, $III$ | $(1,0,0)\ 0$; $(4,-1,1)\ \tfrac12$; $(6,-1,1)\ \tfrac13$; $(3,-1,1)\ \tfrac23$ |
| $(0,\tfrac12)$ | $III$, $I_n$, $I_n$ | $(4,3,3)\ 0$; $(12,7,11)\ \tfrac13$; $(12,5,13)\ \tfrac23$; $(2,0,3)\ \tfrac32$ |
| $(0,-\tfrac12)$ | $III$, $I_n$, $I_n$ | $(4,1,1)\ 0$; $(2,0,1)\ \tfrac12$; $(12,1,5)\ \tfrac13$; $(12,-1,7)\ \tfrac23$ |
| $(\tfrac12,-\tfrac12)$, $(-\tfrac12,\tfrac12)$ | three $III$, two gauges | $(4,3,3)\ 0$; $(4,1,1)\ 0$ |
| $(0,\pm1)$ | control: $\rho_r\in\mathbf Z$, apparent points expected | $(4,3,5)$, $(4,-1,1)$, $\delta_\infty=\tfrac12$ |

Because $\rho_p=\pm\tfrac12\notin\mathbf Z$ and $\rho_r=0$ forces a logarithm, the classes
in the first two rows admit **no apparent singularity at a finite point**; only $\infty$
can be apparent, and only when $\delta_\infty\in\mathbf Z$. They therefore carry genuine
five-point rows automatically.

> **What they are.** The quadratic character ramified at exactly $\{t_2,t_3\}$ twists the
> class $(0;M',j_1',j_2')$ into $(\tfrac12,0;M,j_1,j_2)$: it turns two of the $I_n$ fibres
> into $I_n^*$ and leaves the *projective* local system unchanged. So these classes are not
> a new geometry; they are a **new integral form** of the same geometry, and — exactly as
> `FOUR_TERM_SCAN.md` §5.5 note 3 records for the $\pm\tfrac12$ twist pairs — their Apéry
> limits differ from those of the untwisted rows.

---

## 3. Theorems D1 and D2: no decayer with a rational characteristic root

**Theorem D1.** *Let $(n+1)^2u_{n+1}=P(n)u_n-Q(n)u_{n-1}+R(n)u_{n-2}$ be an integral
four-term row with $g=CM^2\ne0$, and let $|\lambda_1|\ge|\lambda_2|\ge|\lambda_3|$ be its
characteristic roots. If $|\lambda_2|<1$ then $\chi(\lambda)=\lambda^3-a\lambda^2+d\lambda-g$
is irreducible over $\mathbf Q$.*

*Proof.* Suppose $\chi$ has a rational root. $\chi$ is monic with integer coefficients, so
that root is an integer $r$, and $r\ne0$ because $\chi(0)=-g\ne0$. If $r\in\{\lambda_2,\lambda_3\}$
then $|r|\le|\lambda_2|<1$ with $r$ a non-zero integer — impossible. Hence $r=\lambda_1$ and,
by Gauss's lemma, $\chi=(\lambda-r)(\lambda^2-\sigma\lambda+\pi)$ with $\sigma,\pi\in\mathbf Z$.
Then $\pi=\lambda_2\lambda_3$ and $|\pi|\le|\lambda_2|^2<1$, so $\pi=0$ and $g=r\pi=0$ —
contradiction. $\square$

**Corollary D1.1.** *No class with $\rho_1\ne\rho_2=\rho_3$ — in particular none of the
twenty mixed classes of §2.1 — contains a row with $|\lambda_2|<1$.* Indeed Corollary F2.1
makes unequal exponents force a reducible cubic, which D1 forbids. $\square$

**Corollary D1.2.** *Every row of box W of `FOUR_TERM_SCAN.md` §6.1 with $|\lambda_2|<1$ has
equal exponents at all three finite points, so its three finite fibres are Galois-conjugate
and of one Kodaira type.* This explains, rather than merely records, why box W is populated
only by the equal-$\rho$ classes.

**Theorem D2.** *If moreover $\lambda_2\in\mathbf R$ (a real decayer) then $\chi$ is
irreducible and totally real, $\Delta=18adg-4a^3g+a^2d^2-4d^3-27g^2>0$, and
$|g|<|\lambda_1|\le|a|+2$.*

*Proof.* An irreducible cubic over $\mathbf Q$ either has three real roots or one real root
and a complex-conjugate pair. In the second case the pair has equal modulus; if the real root
were $\lambda_2$ then $|\lambda_1|=|\lambda_3|$, contradicting $|\lambda_2|\ge|\lambda_3|$ and
$|\lambda_2|<|\lambda_1|$. So all three roots are real and $\Delta>0$. Finally
$|g|=|\lambda_1||\lambda_2\lambda_3|<|\lambda_1|$. $\square$

**Verification.** `09_prop.gp` runs both statements against the $30$ genuine five-point rows
of the first census. Exactly one has $|\lambda_2|<1$,
$$(0;2,1,1),\ (a,c,d,f,C)=(154,42,-128,-12,8),\qquad |\lambda_2|=0.45462\ldots,$$
and its cubic is irreducible with $\Delta=-81\,920\,000<0$ — so its two small roots are the
complex pair $-0.4140\pm0.1878i$ and it is **not** a real decayer. The only other row with
$|\lambda_2|\le1$, $(87,24,-152,-16,16)$, has $|\lambda_2|=1$ exactly and a **reducible**
cubic $(\lambda+1)(\lambda^2-88\lambda-64)$ — precisely the boundary case D1's proof leaves
open. **[verified]**

Consequently: **in the four-term world, a real decayer must live in an equal-$\rho$ class with
an irreducible totally real cubic.** That is what box DW below searches.


---

## 4. The boxes

Four scans were run, on two machines (`snake`, 8 cores; the dev box, 12 cores — the latter
shared with an unrelated $12$-process job for most of the run, so its effective throughput
was about half).

**Cost law.** The scanner cost is essentially proportional to
$\text{AMAX}\cdot\text{CMAX}\cdot\text{DMAX}\cdot\text{FMAX}\cdot\text{GMAX}$; measured on
the dev box, one value of $a$ in the class $\mathbf S_0=(0;2,1,1)$ with
$(\text{CMAX},\text{DMAX},\text{FMAX},\text{GMAX})=(120,4000,240,600)$ costs $298$ s, and
box G of `FOUR_TERM_SCAN.md` costs $\approx2.0$ s per value of $a$ per class. A literal
$100\times$ extension of box G in all five directions is therefore about
$4\times10^{6}$ core-seconds — an order of magnitude more than an overnight run on
$20$ cores. The boxes below were sized to the budget, with the **accessory** directions
$c=u_1$ and $f=Q(0)$ (where the period lives, and where box G's surviving rows sat right at
the boundary: $|c|$ up to $42$ against $\text{CMAX}=40$, $|f|$ up to $64$ against
$\text{FMAX}=80$) given priority over $a$ and $C$.

**Validation of the scanners.** `02_fscan` was checked to return bit-identical hit lists to
`lattice/four_term/03_fscan` on $14$ classes over a common box. `03_fmix` (mixed classes) was
checked against an independent exact brute force (`03b_checkmix.py`: no congruence shortcuts,
$U_n$ to $n=16$ with $v_p$ tested for $p\le13$) on the class $(-\tfrac12,0;1,0,0)$ over
$|r|\le3$, $|a|\le16$, $|c|\le9$, $|d|\le30$, $|f|\le5$, $|C|\le9$: both return **exactly the
same $4$ rows** — soundness *and* completeness. `04_fscan5` likewise (§7.2). **[verified]**

**Speed.** Two optimisations were needed and are in `02_fscan.c`:
$U_5=P(4)U_4-16Q(4)U_3+144R(4)U_2$ is still **linear in $C$** (the first non-linear one is
$U_6$), so the analytic-congruence chain of `FOUR_TERM_SCAN.md` §4.1 extends by
$14400\mid U_5$, which thins the innermost arithmetic progression by one to two orders of
magnitude; and the modular inverses in the $36$- and $576$-congruences depend only on the
class and on $(a,c)$ respectively, so they are hoisted out of the $(d,f)$ loops.

| box | classes | $a$ | $|c|$ | $|d|$ | $|f|$ | $|C|$ | purpose |
|---|---|---|---|---|---|---|---|
| **DG** | the $15$ "carrier" classes (those with $\rho\notin\mathbf Z$, or $\rho=0$ and no integral exponent difference — the ones that admit no apparent singularity) | ascending | $200$ | $2000$ | $400$ | $150$ | general census; new genuine rows and their periods |
| **DW** | all $27$ equal-$\rho$ classes | ascending | $250$ | $2a+5$ | $500$ | $(a+2)/M^2+1$ | the $|\lambda_2|<1$ window — question (b); complete in $(d,C)$ for every $a$ reached |
| **MIX** | the $20$ mixed classes of §2.1 | $\le400$ | $100$ | $3000$ | $200$ | $400$ | new territory; $r$ (the rational characteristic root) sharded in blocks |
| **F5** | $7$ five-term classes | $\le24$ | $12$ | $90$ | $40$ | $90$ ($|g|$), $60$ ($|j|$), $90$ ($|C|$) | six singular points, §7 |

Box DG and box DW run in ascending $a$-blocks $[1,30],[31,70],[71,120],[121,180],\dots$, so
that a prefix of the queue is always a *complete* scan of a low-$a$ box — "prioritise by $a$
then $d$" — and every block writes its own output file, so partial progress survives.


### 4.1 What the boxes produced

Counts at the point the scans were stopped (the queues are $a$-block ordered, so each is a
*complete* scan of the stated prefix).

| box | reach | distinct integral rows | `DISGUISED` | `RESCALED` | `NOT-INTEGRAL` | **`CANDIDATE`** |
|---|---|---|---|---|---|---|
| **DG** + **DW** (all equal-$\rho$) | DG $a\le12$; DW complete to $a\le180$ | $251$ | $232$ $(92.4\%)$ | $1$ | $0$ | $18$ |
| **MIX** (20 mixed classes, partial) | $|r|\le30$ of $20$ classes | $143$ | $93$ | $18$ | $2$ | $30$ |
| **MIX2** (the two productive classes complete, plus $(0,\mp\tfrac12;2,\cdot)$ to $|r|\le20$) | $|r|\le60$, $|a|\le200$, $|c|\le200$, $|d|\le2000$, $|f|\le400$, $|C|\le600$ | $85$ | $0$ | $30$ | $22$ | $24$ |
| **F5** (five-term) | $a\le24$, census mode | $2576$ | — | — | — | $4$ (the rest have a repeated root) |

**The equal-$\rho$ half produced nothing new.** All $18$ `CANDIDATE` rows of the DG/DW sweep
are already in the $30$-row list of `FOUR_TERM_SCAN.md` §5.5 — **zero new equal-$\rho$
five-point rows** in the region covered — and the Catalan battery returns
$\textbf{0}$ hits on their limits. The disguise fraction is $92\%$ here rather than $98\%$
only because the deep box is weighted towards the carrier classes.

**The mixed half produced everything.** Every genuine five-point row with an identified
period found in this run is in a mixed class; see §6.4–§6.8.

### 4.2 Question (b), settled empirically as well as structurally

`12_window.py` over the whole extended window (all $27$ equal-$\rho$ classes, $|c|\le250$,
$|f|\le500$, complete in $(d,C)$ for every $a\le180$ — about $17\times$ the $(c,f)$ volume of
box W of `FOUR_TERM_SCAN.md`) finds **exactly three rows with $|\lambda_2|<1$**:

| class | $(a,c,d,f,C)$ | $|\lambda_2|$ | $\lambda_2$ real? | $k$ | score | verdict |
|---|---|---|---|---|---|---|
| $(0;2,-1,3)$ | $(144,32,96,-32,4)$ | $0.35064$ | yes | $1$ | $+0.0480$ | `DISGUISED` (the $I_0^*$ Padé row of `FOUR_TERM_SCAN.md` §6.1) |
| $(0;2,-1,3)$ | $(36,8,-48,16,4)$ | $0.65496$ | no | $1$ | $-0.5768$ | `DISGUISED` |
| $(0;2,1,1)$ | $(154,42,-128,-12,8)$ | $0.45462$ | **no** ($\Delta=-81\,920\,000$) | $2$ | $-1.2117$ | `GENUINE` |

So the only genuine five-point row with $|\lambda_2|<1$ anywhere in the scanned region is the
one already known, and its two small roots are a **complex** pair — Theorem D2 in action.
Corollary D1.1 removes the twenty mixed classes from consideration a priori, and the census
of §6.7 confirms it: the smallest $|\lambda_2|$ over all mixed `CANDIDATE` rows is
**exactly $1$**, never below. **Question (b): no.**


---

## 5. The complex-fold periods of the first census

Nine of the thirty genuine five-point rows of `FOUR_TERM_SCAN.md` have a
complex-conjugate dominant pair and therefore **no archimedean Apéry limit**; §6 of that
file left their periods uncomputed. `08_fold.py` generalises the connection-constant
machinery of `lattice/k3_period/10_fold.py` from the one elliptic-$K3$ row to an arbitrary
class $(\rho;M,j_1,j_2)$: it builds $\mathcal R,\mathcal S,\mathcal V$ from the parameters
(the identity $\mathcal S=\mathcal R+(1-\rho)t\mathcal R'$ replaces the self-adjointness
$\mathcal S=(t\mathcal R)'$ that holds only for $\rho=0$), continues $(A,A')$ and $(B,B')$
along a radial polyline from inside $|t|<\min|t_i|$, and matches to the local Frobenius
basis, which is logarithmic for $\rho\in\mathbf Z$ and
$\{1+\cdots,\ s^{\rho}(1+\cdots)\}$ for $\rho\notin\mathbf Z$. The fold constant is
$$\xi(t_i)=\frac{\text{coefficient of the singular basis element in }B}
                 {\text{coefficient of the singular basis element in }A}.$$

**Self-test.** On the $K3$ row it reproduces the three published constants of
`K3_ROW_PERIOD.md` — $\xi(1)=\tfrac{2\sqrt2}3L(g,2)$ and
$\xi\bigl(\tfrac{5\pm i\sqrt2}{27}\bigr)=\tfrac{\sqrt2\pm i}3L(g,2)$ — to
**210, 211, 211 digits**. The non-integer-$\rho$ branch was validated separately on the two
$\rho=\pm\tfrac13$ rows: changing the step fraction and the matching offset leaves $\xi$
invariant to $85$ digits. **[verified]**

**Result.** All nine rows succeed; $27$ fold constants, each self-agreeing to
$154$–$156$ digits between working precisions $130$ and $210$
(`out/fold_existing.txt`, `out/census_log.txt`). In every one of the nine the three finite
singular points are one real $t_1$ and a conjugate pair, so each row contributes **three**
real periods: $\xi(t_1)$, $\xi_++\xi_-$ and $(\xi_+-\xi_-)/i$.

| row (class $\_\_$ parameters) | $\xi(t_1)$ | $\xi_++\xi_-$ | $(\xi_+-\xi_-)/i$ |
|---|---|---|---|
| $(0;3,1,2)\_\_11,4,37,3,3$ ($K3$) | $1.05749891119153366755$ | $1.05749891119153366755$ | $0.74776465120092405006$ |
| $(-\tfrac13;1,0,0)\_\_42,22,480,64,-64$ | $-0.77725527344708206317$ | $0.16564785024183274686$ | $-0.10200607178273016814$ |
| $(0;4,1,3)\_\_14,5,97,8,12$ | $0.54014823244415968785$ | $0.29278925650488138314$ | $-0.61000046061746984392$ |
| $(0;6,-1,7)\_\_13,4,432,-24,-48$ | $-0.26217135798428882630$ | $0.06936248989328532180$ | $0.22439947525641380275$ |
| $(0;6,1,5)\_\_17,6,56,0,-12$ | $-0.33823337021568890069$ | $0.37186283020499194746$ | $0.31608325840320146871$ |
| $(\tfrac13;1,1,1)\_\_42,8,480,64,-64$ | $-0.32280970782578483005$ | $0.41857135889687961806$ | $-0.36864649809279812009$ |
| $(0;2,1,1)\_\_16,6,96,12,48$ | $0.61944643234552988734$ | $0.52799493919332617147$ | $0.71089792549773360321$ |
| $(0;2,1,1)\_\_17,6,112,8,24$ | $0.56118534995005655738$ | $0.35102266910015264163$ | $0.44533799322124492439$ |
| $(0;2,1,1)\_\_6,2,64,4,-40$ | $-0.44601313717591903062$ | $0.18789250407736198146$ | $-0.48295435974924497113$ |

Two remarks.

* **$\xi(t_1)=\xi_++\xi_-$ only for the $K3$ row.** In the other eight the real fold
  constant and the conjugate-pair sum are independent numbers, so those rows contribute
  *two* independent real periods each, not one.
* **A cross-row identity.** $\operatorname{Re}\xi(t_2)$ of $(0;4,1,3)\_\_14,5,97,8,12$
  equals the **archimedean Apéry limit** of the *different* row
  $(0;2,1,1)\_\_42,15,441,24,100$ to all $70$ stored digits:
  $$\tfrac12(\xi_++\xi_-)=0.1463946282524406915716543362566885954366848059601759545954853976790837\ldots$$
  (difference $6.1\times10^{-72}$, i.e. $70$ agreeing digits — far beyond coincidence).
  The two operators have different singular points —
  $\lambda=\{3,\ \tfrac{11\pm i\sqrt{115}}2\}$ against $\lambda=\{1,16,25\}$ — so this is
  neither a twist nor a rescaling. **Unexplained; flagged for follow-up.** **[measured]**

**Mixed classes.** `11_foldmix.py` extends the same machinery to the mixed-exponent rows of
§2: it builds $\mathcal R,\mathcal S,\mathcal V$ from Theorem D3 and uses the *correct*
$\rho$ at each point ($\rho_r$ at $t=1/r$, $\rho_p$ at the other two), and it falls back to a
rotated continuation path when the radial one runs into another singular point on the ray.
Validated on $(\tfrac12,0;4,3,5)$, $r=-2$, $(30,2,-128,-28,8)$: $134$–$135$ agreeing digits,
and the fold constant at the singular point nearest the origin,
$\xi(t_3)=0.316925359210354980978895\ldots$, reproduces that row's ordinary archimedean
Apéry limit exactly, as it must. **[verified]**

**The new mixed-class rows.** `11_foldmix.py` was also run on the mixed-class `CANDIDATE`
rows — $17$ rows, $51$ fold constants, each self-agreeing to $132$–$137$ digits, in
`out/fold_mix_new.txt`. Where the row also has an ordinary Apéry limit, the fold constant at
the singular point nearest the origin reproduces it exactly, as it must. **Every fold
constant that is not one of those Apéry limits is unidentified** against the Catalan battery
(the only other "hits" are the identically-zero imaginary parts at real singular points).
**[verified]**

**Identification.** All $27$ fold constants (and the derived sums and differences) were put
through the $145$-constant Catalan battery of §6.2. The only hits are the $K3$ row's, which
come back as the expected $L(g,2)$ multiples at $124$–$125$ digits — a strong positive
control. **Every one of the $24$ genuinely new real periods is unidentified.**
**[verified]**, `out/catalan_ident_folds.log`.

---

## 6. Question (a): Catalan

### 6.1 The one integral three-term row with a Catalan period, and what the four-term world does to it

The project knows exactly **one** integral three-term row in Zagier's sense whose Apéry
limit involves $G$: **Zagier $\mathbf E$**, $(a,b,c)=(12,4,32)$,
$$(n+1)^2u_{n+1}=(12n^2+12n+4)u_n-32n^2u_{n-1},\qquad
  u_n=1,4,20,112,676,4304,28496,\dots,$$
$\lambda=(8,4)$, $k=2$, $\xi=\tfrac12G$
(`ROW_LEDGER.md` Table 1; `SPORADIC_SCAN2.md`; modular source $\Gamma_0(8)$,
`CATALAN_AL_HOSTS.md` §4.1). Everything else in the project that carries $G$ is *not* an
integral three-term row: the level-$16$ rows of `CATALAN_TWO_CLASSES.md` §2.1–2.2 are
Eichler/Eisenstein periods on $X_0(16)$ with **six** singular points (and
$\lambda_2=2\sqrt2\notin\mathbf Q$); the five $\Gamma_0(12)$ placements of
`CATALAN_AL_HOSTS.md` §4.2 also have six — but neither is a five-term row, because the fold
is a *double* root of $\mathcal R_c$ (§7.4); Zudilin's Catalan row is a non-integral
hypergeometric decayer; Nesterenko's $(4,7)$ is a two-term linear form; $\beta(4)_{24}$ is
$\operatorname{Sym}^3(\mathbf E)$.

So the only route to a Catalan period inside the four-term world that is guaranteed to
exist *a priori* — before the mixed-class scan of §6.4 found the real ones — is the disguised
one: Proposition F4's signed binomial transform of
Zagier $\mathbf E$. That family is
$$v_n=\sum_m\binom nm(-\nu)^{n-m}u_m,\qquad
  (n+1)^2v_{n+1}=P v_n-Qv_{n-1}+Rv_{n-2},$$

| $\nu$ | $P$ | $Q$ | $R$ | $v_n$ | $\lambda$ | $k$ | score |
|---|---|---|---|---|---|---|---|
| $1$ | $9n^2{+}9n{+}3$ | $11n^2$ | $-21n(n{-}1)$ | $1,3,13,63,333,1863,\dots$ | $7,3,-1$ | $2$ | $-3.0986$ |
| $2$ | $6n^2{+}6n{+}2$ | $-4n^2$ | $-24n(n{-}1)$ | $1,2,8,32,148,712,\dots$ | $6,2,-2$ | $2$ | $-2.6931$ |
| $-1$ | $15n^2{+}15n{+}5$ | $59n^2$ | $45n(n{-}1)$ | $1,5,29,185,1261,\dots$ | $9,5,1$ | $2$ | $-3.6094$ |
| $-2$ | $18n^2{+}18n{+}6$ | $92n^2$ | $120n(n{-}1)$ | $1,6,40,288,2196,\dots$ | $10,6,2$ | $2$ | $-3.7918$ |

(class $(0;1,0,1)$; char. roots exactly $\{8-\nu,4-\nu,-\nu\}$; integrality **[verified]**
to $n=459$; the twelve-parameter fit is one-dimensional and reproduces Proposition F4
verbatim; the Frobenius obstruction at $\infty$ is $S=0$ *exactly*, so all four are
`DISGUISED`, and three of them are already in `lattice/four_term/out/analysis_all.json`
with `"inf_apparent": true`).

**Proposition D4 (the disguise costs the period).** *The signed binomial transform $T$ of
the companion is **not** the companion of the transformed row.* Writing $L_0B=t$ for Zagier
$\mathbf E$ and $\tilde y(t)=(1+\nu t)^{-1}y(t/(1+\nu t))$ for the gauge, one has
$L'(Ty)(t)=(L_0y)(s)$ with $s=t/(1+\nu t)$, so $L'(TB)=t/(1+\nu t)\ne t$: $TB$ fails the
four-term recurrence at **every** $n\ge1$. The row's own companion $W$ ($w_1=1$, $L'W=t$)
is $W=TB+Ty$ with $(n+1)^2y_{n+1}=P_0y_n-Q_0y_{n-1}+\nu^n$, $y_0=y_1=0$, and the exact
Casoratian $(m+1)^2(a_mb_{m+1}-a_{m+1}b_m)=32^m$ gives
$$\boxed{\ \xi_\nu=\tfrac12G\cdot A(\nu/32)-B(\nu/32)\ },\qquad
  A=\sum a_nt^n,\ B=\sum b_nt^n .$$
**[proved + verified to 70 digits]** against the directly recursed companion for
$\nu=\pm1,\pm2$.

| $\nu$ | $\xi_\nu$ |
|---|---|
| $1$ | $0.4863068456612505674910777191229431581350\dots$ |
| $2$ | $0.5203715170059685771322550428796926054562\dots$ |
| $-1$ | $0.4339161236197114306836470378677542083383\dots$ |
| $-2$ | $0.4131171653710453594220235058204506785721\dots$ |

**All four are unidentified** against the $145$-constant Catalan battery of §6.2
(T1 heights $10^8$, T2 pairs $10^5$, `algdep` to degree $6$) — and structurally they must
be, since $A(\nu/32)$ and $B(\nu/32)$ are values of the period and of its Eichler integral
at an *interior* point of the disc. The control $\tfrac12G$ itself is recovered by the
battery at $79$ digits. **[verified]**

> **The moral.** The four-term world *does* contain integral rows built on the Catalan
> local system — the whole $\nu$-family — and the cusp move even improves the elementary
> score by up to $+0.69$ nats ($-3.3863\to-2.6931$) with $k$ still $2$. But it is a
> four-point local system wearing a fifth point, so it brings no new $\kappa_p$; and the
> canonical Apéry limit is no longer $\tfrac12 G$ but an interior value of the Eichler
> integral. **A disguised row cannot answer (a) even when its local system carries $G$.**


### 6.2 The Catalan battery

`07_catalan.gp` is a Catalan-focused replacement for `lattice/four_term/07_ident.gp`.
Constants are built at $250$ digits, `lindep` is run at $60$, and **every** candidate
relation is verified by reconstructing $\xi$ from it and reporting the agreement, with a
$55$-digit threshold. The battery has **$145$ constants** after value-deduplication:

* $G$, $\pi G$, $G\sqrt2$, $G\sqrt3$, $G/\sqrt2$, $G/\pi$, $G/\pi^2$;
* $\zeta(2)$, $\pi^2$, $\pi^2\sqrt D^{\pm1}$ for $D\mid24$;
* $L(2,\chi_D)$ for all $17$ fundamental discriminants $|D|\le24$;
* $L(3,\chi_D)$ for $|D|\le8$, $\zeta(3)$, $\zeta(3)/\pi$, $\zeta(3)/\pi^2$;
* $\mathrm{Ti}_2(x)=\operatorname{Im}\mathrm{Li}_2(ix)$ at $11$ arguments, and
  $\mathrm{Ti}_2(2)=G+\tfrac\pi2\log2$;
* the elliptic-$K3$ row's period
  $L(g,2)=\tfrac{\pi}{32}\Gamma(\tfrac18)\Gamma(\tfrac38)/(\Gamma(\tfrac58)\Gamma(\tfrac78))$
  and four twists;
* seven $\Gamma$/CM-period products (levels $3,4,6,8,12$);
* $\log^2m$ and $\pi\log m$ for $13$ values of $m$, $28$ products $\log m_1\log m_2$, the
  five real-quadratic-unit pairs, $\pi\arctan\frac1m$ and $\arctan^2\frac1m$.

Tests: **T1** $\operatorname{lindep}(1,\xi,X)$, height $\le10^8$; **T2**
$\operatorname{lindep}(1,\xi,X,Y)$, height $\le10^5$, over a $24$-element Catalan-relevant
sublist *greedily reduced to an independent basis of $22$* at $220$ digits (this is the
guard of `FOUR_TERM_SCAN.md` §6: without it `lindep` returns relations among the basis
rather than about $\xi$) — $231$ pairs per $\xi$, and every reported relation is checked to
have a non-zero coefficient on $\xi$; **T3** the explicit $(1,\xi,G,\zeta(2))$ vectors,
printed whether or not they are relations; **T4** `algdep` to degree $6$ on
$\xi$, $\xi/G$, $\xi/\zeta(2)$, $\xi/\pi^2$, $\xi/L(g,2)$.

**Self-tests** (run separately): the battery correctly recovers $\tfrac73G$,
$\tfrac52-\tfrac37\pi G$, $\tfrac25G-\tfrac34\zeta(2)+\tfrac16$ (also from a $62$-digit
truncated input), $\sqrt5$ (T4), $\tfrac{2\sqrt2}3L(g,2)$ (T1), and $\tfrac12G$ from
Zagier $\mathbf E$ at $79$ digits. It is not silently blind. **[verified]**

### 6.3 What the battery says about the first census

* the $21$ archimedean Apéry limits of `FOUR_TERM_SCAN.md` §5.5: **no hit**, in T1, T2 and
  T4; `algdep` minimal heights $10^8$–$10^{30}$, i.e. lattice noise. The nearest thing to
  an amusement is $\xi\cdot G=1.000852\ldots$ for
  $(\tfrac12;4,5,5)\_\_16,2,64,20,4$ — agreeing to three digits, hence nothing.
* the $24$ new real fold periods of §5: **no hit** (§5).
* the four disguised Catalan rows of §6.1: **no hit** (§6.1), and structurally so.

---

### 6.4 **The answer to (a): six genuine five-point rows with Catalan / $L$-value periods**

The mixed-exponent classes of §2 — which the first census never scanned — contain them.
All six live in the single class
$$(\rho_p,\rho_r;M,j_1,j_2)=(-\tfrac12,0;1,0,0),\qquad R(n)=Cn^2,$$
i.e. **exponent differences $(0,\tfrac12,\tfrac12)$ at the three finite points and
$\delta_\infty=0$**: one $I_n$ fibre and two points whose projective monodromy has order $2$
(Kodaira $III$/$III^*$), with $I_n$'s at $0$ and $\infty$.
Writing the row as $(n+1)^2u_{n+1}=P(n)u_n-Q(n)u_{n-1}+R(n)u_{n-2}$ and $r$ for the rational
characteristic root:

| # | $P(n)$ | $Q(n)$ | $R(n)$ | $u_n$ | $\lambda$ | $k$ | score | $\xi=\lim b_n/a_n$ |
|---|---|---|---|---|---|---|---|---|
| 1 | $16n^2{+}20n{+}8$ | $48n^2{+}16n$ | $-128n^2$ | $1,8,72,640,5744,51840,\dots$ | $4{+}4\sqrt2,\ 8,\ 4{-}4\sqrt2$ | $2$ | $-4.0794$ | $\boxed{\tfrac14G}$ |
| 2 | $14n^2{+}20n{+}8$ | $28n^2{+}16n{+}4$ | $8n^2$ | $1,8,72,704,7204,75744,\dots$ | $6{+}4\sqrt2,\ 2,\ 6{-}4\sqrt2$ | $2$ | $-2.6931$ | $\tfrac12G-\tfrac3{16}\zeta(2)$ |
| 3 | $6n^2{+}10n{+}4$ | $-32n^2{-}24n{-}8$ | $32n^2$ | $1,4,36,288,2484,21744,\dots$ | $4{+}4\sqrt2,\ -2,\ 4{-}4\sqrt2$ | $2$ | $-2.6931$ | $\tfrac38\zeta(2)-\tfrac12G$ |
| 4 | $16n^2{+}20n{+}8$ | $68n^2{+}36n{+}8$ | $32n^2$ | $1,8,60,448,3350,25104,\dots$ | $8,\ 4{+}2\sqrt3,\ 4{-}2\sqrt3$ | $2$ | $-4.0101$ | $\tfrac1{32}\bigl(2\zeta(2)+15L(2,\chi_{-3})\bigr)$ |
| 5 | $17n^2{+}25n{+}10$ | $32n^2{+}24n{+}8$ | $16n^2$ | $1,10,114,1424,18650,\dots$ | $8{+}4\sqrt3,\ 8{-}4\sqrt3,\ 1$ | $2$ | $-2.0693$ | $\tfrac1{16}\bigl(15L(2,\chi_{-3})-6\zeta(2)\bigr)$ |
| 6 | $13n^2{+}20n{+}8$ | $-13n^2{-}6n{-}1$ | $-n^2$ | $1,8,87,1024,12575,\dots$ | $7{+}4\sqrt3,\ -1,\ 7{-}4\sqrt3$ | $2$ | $-2.0000$ | $\tfrac18\bigl(2\zeta(2)-3L(2,\chi_{-3})\bigr)$ |

(Class parameters $(a,c,d,f,C)$ and $r$: $(16,8,48,0,-128)_{r=8}$, $(14,8,28,4,8)_{r=2}$,
$(6,4,-32,-8,32)_{r=-2}$, $(16,8,68,8,32)_{r=8}$, $(17,10,32,8,16)_{r=1}$,
$(13,8,-13,-1,-1)_{r=-1}$.)

**Verification, independent of the scan pipeline** (`15_verify_cat.py`, `16_catconfirm.gp`).
For each row the coefficients were rebuilt from Theorem D3 and

* $u_n\in\mathbf Z$ by exact rational arithmetic **to $n=200$** ($u_{200}$ has $179$–$234$
  decimal digits);
* the exponents recomputed from the Theorem F1 dictionary
  $\rho_i=-T(t_i)/(t_i\mathcal R'(t_i))$: they come out as $(0,-\tfrac12,-\tfrac12)$ to
  $40$ digits, and the cubic discriminant is $2048$, $32768$, $32768$, $768$, $192$,
  $49152$ — all non-zero, so the five singular points are distinct;
* **no apparent singularity is possible**: $\rho_r=0$ and $\delta_\infty=0$ force
  logarithms, and at the two $III$ points the local monodromy has eigenvalue ratio $-1$,
  i.e. order $2$ in $\mathrm{PSL}_2$ — non-trivial. The projective local system genuinely has
  **five** singular points;
* the sharp denominator exponent of the companion is $k=2$ in every case (measured to
  $n=80$), the four-term echo of Theorem R3;
* $\xi$ recomputed by forward recursion at working precision $380$ with the adaptive length
  ($n=600$–$6211$), then matched: **agreement $139$–$141$ digits**;
* **the sequences satisfy no three-term recurrence.** The nine-parameter ansatz
  $L(n)u_{n+1}=P(n)u_n-Q(n)u_{n-1}$ with $L,P,Q$ quadratic has nullspace dimension
  $\mathbf 0$ over $21$ equations, while the twelve-parameter four-term ansatz has
  nullspace dimension exactly $1$. They are irreducibly four-term. **[verified]**

**What this is.** `CATALAN_OBSTRUCTION.md` §3 lists as evasion **(E1)** *"a host with an
irrational second singularity … three-term rows cannot have it (integer roots); rows with
more singular points can"*, and notes that the mixed class $E\pm rT$ is the one that
descends, *"with period $G\pm r'\zeta(2)$"*. Rows 1–3 are exactly that object, written down
as integral recurrences for the first time:

* the singular points of the Catalan rows are $0,\ \tfrac18\ \text{or}\ \pm\tfrac12,\
  \tfrac1{4\pm4\sqrt2}\ \text{or}\ \tfrac1{6\pm4\sqrt2},\ \infty$ — **irrational, conjugate
  over $\mathbf Q(\sqrt2)$**, the field of the level-$8$/$16$ Catalan hosts;
* the $L(2,\chi_{-3})$ rows 4–6 are the same phenomenon over $\mathbf Q(\sqrt3)$, the
  conductor-$3$ field of `ONE_CLASS_TWO_WORLDS.md`;
* rows 2 and 3 realise the **mixed class** $\tfrac12G\mp\tfrac3{16}\zeta(2)$,
  $\tfrac38\zeta(2)-\tfrac12G$ — CDT's three-period architecture — on a five-point host;
* the best Catalan score here is $-2.6931=\log\tfrac12-2$, i.e. **$+0.693$ nats better than
  Zagier $\mathbf E$** ($-3.3863$) and equal to the best Catalan geometry previously known
  anywhere in the project (the $(4,2)$ placements on $\Gamma_0(12)$,
  `CATALAN_AL_HOSTS.md` §4.2). Row 6 reaches score $-2.0000$ with $|\lambda_2|=1$ exactly.

**Which point is the fold.** For rows 1, 2, 3, 5, 6 the dominant characteristic root — i.e.
the singular point nearest the origin, the one that governs $\lim b_n/a_n$ — is one of the two
$III$ points, an **order-$2$ orbifold point**, not a cusp; for row 4 it is the cusp $t=\tfrac18$.
`CATALAN_AL_HOSTS.md` §5 found that on every proper Atkin–Lehner quotient carrying a row the
fold is an order-$2$ orbifold point *and concluded that there is then no Eichler period and no
conditional function*. Rows 1–3 nonetheless have honest archimedean Apéry limits at such a
fold, equal to $\tfrac14G$ and to $\tfrac12G\mp\tfrac3{16}\zeta(2)$-type combinations. **Whether
the conditional-function machinery of `CATALAN_OBSTRUCTION.md` §2 applies to them is not
settled here** — that is exactly what the $\sigma_p$ computation below has to decide.

**What it is not.** Every one of the six has $|\lambda_2|\ge1$, so **none is a decayer** —
consistent with, and forced by, Theorem D1 (a mixed class has a rational characteristic root,
and D1 forbids $|\lambda_2|<1$ then). They are *engines*. Their $\kappa_p$, $\sigma_p$ and
$p$-adic Apéry limits are computed in **§6.8**, and the verdict there is that the arithmetic
is structurally new but *strictly poorer* than Zagier $\mathbf E$'s as a lattice resource.



### 6.5 The **whole period matrix** of the Catalan rows is classical

The Apéry limit is only the fold constant at the singular point nearest the origin.
`11_foldmix.py` computes the constant at *all three*, and the Catalan battery identifies
**every one of them**, real and imaginary part, at $125$–$129$ digits
(`out/catalan_ident_foldmix.log`). Ordering $t_1,t_2,t_3$ by decreasing $|t|$, so that
$\xi(t_3)$ is the archimedean Apéry limit:

**The $\mathbf Q(\sqrt2)$ / Catalan rows.**

| row | $\xi(t_1)$ | $\xi(t_2)$ | $\xi(t_3)=\lim b_n/a_n$ |
|---|---|---|---|
| $(16,8,48,0,-128)_{r=8}$ | $\tfrac18\bigl(2G-3\zeta(2)\bigr)$ | $\tfrac1{32}\bigl(8G+3\zeta(2)\bigr)+\tfrac3{32}\zeta(2)\,i$ | $\tfrac14G$ |
| $(14,8,28,4,8)_{r=2}$ | $\tfrac1{80}\bigl(40G-39\zeta(2)\bigr)+\tfrac9{10}\zeta(2)\,i$ | $\tfrac1{16}\bigl(8G+3\zeta(2)\bigr)+\tfrac38\zeta(2)\,i$ | $\tfrac1{16}\bigl(8G-3\zeta(2)\bigr)$ |

**The $\mathbf Q(\sqrt3)$ / $\chi_{-3}$ rows.**

| row | $\xi(t_1)$ | $\xi(t_2)$ | $\xi(t_3)$ |
|---|---|---|---|
| $(16,8,68,8,32)_{r=8}$ | $\tfrac1{32}\bigl(15L(2,\chi_{-3})-12\zeta(2)\bigr)+\tfrac5{96}\pi^2\sqrt3\,i$ | $\tfrac{15}{32}L(2,\chi_{-3})+\tfrac1{96}\pi^2\sqrt3\,i$ | $\tfrac1{32}\bigl(2\zeta(2)+15L(2,\chi_{-3})\bigr)$ |
| $(8,4,32,8,64)_{r=4}$ | $\tfrac1{32}\bigl(6\zeta(2)+15L(2,\chi_{-3})\bigr)$ | $\tfrac1{32}\bigl(15L(2,\chi_{-3})-6\zeta(2)\bigr)+\tfrac1{48}\pi^2\sqrt3\,i$ | $\overline{\xi(t_2)}$ |

**[verified]** to $125$–$129$ digits each.

Two consequences.

* **The period lattice is spanned by two classical constants.** For the $\sqrt2$ family every
  fold period lies in $\mathbf QG+\mathbf Q\zeta(2)$; for the $\sqrt3$ family in
  $\mathbf Q\zeta(2)+\mathbf QL(2,\chi_{-3})+\mathbf Q\pi^2\sqrt3$ (and
  $\pi^2\sqrt3=6\sqrt3\,\zeta(2)$, so this is $\mathbf QL(2,\chi_{-3})+\mathbf Q(\sqrt3)\zeta(2)$).
  These are *not* generic five-point periods: `FOUR_TERM_SCAN.md` §6 found all $21$
  equal-$\rho$ Apéry limits and this file finds all $24$ new equal-$\rho$ fold periods
  **unidentified**, over the same battery. The mixed class is arithmetically special.
* **The row $(8,4,32,8,64)_{r=4}$ has no archimedean Apéry limit at all** — its dominant
  characteristic roots are the complex pair $2\pm2\sqrt3\,i$ — and yet all three of its fold
  periods are classical. So the phenomenon is a property of the *local system*, not of the
  recurrence's convergence.


### 6.6 The twist partners: the same class, the other gauge, and the lemniscate period

The class $(-\tfrac12,0;1,0,0)$ has a gauge partner $(+\tfrac12,0;1,1,1)$ — the same
projective local system with $\rho\mapsto-\rho$, which is the $\pm\tfrac12$ twist pairing of
`FOUR_TERM_SCAN.md` §5.5 note 3. (It has $R(n)=C(n-1)^2$, so $R(1)=0$: a Casoratian
degeneracy in the sense of (D4), which is why it was initially left out of the class list.
Including it was worth it.) Its integral rows are the **same** $(a,d,f,C)$ as the Catalan
rows with a different accessory $c$:

| $(a,c,d,f,C)$, $r$ | $P(n)$ | $Q(n)$ | $R(n)$ | $u_n$ | $\lambda$ | $k$ | $\xi$ |
|---|---|---|---|---|---|---|---|
| $(16,4,48,0,-128)_{r=8}$ | $16n^2{+}12n{+}4$ | $48n^2{-}16n$ | $-128(n{-}1)^2$ | $1,4,24,160,1136,8384,\dots$ | $4{+}4\sqrt2,\,8,\,4{-}4\sqrt2$ | $2$ | $\dfrac{\Gamma(\tfrac14)^4}{64\pi}$ |
| $(14,2,28,4,8)_{r=2}$ | $14n^2{+}8n{+}2$ | $28n^2{-}16n{+}4$ | $8(n{-}1)^2$ | $1,2,8,48,356,2952,\dots$ | $6{+}4\sqrt2,\,2,\,6{-}4\sqrt2$ | $2$ | $\dfrac{\Gamma(\tfrac14)^4}{64\pi}$ |
| $(6,0,-32,-8,32)_{r=-2}$ | $6n^2{+}2n$ | $-32n^2{+}24n{-}8$ | $32(n{-}1)^2$ | $1,0,4,16,116,800,\dots$ | $4{+}4\sqrt2,\,-2,\,4{-}4\sqrt2$ | $2$ | $\dfrac{\Gamma(\tfrac14)^4}{64\pi}$ |

$$\frac{\Gamma(\tfrac14)^4}{64\pi}=0.8593982272525466034362619724763196497376070564773951120343109864420004099501\dots$$
Integrality **[verified]** to $n=200$; exponents $(0,\tfrac12,\tfrac12)$; $k=2$; and the limit
matched to **$195$, $241$ and $241$ digits** respectively (`15b_verify_lem.py`,
`out/lem_verify.log`). Three *different* rows, with three different singular configurations,
and **one and the same period** — in sharp contrast to the $-\tfrac12$ gauge, where the three
partners give three different constants $\tfrac14G$, $\tfrac12G-\tfrac3{16}\zeta(2)$,
$\tfrac38\zeta(2)-\tfrac12G$.

> **The pairing is the point.** $\Gamma(\tfrac14)^4/\pi$ is the period of the CM elliptic
> curve $y^2=x^3-x$ (CM by $\mathbf Z[i]$, $j=1728$) — the *algebraic-period* side of the
> Gaussian world — and $G=L(2,\chi_{-4})$ is its *$L$-value* side. The two half-integral
> gauges of one five-point projective local system over $\mathbf Q(\sqrt2)$ produce one each.
> Both are $\mathbf Q(i)$ objects, and $\sqrt2$ is the ramification of $\mathbf Q(i)/\mathbf Q$
> made visible in the singular points. This is the cleanest structural statement the run
> produced, and it is what makes the six rows of §6.4 look like a genuine realisation rather
> than an accident.



### 6.7 The complete census of the two productive classes, and the constant $\xi^*$

The pair of gauge-partner classes was then swept to completion over
$$|r|\le60,\quad |a|\le200,\quad |c|\le200,\quad |d|\le2000,\quad |f|\le400,\quad |C|\le600$$
(`jobs_MIX2.txt`; $69$ distinct integral rows, of which $30$ `RESCALED`, $14$
`NOT-INTEGRAL` at $n=120$, and **$20$ `CANDIDATE`**). The two classes contain **exactly nine
primitive genuine five-point rows each**, and they correspond one-to-one under the twist
$\rho\mapsto-\rho$ — same $(a,d,f,C)$ and same $r$, different accessory $c$:

| $r$ | $(a,d,f,C)$ | $\lambda$ | score | $\xi$ in gauge $\rho_p=-\tfrac12$ ($c$) | $\xi$ in gauge $\rho_p=+\tfrac12$ ($c$) |
|---|---|---|---|---|---|
| $8$ | $(16,48,0,-128)$ | $4{\pm}4\sqrt2,\ 8$ | $-4.0794$ | $\tfrac14G$ $(c{=}8)$ | $\Gamma(\tfrac14)^4/(64\pi)$ $(c{=}4)$ |
| $2$ | $(14,28,4,8)$ | $6{\pm}4\sqrt2,\ 2$ | $-2.6931$ | $\tfrac12G-\tfrac3{16}\zeta(2)$ $(c{=}8)$ | $\Gamma(\tfrac14)^4/(64\pi)$ $(c{=}2)$ |
| $-2$ | $(6,-32,-8,32)$ | $4{\pm}4\sqrt2,\ -2$ | $-2.6931$ | $\tfrac38\zeta(2)-\tfrac12G$ $(c{=}4)$ | $\Gamma(\tfrac14)^4/(64\pi)$ $(c{=}0)$ |
| $8$ | $(16,68,8,32)$ | $8,\ 4{\pm}2\sqrt3$ | $-4.0101$ | $\tfrac1{32}(2\zeta(2){+}15L(2,\chi_{-3}))$ $(c{=}8)$ | $\xi^*$ $(c{=}4)$ |
| $1$ | $(17,32,8,16)$ | $8{\pm}4\sqrt3,\ 1$ | $-2.0693$ | $\tfrac1{16}(15L(2,\chi_{-3}){-}6\zeta(2))$ $(c{=}10)$ | $\xi^*$ $(c{=}2)$ |
| $-1$ | $(13,-13,-1,-1)$ | $7{\pm}4\sqrt3,\ -1$ | $-2.0000$ | $\tfrac18(2\zeta(2){-}3L(2,\chi_{-3}))$ $(c{=}8)$ | $\xi^*$ $(c{=}1)$ |
| $4$ | $(8,32,8,64)$ | $4,\ 2{\pm}2\sqrt3\,i$ | $-3.3863$ | complex pair, no limit | complex pair, no limit |
| $3$ | $(17,123,33,243)$ | $7{\pm}4\sqrt2\,i,\ 3$ | $-4.1972$ | complex pair, no limit | complex pair, no limit |
| $-3$ | $(5,24,12,-144)$ | $4{\pm}4\sqrt2\,i,\ -3$ | $-3.9356$ | complex pair, no limit | complex pair, no limit |

Every row with an archimedean Apéry limit in the $\rho_p=-\tfrac12$ gauge is **identified**;
in the $\rho_p=+\tfrac12$ gauge three give $\Gamma(\tfrac14)^4/(64\pi)$ and three give one
and the same *unidentified* constant. **[verified]**

**The constant $\xi^*$.** The three $\mathbf Q(\sqrt3)$ rows in the $+\tfrac12$ gauge, a
fourth four-term row in a *different* mixed class ($(0,\tfrac12;2,1,2)$, $r=12$,
$(a,c,d,f,C)=(24,6,192,24,144)$, $\lambda=\{12,\,6\pm2\sqrt3\}$), **and**
the five-term six-point row $(0;1,1,1)$, $(a,c,d,f,g,j,C)=(5,2,0,0,-20,-8,-16)$ of §7.4 —
whose characteristic quartic is $(\lambda-4)(\lambda-2)(\lambda+2)(\lambda-1)$, all roots
*rational* — share the value
$$\xi^*=0.7372929961855962401764261978022936979609191432765291852228688851129961062186382\dots$$
**verified to $319$ agreeing digits** between the four-term and the five-term computation
(`20_xistar.py`, `out/xistar.log`), which use different recurrences of different orders and
different singular sets. It is **unidentified** by the $145$-constant battery (T1 height
$10^8$, T2 pairs $10^5$, `algdep` to degree $6$), and a hand `lindep` against
$\Gamma(\tfrac13)^6/\pi^2$, $\Gamma(\tfrac14)^4/\pi$ and rational multiples returns only
height-$10^{55}$–$10^{85}$ noise.

> **Five** rows, in three different classes and two different worlds (five and six singular
> points), over two different fields, with one common transcendental. Whatever $\xi^*$ is, it is not an accident, and
> identifying it is the sharpest single question this run leaves behind.


### 6.8 The arithmetic of the six rows: half-integral slopes, but no new resource

`CATALAN_OBSTRUCTION.md` §3 and `MULTI_PRIME_LATTICE.md` turn on the *arithmetic*
invariants, not the period, so §6.4 is only half an answer. Scripts `17_arith_*.gp`, logs
`out/arith_*.log`; everything exact (rational / $p$-adic), archimedean at $1500$ digits,
measured to $n=800$ ($2$-adic) and $n=400$ (all $p$), with $p$ swept to $71$.

Notation as in `MULTI_PRIME_LATTICE.md` §1:
$\kappa_p=-\lim v_p(\text{den}\,b_n)/n$, $w_p=\lim v_p(\mathrm{Cas}_n)/n$ with
$\mathrm{Cas}_n=a_nb_{n+1}-a_{n+1}b_n$, and $\sigma_p=w_p+2\kappa_p$ the slope of
$v_p\bigl(b_n/a_n-b_{n-1}/a_{n-1}\bigr)$.

| row | $\xi_\infty$ | $g$ | $\kappa_p$ | $w_2$ | $\sigma_2$ | slope primes ($p\le71$) | $\xi_2$ |
|---|---|---|---|---|---|---|---|
| 1 | $\tfrac14G$ | $-128$ | $\kappa_2=-\tfrac12$, else $0$ | $\tfrac92$ | $\tfrac72$ | $\{2\}$ | $\tfrac14\zeta_2(2)$ |
| 2 | $\tfrac12G-\tfrac3{16}\zeta(2)$ | $8$ | $0$ | $2$ | $2$ | $\{2\}$ | $\tfrac12\zeta_2(2)$ |
| 3 | $\tfrac38\zeta(2)-\tfrac12G$ | $32$ | $0$ | $\tfrac72$ | $\tfrac72$ | $\{2\}$ | $-\tfrac12\zeta_2(2)$ |
| 4 | $\tfrac1{32}(2\zeta(2)+15L(2,\chi_{-3}))$ | $32$ | $0$ | $3$ | $3$ | $\{2\}$ | $\tfrac38L_2(2,\chi_{12})$ |
| 5 | $\tfrac1{16}(15L(2,\chi_{-3})-6\zeta(2))$ | $16$ | $0$ | $3$ | $3$ | $\{2\}$ | $\tfrac34L_2(2,\chi_{12})$ |
| 6 | $\tfrac18(2\zeta(2)-3L(2,\chi_{-3}))$ | $-1$ | $0$ | $0$ | $0$ | $\varnothing$ | none |
| **Zagier $\mathbf E$** | $\tfrac12G$ | $32$ | $0$ | $5$ | $5$ | $\{2\}$ | $\tfrac12\zeta_2(2)$ |

The $2$-adic limits are exact to **$1172$–$1200$ $2$-adic digits**, in the project's
`lpgen.gp` normalisation $\Lambda_p(\chi)=L_p(2,\chi\omega^{-1})$, and Theorem F's
**linearity holds exactly** with the transfer
$$G\mapsto\zeta_2(2),\qquad \zeta(2)\mapsto0,\qquad
  L(2,\chi_{-3})\mapsto\tfrac45L_2(2,\chi_{12})$$
(the $\tfrac54$ Euler factor of `EULER_CRITERION.md`): rows 4 and 5 give
$\tfrac{15}{32}\cdot\tfrac45=\tfrac38$ and $\tfrac{15}{16}\cdot\tfrac45=\tfrac34$, exactly as
observed. The control $\xi_2^{\mathbf E}=\tfrac12\zeta_2(2)$ of `EULER_CRITERION.md` §4.1 is
reproduced.

**Three things are genuinely new. [verified]**

1. **Half-integral slopes.** Rows 1 and 3 have $\sigma_2=\tfrac72$ and row 1 has
   $w_2=\tfrac92$. *Every* row in the `MULTI_PRIME_LATTICE.md` §1.2 census has an integer
   $\sigma_p=v_p(c)$. A half-integral Casoratian slope means the two $2$-adic Newton slopes
   of the pair $(a_n,b_n)$ are conjugate over a **ramified quadratic extension** — here
   $\mathbf Q_2(\sqrt2)$: the irrationality of the singular points is visible $2$-adically.
2. **$\sigma_p$ is decoupled from the trailing coefficient.** $\sigma_2=\tfrac72,2,\tfrac72,3,3,0$
   against $v_2(g)=7,3,5,5,4,0$. The three-term law "$\sigma_p=v_p(c)$" **fails for every one
   of the six**, so in the four-term world the slope must be *measured*, not read off. The
   hypothesis of the cancellation identity of `MULTI_PRIME_LATTICE.md` §4.1 is violated here.
   (The $3\times3$ Casoratian still obeys $W_{n+1}=R(n)W_n/(n+1)^2$ exactly, giving slopes
   $v_2(g)$; it is the $2\times2$ one that decouples, and it has **no** closed form —
   $\mathrm{Cas}_n(n+1)^2/(n^2\mathrm{Cas}_{n-1})$ runs $64,144,\tfrac{256}3,\tfrac{271}3,\dots$
   for row 1, against the constant $32$ and the exact $\mathrm{Cas}_n=32^n/(n+1)^2$ for
   Zagier $\mathbf E$.)
3. **A negative $\kappa_2$.** Row 1 is an *integral* row whose numerators are divisible by a
   growing power of $2$: $v_2(a_n)=\tfrac n2+O(1)$ exactly, and in fact
   $v_2(a_n^{(1)})-\tfrac n2=v_2(a_n^{\mathbf E})$ at every checkpoint tested. Then
   $\kappa_2=-\tfrac12$ and $\sigma_2=w_2+2\kappa_2=\tfrac92-1=\tfrac72$ on the nose — the
   invariant $\sigma=w+2\kappa$ confirmed in a regime it was never tested in.

**But the content is old, and the verdict is negative.** The slope *set* is $\{2\}$ for rows
1–5 and empty for row 6 (checked to $p=71$) — exactly as for $\mathbf E$; the $2$-adic values
are the two constants already in the ledger, $\zeta_2(2)$ and $L_2(2,\chi_{12})$, with new
rational scalars. Row 2's $\xi_2=\tfrac12\zeta_2(2)$ is *literally* $\mathbf E$'s value with
$\sigma_2=2$ instead of $5$; rows 4 and 5 sit in $\mathbf F$'s $2$-adic class
($\sigma_2=3=\sigma_2^{\mathbf F}$) but **lack $\mathbf F$'s second slope prime**
$\sigma_3=2$. The engine yield $\mathrm{net}/r=\sum_p\sigma_p\log p-\log|\rho_2|$ is
$$\tfrac12\log2=0.3466,\quad 0.6931,\quad 1.7329,\quad 0.0693,\quad 2.0101,\quad 0$$
for rows 1–6, against $\log\Lambda_{\mathbf E}=2.0794$ for Zagier $\mathbf E$ and
$\log9=2.1972$ for $\mathbf F$: **every one of the six is dominated by the engine it would
replace.**

> **The sharp statement.** The six rows evade the *geometric* hypothesis (E1) of
> `CATALAN_OBSTRUCTION.md` §3 — irrational conjugate singularities, five punctures,
> $+0.693$ nats of archimedean score over Zagier $\mathbf E$ — but they carry **no new
> arithmetic that the two-row machinery can spend**. The Catalan class's deficit is the
> missing $p=3$ partner, and rows 1–3 supply $\{2\}$ again, with a *smaller* $2$-adic
> resource than $\mathbf E$. **§3 of `CATALAN_OBSTRUCTION.md` has been re-decorated, not
> evaded.**

*Not settled:* a closed form for $v_2(\mathrm{Cas}_n)$ (only $w_2n+O(1)$ with the $O(1)$ in
$[-3,3]$ over the tested window, and exactness at multiples of $100$); and the mechanism
behind $v_2(a_n^{(1)})=\tfrac n2+v_2(a_n^{\mathbf E})$, observed at eleven checkpoints, not
proved.

---

## 7. Five-term rows: six singular points

### 7.1 The dictionary, one term further

**Theorem D5 (five-term exponent dictionary; proved, verified).** *For*
$$(n+1)^2u_{n+1}=P(n)u_n-Q(n)u_{n-1}+R(n)u_{n-2}-T(n)u_{n-3},\quad u_0=1,\ u_{-1}=u_{-2}=u_{-3}=0,$$
$P=an^2+bn+c$, $Q=dn^2+en+f$, $R=gn^2+hn+j$, $T=kn^2+ln+m$, *the operator*
$$L=\theta^2-tP(\theta)+t^2Q(\theta+1)-t^3R(\theta+2)+t^4T(\theta+3)
   =t^2\mathcal R\,D^2+t\,\mathcal S\,D+t\,\mathcal V$$
*has*
$$\mathcal R=1-at+dt^2-gt^3+kt^4,\quad
  \mathcal S=1-(a+b)t+(3d+e)t^2-(5g+h)t^3+(7k+l)t^4,$$
$$\mathcal V=-c+(d+e+f)t-(4g+2h+j)t^2+(9k+3l+m)t^3 ,$$
*six singular points $0,t_1,\dots,t_4,\infty$; exponents $(0,0)$ at $0$, $(0,\rho_i)$ at the
roots of $\mathcal R$ with $\rho_i=-T_c(t_i)/(t_i\mathcal R'(t_i))$,
$T_c=\mathcal S-t\mathcal R'=1-bt+(d+e)t^2-(2g+h)t^3+(3k+l)t^4$, and $3-s_1,3-s_2$ at
$\infty$ with $s_i$ the roots of $T$. Fuchs: $\sum_i\rho_i=s_1+s_2-2$. All four exponents
equal $\rho$ iff*
$$b=(1-\rho)a,\quad e=-2\rho d,\quad h=-(1+3\rho)g,\quad l=-(2+4\rho)k,$$
*and then $(j_1+j_2)/M=s_1+s_2=2+4\rho$ for $T(n)=C(Mn-j_1)(Mn-j_2)$.*

The free parameters are $a,c,d,f,g,j,C$: **three** accessory parameters $c=u_1$, $f=Q(0)$,
$j=R(0)$ — exactly $6-3$, as a rank-two Fuchsian equation with six singular points must
have — plus the four position moduli $a,d,g,C$ modulo the scaling $t\mapsto t/\mu$.

The Casoratian is $W_{n+1}=-T(n)W_n/(n+1)^2$ on four solutions, so $T(n_0)=0$ for an
integer $n_0\ge1$ collapses the recurrence's four-dimensional solution space. **This is not
a defect**: the underlying ODE has only a two-dimensional solution space anyway, and — see
§7.3 — the class in which it happens is exactly the all-semistable one.

### 7.2 The scanner and its validation

`04_fscan5.c` scans a class over $(a,c,d,f,g,j,C)$. The chain of analytic congruences goes
one step further than in the four-term case: $U_2,\dots,U_5$ are all **linear in $C$**
(the first non-linear one is $U_6$, which contains $R(5)U_3\ni C\cdot C$), so
$4\mid U_2$, $36\mid U_3$ are filters, $j$ is put on an arithmetic progression modulo $9$
by $36\mid U_3$, and $576\mid U_4$, $14400\mid U_5$ are solved for $C$ by two linear
congruences and a CRT.

**Validation.** On the class $(-\tfrac12;1,0,0)$ and the box
$a\le4,\ |c|\le4,\ |d|\le10,\ |f|\le6,\ |g|\le10,\ |j|\le8,\ |C|\le8$, an independent
exact-integer brute force (`04d_check5.py`, no congruence shortcuts, $U_n$ to $n=16$ with
$v_p$ tested for $p\le13$) finds exactly $17$ integral rows, and `04_fscan5` in census mode
returns **exactly the same $17$** — soundness *and* completeness. **[verified]**
Every one of those $17$ has a **repeated** characteristic root (the quartic factors as
$(\lambda-1)^2(\cdots)$, $(\lambda-2)^2(\lambda+1)^2$, $(\lambda^2-\lambda+1)^2$, …), i.e.
fewer than six singular points.

### 7.3 A six-singular-point row whose period is the elliptic-$K3$ $L$-value

The census of the **all-semistable** class $(\rho;M,j_1,j_2)=(0;1,1,1)$ — all four finite
exponent differences $0$, $\delta_\infty=0$, $T(n)=C(n-1)^2$, so **six $I_n$ fibres** — turned
up

$$\boxed{\ (n+1)^2u_{n+1}=(n^2{+}n)u_n+(14n^2{+}2)u_{n-1}+(20n^2{-}20n{+}8)u_{n-2}
   +8(n{-}1)^2u_{n-3}\ }$$
$$u_0=1:\quad 1,\,0,\,4,\,8,\,40,\,144,\,616,\,2544,\,10936,\,47200,\,206704,\,911840,\dots$$

$(a,c,d,f,g,j,C)=(1,0,-14,-2,20,8,-8)$; integrality **[verified]** by exact rational
arithmetic to $n=400$ ($u_{400}$ has $271$ decimal digits); primitive ($a=1$, so no
rescaling is possible); all $u_n\ge0$. Its characteristic quartic is
$$\lambda^4-\lambda^3-14\lambda^2-20\lambda-8=(\lambda+1)(\lambda+2)(\lambda^2-4\lambda-4),
\qquad \Delta=2048\ne0,$$
$$\lambda=\ 2+2\sqrt2,\ -2,\ -1,\ 2-2\sqrt2\qquad(4.8284,\ 2,\ 1,\ 0.8284),$$
so the six singular points $0,\ -1,\ -\tfrac12,\ \tfrac{\sqrt2-1}2,\ -\tfrac{\sqrt2+1}2,\ \infty$
are **distinct**; the exponents are $(0,0)$ at each of the five finite ones (a logarithm is
forced, so none is apparent) and $(2,2)$ at $\infty$. Sharp $k=2$, and
$$\operatorname{score}=\log\tfrac1{|\lambda_2|}-k=\log\tfrac12-2=-2.6931 .$$

**Its Apéry limit is the elliptic-$K3$ period.** With $b_0=b_{-1}=b_{-2}=b_{-3}=0$,
$b_1=1$, $\xi=\lim b_n/a_n$ converges at rate $|\lambda_2/\lambda_1|=0.4142$ and
$$\xi=0.52874945559576683377595806385986587134115392973305253672438\dots$$
$$\boxed{\ \xi=\frac{\sqrt2}{3}\,L(g,2)=\tfrac12\,\xi_{K3}(t=1)=\operatorname{Re}\xi_{K3}\!\Bigl(\tfrac{5\pm i\sqrt2}{27}\Bigr)\ }$$
where $L(g,2)=\tfrac\pi{32}\Gamma(\tfrac18)\Gamma(\tfrac38)/(\Gamma(\tfrac58)\Gamma(\tfrac78))$
is the critical value of the weight-$3$ CM newform **32.3.d.a** — the period of the
elliptic-$K3$ four-term row of `FOUR_TERM_SCAN.md` §5.3 and `K3_ROW_PERIOD.md`.
**Verified to $110$ digits**, and returned by the Catalan battery as
$\xi=\tfrac13\bigl(L(g,2)\sqrt2\bigr)$ at $110$ digits. **[verified]**

> **Why this matters.** The elliptic-$K3$ four-term row has the complex conjugate dominant
> pair $5\pm i\sqrt2$ and therefore **no archimedean Apéry limit at all**: $L(g,2)$ appears
> there only as a *fold constant*, computed by analytic continuation (§5). In the six-point
> world the same $L$-value appears as an **honest $\lim b_n/a_n$ of an integral row**, and
> with a score $0.955$ nats better ($-2.6931$ against $-3.6479$). This is the first time in
> the project that a critical value of a weight-$3$ cusp form is an ordinary Apéry limit.
>
> The two local systems are not obviously the same: $\{5\pm i\sqrt2,1\}$ over
> $\mathbf Q(\sqrt{-2})$ against $\{2\pm2\sqrt2,-2,-1\}$ over $\mathbf Q(\sqrt2)$, five
> points against six. Since $32.3.d.a$ has CM by $\mathbf Q(\sqrt{-2})$ and $32=2^5$, both
> are $2$-power-level objects, and a correspondence between them is the obvious guess.
> **Identifying it — and running a five-term analogue of the $\mathcal J$-map test on the
> new row — is left open.** A *second* six-point row, $(6,2,0,0,-6,-2,-1)$ with
> $\lambda=\{3\pm2\sqrt2,\pm1\}$, has the **same** period $\tfrac{\sqrt2}3L(g,2)$ and a
> better score, $-2.0000$; see §7.4 and Conjecture D6 in §7.5.

Three further facts about the row. (i) The point $t=\infty$ has exponents $(2,2)$ with a
logarithm, so the operator has a **second MUM point**; the characteristic quartic seen from
there is, after the normalising scaling $t\mapsto t/2$,
$\lambda^4-5\lambda^3+7\lambda^2-\lambda-2=(\lambda-1)(\lambda-2)(\lambda^2-2\lambda-1)$,
roots $1\pm\sqrt2,2,1$ — again $|\lambda_2|=2$, so no decayer there either, and the scan of
the slice $(a,d,g,C)=(5,7,1,-2)$ over $|c|\le40$, $|f|\le60$, $|j|\le80$ finds **no**
integral row: the second MUM point does not carry one. **[verified]**
(ii) $u_n\ge0$ for every $n\le400$ — a positivity phenomenon in the sense of
`POSITIVITY_PROGRAM.md`. (iii) $|\lambda_2|=2$ exactly, so as an *engine* for the two-row
machinery of `MULTI_PRIME_LATTICE.md` its only candidate alignment prime is $2$; its
$\kappa_p$/$\sigma_p$ data is **not** computed here.

### 7.4 The five-term census and the disguise fraction

Seven classes were scanned over $a\le24$, $|c|\le12$, $|d|\le90$, $|f|\le40$, $|g|\le90$,
$|j|\le60$, $|C|\le90$, in **census mode** (the repeated-root rows kept and flagged, so the
fractions can be counted), plus a deeper box $|c|\le20$, $|d|\le150$, $|f|\le60$,
$|g|\le150$, $|j|\le100$, $|C|\le150$ on the productive class.

| class | fibres | integral rows found | with a **repeated** characteristic root | genuine six-point |
|---|---|---|---|---|
| $(0;1,1,1)$, $T=C(n-1)^2$ | six $I_n$ (all-semistable) | $\approx1700$ | $99.6\%$ | $\mathbf 7$ |
| $(-1;1,-1,-1)$, $T=C(n+1)^2$ | six $I_n$, the gauge partner | $1608$ | $1608$ $(100\%)$ | $\mathbf 0$ |

The pattern of `FOUR_TERM_SCAN.md` §5.1 **continues and intensifies**. There the dominant
degeneracy in the four-term world was the *apparent singularity* ($98.2\%$ of hits were
disguised three-term rows) with the Casoratian/constant-coefficient classes empty; here, one
term further, the dominant degeneracy is a **repeated characteristic root** — the singular
points collide and the object is not a six-point one at all — and it accounts for
$99.6\%$–$100\%$. The gauge partner $(-1;1,-1,-1)$ is the exact five-term analogue of
$\mathbf S_{-1}$: **empty**, as Proposition F5 made $\mathbf S_{-1}$ empty in the four-term
world.

The seven genuine six-point rows found, all in $(0;1,1,1)$, parameters $(a,c,d,f,g,j,C)$
(the scan of that class was still running when it was stopped, so this list is not claimed
complete):

| $(a,c,d,f,g,j,C)$ | $\lambda$ | $k$ | score | $\xi$ |
|---|---|---|---|---|
| $(1,0,-14,-2,20,8,-8)$ | $2{+}2\sqrt2,\ -2,\ -1,\ 2{-}2\sqrt2$ | $2$ | $-2.6931$ | $\tfrac{\sqrt2}3L(g,2)$ (§7.3) |
| $(6,2,0,0,-6,-2,-1)$ | $3{+}2\sqrt2,\ -1,\ 1,\ 3{-}2\sqrt2$ | $2$ | $\mathbf{-2.0000}$ | $\tfrac{\sqrt2}3L(g,2)$ again |
| $(8,2,-52,-8,72,24,-16)$ | $12.5918,\ -2.3177,\ -2,\ -0.2741$ | $2$ | $-2.8406$ | $0.2336058844421865619710\dots$ ($=\operatorname{Re}\xi(t_2)$ of a four-term row) |
| $(4,1,-13,-2,9,3,-1)$ | $6.2959,\ -1.1588,\ -1,\ -0.1371$ | $2$ | $-2.1474$ | $0.4672117688843731239421\dots$ unidentified |
| $(9,3,13,2,4,1,-1)$ | $7.2959,\ 1,\ 0.8629,\ -0.1588$ | $2$ | $\mathbf{-2.0000}$ | the same $0.46721176888\dots$ |
| $(5,2,0,0,-20,-8,-16)$ | $4,\ -2,\ 2,\ 1$ | $2$ | $-2.6931$ | $\xi^*$ (§6.7) |
| $(5,2,17,3,31,12,18)$ | $1\pm2\sqrt2\,i,\ 2,\ 1$ | $2$ | $-3.0986$ | none (complex dominant pair) |

$k=2$ for all seven — the same "free integration" phenomenon one term further — and
**every period occurs at least twice**: $\tfrac{\sqrt2}3L(g,2)$ on two rows,
$0.46721176888\dots$ on two, $\xi^*$ on one here and on three four-term rows, and
$0.23360588444\dots$ on one here and as a four-term fold constant (§7.5). Three of the seven
have $|\lambda_2|=1$ exactly, hence score $-2.0000$. The third row
is worth a note: its characteristic set $\{4,2,-2,1\}$ is **exactly** that of the
$\Gamma_0(12)$ hauptmodul placement $c=3$ of `CATALAN_AL_HOSTS.md` §4.2 — but there the roots
$2$ and $-2$ are *doubled* (§7.4) and the minimal row has seven terms, whereas here they are
simple and the row is honestly five-term. Its period is not the $\Gamma_0(12)$ period
$-\tfrac1{27}G$ or $\tfrac49\zeta(2)$, and is unidentified against the $145$-constant
battery. **[verified]**



### 7.5 Conjecture D6: periods repeat across worlds

Five independent numerical identities turned up, in two shapes.

**(A) A six-point Apéry limit equals the real part of a five-point fold constant.**

| # | $\xi=\lim b_n/a_n$ of a **six-point five-term** row | $=\operatorname{Re}\xi(t_2)=\tfrac12(\xi_++\xi_-)$ of a **five-point four-term** row | digits |
|---|---|---|---|
| (i) | $(0;1,1,1)$, $(1,0,-14,-2,20,8,-8)$ | the elliptic $K3$ $(0;3,1,2)$, $(11,4,37,3,3)$: $\tfrac{\sqrt2}3L(g,2)$ | $110$ |
| (ii) | $(0;1,1,1)$, $(6,2,0,0,-6,-2,-1)$ | the same $K3$ row, again — but with score $-2.0000$ against the $K3$ row's $-3.6479$ | $110$ |
| (iii) | $(0;1,1,1)$, $(8,2,-52,-8,72,24,-16)$ | the mixed row $(\tfrac12,0;3,2,4)_{r=27}$, $(23,8,-216,-42,-324)$ | $110$ |

and, inside the four-term world alone,

| (iv) | $\operatorname{Re}\xi(t_2)$ of $(0;4,1,3)$, $(14,5,97,8,12)$ | $=$ the *Apéry limit* of $(0;2,1,1)$, $(42,15,441,24,100)$ | $70$ |

**(B) One Apéry limit shared by five rows in three classes and two worlds** — the constant
$\xi^*$ of §6.7, on the four-term mixed rows $(\tfrac12,0;1,1,1)_{r=8,1,-1}$ and
$(0,\tfrac12;2,1,2)_{r=12}$ and on the five-term row $(0;1,1,1)$, $(5,2,0,0,-20,-8,-16)$,
to $319$ digits.

Two more repetitions inside the five-term list: $\tfrac{\sqrt2}3L(g,2)$ on two rows and
$0.46721176888\ldots$ on two rows (§7.4); and $\Gamma(\tfrac14)^4/(64\pi)$ on three four-term
rows with three different singular configurations (§6.6). **No period found in this run occurs
only once.**

> **Conjecture D6.** *The real part of a complex fold constant of a five-singular-point
> four-term row is, systematically, an honest archimedean Apéry limit $\lim b_n/a_n$ of an
> integral six-singular-point five-term row.*

The mechanism is presumably a quadratic base change: separating a complex-conjugate pair of
singular points adds a puncture and replaces the fold constant by its
$\mathrm{Gal}(\mathbf C/\mathbf R)$-trace $\tfrac12(\xi_++\xi_-)$. Instance (ii) is the
striking one — the same weight-$3$ $L$-value as the elliptic $K3$, at a score $1.65$ nats
better. **[verified]**, not explained.

### 7.6 The modular Catalan hosts are *not* five-term rows: the fold is a double root

`CATALAN_TWO_CLASSES.md` §2 (the level-$16$ host) and `CATALAN_AL_HOSTS.md` §4.2 (the five
$\Gamma_0(12)$ hauptmodul placements) each record an integral power series $A\in\mathbf Z[[t]]$
on a host with **six** cusps and a Catalan period. Since a rank-two operator with six
singular points gives a five-term row, they looked like ready-made members of the §7 census.
They are not, and the reason is uniform and worth recording.

**Level $16$** ($x=\eta(2\tau)\eta(16\tau)^2/(\eta(\tau)^2\eta(8\tau))$,
$A=1-4x^2+16x^3-44x^4+\cdots$, regenerated to $a_{198}$, content $1$, primitive):
the twelve-unknown five-term ansatz with leading $(n+1)^2$ is **inconsistent**
($\operatorname{rank}M=12$, $\operatorname{rank}[M\,|\,b]=13$; the solution fitted from
$n=3..14$ first fails at $n=16$). A scan over (number of terms, degree) finds the minimal
relations to be **six terms at degree $2$** and **five terms at degree $3$**, each with
nullspace dimension exactly $1$. The six-term one has
$$\mathcal R_c=(1+2t)(1+4t)^2(1+4t+8t^2),$$
i.e. the singular *set* is exactly the six predicted cusps and the predicted characteristic
values $\{-2,-4,-2\pm2i\}$ occur — but **the fold $t=-\tfrac14$ is a *double* root of
$\mathcal R_c$**, which is what pushes the term count from five to six. Exponents: $(0,0)$ at
$0,-\tfrac12,\tfrac{-1\pm i}4$ and $(1,1)$ at the fold and at $\infty$;
$\delta_\infty=0$; Fuchs $\sum=4$. **[verified]**

**$\Gamma_0(12)$, placement $3$** ($A=1,2,4,12,36,120,404,1408,\dots$ to $a_{296}$):
**no five-term recurrence at any degree $\le8$**; the minimal one is **seven terms at degree
$2$**, with
$$\mathcal R_c=(t-1)(2t-1)^2(2t+1)^2(4t-1),\qquad
  \chi=(\lambda-4)(\lambda-2)^2(\lambda-1)(\lambda+2)^2 .$$
The predicted characteristic *set* $\{4,2,-2,1\}$ is right; the multiplicities are $1,2,2,1$,
and the two doubled points carry exponents $\{\tfrac12,\tfrac12\}$.
**Placement $6$**: minimal is **six terms at degree $2$**, $\mathcal R_c=-(t+1)(2t-1)(2t+1)(4t+1)^2$,
the fold $t=-\tfrac14$ doubled with exponents $\{\tfrac12,\tfrac12\}$. **[verified]**

**And Proposition D4 holds on all three.** In each case the row's own canonical companion
($w_1=1$, all earlier terms $0$) has a limit that is **not** the Catalan period and is
**unidentified** by the battery, while only the modular/Eichler companion returns the exact
value:

| host | canonical $\xi$ | battery | modular control |
|---|---|---|---|
| level $16$ | $0.6794425027615324994552983059164448566\dots$ ($134$ digits, $k=2$, score $-3.0397$) | **no hit** | $\lim b_n/a_n=-\tfrac12G$, hit at $119$ digits |
| $\Gamma_0(12)$, $c=3$ | $0.8342406042999240582595730095476920289\dots$ ($k=2$, score $-2.6931$) | **no hit** | periods $-\tfrac1{27}G$ and $\tfrac49\zeta(2)$ |
| $\Gamma_0(12)$, $c=6$ | $0.9933\pm0.0002$ only — the fold has exponents $\{\tfrac12,\tfrac12\}$, so $w_n/a_n\to\xi$ like $1/\log n$, not geometrically (Richardson in $1/\log n$ to $n=4\times10^4$) | not identifiable at this precision | periods $\tfrac19G$ and $\tfrac29\zeta(2)$ |

> **The contrast with §6.4 is the point.** On the *modular* Catalan hosts with more than four
> punctures the fold is a repeated root of $\mathcal R_c$, the minimal row is longer than the
> puncture count predicts, and the canonical Apéry limit is not the period. In the
> mixed-exponent four-term class $(-\tfrac12,0;1,0,0)$ the fold is a **simple** root, the row
> is honestly four-term, and the canonical Apéry limit **is** $\tfrac14G$. That is why §6.4
> is a genuinely new realisation and not a re-labelling of a known host.


---

## 8. Status: proved / verified / open

**Proved.**

* Lemma D3.0 and Theorem D3 (the mixed-exponent normalisation forms, $W=\chi[1+\sum\rho_i\lambda_i/(\lambda-\lambda_i)]$), §2.
* Theorem D1 ($|\lambda_2|<1\Rightarrow\chi$ irreducible), Corollary D1.1 (no decayer in any
  unequal-exponent class), Corollary D1.2, Theorem D2 (real $\lambda_2\Rightarrow\Delta>0$), §3.
* Proposition D4 (the signed binomial transform of a companion is not the transformed row's
  companion; the closed form $\xi_\nu=\tfrac12GA(\nu/32)-B(\nu/32)$), §6.1.
* Theorem D5 (the five-term exponent dictionary and the equal-exponent conditions), §7.1.
* The extension of the analytic-congruence chain: $U_2,\dots,U_5$ are all linear in $C$ (the
  first non-linear one is $U_6$), for four-term and five-term rows alike, §4, §7.2.

**Verified (exact finite computation).**

* Theorem D3 against the Theorem F1 dictionary, $15$ shapes $\times$ $6$ parameter points,
  $0$ failures.
* Theorems D1/D2 on the $30$ census rows, and D1 on all $37\,088$ integer cubics with
  $|\lambda_2|<1$ in $|a|\le24$, $|d|\le70$, $|g|\le70$ — $0$ reducible.
* The patched scanner returns bit-identical hit lists to `lattice/four_term/03_fscan` on
  $14$ classes over a common box.
* `04_fscan5` against an independent exact brute force on a small box: both find exactly the
  same $17$ rows, all with a repeated characteristic root.
* The six Catalan / $L$-value rows of §6.4: integrality to $n=200$, exponents to $40$ digits,
  $k=2$ to $n=80$, $\xi$ to $139$–$141$ digits, and **no three-term recurrence**
  (nullspace dimension $0$ against a nine-parameter ansatz over $21$ equations).
* Their three gauge partners of §6.6: $\xi=\Gamma(\tfrac14)^4/(64\pi)$ to $195$–$241$ digits.
* The complete period matrices of §6.5: all three fold constants of four mixed rows, real and
  imaginary parts, identified at $125$–$129$ digits.
* The completeness of the two productive mixed classes over the box of §6.7, and the
  constant $\xi^*$ to $319$ digits from two recurrences of different order.
* The $\sigma_p$ census of §6.8: $\kappa_p$ and $w_p$ to $n=800$ ($p=2$) and $n=400$
  (all $p\le71$), the $2$-adic Apéry limits to $1172$–$1200$ $2$-adic digits.
* The five instances of Conjecture D6 (§7.5), $70$–$319$ digits.
* The six-point row of §7.3: integrality to $n=400$, $\xi=\tfrac{\sqrt2}3L(g,2)$ to $110$
  digits.
* The $27$ fold constants of §5 to $154$–$156$ digits, with the $K3$ row's reproducing the
  published closed forms to $210$–$211$ digits.
* The non-existence of a five-term row for the level-$16$ and $\Gamma_0(12)$ Catalan hosts,
  and the D4 dichotomy on all three, §7.6.
* The five-term census fractions of §7.4: $99.6\%$ repeated-root in the all-semistable class,
  $100\%$ in its gauge partner.

**Open.**

1. **A closed form for $v_2(\mathrm{Cas}_n)$** of the six rows of §6.4 — measured as
   $w_2n+O(1)$ with the $O(1)$ in $[-3,3]$ and exact at multiples of $100$, but no identity;
   and a proof of $v_2(a_n^{(1)})=\tfrac n2+v_2(a_n^{\mathbf E})$, observed at eleven
   checkpoints. The $2\times2$ Casoratian of a four-term row has no analogue of
   $W_{n+1}=R(n)W_n/(n+1)^2$, and finding one would explain the half-integral slopes.
2. **Are the six a family?** The class $(-\tfrac12,0;1,0,0)$ *is* complete over the box of
   §6.6 and contains exactly nine primitive rows, at $r=\pm1,\pm2,\pm3,4,8,8$ — so within that
   box they are sporadic. Whether they are the small members of a one-parameter family living
   at larger $|r|$ or $|a|$ is open, and it would change the picture entirely: the three
   Catalan periods $\tfrac14G$, $\tfrac12G-\tfrac3{16}\zeta(2)$, $\tfrac38\zeta(2)-\tfrac12G$
   are three different points of the two-dimensional space $\mathbf QG+\mathbf Q\zeta(2)$, which
   is what a family would sweep out.
3. **The correspondence behind §7.3 and §7.5**: why does a six-point row over
   $\mathbf Q(\sqrt2)$ have the period of a five-point elliptic $K3$ over
   $\mathbf Q(\sqrt{-2})$ — and why does the second such row do it with a score $1.65$ nats
   better than the $K3$ row's own?
4. **Conjecture D6 (§7.5)**: five independent identities, verified to $70$–$319$ digits, of
   the shape "the real part of a complex fold constant of a five-point four-term row is an
   Apéry limit of a six-point five-term row", plus the constant $\xi^*$ shared by five rows in
   three classes and two worlds. A quadratic base change is the obvious mechanism; nobody has
   written it down. **No period found in this run occurs only once**, which suggests the
   census is seeing families, not sporadic rows. Identifying $\xi^*$ is the sharpest special
   case.
5. **The fold constants of the $+\tfrac12$ gauge** (§6.6) and of the
   $(0,\mp\tfrac12;2,\cdot)$ classes are not computed here — only the $-\tfrac12$ gauge's, in
   §6.5. If the pattern of §6.5 holds, the $+\tfrac12$ periods should all lie in
   $\mathbf Q\Gamma(\tfrac14)^4/\pi+\mathbf Q\cdot(\text{something})$, which would pin
   $\xi^*$.
6. **A $\mathcal J$-map test for mixed classes and for five-term rows.** `06_jtest.gp`
   assumes the equal-$\rho$ parametrisation; none of the new rows has been tested for
   modularity of its projective monodromy.
7. **Completeness.** None of the boxes is complete in the sense box PW of
   `FOUR_TERM_SCAN.md` §4.2 is; the cost law of §4 says a literal $100\times$ extension of
   box G is $\approx4\times10^6$ core-seconds.

---

## 9. Reproduction

```
cd lattice/four_term_deep
gcc -O3 -march=native -o 02_fscan  02_fscan.c      # equal-rho, deep box
gcc -O3 -march=native -o 03_fmix   03_fmix.c       # mixed-rho  (Theorem D3)
gcc -O3 -march=native -o 04_fscan5 04_fscan5.c     # five-term  (Theorem D5)

python3 01_mixcheck.py                 # Theorem D3 against the F1 dictionary
python3 04b_fit5.py                    # the five-term framework via a binomial transform
python3 04d_check5.py -1 2 1 0 0 4 4 10 6 10 8 8   # brute-force check of 04_fscan5

python3 mkjobs.py DG  classes_prod.txt 800 200 2000 400 150 24 > jobs_DG.txt
python3 mkjobs.py DW  classes_eq.txt   250 250    0 500   0 24 > jobs_DW.txt
python3 mkjobs.py MIX classes_mix.txt  400 100 3000 200 400 24 > jobs_MIX.txt
./run_queue.sh jobs_DG.txt  8      # snake
./run_queue.sh jobs_DW.txt  5      # dev box
./run_queue.sh jobs_MIX.txt 4
./run_queue.sh jobs_F5c.txt 2

python3 06_analyse_deep.py --full out/analysis_deep.json out/DG_*.txt out/DW_*.txt out/MIX_*.txt
python3 12_window.py out/DW_*.txt out/DWB_*.txt          # the |lambda_2| < 1 rows
python3 08_fold.py  census 130 210                        # complex folds, equal-rho
python3 11_foldmix.py file out/analysis_deep.json out/fold_new.txt
python3 05b_report5.py < out/F5_all.txt                   # five-term invariants
gp -q -s 2000000000 09_prop.gp                            # Theorems D1/D2
gp -q -s 2000000000 07_catalan.gp                         # the Catalan battery
echo 'XIFILE="/abs/path/limits.txt"; read("07_catalan.gp")' | gp -q -s 2000000000
```

Data files under `lattice/four_term_deep/out/`: `analysis_deep.json` (every hit with its
verdict and invariants), `window_rows.json` (the $|\lambda_2|<1$ rows),
`fold_existing.txt` / `fold_real_periods.txt` (the $27$ fold constants of §5),
`catalan_ident_existing.log`, `catalan_ident_folds.log`, `catsrc_*.log` (§6.1),
`prop_check.log` (§3), `lvl16_*.log` / `lvl12_*.log` (§7.6), `arith_*.log` (§6.8),
`cat_rows_xi.txt`, `cat_verify.log`, `lem_verify.log`, `xistar.log` (§6.4–§6.7).

**Sources.** As in `FOUR_TERM_SCAN.md`, together with `K3_ROW_PERIOD.md`,
`CATALAN_AL_HOSTS.md`, `CATALAN_TWO_CLASSES.md`, `CATALAN_OBSTRUCTION.md`,
`MULTI_PRIME_LATTICE.md`, `ROW_LEDGER.md`.
