# The $K3$ window of the $\mathcal J$-map test ($\deg\mathcal J\le24$)

*Claude (Opus 5), 2026-08-22. Scripts and data: `lattice/herfurtner/15_k3jtest.gp`,
`15_k3rows.py`, `15_k3run.gp`, `out/k3jall.log`.
Closes open item **4** of `HERFURTNER_CLASSIFICATION.md` §8 ("$K3$ and higher").*

---

## 0. Verdict first

`HERFURTNER_CLASSIFICATION.md` §4 ran the $\mathcal J$-map test with the window
$\deg\mathcal J\le12$, the bound appropriate to a **rational** elliptic surface
($\sum e=12$), and §8 item 4 flagged the obvious worry: an elliptic $K3$ with
four singular fibres has $\sum e=24$, so a priori $\deg\mathcal J$ could be as
large as $24$ and a positive answer might have been missed. Two things are
reported here.

* **Proposition K1 (the window was never the constraint).** For *any* elliptic
  surface — rational, $K3$, or higher $\chi$ — whose Picard–Fuchs local system
  is that of a second-order Apéry-like row with four singular points on
  $\mathbf P^1_t$, Riemann–Hurwitz forces
  $$\boxed{\ \deg\mathcal J\ \le\ 6c+4e_3+3e_2-12\ \le\ 12\ }$$
  where $c,e_2,e_3$ count, among the four singular points, those whose
  projective local monodromy is unipotent (a cusp), of order $2$, of order $3$.
  Equality $\deg\mathcal J=12$ happens only in the all-cusp class. So the
  $\deg\le12$ window of §4 is *complete*: the $K3$ case is not excluded
  geometrically, but any such row would already have been found at
  $\deg\mathcal J\le12$. §1.

* **The computation agrees.** The test was re-run with $\deg\mathcal J\le24$,
  $h\le24$, $140$ series terms and $992$ values of $\gamma$ per $h$, on all
  $18$ rows of §6.1 that had a negative verdict (these include all five negative
  root rows), on three Kodaira-inadmissible controls, and on the six known
  positive rows as implementation controls. **Every negative stays negative; every
  positive returns the same $(h,\deg\mathcal J)$ as at $\deg\le12$.** §2–3.

* **Corollary H6 re-checked from the stored data.** Re-reading
  `out/rows_full.json` (149 exactly re-verified rows) and re-running the cheap
  complete pass `13_fastfilter.py`/`14_decay.py` over all $1\,270\,065$ scan hits
  — no re-run of the $75$-minute scan — the rows with positive score and $k\ge2$
  are exactly two: Apéry's $\zeta(2)$ row $(11,3,-1)$ at $+0.40606$ and
  Beukers'/$\sqrt{\text{Apéry}}$ $(136,10,4)$ at $+0.13920$. Nothing else. §4.

---

## 1. Proposition K1: $\deg\mathcal J\le12$, always

Let $\pi:\mathcal E\to\mathbf P^1_t$ be a minimal elliptic surface with section
and non-constant functional invariant $\mathcal J:\mathbf P^1_t\to
\mathbf P^1_j=X(1)$, $\mu:=\deg\mathcal J$, and suppose its singular fibres sit
over the four points $0,t_1,t_2,\infty$ — the singular points of the row's
operator $L$. (Fibres of type $I_0^*$ have trivial image in
$\mathrm{PSL}_2$ and are invisible to $L$; they may occur elsewhere and change
nothing below.)

Write, among the four points,

* $c$ = the number carrying an $I_n$ or $I_n^*$ fibre ($n\ge1$), i.e. a cusp;
* $e_2$ = the number carrying $III$ or $III^*$;
* $e_3$ = the number carrying $II,IV,IV^*,II^*$;

so $c+e_2+e_3\le4$, and $c\ge1$ because $t=0$ is a MUM point.

**Proof.** Kodaira's dictionary for the functional invariant (Miranda,
*The basic theory of elliptic surfaces*, §IV):

* $\mathcal J^{-1}(\infty)$ is exactly the set of $I_n/I_n^*$ fibres, with
  ramification index $n$; hence $\#\mathcal J^{-1}(\infty)=c$ and
  $\sum_{\mathcal J(P)=\infty}(e_P-1)=\mu-c$;
* over $j=0$ the ramification index is $\equiv0\pmod 3$ except at fibres of type
  $II,IV,IV^*,II^*$; so at most $e_3$ preimages have $e_P<3$, giving
  $\#\mathcal J^{-1}(0)\le e_3+\tfrac{\mu-e_3}3$ and
  $\sum_{\mathcal J(P)=0}(e_P-1)\ge\tfrac23(\mu-e_3)$;
* over $j=1728$ the ramification index is even except at $III,III^*$; likewise
  $\sum_{\mathcal J(P)=1728}(e_P-1)\ge\tfrac12(\mu-e_2)$.

Riemann–Hurwitz for $\mathcal J:\mathbf P^1\to\mathbf P^1$ gives
$$2\mu-2\ =\ \sum_P(e_P-1)\ \ge\ (\mu-c)+\tfrac23(\mu-e_3)+\tfrac12(\mu-e_2)
\ =\ \tfrac{13}6\mu-c-\tfrac23e_3-\tfrac12e_2 ,$$
i.e. $\tfrac16\mu\le c+\tfrac23e_3+\tfrac12e_2-2$, which is the boxed bound.
Any further ramification of $\mathcal J$ (over $j\ne0,1728,\infty$, or beyond the
minimum forced above) only strengthens the inequality. Maximising
$6c+4e_3+3e_2-12$ subject to $c+e_2+e_3\le4$ gives $12$, attained only at
$c=4$. $\square$

This is the same computation as the genus formula
$g=1+\tfrac\mu{12}-\tfrac{e_2}4-\tfrac{e_3}3-\tfrac c2$ for a genus-zero
Fuchsian group of index $\mu$ with signature $(0;2^{e_2},3^{e_3};c)$, written so
that it survives extra ramification (i.e. so that it does not assume $t$ is a
Hauptmodul).

### 1.1 The bound class by class

$\rho$ and $\delta_\infty$ are the invariants of `HERFURTNER_CLASSIFICATION.md`
§1; the class $(M;j_1,j_2)$ means $Q(n)=C(Mn-j_1)(Mn-j_2)$,
$P(n)=A\bigl(n^2+\tfrac{2M-j_1-j_2}{2M}n\bigr)+B$.

| class | $(\rho;\delta_\infty)$ | $(c,e_2,e_3)$ | bound on $\deg\mathcal J$ | realised |
|---|---|---|---|---|
| $(1;0,0)$ | $(0;0)$ | $(4,0,0)$ | $\le12$ | $12$ (Zagier $\mathbf A$–$\mathbf F$) |
| $(2;1,1)$ | $(\tfrac12;0)$ | $(2,2,0)$ | $\le6$ | $6$ (#38, #40) |
| $(3;1,1)$ | $(\tfrac13;0)$ | $(2,0,2)$ | $\le8$ | $8$ (#30) |
| $(3;1,2)$ | $(\tfrac12;\tfrac13)$ | $(1,2,1)$ | $\le4$ | — |
| $(4;1,3)$ | $(\tfrac12;\tfrac12)$ | $(1,3,0)$ | $\le3$ | $3$ (#45) |
| $(2;1,3)$ | $(1;1)$ | $\le(4,0,0)$; $(2,0,0)$ if $t_{1,2}$ apparent | $\le12$; **$\le0$, i.e. empty**, if apparent | — |
| $(3;2,4)$ | $(1;\tfrac23)$ | $\le(3,0,1)$; $(1,0,1)$ if apparent | $\le10$; **empty** if apparent | — |
| $(4;3,5)$ | $(1;\tfrac12)$ | $\le(3,1,0)$; $(1,1,0)$ if apparent | $\le9$; **empty** if apparent | — |
| $(8;3,5)$ | $(\tfrac12;\tfrac14)$ | — | **excluded by Theorem H2** | — |
| $(3;2,5)$ | $(\tfrac76;1)$ | — | **excluded by Theorem H2** | — |

Every realised value in the last column meets its bound with equality — as it
must, since equality is exactly the statement that $\mathcal J$ is a Belyi map
(no ramification outside $\{0,1728,\infty\}$ and minimal indices there), which is
Herfurtner's normal form. **The maximum of the bound over all $26$
Kodaira-admissible classes is $12$.**

So Proposition K1 answers §8 item 4 in the negative: raising the window cannot
produce a new positive. An elliptic $K3$ with four singular fibres is perfectly
possible — the fibre data $I_1I_3\,II^*II^*$ has $\sum e=24$, and there
$\deg\mathcal J=1+3=4$, comfortably inside the bound $6\cdot2+4\cdot2-12=8$ for
its signature — but its $\mathcal J$ still has degree $\le12$, so the
$\deg\le12$ test would already have caught it. The difference between $\sum e=12$
and $\sum e=24$ is absorbed by the fibres that contribute to $e$ without
contributing to $\deg\mathcal J$ ($II,III,IV,IV^*,III^*,II^*$ contribute
$2,3,4,8,9,10$ to $\sum e$ and $0$ to $\deg\mathcal J$), not by a larger
$\mathcal J$.

---

## 2. The raised test

`15_k3jtest.gp` is a rewrite of `05_jtest.gp` with the window raised and the
linear algebra made cheap enough to afford it.

**What is tested.** With $q(t)=t\exp(g/y_0)$ the canonical nome of
$L=\theta^2-tP(\theta)+t^2Q(\theta+1)$ and $J=q^{-1}+744+196884q+\cdots$ the
classical modular function, set
$$W\ :=\ t^{\,h}\,J\bigl(\gamma\,q(t)^h\bigr)\ \in\ \mathbf Q[[t]],\qquad
\operatorname{val}W=0 .$$
The row is (projectively) the Picard–Fuchs system of an elliptic surface over its
own $\mathbf P^1_t$, with a cusp of width $h$ at the MUM point, iff for some
$h\ge1$ and some constant $\gamma$
$$V\cdot W\ =\ t^{\,h}\,U,\qquad U,V\in\mathbf Q[t],\quad
\deg U,\deg V\le d_{\max},$$
and then $\mathcal J=U/V$ and $\deg\mathcal J=\max(\deg U,\deg V)$.

**Parameters.**

| | §4 of `HERFURTNER_CLASSIFICATION.md` | here |
|---|---|---|
| $d_{\max}=\deg\mathcal J$ window | $12$ (rational surface) | $\mathbf{24}$ ($K3$) |
| cusp width $h$ | $1\le h\le12$ | $\mathbf{1\le h\le24}$ |
| series length $N$ | $56$ | $\mathbf{140}$ $(\ge2\cdot24+10)$ |
| $\gamma$ candidates per $h$ | $96$ ($\pm m^{\pm h}$, the $24$ bases hard-coded in `05_jtest.gp`) | $\mathbf{992}$ ($\pm m^{\pm h}$, $m$ $7$-smooth $\le4096$: $248$ bases, a **superset** of the old list) |
| linear algebra | `matker` over $\mathbf Q$ | filter over $\mathbf F_p$, $p=2^{61}-1$; certificate over $\mathbf Q$ |
| unknowns / equations | $26$ / $30$ | $\mathbf{25}$ / $\mathbf{115}$ (after eliminating $U$) |

(The §4 text of `HERFURTNER_CLASSIFICATION.md` quotes "$29$ integers $m$,
$116$ candidates per $h$"; the shipped `GAMS` in `05_jtest.gp` has $24$ bases,
so $96$. Either way the new list contains it.)

**Eliminating $U$.** The coefficient of $t^e$ in $VW$ must vanish for every
$e\notin[h,h+d_{\max}]$; that is a linear system in the $d_{\max}+1$ unknowns
$v_0,\dots,v_{d_{\max}}$ alone, with $N-1-d_{\max}=115$ equations. $U$ is then
read off from the coefficients $e\in[h,h+d_{\max}]$.

**Why $\mathbf F_p$ is safe.** All coefficients occurring ($u_n$, $g_n$, the
$q$-expansion of $j$, $\gamma$) are $p$-integral for $p=2^{61}-1$: denominators
are $N$-smooth (they come from $1/n^2$ and $1/m!$ in the Frobenius recursion and
the exponential) and $\gamma=\pm m^{\pm h}$ with $m\le4096$. Reduction mod $p$
is therefore a ring homomorphism, so a kernel over $\mathbf Q$ reduces to a
kernel over $\mathbf F_p$: **the filter has no false negatives.** False positives
are killed by the exact stage.

**The certificate.** Every $(h,\gamma)$ surviving the $\mathbf F_p$ filter is
redone with exact rational arithmetic. $V$ is determined from the **first $50$
equations only** ($50=2d_{\max}+2$, the total number of free coefficients in
$U$ and $V$); the minimal $V$ is extracted as the $\gcd$ of a basis of that
kernel (the kernel is $V_{\min}\cdot\{\deg\le d_{\max}-\deg V_{\min}\}$, so the
$\gcd$ of any basis is $V_{\min}$); $U$ is computed; and then
$$V\,W-t^{\,h}U\ =\ O(t^{140})$$
is verified coefficient by coefficient. That check involves **$65$ further
coefficient identities that were not used to produce the fit** — the "used /
extra" columns of the table below. A positive verdict is thus a $65$-fold
overdetermined certificate; a negative verdict remains, as in §4, only
"no fit within $h\le24$, $\deg\le24$ and the $\gamma$-list", now backed by
Proposition K1.

---

## 3. Results

`out/k3jall.log`, produced by `gp -q 15_k3run.gp`. Runtime $\approx3.6$ min for
the $27$ rows on one core ($\approx10$ s per negative row: $24\times992=23\,808$
$(h,\gamma)$ pairs each).

### 3.1 Implementation controls (rows that were POSITIVE at $\deg\le12$)

| §6.1 | row | class | $(A,B,C)$ | $\deg\le12$ | $\deg\le24$ | $h$ | $\deg\mathcal J$ | $\gamma$ | coeffs fitted / extra verified |
|---|---|---|---|---|---|---|---|---|---|
| 1 | Zagier $\mathbf B$ | $(1;0,0)$ | $(9,3,27)$ | yes $(1,12)$ | **yes** | $1$ | $12$ | $-1$ | $50$ / $65$ |
| 7 | Zagier $\mathbf D$ | $(1;0,0)$ | $(11,3,-1)$ | yes $(1,12)$ | **yes** | $1$ | $12$ | $1$ | $50$ / $65$ |
| 3 | $\sqrt{\mathrm{AZ}(11,5,125)}$ | $(2;1,1)$ | $(88,10,500)$ | yes $(1,6)$ | **yes** | $1$ | $6$ | $-4$ | $50$ / $65$ |
| 18 | $\sqrt{\mathrm{AZ}(9,3,-27)}$ | $(2;1,1)$ | $(72,6,-108)$ | yes $(3,6)$ | **yes** | $3$ | $6$ | $64$ | $50$ / $65$ |
| 4 | NEW, $I_1I_7\,II\,II$ | $(3;1,1)$ | $(117,21,441)$ | yes $(1,8)$ | **yes** | $1$ | $8$ | $-9$ | $50$ / $65$ |
| 6 | NEW, $I_3\,III\,III\,III$ | $(4;1,3)$ | $(72,6,108)$ | yes $(3,3)$ | **yes** | $3$ | $3$ | $-8$ | $50$ / $65$ |

Same $h$, same $\deg\mathcal J$, in every case; the certificate is $65$ coefficient
identities beyond the $50$ used to fit.

### 3.2 (a) + (b): every §6.1 row that was NEGATIVE, and the negative root rows

$\gamma$ ranges over $992$ values for each of $h=1,\dots,24$: $23\,808$
$(h,\gamma)$ pairs per row, each giving a $115\times25$ system over
$\mathbf F_{2^{61}-1}$.

| §6.1 | row | class | $(A,B,C)$ | $\deg\le12$ | $\deg\le24$, $h\le24$ | fit? |
|---|---|---|---|---|---|---|
| 8 | **Beukers $=\sqrt{\text{Apéry}}$** (root row) | $(2;1,1)$ | $(136,10,4)$ | no | **no** | — |
| 9 | **$\sqrt T$** (root row) | $(2;1,1)$ | $(24,2,4)$ | no | **no** | — |
| 16 | **$\sqrt{\text{Domb}}$** (root row) | $(2;1,1)$ | $(20,2,16)$ | no | **no** | — |
| 2 | **$\sqrt{\mathrm{AZ}(7,3,81)}$** (root row) | $(2;1,1)$ | $(56,6,324)$ | no | **no** | — |
| 13 | **$\sqrt{s_7}$** (Cooper, root row) | $(3;1,2)$ | $(26,2,-3)$ | no | **no** | — |
| 10 | level-5 Fricke row | $(4;1,3)$ | $(88,6,-4)$ | no | **no** | — |
| 5 | new, not a surface | $(3;1,2)$ | $(40,4,48)$ | no | **no** | — |
| 14 | new, not a surface ($\mathcal I=-50$) | $(3;1,1)$ | $(180,24,-72)$ | no | **no** | — |
| 17 | new, not a surface | $(4;1,3)$ | $(28,2,-8)$ | no | **no** | — |
| 20 | new, not a surface | $(4;1,3)$ | $(80,6,36)$ | no | **no** | — |
| 23 | new, not a surface | $(4;1,3)$ | $(48,4,32)$ | no | **no** | — |
| 26 | new, not a surface | $(4;1,3)$ | $(68,6,72)$ | no | **no** | — |
| 21 | $A=0$ | $(2;1,3)$ | $(0,-4,-64)$ | no | **no** | — |
| 22 | $A=0$ | $(2;1,3)$ | $(0,4,-64)$ | no | **no** | — |
| 24 | $A=0$ | $(3;2,4)$ | $(0,-3,-81)$ | no | **no** | — |
| 25 | $A=0$ | $(3;2,4)$ | $(0,3,-81)$ | no | **no** | — |
| 27 | $A=0$ | $(4;3,5)$ | $(0,-4,-256)$ | no | **no** | — |
| 28 | $A=0$ | $(4;3,5)$ | $(0,4,-256)$ | no | **no** | — |

All $18$ rows of §6.1 with a negative verdict; the five negative root rows of (b)
are rows 8, 9, 16, 2, 13 and are a subset of this list, as they must be. **No fit
was found for any of them** — no $h\le24$, no $\gamma$ in the list, no degree
$\le24$. There is nothing to certify (a negative verdict has no certificate; see
Proposition K1 for what replaces one).

### 3.3 (c): the Kodaira-inadmissible controls

Expected negative by Theorem H2 — these three rows are not admissible at all, so
a positive here would have indicated a bug in the raised test.

| row | class $(M;j_1,j_2)$ | $P(n)$ | $Q(n)$ | obstruction | $\deg\le24$, $h\le24$ |
|---|---|---|---|---|---|
| $\sqrt{s_{10}}$ | $(8;3,5)$ | $24n^2+12n+2$ | $-4(8n-3)(8n-5)$ | $\delta_\infty=\tfrac14$ | **no** |
| $\sqrt{s_{18}}$ | $(8;3,5)$ | $56n^2+28n+6$ | $12(8n-3)(8n-5)$ | $\delta_\infty=\tfrac14$ | **no** |
| `NONCONGRUENCE_SCAN` §4.4, $\alpha=18$ | $(3;2,5)$ | $18n^2-3n-6$ | $-3(3n-2)(3n-5)$ | $\rho=\tfrac76\equiv\tfrac16$ | **no** |

(The §4.4 row is the $\arctan$/Padé family $(\alpha,-\tfrac\alpha6,-\tfrac\alpha3,
-27,63,-30)$, integral for $\alpha\equiv18\bmod36$; $\alpha=18$ is its first
member, $u_n=1,-6,-15,-108,-1071,-12474,\dots$. The $\sqrt{s_{10}}$ and
$\sqrt{s_{18}}$ parameters reproduce the census sequences
$1,2,34,588,12726,\dots$ and $1,6,90,1716,36774,\dots$ of
`SPORADIC_SCAN2.md` §8, which is the check that the class encoding
$(M;j_1,j_2)=(8;3,5)$ is the right one.)

### 3.4 Reading the table

1. **No new positives.** All $18$ rows of §6.1 with a negative verdict at
   $\deg\le12$ are still negative at $\deg\le24$, with $h$ up to $24$, $140$
   series terms and $10\times$ as many $\gamma$'s. In particular the five
   negative root rows — Beukers'/$\sqrt{\text{Apéry}}$, $\sqrt T$,
   $\sqrt{\text{Domb}}$, $\sqrt{\mathrm{AZ}(7,3,81)}$, $\sqrt{s_7}$ — are not
   Picard–Fuchs systems of an elliptic $K3$ over their own base either. Theorem
   H4 stands verbatim with "rational elliptic surface" replaced by "elliptic
   surface".

2. **The controls behave.** The six known positives return exactly the same
   $(h,\deg\mathcal J)$ as at $\deg\le12$: $(1,12)$ for Zagier $\mathbf B$ and
   $\mathbf D$, $(1,6)$ for $\sqrt{\mathrm{AZ}(11,5,125)}$, $(3,6)$ for
   $\sqrt{\mathrm{AZ}(9,3,-27)}$, $(1,8)$ for the $I_1I_7\,II\,II$ row and
   $(3,3)$ for the $I_3\,III\,III\,III$ row. Widening the window did **not**
   produce a spurious larger-degree fit for any of them, and did not lower any
   degree: the minimal $\deg\mathcal J$ is an invariant, as claimed in the
   Warning of §4. The $\gamma$'s also agree with the $\deg\le12$ run
   ($-1,\,1,\,-4,\,64,\,-9,\,-8$), although $\gamma$ is not in principle unique
   — it absorbs the scaling $t\mapsto t/\lambda$ of Theorem R1, and the enlarged
   $\gamma$-list here ($992$ values against $96$) simply reaches the same one
   first, both lists being enumerated in increasing $m$ with the same
   sign/inverse cycling.

3. **The Kodaira-inadmissible controls are negative**, as Theorem H2 predicts
   and as they must be: $\sqrt{s_{10}}$ and $\sqrt{s_{18}}$ have
   $\delta_\infty=\tfrac14$ (an order-$8$ local monodromy, impossible in
   $\mathrm{PSL}_2(\mathbf Z)$) and the `NONCONGRUENCE_SCAN.md` §4.4 row has
   $\rho=\tfrac76\equiv\tfrac16$. These three are the sanity check that the
   raised test does not manufacture fits out of a larger search space: with
   $23\,808$ $(h,\gamma)$ pairs per row and a $115$-equation system in $25$
   unknowns, a false positive mod $p$ has probability $O(p^{-1})$ per trial and
   none occurred.

4. **The classes with $\rho\in\mathbf Z$ are the ones Proposition K1 predicts to
   be empty.** All six $A=0$ rows live in $(2;1,3)$, $(3;2,4)$, $(4;3,5)$, the
   three classes whose $t_1,t_2$ carry an apparent singularity; for these the
   bound of §1.1 is $\le0$, so no elliptic surface at all, and the computation
   agrees.

---

## 4. Corollary H6 re-checked (no re-run of the scan)

Two independent passes over the **stored** scan data.

* `out/rows_full.json` — the $149$ rows the scan produced that were re-verified
  exactly to $n=300$, with $k$ measured to $n=60$ and the score
  $\log(1/|\lambda_2|)-k$. Rows with a finite score $>0$: exactly two,
  $$(1;0,0)\ (A,B,C)=(11,3,-1),\ k=2,\ \text{score}=+0.40606\quad
    (\text{Apéry's }\zeta(2)\text{ row}),$$
  $$(2;1,1)\ (A,B,C)=(136,10,4),\ k=2,\ \text{score}=+0.13920\quad
    (\text{Beukers}/\!\sqrt{\text{Apéry}}).$$
  Six rows have score `null` (complex-conjugate characteristic roots, no
  archimedean limit); all remaining $141$ have score $\le0$.

* `13_fastfilter.py` + `14_decay.py` re-run on `out/c_*.txt` and `out/a0.txt`
  (the raw scan output, $\approx15$ s; the $75$-minute `02_hscan` scan itself was
  **not** repeated). Of $1\,270\,065$ hits, $69\,434$ lie in Casoratian-non-degenerate
  classes; $1\,374$ of those have $|\lambda_2|<1$ with real characteristic roots;
  and the number with score $>0$ **and** $k\ge2$ is
  $$\#\{\text{score}>0,\ k\ge2\}\ =\ 2,$$
  namely $(1;0,0)$ $(11,3,-1)$ and $(2;1,1)$ $(136,10,4)$ again. Every other
  positive score in the list has $k=1$ and is a member of the classical
  Legendre–Padé $\log$ family in the classes $(2;-1,1)$ / $(2;1,3)$ with
  $C=4$, $B=\pm A/4$ (top of the ranking: $A=600$, score $+2.5115$, $k=1$).

So **Corollary H6 is confirmed from the data as stored**: beyond Apéry's
$\zeta(2)$ row and Beukers' row there is no third row anywhere in the scan with
positive score and $k\ge2$.

---

## 5. What this changes in `HERFURTNER_CLASSIFICATION.md`

* §8 open item **4** ("$K3$ and higher") is **closed**. The correct statement is
  not "the window was not raised" but "the window is irrelevant": Proposition K1
  bounds $\deg\mathcal J$ by $12$ for every row of the census, whatever the
  Kodaira dimension of the surface. The raised computation is the confirmation.
* Theorem H4 and Theorem H5(iii) can drop the word *rational*: the seven
  $\operatorname{Sym}^1$ square roots other than
  $\sqrt{\mathrm{AZ}(11,5,125)}$ and $\sqrt{\mathrm{AZ}(9,3,-27)}$ are not the
  Picard–Fuchs systems of **any** elliptic surface over their own
  $\mathbf P^1_t$.
* §1.1 above adds a per-class a priori bound on $\deg\mathcal J$, which for
  $(2;1,3)$, $(3;2,4)$, $(4;3,5)$ (the classes with an apparent singularity at
  $t_1,t_2$) is vacuous — those classes contain no elliptic-surface row at all,
  independently of any search.
* Open items 1, 2, 3, 5 of §8 are untouched.

---

## 6. Reproduction

```
cd /home/ubuntu/code/math-modular-sources/lattice/herfurtner

# the K3-window J-map test (needs out/jall.log from the deg <= 12 run, which
# supplies the list of rows that had a NEGATIVE verdict)
python3 15_k3rows.py                       # writes 15_k3run.gp, 27 rows
gp -q -s 8000000000 15_k3run.gp </dev/null > out/k3jall.log     # ~3.6 min

# Corollary H6, from the stored scan output only (no re-run of 02_hscan)
python3 -c "import json;R=json.load(open('out/rows_full.json'));\
print([ (r['M'],r['j1'],r['j2'],r['A'],r['B'],r['C'],r['k'],round(r['score'],5)) \
for r in R if r['score'] is not None and r['score']>0 ])"
python3 13_fastfilter.py | head -3
python3 14_decay.py | tail -4
```

Single-row use of the test from `gp`:

```
? read("15_k3jtest.gp");
? k3test(2,1,1,136,10,4)          \\ Beukers: 0  (no fit, deg <= 24, h <= 24)
? k3test(1,0,0,9,3,27)            \\ Zagier B: [1, 12, -1, 50, 65]
                                  \\  = [h, deg J, gamma, #coeffs fitted, #extra verified]
```

Defaults live at the top of `15_k3jtest.gp`: `DMAX=24`, `HMAX=24`, `NTERM=140`,
`GLIM=4096`, `KPRIME=2^61-1`. `k3test(mm,j1,j2,A,B,C,dmax,hmax,N,p)` takes all
of them as optional arguments, so `k3test(1,0,0,9,3,27,12,12,56)` reproduces the
old $\deg\le12$ setting.

**Files.** `lattice/herfurtner/15_k3jtest.gp` (the test),
`15_k3rows.py` (row selection, reads `out/jall.log`), `15_k3run.gp`
(generated driver), `out/k3jall.log` (the log transcribed in §3).

**Sources.** R. Miranda, *The basic theory of elliptic surfaces*, ETS Editrice
Pisa 1989, ch. IV (the functional invariant and its ramification);
K. Kodaira, *On compact analytic surfaces II, III*, Ann. of Math. 77/78
(1963); S. Herfurtner, *Elliptic surfaces with four singular fibres*,
Math. Ann. **291** (1991) 319–342; I. Shimada, D.-Q. Zhang, *Classification of
extremal elliptic K3 surfaces and fundamental groups of open K3 surfaces*,
Nagoya Math. J. **161** (2001) 23–54 (the $\sum e=24$ side of the picture).
