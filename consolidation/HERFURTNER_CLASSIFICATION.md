# Second-order Apéry-like rows and Herfurtner's elliptic surfaces

*Claude (Opus 5), 2026-08-22. Scripts and data: `lattice/herfurtner/`.
Companion to `SPORADIC_SCAN2.md` (the $\operatorname{Sym}^1$ census),
`ROOT_ROWS.md` (Theorems R1–R4), `NONCONGRUENCE_SCAN.md` (Theorems N1–N3 and the
28 host groups), `SQRT_APERY.md`, and `paper/sections/02_sources.tex`.*

---

## 0. Verdict first

The classification of second-order integral Apéry-like rows splits cleanly into a
**local** question (which shapes of $(P,Q)$ are possible) and a **global** one
(which of the possible shapes actually carry an elliptic surface). Both are now
answered, the first completely and by proof, the second completely relative to
Herfurtner's classification and to an exhaustive integrality scan.

* **Theorem H1 (the shape is forced by the exponents).** Write the row as
  $(n+1)^2u_{n+1}=P(n)u_n-Q(n)u_{n-1}$, $P=an^2+bn+c$, $Q=dn^2+en+f$. Its
  Picard–Fuchs operator $L=\theta^2-tP(\theta)+t^2Q(\theta+1)$ has four singular
  points $0,t_1,t_2,\infty$, with exponents $(0,0)$ at $0$, $(0,\rho_i)$ at
  $t_i$ and $(1-r_1,1-r_2)$ at $\infty$, where $r_i$ are the roots of $Q$.
  Then $\rho_1=\rho_2=\rho$ **iff**
  $$\boxed{\ b=(1-\rho)\,a,\qquad e=-2\rho\,d\ }$$
  and in that case $\rho_1+\rho_2=r_1+r_2=2\rho$ (Fuchs). In particular
  * $\rho=0$ forces $P=a(n^2+n)+c$, $e=0$ — and with $\delta_\infty=0$ also
    $f=0$: **Zagier's normalisation $(an^2+an+b,\;cn^2)$ is not a choice, it is
    the statement that all four fibres are semistable**;
  * $\rho=\tfrac12$ forces $P=a(n^2+\tfrac12n)+c$ — the shape
    $An^2+\tfrac A2n+B$ observed empirically for all nine
    $\operatorname{Sym}^1$ square roots in `SPORADIC_SCAN2.md` §7–8 — and with
    $\delta_\infty=0$ also $Q=C(2n-1)^2$: **the root-row shape is the statement
    that $t_1,t_2$ carry fibres of type III/III\* and $\infty$ an $I_m^*$.**

* **Theorem H2 (Kodaira admissibility).** If $L$ is, up to a rational gauge
  factor and a quadratic twist, the Picard–Fuchs operator of an elliptic surface
  over $\mathbf P^1_t$, then every local exponent difference lies in
  $$\{0\}\cup\{\tfrac12\}\cup\{\tfrac13,\tfrac23\}\pmod 1 ,$$
  because the Kodaira monodromies have eigenvalue ratio $1$ ($I_n,I_n^*$),
  $-1$ ($III,III^*$) or $e^{\pm2\pi i/3}$ ($II,IV,IV^*,II^*$). Consequently
  $$\rho\ \equiv\ 0,\tfrac12,\tfrac13,\tfrac23,\qquad
    \delta_\infty:=r_2-r_1\ \equiv\ 0,\tfrac12,\tfrac13,\tfrac23 \pmod 1 .$$
  This **excludes** three objects of the project's census outright: the
  $\operatorname{Sym}^1$ square roots of Cooper's $s_{10}$ and $s_{18}$
  ($Q=-4(8n-3)(8n-5)$ and $12(8n-3)(8n-5)$, so $\delta_\infty=\tfrac14$: an
  order-$8$ local monodromy, impossible in $\mathrm{PSL}_2(\mathbf Z)$), and the
  unexplained positive-score candidate of `NONCONGRUENCE_SCAN.md` §4.4
  ($\rho=\tfrac76\equiv\tfrac16$). These are precisely three of the rows whose
  periods the project's identification batteries have failed to recognise.

* **Theorem H3 (the cross-ratio invariant, and finiteness).** The Möbius
  invariant of the labelled quadruple $(0,\infty;t_1,t_2)$ is
  $$\boxed{\ \mathcal I\ :=\ \frac{(1+z)^2}{z}\ =\ \frac{a^2}{d},\qquad
    z=\frac{t_1}{t_2}=\frac{\lambda_2}{\lambda_1}\ }$$
  ($\lambda_i$ the characteristic roots, $\lambda^2-a\lambda+d=0$). For a
  **rigid** Herfurtner configuration the four base points are known explicitly,
  so $\mathcal I$ is a computable algebraic number; a row can be the
  Picard–Fuchs system of that configuration only if the two agree. Since
  $a,d\in\mathbf Z$, $\mathcal I\in\mathbf Q$, and every configuration whose
  $\mathcal I$ is irrational or non-real is **excluded by integrality alone**.
  The resulting finite lists are Table 3 below.

* **Consequence for Zagier's six.** In the all-cusp class the rational values of
  $\mathcal I$ available from Herfurtner's rigid list are exactly
  $$\mathcal I\in\{\,0,\ 3,\ \tfrac92,\ -\tfrac{49}8,\ \tfrac{100}9,\
    \tfrac{289}{72},\ -121\,\},$$
  and Zagier's rows $\mathbf B,\mathbf E,\mathbf A,\mathbf C,\mathbf F,\mathbf D$
  realise the last six of them, in that order; the two Beauville
  configurations with an irrational or complex $\mathcal I$
  ($I_1I_1I_5I_5$ with the $I_1$'s split, $65+\tfrac{682}{25}\sqrt5$) are
  unreachable over $\mathbf Q$. So the "six" of Beauville–Zagier is, on the
  recurrence side, the count of **rational cross-ratios in Herfurtner's
  four-cusp block**.

* **Theorem H4 (which root rows are elliptic surfaces).** Of the nine
  $\operatorname{Sym}^1$ square roots of `SPORADIC_SCAN2.md`, exactly **two**
  are Picard–Fuchs systems of a rational elliptic surface over the row's own
  $\mathbf P^1_t$:
  | root row | $(A,B,C)$ | $\mathcal I=A^2/D$ | Herfurtner | $\deg\mathcal J$ | $h$ |
  |---|---|---|---|---|---|
  | $\sqrt{\mathrm{AZ}(11,5,125)}=\sqrt\eta$ | $(88,10,500)$ | $484/125$ | **#38 $I_1I_5\,III\,III$** | 6 | 1 |
  | $\sqrt{\mathrm{AZ}(9,3,-27)}$ | $(72,6,-108)$ | $-12$ | **#40 $I_3I_3\,III\,III$** | 6 | 3 |
  Apéry's own square root ($\mathcal I=1156$), the $T$ row ($36$), Domb's
  ($25/4$), $\mathrm{AZ}(7,3,81)$'s ($196/81$) and Cooper's $s_7$'s
  ($-676/27$) match **no** Herfurtner configuration, and an independent
  $\mathcal J$-map test confirms that their projective monodromy is not
  conjugate into $\mathrm{PSL}_2(\mathbf Z)$: they live on Atkin–Lehner
  quotients $X_0(N)/W$, which are commensurable with but not contained in the
  modular tower. This is the precise sense in which Beukers' row is
  "non-congruence", and it is **not** the same dichotomy as the CM/non-CM one of
  `SPORADIC_SCAN2.md` §8 (Domb's and $T$'s roots have identified CM periods yet
  fail here; $\mathrm{AZ}(9,3,-27)$'s passes here yet its period is
  unidentified).

* **The scan.** An exhaustive integrality scan of every Kodaira-admissible
  normalisation class (26 classes; $|a|\le3000$, $|c|\le150$,
  $|d|\le12000$; integrality tested prime-by-prime to $n=30$ and re-verified
  exactly to $n=300$) recovers every second-order row in the project's census
  and produces **no new row that is a Herfurtner surface**. See §6.

---

## 1. The exponent dictionary

Let $\theta=t\,d/dt$ and let the row be
$$(n+1)^2u_{n+1}=P(n)u_n-Q(n)u_{n-1},\qquad u_0=1,\ u_{-1}=0,$$
with $P=an^2+bn+c$, $Q=dn^2+en+f$. Substituting $y=\sum_{n\ge0}u_nt^n$ gives
$$L=\theta^2-t\,P(\theta)+t^2Q(\theta+1),$$
and in the coordinate $D=d/dt$,
$$L=t^2R(t)\,D^2+t\,S(t)\,D+\bigl(-ct+(d+e+f)t^2\bigr),$$
$$R(t)=1-at+dt^2,\qquad S(t)=1-(a+b)t+(3d+e)t^2 .$$

**Singular points.** $t=0$, the two roots $t_1,t_2$ of $R$ (assume $d\ne0$ and
$a^2\ne4d$, i.e. four distinct singular points), and $t=\infty$. The
characteristic roots of the recurrence are $\lambda_i=1/t_i$, the roots of
$\lambda^2-a\lambda+d$.

**Exponents.**
* At $t=0$: $\theta^2$ gives $(0,0)$ with a logarithm — the MUM point. This is a
  fibre with unipotent local monodromy, i.e. of type $I_{n_0}$, $n_0\ge1$.
* At $t=t_i$ ($R$ has a simple zero there): $0$ is always an exponent, and with
  $$T(t):=S(t)-tR'(t)=1-bt+(d+e)t^2$$
  the second exponent is
  $$\rho_i=-\frac{T(t_i)}{t_i\,R'(t_i)} .$$
  $\rho_1+\rho_2=-e/d$ and $\rho_1\rho_2=-\mathrm{Res}(R,T)/\bigl(d^{\deg T}\cdot
  \tfrac1d\cdot(a^2-4d)\bigr)$, so $\rho_1,\rho_2$ are the roots of an explicit
  rational quadratic.
* At $t=\infty$: writing $y\sim t^{\mu}$, the top coefficient $t^2Q(\theta+1)$
  gives $Q(\mu+1)=0$, so the exponents are $1-r_1,\,1-r_2$ with $r_i$ the roots
  of $Q$; the exponent difference is $\delta_\infty=r_2-r_1$.

**Fuchs relation.** $0+0+\rho_1+\rho_2+(1-r_1)+(1-r_2)=2$, i.e.
$\rho_1+\rho_2=r_1+r_2=-e/d$; verified symbolically in `01_exponents.py`.

### Theorem H1 (normalisation classes)

*If $t_1\ne t_2$ then $\rho_1=\rho_2=\rho$ if and only if $b=(1-\rho)a$ and
$e=-2\rho d$.*

**Proof.** $\rho_i$ is a root of $T(t)+\rho\,tR'(t)$ at $t=t_i$ means
$R\mid T+\rho\,tR'$. Both sides are quadratics with constant term $1$, so
$R=T+\rho tR'$, i.e.
$1-at+dt^2=1-(b+\rho a)t+(d+e+2\rho d)t^2$, giving $b+\rho a=a$ and
$e+2\rho d=0$. $\square$

Thus, up to the scaling $t\mapsto t/\mu$ (which sends
$(a,b,c,d,e,f)\mapsto(\mu a,\mu b,\mu c,\mu^2d,\mu^2e,\mu^2f)$ and
$u_n\mapsto\mu^nu_n$), a normalisation class is fixed by the pair
$(\rho,\delta_\infty)$, and inside it the free parameters are $a$, the leading
coefficient $d$ of $Q$, and the **accessory parameter** $c=P(0)=u_1$ (a rank-two
Fuchsian equation with four singular points has exactly one accessory
parameter). Concretely, writing $Q(n)=C\,(Mn-j_1)(Mn-j_2)$ with
$r_i=j_i/M$:
$$\rho=\frac{j_1+j_2}{2M},\qquad \delta_\infty=\frac{j_2-j_1}{M},\qquad
P(n)=A\Bigl(n^2+\frac{2M-j_1-j_2}{2M}\,n\Bigr)+B,\qquad D=CM^2 .$$
This is the parametrisation used by the scanner `02_hscan.c`.

**Worked instances** (`01_exponents.py`; all exact):

| row | $(M;j_1,j_2)$ | $\rho_1,\rho_2$ | exps at $\infty$ | $\delta_\infty$ | types $(t_1,t_2;\infty)$ |
|---|---|---|---|---|---|
| Zagier $\mathbf A$–$\mathbf F$ | $(1;0,0)$ | $0,0$ | $(1,1)$ | $0$ | cusp, cusp; cusp |
| the nine root rows of AZ type | $(2;1,1)$ | $\tfrac12,\tfrac12$ | $(\tfrac12,\tfrac12)$ | $0$ | III, III; $I_m^*$ |
| $\sqrt{s_7}$ (Cooper) | $(3;1,2)$ | $\tfrac12,\tfrac12$ | $(\tfrac13,\tfrac23)$ | $\tfrac13$ | III, III; II/IV |
| $\sqrt{s_{10}}$, $\sqrt{s_{18}}$ | $(8;3,5)$ | $\tfrac12,\tfrac12$ | $(\tfrac58,\tfrac38)$ | $\tfrac14$ | **not Kodaira** |
| NCS §4.4 candidate | — | $\tfrac76,\tfrac76$ | $(-\tfrac23,\tfrac13)$ | $1$ | **not Kodaira** |

### Theorem H2 (Kodaira admissibility)

*Suppose the local system of $L$ is, up to a rational gauge factor $g(t)$ and a
quadratic character twist, $R^1\pi_*\mathbf Q$ of an elliptic surface
$\pi:\mathcal E\to\mathbf P^1_t$. Then $\rho_1,\rho_2,\delta_\infty\in
\{0,\tfrac12,\tfrac13,\tfrac23\}\pmod1$.*

**Proof.** A rational gauge factor shifts exponents by integers, and a quadratic
twist by half-integers, at each singular point; neither changes an exponent
**difference** mod $1$. The exponent difference at a point equals $\frac1{2\pi
i}\log$ of the ratio of the two eigenvalues of the local monodromy, and Kodaira's
list gives the eigenvalue pairs
$(1,1)$ for $I_n$; $(-1,-1)$ for $I_n^*$; $(\pm i)$ for $III,III^*$;
$(e^{\pm2\pi i/3})$ for $IV,IV^*$; $(e^{\pm\pi i/3})$ for $II,II^*$. The
corresponding ratios are $1,1,-1,e^{\mp4\pi i/3},e^{\mp2\pi i/3}$, i.e.
differences $0,0,\tfrac12,\tfrac13$ or $\tfrac23$ mod $1$. $\square$

The three types are exactly the three kinds of point of a Fuchsian group
commensurable with $\mathrm{PSL}_2(\mathbf Z)$: a cusp, an elliptic point of
order two, an elliptic point of order three. The MUM point is always a cusp.

**Corollary H2.1 (Galois).** If $t_1,t_2$ are conjugate irrationals — i.e.
$a^2-4d$ is not a square — then $\rho_1=\rho_2$ (`NONCONGRUENCE_SCAN.md`
Theorem N3), so the two finite singular fibres have the **same** Kodaira class.
A configuration whose two "movable" fibres have different classes therefore
requires $a^2-4d\in\square$, hence $\lambda_1,\lambda_2\in\mathbf Z$ and
$|\lambda_2|\ge1$: such a row can never have a positive score.

---

## 2. Herfurtner's list

S. Herfurtner, *Elliptic surfaces with four singular fibres*, Math. Ann. **291**
(1991) 319–342, classifies all minimal elliptic surfaces over $\mathbf P^1$ with
section, non-constant $\mathcal J$, and exactly four singular fibres. Because
$e(I_0^*)=6$, two fibres in $T^-=\{I_n^*,IV^*,III^*,II^*\}$ already force
$\sum e\ge14>12$, so for a **rational** elliptic surface ($\sum e=12$) at most
one fibre is in $T^-$; Herfurtner's Table 3 lists $50$ twist classes, and
"transfer of $*$" ($I_n\leftrightarrow I_n^*$, $II\leftrightarrow IV^*$,
$III\leftrightarrow III^*$, $IV\leftrightarrow II^*$) produces the full list of
**56** fibre-type quadruples with $\sum e=12$:

* **38 rigid** configurations, all four fibres in $T^+$; the surface is unique up
  to isomorphism (except $I_1I_6\,II\,III$, which has two Galois-conjugate
  models over $\mathbf Q(\sqrt{-3})$), and $\mathcal J$ is one of 38 Belyi maps
  (= Vidunas–Filipuk H1–H38, = Movasati–Reiter Table 1 rows 1–38);
* **18 in one-parameter families**, exactly one fibre in $T^-$. Seven of these
  have an $I_0^*$ whose position $\varrho_4$ is the modulus while $\mathcal J$
  stays rigid; the other eleven have a moving $\mathcal J$.

$\deg\mathcal J=\sum n_i$ over the $I_n$ and $I_n^*$ fibres. The six
**semistable** configurations
$$I_9I_1I_1I_1,\quad I_8I_2I_1I_1,\quad I_6I_3I_2I_1,\quad I_5I_5I_1I_1,\quad
I_4I_4I_2I_2,\quad I_3I_3I_3I_3$$
are Beauville's, all rigid, all with $\deg\mathcal J=12$.

Full source data (base points, Weierstrass models, $\mathcal J$-maps) is
transcribed in `lattice/herfurtner/04_herfurtner.py`; it was read off the GDZ
scan of the original and independently cross-checked against Miranda's account
of Persson's list and against Vidunas–Filipuk's Table 4.

### 2.1 Which configurations can carry a second-order Apéry-like row

The row imposes: (i) a cusp at $t=0$ (MUM); (ii) the two "movable" points
$t_1,t_2$ carry fibres of the **same** Kodaira class (Corollary H2.1), unless
$a^2-4d$ is a square; (iii) the fourth fibre sits at $t=\infty$. Grouping the
$38$ rigid configurations by the pair (class at $t_1,t_2$; class at $\infty$)
gives the nine **normalisation classes** below; the full listing with all
assignments is `lattice/herfurtner/out/classtable.txt`.

| class | $(\rho;\delta_\infty)$ | recurrence shape | rigid Herfurtner configurations |
|---|---|---|---|
| (cusp, cusp) | $(0;0)$, $(0;1)$, $(1;0)$, $(1;1)$ | $P=a(n^2+n)+c$, $Q=dn^2$ (Zagier) etc. | #13–#18 (Beauville's six) |
| (cusp, ord 2) | $(0;\tfrac12)$, $(1;\tfrac12)$ | $Q=C(4n\mp1)(4n\pm1)$ etc. | #23 $I_1I_1I_7III$, #24 $I_1I_2I_6III$, #25 $I_1I_3I_5III$, #26 $I_2I_3I_4III$ |
| (cusp, ord 3) | $(0;\tfrac13),(0;\tfrac23),(1;\cdot)$ | $Q=C(6n\mp1)(6n\pm1)$, $C(3n\mp1)(3n\pm1)$ | #19–#22 ($\cdot\,II$), #27–#29 ($\cdot\,IV$) |
| (ord 2, cusp) | $(\tfrac12;0)$, $(\tfrac12;1)$ | $P=a(n^2+\tfrac n2)+c$, $Q=C(2n-1)^2$ **(root rows)** | #38 $I_1I_5III\,III$, #39 $I_2I_4III\,III$, #40 $I_3I_3III\,III$ |
| (ord 2, ord 2) | $(\tfrac12;\tfrac12)$ | $Q=C(4n-1)(4n-3)$ | #45 $I_3\,III\,III\,III$ |
| (ord 2, ord 3) | $(\tfrac12;\tfrac13),(\tfrac12;\tfrac23)$ | $Q=C(3n-1)(3n-2)$, $C(6n-1)(6n-5)$ | #44 $I_2\,IV\,III\,III$, #48 $I_4\,II\,III\,III$ |
| (ord 3, cusp) | $(\tfrac13;0)$, $(\tfrac23;0)$ | $Q=C(3n-1)^2$, $C(3n-2)^2$ | #30 $I_1I_7II\,II$, #31 $I_2I_6II\,II$, #32 $I_4I_4II\,II$, #36 $I_1I_5II\,IV$, #37 $I_2I_4II\,IV$, #43 $I_2I_2IV\,IV$ |
| (ord 3, ord 2) | $(\tfrac13;\tfrac12)$, $(\tfrac23;\tfrac12)$ | $Q=C(12n-1)(12n-7)$ etc. | #46 $I_3\,II\,III\,IV$, #49 $I_5\,III\,II\,II$ |
| (ord 3, ord 3) | $(\tfrac13;\tfrac13)$ etc. | $Q=C(6n-1)(6n-3)$ etc. | #47 $I_4\,IV\,II\,II$, #50 $I_6\,II\,II\,II$ |

The one-parameter families enter as follows. $I_0^*$ has trivial image in
$\mathrm{PSL}_2$, so at such a point $L$ has an **apparent** singularity
(exponents $(0,m)$, $m\in\mathbf Z_{>0}$, no log): these are the classes with
$\rho\in\mathbf Z_{>0}$. $I_m^*$ has half-integral exponents, hence can only sit
at $\infty$ in our normalisation (where the exponents $1-r_i$ are free to be
half-integers) — and that is exactly what $Q=C(2n-1)^2$ does. So the root-row
class $(\tfrac12;0)$ meets both the rigid block (#38–#40, with $I_m$ at
$\infty$) and, a priori, families with an $I_m^*$ at $\infty$; but
$\sum e=I_a+3+3+(m+6)=12$ forces $a+m=0$, impossible, so **in the root-row class
only #38, #39, #40 survive**.

---

## 3. The cross-ratio invariant

For the row, the four singular points are $0$, $t_1=1/\lambda_1$,
$t_2=1/\lambda_2$, $\infty$; the labels are: cusp at $0$ (MUM), the pair
$\{t_1,t_2\}$ of equal class, and $\infty$. The Möbius invariant of this
labelled quadruple, symmetric under the only two allowed relabellings
($t_1\leftrightarrow t_2$, and $0\leftrightarrow\infty$ when both are cusps
— both of which invert $z$), is
$$\mathcal I=\frac{(1+z)^2}{z},\qquad z=\frac{t_1}{t_2}=\frac{\lambda_2}{\lambda_1},
\qquad\text{i.e.}\qquad \mathcal I=\frac{(\lambda_1+\lambda_2)^2}{\lambda_1\lambda_2}
=\frac{a^2}{d}.$$
For a Herfurtner configuration with base points $(p_1,\dots,p_4)$ and an
assignment $(p_0,p_\infty\mid p',p'')$ the same quantity is
$\mathcal I=(1+z)^2/z$ with
$z=\frac{(p'-p_0)(p''-p_\infty)}{(p'-p_\infty)(p''-p_0)}$.

**Theorem H3.** *A second-order integral row is the Picard–Fuchs system of a
rigid Herfurtner configuration only if its $\mathcal I=a^2/d$ equals the
$\mathcal I$ of that configuration under some admissible assignment. Since
$a,d\in\mathbf Z$, only configurations with $\mathcal I\in\mathbf Q$ are
reachable.*

The complete lists (computed exactly in `04_herfurtner.py`, printed in
`out/classtable.txt`) are:

| class | rational values of $\mathcal I$ available | discarded (irrational / non-real) |
|---|---|---|
| (cusp, cusp) | $0,\;3,\;\tfrac92,\;-\tfrac{49}8,\;\tfrac{100}9,\;\tfrac{289}{72},\;-121$ | $65+\tfrac{682}{25}\sqrt5$ (#16 with the $I_1$'s split) |
| (cusp, ord 2) | $\tfrac{169}{128},\,-\tfrac43,\,\tfrac{49}{12},\,\tfrac{25}4,\,-\tfrac{14884}{375},\,\tfrac{64009}{16000},\,\tfrac{17161}{384}$ | #23 with the $I_1$'s split (complex) |
| (cusp, ord 3) | $\tfrac{196}{81},-\tfrac{289}{1568},\tfrac{16900}{3969},\tfrac{12769}{2592},-\tfrac{6241}{80},\tfrac{25921}{6480},\tfrac{6724}{81},-\tfrac{484}{135},\tfrac{1369}{160},\tfrac{3481}{864},\tfrac92,0,-\tfrac{529}{50},\tfrac{2704}{675},\tfrac{841}{54}$ | #19 with the $I_1$'s split (complex) |
| **(ord 2, cusp)** | $\tfrac{484}{125}$ (#38), $0$ (#39), $-12$ (#40) | — |
| (ord 2, ord 2) | $3$ (#45) | — |
| (ord 2, ord 3) | $0$ (#44), $\tfrac{100}{27}$ (#48) | — |
| (ord 3, cusp) | $\tfrac{169}{49},0,-50,\tfrac{1681}{400},-\tfrac{49}8$ | — |
| (ord 3, ord 2) | $\tfrac{25}4$ (#46), $\tfrac{121}{64}$ (#49) | — |
| (ord 3, ord 3) | $0,\tfrac92,3$ | — |

### 3.1 Zagier's six, read off the table

| row | $(a,b,c)$ | $\mathcal I=a^2/d$ | Herfurtner configuration | fibre at the MUM point |
|---|---|---|---|---|
| $\mathbf A$ | $(7,2,-8)$ | $-49/8$ | #15 $I_1I_2I_3I_6$ | $I_1$ (level 6) |
| $\mathbf B$ | $(9,3,27)$ | $3$ | #13 $I_1I_1I_1I_9$ | $I_1$ (level 9) |
| $\mathbf C$ | $(10,3,9)$ | $100/9$ | #15 $I_1I_2I_3I_6$ | $I_1$ |
| $\mathbf D$ | $(11,3,-1)$ | $-121$ | #16 $I_1I_1I_5I_5$ | $I_1$ (level 5) |
| $\mathbf E$ | $(12,4,32)$ | $9/2$ | #14 $I_1I_1I_2I_8$ | $I_1$ (level 8) |
| $\mathbf F$ | $(17,6,72)$ | $289/72$ | #15 $I_1I_2I_3I_6$ | $I_1$ |

$\mathcal I$ alone does not separate #13 from #18 (both $3$) nor #14 from #17
(both $9/2$); the $\mathcal J$-map test of §4 does, returning a width $h$ of
the MUM cusp: $h=1$ for $\mathbf B$ (so the MUM fibre is $I_1$, hence #13; #18
has all cusps of width $3$) and $h=1$ for $\mathbf E$ (#14, not #17 whose cusp
widths are $2,2,4,4$). $\mathbf A,\mathbf C,\mathbf F$ are three different rows
on the **same** surface $I_6I_3I_2I_1$ with the **same** MUM fibre $I_1$,
differing only in which of the three remaining fibres is placed at $\infty$
(namely $I_2$, $I_3$, $I_6$, giving $\mathcal I=-\tfrac{49}8,\tfrac{100}9,
\tfrac{289}{72}$) —
which is the recurrence-level explanation of the census entries "level $4/6$",
"level $6/9$", "level $6/8$". The two Beauville configurations not realised are
$I_4I_4I_2I_2$ and $I_3I_3I_3I_3$; their only rational invariants are
$\mathcal I=0$ and $\mathcal I=3$, values already taken by $\mathbf E$ and
$\mathbf B$ (and $\mathcal I=0$ means $a=0$, for which the scan of §6 finds no
integral row in this class).

---

## 4. The $\mathcal J$-map test

The cross-ratio test is necessary but not sufficient: it uses the positions of
the singular points, not the accessory parameter. A sufficient test is
available and is decisive.

**Lemma.** *Let $q$ be the canonical nome of $L$ at the MUM point,
$q=t\exp(g/y_0)$ with $y_1=y_0\log t+g$ the second Frobenius solution. If $t$ is
a Hauptmodul of a genus-zero group $\Gamma$ with $t=c\,q_h+O(q_h^2)$ at a cusp
of width $h$, $q_h=e^{2\pi i\tau/h}$, then $q=c\,q_h$; hence
$e^{2\pi i\tau}=q_h^h=(q/c)^h$. Therefore $\Gamma$ is conjugate into
$\mathrm{PSL}_2(\mathbf Z)$ — equivalently $L$ is (projectively) the
Picard–Fuchs system of an elliptic surface over $\mathbf P^1_t$ — if and only if
for some $h\in\mathbf Z_{>0}$ and some $\gamma=c^{-h}$ the composite
$J\bigl(\gamma q(t)^h\bigr)$, $J=j/1$ the classical $q$-expansion
$q^{-1}+744+196884q+\cdots$, is a **rational function of $t$**; its degree is
then $\deg\mathcal J$.*

`05_jtest.gp` implements this: it builds $q(t)$ from the Frobenius recursion
$$n^2u_n=P(n{-}1)u_{n-1}-Q(n{-}1)u_{n-2},\qquad
n^2g_n=P'(n{-}1)u_{n-1}+P(n{-}1)g_{n-1}-Q'(n{-}1)u_{n-2}-Q(n{-}1)g_{n-2}-2nu_n$$
exactly over $\mathbf Q$ to $56$ terms, forms $W=x^h\,J(\gamma q^h)$ as a power
series, and solves the linear system $V\cdot W=x^hU$, $\deg U,\deg V\le12$, by
`matker`, verifying the fit on more coefficients than were used to produce it.
$\gamma$ ranges over $\pm m^{\pm h}$ for $29$ integers $m$ (the scaling
$t\mapsto t/\lambda$ of the root theorems changes $\gamma$ by $\lambda^{h}$;
this is why the naive $\gamma=1$ version of the test fails on rescaled rows).

**Results** (`out/jtest2.log`, `out/jtest3.log`):

| row | in $\mathrm{PSL}_2(\mathbf Z)$? | $h$ | $\deg\mathcal J$ | $\gamma$ |
|---|---|---|---|---|
| Zagier $\mathbf A$ | yes | 1 | 12 | $1$ |
| Zagier $\mathbf B$ | yes | 1 | 12 | $-1$ |
| Zagier $\mathbf C$ | yes | 1 | 12 | $1$ |
| Zagier $\mathbf D$ | yes | 1 | 12 | $1$ |
| Zagier $\mathbf E$ | yes | 1 | 12 | $1$ |
| Zagier $\mathbf F$ | yes | 1 | 12 | $-1$ |
| $\sqrt{\mathrm{AZ}(11,5,125)}$ | **yes** | 1 | **6** | $-4$ |
| $\sqrt{\mathrm{AZ}(9,3,-27)}$ | **yes** | 3 | **6** | $64$ |
| $\sqrt{\text{Apéry}}$ (Beukers) | no | — | — | — |
| $\sqrt{T}$ | no | — | — | — |
| $\sqrt{\text{Domb}}$ | no | — | — | — |
| $\sqrt{\mathrm{AZ}(7,3,81)}$ | no | — | — | — |
| $\sqrt{s_7}$ (Cooper) | no | — | — | — |

The two positive answers agree with the cross-ratio test of §3 on both counts
($\deg\mathcal J=6$ = the degree of #38 and #40; and $h=1$, $h=3$ identify the
MUM fibre as $I_1$ and $I_3$, matching $I_1I_5III\,III$ and $I_3I_3III\,III$
respectively). So:

> **Theorem H4.** $\sqrt{\mathrm{AZ}(11,5,125)}$ is the period of the rational
> elliptic surface with fibres $I_1\,I_5\,III\,III$ (Herfurtner #38,
> $\mathcal J=-\tfrac1{27}(X^2+6XY+4Y^2)^3/[Y^5(2X+11Y)]$), whose projective
> monodromy group is the index-six genus-zero group with two cusps of widths
> $1,5$ and two elliptic points of order two — i.e. $\Gamma_0(5)$;
> $\sqrt{\mathrm{AZ}(9,3,-27)}$ is the period of $I_3\,I_3\,III\,III$
> (#40, $\mathcal J=(X^2+6XY-3Y^2)^3/(X^2-6XY-3Y^2)^3$), index six with two
> cusps of width $3$. The other seven $\operatorname{Sym}^1$ square roots are
> **not** Picard–Fuchs systems of an elliptic surface over their own base: their
> projective monodromy is an Atkin–Lehner quotient $\Gamma_0(N)+W$, commensurable
> with but not contained in $\mathrm{PSL}_2(\mathbf Z)$.

*Why this is consistent with the modular descriptions in `SPORADIC_SCAN2.md`.*
Apéry's parameter $t=(\eta_1\eta_6/\eta_2\eta_3)^{12}$ has degree two on
$X_0(6)$ and is $W_6$-invariant, so $\mathbf P^1_t=X_0(6)/W_6$; the universal
elliptic curve does not descend along that quotient (only its $\operatorname{
Sym}^2$ does, which is why *Apéry's own order-three row* is fine on
$\mathbf P^1_t$ while its square root is not), and $j$ is not $W_6$-invariant,
so no $\mathcal J$-map of finite degree exists. The two fixed points of $W_6$
become the two order-two points of the quotient, i.e. the two $III$-type
exponents $\rho=\tfrac12$ — this is exactly the "fold at the Fricke point
$i/\sqrt6$" of Theorem B. Same story for $T$ ($W_8$ on $X_0(8)$), Domb, and
$\mathrm{AZ}(7,3,81)$.

*Warning.* $(h,\gamma)$ need not be unique: for a $t$ that is simultaneously a
Hauptmodul for two commensurable subgroups of $\mathrm{PSL}_2(\mathbf Z)$
several pairs fit. The invariant statement is the pair
(*is it in $\mathrm{PSL}_2(\mathbf Z)$?*, *the minimal $\deg\mathcal J$*); $h$
is reported as the smallest width found and is used only to break the
$\mathcal I$-ties of §3.1, where the tie-break is by the existence of a cusp of
that width in the configuration.

---

## 5. The scan

`02_hscan.c` scans, class by class, the three-parameter family
$$P(n)=A\Bigl(n^2+\frac{2M-j_1-j_2}{2M}n\Bigr)+B,\qquad
Q(n)=C\,(Mn-j_1)(Mn-j_2),$$
over $1\le A\le3000$ (with the divisibility $2M\mid A(2M-j_1-j_2)$),
$1\le|B|\le150$, $1\le|C|\le\min(4000,12000/M^2)$, plus a separate pass at
$A=0$. Integrality is tested division-free on $U_n=u_n(n!)^2$,
$$U_{n+1}=P(n)U_n-n^2Q(n)U_{n-1},\qquad u_n\in\mathbf Z\iff v_p(U_n)\ge2v_p(n!),$$
prime by prime modulo $p^K$ with $K$ large enough that no precision is lost, for
every $p\le N=30$. Survivors are re-verified with exact rational arithmetic to
$n=300$ in `08_report.py`, which also computes the companion $b_n$
($b_0=0,b_1=1$), the sharp denominator exponent $k$ (least $k$ with
$d_n^kb_n\in\mathbf Z$, $n\le60$), the characteristic roots, the score
$\log(1/|\lambda_2|)-k$, and the cross-ratio match.

The 26 classes scanned are all $(M;j_1,j_2)$ with $r_i=j_i/M\in[-1,2]$ that are
Kodaira-admissible in the sense of Theorem H2:

```
(1;0,0)  (2;1,1)  (3;1,2)  (4;1,3)  (6;1,5)  (3;1,1)  (3;2,2)  (1;0,1)  (1;1,1)
(2;1,3)  (2;-1,1) (1;-1,1) (3;0,2)  (3;1,3)  (3;-1,1) (4;-1,1) (6;1,3)  (6;3,5)
(6;-1,1) (12;1,7) (12;5,11)(2;3,3)  (1;1,2)  (3;2,4)  (4;3,5)  (6;5,7)
```

Two filters are applied before a row is called interesting.

* **Casoratian degeneracy.** $a_nb_{n+1}-a_{n+1}b_n=\pm\prod_{m\le n}Q(m)/((m+1)!)^2$,
  so if $Q(n_0)=0$ for some integer $n_0\ge1$ — i.e. if $M\mid j_i$ with
  $j_i/M\ge1$ — the Casoratian vanishes identically from $n_0$ on, $b_n$ becomes
  a rational multiple of $a_n$, the Apéry limit is rational and the extension
  class is trivial. Such classes ($(1;0,1)$, $(1;1,1)$, $(1;1,2)$, $(1;-1,1)$,
  $(3;0,2)$, $(3;1,3)$, $(3;2,4)$, $(2;1,3)$ at $n=1$ or higher …) are recorded
  but carry no information.
* **Double characteristic root.** $a^2=4d$ gives $\lambda_1=\lambda_2$ and no
  second growth rate; the row cannot separate a limit.
