# Four-term Apéry-like rows: rank two with five singular points

*Claude (Opus 5), 2026-08-23. Scripts and data: `lattice/four_term/`.
Sequel to `HERFURTNER_CLASSIFICATION.md` (the three-term / four-singular-point
case), `HERFURTNER_K3_WINDOW.md` (the Riemann–Hurwitz bound on $\deg\mathcal J$),
`CUSP_MOVE_PROGRAM.md` (the move theorem) and `NONCONGRUENCE_SCAN.md`
(Theorems N1–N3, the Padé families).*

---

## 0. Verdict first

The object is the integral four-term row
$$(n+1)^2u_{n+1}=P(n)u_n-Q(n)u_{n-1}+R(n)u_{n-2},\qquad u_0=1,\ u_{-1}=u_{-2}=0,$$
$P=an^2+bn+c$, $Q=dn^2+en+f$, $R=gn^2+hn+j$ integral, $g\ne0$. Its
Picard–Fuchs operator is **rank two with five singular points**
$0,t_1,t_2,t_3,\infty$ — the first non-rigid case. Seven things are established.

* **Theorem F1 (exponent dictionary; proved, §1).** With
  $\theta=t\,d/dt$,
  $$L=\theta^2-t\,P(\theta)+t^2Q(\theta+1)-t^3R(\theta+2)
   = t^2\mathcal R(t)D^2+t\,\mathcal S(t)D+t\,\mathcal V(t),$$
  $$\mathcal R=1-at+dt^2-gt^3,\quad
    \mathcal S=1-(a+b)t+(3d+e)t^2-(5g+h)t^3,\quad
    \mathcal V=-c+(d+e+f)t-(4g+2h+j)t^2 .$$
  The exponents are $(0,0)$ at $t=0$ (MUM), $(0,\rho_i)$ at the roots $t_i$ of
  $\mathcal R$ with
  $$\boxed{\ \rho_i=-\frac{T(t_i)}{t_i\,\mathcal R'(t_i)},\qquad
    T:=\mathcal S-t\mathcal R'=1-bt+(d+e)t^2-(2g+h)t^3\ }$$
  and $(2-s_1,\,2-s_2)$ at $t=\infty$, where $s_1,s_2$ are the roots of $R$.
  Fuchs reads $\rho_1+\rho_2+\rho_3=s_1+s_2-1$.

* **Theorem F2 (normalisation classes; proved, §2).**
  $\rho_1=\rho_2=\rho_3=\rho$ **iff**
  $$\boxed{\ b=(1-\rho)a,\qquad e=-2\rho\,d,\qquad h=-(1+3\rho)\,g\ }$$
  — three *linear* conditions, the exact analogue of Theorem H1. If the cubic
  $\mathcal R$ is irreducible over $\mathbf Q$ this is automatic. The class is
  then $(\rho;\delta_\infty)$ with $\delta_\infty=s_2-s_1$, and the remaining
  freedom is $a,d,g$ (the positions, two Möbius moduli $\mathcal I_1=a^3/g$ and
  $\mathcal I_2=a^2/d$, plus the scaling $t\mapsto t/\mu$) together with the
  **two accessory parameters** $c=u_1=P(0)$ and $f=Q(0)$ — exactly $5-3=2$, as a
  rank-two Fuchsian equation with five singular points must have.

* **The semistable classes (§2.3).** "All exponent differences $0$" means
  $\rho\in\mathbf Z$ and $\delta_\infty\in\mathbf Z$. Two presentations matter:
  $$\textbf{(S}_0\textbf{)}\quad P=a(n^2+n)+c,\qquad Q=dn^2+f,\qquad R=C(2n-1)^2$$
  (four fibres $I_n$ and one $I_m^*$: exponents $(\tfrac32,\tfrac32)$ at $\infty$),
  and
  $$\textbf{(S}_1\textbf{)}\quad P=an^2+c,\qquad Q=d(n^2-2n)+f,\qquad R=C(n-2)^2$$
  (all five fibres $I_n$: exponents $(0,0)$ at $\infty$), whose gauge partner is
  $$\textbf{(S}_{-1}\textbf{)}\quad P=a(n^2+2n)+c,\ \ Q=d(n^2+2n)+f,\ \ R=C(n+1)^2 .$$
  $(\mathbf S_0)$ is the true analogue of Zagier's $(an^2+an+b,\;cn^2)$ and it is
  the one that carries integral rows. $(\mathbf S_{\pm1})$ does not: see
  Proposition F5 and §5.

* **Theorem F3 (Riemann–Hurwitz; proved, §2.5).** If the local system is that of
  an elliptic surface over $\mathbf P^1_t$ whose singular fibres sit over the
  five singular points, then
  $$\deg\mathcal J\ \le\ 6c+4e_3+3e_2-12\ \le\ 18,$$
  $c,e_2,e_3$ counting the points whose projective monodromy is parabolic, of
  order $2$, of order $3$. The $\mathcal J$-map test is run in that window.

* **The deduplication is the whole story (§3, §5).** Of the $2178$ distinct
  integral rows the scan produces, $2138$ — **98.2 %** — have an *apparent*
  singularity: a point where the exponent difference is a positive integer and
  the Frobenius obstruction vanishes, so the **projective** local system has only
  **four** singular points and the row is a gauge / cusp-move image of a
  three-term row.
  Proposition F4 writes the largest such family down in closed form. What is left
  is $30$ genuinely five-point rows (plus $6$ rescalings of them); $28$ of them
  have $k=2$; and no Apéry limit among them is a classical weight-two constant.

* **One of them is an elliptic $K3$ (§5.3).** The row
  $$(n+1)^2u_{n+1}=(11n^2+11n+4)u_n-(37n^2+3)u_{n-1}+3(3n-1)(3n-2)u_{n-2},$$
  $u_n=1,4,16,64,250,928,3136,8704,11866,-79400,\dots$, passes the
  $\mathcal J$-map test with $h=4$, $\deg\mathcal J=16$, and the certified
  $\mathcal J$ gives the fibre configuration
  $$I_4\ (t=0),\quad I_6\ (t=1),\quad I_3,I_3\ \bigl(t=\tfrac{5\pm i\sqrt2}{27}\bigr),
    \quad IV^*\ (t=\infty),\qquad \textstyle\sum e=24 .$$
  Because $\sum_{I_n}n_i=16>12$ this is **not** a rational elliptic surface: it is
  an elliptic $K3$ with exactly five singular fibres. Proposition K1 of
  `HERFURTNER_K3_WINDOW.md` showed that a *three-term* row can never carry one
  ($\deg\mathcal J\le12$ there, always). **The four-term row is where the
  elliptic $K3$s enter the Apéry world.**

* **Corollary F6 (no new irrationality candidate; verified, §6.1).** A positive
  score $\log(1/|\lambda_2|)-k>0$ needs $|\lambda_2|<e^{-1}$, which forces
  $|d|\le0.74(|a|+1)+2$ and $|g|\le0.136(|a|+1)+1$: a *provably complete* window,
  swept for $a\le800$, $|c|\le60$, $|f|\le120$ over all $27$ classes. It contains
  exactly two rows, both with $k=1$ and both carrying an $I_0^*$ at $\infty$
  (exponents $\tfrac12,\tfrac52$, no logarithm) — so their projective local
  system has four points and they are the four-term avatars of the classical
  Legendre–Padé family. **No genuine five-singular-point row has a positive
  score.**

---

## 1. The exponent dictionary

Let $\theta=t\,d/dt$, $D=d/dt$, and write the row as
$$(n+1)^2u_{n+1}=P(n)u_n-Q(n)u_{n-1}+R(n)u_{n-2},\qquad u_0=1,\ u_{-1}=u_{-2}=0,$$
$$P=an^2+bn+c,\qquad Q=dn^2+en+f,\qquad R=gn^2+hn+j .$$

**Lemma 1.1 (the operator).** $y=\sum_{n\ge0}u_nt^n$ satisfies $Ly=0$ with
$$L=\theta^2-t\,P(\theta)+t^2Q(\theta+1)-t^3R(\theta+2).$$
*Proof.* $[t^{n+1}]\,\theta^2y=(n+1)^2u_{n+1}$; $[t^{n+1}]\,tP(\theta)y=P(n)u_n$;
$[t^{n+1}]\,t^2Q(\theta+1)y=Q(n)u_{n-1}$; $[t^{n+1}]\,t^3R(\theta+2)y=R(n)u_{n-2}$.
$\square$ (Verified symbolically to $n=12$ in `01_exponents.py`.)

Expanding $\theta=tD$, $\theta^2=t^2D^2+tD$:
$$\boxed{\ L=t^2\mathcal R(t)\,D^2+t\,\mathcal S(t)\,D+t\,\mathcal V(t)\ }$$
$$\mathcal R(t)=1-at+dt^2-gt^3,\qquad
  \mathcal S(t)=1-(a+b)t+(3d+e)t^2-(5g+h)t^3,$$
$$\mathcal V(t)=-c+(d+e+f)t-(4g+2h+j)t^2 .$$

**Singular points.** $t=0$; the three roots $t_1,t_2,t_3$ of $\mathcal R$; and
$t=\infty$. Since $t^3(\lambda^3-a\lambda^2+d\lambda-g)=\mathcal R(1/\lambda)\,$
with $\lambda=1/t$, the **characteristic roots** $\lambda_i=1/t_i$ are the roots
of
$$\lambda^3-a\lambda^2+d\lambda-g=0 .$$
Five distinct singular points requires $g\ne0$ and
$$\Delta:=18adg-4a^3g+a^2d^2-4d^3-27g^2\ \ne\ 0 .$$

**Exponents.**

* At $t=0$: $\theta^2$ gives $(0,0)$ with a logarithm — a MUM point, hence a
  fibre of type $I_{n_0}$, $n_0\ge1$.
* At $t=t_i$ (a simple zero of $\mathcal R$): $0$ is always an exponent, and
  with $T:=\mathcal S-t\mathcal R'$,
  $$T(t)=1-bt+(d+e)t^2-(2g+h)t^3,\qquad
    \rho_i=1-\frac{\mathcal S(t_i)}{t_i\mathcal R'(t_i)}=-\frac{T(t_i)}{t_i\mathcal R'(t_i)} .$$
* At $t=\infty$: $Lt^{\mu}$ has leading term $-R(\mu+2)t^{\mu+3}$, so
  $\mu+2\in\{s_1,s_2\}$, the roots of $R$; with the convention
  (exponent) $=-\mu$ the exponents at $\infty$ are
  $$2-s_1,\qquad 2-s_2,\qquad \delta_\infty:=s_2-s_1 .$$

**Proposition 1.2 (Fuchs).** $\rho_1+\rho_2+\rho_3=s_1+s_2-1=-h/g-1$.

*Proof.* $\sum_i\rho_i=-\sum_i\operatorname{Res}_{t_i}\bigl(T/(t\mathcal R)\bigr)$.
The rational function $T/(t\mathcal R)$ has residue $T(0)/\mathcal R(0)=1$ at
$t=0$ and behaves like $\bigl((2g+h)/g\bigr)t^{-1}$ at $\infty$, so the sum of
its finite residues is $(2g+h)/g$; hence $\sum\rho_i=1-(2g+h)/g=-1-h/g$. The sum
of all ten exponents is $0+0+\sum\rho_i+(2-s_1)+(2-s_2)=3=(5-2)\cdot1$, the Fuchs
relation for a rank-two operator with five singular points. $\square$

(Both identities are checked symbolically in `01_exponents.py`.)

**Möbius invariants.** The scaling $t\mapsto t/\mu$, i.e. $u_n\mapsto\mu^nu_n$,
acts by $(a,b,c)\mapsto\mu(a,b,c)$, $(d,e,f)\mapsto\mu^2(d,e,f)$,
$(g,h,j)\mapsto\mu^3(g,h,j)$. The labelled configuration
$(0,\infty;\{t_1,t_2,t_3\})$ therefore has exactly two invariants,
$$\mathcal I_1=\frac{a^3}{g}=\frac{(\lambda_1+\lambda_2+\lambda_3)^3}{\lambda_1\lambda_2\lambda_3},
\qquad
\mathcal I_2=\frac{a^2}{d}=\frac{(\lambda_1+\lambda_2+\lambda_3)^2}{\sum_{i<j}\lambda_i\lambda_j},$$
the analogue of $\mathcal I=a^2/d$ of Theorem H3. Both are rational for an
integral row.

> **Why the cross-ratio test does not transfer.** In the four-point case
> Herfurtner's rigid configurations pin $\mathcal I$ to one of finitely many
> algebraic numbers, and rationality kills most of them. Here the fibre
> configurations with five singular fibres come in **one-parameter families**
> (Persson), so the locus of realisable $(\mathcal I_1,\mathcal I_2)$ is a
> *curve* in the $(\mathcal I_1,\mathcal I_2)$-plane for each fibre type, not a
> point. Rationality is no longer restrictive, and the finiteness mechanism of
> Theorem H3 disappears. This is the precise sense in which the four-term case is
> "non-rigid", and it is why the $\mathcal J$-map test (§6) has to do all the
> work.

---

## 2. Normalisation classes

### 2.1 Equal exponents

**Theorem F2.** *Assume $t_1,t_2,t_3$ are distinct. Then
$\rho_1=\rho_2=\rho_3=\rho$ if and only if*
$$b=(1-\rho)a,\qquad e=-2\rho d,\qquad h=-(1+3\rho)g .$$

*Proof.* $\rho_i=\rho$ for all $i$ says $\mathcal R\mid T+\rho\,t\mathcal R'$.
Both are cubics with constant term $1$, so $T+\rho t\mathcal R'=\mathcal R$:
$$1-(b+\rho a)t+(d+e+2\rho d)t^2-(2g+h+3\rho g)t^3=1-at+dt^2-gt^3 .\qquad\square$$

**Corollary F2.1 (Galois).** *The $\rho_i$ are values at $t_i$ of a fixed
rational function with rational coefficients, hence are permuted by
$\operatorname{Gal}(\overline{\mathbf Q}/\mathbf Q)$ exactly as the $t_i$ are.
Kodaira admissibility (F2.2) forces $\rho_i\in\mathbf Q$; therefore $\rho_i$ is
constant on each Galois orbit of $\{t_1,t_2,t_3\}$. In particular if
$\lambda^3-a\lambda^2+d\lambda-g$ is irreducible then $\rho_1=\rho_2=\rho_3$.*

This is the exact analogue of Theorem N3 / Corollary H2.1: a configuration whose
finite fibres are of different Kodaira classes needs the cubic to be reducible.

**Corollary F2.2 (Kodaira admissibility).** *If the local system of $L$ is, up to
a rational gauge factor and a quadratic twist, $R^1\pi_*\mathbf Q$ of an elliptic
surface over $\mathbf P^1_t$, then*
$$\rho_1,\rho_2,\rho_3,\ \delta_\infty\ \in\ \{0,\tfrac12,\tfrac13,\tfrac23\}
\pmod 1 .$$
The proof is Theorem H2 verbatim: neither a rational gauge nor a quadratic twist
changes an exponent *difference* mod $1$, and Kodaira's monodromy eigenvalue
ratios are $1$ ($I_n,I_n^*$), $-1$ ($III,III^*$), $e^{\pm2\pi i/3}$
($IV,IV^*$), $e^{\mp2\pi i/3}$ ($II,II^*$).

### 2.2 The class parametrisation used by the scanner

Fix $\rho$ and write $R(n)=C(Mn-j_1)(Mn-j_2)$, so $s_i=j_i/M$,
$$g=CM^2,\qquad h=-CM(j_1+j_2),\qquad j=Cj_1j_2,\qquad j_1+j_2=(1+3\rho)M$$
(the last is Proposition 1.2 combined with Theorem F2). A **class** is the tuple
$(\rho;M,j_1,j_2)$, and inside it the scan runs over
$$\boxed{\ P(n)=a\bigl(n^2+(1-\rho)n\bigr)+c,\quad
  Q(n)=d\bigl(n^2-2\rho n\bigr)+f,\quad
  R(n)=C(Mn-j_1)(Mn-j_2)\ }$$
with five integer parameters $a,c,d,f,C$ (four after $t\mapsto t/\mu$).
`04_classes.py` enumerates the Kodaira-admissible classes; there are $156$ with
$\rho,\delta_\infty$ of denominator $\le6$, $\rho\in[-1,2]$,
$\delta_\infty\in[0,4]$ and $M\le12$.

### 2.3 The semistable classes, explicitly

"All five exponent differences $\equiv0$" means $\rho\in\mathbf Z$,
$\delta_\infty\in\mathbf Z$. Since $s_1+s_2=1+3\rho$ and $\delta_\infty=0$ forces
$s_1=s_2=(1+3\rho)/2$, the exponents at $\infty$ are integral iff $\rho$ is odd.

| class | $\rho$ | $\delta_\infty$ | $P$ | $Q$ | $R$ | exps at $\infty$ | fibres |
|---|---|---|---|---|---|---|---|
| $\mathbf S_0$ | $0$ | $0$ | $a(n^2+n)+c$ | $dn^2+f$ | $C(2n-1)^2$ | $(\tfrac32,\tfrac32)$ | $I\,I\,I\,I\;I_m^*$ |
| $\mathbf S_1$ | $1$ | $0$ | $an^2+c$ | $d(n^2-2n)+f$ | $C(n-2)^2$ | $(0,0)$ | $I\,I\,I\,I\,I$ |
| $\mathbf S_{-1}$ | $-1$ | $0$ | $a(n^2+2n)+c$ | $d(n^2+2n)+f$ | $C(n+1)^2$ | $(3,3)$ | $I\,I\,I\,I\,I$ |
| $\mathbf S_0'$ | $0$ | $1$ | $a(n^2+n)+c$ | $dn^2+f$ | $Cn(n-1)$ | $(2,1)$ | apparent at $\infty$ |

$\mathbf S_1$ and $\mathbf S_{-1}$ are the two gauge presentations of the same
"five $I_n$ fibres" local system: the gauge $y\mapsto\prod_i(1-\lambda_it)^{-\rho}y$
sends $\rho\mapsto-\rho$ and $s\mapsto s-3\rho$.

The distinction between $\mathbf S_0$ and $\mathbf S_{\pm1}$ is not a matter of
gauge. A quadratic character twist flips the "starredness" of the fibres at an
**even** number of points, so the parity of the number of $I^*$ fibres is an
invariant of the projective local system: it is $1$ for $\mathbf S_0$ (the
exponents at $\infty$ are half-integral, the others integral) and $0$ for
$\mathbf S_{\pm1}$. $\mathbf S_0$ is therefore "four $I_n$ and one $I_m^*$" up to
twist, with $\sum e=\sum_{i=0}^3 n_i+m+6$; $\mathbf S_{\pm1}$ is "five $I_n$",
with $\sum e=\sum_{i=0}^4 n_i$.

$(\mathbf S_0)$ is the analogue of Zagier's $(an^2+an+b,\;cn^2)$: the *linear*
conditions are $b=a$, $e=0$, $h=-g$, $j=g/4$, and the two accessory parameters
$c=u_1$ and $f=Q(0)$ are free. Note that, unlike the three-term case, $f$ is
**not** forced to vanish — the fibre at $\infty$ is governed by $R$, not by $Q$.

For a *rational* elliptic surface ($\sum e=12$) this gives
$\sum_{i=0}^{3}n_i+m=6$ in $\mathbf S_0$ and $\sum_{i=0}^4n_i=12$ in
$\mathbf S_{\pm1}$ — the latter being Beauville's condition with five instead of
four fibres, and those configurations form one-parameter families (Persson),
which is the geometric source of the non-rigidity.

### 2.4 Apparent singularities

At $t_i$ the exponents are $(0,\rho_i)$. If $\rho_i\in\mathbf Z\setminus\{0\}$
the two exponents differ by an integer and the second solution either carries a
logarithm (a genuine $I_n$ fibre) or does not, in which case the local monodromy
is **trivial**: an *apparent* singularity. The same at $\infty$ when
$\delta_\infty\in\mathbf Z_{>0}$. If $\rho_i=0$ (or $\delta_\infty=0$) the
exponents coincide and a logarithm is unavoidable: such a point is always a
genuine $I_n$.

The obstruction is one Frobenius coefficient. Writing the ODE at the point as
$A(x)y''+B(x)y'+C(x)y=0$ and running the recursion from the smaller exponent
$\nu$, the coefficient $c_m$ at $m=|\rho_i|$ (resp. $m=\delta_\infty$) has
vanishing indicial factor; the point is apparent iff the accumulated right-hand
side vanishes. `05_report.py` evaluates this at $100$ digits.

Two cases must be distinguished. If the two exponents are integers the local
monodromy is then the **identity** and the point is not a singular point of the
local system at all. If they are half-integers (which happens at $\infty$ in the
classes with $\delta_\infty\in\mathbf Z$ but $s_i\notin\mathbf Z$) the monodromy
is $-I$: trivial in $\mathrm{PSL}_2$, i.e. an $I_0^*$ fibre. Both are invisible
to the **projective** local system, and it is the projective one that the
$\mathcal J$-map test and the Kodaira dictionary see; so in both cases the row is
counted as a four-point object. This is exactly the $I_0^*$ mechanism that put
the class $(2;1,3)$ into Herfurtner's one-parameter block in Theorem H5(iv).

**Proposition F5 (worked case $\mathbf S_{-1}$).** *In the class $(\rho;\delta_\infty)=(-1;0)$
the obstruction at $t_i$ is*
$$S_i \;=\; t_i\bigl[(a-c)+(f-d)\,t_i\bigr].$$
*Hence at most one $t_i$ is apparent unless $a=c$ and $d=f$, and in that case
all three are, which happens exactly when*
$$P(n)=a(n+1)^2,\qquad Q(n)=d(n+1)^2,\qquad R(n)=C(n+1)^2,$$
*i.e. exactly when the row degenerates to the constant-coefficient recurrence
$u_{n+1}=au_n-du_{n-1}+Cu_{n-2}$ and $\sum u_nt^n$ is a rational function.*

*Proof.* Substituting $b=2a$, $e=2d$, $h=2g$, $j=g$ into
$S_i=2A_2-B_1+C_0$ with $A=t^2\mathcal R$, $B=t\mathcal S$, $C=t\mathcal V$
expanded at $t_i$, and using $\mathcal R(t_i)=0$ and
$\mathcal S(t_i)=2t_i\mathcal R'(t_i)$ (which is $\rho_i=-1$), gives
$S_i=t_i\bigl[2\mathcal R'+t_i\mathcal R''-\mathcal S'+\mathcal V\bigr](t_i)
=t_i[(a-c)+(f-d)t_i]$ identically in $t_i$. $\square$

This is the four-term analogue of the "$Q(n_0)=0$" degeneracy: it is detected by
the scanner directly (`03_fscan.c`) and is by far the largest single family in
the naive output.

### 2.5 Theorem F3: $\deg\mathcal J\le18$

*Let $\pi:\mathcal E\to\mathbf P^1_t$ be a minimal elliptic surface with section
and non-constant $\mathcal J$ whose singular fibres lie over the five singular
points of $L$, and let $c,e_2,e_3$ count the points whose projective local
monodromy is parabolic, of order $2$, of order $3$ ($c+e_2+e_3\le5$, $c\ge1$
because $t=0$ is MUM). Then $\mu=\deg\mathcal J\le6c+4e_3+3e_2-12\le18$.*

*Proof.* Identical to Proposition K1 of `HERFURTNER_K3_WINDOW.md`; the argument
never used the number of singular points. Riemann–Hurwitz for
$\mathcal J:\mathbf P^1\to\mathbf P^1$ and Kodaira's dictionary give
$2\mu-2\ge(\mu-c)+\tfrac23(\mu-e_3)+\tfrac12(\mu-e_2)$. Maximising subject to
$c+e_2+e_3\le5$ gives $\mu\le6\cdot5-12=18$, attained only in the all-parabolic
case. $\square$

Accordingly `06_jtest.gp` runs the $\mathcal J$-map test with
$\mathrm{DMAX}=\mathrm{HMAX}=18$, $130$ series terms and $760$ values of
$\gamma$ per $h$ ($\pm m^{\pm h}$, $m$ running over the $7$-smooth integers $\le2048$).

---

## 3. Degeneracies and deduplication

A hit of the scan is discarded, or reclassified, by the following tests, in
order.

**(D1) Repeated characteristic root.** $\Delta=0$: two of the $t_i$ collide and
the point is no longer a regular singular point of the normal form. Applied
inside `03_fscan.c`.

**(D2) Terminating rows.** If $u_1=u_2=u_3=0$ then $u_n=0$ for all $n\ge1$; more
generally three consecutive zeros kill the sequence. Applied inside the scanner
(`nontrivial()`); in the class $\mathbf S_1$ this alone removes **all** of the
$6.5\times10^5$ raw hits at $|a|\le40$.

**(D3) Constant-coefficient rows.** $P,Q,R$ all proportional to $(n+1)^2$
(Proposition F5). Applied inside the scanner.

**(D4) Casoratian degeneracy.** For three solutions of the row the Casoratian
$W_n=\det\bigl[u^{(i)}_{n-2+k}\bigr]$ satisfies
$$W_{n+1}=\frac{R(n)}{(n+1)^2}\,W_n ,$$
(expand the determinant along the last row and use the recurrence; checked
numerically on the $K3$ row of §5.3 for $2\le n\le8$), so if $R(n_0)=0$ for an
integer $n_0\ge1$ — i.e. if $M\mid j_i$ and $j_i/M\ge1$ — then $W\equiv0$ from
$n_0$ on and the row's three-dimensional solution space collapses to two
dimensions. These classes ($\mathbf S_1$, $(1;1,1,3)$,
$(\tfrac23;1,1,2)$, $(0;1,0,1)$, …) are recorded but carry no new information.

**(D5) Apparent singularity $\Rightarrow$ disguised three-term row.** If the
Frobenius obstruction vanishes at some $t_i$ or at $\infty$, the **projective**
local system has only **four** singular points (§2.4), and the row is a gauge /
cusp-move image of a three-term row — over $\mathbf Q$ when the characteristic
root to be sent to $\infty$ is rational, and over $\mathbf Q(\lambda_i)$ in
general, since the Möbius map $t\mapsto t/(1-\lambda_it)$ that makes the
presentation three-term is defined over that field. Two mechanisms produce these,
and between them they account for essentially all such hits:

* *the generic cusp move.* Theorem 1 of `CUSP_MOVE_PROGRAM.md` moves a three-term
  row by $s=t/(1-\nu t)$ **with $\nu$ a characteristic root**, and the result is
  again three-term. If $\nu$ is *not* a characteristic root, the ordinary point
  $t=1/\nu$ is sent to $s=\infty$ and $t=\infty$ to $s=-1/\nu$: the new operator
  has the five singular points $0,\ s(t_1),\ s(t_2),\ -1/\nu,\ \infty$, of which
  $s=\infty$ is **apparent**. Equivalently, on the generalised signed binomial
  transform
  $$u_n^{\sharp}=\sum_{m\le n}\binom{n}{n-m}(-\nu)^{n-m}u_m ,$$
  a three-term row becomes a four-term row in a class with
  $\delta_\infty\in\mathbf Z$; Proposition F4 makes this explicit.
  *Worked example.* Zagier $\mathbf A$ ($P=7n^2+7n+2$, $Q=-8n^2$, characteristic
  roots $\{8,-1\}$) with $\nu=1$ gives
  $$(n+1)^2u_{n+1}=(2n+1)^2u_n+19n^2u_{n-1}+14n(n-1)u_{n-2},\qquad
    u_n=1,1,7,31,175,991,\dots$$
  in the class $(\rho;\delta_\infty)=(0;1)$, i.e. $(M;j_1,j_2)=(1;0,1)$,
  $(a,c,d,f,C)=(4,1,-19,0,14)$, with characteristic roots
  $\{7,-2,-1\}=\{8,-1,0\}-1$ (the translation by $-\nu$ in the
  coordinate $v=1/t$, exactly Theorem 3 of `CUSP_MOVE_PROGRAM.md`). Its
  Frobenius obstruction at $\infty$ vanishes, and the $\mathcal J$-map test
  returns $(h,\deg\mathcal J,\gamma)=(1,12,1)$ — Zagier $\mathbf A$'s own
  invariants, as it must, the local system being unchanged.
* *a rational gauge at an ordinary point.* $y\mapsto y/(1-\nu t)$ with
  $\mathcal R(1/\nu)\ne0$ creates a fifth singular point at $t=1/\nu$ with
  exponents $(0,-1)$ and no logarithm. On Zagier $\mathbf A$ with $\nu=1$:
  $$(n+1)^2u_{n+1}=(8n^2+9n+3)u_n-(-n^2+7n+2)u_{n-1}-8n^2u_{n-2},$$
  $u_n=1,3,13,69,415,2667,\dots$, characteristic roots $\{8,1,-1\}$, with
  $\rho=(0,0,-1)$ — a class with *unequal* $\rho_i$, so the cubic must be (and
  is) reducible, in accordance with Corollary F2.1.

**(D6) Rescaling.** $u_n\mapsto\mu^nu_n$ acts by
$(a,c,d,f,C)\mapsto(\mu a,\mu c,\mu^2d,\mu^2f,\mu^3C)$; a hit is marked
`RESCALED` when dividing by some $\mu\ge2$ leaves an integral row. $\mu=-1$ is
$t\mapsto-t$ and is absorbed by normalising $a>0$.

**(D7) Pullbacks $t\mapsto t^2$.** These would require $\mathcal R(-t)=\mathcal R(t)$,
i.e. $a=g=0$; but $g\ne0$ by hypothesis. **There are no $t\mapsto t^2$ pullbacks
in the four-term world** — the odd-coefficient phenomenon of the three-term case
cannot occur with an odd-degree $\mathcal R$.

**(D8) The Padé $\log$/$\arctan$ families** of `NONCONGRUENCE_SCAN.md` §4 are
three-term, with $k=1$ resp. $k=2$; their four-term images are removed by (D5).
Exactly two of the surviving five-point rows have $k=1$ — the diagnostic of
`NONCONGRUENCE_SCAN.md` §4.2 marks them as Padé-type; both have negative score
and neither belongs to a family the scan could detect (§5.5 note 1).

---
## 4. The scan

### 4.1 Design

`03_fscan.c` scans one class $(\rho;M,j_1,j_2)$ at a time over the five integer
parameters $a,c,d,f,C$. Integrality is tested division-free on
$U_n=u_n\,(n!)^2$:
$$U_{n+1}=P(n)\,U_n-n^2Q(n)\,U_{n-1}+n^2(n-1)^2R(n)\,U_{n-2},\qquad
  u_n\in\mathbf Z\iff v_p(U_n)\ge2v_p(n!)\ \ \forall p,$$
prime by prime modulo $p^K$ with $K=2v_p(N!)+2$, $N=24$. Three cheap layers
precede the exact test.

1. **Analytic congruences.** $U_2=P(1)c-Q(1)$ and $U_3=K_3+A_3C$, $U_4=K_4+A_4C$
   are *linear in $C$*, with
   $$A_3=4(2M-j_1)(2M-j_2),\qquad A_4=P(3)A_3+36\,c\,(3M-j_1)(3M-j_2).$$
   The conditions $4\mid U_2$, $36\mid U_3$, $576\mid U_4$ (i.e.
   $v_2\ge2v_2(n!)$, $v_3\ge2v_3(n!)$ for $n=2,3,4$) are solved for $C$ by two
   linear congruences and a CRT, so the innermost loop runs over an arithmetic
   progression instead of an interval — typically a $100$-fold saving.
2. **A mask-arithmetic $2$-adic prefilter** to $n=14$ modulo $2^{24}$ (no
   divisions; the valuation test is a `ctz`).
3. **A 64-bit prefilter** for $p=3,5,7,11$ to $n=12$.

Degeneracies (D1)–(D3) of §3 are applied inside the scanner. Survivors are
re-verified with **exact rational arithmetic to $n=120$** in `05_report.py`,
which also computes the characteristic roots, the Frobenius obstruction at every
point with an integral exponent difference, the companion $b_n$ ($b_0=0$,
$b_1=1$) and its sharp denominator exponent $k$ (least $k$ with
$d_n^kb_n\in\mathbf Z$ for $n\le60$), the Apéry limit $\xi=\lim b_n/a_n$ to
$\ge60$ digits (adaptive length, $n$ up to $2\times10^4$), and the score
$\log(1/|\lambda_2|)-k$ with $\lambda_2$ the second-largest characteristic root
in modulus.

### 4.2 Three boxes

**Box G (general census).** $1\le a\le200$ ($\le140$ for the classes scanned in
the second pass), $|c|\le40$, $|d|\le1200$, $|f|\le80$, $1\le|C|\le120$,
integrality prime-tested to $n=24$ and re-verified exactly to $n=120$.

**Box W (the $|\lambda_2|<1$ window; complete).** A positive score needs
$\log(1/|\lambda_2|)>k\ge0$, hence $|\lambda_2|<1$ and therefore
$|\lambda_3|\le|\lambda_2|<1$. From $a=\sum\lambda_i$,
$d=\sum_{i<j}\lambda_i\lambda_j$, $g=\lambda_1\lambda_2\lambda_3$ this forces
$$|\lambda_1|\le|a|+2,\qquad |d|\le2|a|+5,\qquad |C|=|g|/M^2\le(|a|+2)/M^2 .$$
Box W therefore takes $|c|\le60$, $|f|\le120$ and lets $|d|,|C|$ run over exactly
those ranges; it was swept to $a\le300$ for the $18$ classes of `classes2.txt`
and to $a\le400$ for the semistable class $\mathbf S_0$. **Inside box W the scan
is complete for the question "is there a row with $|\lambda_2|<1$?"** in those
classes — no such row with $|c|\le60$, $|f|\le120$ and $a$ in range can escape
it.

**Box PW (the positive-score window; complete).** Every row of the census has
$k\ge1$, so a positive score needs $|\lambda_2|<e^{-1}$, whence
$$|\lambda_1|\le|a|+2e^{-1},\qquad |d|\le0.74(|a|+1)+2,\qquad
  |C|\le\bigl(0.136(|a|+1)+1\bigr)/M^2 .$$
This is an order of magnitude smaller than box W and lets $a$ run much further:
box PW takes $1\le a\le1500$, $|c|\le60$, $|f|\le120$, all $27$ classes.
**Inside box PW the scan is complete for the question "is there a row with
positive score?"** among rows with $k\ge1$.

### 4.3 The classes scanned

| $(\rho;M,j_1,j_2)$ | $(\rho;\delta_\infty)$ | types $(0;t_1t_2t_3;\infty)$ | $R(n)$ |
|---|---|---|---|
| $(0;2,1,1)$ | $(0;0)$ | cusp; cusp$^3$; cusp | $C(2n-1)^2$ |
| $(0;3,1,2)$ | $(0;\tfrac13)$ | cusp; cusp$^3$; ord $3$ | $C(3n-1)(3n-2)$ |
| $(0;4,1,3)$ | $(0;\tfrac12)$ | cusp; cusp$^3$; ord $2$ | $C(4n-1)(4n-3)$ |
| $(0;6,1,5)$ | $(0;\tfrac23)$ | cusp; cusp$^3$; ord $3$ | $C(6n-1)(6n-5)$ |
| $(0;3,-1,4)$, $(0;4,-1,5)$, $(0;6,-1,7)$ | $(0;\tfrac53),(0;\tfrac32),(0;\tfrac43)$ | as above | $C(3n+1)(3n-4)$ etc. |
| $(\pm\tfrac12;4,\pm(5,5)\ \text{resp.}\ (-1,-1))$ | $(\pm\tfrac12;0)$ | cusp; ord $2^3$; cusp | $C(4n-5)^2$, $C(4n+1)^2$ |
| $(\tfrac13;1,1,1)$, $(-\tfrac13;1,0,0)$ | $(\pm\tfrac13;0)$ | cusp; ord $3^3$; cusp | $C(n-1)^2$, $Cn^2$ |
| $(\tfrac23;2,3,3)$, $(-\tfrac23;2,-1,-1)$ | $(\pm\tfrac23;0)$ | cusp; ord $3^3$; cusp | $C(2n-3)^2$, $C(2n+1)^2$ |
| $(\pm1;M,\dots)$, $\delta_\infty\in\{\tfrac13,\tfrac12,\tfrac23\}$ | | cusp; cusp$^3$; ord $2$/ord $3$ | $C(6n-11)(6n-13)$ etc. |
| $(0;2,1,1)$, $(\pm1;\dots)$, $(0;1,0,1)$, $(1;1,1,3)$, $(1;2,3,5)$, $(\tfrac13;2,1,3)$, $(\tfrac23;1,1,2)$, $(0;2,-1,3)$ | $\delta_\infty$ or $\rho\in\mathbf Z$ | | the degenerate / disguised classes |

Signature list. By Gauss–Bonnet a genus-zero group $\Gamma\le\mathrm{PSL}_2(\mathbf Z)$
with exactly **five** special points ($c$ cusps, $e_2$ elliptic of order $2$,
$e_3$ of order $3$, $c+e_2+e_3=5$, $c\ge1$) has index
$$\mu=12\Bigl(-1+\frac c2+\frac{e_2}4+\frac{e_3}3\Bigr)
   \;=\;6c+4e_3+3e_2-12,$$
which is exactly the Riemann–Hurwitz bound of Theorem F3 — the two agree because
equality in F3 says precisely that $\mathcal J$ is a Belyi map. The fifteen
possible signatures and indices are

| $(c,e_2,e_3)$ | $(5,0,0)$ | $(4,1,0)$ | $(4,0,1)$ | $(3,2,0)$ | $(3,1,1)$ | $(3,0,2)$ | $(2,3,0)$ | $(2,2,1)$ |
|---|---|---|---|---|---|---|---|---|
| $\mu$ | $18$ | $15$ | $16$ | $12$ | $13$ | $14$ | $9$ | $10$ |

| $(c,e_2,e_3)$ | $(2,1,2)$ | $(2,0,3)$ | $(1,4,0)$ | $(1,3,1)$ | $(1,2,2)$ | $(1,1,3)$ | $(1,0,4)$ |
|---|---|---|---|---|---|---|---|
| $\mu$ | $11$ | $12$ | $6$ | $7$ | $8$ | $9$ | $10$ |

so there are **finitely many host groups**, the exact five-point analogue of
Theorem N2 (four special points, $\mu\le12$, $28$ groups). The enumeration of the
conjugacy classes themselves is not carried out here.

---
## 5. What the scan produces

### 5.1 The shape of the output

Raw hit counts are dominated by degenerate and disguised rows. Two classes make
the point on their own: in $\mathbf S_1$ ($\rho=1$, $R=C(n-2)^2$) the box
$a\le40$, $|c|\le30$, $|d|\le600$, $|f|\le50$, $|C|\le80$ contains
$646\,297$ parameter tuples passing integrality mod $p\le24$ to $n=24$ with
$\Delta\ne0$ — and **every one of them** has $u_n=0$ for $n\ge1$, the trivial
solution created by $R(2)=0$; in $\mathbf S_{-1}$ ($\rho=-1$, $R=C(n+1)^2$) the
same box contains $484\,697$, and **every one of them** is the
constant-coefficient row of Proposition F5. After (D2)–(D3) both classes are
empty.

Over the whole census — $27$ classes, box G for $21$ of them and the probe box
for the six purely disguised ones (§4.3) — the filters act as follows.

| stage | rows |
|---|---|
| integral mod $p\le24$ to $n=24$, $\Delta\ne0$, non-terminating, not constant-coefficient | $2178$ |
| of which exact integrality to $n=120$ fails | $4$ |
| of which **(D5)** an apparent singularity: four-point local system, a disguised three-term row | $2138$ $(98.2\%)$ |
| of which **(D6)** a rescaling $u_n\mapsto\mu^nu_n$ of another integral row | $6$ |
| **genuine five-point rows, primitive** | $\mathbf{30}$ |

Two observations, both stable across all boxes tried.

* **The classes producing $10^5$–$10^6$ raw hits are exactly the
  Casoratian-degenerate ones**, i.e. those in which $R$ has a positive integer
  root ($M\mid j_i$, $j_i/M\ge1$): $\mathbf S_1$ ($R=C(n-2)^2$), $(1;1,1,3)$,
  $(\tfrac23;1,1,2)$, and $\mathbf S_{-1}$ through Proposition F5. In each of
  them all the hits are rows whose sequence terminates ($u_n=0$ from $n=1$ or
  $n=2$ on) or which are constant-coefficient, and after (D2)–(D3) the class is
  empty. The one Casoratian-degenerate class that survives (D2)–(D3) is
  $(0;1,0,1)$, $R=Cn(n-1)$ — the class of Proposition F4 — and all of its rows
  are disguised.
* **$98\%$ of the remaining hits are disguised three-term rows.** Every class with
  $\rho\in\mathbf Z\setminus\{0\}$ or $\delta_\infty\in\mathbf Z_{>0}$ produces a
  large family, and the Frobenius obstruction vanishes for all of them. The class
  $(1;2,3,5)$ is the extreme case: it has *two* apparent points, so the local
  system has only **three** singular points and the rows are hypergeometric —
  many different $(a,c,d,f,C)$ give the same sequence
  $1,-4,-12,-80,-700,-7056,\dots$

The classes with $\rho\notin\mathbf Z$ or ($\rho=0$ and $\delta_\infty\notin\mathbf Z$)
admit **no** apparent singularity at all — when $\rho=0$ the two exponents at
$t_i$ coincide and a logarithm is forced — so every hit there is a genuine
five-point row. Those are the classes that carry the census.

### 5.2 The five-singular-point rows

$28$ of the $30$ surviving rows have $k=2$: the four-term world reproduces the
"free integration" phenomenon of Theorem R3 almost exactly (the two exceptions
have $k=1$ and are Padé-type; §5.5 note 1). All of them are sporadic — the
scan produced **no one-parameter integral family of genuine five-point rows**;
every infinite family it produced is a family of disguised three-term rows or a
Casoratian degeneration.

### 5.3 The one elliptic-surface row: an elliptic $K3$ with five singular fibres

$$\boxed{\ (n+1)^2u_{n+1}=(11n^2+11n+4)\,u_n-(37n^2+3)\,u_{n-1}+3(3n-1)(3n-2)\,u_{n-2}\ }$$
$$u_0=1:\qquad 1,\;4,\;16,\;64,\;250,\;928,\;3136,\;8704,\;11866,\;-79400,\;
-975488,\;\dots$$

Class $(\rho;M,j_1,j_2)=(0;3,1,2)$, i.e. $(\rho;\delta_\infty)=(0;\tfrac13)$:
exponents $(0,0)$ at $t=0$ and at each $t_i$, and $(\tfrac53,\tfrac43)$ at
$\infty$. Characteristic polynomial
$$\lambda^3-11\lambda^2+37\lambda-27=(\lambda-1)\bigl(\lambda^2-10\lambda+27\bigr),
\qquad \lambda=1,\ 5\pm i\sqrt2 ,$$
so the singular points are $t=0,\ 1,\ \tfrac{5\pm i\sqrt2}{27},\ \infty$;
$\mathcal I_1=a^3/g=1331/27$, $\mathcal I_2=a^2/d=121/37$. Integrality is verified
exactly to $n=300$ ($u_{300}$ has $213$ decimal digits); $k=2$; the two dominant
characteristic roots are the complex conjugate pair $5\pm i\sqrt2$, so there is
no archimedean Apéry limit.

**The $\mathcal J$-map test is positive.** `06_jtest.gp` returns
$$h=4,\qquad \deg\mathcal J=16,\qquad \gamma=-1,$$
certified: $\mathcal J=U/V$ with $\deg U=15$, $\deg V=16$ fitted from $38$
coefficients of the $q$-expansion and verified on $73$ further ones. Factoring
$V$, $U$ and $U-1728V$ (`13_jdetail.gp`) gives the complete fibre configuration:

| $t$ | $0$ | $1$ | $\tfrac{5\pm i\sqrt2}{27}$ | $\infty$ |
|---|---|---|---|---|
| pole order of $\mathcal J$ | $4$ | $6$ | $3,\,3$ | — |
| fibre | $I_4$ | $I_6$ | $I_3,\ I_3$ | $IV^*$ |

The full ramification data is
$$\mathcal J^{-1}(\infty):\ 4+6+3+3=16,\qquad
  \mathcal J^{-1}(0):\ 3\cdot5+1=16,\qquad
  \mathcal J^{-1}(1728):\ 2\cdot8=16,$$
$U-1728V$ being a perfect square (so no odd ramification over $j=1728$: no
$III/III^*$ fibres, in agreement with $e_2=0$) and $\mathcal J$ having a
**simple** zero at $t=\infty$ (index $1\equiv1\bmod3$ over $j=0$: the $II$/$IV^*$
signature). Riemann–Hurwitz closes exactly,
$$\sum_P(e_P-1)=12+10+8=30=2\cdot16-2,$$
so $\mathcal J$ is a **Belyi map** — no ramification anywhere else, hence no
further singular fibres — and the projective monodromy is the genus-zero group of
signature $(0;3;\,4\ \text{cusps of widths }3,3,4,6)$ and index $16$: exactly the
$(c,e_2,e_3)=(4,0,1)$ line of the signature table of §4.3, and exactly the
Riemann–Hurwitz bound $6\cdot4+4\cdot1-12=16$ of Theorem F3 with **equality**.

**This is a $K3$, not a rational elliptic surface.** $\sum_{I_n}n_i=16>12$, so
$\chi\ne1$; the Euler numbers give
$$\textstyle\sum e = 4+6+3+3+e(\infty),\qquad e(IV^*)=8\ \Longrightarrow\ \sum e=24,$$
the only value $\le26$ divisible by $12$. Hence $\chi=2$: an **elliptic $K3$
surface over $\mathbf P^1_t$ with exactly five singular fibres
$I_4\,I_6\,I_3\,I_3\,IV^*$**, defined over $\mathbf Q$ with the two $I_3$'s
conjugate over $\mathbf Q(\sqrt{-2})$.

> This is the structural pay-off of going from four singular points to five.
> Proposition K1 of `HERFURTNER_K3_WINDOW.md` showed that with **four** singular
> points $\deg\mathcal J\le12$ always, so *no elliptic $K3$ can ever be the
> Picard–Fuchs system of a three-term Apéry-like row*. With **five** the bound is
> $18$, and the first integral row that meets it is a $K3$. The four-term row is
> the recurrence-level home of the elliptic $K3$s with five singular fibres, just
> as the three-term row is the home of Herfurtner's rational elliptic surfaces
> with four.

Its $\mu$-rescalings $(22,8,148,12,24)$ ($\mu=2$) and $(33,12,333,27,81)$
($\mu=3$) are in the scan output and are the same surface.

### 5.4 The one-parameter families are all disguised: the explicit statement

**Proposition F4 (the generic cusp move, four-term form).** *Let $(P_0,Q_0)$ be a
three-term row in Zagier's class, $P_0=a_0(n^2+n)+c_0$, $Q_0=d_0n^2$, with
integral $u_n$, and let $\nu\in\mathbf Z$ be arbitrary (in particular **not** a
characteristic root). Then the signed binomial transform*
$$v_n=\sum_{m=0}^n\binom nm(-\nu)^{\,n-m}u_m$$
*satisfies the four-term row of class $(\rho;M,j_1,j_2)=(0;1,0,1)$*
$$\boxed{\
 P(n)=(a_0-3\nu)(n^2+n)+(c_0-\nu),\quad
 Q(n)=\bigl(d_0-2a_0\nu+3\nu^2\bigr)n^2,\quad
 R(n)=-\nu\bigl(d_0-a_0\nu+\nu^2\bigr)\,n(n-1),}$$
*i.e. $(a,c,d,f,C)=\bigl(a_0-3\nu,\ c_0-\nu,\ d_0-2a_0\nu+3\nu^2,\ 0,\
-\nu(d_0-a_0\nu+\nu^2)\bigr)$, with $f=Q(0)=0$ identically. Its characteristic
roots are $\{\lambda_1-\nu,\ \lambda_2-\nu,\ -\nu\}$ (the translation
$v\mapsto v-\nu$ of $\{\lambda_1,\lambda_2,0\}$ in the coordinate $v=1/t$), and
$t=\infty$ is an **apparent** singularity: the local system is the original
four-point one.*

*Proof.* $\sum_nv_nt^n=(1+\nu t)^{-1}y\bigl(t/(1+\nu t)\bigr)$, so the singular
points of the transform are the preimages $t=1/(\lambda_i-\nu)$ of
$t_i$, the preimage $t=-1/\nu$ of $s=\infty$, and $t=\infty$ (the image of the
*ordinary* point $s=1/\nu$, carrying only the pole of the gauge factor, hence
exponents $(0,1)$ and no logarithm). $\mathcal R(t)=\prod_i(1-(\lambda_i-\nu)t)\cdot(1+\nu t)$
gives $a,d,g$; $c=v_1=c_0-\nu$; and $f$ falls out of
$4v_2=P(1)v_1-Q(1)$. Verified for Zagier $\mathbf A$–$\mathbf F$ and
$\nu\in\{-3,\dots,3\}$ to $n=18$ in `02_fit.py`. $\square$

For a three-term row of another class the exponents $1-r_i$ at $s=\infty$ are not
integral, and the gauge exponent $\alpha$ must be taken in $\{1-r_1,1-r_2\}$
(the *generalised* signed binomial transform of `CUSP_MOVE_PROGRAM.md` §2.1) —
the difference from Theorem 1 there being only that $\nu$ is now allowed to be a
non-characteristic root, which is exactly what makes the output four-term instead
of three-term. The resulting classes are those with $\delta_\infty\in\mathbf Z$,
and every one of them is populated in the scan.

**This is the four-term analogue of the family classes $(2;-1,1)$ and $(2;1,3)$
of Theorem H5(iv).** There, the infinite families were the $I_m^*/I_0^*$
one-parameter Herfurtner families, i.e. the classical Legendre–Padé rows. Here
they are the cusp-move images: a fifth singular point that is not really there.

### 5.5 The complete list of five-singular-point rows

All $30$ primitive rows with five genuine singular points produced by the
census (boxes G and W together), sorted by score. $(\rho;M,j_1,j_2)$ is
the class, $(a,c,d,f,C)$ the parameters of §2.2, $\lambda_i$ the
characteristic roots (roots of $\lambda^3-a\lambda^2+d\lambda-CM^2$), $k$ the
sharp denominator exponent of the companion, score $=\log(1/|\lambda_2|)-k$
with $\lambda_2$ the second largest in modulus, $\mathcal I_1=a^3/g$,
$\mathcal I_2=a^2/d$, and $\mathcal J$ the outcome of the $\deg\le18$
$\mathcal J$-map test.

| class | $(a,c,d,f,C)$ | $\lambda_1,\lambda_2,\lambda_3$ | $k$ | score | $\mathcal I_1$ | $\mathcal I_2$ | $\mathcal J$ |
|---|---|---|---|---|---|---|---|
| $(0;2,1,1)$ | $(154, 42, -128, -12, 8)$ | $154.8, -0.414+0.1878i, -0.414-0.1878i$ | 2 | $-1.2117$ | $456533/4$ | $-5929/32$ | no |
| $(0;2,1,1)$ | $(87, 24, -152, -16, 16)$ | $88.72, -1, -0.7214$ | 2 | $-2.0000$ | $658503/64$ | $-7569/152$ | no |
| $(1/2;4,5,5)$ | $(72, 6, -192, -60, 8)$ | $74.6, -1.298-0.1732i, -1.298+0.1732i$ | 2 | $-2.2700$ | $2916$ | $-27$ | no |
| $(0;3,-1,4)$ | $(22, 6, -72, 12, 36)$ | $25.35, -1.673+3.16i, -1.673-3.16i$ | 1 | $-2.2741$ | $2662/81$ | $-121/18$ | no |
| $(0;2,1,1)$ | $(22, 6, -96, -12, 24)$ | $25.86, -2, -1.856$ | 2 | $-2.6931$ | $1331/12$ | $-121/24$ | no |
| $(0;2,1,1)$ | $(11, 3, -73, -8, 28)$ | $16, -2.5+0.866i, -2.5-0.866i$ | 2 | $-2.9730$ | $1331/112$ | $-121/73$ | no |
| $(-1/2;4,-1,-1)$ | $(16, 10, 64, 20, 4)$ | $10.47, 4, 1.528$ | 2 | $-3.3863$ | $64$ | $4$ | no |
| $(1/2;4,5,5)$ | $(16, 2, 64, 20, 4)$ | $10.47, 4, 1.528$ | 2 | $-3.3863$ | $64$ | $4$ | no |
| $(0;2,1,1)$ | $(60, 18, 240, 12, 16)$ | $55.71, 4, 0.2872$ | 2 | $-3.3863$ | $3375$ | $15$ | no |
| $(0;2,1,1)$ | $(20, 6, -80, -4, -16)$ | $23.31, -4, 0.6863$ | 2 | $-3.3863$ | $-125$ | $-5$ | no |
| $(0;2,1,1)$ | $(3, 0, -104, -16, 80)$ | $12.94, -5, -4.944$ | 2 | $-3.6094$ | $27/320$ | $-9/104$ | no |
| $(0;3,1,2)$ | $(11, 4, 37, 3, 3)$ | $5-1.414i, 5+1.414i, 1$ | 2 | $-3.6479$ | $1331/27$ | $121/37$ | **yes**, $(h,\deg\mathcal J)=(4,16)$ |
| $(0;2,1,1)$ | $(16, 6, 96, 12, 64)$ | $8, 4-4i, 4+4i$ | 2 | $-3.7329$ | $16$ | $8/3$ | no |
| $(0;2,1,1)$ | $(16, 6, 76, 8, 24)$ | $8, 6, 2$ | 2 | $-3.7918$ | $128/3$ | $64/19$ | no |
| $(0;2,1,1)$ | $(16, 6, 96, 12, 48)$ | $6+3.464i, 6-3.464i, 4$ | 2 | $-3.9356$ | $64/3$ | $8/3$ | no |
| $(0;4,1,3)$ | $(14, 5, 97, 8, 12)$ | $5.5-5.809i, 5.5+5.809i, 3$ | 2 | $-4.0794$ | $343/24$ | $196/97$ | no |
| $(0;2,1,1)$ | $(16, 6, 16, -4, -96)$ | $12, 8, -4$ | 2 | $-4.0794$ | $-32/3$ | $16$ | no |
| $(0;6,-1,7)$ | $(13, 4, 432, -24, -48)$ | $8.262-20.54i, 8.262+20.54i, -3.525$ | 1 | $-4.0975$ | $-2197/1728$ | $169/432$ | no |
| $(0;2,1,1)$ | $(6, 2, 64, 4, -40)$ | $4-8i, 4+8i, -2$ | 2 | $-4.1910$ | $-27/20$ | $9/16$ | no |
| $(0;3,1,2)$ | $(3, 0, -216, -24, 108)$ | $18, -9, -6$ | 2 | $-4.1972$ | $1/36$ | $-1/24$ | no |
| $(0;3,1,2)$ | $(33, 12, 324, 36, 108)$ | $18, 9, 6$ | 2 | $-4.1972$ | $1331/36$ | $121/36$ | no |
| $(0;3,-1,4)$ | $(3, 0, -216, 48, 108)$ | $18, -9, -6$ | 2 | $-4.1972$ | $1/36$ | $-1/24$ | no |
| $(0;2,1,1)$ | $(17, 6, 112, 8, 24)$ | $8-5.657i, 8+5.657i, 1$ | 2 | $-4.2822$ | $4913/96$ | $289/112$ | no |
| $(0;4,1,3)$ | $(28, 10, 208, 12, 16)$ | $16, 10.47, 1.528$ | 2 | $-4.3487$ | $343/4$ | $49/13$ | no |
| $(0;6,1,5)$ | $(17, 6, 56, 0, -12)$ | $10.21+4.658i, 10.21-4.658i, -3.428$ | 2 | $-4.4182$ | $-4913/432$ | $289/56$ | no |
| $(0;2,1,1)$ | $(32, 12, 272, 16, 64)$ | $16, 14.93, 1.072$ | 2 | $-4.7033$ | $128$ | $64/17$ | no |
| $(0;2,1,1)$ | $(42, 15, 441, 24, 100)$ | $25, 16, 1$ | 2 | $-4.7726$ | $9261/50$ | $4$ | no |
| $(-1/3;1,0,0)$ | $(42, 22, 480, 64, -64)$ | $21.07-6.464i, 21.07+6.464i, -0.1318$ | 2 | $-5.0926$ | $-9261/8$ | $147/40$ | no |
| $(1/3;1,1,1)$ | $(42, 8, 480, 64, -64)$ | $21.07-6.464i, 21.07+6.464i, -0.1318$ | 2 | $-5.0926$ | $-9261/8$ | $147/40$ | no |
| $(0;2,1,1)$ | $(64, 24, 1030, 26, 48)$ | $32, 31.81, 0.1886$ | 2 | $-5.4598$ | $4096/3$ | $2048/515$ | no |

Notes.

1. **$k=2$ for $28$ of the $30$**, the four-term echo of Theorem R3. The two
   exceptions have $k=1$ and lie in the classes with $\delta_\infty>1$,
   $(0;3,-1,4)$ and $(0;6,-1,7)$; by the diagnostic of `NONCONGRUENCE_SCAN.md`
   §4.2 a row with $k=1$ is a Padé construction, not a modular object.
2. **Exactly one row is an elliptic surface** — the $K3$ of §5.3. The other
   $29$ fail the $\mathcal J$-map test over the whole $\deg\le18$ window, so their
   projective monodromy is not conjugate into $\mathrm{PSL}_2(\mathbf Z)$: they are
   the four-term half of the non-congruence picture.
3. **Quadratic-twist pairs.** $(16,10,64,20,4)$ in class $(-\tfrac12;4,-1,-1)$ and
   $(16,2,64,20,4)$ in class $(\tfrac12;4,5,5)$ have the same $(a,d,f,C)$ and
   differ only in the accessory parameter $c$ ($10$ versus $2$), hence have the
   same five singular points; likewise $(42,22,480,64,-64)$ in
   $(-\tfrac13;1,0,0)$ and $(42,8,480,64,-64)$ in $(\tfrac13;1,1,1)$. In each
   pair the classes are exchanged by the twist
   $\prod_i(1-\lambda_it)^{\mp1/2}$, which sends $\rho\mapsto-\rho$; they are the
   two half-integral gauges of one projective local system, exactly as
   $\rho\mapsto-\rho$ relates $\mathbf S_1$ and $\mathbf S_{-1}$. Their Apéry
   limits differ.
4. **The best score is $-1.21$**, attained by $(154,42,-128,-12,8)$ in the
   semistable class $\mathbf S_0$, whose two small characteristic roots are the
   complex pair $-0.4140\pm0.1878i$ of modulus $0.4546$. Its Apéry limit is
   $\xi=0.0235858030130198228326871528054644790660665933510869454291671\dots$,
   unidentified. The next best, $(87,24,-152,-16,16)$, has
   $\lambda^3-87\lambda^2-152\lambda-64=(\lambda+1)(\lambda^2-88\lambda-64)$, the
   quadratic being that of the level-$5$ Fricke row
   $(88n^2+44n+6,\ -4(4n-1)(4n-3))$ of `HERFURTNER_CLASSIFICATION.md` §6.1 with
   one extra genuine $I_n$ fibre at $\lambda=-1$; similarly
   $(20,6,-80,-4,-16)$ has $(\lambda+4)(\lambda^2-24\lambda+16)$, the quadratic
   being $\sqrt T$'s. Adding a fifth *genuine* singular point to a known
   four-point configuration is a recurring motif in the list.
5. The Apéry limits of the $21$ rows that have one are listed in `out/xi.txt`;
   none was identified (§6).
---

## 6. Periods

For each surviving row the companion $b_n$ ($b_0=0$, $b_1=1$, same recurrence
from $n=1$) satisfies $LB=t$, and $\xi=\lim b_n/a_n$ exists exactly when the
largest characteristic root is real and simple, the error being
$O\bigl((\lambda_2/\lambda_1)^n\bigr)$. Nine of the thirty rows have a complex
conjugate dominant pair and therefore **no archimedean Apéry limit** — the
situation of Zagier $\mathbf B$ and $\mathrm{AZ}(11,5,125)$, and in particular
the situation of the $K3$ row of §5.3, whose dominant roots are $5\pm i\sqrt2$.

The limits that do exist were computed to $\ge60$ digits (`05_report.py`,
adaptive length up to $n=2\times10^4$ at $200$ digits of working precision) and
put through two batteries.

* `07_ident.gp`: `lindep` of $(1,\xi,X)$ and $(1,\xi,X,Y)$ over $103$ classical
  weight-two constants — $\zeta(2)$; $L(2,\chi_D)$ for the fundamental
  discriminants with $|D|\le24$; $\log^2m$, $\pi\log m$,
  $\pi\log\frac{\sqrt m+1}{\sqrt m-1}$, $\log^2\frac{\sqrt m+1}{\sqrt m-1}$;
  $\pi\arctan\frac1m$ and $\arctan^2\frac1m$; $\Gamma(1/4)^4/\pi$,
  $\Gamma(1/3)^6/\pi^2$; $\zeta(2)/\sqrt m$, $\pi\log2/\sqrt m$,
  $\pi\log3/\sqrt m$; $\log m_1\log m_2$; and $\mathrm{Li}_2$ at twenty small
  rational arguments.
* `12_polylog.gp`: for each row separately, the weight-two polylogarithm space
  attached to *its own* singular points — $\zeta(2)$, $\log\alpha\log\beta$ and
  $\mathrm{Li}_2(z)$ for $\alpha,\beta,z$ in the multiplicative group generated by
  the ratios $\lambda_i/\lambda_j$ and $1-\lambda_i/\lambda_j$ — reduced to a
  *greedy independent subset* first (without that reduction the classical
  identities $\mathrm{Li}_2(\tfrac12)=\tfrac{\zeta(2)}2-\tfrac{\log^22}2$,
  $\mathrm{Li}_2(-1)=-\tfrac{\zeta(2)}2$, Landen, … make the basis dependent and
  `lindep` returns them instead of anything about $\xi$).

**All $21$ limits came back unidentified**, in both batteries, at $60$ digits
with coefficient heights up to $10^6$ (single constant) and $10^5$ (pairs) in
`07_ident.gp`, and in spaces of dimension $7$–$11$ in `12_polylog.gp`. This is a
negative result but a fairly sharp one, and it is consistent with the geometry:
a five-point configuration is not a modular curve (§1), so the period of the
associated local system has no reason to be a special value of an $L$-function.

### 6.1 Corollary F6: the score

The score of a row is $\log(1/|\lambda_2|)-k$ with $\lambda_2$ the second largest
characteristic root in modulus; a positive score is the archimedean part of an
irrationality proof. Two complete sweeps settle it.

* **Box W** ($|\lambda_2|<1$; $a\le300$ for the $18$ classes of `classes2.txt`,
  $a\le400$ for $\mathbf S_0$): $79$ rows have $|\lambda_2|<1$, all in the classes
  $(\pm1;3,\ldots)$ and $(\pm1;4,\ldots)$, and **every one of them has a vanishing
  Frobenius obstruction** — a disguised three-term row. In $\mathbf S_0$ the only
  rows with $|\lambda_2|\le1$ are the two genuine five-point rows
  $(87,24,-152,-16,16)$ ($\lambda_2=-1$ exactly) and $(154,42,-128,-12,8)$
  ($|\lambda_2|=0.4546$), of scores $-2$ and $-1.2117$.
* **Box PW** ($|\lambda_2|<e^{-1}$, hence complete for a positive score whenever
  $k\ge1$; $a\le800$, $|c|\le60$, $|f|\le120$, all $27$ classes swept to
  completion): the box contains **exactly two integral rows in total**, both in
  the class $(0;2,-1,3)$ ($\delta_\infty=2$), namely
  $$(a,c,d,f,C)=(144,32,96,-32,4):\ \lambda=143.331,\,0.3506,\,0.3184,\quad
    k=1,\ \text{score}=+0.0480,$$
  $$(a,c,d,f,C)=(225,50,-120,40,4):\ \lambda=225.532,\,-0.2662\pm0.0092i,\quad
    k=1,\ \text{score}=+0.3229,$$
  and **both have an exactly vanishing Frobenius obstruction at $t=\infty$**
  (the value is $0$, not merely small). The exponents there are
  $(\tfrac12,\tfrac52)$, so the local monodromy is $-I$: an $I_0^*$ fibre,
  invisible to the projective local system, which therefore has only four
  singular points. These two rows are the four-term avatars of the classical
  Legendre–Padé family — the same $I_0^*$ mechanism that made the class
  $(2;1,3)$ an infinite family with $k=1$ and unbounded score in Theorem H5(iv),
  and $k=1$ identifies them as Padé constructions in the sense of
  `NONCONGRUENCE_SCAN.md` §4.2.

> **Corollary F6.** Over the whole scanned region — every Kodaira-admissible
> class, both windows — **no genuine five-singular-point four-term row has a
> positive score**. The only positive scores that occur belong to disguised
> three-term rows with $k=1$, exactly as in the three-term case, where every
> positive score outside Apéry's and Beukers' rows has $k=1$ and is a classical
> Legendre–Padé row (Corollary H6). The four-term world produces **no new
> irrationality candidate**.

---

## 7. Status: what is proved, what is verified, what is open

**Proved.**

* Lemma 1.1 and the exponent dictionary (§1): the operator identity, the $D$-form,
  $T=\mathcal S-t\mathcal R'$, $\rho_i=-T(t_i)/(t_i\mathcal R'(t_i))$, the
  exponents $2-s_i$ at $\infty$, and the Fuchs relation
  $\sum\rho_i=s_1+s_2-1$. Elementary; all four checked symbolically in
  `01_exponents.py`.
* Theorem F2 (equal exponents $\Leftrightarrow$ $b=(1-\rho)a$, $e=-2\rho d$,
  $h=-(1+3\rho)g$) and Corollary F2.1 (Galois). Elementary.
* Corollary F2.2 (Kodaira admissibility of the exponent differences), given
  Kodaira's list; the proof of Theorem H2 verbatim.
* Theorem F3 ($\deg\mathcal J\le6c+4e_3+3e_2-12\le18$) and the Gauss–Bonnet
  signature list of §4.3. Both elementary; the two agree, and equality holds iff
  $\mathcal J$ is a Belyi map.
* The identification of the accessory parameters ($c=u_1$, $f=Q(0)$; exactly
  $5-3=2$) and of the Möbius invariants $\mathcal I_1=a^3/g$, $\mathcal I_2=a^2/d$.
* Proposition F4 (the generic cusp move sends a Zagier-class three-term row to a
  four-term row in class $(0;1,0,1)$ with $f=0$ and an apparent point at
  $\infty$) — proved, and checked on Zagier $\mathbf A$–$\mathbf F$ for
  $|\nu|\le3$.
* Proposition F5 (the obstruction $S_i=t_i[(a-c)+(f-d)t_i]$ in class
  $(-1;0)$, and the constant-coefficient degeneration).
* (D7): there are no $t\mapsto t^2$ pullbacks, because $\mathcal R$ has odd
  degree and $g\ne0$.
* The Casoratian identity $W_{n+1}=R(n)W_n/(n+1)^2$, hence the Casoratian
  degeneracy of every class in which $R$ has a positive integer root.

**Verified (exact finite computation, not proof).**

* The census of §5 over the stated boxes; integrality prime-tested to $n=24$ and
  re-verified with exact rational arithmetic to $n=120$ for every listed row.
* The Frobenius obstruction (apparent-singularity) verdicts, computed at $100$
  digits with a $10^{-40}$ relative threshold. In the cases that matter most —
  the two positive-score rows of §6.1 — the obstruction evaluates to exactly $0$,
  not merely to something small; in general the verdict is numerically decisive
  but not a certificate.
* The $\mathcal J$-map verdicts. A *positive* answer is a certificate (an explicit
  rational $\mathcal J$ verified on more $q$-coefficients than were used to fit
  it: $73$ extra for the $K3$ row). A *negative* answer means only "no fit with
  $h\le18$, $\deg\le18$ and $\gamma$ in the $760$-element list".
* The fibre configuration $I_4\,I_6\,I_3\,I_3\,IV^*$ of the $K3$ row, read off
  the certified $\mathcal J$ by factoring $V$, $U$ and $U-1728V$.
* $k=2$ for $28$ of the $30$ surviving rows and $k=1$ for the other two,
  measured to $n=60$.
* Corollary F6: the completeness of box PW rests on the elementary bounds
  $|d|\le0.74(|a|+1)+2$, $|g|\le0.136(|a|+1)+1$ (proved), and on the sweep of that
  box for $a\le800$ over all $27$ classes (computed).
* The unidentifiability of the $21$ Apéry limits, over the stated bases and
  heights.

**Open / not settled here.**

1. *The Euler number of the $K3$.* $\sum_{I_n}n_i=16$ forces $\chi\ge2$; $\chi=2$
   with an $IV^*$ at $\infty$ is the only possibility with $\sum e\le26$, but a
   fibre of type $II$ at $\infty$ together with an $I_0^*$ somewhere (invisible
   to the projective local system) would also give $\sum e=24$. Distinguishing
   them needs the Weierstrass model, not done here.
2. *A five-fibre analogue of Herfurtner's table.* The configurations of elliptic
   surfaces with five singular fibres form one-parameter families; the locus they
   cut out in the $(\mathcal I_1,\mathcal I_2)$-plane is a curve for each fibre
   type. Computing those curves would give a genuine analogue of Theorem H3 and
   would predict the $K3$ row instead of discovering it.
3. *The enumeration of the fifteen five-point signatures* into conjugacy classes
   of genus-zero subgroups of $\mathrm{PSL}_2(\mathbf Z)$ (the analogue of the
   $28$ groups of Theorem N2) is not carried out.
4. *The periods.* All surviving Apéry limits are unidentified. They are periods
   of rank-two local systems on $\mathbf P^1$ minus five points; the natural
   guess is that they are $\mathbf Q$-linear combinations of weight-two
   *multiple* polylogarithms at arguments in the field generated by the
   $\lambda_i$ — the batteries here only tested the classical and the depth-one
   parts of that space.
5. *Unequal $\rho_i$.* Corollary F2.1 confines these to reducible cubics, and the
   examples produced by the gauge $y\mapsto y/(1-\nu t)$ all have an apparent
   point; but the scan only covers equal-$\rho$ classes. A class with
   $\rho_1\ne\rho_2=\rho_3$ has $b,e,h$ determined as *rational* functions of
   $a,d,g$ and would need its own parametrisation.
6. *$p$-adic invariants.* Nothing here computes $p$-adic Apéry limits or the
   Dwork crystal side for the $K3$ row.

## 8. Reproduction

```
cd lattice/four_term
python3 01_exponents.py            # the dictionary, verified symbolically
python3 02_fit.py                  # the gauge and binomial-transform constructions
python3 04_classes.py              # the Kodaira-admissible classes
gcc -O3 -march=native -o 03_fscan 03_fscan.c

# general census (box G)
./run_scan.sh    G   200 40 1200 80 120 24 12 classes.txt
./run_scan.sh    G2  200 40 1200 80 120 24 12 classes2.txt
./run_scan.sh    G3  140 40 1200 80 120 24 12 classes3.txt
# |lambda_2| < 1 window (box W), complete
./run_window.sh  W   300 60      120     24 12 classes_all.txt
# positive-score window (box PW), complete for score > 0 with k >= 1
./run_pwindow.sh PW  800 60      120     24 12 classes_ord.txt

python3 08_analyse.py out/analysis_all.json out/G*_*.txt out/P_*.txt out/W*_*.txt
python3 09_full.py    out/analysis_all.json out/full.json
python3 10_jrun.py    out/analysis_all.json 11_jall.gp
gp -q -s 6000000000 11_jall.gp < /dev/null > out/jall.log
gp -q -s 4000000000 07_ident.gp                     # reads out/xi.txt
```

The $K3$ fibre configuration of §5.3:

```
gp -q -s 6000000000 <<'EOF'
\r 06_jtest.gp
\r 13_jdetail.gp
jdetail(0,1,3,1,2, 11,4,37,3,3);
EOF
```

Data files: `out/analysis_all.json` (every hit with its verdict), `out/full.json`
(the candidate rows with $k$, $\xi$, score, $\mathcal I_1$, $\mathcal I_2$),
`out/jall.log` (the $\mathcal J$-map verdicts), `out/ident.log` (the period battery),
and `out/xi.txt` (the Apéry limits). The filter-stage counts of §5.1 are read off `out/analysis_all.json`.

**Sources.** S. Herfurtner, *Elliptic surfaces with four singular fibres*, Math.
Ann. **291** (1991) 319–342; U. Persson, *Configurations of Kodaira fibers on
rational elliptic surfaces*, Math. Z. **205** (1990) 1–47; R. Miranda, *The basic
theory of elliptic surfaces*, Pisa 1989 (Kodaira's dictionary for the functional
invariant, §IV); A. Beauville, C. R. Acad. Sci. Paris **294** (1982) 657–660;
D. Zagier, *Integral solutions of Apéry-like recurrence equations* (2009);
and, inside this project, `HERFURTNER_CLASSIFICATION.md`,
`HERFURTNER_K3_WINDOW.md`, `CUSP_MOVE_PROGRAM.md`, `NONCONGRUENCE_SCAN.md`,
`ROOT_ROWS.md`.
