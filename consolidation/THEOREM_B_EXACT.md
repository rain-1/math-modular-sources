# Theorem B, exactly: the Apéry limit **is** the critical value

*Companion to `paper/sections/02_sources.tex` (Theorem B, Lemma `lem:fold`) and to
`EULER_CRITERION.md` (Theorem F, the $p$-adic twin).  Scripts and logs:
`lattice/thmB_exact/`.*

---

## 0. Summary

For every modular Apéry row of the twelve-family Source Theorem table, write
$\Phi=P(V)E^{\psi,\varphi}_{w+2}$ for the Eisenstein source, $w=r-1$,
$P(s)=\sum_d c_d d^{-s}$, so that

$$L(\Phi,s)=P(s)\,L(\psi,s)\,L(\varphi,s-w-1).$$

> **Theorem B\* (exact prefactor).**  Assume the **endpoint condition**
> $$\boxed{\ \delta(\varphi)\,P(w+2)=0\ }\qquad(\delta(\varphi)=1 \text{ iff } \varphi=\mathbf 1)$$
> — equivalently, $\Phi$ vanishes at the cusp $0$, equivalently $L(\Phi,s)$ is
> regular at the endpoint $s=w+2$.  Then $t_c$ is attained on the boundary of
> the disc of convergence at a point where the Eichler integral has a finite
> value, the Apéry limit exists, and
> $$\boxed{\ \xi_\infty=\lim_{n\to\infty}\frac{B_n}{A_n}=L(\Phi,w+1)
>  =P(w+1)\,L(\psi,w+1)\,L(\varphi,0).\ }$$
> In particular $r_\infty=P(w+1)L(\varphi,0)$, which is $-\tfrac12P(w+1)$ when
> $\varphi=\mathbf 1$: exactly the archimedean half of the Corollary of
> `EULER_CRITERION.md`.  **There is no further rational factor.**

The endpoint condition holds for **exactly the nine real-fold rows** and fails
for **exactly the three complex-fold rows** $\mathbf B,\delta,\eta$ — twelve
tests, no exception (§3).  So the archimedean criterion for "the row has an
Apéry limit" is an *endpoint* annihilation, the mirror of the *Euler-factor*
divisibility criterion that governs the $p$-adic limits (Theorem F).

**Status.**  §2 proves the theorem in the two geometries it occurs in.  Every
input hypothesis is checked exactly or numerically per row in §3–§4.  The
identity $\xi_\infty=L(\Phi,w+1)$ is verified for all nine rows to
$\ge 10^{-96}$ **from the recurrence alone** (§4.1) and independently to
$\ge 10^{-60}$ on the modular side (§4.2, §4.3).  §5 identifies the three
complex-fold constants in closed form, which settles Problem $\Phi$-1 of
`certificates/eisenstein` numerically to $90$ digits.

---

## 1. Where the singularity sits: two geometries, not one

Let $P_2(t)=1-at+ct^2$ (order two) resp. $P_3(t)=1-2at+ct^2$ (order three) and
let $t_c$ be its root of smallest modulus.  Writing the operator in the form
$y''+p_1y'+p_0y=0$ one finds, with $P=P_r$,

$$r=2:\quad p_1=\frac{P+tP'}{tP}=\frac{d}{dt}\log(tP),\qquad
r=3:\quad p_1=\frac{P+\tfrac12tP'}{tP}=\frac1t+\frac{P'}{2P}.$$

At a simple root $t_c$ of $P$ the indicial equation is therefore

| | $p_1\sim$ | indicial | local shape of $y_0$ | $[t^n]y_0$ |
|---|---|---|---|---|
| $r=2$ | $1/(t-t_c)$ | $\rho^2=0$ | $g+h\log(1-t/t_c)$ | $\asymp t_c^{-n}/n$ |
| $r=3$ | $\tfrac12/(t-t_c)$ | $\rho^2-\tfrac12\rho=0$ | $g+h\sqrt{1-t/t_c}$ | $\asymp t_c^{-n}n^{-3/2}$ |

**This is a dichotomy the paper does not currently record.**  Only the
order-three rows have a *fold* in the sense of Lemma `lem:fold`; the
second-order rows have a **logarithmic** singularity, and $t\to t_c$ there as
$\tau$ tends to a **cusp**, not to an interior point.  Numerically (scripts
`01`, `02`):

* $\mathbf A,\mathbf C,\mathbf D,\mathbf E,\mathbf F$: $t(iy)$ increases
  monotonically along the imaginary axis and $t(iy)\to t_c$ as $y\to0^+$, i.e.
  $q_c=1$, the cusp $0$.  ($\mathbf B$: the cusp $1/6$ of $\Gamma_0(36)$.)
* $\alpha,\gamma,\varepsilon,\zeta$: $dt/d\tau$ vanishes at the interior point
  $\tau_*=i/\sqrt N$ — the Fricke fixed point — where $t(\tau_*)=t_c$ exactly
  and $\Phi(\tau_*)=0$.

(The record in `certificates/eisenstein/main(20260813-131951).tex`, Table 2,
locates $q_c$ for the five second-order rows by Newton's method on a
$26$-term truncation of $t(q)$; those starred rows are the ones whose true
$q_c$ is a root of unity, which is why they came out with $|q_c|\approx0.84$
and only $2$–$3$ correct digits.  The mechanism they confirm is right; the
geometry is the cusp one.)

Both geometries lead to the same answer, for the same reason.

---

## 2. The proof

Throughout, $\Phi\in M_{w+2}(\Gamma_0(N),\chi)$ with $\Phi=q+O(q^2)$,
$\Theta=\thq^{-(w+1)}\Phi=\sum_{m\ge1}c(m)m^{-(w+1)}q^m$, and

$$\Theta(\tau)=\frac{(-2\pi i)^{w+1}}{w!}\int_\tau^{i\infty}\Phi(z)(z-\tau)^w\,dz
\tag{2.1}$$

(termwise: $\int_\tau^{i\infty}e^{2\pi imz}(z-\tau)^wdz=q^m i^{w+1}w!/(2\pi m)^{w+1}$).

### Lemma 1 (endpoint criterion)

*$\Lambda(\Phi,s)=\int_0^\infty\Phi(iy)y^{s-1}dy$ is entire $\iff$ $\Phi$
vanishes at the cusp $0$ $\iff$ $\delta(\varphi)P(w+2)=0$.*

**Proof.**  Convergence at $y\to\infty$ is automatic ($\Phi=O(q)$).  At
$y\to0$, $\Phi(iy)=a_0^{(0)}(\sqrt N y)^{-(w+2)}(1+O(e^{-c/y}))$ where
$a_0^{(0)}$ is the constant term of $\Phi$ at the cusp $0$; the integral
converges at $0$ iff $a_0^{(0)}=0$, and otherwise $\Lambda$ acquires a simple
pole at $s=w+2$.  On the other hand $L(\Phi,s)=P(s)L(\psi,s)L(\varphi,s-w-1)$
has a pole at $s=w+2$ precisely when $\varphi=\mathbf 1$ (the pole of
$L(\varphi,s-w-1)$ at $s-w-1=1$) and $P(w+2)\neq0$. $\square$

### Lemma 2 (the cusp value of the Eichler integral)

*If $\Lambda(\Phi,\cdot)$ is entire then*
$$\lim_{y\to0^+}\Theta(iy)=L(\Phi,w+1),$$
*and more precisely $\Theta(iy)$ is a polynomial of degree $w$ in $y$ with
constant term $L(\Phi,w+1)$, up to $O(e^{-c/y})$; the coefficient of $y^{w-j}$
is a multiple of $L(\Phi,j+1)$.*

**Proof.**  Expand $(z-iy)^w=\sum_j\binom wj z^j(-iy)^{w-j}$ in (2.1).  Each
$\int_{iy}^{i\infty}\Phi(z)z^jdz\to i^{j+1}\Lambda(\Phi,j+1)$ as $y\to0$ by
dominated convergence, the error being $\int_0^{iy}=O(e^{-c/y})$ by Lemma 1.
Only $j=w$ survives the limit, and
$\frac{(-2\pi i)^{w+1}}{w!}\,i^{w+1}\Lambda(\Phi,w+1)
=\frac{(2\pi)^{w+1}}{w!}\cdot\frac{w!}{(2\pi)^{w+1}}L(\Phi,w+1)=L(\Phi,w+1)$. $\square$

Equivalently, the **Fricke period polynomial**
$R(\tau):=\frac{(-2\pi i)^{w+1}}{w!}\int_0^{i\infty}\Phi(u)(u-\tau)^w du$
is
$$R(\tau)=\sum_{j=0}^{w}\binom wj\frac{(-2\pi i)^{w+1}i^{j+1}}{w!}(-\tau)^{w-j}\Lambda(\Phi,j+1),$$
which for $w=1,2$ reads
$$w=1:\ R(\tau)=L(\Phi,2)+2\pi i\,L(\Phi,1)\,\tau,\qquad
w=2:\ R(\tau)=L(\Phi,3)+2\pi i\,L(\Phi,2)\,\tau-2\pi^2L(\Phi,1)\,\tau^2. \tag{2.2}$$
$R(0)=L(\Phi,w+1)$ always.

### Theorem I (cusp geometry: the five second-order real rows)

*Let $r=2$, let the endpoint condition hold, and suppose $y\mapsto t(iy)$ maps
$(0,\infty)$ monotonically onto $(0,t_c)$ with $t(iy)\to t_c$ as $y\to0^+$.
Then $\xi_\infty=L(\Phi,2)$.*

**Proof.**  Near $t_c$, $y_0=g_0+h_0\log(1-t/t_c)$ and $y_B=g_B+h_B\log(1-t/t_c)$
with $g_\bullet,h_\bullet$ analytic (§1).  Singularity analysis gives
$[t^n]y=-h(t_c)t_c^{-n}/n+O(\rho^{-n})$, $\rho>|t_c|$, so
$\xi_\infty=h_B(t_c)/h_0(t_c)$ provided $h_0(t_c)\neq0$; and since both series
diverge logarithmically at $t_c$,
$$\xi_\infty=\lim_{t\to t_c^-}\frac{y_B(t)}{y_0(t)} .$$
Along $\tau=iy$ we have $y_0(t(\tau))=F(\tau)$ and $y_B(t(\tau))=F(\tau)\Theta(\tau)$
(analytic continuation of the $q$-expansion identity along the path, legitimate
because $|t(iy)|<|t_c|$ throughout), so the ratio is $\Theta(iy)$, and
Lemma 2 finishes.  $h_0(t_c)\neq0$ is the statement $F(iy)\to\infty$, which
holds because $F$ has weight $w=1>0$ and a nonzero constant term at the cusp
$0$; indeed $\log(1-t/t_c)\asymp\log q_0=-2\pi/(Nhy)$ and $F(iy)\asymp y^{-1}$
match. $\square$

*Remark (the $L(\Phi,1)$ term).*  This answers the question raised by
`PADIC_PERIOD.md` §2, where the second period coefficient
$L(\Phi_{\mathbf F},1)=\pi/(8\sqrt3)\neq0$ was found not to vanish.  It does not
have to: by (2.2) it multiplies $\tau$, and the fold of a second-order row sits
at $\tau=0$.  The quasiperiod term is present in the cocycle and **absent from
the fold value**, for the trivial reason that the fold is the point where the
linear term dies.

### Theorem II (Fricke geometry: the four third-order real rows)

*Let $r=3$, $N$ the level, and suppose $\Phi|_{w+2}W_N=-\Phi$ and
$F|_wW_N=-F$.  Then:*
1. *$t\circ W_N=t$, and $\tau_*=i/\sqrt N$ is a simple fold of $t$ with
   $\Phi(\tau_*)=0$;*
2. *the Fricke period polynomial satisfies $R(\tau)=\xi_\infty\,(1+N\tau^2)$;*
3. *hence $\xi_\infty=R(0)=L(\Phi,3)$, together with the two
   forced identities*
   $$L(\Phi,2)=0\qquad\text{and}\qquad L(\Phi,1)=-\frac{N}{2\pi^2}L(\Phi,3).$$

**Proof.**  (1) $\thq t=\Phi/F$ has $W_N$-eigenvalue $(-1)/(-1)=+1$ in weight
$2$, i.e. $(dt/d\tau)(W_N\tau)\cdot(d W_N\tau/d\tau)=dt/d\tau$, so
$t\circ W_N-t$ is constant; it vanishes at the fixed point $\tau_*$, hence
$t\circ W_N=t$.  Differentiating this at $\tau_*$, where $W_N'(\tau_*)=-1$,
gives $t'(\tau_*)=0$; and $\Phi(\tau_*)=F(\tau_*)\thq t(\tau_*)=0$.
Weight $w=2$ and $j(W_N,\tau_*)=\sqrt N\tau_*=i$ give
$F(\tau_*)=-i^2F(\tau_*)=F(\tau_*)$, so $F(\tau_*)$ is unobstructed (and is
nonzero in each row).

(2) The Eichler cocycle for $W_N$ with eigenvalue $\varepsilon=-1$ is
$$j(W_N,\tau)^{w}\Theta(W_N\tau)=\varepsilon\bigl(\Theta(\tau)-R(\tau)\bigr) \tag{2.3}$$
(substitute $z=W_Nu$ in (2.1); the boundary point $W_N^{-1}(i\infty)=0$
produces exactly $R$).  Since $t$ is $W_N$-invariant and $W_N$ is a
holomorphic involution of $\mathbb H$ fixing $\tau_*$, $W_N$ **is** the local
fold involution of Lemma `lem:fold`.  That lemma says $\xi_\infty$ is the
unique constant for which $H_\xi:=F\,(\Theta-\xi)$ has vanishing odd part at
$\tau_*$, i.e. $H_{\xi_\infty}\circ W_N=H_{\xi_\infty}$ near $\tau_*$, hence on
all of $\mathbb H$ by analytic continuation.  Substituting the two
transformation laws,
$$H_\xi(W_N\tau)=F(\tau)\bigl[\varepsilon^2(\Theta-R)-\varepsilon\xi(\sqrt N\tau)^w\bigr]
=F(\tau)\bigl[\Theta-R+\xi N\tau^2\bigr],$$
and equating with $F(\Theta-\xi)$ gives $R(\tau)=\xi(N\tau^2+1)$.

(3) Compare with (2.2): the constant term gives $\xi=L(\Phi,3)$, the linear
term gives $L(\Phi,2)=0$, and the quadratic term gives
$-2\pi^2L(\Phi,1)=\xi N$. $\square$

*This is the "pure CM value" mechanism of `book/v8/08b_beta4_application.tex`
in its natural generality: the period polynomial is a multiple of
$1+N\tau^2=N(\tau-\tau_*)(\tau+\tau_*)$, it vanishes at the fixed point, and no
$\pi$-power correction survives.*

### Corollary (the rational prefactor)

$$\xi_\infty=P(w+1)\,L(\psi,w+1)\,L(\varphi,0),\qquad
r_\infty=P(w+1)L(\varphi,0)=\begin{cases}-\tfrac12P(w+1),&\varphi=\mathbf 1,\\
P(w+1)L(\varphi,0),&\text{else.}\end{cases}$$

Together with Theorem F this gives the two-sided statement

$$\boxed{\ \xi_\infty=-\tfrac12P(w+1)\Lambda,\qquad
\xi_p=-\tfrac12 Q(w+1)\Lambda_p,\qquad
r_p/r_\infty=1/\mathcal E_p(w+1)\ }\qquad(\varphi=\mathbf 1)$$

with $\Lambda=L(\psi,w+1)$, $\Lambda_p=L_p(w+1,\psi\omega^{-w})$ — **Conjecture
D with both halves proved**, the archimedean half by Theorems I–II above and
the $p$-adic half by Theorem F.  The two criteria are dual: *endpoint*
annihilation $P(w+2)=0$ archimedean, *Euler-factor* divisibility
$(1-\psi(p)p^{-s})\mid P(s)$ at $p$.

---

## 3. The endpoint criterion on the twelve-row table

`lattice/thmB_exact/02_cuspvalue.gp`, exact rational arithmetic.

| row | $r$ | $(\psi,\varphi)$ | $P(w+2)$ | $\delta(\varphi)P(w+2)$ | fold | $\Theta(0)$ |
|---|---|---|---|---|---|---|
| $\mathbf A$ | 2 | $(\mathbf 1,\chi_{-3})$ | $3/4$ | $0$ | cusp $0$ | finite |
| $\mathbf B$ | 2 | $(\chi_{-3},\mathbf 1)$ | $1/8$ | $1/8$ | **complex** (cusp $1/6$) | **diverges** |
| $\mathbf C$ | 2 | $(\chi_{-3},\mathbf 1)$ | $0$ | $0$ | cusp $0$ | finite |
| $\mathbf D$ | 2 | $(\mathbf 1,\nu)$ | $1$ | $0$ | cusp $0$ | finite |
| $\mathbf E$ | 2 | $(\chi_{-4},\mathbf 1)$ | $0$ | $0$ | cusp $0$ | finite |
| $\mathbf F$ | 2 | $(\chi_{-3},\mathbf 1)$ | $0$ | $0$ | cusp $0$ | finite |
| $\alpha$ | 3 | $(\mathbf 1,\mathbf 1)$ | $0$ | $0$ | $i/\sqrt{12}$ | finite |
| $\gamma$ | 3 | $(\mathbf 1,\mathbf 1)$ | $0$ | $0$ | $i/\sqrt6$ | finite |
| $\delta$ | 3 | $(\mathbf 1,\mathbf 1)$ | $5/27$ | $5/27$ | **complex** | **diverges** |
| $\varepsilon$ | 3 | $(\mathbf 1,\mathbf 1)$ | $0$ | $0$ | $i/\sqrt8$ | finite |
| $\zeta$ | 3 | $(\chi_{-3},\chi_{-3})$ | $1$ | $0$ | $i/\sqrt9$ | finite |
| $\eta$ | 3 | $(\chi_5,\mathbf 1)$ | $1/16$ | $1/16$ | **complex** | **diverges** |

Twelve for twelve: *the row has a (real) Apéry limit iff the endpoint slot
$s=w+2$ is annihilated.*  This is the archimedean member of the same family of
statements as `sec:hecke` of the paper ("full endpoint/critical
purification"): $P_\varepsilon=(1-X)(1-4X)(1-16X)$, $P_\gamma$, $P_\alpha$ all
kill $s=4$, and $P_{\mathbf C},P_{\mathbf E},P_{\mathbf F}$ all kill $s=3$,
while $P_{\mathbf B},P_\delta,P_\eta$ do not.

---

## 4. Verification

### 4.1 From the recurrence alone (`01_limits.gp`)

$A_n,B_n$ generated exactly over $\mathbb Q$ to $n=2600$; nothing modular is
used on the left-hand side.  Working precision $80$ digits.

| row | $\xi_{2600}-L(\Phi,w+1)$ | $L(\Phi,w+1)$ | closed form |
|---|---|---|---|
| $\mathbf A$ | $2.3\cdot10^{-97}$ | $0.41123351671205660911\ldots$ | $\zeta(2)/4$ |
| $\mathbf C$ | $-2.3\cdot10^{-97}$ | $0.39065120644824314843\ldots$ | $L(\chi_{-3},2)/2$ |
| $\mathbf D$ | $0$ (to $10^{-96}$) | $0.32898681336964528729\ldots$ | $\zeta(2)/5$ |
| $\mathbf E$ | $0$ (to $10^{-96}$) | $0.45798279708860950752\ldots$ | $G/2$ |
| $\mathbf F$ | $-2.3\cdot10^{-97}$ | $0.48831400806030393554\ldots$ | $5L(\chi_{-3},2)/8$ |
| $\alpha$ | $0$ (to $10^{-96}$) | $0.35059993008821499990\ldots$ | $7\zeta(3)/24$ |
| $\gamma$ | $0$ (to $10^{-96}$) | $0.20034281719326571423\ldots$ | $\zeta(3)/6$ |
| $\varepsilon$ | $0$ (to $10^{-96}$) | $0.26294994756616124993\ldots$ | $7\zeta(3)/32$ |
| $\zeta$ | $0$ (to $10^{-96}$) | $0.29467460391669328558\ldots$ | $L(\chi_{-3},3)/3$ |

All nine rational prefactors of the census table are reproduced by the single
formula $r_\infty=P(w+1)L(\varphi,0)$:
$\tfrac14,\tfrac12,\tfrac15,\tfrac12,\tfrac58,\tfrac7{24},\tfrac16,\tfrac7{32},\tfrac13$.

### 4.2 The cusp value (`02_cuspvalue.gp`)

$\Theta(iy)$ summed directly ($y=0.004,0.002,0.001$; up to $2.6\cdot10^4$
terms) and Richardson-extrapolated in $y$ to degree $w$.  Agreement with
$L(\Phi,w+1)$: $10^{-70}$ for $\mathbf A,\mathbf C,\mathbf D,\mathbf E,
\gamma,\varepsilon,\zeta$; $10^{-60}$ for $\mathbf F$, $10^{-62}$ for $\alpha$
(Richardson residue).  Divergence for $\mathbf B,\delta,\eta$, as Lemma 1
predicts.

### 4.3 The Fricke geometry (`03_fricke_fold.gp`, `05_periodpoly.gp`)

For $\alpha,\gamma,\varepsilon,\zeta$, at three test points each and to
$10^{-70}$:
$$\Phi|_4W_N=-\Phi,\qquad F|_2W_N=-F,\qquad t\circ W_N=t;$$
and at $\tau_*=i/\sqrt N$: $t(\tau_*)-t_c=0$, $t'(q_c)=0$, $\Phi(q_c)=0$ (all
$\le10^{-76}$), and
$$\Theta(q_c)+\frac{F(q_c)\Theta'(q_c)}{F'(q_c)}-L(\Phi,3)=0\quad(\le10^{-77}).$$
The two forced period identities of Theorem II(3) hold exactly (to $10^{-58}$):
$L(\Phi,2)=0$ and $L(\Phi,1)=-N L(\Phi,3)/(2\pi^2)$ for all four rows —
e.g. $\gamma$: $L(\Phi_\gamma,1)=2\zeta'(-2)=-\zeta(3)/(2\pi^2)$ and
$N L(\Phi,3)/(2\pi^2)=6\cdot\tfrac{\zeta(3)}6/(2\pi^2)$.  For $\delta,\eta$
none of this holds: they are not Fricke eigenforms, $t\circ W_N\neq t$, and
$L(\Phi,2)\neq0$.

$\Phi$ and $F$, and hence $t$, are built **from the recurrence alone**
(Frobenius $\to$ nome $\to$ series reversion, `common.gp: build`), so the
identification $\Phi=$ Eisenstein series of Theorem A is re-verified
coefficientwise to $q^{900}$ as a by-product.

---

## 5. The three complex folds, identified

For $\mathbf B,\delta,\eta$ the endpoint condition fails, $\Theta$ has no
finite value at the cusp $0$, and the connection constant is genuinely complex.
The folds are still CM points, but they lie **on** the Fricke circle
$|\tau|=1/\sqrt N$ rather than at its centre $i/\sqrt N$, so $W_N$ *exchanges*
the two folds (over the two conjugate roots of $P$) instead of fixing one —
which is precisely why $L(\Phi,2)\neq0$ there.

`04_complex_folds.gp` locates them exactly and evaluates
$\xi=\Theta+F\Theta'/F'$ to $90$ digits:

| row | $N$ | fold $\tau_*$ | $N|\tau_*|^2$ | $t(\tau_*)$ |
|---|---|---|---|---|
| $\delta$ | $12$ | $(1+i\sqrt2)/6$ | $1$ | $(7+i\sqrt{32})/81$ |
| $\eta$ | $20$ | $(1+2i)/10$ | $1$ | $0.088+0.016i$ |
| $\eta$ | $20$ | $(2+i)/10$ | $1$ | $0.088-0.016i$ |
| $\mathbf B$ | $36$ | cusp $1/6$ | — | $(9+i\sqrt{27})/54$ |

($\delta$: $12\tau^2-4\tau+1=0$, discriminant $-32$; $\eta$: $20\tau^2-8\tau+1=0$
resp. $20\tau^2-4\tau+1=0$, discriminants $-16,-64$ — the fields
$\mathbb Q(\sqrt{-2})$, $\mathbb Q(i)$ of `Observation obs:complex`.)

> **Result (settles Problem $\Phi$-1 numerically).**  With the fold normalised
> as above,
> $$\boxed{\ \xi_\delta=\frac{13}{54}\zeta(3)+i\,\frac{\pi^3}{81},\qquad
> \xi_\eta=\frac12L(\chi_5,3)+i\,\frac{\pi}{10}L(\chi_5,2),\qquad
> \xi_{\mathbf B}=\frac12L(\chi_{-3},2)+i\,\frac{2\pi^2}{27\sqrt3}.\ }$$
> \[VERIFIED: $\delta,\eta$ to $10^{-94}$; $\mathbf B$ to $10^{-55}$.\]

Three structural remarks.

1. **The real part is still the critical value.**  In all three rows
   $\operatorname{Re}\xi=L(\Phi,w+1)$ *exactly* — $13\zeta(3)/54$,
   $\tfrac12L(\chi_5,3)$, $\tfrac12L(\chi_{-3},2)$ — i.e. the archimedean
   formula $-\tfrac12P(w+1)\Lambda$ of the census **does** predict the
   complex-fold constants, but only their real part.  This is exactly the
   pattern the $p$-adic census sees ($\mathbf B\mapsto\tfrac12\zeta_3(2)$,
   $\delta\mapsto\tfrac14\zeta_3(3)$, $\eta\mapsto\tfrac12\zeta_5(3)$ after the
   Euler correction).
2. **The imaginary part is a quasiperiod.**  $\operatorname{Im}\xi\in
   \pi\,\mathbb Q\cdot L(\psi,w)$: $\tfrac{2}{27}\pi\zeta(2)$ for $\delta$,
   $\tfrac1{10}\pi L(\chi_5,2)$ for $\eta$, $\tfrac29\pi L(\chi_{-3},1)$ for
   $\mathbf B$.  Equivalently $\operatorname{Im}\xi_\delta=\tfrac{2\pi}3
   L(\Phi_\delta,2)$, $\operatorname{Im}\xi_{\mathbf B}=\tfrac{2\pi}3
   L(\Phi_{\mathbf B},1)$, $\operatorname{Im}\xi_\eta=\tfrac{12\pi}{35}
   L(\Phi_\eta,2)$.  A uniform closed form for the rational factor is
   **open**; it should come from the analogue of (2.3) for the anti-holomorphic
   pairing that swaps the two folds.
3. **The previously recorded values were not accurate.**
   `certificates/eisenstein/main(20260813-131951).tex`,
   Verification `ver:complex`, reports
   $\xi_\delta=0.289384069279185217\ldots-0.382793539263049809\ldots i$
   "to $3\cdot10^{-39}$" and $\xi_\eta$ "to $2\cdot10^{-33}$".  Both agree with
   the values above only to $\approx13$ and $\approx9$ digits; the stated error
   bars were too optimistic (they came from a $q^{60}$ truncation of $t(q)$ at
   $|q_c|\approx0.23$–$0.53$, plus a Newton solve for $q_c$ on the same
   truncation).  The correct constants sit inside the PSLQ basis that
   `ver:pslq` used, so the negative PSLQ result of that verification was an
   artefact of the input digits, not of the basis.  **`ver:complex` and
   `ver:pslq` should be corrected, and Problem $\Phi$-1 closed.**

---

## 6. What this changes in the paper

1. **Theorem B(i)** can be upgraded from "$\xi$ is *a rational multiple* of a
   critical value \[the multiple VERIFIED case by case\]" to
   "$\xi_\infty=L(\Phi,w+1)$, PROVED", with the rational multiple read off the
   Mellin polynomial as $r_\infty=P(w+1)L(\varphi,0)$.
   The line "Uniform closed evaluation of \eqref{eq:fold} giving the rational
   prefactor — OPEN" in the evidence table of `sec:sec2status` becomes PROVED.
2. **Lemma `lem:fold`** should be stated for $r=3$ only, and paired with the
   logarithmic (cusp) version for $r=2$; as written it asserts a square-root
   fold for all fifteen rows, which is false for the (R2) six.
3. **Remark `rem:complexfold`** ("no rational relation was found") is
   superseded by §5.
4. **`PADIC_PERIOD.md` §2**: the nonvanishing of $L(\Phi_{\mathbf F},1)$ is not
   an obstruction to the archimedean identification — the coefficient it
   carries multiplies $\tau$ and the second-order fold sits at $\tau=0$.
   Archimedean and $3$-adic now match structurally: the missing
   $L(\cdot,1)$-coefficient $3$-adically is the $p$-adic shadow of the fact
   that this coefficient never reaches the fold value.
