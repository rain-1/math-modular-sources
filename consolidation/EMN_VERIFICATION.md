# The Eskandari–Murty–Nemoto Catalan motive: verification of the ChatGPT note, and its price in our units

*Fable, 2026-08-23.  Scripts: `lattice/emn/` (`emnexact.py`, `emngen.py`, `01`–`25`).
Source under test: `chatgpt-research-archive/A Non-Eisenstein Catalan Programme from the
Eskandari–Murty–Nemoto Motive.md` (1185 lines).  Scoring framework:
`CATALAN_OBSTRUCTION.md`, `CDT_UNPACKED.md` §§1–4, `ADELIC_HOLONOMY.md` §2.6 (Theorem A),
`CATALAN_MU4.md` (pricing template; `lattice/catalan_mu4/calib.py` reproduces CDT's numbers),
`paper/sections/05_two_row.tex` (design rule), `MULTI_PRIME_LATTICE.md` (the conductor-12 slot).*

**No theorem is claimed.**  Everything below is either an exact computation, a proof written out,
or a measurement with its range stated.

---

## 0. Verdict in one paragraph

Every substantive claim of the note is **true**.  Of 13 checkable items, **12 verify exactly**
(several to far higher precision than the note asserts, and three of them are *sharpened* here to
exact laws), **1 is corrected** (§13's "conditional" function is unconditional: the hypothesis
$G\in\mathbf Q$ is not used), and the reference is confirmed with two bibliographic corrections
(arXiv:2510.20648; the Murty is **V. Kumar** Murty, not M. Ram Murty; Nemoto is a real coauthor).
The note's own negatives reproduce.  In our units:

* **the radial Padé family scores $\delta\approx0.34$** (Zudilin $\times\mathbf E$: $0.9025$) and its
  holonomic host is $\mathbf P^1$ minus $\{0,\infty\}$ and a **single irreducible cubic Galois
  orbit**, entry $\approx-24$: dead;
* **the EMN moving-period host is the first Catalan architecture in this project whose CDT bound is
  not vacuous.**  Its singular set is $\{1,\infty\}$ — *two* points — so the host is parabolic and
  the **entry obstruction is absent** (best entry $+0.096$ against $-0.077$ at level 8, $-0.563$ at
  $\mu_4$).  But the fold $z=1$ is the *only* finite singularity, so the CDT hypothesis mechanism has
  nothing to buy: removing the fold removes the host.  Best margin $\mathbf{-4.05}$;
* **the 3-adic pair does not close.**  $\Phi_3(3/2)=L_3(2,\chi_{12})$ is verified to **89** 3-adic
  digits and does fill the slot left open in `MULTI_PRIME_LATTICE.md` §8.7 — but the conductor-6
  decayer's $\sigma_3=3$ class is $\zeta_3(2)$, and $L_3(2,\chi_{12})$ and $\zeta_3(2)$ are
  3-adically $\mathbf Q$-independent (lindep at precision $3^{78}$ returns height $10^{12}$).
  There is no partner.  The EMN row that decays archimedeanly ($z=1/2$) has $\sigma_3=0$;
  the row with $\sigma_3=1$ ($z=3/2$) does not decay.  $\delta=0.19$, $F>0$.

The one genuinely new structural fact is §5 below: **the EMN period is the $\mu_4$ dilogarithm
module read at a regular point rather than at a cusp**, and that is exactly why its entry deficit
disappears.

---

## 1. Claim-by-claim verdict table

| # | claim (note) | verdict | evidence |
|---|---|---|---|
| (a) | $a_n=\frac1{2^{n+1}(n+1)}\sum_k\binom nk\frac1{2k+1}$ equals $\iint_\Delta(x^2+y^2)^n$ | **verified** | exact, $n\le60$ against the direct integral; `01_moments.gp` |
| (a) | $(n{+}3)(2n{+}5)a_{n+2}-(n{+}2)(3n{+}5)a_{n+1}+(n{+}1)^2a_n=0$ | **verified**, and **not minimal** | 0 failures $n\le203$; superseded by (b′) |
| (b) | $\bigl[(\theta{+}1)(2\theta{+}1)-z(\theta{+}1)(3\theta{+}2)+z^2(\theta{+}1)^2\bigr]\Phi=\tfrac12$ | **verified** | 0 failures to $z^{298}$; `02_ode.gp` |
| (b) | $z(z-2)H''+(z-1)H'=\frac1{2(z-1)}$, $H=z\Phi$ | **verified** | 0 failures to $z^{300}$ |
| (b′) | **sharpening (new):** the collapsed form is the *two-term inhomogeneous* recurrence $\ (m{+}1)(2m{+}1)a_m-m^2a_{m-1}=\tfrac12$ | **proved + verified** | 0 failures $m\le300$; the note's three-term relation is its first difference |
| (b) | homogeneous periods $1,\ S(z)=2\arcsin\sqrt{z/2}$ | **verified** | ODE residual $O(w^{28})$ in $z=2w^2$; `23_singset.gp` |
| (c) | $\mathscr H'(S)=\tfrac12\log\tan(\tfrac\pi4+\tfrac S2)$; Fourier form | **verified** | 25 digits at 3 points; `03_values.py` |
| (c) | $H(z)+H(2-z)=2G$ | **verified, and proved** | 40 digits; one line from the Fourier series |
| (c) | $H(3/2)=\tfrac53G$, $G=\tfrac32\Phi(1/2)$ | **verified, and proved** | 50 digits; the $\tfrac{10}9$ falls out as $1-\chi_{-4}(3)3^{-2}$ |
| (c) | $R_3(z)=z(3-2z)^2=1-\cos 3S$; $\ \mathscr H(3S)+3\sum_j\mathscr H(S{+}\tfrac{2\pi j}3)=10G$ | **verified, and proved** | algebra + Fourier; 35 digits |
| (d) | $\operatorname{PV}\Phi(3/2)=\tfrac{10}9G$ | **verified** | **50 digits**, by the radial principal-value integral $-\frac1{2z}\int_0^{\pi/2}\log\bigl|1-\frac z{(\cos t+\sin t)^2}\bigr|dt$ |
| (§6) | $G=\tfrac32\operatorname{Ti}_2(2-\sqrt3)+\tfrac\pi8\log(2+\sqrt3)$ | **verified** | 50 digits |
| (e) | $\sum a_n(3/2)^n$ converges $3$-adically and equals $L_3(2,\chi_{12})$ | **verified** | **89 3-adic digits** (note claims 30); `04_padic3.gp`, `05_kl_sanity.gp` |
| (f) | $b_{m,t}=4^{-m}\binom m{m/2}\binom mt$ (even $t$) | **verified**, general form carries $(-1)^t$ | 140 exact cases $m\le18$, $-4\le t\le m$; `07_check_f.py` |
| (f) | $I_{4,2}=-\tfrac{49}{384}+\tfrac9{64}G$ and the pole-lowering recurrence | **verified** | exact; matches arXiv:2510.20648 Thm 1.1.1 verbatim; 45 cases |
| (g) | $256^nb_n=\binom{4n}{2n}^2\in\mathbf Z$; $|I_n|=64^{-n+o(n)}$ | **verified** | $n\le80$; $-\log|I_n|/n=4.269$ at $n=80$, $\log64=4.1589$ |
| (g) | common 2nd-order system, characteristic roots $256,4$ | **verified exactly** | order 2, coefficient degree 10, nullspace dim 1, char. poly. $21504\lambda^2-5591040\lambda+22020096=21504(\lambda-256)(\lambda-4)$; `12_recpoly.gp` |
| (g) | "large odd denominators" of $a_n$ | **verified, made exact** | $\operatorname{den}(a_n)\mid 2^{8n}d_{4n}d_{2n}$ (0 exceptions $n\le48$); $v_2=8n-O(\log n)$, deviation in $[-9,-1]$ for $n\le48$ (equivalently $256^na_n$ has **odd** denominator); the tightest tested shape ($2^{8n}d_{4n}^2$ also divides but is not tight; $2^{8n}d_{4n}d_n$ and $2^{8n}d_{4n}$ both **fail**); total rate $\approx10.4$ |
| (h) | $256^nb_n\in\mathbf Z$ for the radial Padé family | **verified, and proved** | closed form $\binom{4n}{2n}\sum_k\binom{2n}k(-4)^k\binom{4n}{2n-2k}$; $n\le40$ |
| (h) | $0<256^nL_n\le G(16/27)^n$ | **verified** | $n\le16$ exactly; `13_pade.py` |
| (i) | $-\log M(10,3)=18.9221280514\ldots$; $M(4,1)=1/432$ | **verified** | $18.92212805141808$ |
| (i) | $2^{26n+O(\log n)}b_n\in\mathbf Z$; primary margin $0.9003013569$ | **verified** | $v_2(\operatorname{den}b_n)/n = 22,24,24.3,25,25,25.2$ for $n\le6$; odd part $=1$ exactly |
| (i) | the primitive score is negative by $>20n$ | **verified** | primitive score $\approx-22$ per $n$; gap to $+0.90$ is $\approx23n$ |
| (j) | $\tilde K=H(R_3)+3H$ has rational coefficients for $n\ge1$ and $v_2=n+O(\log n)$ | **verified and sharpened** | $v_2([z^n]\tilde K)=n+\mathbf 1_{2\mid n}$ **exactly**, 0 exceptions $1\le n\le2000$ |
| (j) | nearest singularity at $1-\sqrt3/2$ | **verified** | $c_n/c_{n+1}\to0.134109$ vs $1-\sqrt3/2=0.1339746$ at $n=2000$ |
| (j) | rigidity: the only full-orientation pattern is $(1,3,3,3)$ | **verified, and proved** | local log-coefficients $-cR_3'(z^*)$ with $R_3'(1)=-3$, $R_3'(1\pm\tfrac{\sqrt3}2)=+6$, and $\zeta'=-2$ at the conjugate branches: $\lambda_1=\lambda_2=\lambda_3=3$ forced |
| (§13) | "*if $G\in\mathbf Q$*, $K=\tilde K-4G$ has rational coefficients and is holomorphic across $z=1$" | **CORRECTED** | $\tilde K\in\mathbf Q[[z]]$ and is holomorphic at $z=1$ **unconditionally**; the hypothesis affects only the constant term.  $\tilde K$ is an *unconditional* function, not a conditional generator |
| ref | "Eskandari–Murty–Nemoto 2025, 2-dim'l mixed motive with period $G$" | **verified with corrections** | arXiv:2510.20648, *Mixed motives and linear forms in the Catalan constant*, Payman Eskandari, **V. Kumar** Murty, Yusuke Nemoto; 23 Oct 2025, v2 26 Jan 2026.  Period matrix $\begin{pmatrix}1&G\\0&\pi^2i\end{pmatrix}$.  §1.2 gives $G=\int_\Delta\frac{dx\,dy}{1-x^2-y^2}$, $\Delta=\{x,y,1-x-y\ge0\}$ **verbatim** |

Counts: **12 verified / 0 refuted / 1 corrected**, plus 3 claims sharpened to exact laws
((b′), (g) denominator type, (j) $v_2$) and 5 claims upgraded from "computational evidence" to
proof ((c) three times, (f), (h), (j)-rigidity).

---

## 2. What was proved rather than measured

### 2.1 The two-term recurrence and the collapsed ODE are the same statement

$H=\sum_{m\ge0}a_mz^{m+1}$ gives, coefficientwise,
$$z(z-2)H''+(z-1)H'=\frac1{2(z-1)}\iff m^2a_{m-1}-(m+1)(2m+1)a_m=-\tfrac12\quad(m\ge0).$$
This is *order one* in $m$; the note's three-term relation follows by subtracting consecutive
instances.  It also gives $a_n$ in $O(n)$ operations, which is what makes §§(g)–(j) feasible.

### 2.2 $H(z)+H(2-z)=2G$, $H(3/2)=\tfrac53G$, and the cubic distribution

With $\mathscr H(S)=\sum_{m\ge0}c_m\bigl(1-\cos((2m{+}1)S)\bigr)$, $c_m=\frac{(-1)^m}{(2m+1)^2}$,
$\sum c_m=G$:

* $\cos((2m{+}1)(\pi-S))=-\cos((2m{+}1)S)$ gives $\mathscr H(S)+\mathscr H(\pi-S)=2\sum c_m=2G$, and
  $1-\cos(\pi-S)=2-z$.  **QED**
* At $S=2\pi/3$: $\cos((2m{+}1)\tfrac{2\pi}3)=1$ if $3\mid 2m{+}1$ and $-\tfrac12$ otherwise; for
  $2m{+}1=3(2j{+}1)$ one has $m=3j+1$, $c_m=-\tfrac19c_j$; hence
  $\mathscr H(2\pi/3)=\tfrac32\bigl(G+\tfrac G9\bigr)=\tfrac53G$.  **The $\tfrac{10}9$ is
  $1-\chi_{-4}(3)3^{-2}$, i.e. the Euler factor at 3, read off the cubic geometry.**  QED
* $\mathscr H(3S)+3\sum_{j=0}^2\mathscr H(S+\tfrac{2\pi j}3)$: the $j$-sum kills every frequency
  with $3\nmid 2m{+}1$; what survives is $-3\cdot(-\tfrac19)X$ with
  $X=\sum_mc_m\cos(3(2m{+}1)S)$, and $\mathscr H(3S)=G-X$.  Total $G-X+9G+X=10G$.  QED
* $R_3(z)=1-\cos3S$ with $c=\cos S=1-z$: $1-(4c^3-3c)=4z^3-12z^2+9z=z(3-2z)^2$.  QED

### 2.3 The higher-pole periods in closed form (this is the engine for (f)–(i))

`lattice/emn/emnexact.py`, `emngen.py`.  In polar coordinates with $u=r^2$ and
$p=\sin2\theta$ (so $\int_0^{\pi/2}g(\sin2\theta)\,d\theta=\int_0^{\pi/2}g(\sin\psi)\,d\psi$),
$$T(p,q,t)=\iint_\Delta\frac{x^py^q}{(1-x^2-y^2)^{t+1}}
=\int_0^{\pi/2}\!\cos^p\!\theta\sin^q\!\theta\;\Bigl[\tfrac12\!\int_0^{W}\!u^{\frac{p+q}2}(1-u)^{-t-1}du\Bigr]d\theta,
\quad W=\tfrac1{1+\sin2\theta}.$$
Expanding $u^{mm}=\bigl(1-(1-u)\bigr)^{mm}$, the inner integral is elementary except for the
$i=t$ term, which produces $-\log\frac{s}{1+s}$.  Everything then reduces to four atoms in
$\mathbf Q\oplus\mathbf Q\pi\oplus\mathbf Q\log2\oplus\mathbf Q\pi\log2\oplus\mathbf QG$:
$$C_k=\!\int_0^{\pi/2}\!\!\sin^k,\quad P_r=\!\int_0^{\pi/2}\!\!(1+\sin)^r,\quad
\Lambda_k=\!\int_0^{\pi/2}\!\!\sin^k\log\sin,\quad
\Sigma_k=\!\int_0^{\pi/2}\!\!\sin^k\log(1+\sin),$$
with $P_{-q}\in\mathbf Q$ (the substitution $t=\tan\frac\psi2$ leaves no $u^{-1}$ term), and the
key recurrence — obtained from $\frac d{d\psi}\bigl[\cos\psi\sin^{k-1}\psi\bigr]$ and
$\frac{\cos^2}{1+\sin}=1-\sin$ —
$$\boxed{\ k\,\Sigma_k=(k-1)\Sigma_{k-2}+C_{k-1}-C_k\ },\qquad
\Sigma_0=2G-\tfrac\pi2\log2,\quad \Sigma_1=\tfrac\pi2-1 .$$
Only $\Sigma_0$ carries $G$, so the $G$-part of $\Sigma_{2k}$ is $2\binom{2k}k4^{-k}$, and
$$\boxed{\ b_{m,t}=(-1)^t\,4^{-m}\binom m{m/2}\binom mt\ }$$
**falls out as a theorem**, with $b_{m,t}=0$ for $t<0$.  The note's version is the case $t$ even.
The $\pi$, $\log2$ and $\pi\log2$ components must cancel (EMN Thm 1.1.1) and do — this is asserted
in the code and is an independent check on every value.
`06_poles.py` re-derives 21 of these values by 140-digit PSLQ: **agreement in all cases**.

### 2.4 $256^n b_n^{\rm Padé}\in\mathbf Z$

$(1-4g^2)^{2n}=\sum_k\binom{2n}k(-4)^kg^{2k}$ with $g=1-x^2-y^2$ gives
$L_n=\sum_k\binom{2n}k(-4)^kI_{4n,\,2n-2k}$, hence by §2.3
$$256^nb_n^{\rm Padé}=\binom{4n}{2n}\sum_{k=0}^{2n}\binom{2n}k(-4)^k\binom{4n}{2n-2k}\in\mathbf Z .$$
Verified against the exact integral for $n\le16$ and evaluated to $n=160$.

### 2.5 §13 rigidity

Near $w=1$, $H(w)=\text{hol}+c\,(1-w)\log(1-w)$ with $c=\tfrac12$ (from
$H'(z)\sim-\tfrac12\log(1-z)$).  For $u$ with $u(z^*)=1$ the $(z-z^*)\log(z-z^*)$ coefficient of
$H(u(z))$ is $-c\,u'(z^*)$.  With $R_3'=3(3-2z)(1-2z)$:
$R_3'(1)=-3$, $R_3'(1\pm\tfrac{\sqrt3}2)=+6$; and the conjugate branches
$\zeta_\pm(z)=\frac{3-z}2\pm\frac{\sqrt{3z(2-z)}}2$ of $R_3^{-1}R_3$ have $\zeta'=-2$ at
$z=1\mp\tfrac{\sqrt3}2$.  Cancelling at all three preimages of the fold forces
$\lambda_1=3$ and $\lambda_2=\lambda_3=3$: the pattern $(1,3,3,3)$, i.e. the full cubic trace,
which is the constant $10G$.  **There is no nonconstant full-orientation generator.**  QED

---

## 3. The corrected item: §13 is unconditional

$\tilde K(z)=H(R_3(z))+3H(z)$ has **rational** Taylor coefficients unconditionally
($H\in\mathbf Q[[z]]$, $R_3\in\mathbf Q[z]$) and is holomorphic across $z=1$ unconditionally
(§2.5 with $\lambda_1=3$).  Subtracting $4G=\tilde K(1)$ changes only the constant term.  So the
sentence "*if $G\in\mathbf Q$, $K=\tilde K-4G$ has rational coefficients and is holomorphic across
the fold*" reverses the logic: **the hypothesis is not used**.

This matters, because the note then reads off "$\varsigma_2=1$ for the conditional generator".  The
slope is real and is now known exactly:
$$\boxed{\ v_2\bigl([z^n]\tilde K\bigr)=n+\mathbf 1_{2\mid n}\quad\text{for every }1\le n\le2000\ }$$
(0 exceptions; `15_cubic2000.gp`).  But it is the slope of an **unconditional** function.  In the
adelic bound (`ADELIC_HOLONOMY.md` §2.6) an unconditional slope-1 function is a *good* thing — it
contributes to $\sum_p\log R_p$ — and `CATALAN_OBSTRUCTION.md` (2.3) is untouched: it is the
*conditional* function whose slope is forced to $0$ by Calegari's theorem, and §13 supplies no
conditional function at all.  Background data: $v_2(a_n)=-1$ **exactly** for every $n\le600$, so
$H$ itself has 2-adic slope $0$ in $z$ and $\tfrac12$ in $\zeta=\sqrt{z/2}$
($v_2(2^{n+1}a_n)=n$ exactly, $n\le300$).

---

## 4. Scoring 1 — the radial Padé family

Family: $\displaystyle L_n=\iint_\Delta\frac{x^{4n}y^{4n}(1-4g^2)^{2n}}{g^{2n+1}}$, $g=1-x^2-y^2$.

| quantity | value | source |
|---|---|---|
| archimedean decay of $256^nL_n$ | $\le\log\frac{27}{16}=0.5232481$ (bound verified $n\le16$; measured $0.964$ at $n=16$, still falling) | `13_pade.py` |
| $256^nb_n\in\mathbf Z$, growth | $\log|256^nb_n|/n\to6.19812$ (dominant pair of $729\lambda^3+124064\lambda^2+176304384\lambda-47775744$) | `21_pade_u.gp` |
| $\operatorname{den}(a_n)$ | $\log/n\approx19.7$–$20.4$; $v_2/n\approx7.7$–$7.9$; odd part $\approx14$; primes up to $\approx12n$ | `13_pade.py` |
| cleared form $D_n=\operatorname{lcm}(\operatorname{den}a_n,256^n)$ | $\log D_n/n\approx19.8$, $\log|D_nL_n|/n\approx+13.3$ | |
| $H=\log|q_n|/n$, $F=\log|\ell_n|/n$ | $20.15,\ 13.27$ at $n=16$ | |
| **worthiness $\delta=1-F/H$** | $\mathbf{0.341}$ (values $0.34$–$0.38$ over $n=4$–$16$) | |
| $\delta$ in the naive "decay / denominator rate" reading | $0.5232/19.74=\mathbf{0.0265}$ | |
| Zudilin $\times\mathbf E$ for comparison | $\mathbf{0.9025266029}$ | `05_two_row.tex` |

**The holonomic system.**  The $b$-side $U_n=256^nb_n$ satisfies a minimal recurrence of **order 3**
with degree-16 polynomial coefficients; its characteristic polynomial
$729\lambda^3+124064\lambda^2+176304384\lambda-47775744$ is **irreducible over $\mathbf Q$**, with
roots $0.27093\ldots$ and a complex pair of modulus $491.83$.  So the generating function's
singular set is $\{0,\infty\}$ together with a **single irreducible cubic Galois orbit** at
$|x|=0.002033$ and $x=3.6910$.

That is the `HADAMARD_HOST.md` / Proposition-G configuration verbatim: a singular set that is one
Galois orbit admits **no pure submodule**, so $\tau^\flat\ge\sigma_m$, and the conformal ceiling is
set by the nearest singularity, $\log(4\cdot0.002033)=-4.81$.  With $\sigma\gtrsim19.7$ per unit of
$n$ (of which $7.8\log2=5.4$ is a genuine *geometric* 2-adic denominator, not an lcm layer), the
entry number is
$$\text{entry}\ \approx\ -4.81-19.7\ \approx\ \mathbf{-24.5}.$$
No contour, symmetrisation or adelic twist moves that.  **The radial Padé family is dead by two
independent mechanisms** (Galois-orbit singular set, and a denominator rate 38× the decay).

---

## 5. Scoring 2 — the EMN moving-period host.  *This is where the note is genuinely new.*

### 5.1 The geometry, corrected and completed

The note's ODE has singular points $\{0,1,2,\infty\}$ with the fold at $z=1$.  But the *solution*
$H$ is regular at two of them:

* $z=0$: $H=\sum a_nz^{n+1}\in\mathbf Q[[z]]$ — exponent 1, not $\tfrac12$;
* $z=2$: $H(z)=2G-H(2-z)$ (§2.2), and $H(2-z)$ is the Taylor series — **regular**;
* $z=1$: $H'(z)\sim-\tfrac12\log(1-z)$, a genuine branch point.  $a_n\sim\frac1{2n^2}$ with
  $n^2a_n\uparrow\tfrac12$ monotonically (checked to $n=800$): a *single* dominant singularity, no
  oscillation, local type $(1-z)\log(1-z)$;
* $z=\infty$: $\mathscr H(S)\sim\frac{i\pi}4S$, a log singularity.

$$\boxed{\ \text{Sing}(H)=\{1,\infty\}\ }$$
$H'$ too ($H'(z)=H'(2-z)$), so the whole $\theta$-orbit lives there.  **The host is $\mathbf P^1$
minus two points.**  Every other Catalan host in the ledger has 4–6 singular points.

**Why.**  Put $q=e^{iS}$ in $\mathscr H(S)=\sum\frac{(-1)^m}{(2m+1)^2}(1-\cos((2m{+}1)S))$:
$$\boxed{\ \mathscr H(S)=G-\tfrac12\bigl(\operatorname{Ti}_2(q)+\operatorname{Ti}_2(q^{-1})\bigr),\qquad
z=1-\tfrac{q+q^{-1}}2=-\tfrac{(q-1)^2}{2q}\ }$$
so **the EMN moving period is the $\mu_4$ dilogarithm module of `CATALAN_MU4.md` (host A), pushed
down the involution $q\mapsto1/q$, and expanded at $q=1$ instead of at $q=0$.**  Under $q\mapsto z$:
$q=\pm i\mapsto z=1$ (the two folds merge), $q\in\{0,\infty\}\mapsto z=\infty$,
and the *fixed points* $q=\pm1$ become $z=0,2$ — which is exactly why they are regular for $H$.
`CATALAN_MU4.md` §4 states that host A "admits no involution swapping an outer singularity with
$\infty$"; $q\mapsto1/q$ **is** such an involution, but it does not fix host A's base point $q=0$.
EMN's base point $z=0$ is $q=1$, which it does fix.  *The whole gain of the EMN construction is a
change of base point: from a cusp of the $\mu_4$ host to a regular point of it.*

### 5.2 The denominators — exact

Measured to $n=400$ with **zero exceptions**:
$$\boxed{\ \operatorname{den}(a_n)\ \Bigm|\ 2\,(n+1)\,[1,\dots,2n+1]\ },\qquad v_2(a_n)=-1\ \text{exactly}.$$
So $H$ has CDT type $n\,[1..2n]$: one column with $b=2$, $e=1$, $\sigma=2$ — *the same $\sigma$ as
CDT's own conditional function.*  The true rate is smaller, $\log\operatorname{den}(a_n)/n=1.3946$
at $n=400$ and still falling (against $\log d_{2n+1}/n=1.98$), plausibly $\to\log4$; but the type is
what $\tau$ charges, and $2\,[1..2n+1]$ alone fails for 355 of 400 values, so $b=2$ cannot be
lowered on the $z$-line.  In the coordinate $\zeta=\sqrt{z/2}$ the index doubles and the type
becomes $N[1..N]$, i.e. $b=1$; and $v_2$ of the $\zeta^{2n+2}$-coefficient is exactly $n$, giving a
genuine 2-adic radius $R_2=\sqrt2$.

### 5.3 Conformal ceilings (`18_geom.py`, exact via $\lambda'(\tau)=i\pi\lambda(1-\lambda)\theta_3^4$)

| geometry | $\sup|\varphi'(0)|$ | $\log$ |
|---|---|---|
| $\mathbf C\setminus\{\pm1\}$ at $0$ (= the $S$-plane at $S=0$, since $\sin$ is an unbranched covering $\mathbf C\setminus(\tfrac\pi2+\pi\mathbf Z)\to\mathbf C\setminus\{\pm1\}$) | $4.376879230452953$ | $1.47634$ |
| $\mathbf C\setminus\{\pm1/\sqrt2\}$ at $0$ ($\zeta$-line, folds only) | $3.094920984287841$ | $1.12976$ |
| $\mathbf C\setminus\{1,2\}$ at $0$ | $8.753758460905907$ | $2.16948$ |
| $\mathbf C\setminus\{1\}$ at $0$ | $+\infty$ (**parabolic**) | — |
| Koebe cap for *univalent* $\varphi$, omitting $\{1\}$ | $4$ | $1.38629$ |

The number $4.376879230452953$ is the one `CATALAN_MU4.md` §4 quotes for
$\mathbf C\setminus\{\pm i\}$ (same by rotation) and where it says the extra condition
$\varphi^{-1}(0)=\{0\}$ costs $\log(4.3769/4)=0.090$.  **On the EMN host that condition is not
imposed, because $z=0$ is not a cusp** — the note recovers exactly that $0.090$.

### 5.4 Entry — the obstruction is absent

$\mathbf C\setminus\{1\}$ is not hyperbolic: $\sup|\varphi'(0)|=\infty$ (take
$\varphi_R(z)=1-e^{-Rz}$, $\varphi_R'(0)=R$).  So there is **no conformal ceiling to compare with
$\tau$**, and the entry condition $\log|\varphi'(0)|>\tau$ can always be met.  This is the first
Catalan architecture in the ledger for which the CDT bound is not vacuous.

| architecture | ceiling | $\tau$ | entry |
|---|---|---|---|
| CDT's own $L(2,\chi_{-3})$, symmetrised | $\log256$ | $4.2355$ | $+1.3097$ |
| Catalan level 8, symmetrised (best modular) | $\log64$ | $4.2355$ | $-0.0766$ |
| $\mu_4$ host A | $\log4$ | $1.9490$ | $-0.5627$ |
| $\mu_4$ host B / C | $\log16$ / $\log256$ | $3.898$ / $7.796$ | $-1.125$ / $-2.251$ |
| **EMN $z$-line, $\{1,\mathrm{Li}_1,H,\theta H\}$** | unbounded | $2.3021$ | $\mathbf{+0.0958}$ (at $R=e^{2.398}$) |
| EMN $\zeta$-line, $\{1,H\}$, adelic ($R_2=\sqrt2$, $m=2$) | $\log3.0949$ | $1.5$ | $-0.197$ |
| EMN $z$-line, $\{1,H\}$, Koebe (univalent only) | $\log4$ | $2.25$ | $-0.864$ |

### 5.5 Margin — where it dies, and why the hypothesis cannot help

With $\varphi_R(z)=1-e^{-Rz}$ and the Bost–Charles numerator computed numerically
(`24_bc.py`, `25_scan.py`; the integrator reproduces $\mathrm{BC}(\rho z)=\log\rho$ to 8 decimals):

| inventory | $m$ | $\tau$ | best entry | $\mathrm{BC}$ | bound $m\le$ | **margin** |
|---|---|---|---|---|---|---|
| $\{1,H,\theta H\}$ | 3 | $2.6111$ | $+0.1615$ | $6.46$ | $15.02$ | $-5.98$ |
| $\{1,\mathrm{Li}_1,H,\theta H\}$ | 4 | $2.3021$ | $+0.0958$ | $4.43$ | $11.04$ | $\mathbf{-4.05}$ |
| $\{1,\mathrm{Li}_1,\mathrm{Li}_2,H,\theta H\}$ | 5 | $2.4867$ | $+0.0783$ | $5.24$ | $13.28$ | $-4.85$ |
| $\;+\;\mathrm{Li}_3..\mathrm{Li}_{12}$ | 6–15 | $2.70$–$4.44$ | — | — | $16$–$93$ | $-6.0$ to $-35.4$ |

(margin $=m(\log|\varphi'(0)|-\tau)-\mathrm{BC}$; $>0$ would be a contradiction.)  Adding functions
makes it worse: on a parabolic host $\mathrm{BC}(\varphi_R)$ grows like $0.405R$ while
$\log|\varphi'(0)|=\log R$, so the bound rises faster than $m$.  Best margin $-4.05$
(vs $-14.04$ for $\mu_4$ host A, $-5.8$ for the modular hosts) — **the best Catalan margin in the
project, and the only one that is a real number rather than bookkeeping.**

**But the architecture cannot prove anything, for a structural reason.**  Every function above is
*unconditional*, so a positive margin is impossible a priori (it would say $1,\mathrm{Li}_1,H,\theta H$
are $\mathbf Q(z)$-dependent, which is false).  A CDT irrationality proof needs a **conditional**
function, and on this host there is none to be had:

1. **$H$ itself is not conditional** (the note's own §14.1, confirmed): the log coefficient of
   $H'$ at $z=1$ is $-\tfrac12$, a *rational* number, so no multiple of $G$ removes it.
2. **The conditional function must be manufactured**, exactly as on the $\mu_4$ host:
   $Y(z)=\int_0^z\frac{H(w)}{1-w}dw$ has $\Delta_{z=1}Y=-2\pi iG$ (since $H(1)=G$), so
   $Y_{\rm cond}=bY-aA$, $A=-\log(1-z)$, is fold-regular under $G=a/b$.
3. **And it buys nothing.**  In CDT the gain is that $\varphi$ may avoid *one point fewer*.  Here
   $z=1$ is the **only finite singularity of the whole module**; if every function were
   conditionally regular there, the module would consist of entire $G$-functions, i.e. polynomials.
   So $\varphi$ must still avoid $z=1$, the host is unchanged, and the hypothesis contributes
   exactly $+1$ to $m$ against a bound that rises by $\approx2.2$.

$$\boxed{\ \text{EMN removes the entry obstruction by removing the outer singularity — and with it the fold mechanism.}\ }$$

This is a new entry for `CATALAN_OBSTRUCTION.md` §3: it is neither (E1) (irrational second
singularity) nor (E2) nor (E3) as stated, but a fourth possibility — *no second singularity at
all* — which turns out to be self-defeating.  The width law
$\text{entry}\le\frac1w\log16+\log|t_2|-k$ of (2.2) does **not** bind here: it is a statement about
placements on the Beauville surface $I_8I_2I_1I_1$, and this host is not a modular curve, has no
$t_2$, and its ceiling is $+\infty$ rather than $16|t_2|$.  The trigonometric uniformisation
$z=1-\cos S$ genuinely gives something different — it just gives the wrong different thing.

---

## 6. Scoring 3 — the 3-adic pair

### 6.1 The identity, and why $\chi_{12}$ is the right character

At $p=3$ the Teichmüller character is $\omega_3=\chi_{-3}$, so
$\chi_{12}=\chi_{-4}\chi_{-3}=\chi_{-4}\,\omega_3$ is the **even** character whose Kubota–Leopoldt
$L$-function is the 3-adic avatar of $L(s,\chi_{-4})$.  (Dually, at $p=2$, $\omega_2=\chi_{-4}$ and
$\chi_{12}=\chi_{-3}\omega_2$, which is why `MULTI_PRIME_LATTICE.md` writes the *2*-adic avatar of
the $\chi_{-3}$ class as $L_2(2,\chi_{12})$.  Same $L$-function, two primes, two classes.)

**Verified.**  $\sum_{n\le260}a_n(3/2)^n$ and $L_3(2,\chi_{12})$ computed independently
(Washington Thm 5.11, $F=12$, 140 Bernoulli terms) agree to **$O(3^{89})$**, i.e. 89 3-adic digits;
the Kubota–Leopoldt implementation is validated separately at $s=0,-1,\dots,-5$ against
$-\bigl(1-\chi\omega^{-n}(3)3^{n-1}\bigr)B_{n,\chi\omega^{-n}}/n$ to $O(3^{59})$
(`04_padic3.gp`, `05_kl_sanity.gp`).  Convergence: $\min_nv_3(a_n)=-9$ for $n\le200$, so
$v_3(a_n(3/2)^n)\ge n-O(\log n)$.

This **does** fill the slot `MULTI_PRIME_LATTICE.md` §8 item 7 left open ("*whether a conductor-12
analogue exists for $L(2,\chi_{-4})$ … Building the missing 3-adic Catalan row is the interesting
half*").  The EMN row at $z=3/2$ is a Catalan-class object with $\sigma_3=1$ per unit of $N$
($v_3(s_N-L_3)=21,36,57,79$ at $N=20,40,60,80$) — the census of §1.2 there, "no modular Catalan row
has $\sigma_3>0$", is not contradicted (this is not a modular row) but it is *side-stepped*.

### 6.2 Same 3-adic class as the conductor-6 decayer?  **No.**

The $\sigma_3=3$ side of `MULTI_PRIME_LATTICE.md` is the conductor-6 well-poised row, whose period
is $L(2,\chi_{-3})$ and whose 3-adic limit is $\zeta_3(2)$
(`lattice/multi_prime/02_decayer_padic.gp`: "$p=3$: $\mathcal E_3(2)=1\Rightarrow\xi_3=\zeta_3(2)=2\xi_3(\mathbf C)$").
Computing $\zeta_3(2)=L_3(2,\mathbf 1)$ with the same machinery:

* no $p/q$ with $|p|,q\le40$ gives $v_3\bigl(L_3(2,\chi_{12})-\tfrac pq\zeta_3(2)\bigr)>40$;
* `lindep([L_3(2,\chi_{12}),\zeta_3(2),1])` at precision $3^{78}$ returns
  $(1459115229107,\,-1946436351699,\,-1526619413247)$ — height $\sim10^{12}$, i.e. **no relation**.

$$\boxed{\ L_3(2,\chi_{12})\ \text{and}\ \zeta_3(2)\ \text{are 3-adically }\mathbf Q\text{-independent.}}$$

So the two rows carry different classes at $p=3$ and cannot be bridged: the cross-determinant
divisor $T_n$ of Theorem E exists only when both rows have the *same* $p$-adic limit.

### 6.3 The EMN rows on their own: design-rule numbers

| row | archimedean | $p$-adic | denominators |
|---|---|---|---|
| $t_N=\sum_{n\le N}a_n2^{-n}\to\tfrac23G$ | $-\log|\tfrac23G-t_N|/N\to\log2=0.693$ (measured $0.793$ at $N=100$) | $\sigma_2<0$ (terms have $v_2=-1-n$, no 2-adic limit); $\sigma_3=0$ | $\log\operatorname{den}/N=3.59$, $\kappa_2=1$ |
| $s_N=\sum_{n\le N}a_n(3/2)^n\to L_3(2,\chi_{12})$ | **grows** at $\log\tfrac32=0.405$; no archimedean limit | $\sigma_3=1$ | $\log\operatorname{den}/N=3.55$, $\kappa_2=1$ |

For the archimedean row alone: $H=\log|q_N|/N=3.59$, $F=\log|\ell_N|/N=3.59-0.693=2.90$,
$$\boxed{\ \delta=1-F/H=\mathbf{0.192},\qquad F=+2.90>0\ }$$
— the linear form **grows**; there is no Diophantine content, only a worthiness statement.
(For comparison: Zudilin $\times\mathbf E$ gives $\delta=0.9025$; the conductor-6 multi-prime
best is $0.94$.)  The design rule
$\delta>1\iff\log\Lambda_{\rm dec}>k\max(r,1)+r\log\rho_2^{\rm eng}$ cannot even be applied,
because the two EMN rows are **not a bridgeable pair**: the decayer has $\sigma_3=0$ and the
$\sigma_3=1$ row has no archimedean decay.  The 3-adic resource is real but unreachable.

**What would close it:** a single sequence pair $(Q_n,P_n)$ with $Q_nG-P_n\to0$ archimedeanly
*and* $P_n/Q_n\to\tfrac9{10}L_3(2,\chi_{12})$ 3-adically.  The EMN motive supplies both limits but
on two different arguments of the same generating function, which is not the same thing.

---

## 7. What is new, what is not

**New (and correct):**
1. The moving period $\Phi(z)$, its order-2 inhomogeneous ODE, and the two-term recurrence.
2. The cubic distribution law and the geometric explanation of $\tfrac{10}9$ as $1-\chi_{-4}(3)3^{-2}$.
3. $b_{m,t}=(-1)^t4^{-m}\binom m{m/2}\binom mt$ — a genuinely useful closed form, here promoted to a proof.
4. $\Phi_3(3/2)=L_3(2,\chi_{12})$: the first non-Eisenstein 3-adic realisation of the Catalan class in this project, verified to 89 digits.
5. The central family with characteristic roots exactly $256,4$.
6. $256^nL_n\le G(16/27)^n$: a real free-integration phenomenon.
7. **(not in the note)** the identification $\mathscr H(S)=G-\tfrac12(\operatorname{Ti}_2(q)+\operatorname{Ti}_2(q^{-1}))$, $q=e^{iS}$, which places the EMN host inside the $\mu_4$ tower as the $q\mapsto1/q$ quotient expanded at $q=1$.
8. **(not in the note)** $\operatorname{Sing}(H)=\{1,\infty\}$, the resulting parabolic host, and the fact that this is the first Catalan architecture with a non-vacuous CDT bound.

**Not new / not usable:**
* §13's "conditional" generator is unconditional (§3).
* Every scored family is far from $\delta=1$: radial Padé $0.34$, EMN archimedean row $0.19$, saddle-killer primitive score $\approx-22$.
* The 3-adic class is real but has no partner (§6.2).

## 8. Reproduction

```
cd lattice/emn
gp -q 01_moments.gp 02_ode.gp            # (a),(b)
python3 03_values.py                     # (b),(c),(d),(6)  50 digits
gp -q 04_padic3.gp 05_kl_sanity.gp       # (e)  89 3-adic digits + KL validation
python3 06_poles.py 140                  # (f) by 140-digit PSLQ
python3 07_check_f.py                    # (f) exactly, + pole lowering
python3 08_central.py 26                 # (g)
python3 09_central_analysis.py 80        # (g) decay + denominators
gp -q 11_fitrec.gp 12_recpoly.gp         # (g) the order-2 system, roots 256 and 4
python3 13_pade.py 16                    # (h)
gp -q 14_cubic.gp 15_cubic2000.gp        # (j) to n=2000
python3 16_saddle.py 6                   # (i)
gp -q 17_denoms.gp 23_singset.gp         # denominator type, singular set
python3 18_geom.py 19_tau_entry.py       # conformal radii, tau, entry
gp -q 20_pade_rec.gp 21_pade_u.gp        # Pade holonomic system
gp -q 22_pair3adic.gp                    # the 3-adic pair
python3 24_bc.py 25_scan.py              # Bost-Charles numerators and margins
```
Libraries: `emnexact.py` (exact $I_{m,t}$), `emngen.py` (exact $T(p,q,t)$, $p\ne q$).
Both assert that the $\pi$, $\log2$ and $\pi\log2$ components cancel, which is an independent check
of EMN Theorem 1.1.1 on every value computed.

**Not done.**  No slit or gobble contour was designed for the EMN host (only the exponential family
$\varphi_R=1-e^{-Rz}$ and concentric discs), so the margin $-4.05$ is an upper bound on the
difficulty, not a sharp one; the $\mathbf Q(z)$-independence of $\{1,\mathrm{Li}_k,H,\theta H\}$ is
taken on structural grounds and not machine-checked; the exact asymptotic rate of
$\log\operatorname{den}(a_n)/n$ (measured $1.3946$ at $n=400$, conjecturally $\log4$) is not
identified; and the conjectural exactness of $v_2(\operatorname{den}b_n^{\rm saddle})=26n+O(\log n)$
is measured only to $n=6$.
