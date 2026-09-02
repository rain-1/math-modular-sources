# An explicit modular approximation to $\zeta(7)$: exact sequences, closed-form asymptotics, error law

Working directory: `scratchpad/zeta7/`.  All PARI/GP 2.15.
Every number below is tagged **[exact]** or **[numerical, N digits]**.

---

## 0. The system

**Host.** $h=\eta_1^3\eta_4\eta_6^2/(\eta_2^2\eta_3\eta_{12}^3)=q^{-1}-3+2q+\cdots$, a Hauptmodul of
$\Gamma_0(12)$; $v:=1/h=q+3q^2+\cdots$; the Fricke-invariant coordinate
$$x=\frac{h}{(h+3)(h+4)}=\frac{v}{1+7v+12v^2}=q-4q^2+2q^3+12q^4-17q^5-\cdots \qquad[\text{exact}]$$
is a Hauptmodul of $\Gamma_0(12)^+$.  This is **the same host** as the $\zeta(5)$ showcase.
Fold at $\tau_c=i/\sqrt{12}$ ($q_c=e^{-\pi/\sqrt3}=0.16303\ldots$), $h(\tau_c)=2\sqrt3$,
$x(\tau_c)=x_+=7-4\sqrt3$.  The finite singularities of the rows are $x=-1$ (cusp) and $x=7\mp4\sqrt3$.

**Source.** $\Phi=\frac1{480}\sum_{d\mid12}c_dE_8(d\tau)$, $c=(1,-572,11583,-36608,46332,-20736)$ on
$d=(1,2,3,4,6,12)$; weight $8$ on $\Gamma_0(12)$,
$$\Phi=q-443q^2+13771q^3-93883q^4+78126q^5+571255q^6+823544q^7-12054203q^8+\cdots\qquad[\text{exact}]$$
(the brief in the task said "$1,-571,\dots$"; the true second coefficient is
$g_2=\sigma_7(2)-572=129-572=-443$).

**Rows.**
$$U=\frac{\Phi}{Dx}=\sum_{n\ge0}c_nx^n\ (\text{weight }6),\qquad
\Psi=D^{-7}\Phi=\sum_{m\ge1}g_mm^{-7}q^m,\qquad
V=U\cdot\Psi=\sum_{n\ge0}d_nx^n,\qquad D=q\tfrac{d}{dq}.$$
$$U(q)=1-435q+10285q^2-9041q^3-34947q^4-184658q^5+832527q^6-\cdots\qquad[\text{exact}]$$

**The showcase row** (see §4: $U$ itself has a simple pole at $x=-1$, so the natural object is)
$$U'=(1+x)\,U=\sum_{n\ge0}c'_nx^n,\qquad V'=(1+x)\,V=\sum_{n\ge0}d'_nx^n,\qquad
c'_n=c_n+c_{n-1},\ \ d'_n=d_n+d_{n-1}.$$
$$c'_n:\ 1,\ -434,\ 8110,\ 68734,\ 581830,\ 5357278,\ 53776054,\dots\qquad
d'_1=1,\ d'_2=-\tfrac{55483}{128},\ d'_3=\tfrac{1103599549}{139968}\quad[\text{exact}]$$
$c'_n\in\mathbf Z$ and $D_n^7d'_n\in\mathbf Z$ automatically.  $W'=(1+x)W$ is then genuinely
holomorphic at $x=-1$ with a clean $\log^7$ leading term, and *every* structure below is one step
smaller for the primed row than for the unprimed one.

**Theorem (limit and rate)** — everything below is either exact or verified to the stated precision.
For the showcase row $U'=(1+x)U$, $V'=(1+x)V$:
$$\frac{1728}{209}\,\frac{d'_n}{c'_n}\ \longrightarrow\ \zeta(7),\qquad
\frac{d'_n}{c'_n}-\frac{209}{1728}\zeta(7)\ \sim\ \frac{7\alpha'}{K_7'}\,(-1)^{n+1}\,
\bigl(7-4\sqrt3\bigr)^{n}\,n^{1/2}\,(\log n)^{6},$$
$$c'_n\sim K_7'\,(7+4\sqrt3)^n n^{-3/2},\qquad
K_7'=(8-4\sqrt3)K_7,\qquad
K_7=\frac{3^{27/4}\bigl(739-356\sqrt3\bigr)}{2^{77/6}}\cdot\frac{\Gamma(1/3)^{12}}{\pi^{19/2}}
=72.0651039079814722268\ldots,$$
$$\alpha'=\frac{14161}{5040}=\frac{119^2}{7!},\qquad 7\alpha'=\frac{14161}{720},\qquad
d'_n-\xi c'_n\sim\frac{14161}{720}(-1)^{n+1}\frac{(\log n)^6}{n}.$$
For the unprimed row $U=\Phi/Dx$ (which has a simple pole at $x=-1$, §4) the same constant gives
$d_n-\xi c_n\sim\alpha'(-1)^{n+1}(\log n)^7$ and
$d_n/c_n-\xi\sim(\alpha'/K_7)(-1)^{n+1}x_+^{\,n}n^{3/2}(\log n)^7$.

---

## 1. Exact sequences  (item 1)

Computed by exact rational $q$-series arithmetic (eta products, $E_8$'s, `serreverse`, `subst`),
$q$-precision $930$.

* $c_n\in\mathbf Z$ for **all $n\le900$** **[exact]**.
  $$c_n:\ 1,\ -435,\ 8545,\ 60189,\ 521641,\ 4835637,\ 48940417,\ 521143773,\ 5769030937,\
  65751852165,\ 766674975529,\dots$$
  $\log_{10}|c_{900}|=1026.933\ldots$; $c_n>0$ for **every** $2\le n\le900$ **[exact check]** (only $c_1=-435$ is negative).
* $d_0=0$, $d_1=1$, $d_2=-\tfrac{55611}{128}$, $d_3=\tfrac{2328820355}{279936}$,
  $d_4=\tfrac{243357948415}{3981312}$ **[exact]**.
* $D_n^7d_n\in\mathbf Z$ for **all $n\le900$** ($D_n=\mathrm{lcm}(1,\dots,n)$) **[exact]**, and
  $D_n^6d_n\notin\mathbf Z$ for **every** $2\le n\le900$ (the only $n\le900$ with $D_n^6d_n\in\mathbf Z$
  are $n=0,1$) — so the exponent $7$ is sharp.  For scale: $\mathrm{den}(d_{900})$ has $2726$ digits,
  $D_{900}^7$ has $2728$ digits, and $D_{900}^7/\mathrm{den}(d_{900})=69$.

**Extension to $n=3000$.**  Using the exact order-$19$ recurrence and its exact inhomogeneity
(§6.3), $c_n$ and $d_n$ were continued to $n=3000$ **exactly**; $c_{3000}\in\mathbf Z$
($\log_{10}|c_{3000}|=3428.328$, positive) and $D_{3000}^7d_{3000}\in\mathbf Z$ — a strong
independent check on both the recurrence and the residual polynomial.

Files: `cn.txt` ($c_0\ldots c_{900}$), `cn400.txt`, `dn.txt` ($d_0\ldots d_{400}$),
`dn900.txt` ($d_0\ldots d_{900}$), `cn3000.txt`, `dn3000.txt` ($n\le3000$) — all GP vectors of
exact integers/rationals.

---

## 2. Digit table  (item 2)

$\xi:=\tfrac{209}{1728}\zeta(7)$.  Exact $c_n,d_n$ to $n=3000$ (§1), arithmetic at $3600$ digits.
Leading asymptotics (§4):
$$\text{primed: }\ \frac{d'_n}{c'_n}-\xi\sim\frac{7\alpha'}{K_7'}(-1)^{n+1}x_+^{\,n}n^{1/2}(\log n)^6,
\qquad
\text{unprimed: }\ \frac{d_n}{c_n}-\xi\sim\frac{\alpha'}{K_7}(-1)^{n+1}x_+^{\,n}n^{3/2}(\log n)^7,$$
with $\alpha'=\tfrac{14161}{5040}$, $7\alpha'=\tfrac{14161}{720}$, $x_+=7-4\sqrt3$,
$K_7'=(8-4\sqrt3)K_7$.

**Showcase row $U'=(1+x)U$, $V'=(1+x)V$:**

| $n$ | digits of $\zeta(7)$ in $\frac{1728}{209}\frac{d'_n}{c'_n}$ | $\lvert d'_n/c'_n-\xi\rvert^{1/n}$ | err / leading |
|---:|---:|---:|---:|
| 5 | 2.25 | 0.23222 | 35.77 |
| 10 | 7.10 | 0.15777 | 21.87 |
| 20 | 17.87 | 0.11493 | 14.84 |
| 50 | 51.45 | 0.089654 | 10.31 |
| 100 | 108.17 | 0.081125 | 8.313 |
| 200 | 222.12 | 0.076704 | 6.940 |
| 400 | 450.49 | 0.074383 | 5.954 |
| 900 | 1022.00 | 0.073017 | 5.116 |
| 1500 | 1708.07 | 0.072556 | 4.708 |
| 3000 | 3423.57 | 0.072194 | 4.259 |

**Unprimed row $U,V$** (same $n$; slightly fewer digits, and the ratio converges to $1$ even more
slowly because the log power is $7$ rather than $6$):

| $n$ | digits of $\zeta(7)$ | $\lvert d_n/c_n-\xi\rvert^{1/n}$ | err / leading |
|---:|---:|---:|---:|
| 5 | 1.98 | 0.26323 | 54.33 |
| 20 | 16.76 | 0.13065 | 20.99 |
| 100 | 106.21 | 0.084868 | 10.72 |
| 400 | 447.83 | 0.075530 | 7.409 |
| 900 | 1018.94 | 0.073591 | 6.275 |
| 1500 | 1704.76 | — | 5.729 |
| 3000 | 3419.93 | — | 5.133 |

$\lvert d_n/c_n-\xi\rvert^{1/n}\to x_+=7-4\sqrt3=0.0717967697\ldots$ in both cases; the ratio to the
leading term tends to $1$ only logarithmically (the $(\log n)^{5},\dots$ corrections).

Exact display values **[exact]**:
$$c_n:\ 1,\ -435,\ 8545,\ 60189,\ 521641,\dots\qquad
\frac{1728}{209}\frac{d_n}{c_n}:\ -\frac{576}{30305},\ -\frac{1501497}{3571810},\
\frac{2328820355}{2037879162},\ \frac{243357948415}{251188920576},\dots$$

## 3. The asymptotic constant $K_7$  (item 4)

$$c_n\ \sim\ +K_7\,(7+4\sqrt3)^n\,n^{-3/2}\qquad (\textbf{sign }+;\ c_n>0\text{ for large }n).$$

### 3(0) An exact closed form for the row  [exact — new]

$U$ has weight $6$ on $\Gamma_0(12)$ and $Dv$ weight $2$, so $v^3U/(Dv)^3$ is a modular *function*,
hence rational in the Hauptmodul $v=1/h$.  Exact over-determined fit ($169$ equations, $26$
unknowns; kernel exactly $1$-dimensional, and degree $11$ admits none, so degree $12$ is minimal;
the identity then holds identically to $O(q^{701})$):
$$\boxed{\ \frac{v^3\,U}{(Dv)^3}
=\frac{\mathcal N(v)}{(1+2v)^4(1+3v)^2(1+4v)^2(1+6v)^4}\ }\qquad[\text{exact}]$$
$$\mathcal N(v)=1-398v-3900v^2+95880v^3+2131344v^4+18046944v^5+82569024v^6+216563328v^7$$
$$\qquad\ \ +306913536v^8+165680640v^9-80870400v^{10}-99035136v^{11}+2985984v^{12},$$
$\mathcal N$ is irreducible over $\mathbf Q$; the denominator's four roots $v=-\tfrac12,-\tfrac16$ (where $x=-1$) and
$v=-\tfrac13,-\tfrac14$ (where $x=\infty$) are exactly the four non-trivial cusps.  This is the
$\zeta(7)$ analogue of the $\zeta(5)$ row's $v^2G/(Dv)^2=\frac{(1+3v)(1+4v)(1-12v^2)}{(1+6v)^2(1+2v)^2}$
(REPORT §3.1); note that here the numerator does **not** vanish at $v_c=1/(2\sqrt3)$ — the
Fricke sign is $i^6\varepsilon=+1$, so $U(\tau_c)\ne0$ and this is *not* the CM-vanishing case.

Pleasantly, $\mathcal N$ at the four cusp values of $v$ returns the constant terms $a_0$ of $\Phi$ (§4) up
to explicit rational factors **[exact]**:
$$\mathcal N(-\tfrac12)=119=a_0(1/6),\quad \mathcal N(-\tfrac16)=\tfrac{119}{729}=-\tfrac19a_0(1/2),$$
$$\mathcal N(-\tfrac13)=-\tfrac{697}{81}=\tfrac19a_0(1/4),\quad \mathcal N(-\tfrac14)=-\tfrac{6273}{4096}=-\tfrac1{16}a_0(1/3)$$
(the $v$-value and the cusp are $W_{12}$-partners), a completely independent confirmation of the
constant-term computation that produces $\alpha'$.

At the fold, $v_c=1/(2\sqrt3)$ gives, **entirely by exact algebra**,
$$\frac{U(\tau_c)}{(Dv(\tau_c))^3}=24\sqrt3\cdot\Bigl(\tfrac{475659}{2}\sqrt3-\tfrac{1647621}{4}\Bigr)
=17123724-9885726\sqrt3\qquad[\text{exact}],$$
which is **identically equal** (difference $0$ in $\mathbf Q(\sqrt3)$) to the value
$\dfrac{D\Phi(\tau_c)/(Dv)^5}{144\,x''(h_c)}$ obtained from the census `lindep` product $P_5$.
So the census identification of §3.3 is now *proved* rather than merely `lindep`-guessed, and the
only numerically-identified input left in $K_7$ is the Chowla–Selberg value
$(Dv(\tau_c)/\Omega^2)^3$.

### 3(a) Richardson  [numerical]

$K_n=c_nx_+^{\,n}n^{3/2}$, exact $c_n$ to $n=3000$, nodes $n=3000-45(i-1)$, 600-digit arithmetic.
Relative difference from the closed form: $m=10\to4.8\cdot10^{-33}$, $m=20\to2.6\cdot10^{-60}$,
$m=30\to1.5\cdot10^{-84}$, $m=40\to1.9\cdot10^{-106}$, $m=50\to8.1\cdot10^{-126}$.
**Agreement: $\ge125$ digits.**

For the showcase row, $c'_n\sim K_7'(7+4\sqrt3)^nn^{-3/2}$ with
$$K_7'=(1+x_+)K_7=(8-4\sqrt3)\,K_7=77.23914557843432189271394413662936439303769247609466\ldots,$$
verified by the same Richardson on $c'_n$ ($n\le3000$): $m=22\to2.1\cdot10^{-65}$,
$m=34\to2.8\cdot10^{-92}$.  **Agreement: $\ge91$ digits.**

### 3(b) Fold formula  [numerical, 500 digits]

$U$ has weight $6$ and $U|_6W_{12}=-U$; $i^6\varepsilon=+1$, so at the fixed point
$DU(\tau_c)=\frac{6\sqrt{12}}{4\pi}U(\tau_c)$ **[numerical, $712$ digits: relative error $<10^{-712}$]**, and
$$K_7^2=\frac{x_+\,(DU)(\tau_c)^2}{2\pi\,(-D^2x)(\tau_c)},\qquad
K_7=\frac{6\sqrt{12}}{8\pi^{3/2}}\,|U(\tau_c)|\sqrt{\frac{2x_+}{-D^2x(\tau_c)}} .$$
Numerically at $q_c$ (500 digits):
$$\Phi(\tau_c)=0\ (<10^{-499}),\quad Dx(\tau_c)=0\ (<10^{-500}),$$
$$D\Phi(\tau_c)=-12.876374967680815365097\ldots,\quad
D^2x(\tau_c)=-0.099932959786642016957\ldots,$$
$$Dv(\tau_c)=0.482890362550602961319\ldots,\quad
U(\tau_c)=\frac{D\Phi}{D^2x}(\tau_c)=128.850131079595955598\ldots$$
$$\boxed{K_7=72.065103907981472226808245902065735898065741496059462136043084965084\ldots}$$

### 3(c) Closed form  [exact algebra + numerical verification]

Ingredients (all verified against the $q$-series):

| quantity | value | status |
|---|---|---|
| $x_+$ | $7-4\sqrt3=(2-\sqrt3)^2$ | exact |
| $x''(h_c)$, $h_c=2\sqrt3$ | $\dfrac{168-97\sqrt3}{3}=-\dfrac{x_+^2}{\sqrt3}$ | exact |
| $x_+/(-x'')$ | $12+7\sqrt3$ | exact (matches the $\zeta(5)$ showcase) |
| $x_+^{-5/2}$ | $(2+\sqrt3)^5=362+209\sqrt3$ | exact |
| $Dh(\tau_c)$ | $-12\,Dv(\tau_c)$ | numerical, 712 d |
| $D^2x(\tau_c)$ | $x''(h_c)\,Dh(\tau_c)^2$ | numerical, 710 d |
| $\bigl(Dv(\tau_c)/\Omega^2\bigr)^3$ | $\dfrac{3(45+26\sqrt3)}{2^{17}}=\dfrac{3^{3/2}(2+\sqrt3)^3}{2^{17}}$, i.e. $\dfrac{Dv}{\Omega^2}=\dfrac{\sqrt3\,(2+\sqrt3)}{2^{17/3}}$ | numerical, 712 d (census §3.1) |
| $D\Phi(\tau_c)/(Dv(\tau_c))^5$ | $2025(490730\sqrt3-849969)\cdot\frac{1152\sqrt3-2016}{25}=276169531104-159446553408\sqrt3$ | **verified, 704 digits** |

The census product of REPORT §3.3 is therefore **confirmed**, in the clean form
$$\frac{D\Phi(\tau_c)}{(Dv(\tau_c))^5}=276169531104-159446553408\sqrt3
=-2^5\,3^6\sqrt3\,\bigl(6834986-3946181\sqrt3\bigr)\quad[\text{exact identification, }704\text{ d}]$$
(`lindep` on $[\,D\Phi/(Dv)^5,1,\sqrt3\,]$ at $200$ digits; stable, and the resulting exact value
reproduces the numerical one to $704$ digits).  With $\Omega^2=\Gamma(1/3)^6/\pi^4$, assembling
$$K_7=\frac{\sqrt3}{1152\,\pi^{3/2}}\;\Bigl|\frac{D\Phi}{(Dv)^5}\Bigr|\;(Dv)^2\;
\frac{\sqrt{2x_+}}{|x''|^{3/2}},\qquad
\frac{\sqrt{2x_+}}{|x''|^{3/2}}=\sqrt2\,3^{3/4}(2+\sqrt3)^5,$$
and using $(6834986-3946181\sqrt3)(2+\sqrt3)^7=739-356\sqrt3$ **[exact]** gives the remarkably simple

$$\boxed{\;K_7=\frac{3^{27/4}\,\bigl(739-356\sqrt3\bigr)}{2^{77/6}}\cdot
\frac{\Gamma(1/3)^{12}}{\pi^{19/2}}\;}$$

$$\kappa:=\frac{K_7\,\pi^{19/2}}{\Gamma(1/3)^{12}}
=\frac{3^{27/4}(739-356\sqrt3)}{2^{77/6}}
=27.867359229162067584492644331603077620480318677886350178353849679\ldots$$

*Verification.* Closed form vs. the exact-algebra route of §3(0): **identical**.  Vs. fold formula: **agree to $712$ digits** (at $q$-precision $1100$,
$700$-digit arithmetic); vs. Richardson from the exact $c_n$ to $n=3000$: **$\ge125$ digits**.

*Degree.* $\kappa^6=\dfrac{3^{40}\sqrt3\,(739-356\sqrt3)^6}{2^{77}}\in\mathbf Q(\sqrt3)$ **[exact]**,
so $[\mathbf Q(\kappa):\mathbf Q]\mid12$; $N(739-356\sqrt3)=739^2-3\cdot356^2=165913=11\cdot15083$, and
the minimal polynomial of $\kappa$ over $\mathbf Q$ is exactly of degree $12$:
$$2^{154}\,y^{12}+2^{81}\cdot 8167691949042512779895308753260859635\;y^{6}-3^{81}\cdot165913^{6}=0 .$$
Equivalently $\kappa^{12}-T\kappa^6+N=0$ with
$T=-8167691949042512779895308753260859635/2^{73}$ and $N=-3^{81}165913^6/2^{154}$.

---

## 4. The error law  (item 5)

### The singularity at $x=-1$ is a **pole times $\log^7$**, not a pure $\log^7$

$W:=V-\xi U$ is regular at the fold $x_+=7-4\sqrt3$ (Fricke invariance kills the square-root
branch that $U$ and $V$ each have there), so its radius of convergence jumps from $x_+=0.0718$ to
$1$ and the dominant singularity becomes $x=-1$ — this *is* the mechanism of the approximation.
$x=-1$ is the point $v=-\tfrac12$ / $v=-\tfrac16$, i.e. the $W_{12}$-orbit of cusps
$\{1/2,\ 1/6\}$ of $X_0(12)$ (verified numerically: $h\to-6$ at $1/2$, $h\to-2$ at $1/6$;
the other orbit $\{1/3,1/4\}$ has $h\to-4,-3$, i.e. $x\to\infty$).

Because $x$ is a local coordinate there, $Dx$ vanishes simply, while $\Phi$ does **not** vanish at
that cusp; hence $U=\Phi/Dx$ has a **simple pole** at $x=-1$ and
$$W(x)=\frac{1}{1+x}\sum_{j=0}^{7}G_j(x)\log^j(1+x),\qquad G_j \text{ analytic at }x=-1 .$$
(This is the one structural difference from the $\zeta(5)$ showcase, where the row was
$U_5=7E/Q(x)$ with $Q(-1)=32\ne0$ and hence no pole.)  The cure is the showcase row of §0:
$$W'=(1+x)W=V'-\xi U'=\sum_{j=0}^{7}G_j(x)\log^j(1+x)$$
is holomorphic at $x=-1$ with a genuine $\log^7$ leading term, and
$$\alpha':=G_7(-1)$$
is the single constant governing both rows.

### Route B (modular): the cusp data

$x=-1$ is reached as $\tau\to\tfrac12$ along $\operatorname{Re}\tau=\tfrac12$ (there $q=-e^{-2\pi Y}$
is real negative and $x$ runs down the negative real axis from $0$).  Constant terms of $\Phi$ at
all six cusps, from $\mathrm{const}\bigl(E_k(d\tau)|_k\gamma\bigr)=(\gcd(c,d)/d)^k$ at $a/c$
**[exact]**:
$$a_0(\infty)=0,\quad a_0(0)=\tfrac{P(8)}{480}=0,\quad
a_0(1/2)=-\tfrac{119}{81},\quad a_0(1/3)=\tfrac{6273}{256},\quad
a_0(1/4)=-\tfrac{697}{9},\quad a_0(1/6)=119 ,$$
with widths $w(1/2)=3$, $w(1/3)=4$, $w(1/4)=3$, $w(1/6)=1$.  Note $a_0$ vanishes exactly at the two
cusps where $x=0$, because $P(8)=0$.

**Numerical confirmation** (Richardson in $q'$, 11 nodes $Y'=2.0(0.25)4.5$ at
$\tau=\gamma(iY')$, $\gamma=\binom{1\ 0}{2\ 1}$, $9000$ $q$-terms, 220-digit arithmetic):
$a_0=-1.46913580246913580246913584813\ldots$ against $-119/81$ **[numerical, 26 digits]**;
$\lim (U|_6\gamma)(1+x)=-4.40740740740740740740740742158\ldots$ against $wa_0=-119/27$
**[numerical, 25 digits]** — this last limit is exactly the constant term $\varphi_0$ of the
showcase row $U'=(1+x)U$ at that cusp; and at the other cusp orbit
$a_0(1/3)=24.50390625000000002758\ldots$ against $6273/256$ **[numerical, 16 digits]**.

### The general law

Let $s$ be the cusp with $x(s)=x_c$ finite, of width $w$; let $a_0$ be the constant term of
$\Phi|_k\gamma$ there and $\varphi_0$ the constant term of the row $|_{k-2}\gamma$ there.  Then
$$\boxed{\ W\ \sim\ \varphi_0\,a_0\,\frac{w^{\,k-1}}{(k-1)!}\,\log^{k-1}(x-x_c)\ }
\qquad\text{independently of }\xi .$$
*Proof sketch.*  $W=(\text{row}|_{k-2}\gamma)\bigl[(\Psi|_{2-k}\gamma)-\xi(c\tau'+d)^{k-2}\bigr]$;
Bol's identity $D^{k-1}(\Psi|_{2-k}\gamma)=\Phi|_k\gamma$ forces
$[\tau'^{\,k-1}](\Psi|_{2-k}\gamma)=a_0(2\pi i)^{k-1}/(k-1)!$, while the $\xi$-term has degree
$k-2$; and $\tau'=\frac{w}{2\pi i}\log(x-x_c)+O(1)$.  $\square$

* **Our showcase row.** $U'=(1+x)U$ is holomorphic at $x=-1$ with
  $\varphi_0=\lim(U|_6\gamma)(1+x)=w\,a_0=-\tfrac{119}{27}$
  (**verified numerically to 25 digits**, §Route B above).  With $k=8$, $w=3$, $a_0=-119/81$:
  $$\alpha'=\varphi_0a_0\frac{w^7}{7!}=\frac{a_0^2w^8}{7!}=\frac{(119/81)^2\,3^8}{5040}
  =\boxed{\dfrac{14161}{5040}}=\dfrac{2023}{720}=2.80972\overline{2},\qquad 14161=119^2=(7\cdot17)^2.$$
  The cusp $1/6$ ($w=1$, $a_0=119$) gives the same value, as Fricke symmetry demands.
* **Sanity check on $\zeta(5)$.**  There $k=6$, $w=3$, $a_0(\Phi_5,1/2)=\tfrac{13}{27}$,
  $\varphi_0=7e_0/Q(-1)=\tfrac7{32}\cdot\tfrac8{63}=\tfrac1{36}$, giving
  $\alpha_5=\tfrac1{36}\cdot\tfrac{13}{27}\cdot\tfrac{3^5}{120}=\tfrac{13}{480}$ — **exactly the
  published value**.  [exact]

### Route C (fit): numerical confirmation of $\alpha'$

Fitting $e'_n=d'_n-\xi c'_n=[x^n]W'$ against the exact basis
$[x^n]\bigl((1+x)^i\log^j(1+x)\bigr)$, $0\le i\le I$, $1\le j\le7$ (the $\log^j$ coefficients built
exactly from $(1+x)(L^j)'=jL^{j-1}$), using exact $e'_n$ to $n=3000$ at $3600$-digit precision,
nodes geometrically spread in $[n_0,3000]$:

| $I$ | $n_0=400$ | $n_0=800$ | $n_0=1500$ | $n_0=2000$ |
|---|---|---|---|---|
| 0 | 2.78948 | 2.78633 | 2.78715 | — |
| 1 | 2.81044 | 2.80775 | 2.80758 | — |
| 2 | 2.80609 | 2.80846 | 2.80937 | — |
| 3 | 2.81035 | 2.80989 | 2.80975 | — |
| 4 | 2.80962 | 2.80969 | 2.809719993 | 2.809721772 |
| 5 | — | — | 2.809722444 | **2.8097222335** |

against $14161/5040=2.8097222222\ldots$  Best entry ($I=5$, nodes $2000\ldots3000$): relative error
$4.0\cdot10^{-9}$.  **[numerical, $\ge8$ significant digits; decisive]**

### The laws

$$\boxed{\ d'_n-\xi c'_n\ \sim\ 7\alpha'\,(-1)^{n+1}\frac{(\log n)^6}{n}
=\frac{14161}{720}\,(-1)^{n+1}\frac{(\log n)^6}{n}\ }$$
$$\boxed{\ \frac{d'_n}{c'_n}-\frac{209}{1728}\zeta(7)\ \sim\ \frac{7\alpha'}{K_7'}\,(-1)^{n+1}\,
\bigl(7-4\sqrt3\bigr)^n\,n^{1/2}\,(\log n)^6\ },\qquad
\frac{7\alpha'}{K_7'}=0.25463844023990610004\ldots$$
For the unprimed row (pole present),
$$d_n-\xi c_n\sim\alpha'(-1)^{n+1}(\log n)^7,\qquad
\frac{d_n}{c_n}-\xi\sim\frac{\alpha'}{K_7}(-1)^{n+1}x_+^{\,n}n^{3/2}(\log n)^7,\qquad
\frac{\alpha'}{K_7}=0.03898866538525916525\ldots$$
Signs: $\operatorname{sign}(d_n-\xi c_n)=(-1)^{n+1}$ for **every** $2\le n\le900$ **[exact check]**;
$c_n>0$ for every $2\le n\le900$ **[exact check]** (only $c_1=-435$ is negative).

**On Route A.**  The connection problem at $x=-1$ was not attempted: by §6 the minimal operator
annihilating $U$ has order $7$ and $x$-degree $55$ (with a degree-$42$ apparent-singularity divisor),
so a Frobenius connection computation there is far heavier than in the $\zeta(5)$ case, while Route B
delivers $\alpha'$ *exactly*, with the method independently validated against the published
$\alpha_5=13/480$ and confirmed to $8$ digits by Route C.

---

## 5. Small exact facts  (item 6)

**Dirichlet polynomial** $P(s)=\sum_d c_dd^{-s}$ **[exact]**:
$$P(0)=0,\ P(1)=418,\ P(2)=0,\ P(3)=-11,\ P(4)=0,\ P(5)=\tfrac{11}{12},\ P(6)=0,\
P(7)=-\tfrac{209}{864},\ P(8)=0 .$$
$P$ vanishes at $s=0,2,4,6,8$ — precisely the five points where $\zeta(s-7)$ is non-zero
($\zeta(-7),\zeta(-5),\zeta(-3),\zeta(-1)$ and the pole $\zeta(1)$).  Since
$L(\Phi,s)=P(s)\zeta(s)\zeta(s-7)$ and $\zeta(s-7)=0$ at $s=2,4,6$ anyway, **every period except
$\zeta(7)$ is annihilated**, and
$$L(\Phi,7)=P(7)\zeta(7)\zeta(0)=-\tfrac12\cdot\bigl(-\tfrac{209}{864}\bigr)\zeta(7)
=\tfrac{209}{1728}\zeta(7)=\xi .$$

**Fricke** [numerical, 240 digits, test point $\tau=e^{1.9i}/\sqrt{12}$]:
$$\Phi|_8W_{12}=-\Phi\ (\text{ratio }-1+8\cdot10^{-251}i),\qquad
U|_6W_{12}=-U\ (\text{ratio }-1+3\cdot10^{-251}i),\qquad x\circ W_{12}=x\ (<10^{-251}).$$
$$\Phi(\tau_c)=-3.5\cdot10^{-250}\ (=0),\qquad
U(\tau_c)=128.85013107959595559816133852095848074590517024\ldots$$

**Old coordinate (repository `lattice/zeta7_level60/level12_parent.gp`)**:
$x'=\eta_4^2\eta_{12}^2/(\eta_1^2\eta_3^2)$, $t=x'/(16x'^2+2x'+1)$, $F=\Phi/Dt$,
$A_n=[t^n]F=1,-443,13816,-120335,1042756,-6945002,\dots$

| $n$ | $\lvert B_n/A_n-\xi\rvert^{1/n}$ | digits of $\xi$ |
|---:|---:|---:|
| 50 | 0.75628 | 6.07 |
| 100 | 0.68732 | 16.28 |
| 150 | 0.66179 | 26.89 |
| 200 | 0.64825 | 37.65 |
| 230 | 0.64276 | 44.15 |

so the old rate is $\to 0.6=\tfrac{1/10}{1/6}$ (fold $t_+=1/10$, cusp $t=-1/6$), against
$7-4\sqrt3=0.0718$ for the new host: **$\log(1/0.6)=0.511$ vs $\log(1/(7-4\sqrt3))=2.634$**,
a factor $5.16$ improvement in the exponential rate, and the reason for the new coordinate.

---

## 6. Recurrence and differential operator  (item 3)

### 6.1 The full existence grid  [exact statement, verified mod $p=2^{61}-1$]

Call a relation $\sum_{j=0}^{r}P_j(n)c_{n-j}=0$ ($\deg P_j\le D$, valid for all $n>r$) *of shape
$(r,D)$*.  Under the dictionary $\sum_j x^jQ_j(\theta)$, $\theta=x\,d/dx$, $P_j(n)=Q_j(n-j)$, a shape
$(r,D)$ is exactly an ODE of **order $D$** and **$x$-degree $r$**.  Using the exact $c_n$, $n\le900$,
and mod-$p$ kernels ($p=2^{61}-1$), the dimension of the space of shape-$(r,D)$ relations is
$$\boxed{\ \dim(r,D)=\max\bigl(0,\ (r-12)(D-6)-42\bigr)\ }$$
verified at **every** $(r,D)$ tested: the full grid $6\le r\le30$, $3\le D\le16$ with
$(r+1)(D+1)\le430$; $13\le r\le70$ for each of $D=7,8,9$; and the small-$r$ corners
$(16,16\!-\!18)$, $(15,20\!-\!22)$, $(14,27\!-\!29)$, $(13,48\!-\!50)$, $(12,60)$.
Consequences:

* **There is no relation with $D\le6$ and none with $r\le12$** (the law gives $-42$ there for
  every $D$).  The trade-off is the hyperbola $(r-12)(D-6)\ge43$.
* **Fewest terms: $r=13$, i.e. a $14$-term recurrence** — and then $D=49$ is forced, with the
  solution space exactly $1$-dimensional.  (Compare $\zeta(5)$: $15$ terms of degree $5$.)
* **Minimal ODE order $7$** ($=\mathrm{wt}(U)+1$, the $\mathrm{Sym}^6$ prediction for a weight-$6$
  form as a function of a Hauptmodul), first attained at $x$-degree $r=55$, again with a
  $1$-dimensional space; $\dim(r,7)=r-54$.
* Boundary shapes along the hyperbola, with $\dim$ and the number of unknowns $(r+1)(D+1)$:
  $(13,49){:}1,700$; $(14,28){:}2,435$; $(15,21){:}3,352$; $(16,17){:}2,306$; $(17,15){:}3,288$;
  $(18,14){:}6,285$; $(19,13){:}7,280$; $(20,12){:}6,273$; $(21,11){:}3,\mathbf{264}$;
  $(23,10){:}2,\mathbf{264}$; $(27,9){:}3,280$; $(34,8){:}2,315$; $(55,7){:}1,448$.
* **In the range the task prescribed ($6\le r\le32$, $3\le D\le12$) the first shape that admits a
  relation at all is $(r,D)=(20,12)$, with a $6$-dimensional solution space.**

This is a genuine qualitative difference from the $\zeta(5)$ row, whose minimal operator is
order $5=\mathrm{wt}+1$ with $x$-degree only $14$ and no apparent singularities
(recurrence: $15$ terms, $\deg P_j=5$, $P_0(n)=n^5$, $\chi=(\lambda+1)^4(\lambda^2-14\lambda+1)^5$).
Here the dimension law rewrites as $\dim(r,D)=(D-6)(r-54)+42(D-7)$: the first term counts the
naive left multiples $A\cdot L_{\min}$ of the order-$7$, degree-$55$ operator, and the second says
that **$42$ extra relations are gained per unit of order** by clearing content — the signature of an
apparent-singularity divisor of degree $42$ in the leading coefficient of $L_{\min}$, leaving a true
singular part of degree $55-42=13$ **[observation, not proved here]** — consistent with the fact
that $r=13$ is exactly where the term count bottoms out.

### 6.2 The canonical $14$-term recurrence  [mod-$p$, three independent primes]

$\chi(\lambda)=\sum_j p_j\lambda^{r-j}$, $p_j=[n^D]P_j(n)$, has for its roots the reciprocals of the
singularities, so a *desingularised* operator is one with
$\chi=(\lambda+1)^a(\lambda^2-14\lambda+1)^b$ (roots $-1$ and $7\pm4\sqrt3$ only, $a+2b=r$).

**The unique shape-$(13,49)$ relation is desingularised**, with
$$\boxed{\ \chi(\lambda)=(\lambda+1)^{7}\bigl(\lambda^2-14\lambda+1\bigr)^{3}\ }$$
$$(p_0,\dots,p_{13})=(1,-35,318,462,-8229,-36057,-67052,\ -67052,-36057,-8229,462,318,-35,1)$$
(palindromic), and $n^7\mid P_0(n)$ exactly ($P_0=n^7\cdot(\text{degree }42)$), i.e. $x=0$ is a
**maximally unipotent** point of the order-$49$ operator, exactly as $P_0(n)=n^5$ makes it for
$\zeta(5)$.  All $\deg P_j=49$.  Verified independently at $p=2^{61}-1$, $2^{61}-229$ and
$2^{60}-29$: same $\chi$, same $n^7\|P_0$, kernel always $1$-dimensional.

So the three degrees fit together:  the minimal-*order* operator has order $7$ and $x$-degree $55$;
the minimal-*$x$-degree* (= fewest-terms) operator has $x$-degree $13$ and order $49$; and
$55-13=42$ is exactly the apparent-singularity excess read off the dimension law.

**No palindromic symmetry.**  Unlike $\zeta(5)$ (where $P_{14-j}(n)=-P_j(11-n)$), there is
**no** relation $P_{13-j}(n)=\pm P_j(A-n)$ for any $A\in\tfrac12\mathbf Z\cup\tfrac13\mathbf Z$ with
$|A|\le120$ — only the *leading* coefficients $p_j$ are palindromic (forced by
$\chi(\lambda)=\lambda^{13}\chi(1/\lambda)$, since both $(\lambda+1)$ and $\lambda^2-14\lambda+1$ are).

Clean $\chi$ elsewhere on the boundary of the grid:
$(14,28)$: $(a,b)=(8,3)$; $(15,21)$: $(7,4)$ and $(9,3)$; $(18,14)$: $(8,5),(10,4),(12,3)$;
$(19,13)$: $(7,6),(9,5),(11,4),(13,3)$.  **None** at $(16,17)$, $(17,15)$, $(20,12)$, $(21,11)$,
$(22,11)$, $(23,10)$, $(24,10)$, $(27,9)$, $(34,8)$, $(55,7)$.

### 6.3 Exact recurrences, and the inhomogeneity

Three shapes were lifted to **exact integer polynomials** (mod-$p$ kernels at $67$–$84$ primes near
$2^{61}$; kernel put in reduced row-echelon form, whose pivot set is prime-independent; CRT +
`bestappr`; then exact verification over $\mathbf Q$):

| shape | $\chi$ (content removed) | $\deg P_j$ | max coeff digits | verified |
|---|---|---|---|---|
| $(19,13)$ **canonical** | $(\lambda+1)^7(\lambda^2-14\lambda+1)^6$ — **no apparent factor** | $13$ (all $j$) | $642$ | exact, $n\le\mathbf{3000}$ |
| $(23,10)$ | $4(\lambda+1)^8(\lambda^2-14\lambda+1)^3\,V_9(\lambda)$, $V_9$ irreducible | $10$ | $577$ | exact, $n\le900$ |
| $(17,15)$ | $34(\lambda+1)^7(\lambda^2-14\lambda+1)^4\,V_2(\lambda)$, $V_2$ irreducible | $15$ | $715$ | exact, $n\le900$ |

The $(19,13)$ constrained subspace is **exactly $1$-dimensional**: the desingularised recurrence is
unique up to scale.  Uniform structure at every shape:
$$P_0(n)=n^7\,Q_0(n),\qquad P_r(n)=\bigl(n-(r-1)\bigr)^7Q_r(n),\qquad Q_0,Q_r\ \text{irreducible},\
\deg=D-7,$$
so $x=0$ **and** the reflected end are both maximally unipotent.  Contents of the individual $P_j$
are $1,2$ or $4$ except $P_r$, whose content is $6273=3^2\cdot17\cdot41$ at both $(23,10)$ and
$(19,13)$.  Coefficients are $577$–$741$-digit integers, and this is **intrinsic**: saturating each
solution lattice (`matrixqz(·,-1)`) and LLL-reducing gives shortest vectors of $580/738/642$
digits.  **No compact TeX form exists** — the machine-readable files are the usable object.

**Inhomogeneity.**  For *every* $c$-annihilator $L$ of shape $(r,D)$,
$$\boxed{\ L[d]_n=(-1)^n\,W_L(n),\qquad W_L\in\mathbf Z[n],\ \deg W_L\le D-7\ }$$
and $L\mapsto W_L$ has rank exactly $D-6$, which is precisely why $\dim_d=\dim_c-(D-6)$.
For the canonical $(19,13)$: $\deg W=6$, $\mathrm{cont}(W)=2^{17}\cdot7^2\cdot17$, verified
**exactly for all $n\le3000$**.  For $(23,10)$: $\deg W=2$, $\mathrm{cont}(W)=2^{21}\cdot7^3\cdot17^2$,
exact to $n=900$.  At $(19,13)$, $\dim=7=D-6$, so the residual map is bijective and there is a
*unique* element with $R_n=(-1)^nC$, $C$ constant — a negative $659$-digit integer with smooth part
$2^{17}3^55\cdot7\cdot17$ — but its $\chi$ then acquires a degree-$5$ apparent factor.  **One cannot
have a clean $\chi$ and a constant right-hand side simultaneously**, unlike $\zeta(5)$
(where $13\cdot2^{18}(-1)^{n+1}$ accompanies a clean $\chi$).  The primes $7,17$ in the contents are
those of $a_0=\mp7\cdot17/81$ and $\alpha'=(7\cdot17)^2/7!$.

**Palindromic symmetry: none.**  Tested mod $p$ for $A=m/2$, $|A|\le80$, at all three shapes, and
also as the stronger property that $v\mapsto(P_{r-j}(A-n))_j$ preserves the whole solution space;
the scanner was validated on Apéry's $\zeta(3)$ recurrence (it correctly returns $A=1$,
$\varepsilon=-1$).  Structural reason: $P_0=n^7Q_0$ and $P_r=(n-(r-1))^7Q_r$ force $A=r-1$, but
$\mathrm{lc}(Q_r)\ne\pm\mathrm{lc}(Q_0)$.  Likewise none for the $(13,49)$ and primed $(12,49)$
shapes ($A\in\tfrac12\mathbf Z$, $|A|\le120$).

### 6.4 The showcase row $c'_n$, and the four-fold ladder

$c'=(1+S^{-1})c$ and $d'=(1+S^{-1})d$ (with $S^{-1}y_n=y_{n-1}$).  The four dimension laws are the
*same* function shifted:
$$\dim_{c'}(r,D)=\max\bigl(0,(r-11)(D-6)-42\bigr),\quad
\dim_{c}(r,D)=\max\bigl(0,(r-12)(D-6)-42\bigr),$$
$$\dim_{d'}(r,D)=\max\bigl(0,(r-12)(D-6)-42\bigr),\quad
\dim_{d}(r,D)=\max\bigl(0,(r-13)(D-6)-42\bigr),$$
i.e. $\dim_{c'}(r,D)=\dim_c(r+1,D)=\dim_d(r+2,D)$ and $\dim_{d'}=\dim_c$
**[verified mod $p$ at 20+ shapes each on the exact sequences to $n=3000$]**.  Fewest terms:
$$c'_n:\ \mathbf{13}\ \text{terms}\ (r=12,\ D=49),\qquad c_n:\ 14\ (13,49),\qquad
d'_n:\ 14\ (13,49),\qquad d_n:\ 15\ (14,49),$$
each with a $1$-dimensional solution space; and $(11,D)$ is empty for every $D$, so $13$ terms is
absolutely minimal for the showcase row.

**The minimal shape-$(12,49)$ relation for $c'_n$ is desingularised** (three independent primes):
$$\chi'(\lambda)=(\lambda+1)^{6}\bigl(\lambda^2-14\lambda+1\bigr)^{3}=\chi(\lambda)/(\lambda+1),$$
$$(p_0,\dots,p_{12})=(1,-36,354,108,-8337,-27720,-39332,-27720,-8337,108,354,-36,1)\ \text{(palindromic)},$$
with $n^7\,\|\,P'_0(n)$.  Minimal ODE order is again $7$, now first at $x$-degree $54$ (one less
than $55$).

**An exact $19$-term recurrence for $c'_n$, for free.**  The canonical $(19,13)$ operator satisfies
$\sum_{j=0}^{19}(-1)^jP_j(n)\equiv0$ **[exact]**, i.e. $(1+S^{-1})$ right-divides it:
$$L=M\circ(1+S^{-1}),\qquad M_j(n)=\sum_{k=0}^{j}(-1)^{j-k}P_k(n)\quad(j=0,\dots,18).$$
$M$ has shape $(18,13)$, $\chi_M=(\lambda+1)^6(\lambda^2-14\lambda+1)^6$, $n^7\mid M_0$, and
**annihilates $c'_n$ exactly for all $n\le3000$ (0 failures)**, with **exactly the same**
inhomogeneity on the primed companion:
$$\sum_{j=0}^{18}M_j(n)\,d'_{n-j}=(-1)^n\,W(n),\qquad \deg W=6,\ \mathrm{cont}(W)=2^{17}7^217,$$
verified exactly for all $19\le n\le3000$.  Saved as `rec_prime_18_13.gp` (`PVP`).

---

## 7. What this does and does not do

The exponential rate is $\log\frac1{7-4\sqrt3}=2.6339158$ against the denominator growth
$\log D_n^7\sim7n$: $2.63<7$, so no irrationality statement follows (as for $\zeta(5)$, where
$2.63<5$).  What is explicit here is the whole package: modular source with **every period other than
$\zeta(7)$ annihilated by construction** ($P(0)=P(2)=P(4)=P(6)=P(8)=0$), integral $c'_n$, sharp
$D_n^7$ denominators, closed-form $K_7'$ as a $\Gamma(1/3)^{12}/\pi^{19/2}$ period of degree $12$
over $\mathbf Q$, the exact error constant $\alpha'=14161/5040$ with its cusp-theoretic proof, and
an absolutely minimal $13$-term recurrence with a clean characteristic polynomial.

Compared with the repository's earlier level-12 $\zeta(7)$ parent (host
$t=x'/(16x'^2+2x'+1)$), the new Fricke host improves the rate from $0.6$ to
$7-4\sqrt3=0.0718$ — the exponential gain that motivated the change of coordinate.

---

## 8. Files

| file | content |
|---|---|
| `exact_seq.gp` | exact $c_n$, $n\le900$, from the $q$-series ($q$-precision $930$) $\to$ `cn.txt` |
| `exact_dn.gp` | exact $c_n,d_n$, $n\le400$ $\to$ `cn400.txt`, `dn.txt`; integrality tests |
| `dn900.gp` | exact $d_n$, $n\le900$ $\to$ `dn900.txt` |
| `integ.gp`, `integ2.gp` | $c_n\in\mathbf Z$; $D_n^7d_n\in\mathbf Z$ sharp |
| `fold.gp` | CM point $\tau_c$: all fold data, $K_7$ by the fold formula, census checks $\to$ `foldvals.txt` |
| `closed.gp`, `closed2.gp`, `algsteps.gp`, `minpol.gp` | closed form for $K_7$, every algebraic step, minimal polynomial, Richardson |
| `table_facts.gp`, `table2.gp` | digit table, $P(j)$, exact display fractions |
| `fricke_old.gp` | Fricke signs, $\Phi(\tau_c)=0$, $U(\tau_c)$; old-coordinate comparison |
| `rho.gp`, `rhofac.gp`, `rhoexact.gp`, `rhoverify.gp`, `crosscheck.gp` | exact rational form of $v^3U/(Dv)^3$; $U(\tau_c)$ exactly; cross-check with the census $P_5$ |
| `p5check.gp`, `p5check2.gp` | CM identities at $700$-digit precision |
| `cusps.gp`, `cusp_a0.gp`, `cusp_a0b.gp` | which cusp is $x=-1$; constant terms $a_0$; $\lim(U|\gamma)(1+x)$ |
| `errlaw1.gp` … `errlaw4.gp` | singularity type at $x=-1$; Route-C fit for $\alpha'$ (superseded by `alpha3000*.gp`) |
| `signs.gp` | sign of $c_n$ and of $d_n-\xi c_n$ for all $n\le900$ |
| `search_rec.gp`, `search_rec2.gp`, `search_rec3.gp`, `minr.gp`, `charpoly_search.gp`, `chicorners.gp` | recurrence/ODE existence grid; small-$r$ corners; clean-$\chi$ search |
| `chi1349.gp`, `sym1349.gp`, `sym1349b.gp` | the canonical $14$-term shape $(13,49)$: $\chi$, palindromy, $n^7\|P_0$, symmetry scan |
| `rec_19_13.gp`, `rec_23_10.gp`, `rec_17_15.gp` | **exact** recurrences (GP vectors `PV` of integer polynomials in `n`) |
| `verify_rec.gp`, `verify_all.gp`, `p0check.gp`, `maxdig.gp` | independent exact verification of all three recurrences on $n\le900$ and of $R_n=(-1)^nW(n)$; coefficient sizes |
| `extend.gp` | exact continuation of $c_n,d_n$ to $n=3000$ via `rec_19_13.gp` $\to$ `cn3000.txt`, `dn3000.txt` |
| `alpha3000.gp`, `alpha3000b.gp` | Route-C fit for $\alpha'$ on the primed row to $n=3000$ |
| `primed.gp`, `rich3000.gp` | $K_7'$ by Richardson; primed and unprimed digit tables; symmetry scan |
| `grid_cp.gp`, `cp1249.gp`, `dgrid_min.gp` | dimension grids for $c'_n$, $d_n$, $d'_n$; minimal shapes |
| `divide.gp`, `verify_prime.gp`, `rec_prime_18_13.gp` | exact right division $L=M(1+S^{-1})$; the primed recurrence and its verification to $n=3000$ |
| `rec_19_13_constres.gp`, `constres_C.txt`, `resid_W.txt`, `P0P1_*.txt` | constant-residual operator; residual polynomials; $P_0,P_1$ in full |
| `cn.txt`, `cn400.txt`, `dn.txt`, `dn900.txt`, `cn3000.txt`, `dn3000.txt`, `foldvals.txt` | data |
