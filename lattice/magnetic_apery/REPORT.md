# Magnetic sources on Apéry's host: the free integration cannot be transplanted

*Working directory `lattice/magnetic_apery/`.  All scripts and raw outputs listed in §10.
PARI/GP 2.15.4, exact over **Q** unless a precision is quoted.
Conventions: `lattice/cooper_sources/REPORT.md` §§2, 5 (magnetic sources, $m\mid c(m)$),
`lattice/cuspform_sources/REPORT.md` §§1, 9 (Fricke hosts, Theorem A, the twelve rows),
`lattice/hostscan/REPORT.md` §4.4 (the twelve CDT-shape Fricke rows),
`consolidation/CDT_FINDER.md` and `lattice/cdt_finder/cdt_bound.py` (entry / margin).
Tags: **[exact]** = closed form or exact rational arithmetic; **[verified, d digits]** =
high-precision computation performed here; **[estimated]**.
Nothing outside this directory was modified.*

---

## 0. Verdict

| claim | verdict |
|---|---|
| **The parametrisation is complete.** On Apéry's ζ(3) host ($N=6$, $C=72$, $B=17$) the eta quotient $u=\eta_2\eta_6^5/(\eta_1^5\eta_3)$ is a **degree-one Hauptmodul** of $X_0(6)$ and $F=D\log u\in M_2(\Gamma_0(6))$ has divisor $(\text{cusp }\tfrac12)+(\text{cusp }\tfrac13)$ and **no zero in $\mathbf H$**. Hence $\Phi=F^2\rho(u)$, $\rho\in\mathbf Q(u)$, $\rho(0)=0$, is *every* weight-4 meromorphic form on $\Gamma_0(6)$ with $a_0=0$ at $\infty$. The search below is therefore a search over the whole space, cut only by the height/degree of the polar divisor. | **[exact]** §1 |
| **The dichotomy on the twelve Fricke rows.** The canonical source $\Phi_0=F\,Dx$ is **magnetic $\iff$ it is meromorphic on $\mathbf H$ $\iff$ the point $x=\infty$ (i.e. $W:=Cu+1/u=-B$) is *not* the cusp pair.** Measured by $R=\limsup|c(m)|^{1/m}$: $R>1.87$ in the eight rows with free integration, $R=1$ (polynomial growth) in the four without — Apéry, AZ $\varepsilon$, AZ $\zeta$, Domb. | **[verified, $m\le400$]** §2 |
| **A complete census of magnetic sources on $\Gamma_0(6)$** (exact lattice computation of $\{\rho:\ D^{-1}\Phi\in\mathbf Z[[q]]\}$, dimensions up to **60**, polar divisors running over every CM orbit of $X_0(6)$ of discriminant $\ge-96$, all pairs and triples, all single poles $W+b$, $\lvert b\rvert\le800$, orders 1–4). Exactly **three** families survive: the two folds, the cusp pair, and the disc$(-8)$/disc$(-12)$ points. | **[exact]** §§3–5 |
| **The unit obstruction (the reason the idea fails).** $1/x = W+B$, so the value of the host coordinate at a pole $\tau_0$ is $x(\tau_0)=1/(W_0+B)$. $W$ is a Hauptmodul of $\Gamma_0(6)+6$ with $q$-expansion $q^{-1}+\mathbf Z[[q]]$, so every singular value $W_0$ is an **algebraic integer**; the polar divisor is Galois-stable, so $\prod_i\lvert W_i+17\rvert=\lvert N(W_0+17)\rvert\ge1$ and therefore $\max_i\lvert W_i+17\rvert\ge1$. The pole contributes a term $\max_i\lvert W_i+17\rvert^{\,n}n^{-2}$ to $A_n\xi-B_n$: **no geometric decay, whatever $k$ is.** In fact $N(W_0+17)=\pm1$ on *every* CM orbit tested. | **[exact]** + **[verified, 14 orbits]** §4 |
| **The only escape is $W=-17$, i.e. $x=\infty$, which on Apéry's host is the cusp pair $\{1/2,1/3\}$** — and every $W_6$-antiinvariant magnetic source with poles only there is **Bol-trivial**: $B_n/A_n$ is *exactly* constant for large $n$ and the limit is rational ($-1,\ \tfrac14,\ 0,\ \tfrac1{36},\ -\tfrac5{36},\dots$). A $(d-2)$-dimensional family exists for each pole order $d\ge3$; all of it is trivial. | **[verified, $n\le380$, exact]** §5 |
| **The one genuine new row.** $\Phi_{-8}=F^2\cdot\dfrac{-u(1+18u+72u^2)(72u^2-1)}{(1+16u+72u^2)^2(1+17u+72u^2)}$ is magnetic ($m\mid c(m)$ for all $m\le1000$), $W_6$-**antiinvariant**, and has a single double pole at the disc$(-8)$ CM point $\tau_0=-\tfrac13+\tfrac{i\sqrt2}{6}$ (height $0.2357<1/\sqrt6$). It equals the canonical source of the *other* level-6 host $(C,B)=(81,14)$. On **Apéry's** host it gives $k=2$ (free integration survives the transplant!) and $B_n/A_n\to\zeta(2)/8=\pi^2/48$. | **[exact]** + **[verified, 185 digits]** §6 |
| …**but it is killed by the unit obstruction**: $x(\tau_0)=1$ exactly, so $A_n\zeta(2)/8-B_n\asymp n^{-2}$ and $d_n^2\lvert A_n\xi-B_n\rvert\to\infty$. The Apéry-perfect rate $\lambda_2=(\sqrt2-1)^4$ never appears in the error term. | **[verified, $n\le250$]** §6 |
| **What was at stake.** A magnetic source on Apéry's host with a *geometric* error would score, in the CDT model at $\lambda_2^{\rm norm}=\lvert c\rvert^{1/2}=1$: $\tau=4.2355$, entry $=+0.8464$, margin $=+0.0053$ at $p_0=7$ — **exactly CDT's own passing score** — against entry $=-0.8984$, margin $=-24.42$ for the real $k=3$ Apéry row. Classically it would have given $\mu(\zeta(2))\le4.6221$, better than the record $5.0954$. | **[exact]** §8 |
| **Task 4** (level-5 and the other two level-6 hosts, plus $s_7$, $s_{10}$, Domb, $N=12$, $s_{18}$): in the family $Q=(W+b)^e$, $\lvert b\rvert\le300$, $e\le3$, the *only* magnetic sources sit at $b=B$ — the host's own configuration (plus its cube, and one extra at $(N,C,B)=(6,81,14)$, $b=18$). $N(W+B)=0$ exactly at the CM discriminant carrying the host's own pole ($-4$ for $N=5$; $-8$ for $C=81$; $-12$ for $C=64$), and $\ge1$ elsewhere. | **[verified]** §7 |
| **Task 5 (Paşol–Zudilin control).** $\Delta/E_4^2$ and $E_4\Delta/E_6^2$ are *strongly* magnetic ($\delta^{-1}F\in q\mathbf Z[[q]]$, no scaling needed, checked to $m=400$). $L(\Xi,2)=4\pi^3\Lambda(\Phi,3)$ exists for $\Delta/E_4^2$ and equals $0.06103925100753799194274901793938221534102612004260042345485\ldots$ — **not** a rational multiple of $1,\zeta(2),\zeta(3),\pi^2,\pi^3$ and not algebraic of degree $\le4$: **an apparently new constant**. For $E_4\Delta/E_6^2$ the double pole sits at $\tau=i$, i.e. *on* the Fricke geodesic, and $\Lambda$ **diverges**. | **[verified, 60+ digits]** §9 |

**One sentence.**  *The idea is exactly right in its arithmetic and exactly wrong in its
geometry: a magnetic weight-4 source really does exist on $\Gamma_0(6)$, really is
$W_6$-antiinvariant, really keeps the free integration $k=2$ when transplanted onto Apéry's
depth-3 host, and really does produce a linear form in $1$ and $\zeta(2)/8$ there — but the
double pole of a magnetic form sits at a CM point, the value of the host coordinate at a CM
point is an algebraic **unit** ($x(\tau_0)=1$ in this case), and a pole at $\lvert x\rvert=1$
contributes an $n^{-2}$ tail that swamps the Apéry-perfect rate $(\sqrt2-1)^4$; the only pole
locus on Apéry's host that is invisible to the coordinate is $x=\infty$, which there is the
**cusp** pair, where the magnetic forms are Bol images $D^3g$ of weight $-2$ modular functions
and their Apéry limits are rational.*

---

## 1. The complete parametrisation on Apéry's host  **[exact]** (`33_verify.gp`)

Host: $N=6$, $u=q\prod_d(q^d;q^d)^{r_d}$ with $r=(-5,1,-1,5)$ on $d=(1,2,3,6)$, $C=72$,
$B=17$, $x=u/(1+17u+72u^2)$, $F=D\log u=1+5q+13q^2+23q^3+\cdots\in M_2(\Gamma_0(6))$,
$A_n=[x^n]F=1,5,73,1445,33001,\dots$ (Apéry's $\zeta(3)$ numbers),
$\lambda_{1,2}=17\pm12\sqrt2$, $\lambda_1\lambda_2=1$, $\lambda_2=(\sqrt2-1)^4=0.0294372515\ldots$

*Facts used throughout.*
$\Gamma_0(6)$ has genus $0$, **no elliptic points**, four cusps
$\infty,0,\tfrac12,\tfrac13$ with $u=0,\infty,-\tfrac19,-\tfrac18$.
$\deg u=1$: $u$ is a **Hauptmodul** of $X_0(6)$.  $\operatorname{div}(F)=(\tfrac12)+(\tfrac13)$
(total degree $2=\tfrac{2}{12}[\mathrm{SL}_2(\mathbf Z):\Gamma_0(6)]$), so $F$ has **no zero in
$\mathbf H$**.  Consequently, for any weight-4 meromorphic $\Phi$ on $\Gamma_0(6)$ the quotient
$\Phi/F^2$ is a modular function, i.e. a rational function of $u$; and if $\Phi$ has rational
$q$-coefficients then $\rho:=\Phi/F^2\in\mathbf Q(u)$.

> **Lemma 1.** $\{\text{weight-4 meromorphic }\Phi\text{ on }\Gamma_0(6),\ a_0(\Phi,\infty)=0\}
> =\{F^2\rho(u):\rho\in\mathbf Q(u),\ \rho(0)=0\}$, and
> $\Phi$ is holomorphic at the cusps $\tfrac12,\tfrac13$ iff $\operatorname{ord}\rho\ge-2$ there,
> holomorphic at $0$ iff $\rho(\infty)=0$.
> $\Phi|_4W_6=\varepsilon\Phi\iff\rho(1/(72u))=\varepsilon\,\rho(u)$.

The last equivalence uses $u|W_6=1/(72u)$ and $F|_2W_6=-F$.  Writing
$w=u+1/(72u)$ (a Hauptmodul of $X_0(6)+6$) and $v=u-1/(72u)$, the antiinvariant $\rho$ are
exactly $\rho=v\,P(w)/Q(w)$ with $\deg P\le\deg Q-2$.

**The normalised pole coordinate.**  Throughout we use
$$W:=72u+\frac1u=72w,\qquad \boxed{\ \frac1x=W+17\ }$$
($1/x=(1+17u+72u^2)/u$).  $W=q^{-1}+O(1)$ has integer $q$-coefficients and is a Hauptmodul of
$\Gamma_0(6)+6$.  The locus $1+B'u+72u^2=0$ is exactly $W=-B'$; the cusps $\tfrac12,\tfrac13$
are $W=-17$ ($x=\infty$); the two Fricke fixed points (folds, disc $-24$) are $W^2=288$,
$x=1/\lambda_1$ and $1/\lambda_2$.

---

## 2. Magnetic $\iff$ meromorphic, on all twelve Fricke rows  **[verified]** (`11_growth.gp`)

For each of the twelve CDT-shape hosts we computed $\Phi_0=F\,Dx=F^2u(1-Cu^2)/(1+Bu+Cu^2)^2$
to $q^{400}$, its exact magnetic lattice, and $R=\max_{380\le m\le400}|c(m)|^{1/m}$.

| host | $C$ | $B$ | roots of $1+Bu+Cu^2$ | $R$ | $\operatorname{Im}\tau_0$ | magnetic? |
|---|---|---|---|---|---|---|
| $N5C125B22$ | 125 | 22 | complex | 3.5747 | 0.2027 | **yes** |
| $N6C81B14$ | 81 | 14 | complex | 4.4607 | 0.2380 | **yes** |
| **$N6C72B17$ (Apéry)** | 72 | 17 | $-\tfrac19,-\tfrac18$ = **cusps** | **1.051** (→1) | — | **no** |
| $N6C64B20$ | 64 | 20 | $-\tfrac1{16},-\tfrac14$ | 6.2123 | 0.2907 | **yes** |
| $N7C49B13$ ($s_7$) | 49 | 13 | complex | 2.2170 | 0.1267 | **yes** |
| $N8C32B12$ (AZ $\varepsilon$) | 32 | 12 | $-\tfrac18,-\tfrac14$ = cusps | **1.049** | — | **no** |
| $N8C16B24$ (Catalan) | 16 | 24 | real | 9.3232 | 0.3553 | **yes** |
| $N9C27B9$ (AZ $\zeta$) | 27 | 9 | complex = cusps | **1.048** | — | **no** |
| $N10C25B6$ ($s_{10}$) | 25 | 6 | complex | 1.9106 | 0.1030 | **yes** |
| $N12C9B10$ (Domb) | 9 | 10 | $-\tfrac19,-1$ = cusps | **1.048** | — | **no** |
| $N12C1B34$ | 1 | 34 | real | 13.133 | 0.4098 | **yes** |
| $N18C1B14$ ($s_{18}$) | 1 | 14 | real | 2.8956 | 0.1692 | **yes** |

The four "$R=1$" values are the finite-$m$ artefact of polynomial growth
($c(m)\sim m^3$ gives $m^{3/m}=1.046$ at $m=400$).

> **Observation 2.**  Across the whole family, *free integration $\iff$ the canonical source is
> meromorphic on $\mathbf H$ $\iff$ $x=\infty$ is not the cusp pair.*  This is Paşol–Zudilin's
> folklore conjecture ("no holomorphic form is magnetic") **and its converse**, verified 12/12.
> Apéry's host fails because $\Gamma_0(6)$ has no elliptic points, so the only $W_6$-stable
> pair of special points available to $B\in\mathbf Z$ is the cusp pair, and $F^2$ kills the pole.

---

## 3. The census of magnetic sources on $\Gamma_0(6)$  **[exact]**

**Method** (`lib.gp`).  For a prescribed denominator $Q$ the admissible $\rho$ form a
$\mathbf Q$-vector space $V$ with an explicit basis; $\Xi_j=D^{-1}(F^2\rho_j)$ are computed
exactly to $q^{M}$; then
$$L(M)=\{\alpha\in\mathbf Q^{\dim V}:\textstyle\sum_j\alpha_j\Xi_j\in\mathbf Z[[q]]\bmod q^{M+1}\}$$
is a decreasing family of full-rank lattices, computed exactly by an HNF/`matsolvemod`
descent over the $M$ coefficient conditions.  **A magnetic element is precisely a vector whose
norm stays bounded as $M\to\infty$**, so the test is: LLL-reduce $L(60)$ and $L(150)$, require
the minimal norm to agree, then verify integrality to $m=400$ (and $1000$ for the survivors).
*Calibration* (`t4`, `03_calib.gp`): on Cooper's $s_7$ host the machine returns
$\alpha=(1,0,-49)$, i.e. exactly $\rho=u(1-49u^2)/(1+13u+49u^2)^2$, with minimal norm $2402$
stable at $M=100,\dots,600$; and on all twelve hosts the own-$Q$ lattice generator is
$-C$ (magnetic) or blows up (not), reproducing `hostscan` §4.4's free-integration column.

**What was searched on Apéry's host** (`20_apery_W.gp`, `21_apery_CM.gp`, `27`, `28`):

| family | space | range | hits |
|---|---|---|---|
| $Q=(W+b)^2$, $\rho(\infty)=0$ | dim 3 | $\lvert b\rvert\le800$ | none |
| $Q=(W+b)^2$, widest | dim 4 | $\lvert b\rvert\le800$ | none |
| $Q=(W+b)^{1,3}$, widest | dim 2, 6 | $\lvert b\rvert\le300$ | $b=17$, order 3 (the cusp tower) |
| $Q=(aW+b)^2$, widest | dim 4 | $2\le a\le40$, $\lvert b\rvert\le40a$ | none |
| $Q=h(W)^{1,2,3}$, $h$ = **CM minimal polynomial**, both spaces | dim 1–36 | all 13 orbits of disc $\ge-96$ + the cusp | fold$^2$, cusp$^3$ |
| $Q=(h_ih_j)^2$, all pairs | dim 3–20 | $\deg\le4$ | $(-8)\cdot$cusp, $(-12)\cdot$cusp |
| $Q=(h_ih_jh_k)^2$, triples | dim 5–20 | $\deg\le5$ | none new |
| $Q=\bigl(\prod_{\deg\le2}h_i\bigr)^2$ (all nine orbits at once) | **dim 59 / 60** | — | none new |
| $Q=(72w+17)^d$, antiinvariant | dim $d-1$, $d\le12$ | — | a $(d-2)$-dim magnetic space, all Bol-trivial (§5) |

("$\rho(\infty)=0$" is $\{uP(u)/Q(u):\deg P\le\deg Q-2\}$, i.e. $a_0(\Phi)=0$ at both cusps
$\infty$ and $0$; "widest" is $\deg P\le\deg Q-1$, $a_0=0$ at $\infty$ only.  **Both contain the
$W_6$-antiinvariant subspace** $\{v\,P(w)/Q(w)\}$, so no antiinvariant magnetic form can escape
them; the last row uses the antiinvariant basis directly.)

CM data: `02_cm.gp`, `02b_cmpoly.gp`, `14_Wpoly.gp` enumerate the $\Gamma_0(6)$-orbits of CM
points for $D=-8,\dots,-163$ by primitive forms $[6k,b,c]$, dedupe by $u(\tau_0)$ at 150
digits, and recognise the minimal polynomial of $W$ exactly.  **All are monic over $\mathbf Z$.**
Only the discriminants with $\max\operatorname{Im}\tau_0<1/\sqrt6$ can ever be used (else the
Fricke geodesic and $\Theta(q_c)$ diverge): $-8,-12,-15,-20,-23,-32,-36,-44,-48,-72,-80$;
$-24$ *is* the fold and $-39,-56,-96,-120$ sit above it.

**Result.** Exactly three families of magnetic $\Phi$ on $\Gamma_0(6)$ were found:

| # | $\rho$ | polar divisor of $\Phi$ | $W_6$ | usable? |
|---|---|---|---|---|
| F | $\dfrac{-2u^2(1+16u+72u^2)(1+18u+72u^2)}{(72u^2-1)^4}$ | order 4 at the **two folds** (disc $-24$) | $+1$ | no: the pole *is* the fold, $B_n/A_n$ diverges |
| C | $v\,P(w)/(72w+17)^d$, $d\ge3$ | poles $\le d-2$ at the **cusps** $\tfrac12,\tfrac13$ | $-1$ | no: **Bol-trivial**, $\xi\in\mathbf Q$ (§5) |
| M | $\dfrac{-u(1+18u+72u^2)(72u^2-1)}{(1+16u+72u^2)^2(1+17u+72u^2)}$ | double pole at the **disc $-8$** CM pair | $-1$ | $k=2$, $\xi=\zeta(2)/8$, but $n^{-2}$ tail (§6) |

(plus $W_6$-*mixed* magnetic forms with cusp poles, dim up to 7 at $d=6$, all with $O(1/n)$
convergence, and one with a double pole at a single disc $-12$ point.)

---

## 4. The unit obstruction  **[exact]** (`25_obstruction.gp`)

**The identity.**  On the polar locus $1+B'u+Cu^2=0$ one has, modulo that polynomial,
$1+Bu+Cu^2\equiv(B-B')u$, hence
$$\boxed{\ x(\tau_0)=\frac{u_0}{1+Bu_0+Cu_0^2}=\frac1{B-B'}=\frac1{W_0+B}\ }$$
— the host coordinate at a pole is the reciprocal of the *shifted singular modulus*.

**The tail.**  $\Theta=D^{-3}\Phi$ is a triple $q$-antiderivative of a double pole at $q_0$, so
$\Theta$ has a dilogarithmic singularity at $x=x(\tau_0)$ and
$$B_n-\xi A_n=[x^n]\bigl(F(\Theta-\xi)\bigr)=\underbrace{O(\lambda_2^n n^{-3/2})}_{\text{fold, killed by }W_6=-1}
+\sum_{i}C_i\,x(\tau_i)^{-n}n^{-2}.$$

**The arithmetic.**  $W$ is a Hauptmodul of $\Gamma_0(6)+6$ with $q$-expansion $q^{-1}+\mathbf Z[[q]]$,
so its singular values are algebraic integers (verified: monic minimal polynomials for all 14
orbits computed).  The polar divisor of a $\mathbf Q$-rational $\rho$ is Galois-stable, so
$$\prod_i\bigl|W_i+17\bigr|=\bigl|N_{K/\mathbf Q}(W_0+17)\bigr|\in\mathbf Z_{\ge1}
\quad\Longrightarrow\quad\max_i\frac1{|x(\tau_i)|}=\max_i|W_i+17|\ \ge\ 1 .$$

> **Theorem U.**  Let $\Phi$ be a weight-4 meromorphic form on $\Gamma_0(6)$ whose poles in
> $\mathbf H$ lie at CM points (equivalently, at points with $W_0$ an algebraic integer), not all
> of them at $W_0=-17$.  Then the Apéry companion on the host $(6,72,17)$ satisfies
> $|A_n\xi-B_n|\gg n^{-2}$, so $d_n^k|A_n\xi-B_n|\to\infty$ for every $k\ge0$ and **no
> irrationality can be extracted, regardless of the free integration.**
> Irrationality would need $\max_i|W_i+17|<e^{-2}=0.13534$.

**[verified]** the norms, for every CM orbit of $X_0(6)$ ($Y:=W+17$, $1/|x_i|=|Y_i|$):

| $D$ | minpoly of $Y$ | $N(Y)$ | $\max_i\lvert Y_i\rvert$ |
|---|---|---|---|
| $-8$ | $Y-1$ | $1$ | **1** |
| $-12$ | $Y+1$ | $-1$ | **1** |
| $-15$ | $Y^2-7Y+1$ | $1$ | 6.854 |
| $-20$ | $Y^2-18Y+1$ | $1$ | 17.944 |
| $-23$ | $Y^6-23Y^5-151Y^4-879Y^3-151Y^2-23Y+1$ | $1$ | 29.207 |
| $-24$ (folds) | $Y^2-34Y+1$ | $1$ | $\lambda_1=33.971$ |
| $-32$ | $Y^4+198Y^2+1$ | $1$ | 14.071 |
| $-36$ | $Y^2+14Y+1$ | $1$ | 13.928 |
| $-44$ | $Y^6-58Y^5+1487Y^4+276Y^3+1487Y^2-58Y+1$ | $1$ | 38.717 |
| $-48$ | $Y^2-52Y+1$ | $1$ | 51.981 |
| $-72$ | $Y^2-98Y+1$ | $1$ | 97.990 |
| $-80$ | $Y^8+80Y^7+10076Y^6-48080Y^5+346246Y^4-48080Y^3+10076Y^2+80Y+1$ | $1$ | 102.254 |
| $-96$ | $Y^4-176Y^3-1002Y^2-176Y+1$ | $1$ | 181.525 |
| $-120$ | $Y^2-322Y+1$ | $1$ | 321.997 |

Every norm is $\pm1$: **the singular values of the Apéry coordinate $x$ are algebraic units.**
(For Theorem U only $|N|\ge1$ is needed, which is the integrality of $W_0$; the unit property
is an extra empirical fact, a Gross–Zagier-type statement for the level-6 Hauptmodul.)
The two "best" orbits, disc $-8$ and $-12$, achieve $\max|Y_i|=1$ exactly — and give the
$n^{-2}$ tail observed in §6, the least bad possible.

---

## 5. The only pole locus the coordinate cannot see: the cusps  **[verified, exact]** (`25`, `27`, `28`)

$x=\infty$ exactly at $W=-17$, i.e. at the cusps $\tfrac12,\tfrac13$ — a pole there contributes
**no** $|x|^{-n}$ tail.  What lives there?

* $Q=(72u^2+17u+1)^d$, $d=1,2$: **no magnetic element** (dim 2, 4).
* $d\ge3$: the $W_6$-antiinvariant space has dimension $d-1$ and its magnetic sublattice has
  rank exactly $d-2$ — verified for $d=3,\dots,12$ (`27_cuspspace.gp`).
* **Every one of them is Bol-trivial.**  For each generator, $B_n/A_n$ is *exactly* constant
  from some $n$ on (agreement between $n=300$ and $n=380$ is exact, not merely 300-digit), and
  $$\xi=-1\ (d=3);\quad \tfrac14,\ 0\ (d=4);\quad 0,\ \tfrac1{36},\ -\tfrac5{36}\ (d=5);\quad 0,\dots\ (d=6).$$
  Equivalently $F\cdot(\Theta-\xi)$ is a *polynomial* in $x$, i.e. $\Theta-\xi$ is a weight-$(-2)$
  meromorphic modular form with integral $q$-expansion and $\Phi=D^3(\Theta)$ is a **Bol image**.
  For $d=3$: $\Theta=1/F-1$ and $B_n=-A_n$ for all $n\ge1$ **[exact]**.

Magnetism of a Bol image is a triviality ($D^{-1}D^3g=D^2g\in\mathbf Z[[q]]$ whenever
$g\in\mathbf Z[[q]]$); these are the "$W_4^0$ with $a>2$" forms of Paşol–Zudilin's Remark 1.
In the *mixed*-parity space one also finds magnetic forms with $k=1$ (e.g.
$\rho=(16416u^7+14148u^6-15914u^5-5945u^4-591u^3-7u^2+u)/(1+17u+72u^2)^4$, $c=1,-60,1107,\dots$,
$m\mid c(m)$ to $m=400$), but they are neither $+1$ nor $-1$ under $W_6$ and their companions
converge like $1/n$ (4.77 digits at $n=380$) — no usable limit.

> **Corollary.**  On Apéry's host, the intersection of "magnetic" with "no polynomial tail" is
> exactly the Bol-trivial locus.  The free integration cannot be transplanted.

---

## 6. The one real row: $\Phi_{-8}$ on Apéry's host  **[exact]** + **[verified, 185 digits]** (`22`, `23`, `24`, `31`, `33`)

$$\Phi_{-8}=F^2\cdot\frac{-u\,(1+18u+72u^2)(72u^2-1)}{(1+16u+72u^2)^2\,(1+17u+72u^2)} .$$

* **[exact]** $\rho(1/(72u))+\rho(u)=0$, so $\Phi_{-8}|_4W_6=-\Phi_{-8}$; $\rho(0)=0$.
* **[exact]** polar divisor: a **double pole at the two disc $-8$ CM points**
  $u_0=-\tfrac19\pm\tfrac{i\sqrt{32}}{144}$, i.e. $\tau_0=-\tfrac13+\tfrac{i\sqrt2}{6}$
  (form $[6,4,1]$) and $W_6\tau_0$; simple zeros at the four cusps (the simple pole of $\rho$
  at $u=-\tfrac19,-\tfrac18$ is over-killed by the double zero of $F^2$); zeros at the two
  folds and at the disc $-12$ points.  $\operatorname{Im}\tau_0=\sqrt2/6=0.2357<1/\sqrt6$.
* **[exact]** $c(m)=1,-16,3,768,-4370,-48,123144,-626688,9,15429920,\dots$, integral and
  **$m\mid c(m)$ for all $m\le1000$**.  $|c(m)|^{1/m}=4.458$ at $m=400$ vs
  $e^{2\pi\sqrt2/6}=4.3972$.
* **[exact]** $\Phi_{-8}$ *is* the canonical source $F'Dx'$ of the level-6 host
  $(C,B)=(81,14)$, $u'=(\eta_3\eta_6/\eta_1\eta_2)^4$ — identical $q$-expansions.  So this row
  is a genuine **transplant**: a magnetic source native to one Fricke host, run on another.
* **[exact]** $U_p$ data (`31`, `32`): $\Xi|U_3=\Xi$ exactly ($c'(3m)=c'(m)$, $m\le333$);
  $c'(p^im)\equiv c'(m)\bmod p^2$ sharply for $7\le p\le43$; $e=3$ at $p=5$; **$e=0$ at $p=2$**
  (no mod-4 congruence with either sign — the only cell of Cooper's table that fails here).
  Strong $p$-magnetic ($p^n\mid m\Rightarrow p^n\mid c(m)$) verified for $p\le13$, $m\le1000$.
* **[verified, $n\le250$]** on Apéry's host the free integration **survives**:
  $B=0,1,15,\tfrac{2674}9,\tfrac{61070}9,\dots$ and $k=\min\{k:d_n^kB_n\in\mathbf Z\}=\boxed{2}$.
* **[verified, 185 digits]** $\xi=\lim B_n/A_n$; `lindep([xi,zeta(2)]) = [-8,1]`, i.e.
  $$\xi=\frac{\zeta(2)}8=\frac{\pi^2}{48}=0.205616758356028304559051895830753148652368737650849804716944\ldots$$
  (This is the same period the $(81,14)$ host carries on its own — consistent with the
  host-independence lemma of `cuspform_sources` §6.)
* **[verified]** but $x(\tau_0)=1$ exactly, and the measured error is
  $$B_n/A_n-\xi\ \sim\ \lambda_2^{\,n}\quad(\text{ratio }0.029377\to\lambda_2=0.0294372515),$$
  i.e. $A_n\xi-B_n=O(n^{-2})$, *not* $O(\lambda_2^n)$: $L_{230},L_{240},L_{250}=
  1.176\cdot10^{-6},\,1.081\cdot10^{-6},\,9.960\cdot10^{-7}$, ratio $\approx(n/n')^{-2}$.
  $d_n^2|L_n|$ grows without bound.  **No irrationality.**

The same computation for the $\operatorname{disc}(-12)$ magnetic source
$\rho=-u/((1+12u)^2(1+9u))$ gives $k=2$ but $W_6$-mixed parity and only 3.96 digits of
agreement at $n=250$ — the $O(1/n)$ regime.  The fold source (family F) has $k=2$ and
$B_n/A_n$ **divergent** (its pole is the evaluation point itself, $R\,q_c=1$).

---

## 7. Task 4 — the other hosts  **[verified]** (`30_hosts.gp`)

For eight hosts we tabulated the singular values of $W=Cu+1/u$ (all minimal polynomials monic
over $\mathbf Z$) and $N(W+B)$, and scanned $Q=(W+b)^e$, $|b|\le300$, $e\le3$, in the widest space.

| host | disc with $N(W+B)=0$ (= the host's own pole) | magnetic $b$ in $Q=(W+b)^e$ |
|---|---|---|
| $N5C125B22$ | $D=-4$ | $b=22$ ($e=2$: $\rho=u-125u^3$; $e=3$) |
| $N6C81B14$ | $D=-8$ | $b=14$ ($e=2,3$), **plus $b=18$, $e=3$** |
| $N6C64B20$ | $D=-12$ | $b=20$ ($e=2,3$) |
| $N7C49B13$ ($s_7$) | $D=-3$ (elliptic points of order 3) | $b=13$ ($e=2,3$) |
| $N10C25B6$ ($s_{10}$) | $D=-4$ (elliptic points of order 2) | $b=6$ ($e=2,3$) |
| $N12C9B10$ (Domb) | **none** ($W=-10$ is the cusp pair) | $b=10$, $e=3$ only (Bol-trivial); plus $b=\pm6$, i.e. $u=\mp1/3$ |
| $N12C1B34$ | $D=-96$ | $b=34$ ($e=2,3$), plus $b=-14$ ($e=2,3$) and $b=-2$ ($e=3$) |
| $N18C1B14$ ($s_{18}$) | $D=-36$ | $b=14$ ($e=2,3$), plus $b=-4$ ($e=2,3$) |
| **$N6C72B17$ (Apéry)** | **none** — $N(W+17)=\pm1$ always | $b=17$, $e=3$ only (Bol-trivial) |

The disc column reproduces `cooper_sources` §2.1 independently ($s_7\mapsto-3$, $s_{10}\mapsto-4$,
$s_{18}\mapsto-36$) and extends it to the five new rows.  The extra hits ($b=18$ at $C=81$;
$b=\pm6$ at Domb; $b=-14,-2$ at $N12C1B34$; $b=-4$ at $s_{18}$) are further magnetic weight-4
forms, apparently not in the literature, but all have $|W+B|\ge4$ and so are useless as sources
by Theorem U.

Away from $b=B$ every $N(W+B)$ is a non-zero integer, so Theorem U applies verbatim on every
host: **a magnetic source can only help the host it is native to.**  Since the level-5 and the
two other level-6 hosts have $\lambda_2^{\rm norm}=4$ they lose $\log4=1.386$ nats to Apéry's
host and gain nothing from a $k$-drop they already have; and Apéry's host, which would gain a
full nat, admits no usable magnetic source at all.  (Level-10/12/18 rows were included in the
scan for completeness; same picture.)

*Caveat.*  For hosts with $\deg u=2$ or $4$ the parametrisation $\Phi=F^2\rho(u)$ is **not**
complete (it only sees forms pulled back from the degree-$\deg u$ cover), so §7's search is a
census of that subspace, not of all of $M_4^{\rm mero}(\Gamma_0(N))$.  On Apéry's host
($\deg u=1$) it *is* complete — that is what makes §§3–5 a genuine census.

---

## 8. What it would have been worth  **[exact]** (`29_score.py`)

CDT model of `lattice/cdt_finder/hosts.py`, number-field normalisation
$\lambda_2^{\rm norm}=|c|^{1/2}=1$ (Apéry-perfect), CDT's own inventory $p_0=7$, $m=14$:

| row | $k$ | $\tau$ | $\log|\varphi'(0)|$ | entry | BC/entry | margin |
|---|---|---|---|---|---|---|
| CDT's own $L(2,\chi_{-3})$ proof | 2 | 4.23546 | 5.08191 | 0.84645 | 13.9938 | **+0.0053** |
| Apéry host, real row | **3** | 5.98036 | 5.08191 | $-0.89845$ | — | $-24.42$ |
| Apéry host, hypothetical magnetic row | **2** | 4.23546 | 5.08191 | **0.84645** | 13.9938 | **+0.0053** |

i.e. a magnetic source on Apéry's host would have reproduced CDT's score *exactly* (their
$\lambda_2=1$ is our $|c|^{1/2}=1$), and would have improved with the inventory
($p_0=10\Rightarrow$ margin $+5.98$; $p_0=14\Rightarrow+15.21$).  Classically:

| $k$ | $|q_n|$ | $|q_n\xi-p_n|$ | $\mu\le$ |
|---|---|---|---|
| 3 (real) | $e^{6.5255n}$ | $e^{-0.5255n}$ | 13.4178 (Apéry) |
| 2 (hypothetical) | $e^{5.5255n}$ | $e^{-1.5255n}$ | **4.6221** (record for $\zeta(2)$: 5.0954) |

The gap between "$k=2$ holds" (it does, §6) and "the error is geometric" (it is not, §4) is
the whole content of this note.

---

## 9. Task 5 — the Paşol–Zudilin level-1 controls  **[verified, 60+ digits]** (`12_pz.gp`)

* $F_{4a}=\Delta/E_4^2$ and $F_{4b}=E_4\Delta/E_6^2$ are **strongly magnetic**:
  $\delta^{-1}F\in q\mathbf Z[[q]]$ with **no scaling**, checked for $m\le400$ (PZ prove
  magnetism; Li–Neururer's $64\Delta/E_4^2$ is not needed).
* Growth: $|c(m)|^{1/m}=231.19$ at $m=400$ vs $e^{\pi\sqrt3}=230.765$ (pole at $\rho$, disc $-3$);
  $535.39$ vs $e^{2\pi}=535.492$ (pole at $i$, disc $-4$).
* $\Lambda(\Phi,s)=\sum_m c(m)\bigl[(2\pi m)^{-s}\Gamma(s,2\pi m)+(2\pi m)^{s-4}\Gamma(4-s,2\pi m)\bigr]$
  converges for $F_{4a}$ (ratio $e^{\pi\sqrt3-2\pi}=0.431$) and reproduces $\Lambda(1)=\Lambda(3)$
  to 78 digits (the level-1 $\varepsilon=+1$ functional equation).  Then
  $$L(\Xi_a,2)=L(F_{4a},3)=4\pi^3\Lambda(3)
  =0.061039251007537991942749017939382215341026120042600423454851567419288\ldots$$
  `lindep` against $\{1,\zeta(2),\zeta(3)\}$, $\{1,\pi^2,\pi^3\}$ and `algdep` to degree 4 all
  return junk (coefficients $\ge10^{19}$).  **It is a new, unidentified constant** — in
  particular the "critical slot is Eisenstein" phenomenon of Cooper's three rows
  ($\zeta(2)/7$, $\zeta(2)/5$, $L(2,\chi_{-3})/2$) does **not** persist to level 1.
* For $F_{4b}$ the double pole sits at $\tau=i$, i.e. **on** the Fricke geodesic $(0,i\infty)$:
  $\Phi_b(iy)\sim\text{const}/(y-1)^2$ and $\Lambda$ diverges (the partial sums grow linearly:
  $0.0371,\,0.0740,\,0.1109,\,0.1477$ at $m\le100,200,300,400$).  There is no critical-slot
  period for $E_4\Delta/E_6^2$ by this route.  This is the level-1 shadow of the condition
  $\operatorname{Im}\tau_0<1/\sqrt N$ that Cooper's three rows satisfy.

---

## 10. Honest ledger

* **Proved here:** Lemma 1 (completeness of the parametrisation on Apéry's host, from
  $\deg u=1$, $\operatorname{div}F$, and no elliptic points); the identity $1/x=W+B$ and hence
  $x(\tau_0)=1/(W_0+B)$; $\max_i|W_i+B|\ge1$ from integrality of singular values and Galois
  stability.  Theorem U's *conclusion* $|A_n\xi-B_n|\gg n^{-2}$ additionally assumes the
  residue $C_i$ at the dominant pole is non-zero — verified in the one case where it matters
  (§6, $C_i\ne0$ measured) but not proved in general.
* **Verified, not proved:** the census of §3 is exhaustive only within the stated
  height/degree boxes.  A magnetic source with a pole at a *non-CM* algebraic point of very
  large height, or with a polar divisor of degree $>36$, is not excluded — although
  Paşol–Zudilin's Shimura–Borcherds heuristic ("magnetic $\Rightarrow$ poles at quadratic
  irrationalities") says CM is the only place to look, and the non-monic scan
  $Q=(aW+b)^2$, $a\le40$, found nothing.
* **Verified, not proved:** $N(W_0+17)=\pm1$ for every CM orbit (14 orbits).  Only
  $|N|\ge1$ is used.
* **Not attempted:** a proof that the cusp-supported magnetic space is *exactly* the Bol
  images (we verified triviality generator by generator up to $d=12$); the identification of
  $L(\Delta/E_4^2,3)$; a Shimura–Borcherds construction of $\Phi_{-8}$.
* **Sharpness of $k$:** $k=2$ for $\Phi_{-8}$ on Apéry's host is minimal over $n\le250$;
  $k=0$ for the Bol-trivial rows; $k=1$ for the mixed cusp$^4$ source.
* **Numerics:** $\xi=\zeta(2)/8$ to 185 digits at $n=250$ ($\lambda_2^n$ convergence);
  all $q$-series arithmetic exact over $\mathbf Q$; CM points at 120–150 digits;
  $L(\Delta/E_4^2,3)$ at 78 digits (truncation error $<10^{-140}$, `realprecision` 80).

**Open.**  (i) Is there a magnetic weight-4 form on $\Gamma_0(6)$ with a pole at a point where
$|W+17|<e^{-2}$?  By Theorem U it cannot be CM, so by the PZ heuristic it should not exist —
a proof would close the question completely.  (ii) The same analysis on the three remaining
$k=3$ hosts (AZ $\varepsilon$, AZ $\zeta$, Domb): their $x=\infty$ locus is also the cusp pair,
so the same dichotomy should hold, but their $\lambda_2^{\rm norm}=4,3\sqrt3,4$ makes them
uninteresting even if it failed.  (iii) Identify $L(\Delta/E_4^2,3)$.

---

## 11. Scripts

| file | contents |
|---|---|
| `lib.gp` | hosts, $u$, $F$, the three parametrisations ($W_6$-antiinvariant / $a_0=0$ at both cusps / widest), the exact magnetic-lattice descent `latmag`, the LLL detector `magmin`/`magsub`/`scanT` |
| `01_scan1.gp/.out` | first pass: the Cooper-shape family $\rho=u(1-Cu^2)/(1+B'u+Cu^2)^2$, $|B'|\le60$, on the level-5 and level-6 hosts |
| `02_cm.gp/.out`, `02b_cmpoly.gp/.out` | CM points of $X_0(6)$: orbits, heights, minimal polynomials of $w=u+1/(72u)$ |
| `03_calib.gp/.out` | calibration: each host's own $Q$; reproduces `hostscan` §4.4's free-integration column |
| `08_diag.gp/.out` | where magnetism dies in the $B'$ family; $(n+1)\nmid A_n$ for Apéry |
| `11_growth.gp/.out` | §2: $R=\limsup|c(m)|^{1/m}$ and the magnetic/meromorphic dichotomy on all twelve hosts |
| `12_pz.gp/.out` | §9: Paşol–Zudilin's $\Delta/E_4^2$, $E_4\Delta/E_6^2$ |
| `14_Wpoly.gp/.out` | the CM singular values in the normalised coordinate $W$ (monic, integral) |
| `20_apery_W.gp/.out` | §3: systematic scan in $W$ ($(W+b)^{1,2,3}$, $(aW+b)^2$) |
| `21_apery_CM.gp/.out` | §3: CM-guided scan (singles, pairs, triples, the full degree-$\le2$ product, dim 59/60) |
| `22_hits.gp/.out` | exact data for the four hits: factorisations, $c(m)$, magnetism to $m=1000$, growth |
| `23_apery_row.gp/.out` | companions $B_n$ on Apéry's host for all four, $k$, limits |
| `24_D8row.gp/.out` | §6: the disc$(-8)$ row in detail — $\xi=\zeta(2)/8$, the $n^{-2}$ tail, $d_n^2|L_n|$ |
| `25_obstruction.gp/.out` | §4: the identity $x=1/(W+B)$, the norm table, the cusp$^d$ magnetic spaces |
| `26_cusp4.gp/.out` | the mixed-parity cusp$^4$ magnetic source ($k=1$, $O(1/n)$) |
| `27_cuspspace.gp/.out` | §5: magnetic sublattices for cusp poles, both parities, $d\le12$ |
| `28_cuspanti.gp/.out` | §5: the $W_6=-1$ cusp tower — all limits rational, convergence exact |
| `29_score.py` | §8: CDT entry/margin for $k=3$ vs $k=2$ at $\lambda_2^{\rm norm}=1$, and $\mu$ |
| `30_hosts.gp/.out` | §7: the other hosts — singular values, $N(W+B)$, magnetic scan |
| `31_upcong.gp/.out`, `32_psi2.gp` | $U_p$ congruences and strong $p$-magnetism for $\Phi_{-8}$ |
| `34_data.gp`, `data_D8_c.txt`, `data_D8_cp.txt`, `data_D8_row.txt` | $c(m)$, $c'(m)=c(m)/m$ for $m\le400$ and the Apéry-host row $(A_n,B_n)$ for $n\le200$ |
| `33_verify.gp/.out` | §6: consolidated verification of $\Phi_{-8}$ (Fricke sign, polar divisor, identity with the $(81,14)$ source) |

Run with `gp -q <file>`; the systematic scans take minutes, everything else seconds.
