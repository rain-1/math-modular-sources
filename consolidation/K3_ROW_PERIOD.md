# The period of the elliptic-$K3$ four-term row

*The row of `FOUR_TERM_SCAN.md` §5.3 — the unique row of the four-term census
whose Picard–Fuchs system is an elliptic surface, and the first elliptic $K3$
to appear in the Apéry world.  Scripts and logs: `lattice/k3_period/`.*

---

## 0. Summary

$$\boxed{\ (n+1)^2u_{n+1}=(11n^2+11n+4)u_n-(37n^2+3)u_{n-1}+3(3n-1)(3n-2)u_{n-2}\ }$$
$$u_0=1:\quad 1,4,16,64,250,928,3136,8704,11866,-79400,-975488,\dots$$

| item | value |
|---|---|
| Weierstrass model | $y^2=x^3-27c_4(t)x-54c_6(t)$, $c_4,c_6$ of §1 |
| fibres | $I_4(0),\ I_6(1),\ I_3,I_3\bigl(\tfrac{5\pm i\sqrt2}{27}\bigr),\ IV^*(\infty)$ |
| Picard number | $\rho=20$ (singular $K3$), $\operatorname{MW}$ rank $0$ |
| $\operatorname{MW}_{\rm tors}$ | $\mathbf Z/3$ (a section defined over $\mathbf Q(\sqrt3)$) |
| transcendental lattice | $\det T=72$, $T\cong\begin{pmatrix}6&0\\0&12\end{pmatrix}$, CM by $\mathbf Q(\sqrt{-2})$ |
| newform | **`32.3.d.a`** $=$ `8.3.d.a` $\otimes\chi_{-4}$, weight $3$, level $32$, nebentypus $\chi_{-8}$, CM by $\mathbf Q(\sqrt{-2})$ |
| period | $\displaystyle L(g,2)=\frac{\pi}{32}\,
  \frac{\Gamma(\tfrac18)\Gamma(\tfrac38)}{\Gamma(\tfrac58)\Gamma(\tfrac78)}$ |
| complex fold constant | $\displaystyle \xi\Bigl(\tfrac{5\mp i\sqrt2}{27}\Bigr)=\frac{\sqrt2\mp i}{3}\,L(g,2)$ |
| real ($I_6$) fold constant | $\displaystyle \xi(1)=\frac{2\sqrt2}{3}\,L(g,2)=\xi_++\xi_-$ |
| $p$-adic | $\xi_2=0$ (logarithmic rate); no $3$-adic limit |

Everything below is verified to $\ge165$ decimal digits, or exactly over
$\mathbf Q$.  **This is the first Apéry row in the census whose Apéry constant
is a critical value of a weight-$3$ cusp form rather than of an Eisenstein
series** — the archimedean signature of the jump from rational elliptic surfaces
to $K3$s.

---

## 1. The Weierstrass model

### 1.1 The $\mathcal J$-map, made explicit

`lattice/k3_period/01_jmap.gp` re-runs the certification of
`lattice/four_term/06_jtest.gp` on this row and prints the rational function
$\mathcal J=U/V$ with $\mathcal J(t)=j\bigl(-q(t)^4\bigr)$, $q(t)=t\exp(g/y_0)$
the canonical nome ($h=4$, $\gamma=-1$, $\deg\mathcal J=16$; fitted from $38$
coefficients, checked on $73$ more).  Factoring (`02_fac.gp`):

$$V=t^4(t-1)^6(27t^2-10t+1)^3,\qquad
  U=-(4t-1)^3\bigl(216t^4+64t^3-48t^2+12t-1\bigr)^3,$$
$$U-1728V=-\bigl(5832t^8+34560t^7-30016t^6+12624t^5-4380t^4+1280t^3-240t^2+24t-1\bigr)^2 .$$

The two conjugate poles are the roots of $27t^2-10t+1$, i.e.
$t=\tfrac{5\pm i\sqrt2}{27}$, and $\mathcal J$ has a simple zero at $t=\infty$.
Everything in the fibre table of `FOUR_TERM_SCAN.md` §5.3 follows.

### 1.2 Branch reconstruction

Following the method of `chatgpt-research-archive/The Conductor-3 Apéry
Surface_ Picard Maximality and Explicit Branch Reconstruction.md`: with
$j=c_4^3/\Delta$ and $c_4^3-c_6^2=1728\Delta$ the shape of $U,V$ forces
$$c_4=w^2\,(4t-1)P_4(t),\qquad c_6=\pm w^3Q_8(t),\qquad \Delta=-w^6V(t)$$
for a single constant $w$ (a quadratic twist), where $P_4,Q_8$ are the quartic
and the octic above.  **The twist is fixed by the period.**  With
$y_0(t)=\sum u_nt^n$ and $Q(t)=-q(t)^4$, `03_model.gp` verifies as power series
in $t$ (to $O(t^{41})$, exact over $\mathbf Q$)

$$\boxed{\ c_4(t)=\frac{E_4\bigl(Q(t)\bigr)}{y_0(t)^4},\qquad
          c_6(t)=\frac{E_6\bigl(Q(t)\bigr)}{y_0(t)^6}\ }$$

with **both** ratios identically $1$ — so $w^2=w^3=1$, $w=1$, and the classical
$c_4=(2\pi/\omega_1)^4E_4(\tau)$, $c_6=(2\pi/\omega_1)^6E_6(\tau)$ identify the
period of the model as $\omega_1(t)=2\pi\,y_0(t)$.  Hence:

$$\boxed{\begin{aligned}
c_4(t)&=864t^5+40t^4-256t^3+96t^2-16t+1=(4t-1)\bigl(216t^4+64t^3-48t^2+12t-1\bigr),\\
c_6(t)&=-\bigl(5832t^8+34560t^7-30016t^6+12624t^5-4380t^4+1280t^3-240t^2+24t-1\bigr),\\
\Delta(t)&=\tfrac1{1728}\bigl(c_4^3-c_6^2\bigr)=-\,t^4(t-1)^6(27t^2-10t+1)^3,
\end{aligned}}$$
$$\mathcal X:\qquad y^2=x^3-27c_4(t)\,x-54c_6(t)\qquad\text{over }\mathbf Q(t).$$

(The exact identity $c_4^3-c_6^2=1728\Delta$ is checked over $\mathbf Q[t]$ in
`03_model.gp`.)

### 1.3 The fibres, from the model

$\deg c_4=5\le8$, $\deg c_6=8\le12$, $\deg\Delta=16\le24$: $\mathcal X$ is the
Weierstrass model of an elliptic surface with $\chi=2$.

| $t$ | $v(c_4)$ | $v(c_6)$ | $v(\Delta)$ | Kodaira |
|---|---|---|---|---|
| $0$ | $0$ | $0$ | $4$ | $I_4$ |
| $1$ | $0$ ($c_4(1)=729$) | $0$ | $6$ | $I_6$ |
| $27t^2-10t+1$ | $0$ ($\gcd(c_4,27t^2{-}10t{+}1)=1$) | $0$ | $3$ each | $I_3,I_3$ |
| $\infty$ | $8-5=3$ | $12-8=4$ | $24-16=8$ | $IV^*$ |

$\sum e=4+6+3+3+8=24$: an elliptic $K3$.  This *derives* the fibre table that
`FOUR_TERM_SCAN.md` §5.3 read off the $\mathcal J$-map.

---

## 2. $\rho=20$, the torsion, and $\det T$

**Shioda–Tate.**  $\operatorname{Triv}=U\oplus A_3\oplus A_5\oplus A_2\oplus
A_2\oplus E_6$ has rank $2+3+5+2+2+6=20$, so $\rho\ge20$; for a complex $K3$
$\rho\le h^{1,1}=20$.  Hence

$$\boxed{\rho(\mathcal X)=20\quad\text{(a singular $K3$)},\qquad
\operatorname{rank}\operatorname{MW}=\rho-\operatorname{rank}\operatorname{Triv}=0.}$$

**Torsion.**  $|\!\det\operatorname{Triv}|=4\cdot6\cdot3\cdot3\cdot3=648$ and
$[\operatorname{NS}:\operatorname{Triv}]=|\operatorname{MW}_{\rm tors}|=:m$, so
$\det T=|\!\det\operatorname{NS}|=648/m^2$.

* *A $3$-torsion section exists.*  With $A=-27c_4$, $B=-54c_6$ the $3$-division
  polynomial $\psi_3=3x^4+6Ax^2+12Bx-A^2$ factors over $\mathbf Q(t)$ as
  $$\psi_3=3\bigl(x-9(4t-1)^2\bigr)\cdot(\text{irreducible cubic}),$$
  and at $x_0=9(4t-1)^2$
  $$y_0^2=x_0^3+Ax_0+B=432\,(t-1)^4(27t^2-10t+1)^2 ,$$
  so $y_0=12\sqrt3\,(t-1)^2(27t^2-10t+1)$: **a $3$-torsion section, defined over
  $\mathbf Q(\sqrt3)$, geometrically over $\overline{\mathbf Q}(t)$**
  (`04_torsion.gp`, `04b.gp`).  Hence $3\mid m$.
* *No more.*  $T$ is an even positive-definite rank-$2$ lattice, so
  $\det T\equiv0,3\pmod4$; $648/m^2$ with $m=2,6$ gives $162,18\equiv2\pmod4$ —
  excluded.  $m=9$ needs $\mathbf Z/9$ or $(\mathbf Z/3)^2$ inside
  $\prod_vG_v=\mathbf Z/4\times\mathbf Z/6\times(\mathbf Z/3)^3$ (Shioda: the map
  $\operatorname{MW}_{\rm tors}\to\prod G_v$ is injective); that group has
  exponent $12$, killing $\mathbf Z/9$, and full $3$-torsion over
  $\overline{\mathbf Q}(t)$ would force $3\mid4$ at the $I_4$ fibre (Tate
  parametrisation), killing $(\mathbf Z/3)^2$.  So

$$\boxed{\operatorname{MW}(\overline{\mathbf Q}(t))=\mathbf Z/3,\qquad
\det T=648/9=72 .}$$

**The lattice $T$.**  Even positive-definite binary lattices of determinant $72$
come in three classes: $\operatorname{diag}(2,36)$, $\operatorname{diag}(4,18)$
(the two primitive forms $x^2+18y^2$, $2x^2+9y^2$ of discriminant $-72$) and
$\operatorname{diag}(6,12)$ (the imprimitive $3(x^2+2y^2)$).  The first two have
$3$-part $\mathbf Z/9$ in the discriminant group, the third has
$(\mathbf Z/3)^2$.  Computing the discriminant form of $\operatorname{NS}$: the
$3$-torsion section is an isotropic element $\ell$ of
$q_{\operatorname{Triv}}$ on
$\mathbf Z/4\oplus\mathbf Z/6\oplus\mathbf Z/3\oplus\mathbf Z/3\oplus\mathbf Z/3$
with $q$-values $-3/4,-5/6,-2/3,-2/3,-4/3$ on the generators; $\ell$ is
$0$ in $\mathbf Z/4$ ($\mathbf Z/4$ has no $3$-torsion), and isotropy
$\tfrac23a+\tfrac43c_1+\tfrac43c_2+\tfrac23c_3\equiv0\pmod{2}$ together with
$c_1=c_2$ (the two $I_3$'s are conjugate over $\mathbf Q(\sqrt3,\sqrt{-2})$,
which fixes the section) forces $\ell$ to be **non-trivial at all four of
$I_6,I_3,I_3,IV^*$**.  Then the $3$-part of
$\ell^{\perp}/\langle\ell\rangle$ is a quotient of a subgroup of
$(\mathbf Z/3)^4$, hence elementary abelian of order $9$.  Therefore

$$\boxed{\ T(\mathcal X)\cong\begin{pmatrix}6&0\\0&12\end{pmatrix}
 \ \ \bigl(=3\cdot\langle x^2+2y^2\rangle,\ \operatorname{disc}=-72\bigr),
 \qquad \text{CM field }\mathbf Q(\sqrt{-2}).}$$

The $2$-parts agree as a check: $\mathbf Z/4\oplus\mathbf Z/2$ on both sides.

---

## 3. The motive: the newform `32.3.d.a`

### 3.1 The trace identity used

For $\pi:\mathcal X\to\mathbf P^1$ with section and $\mathcal F=R^1\pi_*\mathbf Q_\ell$
on the good locus $U$, $H^0_c(U,\mathcal F)=H^2_c(U,\mathcal F)=0$ (the monodromy
is irreducible and non-trivial), so
$\sum_{t\in U(\mathbf F_p)}a_t=-\operatorname{tr}(\mathrm{Frob}\mid H^1_c(U,\mathcal F))$;
and in
$$0\to\bigoplus_{v\in S}\mathcal F^{I_v}\to H^1_c(U,\mathcal F)\to
H^1(\mathbf P^1,j_*\mathcal F)\to0$$
the skyscraper term contributes exactly $\sum_{v\ \mathrm{rational}}a_v$ with
$a_v=\pm1$ (split/non-split $I_n$) resp. $0$ (additive).  The two cancel, so

$$\boxed{\ a_p:=\operatorname{tr}\bigl(\mathrm{Frob}_p\mid T\otimes\mathbf Q_\ell\bigr)
=-\sum_{t\in\mathbf P^1(\mathbf F_p)}a_t
=\sum_{t\in\mathbf F_p}\ \sum_{x\in\mathbf F_p}
\chi\bigl(x^3-27c_4(t)x-54c_6(t)\bigr)\ }$$

($\chi$ the quadratic character, $\chi(0)=0$; the fibre at $t=\infty$ is $IV^*$,
$a_\infty=0$; here $H^1(\mathbf P^1,j_*\mathcal F)=T\otimes\mathbf Q_\ell$ because
$\rho=20$).  The bad points $0,1,\tfrac{5\pm i\sqrt2}{27}$ stay distinct and
distinct from each other modulo every $p\ge5$: the relevant resultants are
$\operatorname{Res}(t,t-1)=-1$, $(27t^2-10t+1)|_{t=0}=1$,
$(27t^2-10t+1)|_{t=1}=18$, $\operatorname{disc}(27t^2-10t+1)=-8$, and the
leading coefficient $27$ — all supported on $\{2,3\}$.  So the fibre
configuration is constant for $p\ge5$.
Script `05_ap.gp` / `21_ap_ext.gp`.

### 3.2 The match

`06_newform.gp`, `07_label.gp`, `23_apcheck.gp`: the resulting $a_p$ agree with
the unique newform of $S_3^{\rm new}(\Gamma_0(32),\chi_{-8})$ for **all $76$
primes $5\le p\le397$**, with no free parameter.  LMFDB label

$$\boxed{\ g=\texttt{32.3.d.a}\ =\ \bigl(\eta(z)^2\eta(2z)\eta(4z)\eta(8z)^2\bigr)\otimes\chi_{-4}
 \ =\ \texttt{8.3.d.a}\otimes\chi_{-4},}$$

weight $3$, level $32$, nebentypus $\chi_{-8}$, CM by $\mathbf Q(\sqrt{-2})$,
$q$-expansion $q+2q^3-5q^9-14q^{11}+2q^{17}+34q^{19}+\dots$ (all $a_n$ with
$n$ even vanish).  On the split primes the twisting character is
$\chi_{-4}\equiv\chi_8$ (they agree on $p\equiv1,3\bmod 8$, and
$g\otimes\chi_{-8}=g$).

| $p$ | $5$ | $7$ | $11$ | $13$ | $17$ | $19$ | $23$ | $29$ | $31$ | $37$ | $41$ | $43$ | $47$ | $53$ | $59$ | $61$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $a_p$ (surface) | $0$ | $0$ | $-14$ | $0$ | $2$ | $34$ | $0$ | $0$ | $0$ | $0$ | $-46$ | $-14$ | $0$ | $0$ | $82$ | $0$ |
| `32.3.d.a` | $0$ | $0$ | $-14$ | $0$ | $2$ | $34$ | $0$ | $0$ | $0$ | $0$ | $-46$ | $-14$ | $0$ | $0$ | $82$ | $0$ |

| $p$ | $67$ | $71$ | $73$ | $79$ | $83$ | $89$ | $97$ | $101$ | $103$ | $107$ | $109$ | $113$ | $127$ | $131$ | $137$ | $139$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $a_p$ | $-62$ | $0$ | $-142$ | $0$ | $-158$ | $146$ | $-94$ | $0$ | $0$ | $178$ | $0$ | $98$ | $0$ | $-62$ | $-238$ | $-206$ |

$a_p=0$ exactly at the primes inert in $\mathbf Q(\sqrt{-2})$
($p\equiv5,7\bmod8$), as a CM form by $\mathbf Q(\sqrt{-2})$ must.
By Livné/Shioda–Inose this **identifies the motive**: $T(\mathcal X)$ is the
motive of $g$.

---

## 4. The periods

### 4.1 What "the Apéry limit" means here

The operator (dividing $L=\theta^2-tP(\theta)+t^2Q(\theta+1)-t^3R(\theta+2)$ by
$t$) is **self-adjoint**,
$$\bigl(t\,R(t)\,y'\bigr)'-4(1-10t+15t^2)\,y=g,\qquad
R(t)=1-11t+37t^2-27t^3=(1-t)(1-10t+27t^2),$$
with $g=0$ for $A(t)=\sum u_nt^n$ and $g=1$ for the companion
$B(t)=\sum b_nt^n$ ($b_0=0$, $b_1=1$, $b_2=\tfrac{13}2$; $d_n^2b_n\in\mathbf Z$,
$k=2$).  Every finite singular point $0,\ \tfrac{5\pm i\sqrt2}{27},\ 1$ has
exponents $(0,0)$: a **logarithmic (cusp) fold** in the sense of
`THEOREM_B_EXACT.md` §1, not a square-root fold.  For a fold $t_c$ put

$$A=c_0w_0+c_1\bigl(w_0\log(t-t_c)+h\bigr),\qquad
  B=p+d_0w_0+d_1\bigl(w_0\log(t-t_c)+h\bigr),\qquad
  \xi(t_c):=\frac{d_1}{c_1},$$

$w_0$ the analytic Frobenius solution, $p$ the analytic particular solution.
$\xi$ is independent of the branch and of the normalisation of $w_1$.  When
$t_c$ is the unique nearest singularity this is $\lim b_n/a_n$; here the nearest
singularities are the **conjugate pair** $\tfrac{5\mp i\sqrt2}{27}$
($|t|=1/\sqrt{27}$), the archimedean limit does not exist, and $\xi$ is the
complex fold constant of `THEOREM_B_EXACT.md` §5.

### 4.2 Method (b): analytic continuation of the ODE

`10_fold.py`, `11_fold2.py` (mpmath, $230$ digits).  $A,A',B,B'$ are summed at
$t_1$ inside $|t|<1/\sqrt{27}$ from $2200$ exact Taylor coefficients, then
continued by adaptive Taylor stepping (step $=\tfrac12\,\mathrm{dist}$ to the
nearest singularity) to a point $t_*$ inside the disc of the local Frobenius
basis at $t_c$, and matched.  Residual check of the ODE at $t_1$: $10^{-220}$.
The two $I_3$ folds are reached through the lower/upper half plane; $t=1$ along
the **real axis** (the segment $(0,1)$ is free of singularities, so $A,B$ stay
real there).

$$c_1\bigl(\tfrac{5-i\sqrt2}{27}\bigr)=\frac1\pi\Bigl(-\frac{7}{2\sqrt3}
+\frac{4i}{\sqrt6}\Bigr),\qquad
c_1(1)=-\frac{2}{\pi\sqrt3}$$
(both verified to $>170$ digits, `20_c1.gp`) — the log-coefficients of $A$ are
algebraic multiples of $1/\pi$; all the transcendence sits in the ratio.

**Results** ($\ge165$ digits; first $100$ shown):

$$\operatorname{Re}\xi_\mp=0.5287494555957668337759580638598658713411539297330525367243832947211294834920535857738431390794928\ldots$$
$$\bigl|\operatorname{Im}\xi_\mp\bigr|=0.3738823256004620250275642676907380977439755023728693417149634112084312184972757926014709015389678\ldots$$
$$\xi(1)=1.0574989111915336675519161277197317426823078594661050734487665894422589669841071715476862781589857\ldots$$

where $\xi_\mp:=\xi\bigl(\tfrac{5\mp i\sqrt2}{27}\bigr)$ are complex
conjugates, $\operatorname{Im}\xi_-<0$, and exactly
$\xi(1)=\xi_++\xi_-=2\operatorname{Re}\xi$.

### 4.3 Method (a): from the recurrence alone

`15_check.py`.  If $\xi=X-iY$ then $a_n\simeq2\operatorname{Re}(\alpha\varphi_n)$,
$b_n\simeq2\operatorname{Re}(\xi\alpha\varphi_n)$ with $\varphi$ the exact
$\mu^n/n$-solution, $\mu=5+i\sqrt2$; hence
$z_n:=a_n+i\,(b_n-Xa_n)/Y$ must satisfy $z_{n+1}/z_n\to\mu$.  Measured:

| $n$ | $100$ | $200$ | $300$ | $400$ | $500$ | $590$ |
|---|---|---|---|---|---|---|
| $\bigl|z_{n+1}/z_n-\mu\bigr|$ | $5.12\cdot10^{-2}$ | $2.58\cdot10^{-2}$ | $1.72\cdot10^{-2}$ | $1.29\cdot10^{-2}$ | $1.04\cdot10^{-2}$ | $8.79\cdot10^{-3}$ |
| $|\mu|/n$ | $5.20\cdot10^{-2}$ | $2.60\cdot10^{-2}$ | $1.73\cdot10^{-2}$ | $1.30\cdot10^{-2}$ | $1.04\cdot10^{-2}$ | $8.81\cdot10^{-3}$ |

exactly the $1/n$ of $\varphi_{n+1}/\varphi_n=\mu\,n/(n+1)$; the opposite sign
of $i$ gives $2.82$ and does not converge.  This is an independent, purely
recurrence-side confirmation that $(X,Y)$ is right (to the $\sim3$ digits the
$1/n$ rate allows); the $165$-digit confirmation is §5.

---

## 5. Identification

`12_lvalues.gp`, `13_verify.gp`, `17_cs.gp`, `22_final.gp`.  With $g$ the
newform `32.3.d.a` and $L(g,s)$ computed independently by `lfunmf`/`lfun`:

$$\boxed{\
\xi\Bigl(\tfrac{5\mp i\sqrt2}{27}\Bigr)=\frac{\sqrt2\mp i}{3}\,L(g,2),
\qquad
\xi(1)=\frac{2\sqrt2}{3}\,L(g,2),
\qquad
\bigl|\xi_\mp\bigr|=\frac{L(g,2)}{\sqrt3}. }$$

Residuals (`13_verify.gp`):
$$\Bigl|\operatorname{Im}\xi-\tfrac13L(g,2)\Bigr|=3.5\cdot10^{-166},\quad
\Bigl|\operatorname{Re}\xi-\tfrac{\sqrt2}{3}L(g,2)\Bigr|=2.5\cdot10^{-167},\quad
\Bigl|\xi(1)-\tfrac{2\sqrt2}{3}L(g,2)\Bigr|=2.0\cdot10^{-165}.$$

The two sides are computed by disjoint machinery (ODE continuation of the
recurrence vs. the $L$-function of a modular form), so this is a genuine
$165$-digit verification, not a fit.

**Chowla–Selberg.**  $g$ has CM by $\mathbf Q(\sqrt{-2})$ ($d=-8$, $h=1$,
$w=2$, $\chi_{-8}(1)=\chi_{-8}(3)=1$, $\chi_{-8}(5)=\chi_{-8}(7)=-1$), and

$$\boxed{\ L(g,2)=\frac{\pi}{32}\cdot
\frac{\Gamma\bigl(\tfrac18\bigr)\Gamma\bigl(\tfrac38\bigr)}
     {\Gamma\bigl(\tfrac58\bigr)\Gamma\bigl(\tfrac78\bigr)}\ }
\qquad\bigl[\text{residual }0\text{ at }230\text{ digits}\bigr],$$

$\Gamma(\tfrac18)\Gamma(\tfrac38)/\Gamma(\tfrac58)\Gamma(\tfrac78)
=11.42500228876934711347817240908676707711697799870160023651283605559941\ldots$
Equivalently, in closed form,

$$\operatorname{Re}\xi=\frac{\sqrt2\,\pi}{96}\,
\frac{\Gamma(\tfrac18)\Gamma(\tfrac38)}{\Gamma(\tfrac58)\Gamma(\tfrac78)},
\qquad
\operatorname{Im}\xi=\mp\frac{\pi}{96}\,
\frac{\Gamma(\tfrac18)\Gamma(\tfrac38)}{\Gamma(\tfrac58)\Gamma(\tfrac78)},
\qquad
\xi(1)=\frac{\sqrt2\,\pi}{48}\,
\frac{\Gamma(\tfrac18)\Gamma(\tfrac38)}{\Gamma(\tfrac58)\Gamma(\tfrac78)} .$$

The functional equation ($N=32$, $k=3$, root number $+1$) gives the other
critical value $L(g,1)=\tfrac{2\sqrt2}{\pi}L(g,2)=\tfrac{\sqrt2}{16}\,
\Gamma(\tfrac18)\Gamma(\tfrac38)/\Gamma(\tfrac58)\Gamma(\tfrac78)$ (verified to
$160$ digits), so equivalently $\xi(1)=\tfrac{\pi}{3}L(g,1)$ and
$\xi\bigl(\tfrac{5\mp i\sqrt2}{27}\bigr)=\tfrac{\pi}{6}(1\mp\tfrac{i}{\sqrt2})L(g,1)$.

**Bases that were tried and failed** (`17_cs.gp`, `12_lvalues.gp`): $\zeta(2)$,
$\pi^2$, $\pi^3$, $\pi^2\sqrt D$ for $D\mid24$, $\pi\,L(\chi,1)$,
$L(g,2)/\pi^e$ for $e\le4$, $\omega^4$ and $\pi^2\omega^2$ with
$\omega^2=\Gamma(\tfrac18)\Gamma(\tfrac38)/\Gamma(\tfrac58)\Gamma(\tfrac78)$ —
`lindep` returns only garbage vectors for these; the single hit is
$L(g,2)=\pi\omega^2/32$.

**Structural reading** (cf. `THEOREM_B_EXACT.md` §5, remarks 1–2).  For the
nine real-fold rows $\xi_\infty=L(\Phi,w+1)$ with $\Phi$ an *Eisenstein* series
of weight $3$ or $4$; for the three complex-fold rows
$\operatorname{Re}\xi=L(\Phi,w+1)$ and $\operatorname{Im}\xi\in\pi\mathbf Q\cdot
L(\psi,w)$.  Here $w=1$, the source is the weight-$3$ **cusp** form $g$, and the
pattern is cleaner than in either case:
$$\operatorname{Re}\xi=\sqrt2\cdot\operatorname{Im}|\xi|,\qquad
\xi_\mp=\tfrac13(\sqrt2\mp i)L(g,2),$$
i.e. $\xi_\mp$ is $L(g,2)/3$ times a generator of the order
$\mathbf Z[\sqrt{-2}]$ of norm $3$ — and $3$ is exactly the order of the
Mordell–Weil torsion, and $\tfrac{5\mp i\sqrt2}{27}$ has $27=3^3$ in the
denominator.  (Not explained here; recorded.)

---

## 6. $p$-adic

`16_padic.py`, `18_padic2.py`, `19_padic3.py`, exact rational $a_n,b_n$ to
$n=5000$.

* **$p=2$:** along $n=a\cdot2^s$ for $a=1,3,5,7,9,11$,
  $v_2(b_n/a_n)$ increases by **exactly $1$** per $s$ once $s$ is large, and
  $v_2\bigl(r_{a2^{s+1}}-r_{a2^s}\bigr)$ likewise.  Hence $b_n/a_n\to0$
  $2$-adically for every $a$:
  $$\boxed{\ \xi_2=0,\qquad v_2\bigl(b_n/a_n\bigr)=\log_2 n-c_a+O(1).}$$
  This is the **cuspidal clause** of the Euler-factor criterion
  (`EULER_CRITERION.md` §0: a cuspidal source has $\kappa_p=0$, hence
  $\xi_p=0$ whenever the limit exists) — confirmed, and it is the first row of
  the programme where that clause is actually tested.
  Note the rate is *logarithmic*, not linear: the slope
  $\sigma_2=\lim v_2(b_n/a_n)/n$ is $0$, so the row does **not** satisfy the
  positive-slope hypothesis (H1)/(d) of Theorem F, and Proposition C
  ($\sigma_p=v_p(c)+2\kappa_p$, which would predict $\sigma_2=v_2(4)=2$) fails
  for it.  Recorded as an honest exception: Proposition C is a statement about
  the Eisenstein census, not about cuspidal (K3) rows.
* **$p=3$:** $v_3(b_n)=-2\lfloor\log_3n\rfloor$ exactly (the $d_n^2$ denominator
  is sharp at $3$ for every $n$ tested), $v_3(a_n)=0$, so
  $v_3(b_n/a_n)\to-\infty$ and $b_n/a_n$ **diverges** $3$-adically.  No
  $3$-adic limit.
* **$p=5,7,11,13$:** $v_p(r_{ap^{s+1}}-r_{ap^s})$ stays $\le0$ and does not
  grow; no $p$-adic limit.  (Expected: $p\nmid$ the level $32$ and $p\nmid c=4$.)

So the adelic picture is one-sided: an archimedean fold constant that is a
critical value of a weight-$3$ cusp form, and a single $p$-adic limit, $\xi_2=0$,
at logarithmic rate.

---

## 7. Honest labels

| statement | status |
|---|---|
| $c_4,c_6,\Delta$ and $c_4^3-c_6^2=1728\Delta$ | **exact** over $\mathbf Q[t]$ |
| $c_4=E_4(Q)/y_0^4$, $c_6=E_6(Q)/y_0^6$ (twist fixed) | **verified** to $O(t^{41})$ exactly over $\mathbf Q$ |
| fibre types $I_4I_6I_3I_3IV^*$ | **proved** (Tate valuations from the model) |
| $\rho=20$, $\operatorname{MW}$ rank $0$ | **proved** (Shioda–Tate + $\rho\le h^{1,1}$) |
| $\operatorname{MW}_{\rm tors}=\mathbf Z/3$ | **proved** (explicit section; the exclusions use only standard lattice/Tate facts) |
| $\det T=72$ | **proved** given the above |
| $T\cong\operatorname{diag}(6,12)$ | **proved modulo** the discriminant-form bookkeeping of §2 (finite, but done by hand here, not by machine) |
| $a_p$ formula of §3.1 | standard (Grothendieck–Lefschetz); the two skyscraper contributions cancel |
| $a_p=a_p(\texttt{32.3.d.a})$, $5\le p\le397$ | **verified**, $76$ primes, no free parameter.  Modularity itself is Livné's theorem |
| $\xi_\mp=\tfrac13(\sqrt2\mp i)L(g,2)$, $\xi(1)=\tfrac{2\sqrt2}3L(g,2)$ | **verified** to $1.7\cdot10^{-165}$ (two disjoint computations) |
| $L(g,2)=\tfrac{\pi}{32}\Gamma(\tfrac18)\Gamma(\tfrac38)/\Gamma(\tfrac58)\Gamma(\tfrac78)$ | **verified** to $230$ digits; it is Chowla–Selberg + Deligne for a CM form, so it should be provable outright |
| $\xi_2=0$ | **verified** to $v_2=10$ along six residue classes; it is the predicted cuspidal value |
| no $3$-adic limit | **verified** ($v_3\to-\infty$ exactly) |

**Not done here.** (i) A direct point count $\#\mathcal X(\mathbf F_p)=1+p^2+p\tau_p+a_p$
with the Kodaira fibres resolved and the Galois action on the $20$ $\operatorname{NS}$
classes tracked — the trace identity of §3.1 makes it unnecessary, but it would
be an independent check of $\rho=20$ over $\mathbf F_p$.  (ii) The modular
source $\Phi$ itself: $\Theta=B(t(q))/y_0(t(q))$ should be the Eichler integral
of a weight-$3$ form on the index-$16$ genus-$0$ group of signature
$(0;3;\text{cusps }3,3,4,6)$, and §5 says its $L$-value at $s=2$ is $L(g,2)$ up
to $\tfrac13$; the direct identification $\Phi\leftrightarrow g$ (which
necessarily involves the width-$4$ cusp and the $\gamma=-1$ half-shift
$q^4=-e^{2\pi i\tau}$) is **open**.  (iii) Whether the group is congruence.

---

## 8. Scripts

`lattice/k3_period/`

| file | what |
|---|---|
| `01_jmap.gp` | recovers $U,V$ (uses `../four_term/06_jtest.gp`, `13_jdetail.gp`) |
| `02_fac.gp` | factors $V$, $U$, $U-1728V$ |
| `03_model.gp` | the model; $c_4=E_4(Q)/y_0^4$, $c_6=E_6(Q)/y_0^6$ |
| `04_torsion.gp`, `04b.gp` | $\psi_3$, the $3$-torsion section |
| `05_ap.gp`, `21_ap_ext.gp` | $a_p$ by the character sum, $p\le200$ / $p\le400$ |
| `06_newform.gp`, `07_label.gp`, `13_verify.gp`, `23_apcheck.gp` | newform search and match |
| `10_fold.py`, `11_fold2.py` | ODE analytic continuation, fold constants |
| `15_check.py` | recurrence-side confirmation |
| `12_lvalues.gp`, `17_cs.gp`, `20_c1.gp`, `22_final.gp` | $L$-values, Chowla–Selberg, `lindep` |
| `16_padic.py`, `18_padic2.py`, `19_padic3.py` | $p$-adic |
| `out/` | logs, `xi_dps230.txt`, `ap_ext.log`, `jmap.txt` |
