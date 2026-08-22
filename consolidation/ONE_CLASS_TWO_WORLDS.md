# One class, two worlds: a hypergeometric realisation of the $L(2,\chi_{-3})$ Eisenstein class at $p=3$

*Claude (Fable), 2026-08-22.  Scripts: `lattice/two_worlds/`.  Tags: **[proved]** = exact
identity or exact closed form verified symbolically; **[verified]** = exact rational/$p$-adic
computation over a stated range, no floating point in the arithmetic; **[measured]** = numerical
rate, not converged to a proof; **[conj]**/**[open]** as usual.  No irrationality is claimed
anywhere; see `CATALAN_AUDIT.md` and §7.*

---

## 0. Verdict

1. **The missing decayer for $L(2,\chi_{-3})$ exists, and it is hypergeometric.**  There is an
   explicit one-parameter family of rational functions with double poles at $\mathbb Z+\tfrac13$
   and $\mathbb Z+\tfrac23$, antisymmetric under $t\mapsto-t-b-1$, whose residue sums are **exact
   two-term linear forms** $Q_m\,L(2,\chi_{-3})-P_m$ with $Q_m,P_m\in\mathbb Q$ — no $\zeta(2)$,
   no $\pi/\sqrt3$, no $\log3$.  **[proved]** (Theorem 1, §2).  This is the $\chi_{-3}$ analogue
   of Zudilin's Catalan integral: conductor $4\leftrightarrow$ half-integer poles with weight
   $(-1)^t$; conductor $3\leftrightarrow$ third-integer poles with weight $1$ (§2.3).

2. **The rows realise the same Eisenstein extension class as the modular rows $\mathbf B$,
   $\mathbf C$, $\mathbf F$, $s_{18}$.**  Theorem F predicts $\xi_3=\Lambda_3=\zeta_3(2)$ (the
   $3$-Euler factor is $1$, so $r_3=r_\infty=1$).  Measured:
   $v_3\!\big(P_m/Q_m-\zeta_3(2)\big)=3m+O(1)$ for the $\alpha=2$ member,
   $\to\infty$ for **every** member tested.  **[verified, exact, $m\le70$, $\ge3^{199}$]** (§4).
   The cross-determinant against row $\mathbf C$ at the balanced sampling has the predicted
   divisibility rate to $1\%$ **[verified, $s\le12$]** (§4.3).

3. **$\kappa_3>0$, from factorial integrality, exactly as `THEORY_NOTES_03` §5 demands.**  For
   the distinguished member ($a=2m$, $b=m$) the coefficient row is a closed hypergeometric term
   $$Q_m=27\,(-256)^m\,\frac{(2m)!\;(1/6)_m^2\,(5/6)_m^2}{(m!)^4\,(2/3)_m\,(4/3)_m},\qquad
     v_3(Q_m)=-3m+2s_3(m)-\tfrac12 s_3(2m)+3,$$
   with $\operatorname{den}(Q_m)$ a **pure power of $3$**.  **[proved]** closed form,
   **[verified]** $m\le25$ for the valuation law, $m\le100$ for purity.  So $\kappa_3=3$, i.e.
   $\kappa_p\log p=3.296$ against Zudilin's $\kappa_2\log2=2.773$: the largest denominator
   resource in the project.

4. **But $\delta<1$.**  Evaluating the master formula (`05_two_row.tex` §5.2, $\kappa$-corrected)
   over the whole family against engines $\mathbf C,\mathbf B,\mathbf F$ gives a best value
   $$\boxed{\ \delta_{\max}=0.9191\ \ (\alpha=5/3\ \text{decayer}\times\mathbf C\ \text{engine},\ r=4.15)\ }$$
   **[measured]**, marginally above the Catalan record $0.9025$ and below $1$.  The obstruction is
   *not* the one the census predicted (no decayer): it is that the $D$-type denominators of
   $P_m$ grow faster than the linear form decays, and the net $3$-adic resource
   $(w+\kappa_3)\log3$ — where $w=\lim v_3(\mathrm{Cas}_m)/m$ — is only $\approx2$ per index,
   because *the Casoratian inherits the $3$-power denominators*.  §5, §6.

5. **A sharp structural statement, and it is the interesting one.**  For the distinguished
   member $\sigma_3=\kappa_3=3$ exactly, so the *net* $3$-adic resource is $(\sigma_3-\kappa_3)\log3=0$:
   the row pays for its own denominators.  Zudilin's Catalan row escapes this because
   $\sigma_2=2\kappa_2$ (its Casoratian has bounded $v_2$, so $w=0$).  **The design
   criterion of `THEORY_NOTES_03` §5 must be strengthened**: not "$\kappa_p>0$" but
   $$\boxed{\ w+\kappa_p>0,\qquad w:=\lim_m \tfrac1m v_p(\mathrm{Cas}_m)\ }$$
   equivalently $\sigma_p>\kappa_p$.  §6.1.  This is a *new* invariant of the realisation and it
   is the one that decides.

6. **Target (2), the weight-three cusp row $L(f,2)$: structurally out of reach for this method.**
   The residue mechanism of Theorem 1 produces exactly the $L$-values of *periodic* arithmetic
   functions, i.e. Dirichlet $L$-values = Eisenstein regulators.  A cusp-form $L$-value is not of
   that shape, and consistently Theorem F's cuspidal clause gives $\xi_p=0$.  §8.  **[argued, not
   a theorem]**  Target (3), $\beta(4)$ at $p=2$: the same machine works with **order-four** poles
   at $\mathbb Z+\tfrac12$ and three linear conditions (Rivoal–Zudilin territory); not attempted
   here.  §9.

---

## 1. What was asked and what is new

`THEORY_NOTES_06` §5(3) asks for the cellular-world realisation of a Dirichlet Eisenstein class
whose $p$-adic limit is then *forced* by regulator functoriality, and which decays.
`THEORY_NOTES_03` §5 identifies the arithmetic requirement: $\kappa_p>0$, which only
factorial/hypergeometric normalisation supplies.  `05_two_row.tex` §5.4(c) records the negative:
for $L(2,\chi_{-3})$ four modular rows align at $3$ and **none decays**.

What is new here: the decayer, constructed; the alignment, verified; and the reason the resulting
$\delta$ is still $<1$, which is a different reason from the one anticipated.

---

## 2. Theorem 1: two-term $\chi_{-3}$ linear forms from well-poised residues

### 2.1 Statement

Fix $b\ge0$ and set
$$D_b(t)=\prod_{j=0}^{b}\Bigl(t+j+\tfrac13\Bigr)^2\Bigl(t+j+\tfrac23\Bigr)^2\in\mathbb Q[t],
\qquad \deg D_b=4b+4 .$$

> **Theorem 1.** **[proved]**  Let $N\in\mathbb Q[t]$ satisfy
> $$\deg N\le 4b+2\qquad\text{and}\qquad N(-t-b-1)=-N(t),$$
> and put $R=N/D_b$.  Then $\sum_{t\ge0}R(t)$ converges and
> $$\sum_{t\ge0}R(t)\;=\;Q\,L(2,\chi_{-3})\;-\;P,\qquad Q,P\in\mathbb Q,$$
> with $Q=9\sum_{j=0}^{b}A_j$, $A_j=\lim_{t\to-(j+1/3)}\bigl(t+j+\tfrac13\bigr)^2R(t)$.
> No $\zeta(2)$, no $\pi/\sqrt3$ and no $\log3$ occurs.

### 2.2 Proof

Write the partial fraction expansion $R=\sum_i\bigl[A_i(t-p_i)^{-2}+B_i(t-p_i)^{-1}\bigr]$ over
the $2b+2$ poles $p_i\in\{-(j+\tfrac13)\}\cup\{-(j+\tfrac23)\}$.  Since $\deg N\le\deg D_b-2$ we
have $\sum_iB_i=0$, so $\sum_{t\ge0}R(t)$ converges and equals
$$\sum_i A_i\,\zeta(2,-p_i)\;-\;\sum_i B_i\,\psi(-p_i).$$
For $x=j+f$ with $f\in\{\tfrac13,\tfrac23\}$, $\zeta(2,x)=\zeta(2,f)-\sum_{l<j}(l+f)^{-2}$ and
$\psi(x)=\psi(f)+\sum_{l<j}(l+f)^{-1}$, the correction terms being rational.  Using
$$\zeta(2,\tfrac13)=4\zeta(2)+\tfrac92L,\quad \zeta(2,\tfrac23)=4\zeta(2)-\tfrac92L,\quad
\psi(\tfrac13)-\psi(\tfrac23)=-\tfrac{\pi}{\sqrt3},\qquad L:=L(2,\chi_{-3}),$$
and $\sum_iB_i=0$ (which kills $\gamma$ and the common $-\tfrac32\log3$ in both $\psi$'s), one gets
$$\sum_{t\ge0}R(t)=4\bigl(\Sigma A^{1/3}+\Sigma A^{2/3}\bigr)\zeta(2)
 +\tfrac92\bigl(\Sigma A^{1/3}-\Sigma A^{2/3}\bigr)L
 +\tfrac{\pi}{\sqrt3}\,\Sigma B^{1/3}\;+\;(\text{rational}).$$
The involution $\sigma(t)=-t-b-1$ maps $-(j+\tfrac13)\mapsto -((b-j)+\tfrac23)$ and so permutes
the pole set, exchanging the two classes.  From $R\circ\sigma=-R$ one reads off, comparing Laurent
tails at $p$ and $\sigma(p)$ (with $s\mapsto-s$),
$$A_{\sigma(i)}=-A_i,\qquad B_{\sigma(i)}=+B_i .$$
Hence $\Sigma A^{2/3}=-\Sigma A^{1/3}$ (the $\zeta(2)$ coefficient vanishes) and
$\Sigma B^{1/3}=\Sigma B^{2/3}=\tfrac12\sum_iB_i=0$ (the $\pi/\sqrt3$ coefficient vanishes).
What survives is $9\,\Sigma A^{1/3}\cdot L$ minus a rational.  $\square$

Implementation: `lattice/two_worlds/chi3_forms.py`, `chi3_row.py`.  The two vanishings are
*asserted at runtime* in the code and hold for every instance computed (thousands).

### 2.3 The dictionary — why this is Zudilin's construction for another conductor

Zudilin's Catalan row (`papers/0201024v3.pdf`, Thm. 1 and Thm. 3) is
$\sum_{t\ge0}(-1)^tR_n(t)$ with $R_n$ having order-three poles at $t\in-\tfrac12-\mathbb Z_{\ge0}$,
equivalently $\iint x^{n-1/2}(1-x)^ny^n(1-y)^{n-1/2}(1-xy)^{-n-1}$.  Half-integer poles plus the
weight $(-1)^t$ is exactly $\chi_{-4}$: with $n=2t+2k+1$, $(-1)^t=\pm\chi_{-4}(n)$.  The general
rule the two cases share:

| conductor $F$ | pole classes | weight $\omega(t)$ | value produced |
|---|---|---|---|
| $4$ (Zudilin) | $\tfrac12+\mathbb Z$ | $(-1)^t$ | $\beta(j)=L(j,\chi_{-4})$ |
| $3$ (here) | $\tfrac13+\mathbb Z$ and $\tfrac23+\mathbb Z$ | $1$ | $L(j,\chi_{-3})$ |
| general | $\tfrac aF+\mathbb Z$, $(a,F)=1$ | any $\omega$ with $\omega(t)\chi_{\rm class}=\chi$ | $L(j,\chi)$ |

and the **well-poised involution** $t\mapsto-t-c$ is the device that kills the wrong-parity
periods — $\beta(1),\beta(3)$ in Zudilin's case, $\zeta(2)$ and $\pi/\sqrt3$ in ours.  In the
language of `THEORY_NOTES_06` §1 this is the $\chi$-twisted local system on
$\mathcal M_{0,5}$: Beukers' $\zeta(2)$ kernel $1/(1-xy)$ becomes
$$L(2,\chi_{-3})=\int_0^1\!\!\int_0^1\frac{dx\,dy}{1+xy+x^2y^2},$$
the cyclotomic twist of Brown's $n=5$ cellular integral.  **[proved]** (expand the kernel).

### 2.4 The concrete family

For $0\le a\le 2b$ put
$$N_{a,b}(t)=(2t+b+1)\prod_{j=1}^{a}(t-j+1)\prod_{j=1}^{a}(t+b+j),\qquad \deg N=2a+1 ,$$
which is $\sigma$-antisymmetric (each product is mapped to $(-1)^a$ times the other).  Write
$\alpha=a/b$.  All members satisfy Theorem 1.  Note $R\ge0$ on $t\ge a$ and $R=0$ on
$0\le t<a$, so
$$\boxed{\;Q_m\,L(2,\chi_{-3})-P_m=\sum_{t\ge a}R(t)>0\ \text{strictly}\;}$$
**[proved]** — the individual linear forms never vanish (which is *not* the nonvanishing that
`CATALAN_AUDIT.md` asks for; see §7).

---

## 3. The distinguished member $\alpha=2$: closed form, recurrence, arithmetic

Take $b=m$, $a=2m$ (the maximal numerator degree).  Then **[proved]**:

$$Q_m=27\,(-256)^m\,\frac{(2m)!\,(1/6)_m^2\,(5/6)_m^2}{(m!)^4\,(2/3)_m\,(4/3)_m},\qquad
\frac{Q_{m+1}}{Q_m}=\frac{-32(2m+1)(6m+1)^2(6m+5)^2}{9(m+1)^3(3m+2)(3m+4)} .$$

$Q_0=27$, $Q_1=-300$, $Q_2=67760$, $Q_3=-6618932320/243$;
$P_0=0$, $P_1=-469/2$, $P_2=1588231/30$, $P_3=-868793149369/40824$.

* **Joint recurrence** (fitted exactly over $\mathbb Q$ from $m\le39$, `chi3_recur2.py`;
  order $2$, polynomial degree $9$, nullity $1$ — so it is *the* recurrence):
  $$c_2(m)u_{m+2}+c_1(m)u_{m+1}+c_0(m)u_m=0$$
  with $c_0=-\tfrac{32}{1345005}(2m{+}1)(3m{+}1)(6m{+}1)^2(6m{+}5)^2 g(m{+}1)$,
  $c_2=\tfrac1{149445}(m{+}2)^3(3m{+}4)(3m{+}5)(3m{+}7)g(m)$,
  $g(m)=5535m^3+11916m^2+8337m+1876$, and $c_1$ of degree $9$.
  **Characteristic polynomial $(x-1)(x+1024)$**: $\lambda_1=-1024$, $\lambda_2=1$, $c=-1024=-2^{10}$.
  **[verified]** exactly on $Q$ and $P$ for $m\le39$.
* $|Q_mL-P_m|=0.1979\,m^{-3}(1+o(1))$: **polynomial-size** linear form, $\lambda_2=1$ exactly —
  the row is the hypergeometric twin of Zagier's row $\mathbf C$ (which also has $\lambda_2=1$).
  **[verified]** $m\le64$.
* $|L-P_m/Q_m|\approx 1024^{-m}m^{-3}$; e.g. $|L-P_{64}/Q_{64}|=1.2\cdot10^{-194}$.
* $v_3(Q_m)=-3m+2s_3(m)-\tfrac12s_3(2m)+3$ **[proved from the closed form, verified $m\le25$]**,
  and $\operatorname{den}(Q_m)=3^{3m-\cdots}$ exactly (no other prime) **[verified $m\le100$]**.
  So $\kappa_3=3$, $\nu=0$.
* $\operatorname{den}(P_m)$ carries, besides the same $3$-power, a prime-to-$3$ part of measured
  rate $k_d\approx1.9$–$2.1$ **[measured $m\le100$]**; the provable bound from the tails is
  $D_{3m+2}^2$, rate $6$, so the truth is $3\times$ better than anything currently proved.
* $v_2(Q_m)>0$; **no other prime has a slope**: increments of $v_p(P_m/Q_m-P_{m-1}/Q_{m-1})$ are
  bounded and negative for $p=2,5,7,11,13$ **[verified $m\le25$]** (`chi3_other_primes.py`).

---

## 4. The alignment: one class, two worlds — verified

### 4.1 What Theorem F predicts

Row $\mathbf C$ has $\xi_\infty=\tfrac12L(2,\chi_{-3})$, source $E^{\chi_{-3},\mathbf 1}_3$,
$\xi_3=\tfrac12\zeta_3(2)$ (`EULER_CRITERION.md` §4.1).  Our row has $\xi_\infty=L(2,\chi_{-3})$,
i.e. $r_\infty=1$.  The $3$-Euler factor is $\mathcal E_3(2)=1-\chi_{-3}(3)3^{-2}=1$, so
$r_3=r_\infty/\mathcal E_3(2)=1$ and the prediction is
$$\xi_3=\Lambda_3=L_3(2,\chi_{-3}\omega^{-1})=\zeta_3(2)=2\,\xi_3^{\mathbf C}.$$

### 4.2 What is measured

`chi3_padic.py`, exact rational arithmetic, $\zeta_3(2)$ taken from row $\mathbf C$ at $n=199$
(good to $3^{380}$):

| $m$ | 10 | 20 | 30 | 40 | 50 | 60 | 70 |
|---|---|---|---|---|---|---|---|
| $v_3\bigl(P_m/Q_m-\zeta_3(2)\bigr)$ | 26 | 53 | 88 | 111 | 139 | 173 | 199 |

i.e. $\sigma_3=3$ **[verified]**.  The same test is positive for **all ten** members
$\alpha\in\{1/2,3/4,1,5/4,4/3,3/2,5/3,7/4,9/5,2\}$; at $k=14$ the valuations are

| $\alpha$ | $1/2$ | $3/4$ | $1$ | $5/4$ | $4/3$ | $3/2$ | $5/3$ | $7/4$ | $9/5$ | $2$ |
|---|---|---|---|---|---|---|---|---|---|---|
| $b$ | 28 | 56 | 14 | 56 | 42 | 28 | 42 | 56 | 70 | 14 |
| $v_3$ | 78 | 158 | 32 | 158 | 116 | 78 | 114 | 158 | 198 | 34 |

— i.e. $\approx 2.8\,b$ in every case, converging to $3b$ (the $\alpha=2$ column reaches
$199$ at $b=70$).
**This is the "one class, two worlds" statement made computational**: a hypergeometric row whose
$p$-adic limit is a Kubota–Leopoldt value it never mentions, forced by the class it realises.

### 4.3 The bridge (cross-determinant), balanced sampling

$\Delta(n,m)=2Q_mb_n-a_nP_m$ with $(a_n,b_n)$ row $\mathbf C$.  For the $\alpha=3/2$ member at
$n=27s$, $m=10s$ (which balances $\sigma_3^{\mathbf C}n=2\cdot27s$ against
$\sigma_3^{\rm dec}m=5.4\cdot10s$), `chi3_bridge.py` gives

| $s$ | 2 | 4 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|
| $v_3(\Delta)-v_3(Q_m)-v_3(a_n)$ | 105 | 211 | 321 | 423 | 535 | 643 |

i.e. $53.6\,s$ against the predicted $54\,s$ **[verified $s\le12$]**.  The hidden divisor the
correlated lattice needs is therefore there, at the predicted rate.

---

## 5. The family profile

Per unit of the row's own index $k$ (with $a=pk$, $b=qk$, $\alpha=p/q$).
$\Lambda=|\lambda_1|$, $\lambda=|\lambda_2|$ from successive ratios; $\kappa_3=-v_3(Q_k)/k$;
$w=v_3(\mathrm{Cas}_k)/k$; $\nu,k_d$ the prime-to-$3$ denominator rates of $Q_k,P_k$;
$\sigma_3=w+2\kappa_3$.  All **[measured]** at the largest $k$ reached, averaged over the last
four indices (`chi3_design2.py`).

| $\alpha$ | $k$ | $\log\Lambda$ | $\log\lambda$ | $\kappa_3$ | $w$ | $\nu$ | $k_d$ | $\sigma_3$ |
|---|---|---|---|---|---|---|---|---|
| $2$ | 29 | 6.839 | $-0.110$ | 2.743 | $-2.688$ | 0.000 | 1.880 | 2.797 |
| $9/5$ | 12 | 22.818 | $-14.104$ | 12.106 | $-10.259$ | 1.941 | 21.246 | 13.953 |
| $7/4$ | 15 | 15.883 | $-13.858$ | 9.271 | $-7.401$ | 2.456 | 18.926 | 11.141 |
| $5/3$ | 18 | 9.180 | $-13.296$ | 6.481 | $-4.681$ | 3.088 | 16.140 | 8.281 |
| $3/2$ | 25 | 2.452 | $-12.590$ | 3.531 | $-1.637$ | 3.694 | 13.831 | 5.425 |
| $4/3$ | 18 | $-1.979$ | $-24.340$ | 4.481 | $-0.738$ | 8.845 | 24.745 | 8.224 |
| $5/4$ | 15 | $-6.746$ | $-36.292$ | 5.271 | $+0.305$ | 14.860 | 36.684 | 10.847 |
| $1$ | 29 | $-3.156$ | $-10.207$ | 0.743 | $+1.192$ | 4.545 | 10.126 | 2.677 |

**A uniformity worth recording.**  Re-scaled to the index $b$ (the number of pole pairs) rather
than $k$, the $3$-adic slope is the *same for the whole family*:
$$\sigma_3 = 3\ \text{per unit of }b\qquad\text{[verified, all ten members, } b\le56].$$
Indeed $v_3\bigl(P/Q-\zeta_3(2)\bigr)$ depends only on $b$, not on $a$: the members
$\alpha=1/2$ and $3/2$ (both $b=2k$), and $3/4$ and $5/4$ (both $b=4k$), give *identical*
valuation sequences.  Consequently
$$w=3-2\kappa_3,\qquad \text{net }3\text{-adic resource}=(w+\kappa_3)\log3=(3-\kappa_3)\log3
\quad\text{per unit of }b,$$
so the whole $p$-adic side of the design problem collapses to a single number, $\kappa_3(\alpha)$,
which runs from $2.74$ at $\alpha=2$ (resource $0$) down to $0.74$ at $\alpha=1$ (resource $2.5$).
The trade-off is that $k_d$ rises from $1.9$ to $10.1$ over the same range.

(The $\alpha=2$ row is the only one with $\nu=0$, i.e. with $\operatorname{den}(Q)$ a pure
$3$-power; it is also the only one with $\lambda=1$.  For $\alpha<2$ the linear form decays
genuinely — $\lambda$ as small as $10^{-13}$ per index at $\alpha=1/2$ — at the price of
non-$3$ denominators.)

---

## 6. The master formula, and where it stops

With decayer at index $\gamma n$, engine at $\alpha_s n$, $\gamma=1$, $r=\alpha_s$:
$$F=\tfrac12\Bigl[\underbrace{\max(k_{\rm eng}r,\;k_d)}_{S}
+\underbrace{\kappa_3\log3+\nu}_{\eta}+r\log\rho_2^{\rm eng}+\log\lambda
-\underbrace{\min(\sigma_3^{\rm eng}r,\ \sigma_3^{\rm dec})\log3}_{G}\Bigr],\quad
H=F+\log\tfrac\Lambda\lambda .$$
Engines: $\mathbf C$ $(k{=}2,\sigma_3{=}2,\log\rho_2{=}0)$, $\mathbf B$
$(2,3,\tfrac12\log27)$, $\mathbf F$ $(2,2,\log8)$.  Optimising $r$:

| $\alpha$ | $\times\mathbf C$ | $\times\mathbf B$ | $\times\mathbf F$ |
|---|---|---|---|
| $2$ | $0.8410$ ($r{=}1.40$) | $0.8105$ | $0.7484$ |
| $9/5$ | $0.9128$ ($r{=}6.98$) | $0.8338$ | $0.7740$ |
| $7/4$ | $0.9158$ ($r{=}5.58$) | $0.8369$ | $0.7772$ |
| $5/3$ | $\mathbf{0.9191}$ ($r{=}4.15$) | $0.8409$ | $0.7815$ |
| $3/2$ | $0.9133$ ($r{=}2.71$) | $0.8375$ | $0.7798$ |
| $4/3$ | $0.8969$ | $0.8224$ | $0.7656$ |
| $1$ | $0.8577$ | $0.7871$ | $0.7335$ |

Row $\mathbf C$ is the best engine, exactly as `05_two_row.tex` §5.5 predicts it should be
($\log\rho_2^{\rm eng}=0$ is the cheapest lever, and $\mathbf C$ is the only row in the class
whose linear form is polynomial-size).  Best overall $\delta=0.9191<1$.

### 6.1 Why — and a correction to the design criterion

At the balanced ratio the condition $F<0$ reads
$$\log\tfrac1\lambda\;>\;\underbrace{k_d+\nu}_{\text{denominators}}\;-\;
\underbrace{(w+\kappa_3)\log3}_{\text{net }3\text{-adic resource}} .$$
Two things go wrong, and only one of them was anticipated.

**(i) The Casoratian carries the denominators too.**  `THEORY_NOTES_03` §5 and
`05_two_row.tex` eq. (kappaslope) derive $\sigma_p=v_p(c)+2\kappa_p$ from the *Zagier-normalised*
Casoratian $a_nb_{n-1}-a_{n-1}b_n=-c^{n-1}/n^2$, i.e. assuming $v_p(\mathrm{Cas}_n)=O(\log n)$.
Our rows are not in that normalisation: writing $w=\lim\frac1m v_p(\mathrm{Cas}_m)$, the correct
statement is $\sigma_p=w+2\kappa_p$, and the *net* resource after paying $\eta=\kappa_p\log p$ is
$(\sigma_p-\kappa_p)\log p=(w+\kappa_p)\log p$.  For the $\alpha=2$ member $w=-3$, $\kappa_3=3$:
**net resource exactly zero**.  The row pays for its own denominators.  Zudilin's Catalan row
has $w=0$ (its Casoratian $\mathrm{Cas}_m=(-1)^{m-1}(20m^2-8m+1)/(8m^2(2m-1)^2)$ has bounded
$v_2$), hence $\sigma_2=2\kappa_2=8$ and net resource $\kappa_2\log2=2.77$.  So:

> **Corrected design criterion.**  A decaying partner must have $w+\kappa_p>0$, i.e.
> $\sigma_p>\kappa_p$ — not merely $\kappa_p>0$.  The invariant to shop for is the $p$-adic
> valuation of the **Casoratian**, not of the coefficients.

In our family $w+\kappa_3$ rises monotonically as $\alpha$ falls ($0.0$ at $\alpha=2$, $0.9$ at
$3/2$, $2.0$ at $1$, $3.0$ at $1/2$, per unit of $b$), which is the good direction — but:

**(ii) the $D$-type denominators of $P_m$ outrun the decay.**  $k_d+\nu>\log\frac1\lambda$ for
every member: e.g. $\alpha=5/3$ has $k_d+\nu=19.2$ against $\log\frac1\lambda=13.3$.  The single-row
score is $-5.9$ and the $3$-adic mechanism returns only $(w+\kappa_3)\log3=1.98$.  This is the
classical Apéry barrier, not a $p$-adic one.  Note that $k_d$ is *measured*: the provable bound
from the partial-fraction tails is $D_{3b+2}^2$, rate $6$ per unit of $b$, i.e. $3\times$ worse
than the truth at $\alpha=2$.  **Closing that gap is the single most valuable next step**
(Rhin–Viola's group actions on parameters are exactly the technology that removes such primes).

---

## 7. What is proved / verified / conjectural

* **[proved]** Theorem 1 (§2.2); the exact closed form of $Q_m$ and its ratio; the $3$-adic
  valuation law $v_3(Q_m)=-3m+2s_3(m)-\frac12s_3(2m)+3$; the kernel identity
  $L(2,\chi_{-3})=\iint(1+xy+x^2y^2)^{-1}$; strict positivity $Q_mL-P_m>0$.
* **[verified]** (exact arithmetic, stated ranges) the order-$2$ recurrence with roots
  $-1024,1$ ($m\le39$); $\operatorname{den}(Q_m)$ a pure $3$-power ($m\le100$);
  the $3$-adic limit $\xi_3=\zeta_3(2)$ ($m\le70$, to $3^{199}$); the absence of slopes at
  $p=2,5,7,11,13$; the bridge divisibility ($s\le12$).
* **[measured]** every rate in §5, and hence the whole of §6.  $\lambda$, $\kappa$, $w$, $\nu$,
  $k_d$ are read off at finite index with visible $O(\log m/m)$ drift; the $\delta$ column should
  be read to two decimals at best.
* **[open]** the recurrences for $\alpha\ne2$ (the exact joint $(Q,P)$ nullspace fit over
  $\mathbb Q$ was cut off on time before reaching order $>3$; the $\alpha=2$ case is settled); a *proved* denominator
  bound for $P_m$; the exact asymptotics $\Lambda(\alpha),\lambda(\alpha)$ (the naive saddle-point
  of `chi3_asym.py` is wrong by a constant because $Q_m=9\sum_jA_j$ has a sign change at
  $j\approx b/2$ and the residues cancel).
* **No irrationality is claimed.**  $\delta<1$ everywhere, so the question does not arise; and
  even $\delta>1$ would be void without the nonvanishing of the *selected* linear form
  (`CATALAN_AUDIT.md` §6, `05_two_row.tex` Remark on $F>0$).  The strict positivity of §2.4 is
  nonvanishing for the *rows*, which is not the same statement.

---

## 8. Target (2): the weight-three cusp row $L(f,2)$, $f=\eta_2^3\eta_6^3$

The engine parameters are $k=2$, $\sigma_2=2$, $\log\rho_2^{\rm eng}=\log4$, and
`05_two_row.tex` §5.5 asks for a decayer with $\Lambda_{\rm dec}>14.8,\,29.6,\,160$ according as
$\sigma_{\rm dec}=1,2,3$.  I did not find one, and I believe the reason is structural, not
combinatorial:

> The residue mechanism of Theorem 1 turns $\sum_{t\ge0}R(t)$ into $\sum_i A_i\zeta(s,x_i)$, i.e.
> into $L$-values of **periodic** arithmetic functions.  Every value it can produce is a
> $\mathbb Q$-combination of $L(s,\chi)$ for Dirichlet $\chi$ — precisely the Eisenstein
> regulators.  A cusp-form $L$-value is not of that shape: $a_n(f)$ is not periodic, and no
> rational function's residues can generate it.

This matches Theorem F's cuspidal clause ($\xi_p=0$) and the entry "$\mathcal M_{0,n}$ cellular
integrals $\to$ $\zeta$-values" of `THEORY_NOTES_06` §4(a): the cellular/hypergeometric world is
coextensive with the Eisenstein world.  A hypergeometric row for $L(f,2)$ would have to come from
a different mechanism entirely (the CM structure of $f$ — $f$ has CM by $\mathbb Q(\sqrt{-3})$, so
$L(f,2)$ is a Chowla–Selberg-type period, and the natural approximations are Chudnovsky's
$\Gamma(1/3)$ constructions, which are *not* two-term forms in $1$ and $L(f,2)$).
**[argued; not a theorem, and worth a proper look before it is written into the paper]**

## 9. Target (3): $\beta(4)=L(4,\chi_{-4})$ at $p=2$

Theorem 1 generalises verbatim to pole order $\rho$: with poles of order $\rho$ at
$\tfrac aF+\mathbb Z$ and the antisymmetry $R\circ\sigma=-R$, the Laurent coefficients satisfy
$A^{(j)}_{\sigma(i)}=(-1)^{j+1}A^{(j)}_i$, so the *even*-order coefficients pair with a sign and
the *odd*-order ones without.  For $F=4$, $\rho=4$ this leaves $\beta(4)$ and $\beta(2)$ together
with $\beta(3),\beta(1)$; three linear conditions on the numerator are needed to isolate
$\beta(4)$ — one more than Zudilin needed for $\beta(2)$.  This is Rivoal–Zudilin's territory
(`[RZ]` in `papers/0201024v3.pdf`) and is the natural next construction; the level-$24$ rank-$16$
engine then has to be matched against it.  **Not attempted here.**

---

## 10. Scripts

`lattice/two_worlds/`

| file | what |
|---|---|
| `chi3_forms.py` | Theorem 1: exact partial fractions, the two vanishing assertions |
| `chi3_row.py` | fast exact generator of $(Q_m,P_m)$ for $N_{a,b}$ |
| `chi3_family.py`, `chi3_scan.py` | first scan; verification that $\sum_tR(t)$ *is* the form |
| `chi3_denoms.py` | closed form of $Q_m$ and the $v_3$ law, checked |
| `chi3_recur.py`, `chi3_recur2.py`, `chi3_recur3.py` | exact recurrence fitting over $\mathbb Q$ |
| `chi3_padic.py` | the $3$-adic alignment test against row $\mathbf C$ |
| `chi3_other_primes.py` | no slopes at $p\ne3$ |
| `chi3_bridge.py` | cross-determinant divisibility at balanced sampling |
| `chi3_full.py`, `chi3_design.py`, `chi3_design2.py` | rate profile and the master-formula table |
| `chi3_asym.py` | saddle-point rates (**wrong for $\Lambda$**, kept as a record; see §7) |
