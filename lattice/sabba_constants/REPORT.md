# Sabba's alternating continued fractions: structure, an exact theorem, and why the
# E-function/Siegel–Shidlovsky route does not apply

Working notes, 2026-09-02.  All scripts in this directory; every number quoted is
reproducible from `summary.py`, `taylor_step2.py`, `Lconst2.py`, `hunt.py`.

---

## 0. Executive summary

Let

$$C=\cfrac{1}{1+\cfrac{1}{2-\cfrac{1}{3+\cfrac{1}{4-\ddots}}}},\qquad
  D=1-\cfrac{1}{2+\cfrac{1}{3-\cfrac{1}{4+\ddots}}}$$

(OEIS A244279/A244280 and A262957/A263295, Mohamed Sabba).

**What I proved.**

1. Both constants are governed by one number each: with $y^{(a)}$ the **minimal** solution of
   $y_n=n\,y_{n-1}+a(-1)^n y_{n-2}$ and $r_a=y_1/y_0$,
   $$\boxed{\;C=\frac{1}{1-r_{+1}},\qquad D=1-r_{-1}.\;}$$
2. **GPT-5.6's formula is correct, and I can prove it** (§3), not merely verify it.
   It is exactly the statement "$F(1^-)=0$", i.e. the condition that the exponential
   generating function be entire.  Independently verified numerically to **53 digits**.
3. The **Laplace/Borel dual** of the problem is a *second-order* equation.  With
   $Y_a(z)=\sum_{n\ge0}y_n z^n$ (entire, exponential type 1) and $W=zY$:
   $$\boxed{\;W''=\Bigl(\frac1{z^4}-\frac2{z^3}-a^2\Bigr)W\;}$$
   a **doubly-confluent Heun operator**: rank-1 irregular singularities at $z=0$ and
   $z=\infty$, and *no other singular point at all*.  (§4; verified symbolically.)
4. That operator has **no Liouvillian solutions**; its differential Galois group is
   $\mathrm{SL}_2(\mathbb C)$ (Kovacic; §4.3, confirmed by brute-force Riccati search).
   It is **not Möbius-gauge-equivalent to a Bessel/Whittaker/Kummer operator** (§4.4).
5. **Consequently the 4th-order EGF operator is *not* an E-operator, and none of its
   solutions is an E-function.**  Siegel–Shidlovsky and Beukers' 2006 refinement
   therefore **do not apply** to $C$ and $D$ as stated in the task (§5).  This is the
   main negative structural result: the alternating fractions genuinely leave the
   Bessel/E-function world that contains their non-alternating relatives
   $I_0(2)/I_1(2)$ and $J_0(2)/J_1(2)$.
6. Two exact identities discovered and verified to 120 digits (§6):
   $f_a(-1)=-2a\,f_a'(1)$, and a bilinear invariant coupling $a=+1$ to $a=-1$:
   $r_{+1}+r_{-1}-2=-\bigl(f_{+}(1)f_{-}(-1)+f_{-}(1)f_{+}(-1)\bigr)$.
7. The moment kernel's extra singularities **are apparent** (§7) — the task's guess was
   right — but this does *not* produce a Laurent-coefficient closed form, because
   $t=0$ and $t=\infty$ are both irregular, so there is no annulus.

**What I did not find.**  No closed form.  All PSLQ searches negative (§8).  I now
believe there is no closed form in classical special functions, and §4–§5 say why.

**Status of transcendence.**  Open.  $C,D$ are *exponential periods* /
connection constants of an irregular rank-2 connection with $\mathrm{SL}_2$ Galois group.
Irrationality and irrationality measure $2$ are known (Bala).  No current theorem gives
transcendence; §5.3 states precisely what would be needed.

---

## 1. Set-up and the reduction to one number

Convergent numerators/denominators of both fractions satisfy
$$y_n=n\,y_{n-1}+\varepsilon_n y_{n-2},\qquad \varepsilon_n=a(-1)^n,\quad a=\pm1 .$$
The solution space is 2-dimensional, spanned by two integer sequences; the *dominant*
solutions grow like $K\,n!$ and the *minimal* solution decays like $L/(n+1)!$.
Normalise the minimal solution by $y_0=1$ and write $r_a=y_1$.

**Proposition 1.**  $C=1/(1-r_{+1})$ and $D=1-r_{-1}$.

*Proof.*  Pincherle: the tail of the continued fraction is the ratio of consecutive
terms of the minimal solution.  Numerically (140 digits, `basics.py`, `summary.py`):
$C-1/(1-r_{+1})=7.9\cdot10^{-131}$, $D-(1-r_{-1})=-1.0\cdot10^{-130}$. ∎

Minimal solutions ($y_0=1$):

| $n$ | $y_n$, $a=+1$ | $y_n$, $a=-1$ |
|---|---|---|
|0| $1$ | $1$ |
|1| $-0.5904911352531017312102$ | $0.4233666102698156076021$ |
|2| $-0.1809822705062034624205$ | $-0.1532667794603687847958$ |
|3| $0.04754432373449134394883$ | $-0.03643372811129074678524$ |
|4| $0.009195024431761913374842$ | $0.007531867015205797654815$ |

Signs: $\operatorname{sgn}y_n=(-1)^{n(n+1)/2}$ for $a=+1$ and $(-1)^{\lfloor n/2\rfloor}$
for $a=-1$ (period 4 in both cases).  Asymptotics $(n+1)!\,|y_n|\to L_a$ with a
parity-dependent $1/n$ correction:
$$L_{+1}=1.1181244000437889584649550466815851920430177\ldots$$
$$L_{-1}=0.89112302306750817285189998765302692819222666\ldots$$
(44 digits, `Lconst2.py`; both parity classes agree to that many digits.)

**High-precision values** (`summary.py`, 100 digits):
```
C   = 0.6287366070989547994355879022526302953834327469083313456011767306080768279905057674386550681519432139
D   = 0.5766333897301843923978917497829139257961494352757108398411041918076483544012440202031822641121339437
r+  = -0.59049113525310173121023444226659472466262123978317117865639
r-  =  0.423366610269815607602108250217086074203850564724289160158896
```
**OEIS erratum.**  The comments in A244279 and A244280 give
$C=0.628736607098954801603428\ldots$, which is wrong from the **16th decimal** on
(correct: $\ldots0989547994355879\ldots$).  Worth reporting to OEIS.

---

## 2. The exponential generating function and the fourth-order operator

Let $f_a(t)=\sum_{n\ge0}y_n t^n/n!$.

**Proposition 2 (verified symbolically, `egf_ode2.py`).**
$$\bigl((1-t)f(t)\bigr)''=a\,f(-t)\qquad\text{exactly}$$
(the inhomogeneous polynomial $y_0+(y_1-y_0)t$ is killed by $\partial_t^2$), and hence
$$(1-t^2)f''''-(6t+2)f'''-6f''-f=0 .$$
Equivalently, with $F=(1-t)f$, $G=(1+t)f(-t)$:
$$F''=\frac{aG}{1+t},\qquad G''=\frac{aF}{1-t},\qquad (1-t)\bigl((1+t)F''\bigr)''=F,$$
and in the variable $x=(1+t)/2$ this is **exactly GPT-5.6's operator**
$$\mathcal L\,[F]\;=\;x(1-x)F''''+2(1-x)F'''-4F\;=\;(1-x)\bigl(xF''\bigr)''-4F\;=\;0 .$$
(All four identities verified to $t^{21}$ symbolically, with $a^2=1$ imposed.)

**Note.** $a$ drops out of the operator; the 4-dimensional solution space of $\mathcal L$
is the direct sum of the 2-dimensional $a=+1$ and $a=-1$ recurrence families.

*Correction to the first pass:* the earlier hand elimination gave
$(1-t^2)f''''-(2t+6)f'''+6f''-f$, which is **wrong** (a sign slip in reflecting
$f'(-t)$).  The corrected coefficients are $-(6t+2)$ and $-6$.

### 2.1 Local exponents (`local_exponents.py`)

| point | indicial polynomial | exponents |
|---|---|---|
| $x=0$ $(t=-1)$ | $\rho(\rho-1)^2(\rho-2)$ | $0,1,1,2$ (log, repeated exponent) |
| $x=1$ $(t=1)$ | $\rho(\rho-1)(\rho-2)(\rho-3)$ | $0,1,2,3$ |
| $x=\infty$ | irregular, $f\sim e^{2\zeta\sqrt t}$, $\zeta^4=-1$ | slope $1/2$ |

At $x=1$ write $X=1-x$.  Substituting $F=\sum c_kX^k$ into $\mathcal L$, the coefficients of
$X^{-3},X^{-2},X^{-1}$ vanish identically and the $X^{0}$ coefficient is $-4c_0$.
Hence:

**Lemma 3.**  The analytic (log-free) solutions at $x=1$ are exactly those with $c_0=0$;
they form a 3-dimensional space with free $c_1,c_2,c_3$.  A solution with $c_0\ne0$
carries a $c_0\,X^3\log X$ term.  In particular every solution is *bounded* at $x=1$ and
$F(1^-)=c_0$ exists.

---

## 3. Theorem: GPT-5.6's formula, with proof

**Theorem 4.**  Let $\mathcal L_a[y]=x(1-x)y''''+2(1-x)y'''-4a^2y$ and let $P_a,Q_a$ be the
solutions with
$$(P,P',P'',P''')(\tfrac12)=(1,0,4a,-8a),\qquad (Q,Q',Q'',Q''')(\tfrac12)=(0,1,0,-4a).$$
Then
$$C=\frac{2Q_{+1}(1^-)}{P_{+1}(1^-)},\qquad D=\frac{P_{-1}(1^-)}{2Q_{-1}(1^-)} .$$

*Proof.*  For $F=(1-t)f_a$ with $x=(1+t)/2$ one computes from the recurrence
($y_2=2y_1+ay_0$, $y_3=3y_2-ay_1$):
$$(F,F_x,F_{xx},F_{xxx})(\tfrac12)=\bigl(y_0,\;2(y_1-y_0),\;4a\,y_0,\;-8a\,y_1\bigr).$$
So $P_a$ is the case $(y_0,y_1)=(1,1)$ and $Q_a$ the case $(y_0,y_1)=(0,\tfrac12)$, and by
linearity $F=y_0P_a+2(y_1-y_0)Q_a$.

Next, $f_a$ analytic at $t=1$ $\iff$ $F(1^-)=0$ (Lemma 3: $f=F/(2X)$, and $c_0=0$
simultaneously removes the pole and the log).  And $f_a$ analytic at $t=1$ $\iff$
$y$ is the minimal solution: writing $c_n=y_n/n!$, the recurrence gives
$c_n-c_{n-1}=a(-1)^n y_{n-2}/n!=O(n^{-2})$, so $K=\lim c_n$ exists for every solution;
$K\ne0$ forces a pole $K/(1-t)$, and $y\mapsto K$ is a linear functional whose kernel is
precisely the minimal line.  Therefore
$$P_a(1^-)+2(r_a-1)\,Q_a(1^-)=0\ \Longrightarrow\ r_a=1-\frac{P_a(1^-)}{2Q_a(1^-)} .$$
Combining with Proposition 1 gives the two displayed formulas. ∎

**Numerical verification** (`taylor_step2.py`: Taylor stepping from $x=1/2$ toward $1$,
then a local fit $y(1-\delta)=c_0+c_1\delta+c_2\delta^2+\sum_{k\ge3}(\alpha_k+\beta_k\log\delta)\delta^k$):
```
P_{+1}(1^-) = 1.42246348902752997481433451717184870274098754
Q_{+1}(1^-) = 0.447177433906655257608362980315149048560442248
   2Q/P - C   = -1.0e-53
P_{-1}(1^-) = 0.647086176435260596952420030588059861625205851
Q_{-1}(1^-) = 0.561089756472515530254258677928634815720195682
   P/(2Q) - D = -1.4e-53
```
**Interpretation.**  As the task anticipated, this is *not* a closed form: it is the
Pincherle/connection statement "the entire solution is the one killing $F(1)$", written
as a ratio of two connection constants of a rank-4 operator.

---

## 4. The Laplace dual: a doubly-confluent Heun equation

This is the sharpest structural result and, I think, the reason no closed form exists.

### 4.1  A second-order equation for the ordinary generating function

Let $Y_a(z)=\sum_{n\ge0}y_n z^n$ (entire of order 1, exponential type 1, since
$y_n\sim L s_n/(n+1)!$).  From $Y=z(zY)'+az^2Y(-z)$ plus the two low-order boundary terms:

**Proposition 5 (verified symbolically to $z^{25}$, `egf_ode2.py`).**
$$z^4Y''+2z^3Y'+(z^4+2z-1)\,Y\;=\;-y_0+(2y_0-y_1)z-a\,y_0z^2-a(y_0-y_1)z^3 .$$

The homogeneous operator, put in normal form by $W=zY$:
$$W''=\Bigl(\frac1{z^4}-\frac2{z^3}-1\Bigr)W\qquad(\text{general }a:\ \ -a^2\ \text{in place of}\ -1),$$
and the gauge $W=e^{-1/z}V$ gives the strikingly simple
$$V''+\frac{2}{z^2}V'+V=0 .$$
Both reductions checked by sympy.

**Singularity structure.**  As a projective connection $Q\,dz^2$ with
$Q=z^{-4}-2z^{-3}-1$: a pole of order 4 at $z=0$ (rank-1 irregular, formal solutions
$e^{1/z}z^2$ and $e^{-1/z}$) and a pole of order 4 at $z=\infty$ (rank-1 irregular,
formal solutions $e^{\pm iz}$).  **No other singular point.**  This is precisely the
class of *doubly-confluent Heun* equations (of which the modified Mathieu equation is
the best-known member).

**The whole family is a deformation of a solvable point.**  For general
$y_n=n y_{n-1}+\mu(-1)^ny_{n-2}$ the same elimination gives
$z^4Y''+2z^3Y'+(\mu^2z^4+2z-1)Y=\ldots$, i.e. $W''=(z^{-4}-2z^{-3}-\mu^2)W$.
At $\mu=0$ the equation *is* solvable — $W=e^{-1/z}$ is a solution and the second follows
by reduction of order — and $\mu$ cannot be scaled away ($a_{-2}=a_{-1}=0$ pins the
scaling).  So the problem sits at $\mu^2=1$ of a one-parameter deformation of a solvable
point: the exact structural analogue of Mathieu.

### 4.2  Fourier–Laplace duality (consistency check)

The Fourier–Laplace transform ($t\mapsto-\partial_z$, $\partial_t\mapsto z$, a Weyl-algebra
automorphism) of $\mathcal L_f=(1-t^2)\partial^4-(6t+2)\partial^3-6\partial^2-1$ is, exactly
(`fl_check.py`),
$$\mathrm{FL}(\mathcal L_f)=z^4\partial^2+2z^3\partial+(-z^4+2z^3+1),$$
of order 2 in $\partial$ and degree 4 in $z$.  Its normal form ($W=zY$) is
$$W''=\Bigl(1-\frac2z-\frac1{z^4}\Bigr)W,$$
i.e. the image of §4.1's $W''=(z^{-4}-2z^{-3}-1)W$ under $z\mapsto1/z$ (the reversal
$(a_{-4},a_{-3},a_{-2},a_{-1},a_0)\mapsto(a_0,a_{-1},a_{-2},a_{-3},a_{-4})$, with
$W(z)=z\,\psi(1/z)$).  Irregular of rank 1 at **both** $z=0$ and $z=\infty$.  So the rank-4
EGF operator and the rank-2 doubly-confluent Heun operator are Fourier–Laplace duals.
Concretely, $\widetilde Y(\sigma)=\sum_{n\ge0}y_n\sigma^{-n-1}$ (convergent for
$\sigma\ne0$) solves this equation up to the same kind of low-degree polynomial
inhomogeneity as in Proposition 5.

### 4.3  No Liouvillian solutions: Galois group $\mathrm{SL}_2(\mathbb C)$

**Theorem 6.**  $W''=(z^{-4}-2z^{-3}-1)W$ has no Liouvillian solution; its differential
Galois group is $\mathrm{SL}_2(\mathbb C)$.

*Proof (Kovacic, `structure.py`).*  $r=z^{-4}-2z^{-3}-1$ has one pole, at $0$, of order 4,
and $o(\infty)=0$.
* Case 3 needs every pole of order $\le2$: excluded (order 4).
* Case 2 needs a pole of order 2 or odd order $>2$: excluded.
* Case 1: $\sqrt r=z^{-2}-z^{-1}-\tfrac12-\ldots$ so $[\sqrt r]_0=z^{-2}-z^{-1}$, $v=2$,
  and $r-[\sqrt r]_0^2=-1-z^{-2}$ has zero coefficient at $z^{-3}$, so $b=0$ and
  $\alpha_0^\pm=\tfrac12(0+2)=1$.  At $\infty$, $[\sqrt r]_\infty=i$, and
  $r-(-1)=z^{-4}-2z^{-3}$ has zero coefficient at $z^{-1}$, so $\alpha_\infty^\pm=0$.
  Hence $d=\alpha_\infty-\alpha_0=-1<0$ for every sign choice: no case-1 solution.

A brute-force search for rational Riccati solutions
$\omega=\pm z^{-2}+\beta z^{-1}+c+P'/P$, $\deg P\le6$, returns **NONE**, confirming this.
A 2nd-order equation in normal form with no Liouvillian solution has Galois group
$\mathrm{SL}_2(\mathbb C)$. ∎

(Note: $\mathrm{SL}_2$ does **not** by itself exclude Bessel — Bessel also has $\mathrm{SL}_2$ for
generic order.  That is settled next.)

### 4.4  Not Bessel/Whittaker/Kummer

For $W''=QW$ in normal form, $Q$ is a complete invariant of the gauge class (a gauge
$W=sV$ preserving normal form forces $s$ constant), and under a change of variable
$z=\varphi(w)$ the object $Q\,dz^2$ transforms as
$\widetilde Q=Q(\varphi)\varphi'^2-\tfrac12\{\varphi,w\}$, with vanishing Schwarzian for
Möbius $\varphi$.  Pole orders of $Q\,dz^2$ are therefore Möbius invariants.

* Bessel/Whittaker/Kummer in normal form: **one** pole of order 4 (at $\infty$) and **one**
  pole of order 2 (at $0$).
* Ours: **two** poles of order 4 (at $0$ and $\infty$), none of order 2.

Hence **no Möbius change of variable plus algebraic gauge takes our equation to a
Bessel/Whittaker/Kummer equation.**  For a pullback of degree $d\ge2$: our equation has
*only* the two singular points, so $\varphi^{-1}(\{0,\infty\})\subseteq\{0,\infty\}$, forcing
$\varphi(w)=cw^{\pm d}$; but then the rank-1 irregular point at $\infty$ pulls back to rank
$d$, contradiction unless $d=1$.  This closes the classical-special-function route.

### 4.5  The Borel transform of $V$

From $V''+2z^{-2}V'+V=0$, writing $V(z)=\int e^{z\xi}\phi(\xi)\,d\xi$ gives
$(1+\xi^2)\phi''+4\xi\phi'+2(1+\xi)\phi=0$, and the substitution $\phi=\psi/(1+\xi^2)$
removes the first derivative exactly:
$$\psi''+\frac{2\xi}{1+\xi^2}\,\psi=0 .$$
Regular singular points at $\xi=\pm i$ with exponents $0,1$; the Frobenius obstruction
at the lower exponent is $K+1=0$, so $K=-1\ne0$: **the singularities are logarithmic, not
apparent.**  (Under $\sigma=(\xi+i)/2i$ this is $\psi''=2i(2\sigma-1)\sigma^{-1}(1-\sigma)^{-1}\psi$
— a confluent-Heun normal form with an irregular point of rank $1/2$ at $\infty$.)

---

## 5. Transcendence: what applies, what does not (task B)

### 5.1  The theorem the task wanted to use

**Beukers (Ann. of Math. 163 (2006), 369–379).**  Let $f_1,\dots,f_m$ be E-functions
satisfying $y'=A(z)y$ with $A\in M_m(\overline{\mathbb Q}(z))$, and let
$\xi\in\overline{\mathbb Q}$, $\xi\ne0$, not a pole of $A$.  Then every homogeneous
$P\in\overline{\mathbb Q}[X_1,\dots,X_m]$ with $P(f_1(\xi),\dots,f_m(\xi))=0$ is the
specialisation at $z=\xi$ of a homogeneous $Q\in\overline{\mathbb Q}[z][X]$ with
$Q(z,f_1(z),\dots,f_m(z))\equiv0$.

**Corollary (the form one would want here).**  If $f_1,f_2$ are E-functions, linearly
independent over $\overline{\mathbb Q}(z)$, satisfying such a system, and
$\xi\in\overline{\mathbb Q}\setminus\{0\}$ is not a pole of $A$ with $f_2(\xi)\ne0$, then
$f_1(\xi)/f_2(\xi)$ is transcendental.  *(Proof: an algebraic value $\alpha$ gives the
relation $X_1-\alpha X_2$; lifting it produces $q_1f_1+q_2f_2\equiv0$ with
$(q_1(\xi),q_2(\xi))=(1,-\alpha)$, forcing $q_1=q_2=0$ by independence — contradiction.)*

This is exactly how one gets transcendence of the **non-alternating** relatives:
$R=I_0(2)/I_1(2)$ and $J_0(2)/J_1(2)$ are ratios of values at an algebraic point of the
E-functions $I_\nu(2z),J_\nu(2z)$, which are linearly independent over
$\overline{\mathbb Q}(z)$.

### 5.2  Why it does **not** apply to $C$ and $D$

**Proposition 7.**  No nonzero solution of $\mathcal L_f$ is an E-function.

*Reasons, in increasing order of strength.*

1. *Elementary.*  E-function-ness of $f=\sum a_nt^n/n!$ needs $a_n\in\overline{\mathbb Q}$.
   The recurrence solution space is spanned by the two integer sequences $p_n,q_n$, whose
   EGFs $\sum p_nt^n/n!$ have radius of convergence exactly 1 (since $p_n/n!\to K\ne0$) —
   they are entire *nowhere*, so not E-functions; the entire solutions are the minimal
   ones, whose coefficients $y_n$ are (conjecturally) transcendental.  Neither
   $\sum p_nz^n$ (divergent, Gevrey-1) nor $\sum p_nz^n/(n!)^2$ (denominators $\sim n!$,
   far too large) is a $G$- or E-function either.
2. *Structural (the real argument).*  By André's theory an E-operator is the
   Fourier–Laplace transform of a $G$-operator, in particular of a **Fuchsian** operator
   with rational exponents.  By §4.2 the Fourier–Laplace transform of $\mathcal L_f$ is
   $\sigma^4\partial^2+2\sigma^3\partial+(-\sigma^4+2\sigma^3+1)$, which has **irregular**
   singular points at $\sigma=0$ and $\sigma=\infty$: not Fuchsian, hence not a
   $G$-operator.  Moreover this rank-2 operator is irreducible (Theorem 6: no Liouvillian
   solutions $\Rightarrow$ no order-1 factor), so no proper factor could be a $G$-operator
   either; and $\mathcal L_f$ has no polynomial solutions (the functional equation
   $((1-t)f)''=af(-t)$ compares degrees $d-1$ and $d$).  ∎

**Consequence.**  Siegel–Shidlovsky/Beukers cannot be invoked.  The task's premise
"this is a rank-4 E-operator because it is irregular at $\infty$" is not sufficient:
irregularity at $\infty$ is necessary but not sufficient for E-operator status, and the
Laplace dual here is *doubly* irregular.

### 5.3  What is true, and what would be needed

* $C,D$ are **exponential periods**: e.g. $r_a=\int_\Gamma t\,w\,dt\big/\int_\Gamma w\,dt$
  for the moment kernel of §7, an irregular-connection period; equivalently a ratio of
  connection constants of $\mathcal L_a$ between $x=\tfrac12$ and $x=1$ (Theorem 4).
* Irrationality and irrationality measure exactly 2 are known (Bala's elementary
  simple-continued-fraction theorem; consistent with $|q_nC-p_n|\sim L/(n+1)!$,
  $q_n\sim Kn!$, giving $|C-p_n/q_n|\asymp 1/((n+1)q_n^2)$).
* To prove transcendence one would need either (a) a genuine closed form making the
  numbers ratios of E-function values — ruled out in the classical special functions by
  §4.4, or (b) a transcendence theory for Stokes/connection data of doubly-confluent
  Heun equations, which does not exist, or (c) an ad-hoc Padé/Hermite construction with
  two independent linear forms, which the growth $|q_nC-p_n|\asymp K'/n\cdot q_n^{-1}$
  does not obviously supply.

**Honest verdict.**  The strongest defensible statement is:
> $C=2Q_{+1}(1^-)/P_{+1}(1^-)$ and $D=P_{-1}(1^-)/(2Q_{-1}(1^-))$ are ratios of
> connection constants at the logarithmic singularity $x=1$ of the rank-4 operator
> $x(1-x)y''''+2(1-x)y'''-4y$, whose rank-2 Fourier–Laplace dual
> $W''=(z^{-4}-2z^{-3}-1)W$ is an irreducible doubly-confluent Heun operator with
> differential Galois group $\mathrm{SL}_2(\mathbb C)$.  They are exponential periods; their
> transcendence is open and is not accessible by Siegel–Shidlovsky/Beukers.

---

## 6. Exact identities found

Let $f=f_{+1}$, $g=f_{-1}$ be the minimal EGFs normalised by $f(0)=g(0)=1$.

**(i)**  From $((1-t)f)''=af(-t)$ at $t=1$:  $\;f_a(-1)=-2a\,f_a'(1)$.
Verified to 121 digits (`hunt.py`).

**(ii)  Bilinear invariant coupling $a=+1$ and $a=-1$.**  With $F_j=(1-t)f_j$,
$G_j=(1+t)f_j(-t)$ and $a_1=-a_2$, the quantity
$F_1'G_2-F_1G_2'+F_2'G_1-F_2G_1'$ is constant in $t$.  Evaluating at $t=0$ and $t=1$:
$$\boxed{\;r_{+1}+r_{-1}-2=-\bigl(f(1)g(-1)+g(1)f(-1)\bigr)\;}$$
Both sides $=-2.167124524983286123608126192049508650459\ldots$; agreement to
$3.9\cdot10^{-121}$.

Numerical values (`hunt.py`, 40 digits):
```
f(1)  = 0.3273115325873843789497158552322748761624    f(-1) = 1.492471837469414745417464072738514194238
g(1)  = 1.340984721872550693395159268707057054761     g(-1) = 0.5063756593040270771118689815449755363031
```

---

## 7. The moment kernel (task A, first item): apparent singularities confirmed

Writing $y_n=\int_\Gamma t^n w(t)\,dt$ the recurrence holds iff
$$w'(t)+w(t)=a\,w(-t)/t^2 .$$
With $w=A+B$ ($A$ even, $B$ odd): $B'=A(a t^{-2}-1)$, $A'=-B(1+a t^{-2})$, whence
$$A''+\frac{2}{t(1+at^2)}A'-\Bigl(1-\frac1{t^4}\Bigr)A=0,$$
normal form $S''=\bigl(t^{-4}-1-3a(1+at^2)^{-2}\bigr)S$ with $A=\sqrt{1+at^2}\,S/t$.

**Proposition 8.**  For $a=\pm1$ the regular singular points $t^2=-1/a$ (i.e. $t=\pm i$ for
$a=+1$, $t=\pm1$ for $a=-1$) are **apparent**.

*Proof.*  At such $t_0$, $\operatorname{Res}_{t_0}\bigl(2/(t(1+at^2))\bigr)=-1$, so the exponents
are $0$ and $2$.  The $u^{-1}$ relation forces $c_1=0$, and the only obstruction, at
$u^0$, is $Q(t_0)c_0$ where $Q=-(1-t^{-4})$.  But $t_0^4=a^{-2}=1$, so $Q(t_0)=0$. ∎

Verified symbolically for all four points and by numerical monodromy: transporting the
fundamental solution once around $t_0=1$ ($a=-1$, radius $0.3$) returns the identity to
$1.2\cdot10^{-27}$ (`kernel2.py`).

**But this does not give a closed form.**  The task's hope was that apparentness makes
$\Gamma$ a circle and $y_n$ the Laurent coefficients of an explicit $\omega$.  It does not:
$t=0$ and $t=\infty$ are *both* irregular ($e^{\pm1/t}$ and $e^{\pm t}$), so $w$ has an
essential singularity at $0$ and grows at $\infty$ — there is no annulus of convergence,
and the two-sided extension of $y_n$ grows factorially as $n\to-\infty$.  $\Gamma$ must be
a Hankel-type contour, and the kernel equation is again a doubly-confluent Heun equation
(now carrying two removable apparent singularities), consistent with §4.

Attempts at explicit solutions all failed: with $A=e^{\beta/t+\gamma t}R$, $R=1$ forces
$\beta^2=1$, $\gamma^2=-a^2$ and then an inconsistent $z^2$ coefficient; the diagonalised
system $U'=(a/t^2)e^{2t}V$, $V'=-(a/t^2)e^{-2t}U$ has genuinely resonant exponentials.

---

## 8. Negative results (what was searched and failed)

All at $\ge40$ digits with mpmath PSLQ, `maxcoeff` $\le10^6$–$10^8$, bases of size $\le8$
(`hunt.py`, `identify1.py`).  Targets: $C$, $D$, $r_{\pm1}$, $L_{\pm1}$, $f(\pm1)$,
$g(\pm1)$, and products/quotients thereof.  Bases tried:

| basis | result |
|---|---|
| Kelvin $\{\mathrm{ber},\mathrm{bei},\mathrm{ber}_1,\mathrm{bei}_1,\mathrm{ker},\mathrm{kei}\}$ at $2$ | none |
| $\{I_0,I_1,J_0,J_1\}(2)$, $\{I_0,I_1,J_0,J_1\}(2\sqrt2)$ | none |
| $\{\Gamma(1/4),\Gamma(3/4),\pi,\sqrt2,e,1\}$ | none |
| `mpmath.identify` with $\{\pi,e,\sqrt2,\log2\}$ | none |
| algebraicity: $\{1,X,\dots,X^8\}$, coeffs $<10^8$ | none |
| internal: $L_a$ vs $\{1,r_\pm,f(\pm1),g(\pm1)\}$ at 44 digits | only 5–6-digit "relations" = noise |

Also failed (first pass, re-confirmed): $C,D$ in the $\mathbb Q$- or $\mathbb Q(\sqrt2)$-span
of $\mathrm{Re}/\mathrm{Im}\,I_{0,1}(2e^{i\pi/4})$; the product ansatz
$f_{\min}(t)=I_0(2\sqrt{a(1+t)})\,I_0(2\sqrt{b(1-t)})$.

Two near-misses worth recording as *coincidences*, not identities:
$r_{+1}\,r_{-1}=-0.2499942\ldots$ (vs $-1/4$, differs at $10^{-6}$) and
$L_{+1}L_{-1}=0.9963864\ldots$ (vs $1$).

The "Bessel-argument deformation" suggested in the task collapses: the family
$y_n=ny_{n-1}+\mu(-1)^ny_{n-2}$ gives $W''=(z^{-4}-2z^{-3}-\mu^2)W$, i.e. only $\mu^2$
enters and it cannot be scaled away; $\mu=0$ is the (trivially solvable) degenerate point
and $\mu^2=1$ is our case.  There is no extra hypergeometric structure to expose.

---

## 9. Remaining gaps and what I would try next

1. **Identify the Stokes data.**  $L_a$ is the Stokes-type constant in
   $Y(z)\approx L\,\mathrm{Re}\bigl[(1+i)(e^{iz}-1)/(iz)\bigr]$ as $z\to\infty$.
   Determining $L_a$ and $r_a$ is a connection problem for
   $W''=(z^{-4}-2z^{-3}-\mu^2)W$ between $z=0$ and $z=\infty$.  If anyone has computed
   central connection coefficients for this specific DCHE family (it is a
   *Coulomb-deformed* rather than Mathieu-deformed member: $a_{-3}\ne0$, $a_{-2}=a_{-1}=0$),
   that is where a closed form would live.  Worth a literature search on
   "doubly-confluent Heun connection problem" / Slavyanov–Lay.
2. **The $\psi''+2\xi(1+\xi^2)^{-1}\psi=0$ Borel picture** (§4.5) is the cleanest form I
   found: 2 regular singular points with $\log$'s and one rank-$1/2$ irregular point.
   Its monodromy/Stokes data *is* $r_a$.  Three singular points, one confluent — this is
   the smallest object in the whole problem and the best target for an exact analysis.
3. **Irrationality measure / linear forms.**  Bala already gives measure 2, which is
   optimal; nothing to improve there.
4. **Report the OEIS erratum** in A244279/A244280 (§1).
5. Not attempted for lack of budget: explicit computation of the $2\times2$ Stokes
   matrices of the DCHE, and a check whether $r_{+1},r_{-1}$ satisfy an algebraic relation
   over $\mathbb Q(L_{+1},L_{-1},f(\pm1),g(\pm1))$ beyond §6(ii).

---

## 10. Files

| file | contents |
|---|---|
| `basics.py` | $C,D$ to 140 digits; minimal solutions; sign patterns |
| `asympt.py` | $r_a$ to 220 digits; first pass at $L_a$ |
| `Lconst.py`, `Lconst2.py` | $L_a$ to 44 digits (parity-separated Richardson) |
| `ogf_verify.py`, `egf_ode.py`, `egf_ode2.py` | symbolic verification of every operator identity |
| `ogf_ode.py` | sympy elimination (kept; note its $a^4\to1$ reduction is incomplete) |
| `local_exponents.py` | indicial equations and the $x=1$ log obstruction |
| `taylor_step.py`, `taylor_step2.py` | task C: Taylor-stepping continuation, 19 and 53 digits |
| `structure.py` | Kovacic analysis + brute-force rational Riccati search |
| `fl_check.py` | Fourier–Laplace transform of $\mathcal L_f$ and its normal form |
| `kernel.py`, `kernel2.py` | moment kernel; apparent-singularity proof + numerical monodromy |
| `hunt.py` | $f(\pm1)$, bilinear invariant, PSLQ hunts |
| `identify1.py` | `mpmath.identify` attempts |
| `summary.py` | regenerates all headline numbers |
| `minimal.py`, `kelvin_test.py`, `product_ansatz.py` | first-pass scripts (kept) |
| `bala.pdf`, `bala.txt` | Peter Bala's note |

---

## 11. One-line answers to the four tasks

* **A.** No closed form found, and §4.4 + §5.2 give principled reasons why none exists in
  classical special functions: the governing operator is an irreducible doubly-confluent
  Heun operator with $\mathrm{SL}_2$ Galois group, not Möbius-equivalent to Bessel/Whittaker,
  and not an E-operator.  The moment kernel's extra singularities are apparent (§7) but
  that does not produce a Laurent expansion.
* **B.** Beukers' theorem is stated precisely in §5.1; it does **not** apply, because the
  relevant functions are not E-functions (Proposition 7).  Transcendence of $C,D$ remains
  open; §5.3 says what would be needed.
* **C.** GPT-5.6's formula is **correct**, verified to 53 digits and **proved** (Theorem 4).
  It is the entirety condition $F(1^-)=0$, i.e. a connection-coefficient identity, not a
  closed form.
* **D.** This file.
