# The weight drop of Cooper's rows: the free integration, not the square root

*Claude (Opus 5), 2026-08-22.  Scripts and logs: `lattice/weight_drop/`
(`01_rows.gp`, `02_qseries.gp`, `03_setup.gp`, `04_growth.gp`, `05_period.gp`,
`06_ident.gp`, `07_roots.gp`, `08_slots.gp`, `09_bogner.gp`).  All rational
computations are exact over $\Q$; all floating computations are at $50$–$150$
digits and the achieved agreement is stated per item.
Answers the task "weight drop = square-root phenomenon?" and repairs
`SOURCES_S18_ZUDILIN.md` §0.2, §3.4 and `EULER_CRITERION.md` §4.1 note 8.*

---

## 0. Verdicts

**(V1) The hypothesis is false in its specific form, and true in its slogan.**
The slogan — *"the companion of the order-three row is built with one integration
fewer, and the operative source is a weight-three object, not the meromorphic
weight-four $\Phi$"* — is **exactly right and is now a theorem** (§2, §3).  The
identification of that weight-three object with the square-root row's source
$\Psi_{\rm root}=g^3u$ is **wrong** (§1).  The correct object is

$$\boxed{\ \Xi:=\thq^{-1}\Phi=\int_0^{t}F\,dt=\sum_{n\ge0}\frac{A_n}{n+1}t^{n+1}
\ \in\ \Z[[q]]\ }$$

— Cooper's free integration itself, promoted from a divisibility statement
($(n+1)\mid A_n$, Bogner 2013) to *the source of the row*.

**(V2) The weight drop is exact, not numerical.**  For each of $s_7,s_{10},s_{18}$:
$$\boxed{\ \Theta=\thq^{-3}\Phi=\thq^{-2}\Xi,\qquad B_n=[t^n]\bigl(F\,\thq^{-2}\Xi\bigr),\qquad
\xi_\infty=L(\Xi,2)=4\pi^2\!\int_0^\infty\!\Xi(iy)\,y\,dy,\qquad L(\Xi,1)=0 .}$$
So the companion is a **double** Eichler integral of a source occupying the
weight-three slot — the analytic shape of a **weight-one** row — even though
$F$ has weight $2$ and $\Phi$ weight $4$.  The archimedean period is therefore an
$s=2$ critical slot, not an $s=3$ one.  \[Verified to $10^{-50}$, §3.3.\]

**(V3) The three weight-three-slot data are identified, exactly.**  With the
Theorem B$^*$ shape $\xi_\infty=P(2)L(\psi,2)L(\varphi,0)$ for $w=1$:

| row | $N$ | $\psi$ | $\varphi$ | $P(2)$ | $\xi_\infty$ | agreement |
|---|---|---|---|---|---|---|
| $s_7$    | $7$  | $\mathbf 1$ | odd, $L(\varphi,0)=1$ (fits $\chi_{-7}$) | $1/7$ | $\zeta(2)/7$ | $1/7$ to $55$ digits |
| $s_{10}$ | $10$ | $\mathbf 1$ | odd | $P(2)L(\varphi,0)=1/5$ | $\zeta(2)/5$ | $1/5$ to $55$ digits |
| $s_{18}$ | $18$ | $\chi_{-3}$ | $\mathbf 1$ | $\mathbf{-1}$ | $\tfrac12L(2,\chi_{-3})$ | $-1$ to $55$ digits |

**$s_{18}$'s census entry $(\psi,\varphi)=(\chi_{-3},\mathbf 1)$, $Q(2)=-1$ was
right all along; what was wrong was the weight.**  The row is $w=1$ *in the slot
that matters*, not $w=2$.

**(V4) Consequently Theorem F applies to $s_{18}$, with $w=1$, and gives the
measured $3$-adic value.**  $\psi=\chi_{-3}$, $p=3\mid\operatorname{cond}\psi$, so
criterion (a) is vacuous, $Q=P$, $Q(2)=-1$, and
$$\kappa_3=\tfrac12L_3\bigl(2,\psi\omega^{-1}\bigr)=\tfrac12L_3(2,\mathbf 1)=\tfrac12\zeta_3(2),
\qquad \xi_3=-Q(2)\kappa_3=\tfrac12\zeta_3(2),$$
matching the census measurement to $3^{3025}$ and `SOURCES_S18_ZUDILIN.md` §4 to
$3^{219}$.  Conjecture D holds with $r_3=r_\infty=\tfrac12$ (because
$\mathcal E_3(2)=1$), which is exactly why `SLOPE_CENSUS.md` §1 finds $s_{18}$
$3$-adically **aligned** with Zagier $\mathbf B$ and $\mathbf C$: it carries the
*same* $(\psi,\varphi,P)$ data as $\mathbf C$.  **`SOURCES_S18_ZUDILIN.md` §0
verdict 2 ("structurally wrong … Open") is withdrawn** (§5).

**(V5) For $s_7,s_{10}$ the same theory predicts $\xi_p=0$ whenever a $p$-adic
limit exists** ($\psi=\mathbf 1$ and $w=1$ force $\varphi$ odd, hence
$\delta(\varphi)=0$, hence $\kappa_p=0$).  `SLOPE_CENSUS.md` §1 finds **no**
$p$-adic limit at $p=2,3,5,7$ for either row, so the prediction is not
contradicted; it is also not tested.  Note that this *repairs* an apparent
counterexample to Proposition C: $s_7$ has $3\mid c$ and $s_{10}$ has $2\mid c$
yet no slope — under the $w=1$ reading this is the ordinary "Euler factor does not
divide $P$" failure of criterion (a), not a breakdown.

**(V6) What the square root really contributes.**  The $\Sym^2$ half of the
hypothesis is *true and exact*: $\lambda^nA_n=\sum_ia_ia_{n-i}$ for all three rows
($n\le400$, exact).  The companion half is *false*: $\lambda^nB_n-(a*b)_n\neq0$,
and
$$\frac{\xi_\infty^{\rm par}}{\xi_\infty^{\rm root}}=0.502963744983\ldots,\quad
1.038057712356\ldots,\quad 0.806734635801\ldots$$
admit no algebraic relation of degree $\le6$ and no two- or three-term $\Q$-linear
relation at $120$ digits (§1.2).  Moreover the root rows of $s_{10}$ and $s_{18}$
are **non-congruence**: $v_2\bigl(\operatorname{den}[q^m]\sqrt F\bigr)=m-s_2(m)$
exactly for $m\le60$ (unbounded), so by Calegari–Dimitrov–Tang neither $g=\sqrt F$
nor $\Psi_{\rm root}=g^3u$ lies in *any* congruence space $M_k(\Gamma_0(M),\chi)$,
and `mftobasis` cannot be run at all.  The task's conjecture "$s_{18}$'s root is a
congruence CM theta" is **refuted** (§1.3).  Only $s_7$'s root is congruence, and
there it is beautiful and exact:
$$\sqrt{F_7}=E_1(\chi_{-7}),\qquad \Psi^{s_7}_{\rm root}=\bigl(\sqrt{F_7}\bigr)^3t_7=\eta(\tau)^3\eta(7\tau)^3$$
(both to $q^{90}$) — but its period $L(\eta_1^3\eta_7^3,2)=0.46721176888\ldots$ is
**not** $\zeta(2)/7$, nor twice, nor half it.

**(V7) One genuinely new exact evaluation.**  For all three rows
$\xi_\infty=-\tfrac{2\pi^2}{N}L(\Phi,1)$, i.e. $L(\Phi,1)=-\Xi(\text{cusp }0)$ is
the *conifold period* $\int_{\rm path}F\,dt$, and it is **rational** for the two
$\zeta(2)$ rows:
$$L(\Phi_{s_7},1)=-\tfrac1{12},\qquad L(\Phi_{s_{10}},1)=-\tfrac16
\qquad(\text{to }60\text{ digits}).$$

**(V8) Open, and sharply posed.**  The one step not proved is *why*
$L(\Xi,2)$ equals the Eisenstein critical value: $\Xi$ is meromorphic and its
higher slots are **not** Eisenstein ($P(3),P(4)$ come out irrational, §4.2), so
the statement is a **critical-slot purification**: the principal part of $\Xi$ is
invisible exactly at $s=2$, and nowhere else.  §7 states the clause.

---

## 1. The square root is not the explanation

### 1.1 What is true: the $\Sym^2$ Cauchy square
With $a_n=\lambda^n[t^n]\sqrt{\sum A_nt^n}$ (the $\Sym^1$ rows of
`SPORADIC_SCAN2.md` §8, $\lambda=1,2,2$),
$$\lambda^nA_n=\sum_{i=0}^na_ia_{n-i}\qquad\text{exactly, }n\le400,\ \text{all three rows}$$
(`01_rows.gp`; zero error).  This is `SQRT_APERY.md` Theorem 2/§7 restated at the
level of coefficients.

### 1.2 What is false: any bilinear formula for $B$
With $b$ the root companion ($b_0=0,b_1=1$, same recurrence) and
$\delta_n:=\lambda^nB_n-(a*b)_n$:

| row | $\delta_1,\dots,\delta_8$ |
|---|---|
| $s_7$ | $0,\,0,\,-1,\,-\tfrac{188}9,\,-\tfrac{15455}{36},\,-\tfrac{680372}{75},\,-\tfrac{4934776}{25},\,-\tfrac{5393662648}{1225}$ |
| $s_{10}$ | $0,\,1,\,\tfrac{19}2,\,\tfrac{1796}9,\,\tfrac{149663}{36},\,\tfrac{7332451}{75},\,\tfrac{180362119}{75},\,\tfrac{32357835424}{525}$ |
| $s_{18}$ | $0,\,1,\,\tfrac{45}2,\,\tfrac{1508}3,\,\tfrac{46501}4,\,\tfrac{6951633}{25},\,\tfrac{171209311}{25},\,\tfrac{211890952512}{1225}$ |

$\delta$ is of the same exponential order as $A$, so it *does* move the limit, and
$$\lim\frac{(a*b)_n}{(a*a)_n}=\tfrac12\xi_{\rm root}+\tfrac12\frac{g_B}{g_0}(t_c)\neq\xi_{\rm root},$$
($g_\bullet$ the analytic parts at the conifold; the naive "$\lim=\xi_{\rm root}$"
of the hypothesis ignores the analytic part of the local expansion).
At $n=900$ and $120$ digits (`06_ident.gp`):
$\xi^{\rm par}-\text{target}=3.4\cdot10^{-136},\,-6.9\cdot10^{-136},\,-1.4\cdot10^{-113}$,
$\xi^{\rm root}=0.4672117688843731239\ldots,\;0.3169253592103549809\ldots,\;
0.4842375536043764245\ldots$; `algdep` of the ratio at degrees $4,6$ returns
minimal heights $10^{15}$–$10^{21}$, `lindep` on $(1,\xi^{\rm par},\xi^{\rm root})$
returns coefficients $\approx10^{35}$ and on $(\xi^{\rm par},\xi^{\rm root})$
$\approx10^{53}$.  **No relation.**

### 1.3 The root rows of $s_{10}$, $s_{18}$ are non-congruence
`07_roots.gp`, $m\le60$:
$$v_2\bigl(\operatorname{den}[q^m]\sqrt{F}\bigr)=m-s_2(m)\quad (s_{10}\text{ and }s_{18},\ \text{identical vectors}),
\qquad =0\quad (s_7),$$
and for $\Psi_{\rm root}=g^3t/\lambda$ the denominators are likewise unbounded.
Unbounded denominators $\Rightarrow$ non-congruence (Calegari–Dimitrov–Tang),
so **there is no $(M,\chi)$ for which `mftobasis` could place $\Psi_{\rm root}$**,
and the hypothesis's step (1) cannot even be formulated for $s_{18}$.
(Consistent with `SQRT_APERY.md` §7's dichotomy: $\lambda=1\iff e\ge2\iff$
integral CM theta $\iff$ period identified.  $s_{10},s_{18}$ have $\lambda=2$.)

For $s_7$ ($\lambda=1$) the hypothesis's step (1) *does* hold and is now exact:
$$\sqrt{F_7}=E_1(\chi_{-7})=1+2\sum_{n\ge1}\Bigl(\sum_{d\mid n}\chi_{-7}(d)\Bigr)q^n
\in M_1(\Gamma_0(7),\chi_{-7}),\qquad
\Psi_{\rm root}=\eta_1^3\eta_7^3\in S_3(\Gamma_0(7),\chi_{-7}),$$
both verified as exact $q$-series identities to $q^{90}$.  So $s_7$'s root row is a
*congruence weight-one row with a CM cuspidal weight-three source*, period
$L(\eta_1^3\eta_7^3,2)$ — a perfectly good Theorem B$^*$ row, and a **different**
one from its parent.  Two rows on the same curve, two sources, two periods.

---

## 2. The free integration is the source (exact)

### 2.1 Statement

> **Theorem WD1.**  Let $(\Gamma_0(N),t,F)$ be an order-three row, $w=2$,
> $\Phi=F\thq t$, $A_n=[t^n]F$, $B_n=[t^n](F\thq^{-3}\Phi)$.  Put
> $$\Xi:=\thq^{-1}\Phi=\sum_{m\ge1}\frac{c(m)}m q^m .$$
> Then, as $t$-series, $\Xi=\sum_{n\ge0}\frac{A_n}{n+1}t^{n+1}=\int_0^tF\,dt$; and
> $$\Theta=\thq^{-3}\Phi=\thq^{-2}\Xi,\qquad B_n=[t^n]\bigl(F\,\thq^{-2}\Xi\bigr).$$
> For Cooper's $s_7,s_{10},s_{18}$, $\Xi\in\Z[[q]]$ — equivalently
> $(n+1)\mid A_n$ (Bogner 2013).

*Proof.*  $\thq\Xi=\Phi=F\thq t$ gives $d\Xi/dt=F$ and $\Xi(0)=0$, i.e.
$\Xi=\int_0^tF\,dt$; the two displayed identities are then the definition of
$\thq^{-1}$.  $\square$

\[**Verified** (`02_qseries.gp`): $\Xi\in\Z[[q]]$ to $q^{58}$ for all three rows;
$\thq^{-3}\Phi-\thq^{-2}\Xi=0$; $B_n=[t^n](F\thq^{-2}\Xi)$ exactly for $n\le22$;
$(n+1)\mid A_n$ and $A_n\in\Z$ for $n\le600$ (`09_bogner.gp`).\]

### 2.2 Why this is the $\Sym^2$ statement the hypothesis was reaching for

The hypothesis asked for "the inhomogeneous-equation structure: $\Sym^2$ of the
root system should relate the two companions".  It does, and the free integration
is *precisely what comes out*.  Write $L_1=u^2P\,\widetilde L_1$,
$\widetilde L_1=D^2+pD+\varrho$ ($D=d/du$, $u=t/\lambda$) for the root operator,
$A=\sqrt F$ its analytic solution, $W=e^{-\int p}=A^2\,d\log q/du$ its Wronskian,
and $S=\Sym^2\widetilde L_1$, so that $L_{\rm par}=\ell\cdot S$ (`SQRT_APERY.md`
§7: the $\Sym^2$ identity is exact for all nine rows, Cooper's included).

**Lemma WD2.**  For any $z$, $\;S(Az)=A\,(\widetilde L_1z)'+(3A'+2pA)(\widetilde L_1z)$.

*(Direct expansion, using $A''=-pA'-\varrho A$; the $z''$, $z'$, $z^{(0)}$
coefficients match term by term.)*

Now $L_{\rm par}(B_{\rm par})=\kappa\,t$ (the inhomogeneity created by
$B_1=1\neq0$), so $S(B_{\rm par})=\kappa\lambda u/\ell$.  Setting $z=B_{\rm par}/A$
and $h=\widetilde L_1 z$, Lemma WD2 is a **first-order** equation for $h$,
$$A h'+(3A'+2pA)h=\kappa\lambda u/\ell
\;\Longrightarrow\;
\frac{d}{du}\Bigl(\frac{A^3h}{W^2}\Bigr)=\frac{A^2}{W^2}\cdot\frac{\kappa\lambda u}{\ell},$$
and with $\ell=u^3P$, $W^2=1/(u^2P)$ the right-hand side collapses to
$\kappa\lambda A^2/u\cdot u = \kappa\lambda A^2$, whence
$$\boxed{\ \widetilde L_1\Bigl(\frac{B_{\rm par}}{A}\Bigr)
=\frac{\kappa\lambda\,W^2}{A^3}\int_0^uA^2\,du
=\frac{\kappa\lambda\,W^2}{A^3}\,\Xi\ }$$
because $A^2=F$.  **The free integration $\Xi=\int F\,du$ is exactly the
inhomogeneity of the root system that produces the parent's companion.**
Variation of parameters in the basis $\{A,A\log q\}$ (Wronskian $W=A^2d\log q/du$)
then gives $B_{\rm par}=\kappa\lambda\,A^2\thq^{-2}\Xi=\kappa\lambda F\thq^{-2}\Xi$,
recovering Theorem WD1 with $\kappa\lambda=1$ in the census normalisation.
Compare the root's own companion (`SQRT_APERY.md` Thm 3):
$$b_n=\lambda^n[t^n]\bigl(g\,\thq^{-2}\Psi_{\rm root}\bigr),\qquad
B_n=[t^n]\bigl(g^2\,\thq^{-2}\Xi\bigr),\qquad g=\sqrt F .$$
*Same double Eichler integral, different sources* — which is exactly why the two
periods are unrelated.

---

## 3. The archimedean weight drop

### 3.1 The geometry: all three Cooper rows are Fricke-fold rows
For $F=\tfrac1{c}\sum_{d\mid N}a_d\,dE_2(d\tau)$ one has
$dE_2(d\tau)\big|_2W_N=(N/d)E_2((N/d)\tau)$, i.e. $a_d\mapsto a_{N/d}$.  Cooper's
three all satisfy $a_{N/d}=-a_d$:
$$(a_1,a_7)=(-1,1);\quad (a_1,a_2,a_5,a_{10})=(-1,-1,1,1);\quad
(a_1,a_2,a_3,a_6,a_9,a_{18})=(-1,1,2,-2,-1,1),$$
hence $F|_2W_N=-F$.  Also $x_7|W_7=49/x_7$, $x_{10}|W_{10}=x_{10}$,
$x_{18}|W_{18}=x_{18}$, and in each case $t\circ W_N=t$; therefore
$\thq t|_2W_N=+\thq t$ and
$$\boxed{\ t\circ W_N=t,\qquad F|_2W_N=-F,\qquad \Phi|_4W_N=-\Phi\ }$$
— **exactly the hypotheses of Theorem II of `THEOREM_B_EXACT.md`**.  And the fold
sits where Theorem II says it does (`05_period.gp`, $50$ digits):
$$t_7(i/\sqrt7)=\tfrac1{27},\qquad t_{10}(i/\sqrt{10})=\tfrac1{16},\qquad
t_{18}(i/\sqrt{18})=\tfrac1{16}$$
— the dominant singularity, to $10^{-55}$ in each case.

The **only** failure of Theorem II's hypotheses is that $\Phi$ is meromorphic.
Its poles lie strictly *inside* the Fricke circle: measuring
$|c(m)|^{1/m}\to1/|q_0|$ (`04_growth.gp`, $m\le200$) gives
$$1/|q_0|=2.2586\ (s_7),\quad 1.1776\ (s_{10}),\quad 2.9282\ (s_{18}),$$
against $e^{2\pi/\sqrt N}=10.749,\;7.293,\;4.397$; so
$\operatorname{Im}\tau_{\rm pole}<1/\sqrt N$ in all three cases and the whole
geodesic $(0,i\infty)$ — the contour of the fold argument — avoids the poles.
(For $s_{18}$, $\operatorname{Im}\tau_0=1/6<1/\sqrt{18}$ recovers
`SOURCES_S18_ZUDILIN.md` §3.1's $\tau_0=(3+i)/6$.)

> **Theorem WD3 (Theorem II for meromorphic Fricke sources).**  Let $r=3$,
> $t\circ W_N=t$, $F|_2W_N=-F$, and let $\Phi=F\thq t$ be meromorphic with all
> poles off the geodesic $(0,i\infty)$ and of imaginary part $<1/\sqrt N$.  Then
> $\Lambda(\Phi,s):=\int_0^\infty\Phi(iy)y^{s-1}dy$ is entire, satisfies
> $\Lambda(\Phi,s)=-N^{2-s}\Lambda(\Phi,4-s)$, and
> $$\xi_\infty=4\pi^3\Lambda(\Phi,3),\qquad \Lambda(\Phi,2)=0,\qquad
> \Lambda(\Phi,1)=-N\Lambda(\Phi,3).$$

*Proof sketch.*  Theorem II's proof uses only (i) the Eichler representation
$\Theta(\tau)=\frac{(-2\pi i)^3}{2}\int_\tau^{i\infty}\Phi(z)(z-\tau)^2dz$, valid
by termwise integration for $\operatorname{Im}\tau>\operatorname{Im}\tau_{\rm pole}$
and by continuation along the (pole-free) imaginary axis elsewhere; (ii) the
substitution $z=W_Nu$, which maps the ray $[\tau_*,i\infty)$ onto the segment
$(0,\tau_*]$ — no contour deformation, hence no residues; (iii) the fold lemma at
$\tau_*=i/\sqrt N$.  All three survive verbatim.  $\square$
\[Class: **proved modulo the transposition**; **verified** to $10^{-50}$, §3.3.\]

### 3.2 The drop, in one integration by parts
$\Phi=\thq\Xi$, so on the imaginary axis $\Phi(iy)=-\frac1{2\pi}\frac{d}{dy}\Xi(iy)$
and, the boundary terms vanishing ($\Xi=O(q)$ at $\infty$; $\Xi(iy)\to\Xi(0)$
finite),
$$\Lambda(\Phi,s)=\frac{s-1}{2\pi}\Lambda(\Xi,s-1),\qquad
\Lambda(\Xi,s)=\int_0^\infty\Xi(iy)y^{s-1}dy .$$
Writing $L(f,s)=(2\pi)^s\Lambda(f,s)/\Gamma(s)$, Theorem WD3 becomes

$$\boxed{\ \xi_\infty=L(\Xi,2),\qquad L(\Xi,1)=0,\qquad
\Xi(\text{cusp }0)=-L(\Phi,1)=\frac{N}{2\pi^2}\,\xi_\infty\ }$$

**This is the weight drop.**  The row's period is the $s=2$ slot of a
weight-three-slot source, i.e. the slot of a *weight-one* row, and the
"endpoint" identity $L(\Xi,1)=0$ is the $w=1$ analogue of Theorem II's forced
$L(\Phi,2)=0$.

### 3.3 Numerical verification (`05_period.gp`, `08_slots.gp`)
$\Lambda(\Phi,s)=I(s)-N^{2-s}I(4-s)$ with
$I(s)=\sum_mc(m)(2\pi m)^{-s}\Gamma(s,2\pi m/\sqrt N)$, $q$-precision $320$–$420$,
$50$–$60$ digits:

| row | $\Lambda(\Phi,2)$ | $4\pi^3\Lambda(\Phi,3)$ vs. recurrence limit | $L(\Phi,1)$ |
|---|---|---|---|
| $s_7$ | $0$ ($<10^{-59}$) | agree to $\mathbf{10^{-59}}$ ($=\zeta(2)/7$) | $-\tfrac1{12}$ (60 digits) |
| $s_{10}$ | $0$ ($<10^{-59}$) | agree to $\mathbf{10^{-59}}$ ($=\zeta(2)/5$) | $-\tfrac16$ (60 digits) |
| $s_{18}$ | $0$ ($<10^{-59}$) | agree to $\mathbf{10^{-50}}$ | $-0.3562311836577892449\ldots$ |

(The $s_{18}$ line is compared against the $150$-digit value produced by the
recurrence in `01_rows.gp`, $0.3906512064482431484335937148120461781825671682726\ldots$,
not against a truncated series for $L(2,\chi_{-3})$.)

---

## 4. The weight-three-slot data

### 4.1 $s=2$: exact
`08_slots.gp`, $55$ digits:
$$\frac{L(\Xi_{s_7},2)}{\zeta(2)}=\frac17,\qquad
\frac{L(\Xi_{s_{10}},2)}{\zeta(2)}=\frac15,\qquad
\frac{L(\Xi_{s_{18}},2)}{L(\chi_{-3},2)\,\zeta(0)}=-1 .$$
Each is the Theorem B$^*$ ($w=1$) shape $\xi_\infty=P(2)L(\psi,2)L(\varphi,0)$:
$\psi=\mathbf 1$ forces $\varphi$ odd (parity: $\psi\varphi(-1)=(-1)^w=-1$) and
then $L(\varphi,0)\in\Q$; $\psi=\chi_{-3},\varphi=\mathbf 1$ gives
$L(\varphi,0)=\zeta(0)=-\tfrac12$ and $P(2)=-1$.  Since only $\psi=\mathbf 1$
produces $L(\psi,2)\in\pi^2\Q$, the characters are forced: $s_7,s_{10}$ have
$\psi=\mathbf 1$ and $s_{18}$ has $\psi=\chi_{-3}$.  For $s_7$, taking
$\varphi=\chi_{-7}$ (the only odd character of conductor dividing $7$;
$L(\chi_{-7},0)=1$) gives $P(2)=1/7$ **exactly**.

### 4.2 $s=3,4$: not Eisenstein — the purification is critical-slot only
$$L(\Xi,3)=0.4916851525975040346\ldots,\;0.6116606136143467045\ldots,\;
0.6449077699720922369\ldots$$
$$L(\Xi,4)=0.6932931505629112763\ldots,\;0.7935476058274871165\ldots,\;
0.7968813054872598556\ldots$$
Dividing out the Eisenstein shape gives $P(3)=0.34447778\ldots$,
$P(4)=0.55607778\ldots$ for $s_7$ — `bestappr` at $10^{12}$ returns nothing
rational, and for $s_{18}$, $L(\Xi,3)/L(\chi_{-3},3)=0.72951402598\ldots$ has no
relation to $1,\log2,\log3$ (`lindep` coefficients $\approx10^{11}$).
So $L(\Xi,s)$ is **not** $P(s)L(\psi,s)L(\varphi,s-2)$ as a function of $s$; the
Eisenstein identity holds at $s=2$ and there only.  This is the precise sense in
which the meromorphic principal part is "controlled": it is annihilated by the
critical slot, not by the whole $L$-function.

### 4.3 A structural remark on the residues
At a pole $\tau_0$ that is a fixed point of an involution $\sigma$ with
$t\circ\sigma=t$ and $F|_2\sigma=-F$ (this is $s_{18}$'s $W_9$-fixed pole
$\tau_0=(3+i)/6$), $t$ is even and $F$ is even in a local coordinate in which
$\sigma$ is $u\mapsto-u$, so $\Phi\,d\tau=F\,dt/2\pi i$ is **odd**, and its
residue vanishes.  Hence $\Phi\,d\tau$ is exact and $\Xi$ single-valued near
$\tau_0$: Cooper's notes §12 hoped the principal part is "a total derivative",
and at fold-type poles it provably is.  (For $s_7$ the two poles are exchanged,
not fixed, by $W_7$; the residue question there is untested and is not needed —
the fold argument only ever integrates along the pole-free imaginary axis.)

---

## 5. The $p$-adic side: Theorem F, at $w=1$

Theorem F reads off $(\psi,\varphi,P)$ from the source and $w$ from the number of
Eichler integrations.  With Theorem WD1 the number of integrations acting on the
*operative* source is **two**, i.e. $w=1$, and the source data are those of §4.1.
For $s_{18}$ at $p=3$:

* criterion (a): $p=3\mid\operatorname{cond}\chi_{-3}$, so
  $\mathcal E_3(s)=1-\chi_{-3}(3)3^{-s}=1$ and divisibility is vacuous; $Q=P$,
  $Q(2)=-1$;
* $\omega=\chi_{-3}$ at $p=3$, so $\psi\omega^{-w}=\chi_{-3}\omega^{-1}=\mathbf 1$
  and $\kappa_3=\tfrac12L_3(2,\mathbf 1)=\tfrac12\zeta_3(2)$;
* hence $\xi_3=-Q(2)\kappa_3=\tfrac12\zeta_3(2)$ — the measured value
  ($3^{3025}$ in `EULER_CRITERION.md` §4.1, $3^{219}$ in
  `SOURCES_S18_ZUDILIN.md` §4);
* Conjecture D: $r_3=r_\infty/\mathcal E_3(2)=r_\infty=\tfrac12$, so $s_{18}$
  carries **the same $(\psi,\varphi,P)$ data as Zagier $\mathbf C$** — which is
  exactly the $3$-adic alignment `SLOPE_CENSUS.md` §1 measures ($v_3\approx n$ to
  $n=300$ against both $\mathbf B$ and $\mathbf C$).  The alignment is no longer a
  rigidity transfer; it is an identity of source data.

**Corrections to earlier files.**

1. `SOURCES_S18_ZUDILIN.md` §0 verdict 2 ("the census entry is numerically right
   and structurally wrong … **Open**") is **withdrawn**: the entry is right, and
   the structural error was the weight bookkeeping ($w=2$ instead of $w=1$).
   §3.4 ("What is *not* explained: the weight drop") is superseded by §2–§3 here.
   §0 verdict 4 ("the archimedean Eisenstein-source proof is not available")
   stands only for the *weight-four* factorisation $\Phi=W\cdot x(x-3)/(x+3)^3$;
   the operative source is $\Xi=\thq^{-1}\Phi$, one integration down.
2. `EULER_CRITERION.md` §4.1: the $s_{18}$ row's "(conj.)" may be dropped and its
   $w$ recorded as $1$; note 8's parenthetical "the fit $Q(2)=-1$ is an
   observation, not a derivation" is upgraded to a derivation modulo §7(a).
3. `SOURCES_S18_ZUDILIN.md` §7 item 1 (drop "operator order two" for Cooper's
   rows) stands: the *operator* order is three.  What is two is the number of
   Eichler integrations that act on the operative source.

---

## 6. Why the two mechanisms are cousins but not the same

Both the square root and the free integration lower the integration count of the
order-three row from three to two; that is why `SPORADIC_SCAN2.md` §8 finds
$k=2$ for Cooper's rows *and* for their square roots.  They do it differently:

| | square root ($\Sym^2\to\Sym^1$) | free integration ($(n+1)\mid A_n$) |
|---|---|---|
| stays on the same curve | yes | yes |
| new operator | order $2$, $L_1$ with $\Sym^2L_1=L_{\rm par}$ | none; same $L_{\rm par}$ |
| new row | $a_n=\lambda^n[t^n]\sqrt F$ | none; same $A_n$ |
| source | $\Psi_{\rm root}=g^3u$, weight $3$, holomorphic | $\Xi=\thq^{-1}\Phi$, weight-$3$ slot, meromorphic |
| congruence? | $s_7$ yes, $s_{10},s_{18}$ **no** | all three integral (Bogner) |
| period | $L(\Psi_{\rm root},2)$ | $L(\Xi,2)=\xi_\infty$ |
| cost $k$ | $2$ | $2$ |

The row whose limits the census records is the **parent**, so the operative
mechanism is the free integration.  The square root is a genuine second row on the
same curve with its own, different, period.

---

## 7. The general clause, and what remains

> **Clause (weight drop).**  Let $(\Gamma_0(N),t,F)$ be a row of weight $w=2$ with
> $t\circ W_N=t$, $F|_2W_N=-F$, $\Phi=F\thq t$, $\Xi=\thq^{-1}\Phi$.
> **(i)** If $\Phi$ is holomorphic, $\Lambda(\Phi,3)=\Gamma(3)(2\pi)^{-3}L(\Phi,3)$
> with $L(\Phi,s)=\sum c(m)m^{-s}$ convergent, and $\xi_\infty=L(\Phi,3)$ is a
> genuine weight-two ($s=w+1$) critical value: $\zeta(3)$-type.  **No drop.**
> **(ii)** If $\Phi$ is meromorphic with all poles at points of imaginary part
> $<1/\sqrt N$ off the geodesic $(0,i\infty)$, then $\sum c(m)m^{-3}$ diverges,
> $\xi_\infty=4\pi^3\Lambda(\Phi,3)=L(\Xi,2)$, and the period occupies the
> **$s=2$** slot of the weight-three-slot source $\Xi$.  **Drop by one.**
> In case (ii) $\Xi\in\Z[[q]]$ iff $(n+1)\mid A_n$, and then the row costs
> $k=2$ and Theorem B$^*$/Theorem F are to be applied with $w=1$ and the data
> $(\psi,\varphi,P)$ of $\Xi$.

Cooper's $s_7,s_{10},s_{18}$ are case (ii); Apéry $\gamma$, Domb $\alpha$,
$\varepsilon$, $\zeta$ are case (i) (and indeed $(n+1)\nmid A_n$ there: Apéry
$A_1=5$, Domb $A_2=28$).

**Open, in order of value.**

* **(a) The critical-slot purification.**  Prove that in case (ii)
  $L(\Xi,2)=P(2)L(\psi,2)L(\varphi,0)$ for a holomorphic weight-three Eisenstein
  datum $(\psi,\varphi,P)$, i.e. that the meromorphic part of $\Xi$ contributes
  $0$ at $s=2$.  §4.2 shows this *fails* at $s=3,4$, so any proof must be
  slot-specific — presumably: the principal part of $\Xi$ is (locally at each
  CM pole) a $W$-anti-invariant object whose $s=2$ Mellin slot is the one killed
  by the fold, exactly as $L(\Phi,2)=0$ is killed in Theorem II.  This is the
  whole remaining content of "Cooper's limits as theorems".
* **(b) $P$ itself.**  Determine the Mellin polynomials: $P(2)=1/7$ for $s_7$
  (with $\varphi=\chi_{-7}$), $P(2)L(\varphi,0)=1/5$ for $s_{10}$ (the value is
  Zagier row $\mathbf D$'s, suggesting the same quartic $\varphi$ mod $5$ and a
  Conjecture-D pair), $P(2)=-1$ for $s_{18}$ (Zagier row $\mathbf C$'s data).
  A determination of $P(3)$ would decide the endpoint condition
  $\delta(\varphi)P(3)=0$; §4.2 shows $L(\Xi,3)$ finite, consistent with
  $P(3)=0$ when $\varphi=\mathbf 1$ ($s_{18}$).
* **(c) The $p$-adic prediction for $s_7,s_{10}$.**  $\varphi\neq\mathbf 1$ gives
  $\kappa_p=0$, hence $\xi_p=0$ whenever the limit exists; no limit exists at
  $p\le7$.  Testing at a larger prime, or proving criterion (a) fails, would close
  the loop.
* **(d) Meromorphic Theorem F.**  The $p$-adic argument still runs through
  $\Theta=\thq^{-2}\Xi$ with $\Xi$ meromorphic.  `SOURCES_S18_ZUDILIN.md` §3.3's
  supersingular-CM-pole lemma is the right first step and is unaffected; what is
  needed is that the *depleted* $\Xi^{[p]}$ is overconvergent, i.e. the analogue
  of (a) $p$-adically.  The measurement $\xi_3=\tfrac12\zeta_3(2)$ to $3^{3025}$
  says it is true.
* **(e) The conifold periods.**  $\int_{\rm path}F\,dt=\tfrac1{12}$ ($s_7$),
  $\tfrac16$ ($s_{10}$) are exact rationals (60 digits); a direct proof of either
  *is* a proof of $\xi_\infty=\zeta(2)/7$, $\zeta(2)/5$, by
  $\xi_\infty=\tfrac{2\pi^2}N\Xi(0)$.  This may be the shortest route to (a) in
  those two cases, since it asks only for the monodromy period of $F$ at the
  conifold and not for any $L$-function.
* **(f) The three unidentified $\Sym^1$ periods** ($0.14551\ldots$,
  $0.31693\ldots$, $0.48424\ldots$) remain open, and §1.3 explains why the
  $1792$-element congruence sweep had to fail for the last two: the forms are
  non-congruence.  This is the same explanation `SQRT_APERY.md` §6 gives for
  Apéry's own $\Sym^1$ period.

---

## 8. Reproduction

```
lattice/weight_drop/01_rows.gp     # parent/root rows; A = a*a; delta_n; xi's
lattice/weight_drop/02_qseries.gp  # x,t,F,Phi,Xi,g,Psi_root; Theorem WD1 checks
lattice/weight_drop/03_setup.gp    # shared q-series builder (Cooper's eta data)
lattice/weight_drop/04_growth.gp   # |c(m)|^{1/m} vs e^{2 pi/sqrt N}: poles inside the circle
lattice/weight_drop/05_period.gp   # fold t(i/sqrt N)=t_c; Lambda(Phi,s); xi = 4 pi^3 Lambda(Phi,3)
lattice/weight_drop/06_ident.gp    # n=900 limits; algdep/lindep refutations
lattice/weight_drop/07_roots.gp    # sqrt(F) denominators; sqrt(F_7)=E_1(chi_-7); Psi=eta_1^3 eta_7^3
lattice/weight_drop/08_slots.gp    # L(Xi,2),L(Xi,3),L(Xi,4); P(2)=1/7,1/5,-1
lattice/weight_drop/09_bogner.gp   # A_n in Z and (n+1)|A_n to n=600
```
