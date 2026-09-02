# Cooper's magnetic congruence, round 2: the kernel $\gamma$ identified as a twisted CM
# trace, and the master conjecture as a divisibility of traces of singular moduli

*Working note, 2026-09-02.  PARI/GP 2.15.4; all divisibility claims are exact integer
arithmetic, all CM evaluations are floating point with the precision quoted.  Scripts
`01_*`–`59_*` with matching `.log` files in this directory; three companion verdict files
`FINDINGS_MODP.md`, `FINDINGS_XSIDE.md`, `FINDINGS_ROWS.md`; two working notes
`NOTES_IDENT.md`, `NOTES_TRACE.md`.  Background: `consolidation/COOPER_CONGRUENCE.md`,
`lattice/cooper_congruence/REPORT.md` (round 1), `lattice/cooper_sources/REPORT.md` §2,
`lattice/cooper_sources/pz.txt` (Paşol–Zudilin).  Nothing outside
`lattice/cooper_congruence/round2/` was created or modified.*

Claims are tagged **[proved]**, **[verified, range]** (exact computation over the stated
range), **[num, digits]**, **[refuted, range]**, **[conjectural]**.

---

## 0. Verdict table

| # | statement | status |
|---|---|---|
| **R1** | **THE IDENTIFICATION (row $s_7$).**  Let $f=1/(xF)=q^{-1}+5+19q+52q^2+\cdots$, the weakly holomorphic modular form of **weight $-2$** on $\Gamma_0(7)$ with a simple pole at the cusps, $f|_{-2}W_7=-f$; let $\widehat f:=Df+\dfrac{f}{2\pi y}=-\dfrac{R_{-2}f}{4\pi}$ ($R_{-2}$ = Maass raising), a real-analytic $\Gamma_0(7)$-invariant weight-$0$ function.  Then $$\beta_{s_7}(m)=\sqrt{-3}\sum_{Q}\chi_{-3}(Q)\,\frac{\widehat f(\alpha_Q)}{\omega_Q},$$ the sum over the $\Gamma_0(7)$-classes of Heegner forms $[A,B,C]$ of discriminant $-3m^2$ with $7\mid A$, $B\equiv5m\ (14)$ | **[verified, every $m\le45$ with $7\nmid m$, $\ge70$ significant digits]** §2 |
| **R2** | $f$ is weight $-2$, weakly holomorphic, holomorphic on $\mathbf H$ (vanishing to order $1$ at the polar CM points), $f|W_7=-f$; it is the unique such form with pole order $\le1$ at both cusps | **[proved]** §2.1 |
| **R3** | $\widehat f$ is $\Gamma_0(7)$-invariant of weight $0$; computed two independent ways (closed form via $E_2^*,E_4$ and $\mathrm{SL}_2(\mathbf Z)$-reduction; $q$-series with $\Gamma_0(7)+7$-reduction) agreeing to $50$ digits | **[proved]** (invariance) + **[num, 50]** §2.2 |
| **R4** | **The formula is general, not special to Cooper.**  With $f=E_4E_6/\Delta=q^{-1}-240-\cdots$ on $\mathrm{SL}_2(\mathbf Z)$ and the same $\widehat f$: $\beta_{F_{4a}}(m)=\mathrm{Tr}_{-3m^2}(\widehat f)/192$ and $\beta_{F_{4b}}(m)=-\mathrm{Tr}_{-4m^2}(\widehat f)/432$ for Paşol–Zudilin's $F_{4a}=\Delta/E_4^2$, $F_{4b}=E_4\Delta/E_6^2$ | **[verified, $m\le12$, 50 digits]** §2.5 |
| **R5** | **The constant is the polar coefficient.**  If $\Phi=A_2(\tau-\tau_0)^{-2}+\cdots$ then $\lambda=4\pi^2A_2$ up to sign: $4\pi^2A_2=1/192,\ -1/432$ for $F_{4a},F_{4b}$ (exactly the observed constants), and $|4\pi^2A_2|=\nu^2/|g'(u_0)|=\sqrt3,\,2,\,2/\sqrt3$ for $s_7,s_{10},s_{18}$ — matching $|\sqrt{-3}|=\sqrt3$ for $s_7$ | **[proved]** (the level-one values) + **[verified]** §2.5 |
| **R6** | **Consequence: the master conjecture is a divisibility of traces of singular moduli.**  $n^2\mid\beta(n)$ $\iff$ $n^2$ divides the twisted CM trace of the first Shimura–Maass derivative of a weight $-2$ form — the Ahlgren–Ono/Edixhoven setting.  Route (B) of the task is now a precisely stated problem | **[proved, given R1]** §2.6 |
| **R7** | **Why the naive route (A)(1) fails.**  $\beta(m)$ is **not** a $\mathbf Q$-linear combination of twisted CM traces of modular *functions* at the discriminants $D_0m^2$: those are exponential polynomials $\sum_j2\operatorname{Re}(c_jq_0^{-jm})+O(R^{\theta m})$ with constant $c_j$, whereas $\beta(m)$ carries a $1/m$ tail, $\beta(m)=\lambda_0(1-\kappa/m)R^m\sin(\cdot)+O(R^{\theta m})$, $\kappa=1/(2\pi\operatorname{Im}\tau_0)\notin\mathbf Q$.  Numerically $m\beta(m)/((m-\kappa)\tau_H(m))\to3$ to $11$ digits at $m=33$ | **[proved]** + **[verified, $m\le34$]** §3.1 |
| **R8** | Route (A)(3) is closed: $K=\sum\gamma(n)q^n$ is **not** a weakly holomorphic (or holomorphic, or quasi-, or meromorphic) modular form of any weight and level, and $K$ is **not holonomic** over $\mathbf Q(q)$ — its radius of convergence $e^{-\pi\sqrt{|D_0|}/(2N)}$ is transcendental by Gelfond–Schneider, while a holonomic series has algebraic singularities.  The same for $DK$ and $D^2K$ | **[proved]** §3.2 |
| **R9** | $\Phi^\flat:=\sum_e\mu(e)\psi(e)\,e\,\Phi(e\tau)=D^3K$, whose coefficients are $n\beta(n)=a(n^2)$, is **not** a modular form: it has poles at $\tau_0/e$ for every squarefree $e$, infinitely many $\Gamma_0(N)$-inequivalent points | **[proved]**; the ratio $\Phi^\flat/\Phi$ is not a rational function of $u$ of degree $\le6$ **[refuted, 60 coefficients]** §3.3 |
| **R10** | **Dictionary to the half-integral-weight side, calibrated.**  If $\Phi=\Psi(f_{5/2})$ is the Shimura–Borcherds lift (Paşol–Zudilin (5), $k=2$) of $f_{5/2}=\sum a(n)q^n$, then $\beta(m)=a(m^2)/m$ exactly, so the master conjecture is $m^3\mid a(m^2)$.  Calibrated against the literature: for $F_{4a}$ one gets $a(1)=1$, $a(4)=-506$, matching Paşol–Zudilin's $64f_{4a}=q^{-3}+q-506q^4+\cdots$ | **[proved]** + **[verified]** §4.1 |
| **R11** | Equivalent shapes of the master conjecture: $\;\Xi=\sum_e\psi(e)(D^2K)(q^e)$; $\;D^{-3}\Phi=\sum_n\gamma(n)\mathrm{Li}_2^\psi(q^n)$ (a **dilogarithm ladder**); $\;\mathfrak g=\prod_n\bigl((1-q^n)^n\bigr)^{-\gamma(n)}$; $\;\Phi^\flat$ is **triply magnetic**; $\;c(pn)\equiv\psi(p)p\,c(n)\pmod{p^{3(1+v_p(n))}}$ | **[proved]** (equivalences) §4.2 |
| **R12** | $n^2\mid\beta(n)$, all three rows | **[verified, $n\le12000$]** (round 1 had $1500$) `FINDINGS_MODP.md` |
| **M1–M10** | **A genuine mod-$p$ structure exists, and only at $p=2$ (all rows) and $p=3$ ($s_{18}$).**  $\gamma(n)$ is odd $\iff$ $2^{1+v_2(N)}\nmid n$ and $P\nmid n$ ($P=7,5,3$); so $\gamma\bmod2$ is multiplicative and $K\bmod2$ is a **rational** function of $q$ (periods $14,20,12$).  For $s_{18}$, $K\equiv(q+2q^2+2q^3)/(1+q)^3\pmod3$.  Mod $2$ the whole conjecture degenerates to a weight-$3$ Eisenstein/theta identity: $c'(m)$ is odd $\iff$ the $2P$-free part of $m$ is a perfect square | **[verified, $n\le12000$]** `FINDINGS_MODP.md` |
| **M4, L1–L6, N1** | Outside those four cells $K\bmod p$ is neither rational nor algebraic; all Lucas/Dwork laws for $\gamma$ fail massively; $\gamma(p)\bmod p$ matches nothing in a sweep of several hundred arithmetic candidates over $42$ primes (best score $5/42$ = noise) and is equidistributed | **[refuted, stated ranges]** `FINDINGS_MODP.md` |
| **M9** | New unbounded divisibilities: $v_3(\gamma_{s_7}(n))\ge v_3(n)-1$ and $v_2(\gamma_{s_{10}}(n))\ge v_2(n)-1$, both sharp — so $n^2\mid\beta(n)$ is **not** sharp at $(s_7,3)$ and $(s_{10},2)$, neither of which is an Atkin–Lehner cell | **[verified, $n\le12000$]** `FINDINGS_MODP.md` |
| **W1** | The task's route (A) on the $x$-line is **empty as posed**: in weight $0$ the Hecke correspondence is $T_p=p\,U_p+V_p$, so the Kronecker congruence mod $p$ reads $p\,U_pf\equiv0$, a tautology | **[proved]** `FINDINGS_XSIDE.md` |
| **W3–W7** | Its correct form is a real reduction.  With $H=(\sqrt PF)^{p-1}$, $f|U_p\equiv f_0+\sum_{n\ge1}[x^{pn}](fH)x^n \pmod p$, and $H\equiv P^{(p-1)/2}/F_{<p}$ where $F_{<p}=\sum_{n<p}A_nx^n$ is the **supersingular polynomial** of $X_0(N)/\mathbf F_p$ in the coordinate $x$ ($\deg F_{<p}=\lfloor\mu(N)(p-1)/(6\deg x)\rfloor$ exactly, $20$ primes $\times$ $3$ rows).  So round 1's smallest missing brick becomes $$\bigl[x^{pn}\bigr]\Bigl(l(x)\,P(x)^{(p-1)/2}/F_{<p}(x)\Bigr)\equiv\psi(p)L_n\pmod p$$ — an ordinarity statement at the supersingular divisor | **[proved]** (equivalence) + **[verified, $p\le53$]** `FINDINGS_XSIDE.md` |
| **W11** | **New supercongruence:** $A_{p-1}\equiv\kappa\,\psi(p)\,p\pmod{p^2}$ with $\kappa=-2,-3,-3$; for $s_{18}$ this manufactures $\chi_{-3}$ from the Apéry-like recurrence alone | **[verified, $5\le p\le300$]**, **[proved for $s_{10}$]** `FINDINGS_XSIDE.md` |
| **W20** | The **$\psi$-twisted Dieudonné–Dwork tower** $a_{pn-1}\equiv\psi(p)a_{n-1}\pmod{p^{1+v_p(n)}}$ holds for all three rows, extending round 1's V19 to the character row | **[verified, $p\le31$]** `FINDINGS_XSIDE.md` |
| **GAP** | $n^2\mid\beta(n)$: still open.  It is now exactly the statement that the twisted CM traces $\mathrm{Tr}_{D_0n^2}(\widehat f)$ are divisible by $n^2$ | **open** §7 |

**Headline.**  The kernel $K$ of Conjecture 4.1 is identified: $\gamma(m)=\lambda\,m^{-2}
\mathrm{Tr}_{D_0m^2}(\widehat f)$, where $\widehat f$ is the first Shimura–Maass derivative
of the weight-$(-2)$ weakly holomorphic form $f=1/(xF)$ and $\mathrm{Tr}$ is the
genus-character-twisted Heegner trace.  This is *not* a peculiarity of Cooper's rows — the
same formula computes $\beta$ for Paşol–Zudilin's level-one magnetic forms, with
$\lambda=4\pi^2A_2$ the polar coefficient in both cases — so the identity is a theta lift,
and **what is special about Cooper's rows is purely the arithmetic of the traces**:
$m^2\mid\beta(m)$ holds for them and fails already at $m=2$ for Paşol–Zudilin's.  The master
conjecture is therefore an Ahlgren–Ono/Edixhoven divisibility of traces of singular moduli,
for a Maass-raised weight $-2$ input rather than for a modular function.  Round 1's
"smallest missing brick" is independently reduced (`FINDINGS_XSIDE.md`) to a finite,
explicit congruence whose kernel is the supersingular polynomial.

---

## 1. Setting, and what round 1 left

Cooper's rows $s_7,s_{10},s_{18}$: $N=7,10,18$; $u$ the eta quotient with $u|W_N=1/(Cu)$;
$g(u)=1+Bu+Cu^2$ with $(B,C)=(13,49),(6,25),(14,1)$; $x=u/g(u)$; $F=D\log u\in
M_2(\Gamma_0(N))$; $\Phi=F\,Dx=F^2u(1-Cu^2)/g(u)^2=\sum_{m\ge1}c(m)q^m$, a meromorphic
weight-$4$ form, holomorphic at all cusps, with a **double** pole at the CM points
$\tau_0$ where $g(u)=0$, of discriminants $D_0=-3,-4,-36$; $\psi=\mathbf1,\mathbf1,\chi_{-3}$;
$c'(m)=c(m)/m\in\mathbf Z$ (magnetism); $\beta=c'\star(\mu\psi)$; $\gamma(n)=\beta(n)/n^2$.

Round 1 proved eq:magnetic $\iff$ $(S)$ $\iff$ $\operatorname{rad}(n)^2\mid\beta(n)$, proved
two Atkin–Lehner cells outright, identified the mechanism as a Cartier-eigenvector statement
on the $x$-line, and left two things open: **(i)** what $\gamma$ *is*, and **(ii)** the
mod-$p$ brick $\mathcal C(\eta)\equiv\psi(p)\eta$.  This note answers (i) and reduces (ii).

Numerical baseline recomputed from scratch: $c',\beta,\gamma$ for $n\le12000$, all rows,
exact; $n^2\mid\beta(n)$ throughout (`01_data.gp`, `20_extend12000.gp`).

---

## 2. The identification of $\gamma$

### 2.1 The weight $-2$ form **[proved]**

**Lemma 2.1.**  $f:=\dfrac1{x\,F}=\dfrac{g(u)}{u\,F}$ is a weakly holomorphic modular form
of weight $-2$ on $\Gamma_0(N)$, holomorphic on $\mathbf H$, vanishing to order $1$ (in
$\tau$) at each polar CM point, with $f|_{-2}W_N=-f$ and integral $q$-expansion
$$f_{s_7}=q^{-1}+5+19q+52q^2+137q^3+316q^4+695q^5+1440q^6+2887q^7+5564q^8+\cdots$$

*Proof.*  $x|W_N=x$ and $F|_2W_N=-F$ (round 1, Lemma 6.1) give the weight and the sign.
$x$ has a pole of order $\nu=\operatorname{ord}_\tau(u-u_0)$ at $\tau_0$ ($\nu=3,4,4$), because
$g(u_0)=0$ simply.  $F=Du/u$ vanishes there to order $\nu-1$.  Hence $f=(1/x)(1/F)$ has
order $\nu-(\nu-1)=1$.  Away from $\tau_0$, $F$ does not vanish: the divisor of $F$ on
$X_0(N)$ has degree $2\mu(N)/12=4/3,\,3,\,6$, is supported off the cusps (constant term $1$
at $\infty$, and $F|W_N=-F$), and the polar CM points already account for all of it
($2\cdot\frac13\cdot 2$ for $s_7$; $3\cdot\frac12\cdot2$ for $s_{10}$; $3+3$ for $s_{18}$).
Zeros of $x$ occur only at the zeros of $u$, which are cusps. $\square$

For weight $-2$ with pole order $\le1$ at both cusps the space is one-dimensional
(any such form is $G/F$ with $G$ a modular function vanishing at $u=u_\pm$, i.e. divisible
by $g(u)$, and the pole condition forces $G=cg(u)/u$), so $f$ is canonical.

### 2.2 The Shimura–Maass derivative **[proved + num, 50 digits]**

Let $R_k=2i\,\partial_\tau+k\,y^{-1}$ be the Maass raising operator ($y=\operatorname{Im}\tau$),
so that $R_{-2}$ carries weight $-2$ to weight $0$, and set
$$\widehat f:=-\frac{R_{-2}f}{4\pi}=Df+\frac{f}{2\pi y},\qquad D=q\frac d{dq}.$$
$\widehat f$ is a **real-analytic $\Gamma_0(N)$-invariant function**, with
$\widehat f|W_N=-\widehat f$.  Its CM values are algebraic (Shimura–Maass); e.g.
$\widehat f\bigl(\frac{-5+\sqrt{-3}}{14}\bigr)=-\sqrt{-3}$.

Two independent evaluations agree to $50$ digits (`maass.gp`, `maass2.gp`, `t3.gp`, `t5.gp`):

* closed form $\widehat f=-f^2\bigl(F\,Dx+x\,DF\bigr)+\dfrac f{2\pi y}$ with
  $Dx=Fu(1-Cu^2)/g(u)^2$ and $DF=\sum_d\lambda_d\,d\,(E_2(d\tau)^2-E_4(d\tau))/12$, the
  $E_2^*,E_4$ evaluated by reduction to the $\mathrm{SL}_2(\mathbf Z)$ fundamental domain;
* the $q$-expansions of $f$ and $Df$ after reduction in $\Gamma_0(N)+N$ (used at the zeros
  of $F$, i.e. exactly at the polar CM points, where the closed form is $0\cdot\infty$).

### 2.3 The trace identity **[verified, $m\le45$, $\ge70$ significant digits]**

For $d=D_0m^2$ let $\mathcal Q_m$ be the set of $\Gamma_0(N)$-classes of Heegner forms
$[A,B,C]$ with $N\mid A$, $B\equiv\delta_0m\pmod{2N}$ ($\delta_0^2\equiv D_0\bmod 4N$;
$\delta_0=5$ for $s_7$), $B^2-4AC=D_0m^2$; let $\chi_{D_0}$ be the genus character and
$\omega_Q=|\overline{\Gamma_0(N)}_{\alpha_Q}|$.  Put
$$\operatorname{Tr}_m(\widehat f)=\sum_{Q\in\mathcal Q_m}\chi_{D_0}(Q)\,
\frac{\widehat f(\alpha_Q)}{\omega_Q},\qquad \alpha_Q=\frac{-B+\sqrt{D_0}\,m}{2A}.$$

$$\boxed{\ \beta_{s_7}(m)=\sqrt{-3}\;\operatorname{Tr}_m(\widehat f)\ }\qquad(7\nmid m)$$

**[verified]** for every $m\le45$ with $7\nmid m$: the worst absolute error is
$4\cdot10^{-61}$ against values up to $2.3\cdot10^{15}$, and the number of $\Gamma_0(7)$
Heegner classes found equals the number of $\mathrm{SL}_2(\mathbf Z)$-classes of
discriminant $-3m^2$ in every single case (`13_ident7.gp/.log`, `17_ident7wide.gp/.log`).
Equivalently
$$\gamma_{s_7}(m)=\frac{\sqrt{-3}}{m^2}\operatorname{Tr}_m(\widehat f),\qquad
a(m^2)=m\,\sqrt{-3}\operatorname{Tr}_m(\widehat f).$$
For $7\mid m$ the two $\beta$-classes $\pm5m$ merge, the trace becomes real, and
$\beta(m)=0$ by round 1's Theorem 6.3.

### 2.4 The other two rows

See `FINDINGS_ROWS.md`.  The same law holds for $s_{10}$ ($N=10$, $D_0=-4$, $\beta$-class
$6m\bmod20$, genus character $\chi_{-4}$) with the trace producing clean integers
$2i\operatorname{Tr}(m)=2,-4,-18,32,\ast,-108,294,256,-1782,\ast,6534,-3744,-14534,22148,
\ast,-79872$ against $\beta=1,-4,-9,32,0,-108,147,256,-891,0,3267,-3744,-7267,22148,0,-79872$
($\ast$ at $5\mid m$, where $\beta=0$ by the second Atkin–Lehner cell) — i.e. the identity
with $|\lambda|=2$ once the Heegner-set convention at even/odd $m$ is normalised (§2.5
predicts $|\lambda|=2$; see `FINDINGS_ROWS.md` for the resolution).  $s_{18}$ is more
delicate because $D_0=-36$ is not fundamental.

### 2.5 The formula is general; the constant is the polar coefficient **[verified]**

The same construction at level one (`19_pzcontrol2.gp/.log`): $f=E_4E_6/\Delta=
q^{-1}-240-141444q-\cdots$ is *the* weight-$(-2)$ weakly holomorphic form on
$\mathrm{SL}_2(\mathbf Z)$ with pole order $1$, and with the same $\widehat f$ and the plain
twisted traces,
$$\beta_{F_{4a}}(m)=\frac{\operatorname{Tr}_{-3m^2}(\widehat f)}{192},\qquad
\beta_{F_{4b}}(m)=-\frac{\operatorname{Tr}_{-4m^2}(\widehat f)}{432}
\qquad\textbf{[verified, }m\le12\textbf{, 50 digits]}$$
for Paşol–Zudilin's $F_{4a}=\Delta/E_4^2$ (double pole at $\rho$) and $F_{4b}=E_4\Delta/E_6^2$
(at $i$).  Writing $\Phi=A_2(\tau-\tau_0)^{-2}+\cdots$, one computes
$A_2(F_{4a})=\Delta(\rho)/E_4'(\rho)^2=1/(768\pi^2)$ and
$A_2(F_{4b})=E_4(i)\Delta(i)/E_6'(i)^2=-1/(1728\pi^2)$, so in both cases
$$\lambda=4\pi^2A_2 .$$
For Cooper's rows `cooper_sources` §2.3 gives $A_2=\nu^2/(4\pi^2g'(u_0))$, hence
$$|\lambda|=\frac{\nu^2}{|g'(u_0)|}=\sqrt3\ (s_7),\qquad 2\ (s_{10}),\qquad
\tfrac2{\sqrt3}\ (s_{18}),$$
and $|\sqrt{-3}|=\sqrt3$ confirms $s_7$.  (The sign, and a possible factor $2$, depend on
which of the two Fricke-conjugate CM points the $\beta$-class selects.)

**This is the decisive structural point.**  The trace identity is a *theta lift* — it is the
Bruinier–Funke/Katok–Sarnak correspondence "integral weight $2-2k$ $\to$ half-integral weight
$k+1/2$" at $k=2$, whose input weight is $-2$ and output weight $5/2$, exactly the
Shimura–Borcherds input of Paşol–Zudilin.  It holds for Cooper's rows and for
Paşol–Zudilin's level-one forms alike.  What is special about Cooper's rows is therefore
**not** the shape of $\gamma$ but the *arithmetic of the traces*: $m^2\mid\beta(m)$ holds
for them and fails at once for the level-one forms ($\beta_{F_{4a}}(2)=-253$).

### 2.6 The master conjecture as a divisibility of traces of singular moduli

$$\boxed{\ n^2\mid\beta(n)\ \iff\ n^2\ \Bigm|\ \lambda\operatorname{Tr}_{D_0n^2}(\widehat f)\ }$$
i.e. exactly the Ahlgren–Ono/Edixhoven phenomenon — divisibility of traces of singular
moduli at discriminants $D_0p^{2k}$ by powers of $p$ — but for the **first Shimura–Maass
derivative of a weight $-2$ form** in place of a modular function.  Route (B) of the task is
now a precisely stated problem, and its $p$-adic ingredients are available: the disc-$D_0n^2$
CM points lie in the $p$-adic residue discs of the disc-$D_0$ points, and on the ordinary
locus $R_{-2}$ agrees at CM points with the Katz $\theta$-operator under the Katz $p$-adic
period comparison, so the "one derivative" of $R_{-2}$ is expected to contribute the extra
factor $p$ on top of the factor $p$ that Edixhoven's argument gives for a weight-$0$ input —
which is precisely the observed $p^2$.

---

## 3. What $\gamma$ is **not**: three routes closed

### 3.1 Not a CM trace of a modular function **[proved]**

`NOTES_TRACE.md`, `02_trace7.gp`–`07_obstr7.gp`.  For a modular function $G$ on $X_0(N)$
with poles only at the cusps, the Heegner point of largest height at discriminant $D_0m^2$
is $\alpha=-m\bar\tau_0$ (height $m\operatorname{Im}\tau_0$), together with its Fricke
partner, all other classes having height $\le\frac12m\operatorname{Im}\tau_0$; hence
$$\operatorname{Tr}^{\mathrm{tw}}_m(G)=\sum_{j\ge1}2\operatorname{Re}\bigl(c_jq_0^{-jm}\bigr)
+O(R^{\theta m}),\qquad \theta<1,$$
with **constant** $c_j$ (the principal-part coefficients of $G$).  By `cooper_sources` §2.4
the polar data of $\Phi$ give
$\beta(m)=-\nu^2\,2\operatorname{Re}\bigl(q_0^{-m}/g'(u_0)\bigr)\bigl(1-\kappa/m\bigr)
+O(R^{\theta m})$ with $\kappa=1/(2\pi\operatorname{Im}\tau_0)$.  The $1/m$ tail has size
$R^m/m\gg R^{\theta m}$ and is not an exponential polynomial with constant coefficients, so
$\beta$ is not a finite $\mathbf Q$-linear combination of such traces.  $\square$

Numerically, with $H=1/u-49u$ the unique Fricke-anti-invariant Hauptmodul combination with
poles of order $\le1$, $\operatorname{Tr}^{\mathrm{tw}}_m(H)=i\sqrt3\,\tau_H(m)$ with
$\tau_H(m)\in\mathbf Z$ ($1,-6,5,18,-57,90,\ast,-458,1227,-1138,\dots$) and
$$\frac{m\,\beta(m)}{(m-\kappa)\,\tau_H(m)}\longrightarrow3
\quad(3.000000000012\ \text{at }m=33),$$
while an exact three-parameter fit against the full basis $\{1,1/u-49u,1/u+49u\}$ solved at
$m=22,23,24$ misses $m=25,\dots,34$ by $\approx1\%$ — a relative error of size $\kappa/m$, not
$R^{-m/2}$.  The $E_2^*$-completed control also fails: numerically
$\sum_Q\chi\omega^{-1}(E_2^*H/F)(\alpha_Q)=-3\operatorname{Tr}^{\mathrm{tw}}_m(H)+O(1)$,
producing no new sequence (`06_e2trace7.log`).

The moral is exactly §2: the missing $1/m$ is supplied by the $f/(2\pi y)$ term of the Maass
raising operator, i.e. by the fact that $\Phi$ has a **double** pole and $\Xi=D^{-1}\Phi$ is
an Eichler integral.

### 3.2 $K$ is not modular and not holonomic **[proved]**

$|\gamma(n)|\sim c\,R^n/n^2$ with $R=e^{2\pi\operatorname{Im}\tau_0}=e^{\pi\sqrt{|D_0|}/(2N)}
>1$ (`09_gamma_probe.log`, and $|\gamma(n)|^{1/n}\to R$ from below).  Hence:

* $K$ is not a weakly holomorphic modular form of any weight and level: those have
  coefficients $O(e^{c\sqrt n})$.  A fortiori not holomorphic, not quasi-modular.
* $K$ is not a *meromorphic* modular form either: those have coefficients
  $\mathrm{poly}(n)\cdot R^n$, with no $1/n^2$ factor.  Same for $DK$ ($\sim R^n/n$) and
  $D^2K=B=\sum\beta(n)q^n$ ($\sim R^n(1-\kappa/n)$).
* $K$ is not holonomic over $\mathbf Q(q)$: a $D$-finite series has its singularities among
  the roots of the leading polynomial, hence an algebraic radius of convergence, but
  $R^{-1}=e^{-\pi\sqrt{|D_0|}/(2N)}=e^{-\pi\sqrt d}$ with $d=3/49,\,1/25,\,1/9$ rational
  positive, which is transcendental by Gelfond–Schneider ($e^{\pi\sqrt d}=(-1)^{-i\sqrt d}$).

### 3.3 $\Phi^\flat$ is not modular **[proved]**

$\Phi^\flat:=\sum_{e\ge1}\mu(e)\psi(e)\,e\,\Phi(e\tau)=D^3K=DB$ has $[q^n]\Phi^\flat=n\beta(n)
=a(n^2)$ **[verified, $n\le400$, `08_struct.gp`]**, and $\Phi^\flat\approx\Phi$
asymptotically.  But $\Phi(e\tau)$ has poles at $\tau_0/e$, so $\Phi^\flat$ is meromorphic on
$\mathbf H$ with poles at infinitely many $\Gamma_0(N)$-inequivalent points (accumulating at
the cusp $0$), hence is not a modular form.  Consistently, $\Phi^\flat/\Phi$ is not a
rational function of $u$ of degree $\le6$ (`10b_flat.log`, $60$ exact equations, kernel $0$).

---

## 4. The dictionary and the equivalent shapes of the conjecture

### 4.1 $\beta(m)=a(m^2)/m$ **[proved; calibrated]**

If $\Phi=\Psi(f_{5/2})$ is the Shimura–Borcherds lift (Paşol–Zudilin (5), $k=2$, $D=1$) of
$f_{5/2}=\sum a(n)q^n$, i.e. $c(n)=\sum_{d\mid n}\psi(d)\,d\,a(n^2/d^2)$, then
$$c'(n)=\sum_{m\mid n}\frac{a(m^2)}m\quad\Longrightarrow\quad \beta(m)=\frac{a(m^2)}m,$$
so **the master conjecture is $m^3\mid a(m^2)$**, three powers more than Paşol–Zudilin's
theorem $m\mid a(m^2)$ delivers.  Calibration against the literature (`18_pzcontrol.gp`):
for $F_{4a}$ one gets $\beta=1,-253,60083,-14090996,\dots$, hence
$a(m^2)=m\beta(m)=1,-506,180249,\dots$, and $a(1)=1$, $a(4)=-506$ agree with Paşol–Zudilin's
$64f_{4a}=q^{-3}+q-506q^4+\cdots$ (Lemma 1(a)).  See also `FINDINGS_PZ.md`.

### 4.2 Equivalent shapes **[proved]**

With $K=\sum\gamma(d)q^d$ and $\Lambda_\psi$ the $\psi$-Lambert kernel, the following are
equivalent to $\gamma\in\mathbf Z$ (i.e. to $n^2\mid\beta(n)$):

1. $\Xi=\sum_{e\ge1}\psi(e)\,(D^2K)(q^e)$, i.e. $c(n)=n\sum_{d\mid n}\psi(n/d)d^2\gamma(d)$
   — $\Xi$ is a **$\gamma$-twisted weight-3 Eisenstein series** $E_3^{\psi,\mathbf1}$ with the
   Lambert kernel replaced by $K$;
2. the **dilogarithm ladder** $\;D^{-3}\Phi=\sum_{n\ge1}\gamma(n)\,\mathrm{Li}_2^\psi(q^n)$,
   $\mathrm{Li}_2^\psi(t)=\sum_k\psi(k)t^k/k^2$; one derivative down,
   $D^{-2}\Phi=\log\mathfrak g$ with $\mathfrak g=\prod_n\bigl((1-q^n)^{n}\bigr)^{-\gamma(n)}$;
3. $\Phi^\flat=\sum_e\mu(e)\psi(e)e\Phi(e\tau)$ is **triply magnetic**
   ($D^{-3}\Phi^\flat\in\mathbf Z[[q]]$);
4. $c(pn)\equiv\psi(p)\,p\,c(n)\pmod{p^{3(1+v_p(n))}}$ for every prime $p$ and every $n$
   — the "eq:magnetic with multiplicity" form;
5. $m^3\mid a(m^2)$ for the weight-$5/2$ Shimura–Borcherds input;
6. $n^2\mid\lambda\operatorname{Tr}_{D_0n^2}(\widehat f)$ (§2.6).

(All six are formal rewritings; (1)–(3) verified to $n\le400$ in `08_struct.gp`, (4)
follows from (1) by the computation of round 1 Prop. 2.2 applied at each prime power.)
