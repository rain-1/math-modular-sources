# Three ChatGPT research notes: verification, scoring, and what to absorb

*Fable, 2026-08-23.  Scripts: `lattice/sol_notes/` (`01_cellular.gp`, `02_recur.gp`,
`03_moving.py`, `03b_branch.py`, `04_sig4.gp`, `04b_sig4num.py`, `05_emn_*`,
`06_c3_*`, `07_types.gp`, `08_score.gp`).  Sources reviewed:
`chatgpt-research-archive/Research Program_ A Mixed-Motive–Apéry–Holonomy Attack…md` (**N1**),
`…/A New Mixed-Motive and Cyclotomic Apéry Programme for Catalan's Constant.md` (**N2**),
`…/The Conductor-3 Apéry Surface_ Picard Maximality and Explicit Branch Reconstruction.md` (**N3**).
Judged against `CATALAN_OBSTRUCTION.md`, `CATALAN_MU4.md`, `CDT_UNPACKED.md`,
`ADELIC_HOLONOMY.md` §2.6, `CDT_FINDER.md` §6.
The fourth note (Eskandari–Murty–Nemoto, "A Non-Eisenstein Catalan Programme…") is
being verified separately; §2 below overlaps it heavily and says so where it does.*

---

## 0. Verdicts in one line each

| note | verdict | best new item | best entry achieved | wall it hits |
|---|---|---|---|---|
| **N1** (programme) | **programme, not results**; its only concrete object is N2's cellular row; two factual errors | the cellular row (see N2) | $-1.59$ (that row) | width law (2.2): $16^{1/w}\lvert t_2\rvert\le4$ becomes $\le\varphi^5$, but $k$ doubles |
| **N2** (three realisations) | **largely verified; one central quantitative no-go refuted, one closed form corrected**; the cellular Apéry system is real and new to us | $C_n$: an exact order-2 Catalan Apéry system with $\pi,\log2,\pi^2$ all absent, on the *golden* host | $-1.59$ (cellular row); $+1.27$ **counterfactual** for the moving EMN period | cellular row: width law; moving EMN: **no fold exists** (F3), not entry; EMN Padé: denominators beat the $(16/27)^n$ gain by $15.4$ nats/$n$ |
| **N3** (conductor 3) | **arithmetic all verified; one stated isomorphism is false as written** | $H^1(E,j_*W)\simeq T(X)\oplus\mathbf Q(-1)^3$, $\operatorname{rank}T=2$; the $C_6$ dimension vector $(2,0,1,1,1,0)$ | n/a (not a Catalan host) | n/a; but $\Lambda^2_0W$ needs a $\chi$-twist, and $\lambda=1$ is the one unchecked input |

**Single most useful thing in the three notes:** the moving EMN period
$H(z)=z\iint_\Delta\frac{dxdy}{1-z(x^2+y^2)}$ is the *first* Catalan host in this
programme whose conformal/denominator budget would **pass** the CDT entry test
($+1.27$), and it dies for a completely different reason — the logarithmic branch at
the Catalan point has coefficient $\tfrac12$, independent of $G$, so no conditional
function exists at all.  That relocates the Catalan obstruction from "wrong conformal
size" to "no rationality-sensitive monodromy" for this class of host, and tells us
exactly what to look for next: a host with $G$ as an *Apéry limit* and $\lvert t_2\rvert=O(1)$.

**Second headline.**  For the EMN Padé families the arithmetic is now measured rather than
conjectured: the archimedean side really is solved ($(16/27)^n$, verified exactly), and the
companion costs $15.886$ nats per $n$ against a $0.523$-nat gain.  The "$2$-adic slope
$\approx6$" that N2 §22 hopes for is exactly $\log_264$ — the archimedean saddle counted a
second time — the cleanest worked example of (F5)/(2.3) we have.

---

# 1. N1 — "Research Program: A Mixed-Motive–Apéry–Holonomy Attack…"

## 1.1 Content (15 lines)

A research *plan*, not a set of results.  Premise: the Catalan obstruction is the
arithmetic cost of the realisation, not the difficulty of producing $G$; so replace the
weight-3 Eisenstein host by the two-dimensional mixed motive $C$ of
Eskandari–Murty–Nemoto (arXiv:2510.20648), whose period space is intrinsically
$\langle1,G\rangle_{\mathbf Q}$, with basic period
$G=\iint_\Delta\frac{dxdy}{1-x^2-y^2}$ on $\Delta=\{x,y\ge0,x+y\le1\}$ and order-four
symmetry $\sigma(x,y)=(-y,x)$.  The plan is to *move* $C$ in a one-parameter family
$C_z$ keeping the rank-two $\sigma$-quotient, compute its Gauss–Manin connection and
minimal Picard–Fuchs operator, find an integral cusp, extract an Apéry pair
$B_n/A_n\to G$, and run adelic arithmetic holonomy on the rank-2 quotient rather than
on the full $\mu_4$ polylogarithmic module.  Work packages A–H, three candidate
rationality-dependent mechanisms (conditional splitting / connection cancellation /
conditional linear dependence — the last preferred), five failure criteria F1–F5, five
conjectures A–E.  The one concrete mathematical object in the note is the deformation
$\Phi(z)=\iint_\Delta\frac{dxdy}{1-z(x^2+y^2)}$ with $H=z\Phi$, and the cyclotomic
cellular benchmark row $C_n$ of §11 — both of which are also in N2 and are treated
there.  Everything else is scaffolding, and the scaffolding is competent: F1–F5 are
exactly the right stopping conditions and F3/F4 are the ones that actually fire.

## 1.2 Checkable claims

1. $G=\iint_\Delta\frac{dx\,dy}{1-x^2-y^2}$, $\Delta=\{x,y\ge0,\ x+y\le1\}$.
2. $\Phi(z)=\iint_\Delta\frac{dxdy}{1-z(x^2+y^2)}$ satisfies $\Phi(1)=G$.
3. $H=z\Phi$ satisfies $H'(z)=\operatorname{artanh}\sqrt{z/(2-z)}\big/\sqrt{z(2-z)}$.
4. With $q=\sqrt{z/(2-z)}$: $H=2\int_0^q\frac{\operatorname{artanh}u}{1+u^2}du$; at $q=1$, $G=2\int_0^1\frac{\operatorname{artanh}u}{1+u^2}du$.
5. §11: $C_n=\iint_{[0,1]^2}\frac{x^{2n}(1-x^4)^ny^{2n}(1-y^4)^n}{(1+x^2y^2)^{2n+1}}dxdy = A_nG+B_n$ with $A_n,B_n\in\mathbf Q$ and no $\pi,\log2,\pi^2$.
6. $U_n=(-1)^n16^nA_n$: $U_0=1,U_1=28,U_2=2596,U_3=311536,U_4=41759524$, all integral.
7. $U_n,V_n$ satisfy $(n+2)^2(2n+3)^2(20n^2+32n+13)X_{n+2}-4P(n)X_{n+1}-256(n+1)^2(2n+1)^2(20n^2+72n+65)X_n=0$, $P(n)=3520n^6+26752n^5+83024n^4+134592n^3+120196n^2+56088n+10699$.
8. $\lambda_\pm=88\pm40\sqrt5$; $|\lambda_-|=40\sqrt5-88=1.442719\ldots>1$, so no irrationality.
9. $\operatorname{den}(V_n)\mid d_{2n}^2$.
10. §17: $F_a(t)={}_2F_1(a,1-a;1;t)$, $D=\partial_aF_a|_{a=1/4}$; $\psi_1(1/4)-\psi_1(3/4)=16G$; $d_n[x^n]D(64x)\in\mathbf Z$ ("one LCM layer").
11. §2: "Viola and Marcovecchio obtain effective approximants with exponent approximately $0.6293$".
12. §16: $I_{F,t}=\iint_\Delta\frac{F}{(1-x^2-y^2)^{t+1}}$ reduces to linear forms in $1,G$.

## 1.3 Verification

| # | verdict | evidence |
|---|---|---|
| 1 | **verified, exactly** | $(s,d)\!=\!(x{+}y,x{-}y)$, $d=sw$ reduces it to $\int_0^1\frac{\log(2/(1-w^2))}{1+w^2}dw=\frac\pi4\log2-\bigl(\frac\pi4\log2-G\bigr)=G$, using $\int_0^1\frac{\log(1\pm w)}{1+w^2}dw=\frac\pi8\log2\mp G$. |
| 2 | **verified** | via 3–4 below ($H(1)=2\int_0^1\!\operatorname{artanh}u/(1+u^2)du=G$ to 40 digits); direct partial sums $\sum a_n$ with $1/N$-Richardson agree to $1.6\cdot10^{-8}$ (`03b_branch.py`). |
| 3 | **verified**, 40 digits at $z=0.3,0.7,0.95$ | `03_moving.py` A5 |
| 4 | **verified**, 40 digits ($q=0.4,0.8,1$) | `03_moving.py` A6 |
| 5 | **verified, exactly, to $n=34$** | `01_cellular.gp`: the $\pi$-coordinate of $C_n$ is $0$ for every $n\le34$; the $\log2$ and $\pi^2$ coordinates never arise (all exponents even). |
| 6 | **verified, exactly** | $U_{0..4}=1,28,2596,311536,41759524$; $U_n\in\mathbf Z$ for all $n\le34$; also $V_0=0$, $V_1=-26$ as in N2. |
| 7 | **verified, exactly**, residual $0$ for $n=0..32$ (35 exact terms; the note claims 61) | `02_recur.gp` |
| 8 | **verified** | $\lambda^2-176\lambda-256=0$ from the degree-6 leading behaviour $80\lambda^2-14080\lambda-20480$; roots $-1.4427190999\ldots,\,177.4427190999\ldots$; and the measured error $-V_n/U_n-G$ falls by the factor $\lambda_+/|\lambda_-|=122.9919$ per step, matching to 4 digits from $n=4$ on. |
| 9 | **verified** to $n=34$, and **sharp in the exponent**: $d_n^2$ is *not* enough (fails from $n=2$), and $(1/n)\log\operatorname{den}(V_n)$ reaches $3.55$ at $n=34$ against $(1/n)\log d_{2n}^2=3.91$, ratio $0.907$ and drifting up. So the type is $[1..2n]^2$, $\sigma=4$. | `08_score.gp` |
| 10 | **verified** | $[t^n]{}_2F_1(\tfrac14,\tfrac34;1;t)=\binom{4n}{2n}\binom{2n}n/64^n$ **exactly** (Gauss multiplication: $(\tfrac14)_n(\tfrac34)_n=(4n)!/(64^n(2n)!)$); $\psi_1(\tfrac14)-\psi_1(\tfrac34)=16G$ is classical; and $d_n\,[x^n]D(64x)\in\mathbf Z$ for **every** $n\le40$ — moreover the minimal $k$ with $\operatorname{den}\mid\operatorname{lcm}(1..k)$ satisfies $k\le n$ always. The one-LCM-layer claim is **true**, and it is the only Catalan-carrying object in the archive with a single layer. |
| 11 | **unsupported; almost certainly a conflation** | $0.6293$ is not an irrationality exponent for $G$ in any source we hold; $0.62922$ is CDT's *contour-loss factor* $\lvert\psi'(0)\rvert$ (`CDT_FINDER.md` line 94, `CATALAN_MU4.md` §6). Treat as an error. |
| 12 | see §2 (N2 §10) and the separate EMN review | — |

---

# 2. N2 — "A New Mixed-Motive and Cyclotomic Apéry Programme for Catalan's Constant"

## 2.1 Content (15 lines)

Three realisations, honestly scored by the author.  **(a) Signature four.** With
$F={}_2F_1(\tfrac14,\tfrac34;1;t)$, $\int_0^1F^2=16G/\pi^2$, so $G$ sits in the symmetric
square of a signature-4 elliptic system; the $\mathrm{Sym}^2$ ODE and the moment
recurrence with inhomogeneity $8/\pi^2$ put $G/\pi^2$ and $1/\pi^2$ in one rank-3 period
system.  The parameter derivative $D=\partial_aF_a|_{a=1/4}$ satisfies
$t(1-t)D''+(1-2t)D'-\tfrac3{16}D=\tfrac12F$, has $d_n[x^n]D(64x)\in\mathbf Z$ (one LCM
layer), and carries $16G=\psi_1(\tfrac14)-\psi_1(\tfrac34)$ as a *connection coefficient*
— hence no fold, and the note correctly demotes it to a local model.  **(b) A cyclotomic
cellular Apéry system.** $G=\int_0^1\!\!\int_0^1\frac{dxdy}{1+x^2y^2}$; the family
$C_n$ of §5 has period space exactly $\mathbf Q+\mathbf QG$ and satisfies a second-order
recurrence with $\lambda_\pm=88\pm40\sqrt5$; $|\lambda_-|>1$ and
$\operatorname{den}(V_n)\mid d_{2n}^2$, so no irrationality.  **(c) EMN.** The moving
period $\Phi(z)$ has the two-line Picard–Fuchs equation
$z(z-2)H''+(z-1)H'=\frac1{2(z-1)}$ but a $G$-independent $\log$ branch at $z=1$, which
kills the conditional-holonomy strategy; the fixed motive's higher-pole reduction is made
explicit ($b_{m,t}=4^{-m}\binom m{m/2}\binom mt$), a radial Padé family is claimed to give
a normalised real form decaying like $(16/27)^n$ ("free integration"), and the remaining
obstruction is asserted to be the rational companion's denominator ($\approx L_{12n}^2$).
Jacobi and Racah/dual-Hahn weights are exhibited that annihilate the pole-lowering
boundary functionals; a $2$-adic slope $\approx5.9$ and no $3$-adic slope are measured.

## 2.2 Checkable claims

**Signature four (§§2–3).**

1. $\int_0^1{}_2F_1(\tfrac14,\tfrac34;1;t)^2\,dt=16G/\pi^2$.
2. $8t^2(1-t)^2Y'''+24t(1-t)(1-2t)Y''+2(27t^2-27t+4)Y'+3(2t-1)Y=0$ for $Y=F^2$.
3. $8n^3I_{n-1}-(2n+1)(8n^2+8n+3)I_n+2(n+1)(2n+1)(2n+3)I_{n+1}=8/\pi^2$, $I_n=\int_0^1t^nF^2$.
4. $I_0=16G/\pi^2$, $I_1=(8G+\tfrac43)/\pi^2$, hence $1/\pi^2=\tfrac34(I_1-\tfrac12I_0)$.
5. $t(1-t)D''+(1-2t)D'-\tfrac3{16}D=\tfrac12F$.
6. $F(64x)=\sum\binom{4n}{2n}\binom{2n}nx^n$ and $[x^n]D(64x)=\binom{4n}{2n}\binom{2n}n\cdot4\sum_{k<n}\bigl(\tfrac1{4k+1}-\tfrac1{4k+3}\bigr)$, with $d_n[x^n]D(64x)\in\mathbf Z$.

**Cellular (§§4–6).**  7.–11. are N1's claims 5–9 on the same objects, plus $V_0=0$, $V_1=-26$, $-V_n/U_n\to G$, and the exact period reduction of $\iint x^Ay^B/(1+x^2y^2)^r$ in the basis $1,\pi,\log2,G,\pi^2$.

**Moving EMN (§§8–9).**

12. $a_n=\frac1{2^{n+1}(n+1)}\sum_k\binom nk/(2k+1)$.
13. $(n+3)(2n+5)a_{n+2}-(n+2)(3n+5)a_{n+1}+(n+1)^2a_n=0$.
14. $z(z-2)H''+(z-1)H'=\frac1{2(z-1)}$; homogeneous solutions $1$ and $2\arcsin\sqrt{z/2}$.
15. $H(z)=G+\frac{1-z}2\log(1-z)+O(1-z)$ near $z=1$.
16. After $z=2x$: denominator type $e=1$, $b=2$; nearest log singularity at $x=\tfrac12$; conformal radius $\le2$; $\tau^\flat=\tfrac32$; $\log2<\tfrac32$, so "not a viable CDT host".

**Fixed EMN (§§10–22).**  17.–27.: the pole-lowering recurrence; $I_{4,2}=-\frac{49}{384}+\frac9{64}G$; $b_{m,t}=4^{-m}\binom m{m/2}\binom mt$; $b_{4n,2n}\sim\frac1{2\pi n}$; $\max_\Delta\frac{x^4y^4}{g^2}=\frac1{64}$; $\max_\Delta x^4y^4=\frac{\min(g,1-g)^4}{16}$ at fixed $g$; $W(g)=\frac{\min(g,1-g)^4}{16g^2}$; $256\max_gW(g)(1-4g^2)^2=\frac{16}{27}$; $256^nb_n\in\mathbf Z$ with the stated sum; $c\approx4.6240$ / $0.865$ and $q(h)=1-10h+23h^2$ / $0.50119$; $\operatorname{den}(256^na_n)\mid L_{12n}^2$; $\iint_\Delta\frac{1+569x^2y^2-2800x^4y^4}{1-x^2-y^2}=\frac{13}2G$; Racah weights $(-1)^k\binom{2n}k\binom{4n}{2k}$; $v_2(r_n-r_{n-1})\sim5.9n$; $v_3=O(1)$.  *These overlap the separately-reviewed EMN note; see §2.6.*

## 2.3 Verification — signature four

| # | verdict | evidence |
|---|---|---|
| 1 | **verified**, 40 digits: $I_0=1.4849074908430886518158559216=16G/\pi^2$ to the last digit | `04b_sig4num.py` |
| 2 | **verified**: the series of $F^2$ to $O(t^{40})$ is annihilated exactly | `04_sig4.gp` |
| 3 | **verified**, $n=1..6$: both sides $=0.8105694691387021715510357=8/\pi^2$, residual $<10^{-38}$ | `04b_sig4num.py` |
| 4 | **verified**, 40 digits, with the corollary $1/\pi^2=\tfrac34(I_1-\tfrac12I_0)$ | |
| 5 | **verified exactly** (differentiate Gauss at $a=\tfrac14$: $a(1-a)=\tfrac3{16}$, $1-2a=\tfrac12$) | by hand |
| 6 | **verified exactly**: $(\tfrac14)_n(\tfrac34)_n=(4n)!/(64^n(2n)!)$ by Gauss multiplication, so $[t^n]F=\binom{4n}{2n}\binom{2n}n/64^n$; and $d_n[x^n]D(64x)\in\mathbf Z$ for every $n\le40$ — the minimal $k$ with $\operatorname{den}\mid\mathrm{lcm}(1..k)$ satisfies $k\le n$ throughout.  Mechanism: each prime $p\in(n,4n)$ divides $\binom{4n}{2n}\binom{2n}n$ exactly once (Kummer), and for $p\le n$ the prefactor supplies the excess of $v_p$ over $\lfloor\log_pn\rfloor$.  A proof looks routine but is not written out here | `04_sig4.gp` |

**What the note does not say, and should.**  The moment recurrence has characteristic
polynomial $8\lambda^2-16\lambda+8=8(\lambda-1)^2$: a **double root at $1$**.  Numerically
$I_n\asymp(\log n)^2/n$ (measured $I_7/I_6=0.905$), so there is *no recessive solution and
no geometric decay at all*.  The signature-4 moment system is not an Apéry system in any
useful sense — it is a period-relation statement.  N2 never computes this.

## 2.4 Verification — the cellular row, and its host

Claims 7–11 are N1's 5–9: **all verified exactly** (§1.3), to $n=34$, recurrence residual
$0$ for $n=0..32$, $V_0=0$, $V_1=-26$, $-V_n/U_n\to G$ with the error falling by
$\lambda_+/|\lambda_-|=122.9919$ per step.  The exact period reduction used here is
independent of the note's: with $u=xy$,
$$\iint_{[0,1]^2}\frac{x^Ay^B}{(1+x^2y^2)^r}\,dx\,dy=\frac{J^{(r)}_B-J^{(r)}_A}{A-B}\ (A\ne B),\qquad =K^{(r)}_A\ (A=B),$$
$J^{(r)}_m=\int_0^1\frac{u^m\,du}{(1+u^2)^r}$, $K^{(r)}_m=-\int_0^1\frac{u^m\log u\,du}{(1+u^2)^r}$,
with $J^{(r)}_0=\frac1{2^r(r-1)}+\frac{2r-3}{2r-2}J^{(r-1)}_0$,
$K^{(r)}_0=\frac{J^{(r-1)}_0+(2r-3)K^{(r-1)}_0}{2r-2}$, $J^{(1)}_0=\pi/4$, $K^{(1)}_0=G$,
and $J^{(0)}_m=\frac1{m+1}$, $K^{(0)}_m=\frac1{(m+1)^2}$.  All exponents in $C_n$ are even,
so only $1,\pi,G$ can occur; the $\pi$-coordinate is $0$ for every $n\le34$ — "no $\pi$ or
$\log2$ contamination" confirmed, the mechanism being antisymmetry of the $J$-terms in
$i\leftrightarrow j$.

### The host is Zagier's row D / the $X_1(5)$ golden host, rescaled by 16

The row's singular points are the roots of $256x^2+176x-1$:
$$t_1=\frac{-11+5\sqrt5}{32}=\frac{\varphi^{-5}}{16}=0.00563562\ldots,\qquad
t_2=\frac{-11-5\sqrt5}{32}=-\frac{\varphi^{5}}{16}=-0.69313562\ldots,$$
so $t_1t_2=-\tfrac1{256}$ and $16|t_2|=\varphi^5$, $16|t_1|=\varphi^{-5}$ **exactly**
(38 digits, `08_score.gp`).  These are precisely the singular values of Zagier's sporadic
row $\mathbf D$ / the $X_1(5)$ $\mathrm{Sym}^2$ host of `CDT_FINDER.md` §6
($t_1=\varphi^{-5}$, $t_2=-\varphi^5$, $t_1t_2=-1$, $|N(t_2)|=1$), divided by 16.  The
numbers $177.4427$ and $1.4427$ that `CDT_FINDER.md` records there as "modular ceiling
$16|s|$" **are** $\lambda_+$ and $|\lambda_-|$ of this recurrence.  The cellular
construction is therefore not a new host: it is the golden host carrying a conductor-4
twist, and the twist costs exactly $\log16=4\log2$.

**Scoring** (convention in §4).  Host $\mathbf P^1\setminus\{0,t_1,t_2,\infty\}$ — and the
singular set really is only that: writing the recurrence as
$p_2(\theta-2)+x\,p_1(\theta-1)+x^2p_0(\theta)$ with $\theta=x\,d/dx$, the $\theta^6$
coefficient is $80-14080x-20480x^2=80(1-176x-256x^2)$, so the operator's leading
coefficient is $80x^6(1-176x-256x^2)$ and there are **no apparent singularities**; the
factors $20n^2+32n+13$ and $20n^2+72n+65$ (which are each other's shift) cancel out of the
top $\theta$-degree.  Two of the six solutions at $x=0$ are power series (the recurrence has
order 2), namely $A=\sum U_nx^n$ (radius $|t_1|$) and the recessive one (radius $|t_2|$).
Fold $=t_1$: the conditional
function $bB-aA$ has coefficients $\asymp\lambda_-^n$, hence radius $|t_2|$ and regularity
at $t_1$.  $\tau$ from **measured** denominators: $\operatorname{den}(V_n)\mid d_{2n}^2$,
type $[1..2n]^2$, $\sigma=4$; measured $(1/n)\log\operatorname{den}(V_n)=3.55$ at $n=34$
against $3.91$ for $d_{2n}^2$ (ratio $0.907$, rising).  Conformal ceiling
$\log(16|t_2|)=5\log\varphi=2.40606$.  Hence
$$\boxed{\ \text{entry}=5\log\varphi-4=-1.594\ }\qquad(-1.14\text{ at the measured }\sigma=3.55).$$
Raw Apéry score $\log|\lambda_-|+\sigma=0.3665+4=+4.37$: the cleared linear forms *grow*.
Over $\mathbf Q(\sqrt5)$, where the pure module actually has to live (see §4, (2.5)), the
normalised ceiling is $\log16+\tfrac12\log|N(t_2)|=\log16-\tfrac12\log256=0$ and the entry
is $-4$; $-1.594$ is the optimistic single-place figure.

**Which wall.**  It **does** negate (E1) — the second singularity is irrational (quadratic,
$\mathbf Q(\sqrt5)$) with $|t_2|=\varphi^5/16=0.693>\tfrac14$, worth $+1.02$ nats over the
Beauville/level-8 value $\tfrac14$.  But the Hodge depth is realised with **width $2n$, not
$n$**, so $k$ in the width law is $4$, not $2$, and the $+1.02$ is swallowed by $-2$.  It
falls to (2.2) in its width-law form $\mathrm{entry}\le\frac1w\log16+\log|t_2|-k$.  *This
sharpens (E1): an irrational $t_2$ obtained by a conductor twist is free only if the twist
does not enter the denominators, and here the same $\log16$ is paid twice — once in the
coordinate and once in the doubled lcm width.*

## 2.5 Verification — the moving EMN period, and a refuted no-go

| # | verdict | evidence |
|---|---|---|
| 12 | **verified exactly** for $n\le8$; independently $a_n=\sum_j\binom nj\frac{(2j)!(2n-2j)!}{(2n+2)!}$ and the two agree identically | `03_moving.py` A1 |
| 13 | **verified exactly**, residual $0$ for $n=0..27$ | `03_moving.py` A2 |
| 14 | **verified**: the series of $H=z\Phi$ annihilates $z(z-2)H''+(z-1)H'-\frac1{2(z-1)}$ to $O(z^{20})$ exactly; $2\arcsin\sqrt{z/2}$ solves the homogeneous equation | `03_moving.py` A4, A8 |
| 15 | **verified**, and refined: $(H-G)/((1-z)\log(1-z))\to0.4998$ as $1-z\to0$ (from $1-z=10^{-9},10^{-11}$), so the coefficient is $\tfrac12$; the next term is $-(\tfrac\pi8+\tfrac G2)(1-z)$ (fitted $-0.8507$ against $\pi/8+G/2=0.85069$) | `03b_branch.py` |
| 16 | **denominator type verified; the quantitative no-go REFUTED as stated; conclusion survives for a different reason** | below |

**Denominators.**  In $x=z/2$, $[x^{n+1}]H=T_n/(n+1)$ with $T_n=\sum_k\binom nk/(2k+1)$,
and $\operatorname{den}\mid(n+1)\operatorname{lcm}(1,\dots,2n+1)$ for every $n\le60$,
while $(n+1)\operatorname{lcm}(1..n)$ already fails at $n=1$.  So the claimed type
$n[1..2n]$ ($e=1$, $b=2$) is correct as an upper bound.  **Measured**,
$(1/n)\log\operatorname{den}$ is $1.44$–$1.51$ over $n=10..60$: the effective width is
$\approx0.73$, not $2$, because the binomials cancel most of the odd lcm (`07_types.gp`).

**The refuted step.**  N2 §9 concludes "not a viable CDT host" from: nearest singularity at
$x=\tfrac12$ $\Rightarrow$ conformal radius $\le2$ $\Rightarrow$ $\log2=0.693<\tfrac32=\tau^\flat$.
The arithmetic is right and $\tau^\flat=\tfrac32$ is the correct CDT value for $m=2$, $u=1$,
$b=2$ (`ADELIC_HOLONOMY.md` §2.6.3: $\tau=b(1-u^2/m^2)$).  But *conformal radius
$\le4\,\mathrm{dist}(0,\Sigma)=2$ is the **univalent** (Koebe) bound*, and CDT's method is
explicitly multivalent — the $\lambda$-template on $\mathbf P^1\setminus\{0,t_2,\infty\}$
gives $|\varphi'(0)|=16|t_2|$, which is the entire content of the width law.  Here, in
$x=z/2$, the singular set is $\{0,\tfrac12,1,\infty\}$, fold candidate $x=\tfrac12$,
$t_2=1$; **if** the fold could be removed the ceiling would be $\log16=2.7726$ and
$$\text{entry}=\log16-\tau^\flat=2.7726-1.5=+1.27\qquad(+1.68\text{ at the measured }\sigma),$$
comfortably positive — *the first Catalan host in this programme with a positive entry
number*.  The note's quantitative no-go is therefore wrong.

**What actually kills it, and it is fatal.**  By verified claim 15 the $\log(1-z)$ branch at
the Catalan point has coefficient $\tfrac12$ with no $G$ in it, and the other two solutions
($1$, and $2\arcsin\sqrt{z/2}$, analytic at $z=1$ with value $\pi/2$) cannot cancel it.
Hence **under $G\in\mathbf Q$ no combination of the available functions is regular at
$z=1$: there is no conditional function at all, and CDT has nothing to count**
(`CDT_UNPACKED.md` §2: "the fold is the only place the hypothesis enters").  In the note's
own vocabulary this is failure criterion **F3**; in `CATALAN_OBSTRUCTION.md` terms it is not
(2.1)–(2.6) but a prior failure — the mechanism, not the budget.  Without the fold the host
is genuinely 4-punctured and the ceiling falls below $16|t_2|$ anyway.

**Why this matters.**  $G$ here is a *regular value* of the generating function, exactly as
in the $\mu_4$ architecture (`CATALAN_MU4.md` §2: "the fold has to be created, not found")
and exactly as for the signature-4 connection coefficient of §2.3.  All three of N2's
realisations share this defect.  The moving EMN period shows the *budget* problem is
solvable.  The ledger entry should therefore read: **Catalan now needs a host where $G$ is
an Apéry limit at a fold *and* $|t_2|=O(1)$; the two have never yet co-occurred.**

## 2.6 Verification — the fixed EMN motive (§§10–22)

*This material overlaps the separately-reviewed EMN note ("A Non-Eisenstein Catalan
Programme from the Eskandari–Murty–Nemoto Motive", §§8–13); the results below were obtained
independently here and should be reconciled with that review rather than duplicated.*

**Exact engine** (`05_emn_core.py`, `05_emn_fast.py`).  With $s=x+y$, $d=sw$, then
$w=\tan\theta$, $\phi=2\theta$, every $I(m,t)=\iint_\Delta\frac{x^my^m}{(1-x^2-y^2)^{t+1}}$
collapses to closed form.  The decisive structural fact, which the note does not have:
$$\textbf{the entire }G\textbf{-content of the family is the single number }A_0=\int_0^{\pi/4}\log\tan u\,du=-G,$$
because $A_k=-\frac1{2k}\sum_{i<k}\frac{(-1)^i}{2i+1}\in\mathbf Q$ for $k\ge1$ and every other
ingredient lies in $\mathbf Q+\mathbf Q\pi$ with the $\pi$-parts cancelling identically.  So
$I(m,t)\in\mathbf Q+\mathbf QG$ with **no $\pi$ and no $\log2$ ever present** — a stronger
statement than "they cancel".

| # | claim | verdict | value |
|---|---|---|---|
| E0 | $\iint_\Delta\frac{dxdy}{1-x^2-y^2}=G$ | **verified** exactly | also analytically: $\frac\pi4\log2-\int_0^1\frac{\log(1-w^2)}{1+w^2}dw=G$ |
| E1 | $I_{4,2}=-\frac{49}{384}+\frac9{64}G$ | **verified** exactly | |
| E2 | $b_{m,t}=4^{-m}\binom m{m/2}\binom mt$ | **corrected** | $b_{m,t}=(-1)^t4^{-m}\binom m{m/2}\binom mt$.  Magnitudes right for all 25 pairs $m\in\{0,2,4,6,8\}$; **the sign $(-1)^t$ is missing from the note** (all 10 odd-$t$ cases are off by $-1$).  Harmless for the radial family (there $t=2n-2k$ is even) but not in general. |
| E3 | pole-lowering recurrence | **verified** exactly | residual $0$ for $(4,2),(6,2),(6,4),(8,2),(8,4),(8,6),(10,2),(10,6),(10,10),(12,8)$.  *But:* applied twice from $(4n,2n)$ it lands on $(4n-4,2n-4)$, not $(4n-4,2n-2)$ — **it does not close on the central family**, so there is no $O(1)$ recurrence for $I_n$ along this route. |
| E4 | $b_{4n,2n}=\binom{4n}{2n}^2/4^{4n}\sim\frac1{2\pi n}$ | **verified** | $n=1..4$: $\frac9{64},\frac{1225}{16384},\frac{53361}{1048576},\frac{41409225}{1073741824}$; $2\pi n\,b_n=0.99938$ at $n=200$ |
| E5 | $\max_\Delta\frac{x^4y^4}{g^2}=\frac1{64}$ at $x=y=\tfrac12$ | **verified** exactly | the max is on the edge $x+y=1$, where $g=2xy$ and the ratio is $(xy)^2/4$ |
| E6 | $\max x^4y^4=\frac{\min(g,1-g)^4}{16}$ at fixed $g$ | **verified** exactly | $g\le\tfrac12$: boundary, $2xy=g$; $g\ge\tfrac12$: interior, $x=y=\sqrt{(1-g)/2}$ |
| E7 | $256\max_gW(g)(1-4g^2)^2=\frac{16}{27}$ | **verified** exactly | $256W(g)(1-4g^2)^2=16g^2(2g-1)^2(2g+1)^2$ on $g\le\tfrac12$, max $\tfrac{16}{27}$ at $g=\tfrac1{2\sqrt3}$; the $g\ge\tfrac12$ branch only reaches $0.27093$ |
| E8 | $c\approx4.6240$/$0.865$; $q=1-10h+23h^2$/$0.50119$ | **verified**, labelling caveat | $c^*=4.62402313465907692\ldots$, the unique root $>1$ of $16c^6-80c^5+148c^4-480c^3-335c^2-70c-1$, value $4/c^*=0.865047575$; $q$ gives $0.501188684$.  The apparent inconsistency with E7 is a *labelling* one: in $h=g^2$, $(1-4g^2)^2$ is degree **2** and $(1-cg^2)$ is degree **1**.  By degree: $\deg1\to0.86505$, $\deg2\to\tfrac{16}{27}=0.59259$ for $(1-4h)^2$, $0.50119$ for $q$, **true minimax $0.45488$ at $1-10.15385h+22.91127h^2$**; $\deg3\to0.37503$; $\deg4\to0.27387$. |
| E9 | $256^nb_n=\binom{4n}{2n}\sum_j(-4)^j\binom{2n}j\binom{4n}{2n-2j}$ | **verified** exactly, $n=0..5$ | $1,-12,-19740,3330096,2056806180,-775992197552$ |
| E10 | $\operatorname{den}(256^na_n)\mid L_{12n}^2$ | **verified and sharp**, $n=0..30$ | $\operatorname{den}$ is always **odd** (the $256^n$ exactly clears the $2$-part) and for every $n$ the minimal $k$ with $\operatorname{den}\mid\mathrm{lcm}(1..k)^2$ is **exactly the largest prime $\le12n$** ($11,23,31,47,59,71,83,89,107,\dots$).  The square is genuinely needed: a single lcm requires $k=p^2$ with $p\approx6n$. |
| E11 | $\iint_\Delta\frac{1+569x^2y^2-2800x^4y^4}{1-x^2-y^2}=\frac{13}2G$ | **verified** exactly | $\widehat I(0)=G$, $\widehat I(2)=-\frac5{48}+\frac G8$, $\widehat I(4)=-\frac{569}{26880}+\frac{3G}{128}$; the $\pi$ and $\log2$ coordinates are identically zero in each term separately |
| E12 | $v_2(r_n-r_{n-1})\sim5.9n$, $v_3=O(1)$ | **verified / verified in substance** | $v_2$: LSQ slope $5.931$ on $n=2..80$, local slope over the last 20 terms exactly $\mathbf6$.  $v_3$: bounded in $[-15,-2]$, LSQ slope $-0.079\approx0$, but the envelope is $-1.95\log n$, i.e. $O(\log n)$ **not** $O(1)$.  The substantive claim (no $3$-adic slope) stands. |

### The "free integration" claim, priced

N2 §13 concludes "the archimedean side of an Apéry argument has been solved for this
realization".  The archimedean statement is **true and now exact**:
$$\log|256^nL_n|=-0.52397\,n-1.465\log n-2.977,\qquad \log\tfrac{27}{16}=0.52325,$$
so $256^nL_n\ll(16/27)^n$ with the Laplace $n^{-3/2}$ correction visible.  But the
companion is measured to cost
$$\log\operatorname{den}(256^na_n)=15.886\,n-18.64\log n+28.2$$
(the theoretical ceiling from $L_{12n}^2$ being $24n$; so the true denominator is $e^{8.1n}$
*better* than the note's bound, and still).  **Net balance $+15.36$ nats per $n$: the
denominator exceeds the archimedean gain by a factor $\approx30$.**  There is no regime in
this family in which the Padé construction beats its denominators, and E8's degree-2
minimax improvement is worth $0.264$ nats/$n$ against a $15$-nat deficit.

For contrast the **bare** central family $I_n=I(4n,2n)$ is far better balanced and still
loses:
$$\log|I_n|=-4.15963\,n-1.467\log n-2.36\quad(\log64=4.15888,\ \text{i.e. E5}),\qquad
\log\operatorname{den}(a_n)=+10.4825\,n+\ldots,$$
net $+6.32$ per $n$.  Of the $10.48$, exactly $6\log2=4.1589$ is the power of $2$ — which is
why $256^n=2^{8n}$ clears it and why the measured $2$-adic slope is exactly
$\mathbf{6=\log_264}$: **the $2$-adic gain and the archimedean gain are the same gain counted
twice.**  That is the sharpest available statement of why N2 §22's $\sigma_2\approx6$ is not a
resource — it is the archimedean saddle in $2$-adic clothing, exactly the "always score the
net quantity" of `CATALAN_OBSTRUCTION.md` (2.3)/(F5).  The remaining $\approx6.32$ per $n$ is
the odd lcm part, and it is what kills both families.

**Scoring.**  The fixed-EMN Padé families are not CDT hosts at all: no holonomic host, no
fold, no conditional function — they are a lattice/Apéry construction.  On the archive's
worthiness scale (gain/cost) they score
$$\text{radial }(1-4g^2)^{2n}:\ \tfrac{0.523}{15.886}=0.033,\qquad
\text{bare }I(4n,2n):\ \tfrac{4.159}{10.483}=0.397,$$
against $1-\varepsilon$ for the Zudilin–Nesterenko two-row construction
(`CATALAN_AUDIT.md`, Lean).  They are two orders of magnitude and one order of magnitude
away respectively, and the note's reduction of the problem to "find a low-height integral
lattice for the EMN extension" is correct but understates the size of the gap.

---

# 3. N3 — "The Conductor-3 Apéry Surface: Picard Maximality and Explicit Branch Reconstruction"

*(placed here out of order because it is independent of the Catalan material)*

## 3.1 Content (15 lines)

Not about Catalan.  Take the conductor-3 rank-5 Apéry hypergeometric system
$V={}_5F_4\bigl(\tfrac12,\tfrac16,\tfrac16,\tfrac56,\tfrac56;1,1,\tfrac23,\tfrac43;z\bigr)$
on $\mathbf P^1\setminus\{0,1,\infty\}$, and its rank-4 symplectic spin/Kuga–Satake lift
$W$ with $V\simeq\Lambda^2_0W$.  The minimal cover killing the finite monodromy parts has
ramification $(3,2,6)$ and is the equianharmonic elliptic curve $E:v^2=u^3-1$, $z=u^3$;
over $E$, $W$ is $R^1f_*$ of a semistable genus-2 fibration $f:X\to E$.  From explicit
integral symplectic monodromy matrices the note reads off the stable degeneration:
$N_0=M_0^3-I$ has rank 2 with monodromy pairing $\cong A_2$ (theta graph, 3 nodes) at each
of the two points over $z=0$, and $N_\infty=M_\infty^6-I$ has rank 1 with Smith form
$\mathrm{diag}(4,\ldots)$ (a node of thickness 4, i.e. a 4-cycle), giving $\delta_0=10$,
$\delta_1=0$.  Fougeron's parabolic-degree formula gives $\deg_{\rm par}F^{2,0}_V=1/6$,
hence $\lambda=\deg\det f_*\omega=1$; Mumford's $10\lambda=\delta_0+2\delta_1$ closes.
Numerics: $(p_g,q,K^2,c_2)=(1,1,2,10)$, $b_2=12$, $h^{1,1}=10$, and
$H^2=\mathbf Q(-1)^7\oplus H_{\rm par}$ with $\dim H^1(E,j_*W)=5$.  The $C_6$-isotypic
decomposition of $H_{\rm par}$ is $(2,0,1,1,1,0)$, i.e. $2_{\rm triv}+1_{\rm sgn}+2_{\Phi_3}$
over $\mathbf Q$; Hodge type $(1,3,1)$ then forces three of the five directions to be
$(1,1)$, so $\rho(X)=10=h^{1,1}$ (Picard-maximal) and $\operatorname{rank}T(X)=2$ —
correcting the earlier guesses $\rho=7$ and $T(X)=H_{\rm par}$.  §§14–17 attempt an
explicit branch divisor via a Belyi cover and a $j=0$ specialisation $\tau^2=-1/3$ of
Lyons's family, and honestly label it a candidate, not a theorem.

## 3.2 Checkable claims

1. $M_{0,\mathbf Z},M_{1,\mathbf Z},M_{\infty,\mathbf Z}$ preserve $J_{\mathbf Z}$, $\det J_{\mathbf Z}=1$.
2. Local spectra: $T_0=(\omega J_2,\omega^{-1}J_2)$; $T_1=(-1,-1,1,1)$; $T_\infty=(-J_2,\zeta_6,\zeta_6^{-1})$.
3. A braid relation among the three matrices (the note states none explicitly).
4. $N_0=M_0^3-I$ square-zero of rank 2; $S_0=J_{\mathbf Z}N_0$ psd of rank 2 with $\mathrm{SNF}=\mathrm{diag}(1,3,0,0)$; induced rank-2 form $\cong A_2$.
5. $N_\infty=M_\infty^6-I$ equals the printed matrix, rank 1, $\mathrm{SNF}=\mathrm{diag}(4,0,0,0)$.
6. $E:v^2=u^3-1$, $z=u^3$ is degree 6 with ramification $3,3$ over $0$; $2,2,2$ over $1$; $6$ over $\infty$; $g(E)=1$, $j(E)=0$.
7. $\deg_{\rm par}F^{2,0}_V=1/6$ (Fougeron), hence $\lambda=1$.
8. $\delta_0=10$, $\delta_1=0$, $10\lambda=\delta_0+2\delta_1$.
9. $(p_g,q,K^2,c_2)=(1,1,2,10)$; $b_2=12$; $h^{1,1}=10$.
10. $\dim H^1(E,j_*W)=5$; $C_6$-dimensions $(\dim H_0,\ldots,\dim H_5)=(2,0,1,1,1,0)$.
11. $\rho(X)=10$; $H^1(E,j_*W)\simeq T(X)\oplus\mathbf Q(-1)^3$; $\operatorname{rank}T(X)=2$.
12. Belyi $\phi(x)=-x^6/(2x^3+1)$ has cycle types $(6),(2,2,2),(3,1,1,1)$.
13. $C:y^6=-(2x^3+1)$, $u=x^2/y^2$, $v=(x^3+1)/y^3$ satisfies $v^2=u^3-1$.
14. $A^2+192\tau^6=(3\tau^2+1)(3\tau^6+75\tau^4-15\tau^2+1)$, $A=1-6\tau^2-3\tau^4$; $\tau^2=-1/3\Rightarrow j=0$.
15. Branch coefficients (16.1).

## 3.3 Verification

| # | verdict | evidence (`06_c3_*`) |
|---|---|---|
| 1 | **verified** exactly: $m^{\mathsf T}J_{\mathbf Z}m=J_{\mathbf Z}$ for all three, $\det=1$ each, $\det J_{\mathbf Z}=1$ | |
| 2 | **verified** exactly: charpolys $(x^2{+}x{+}1)^2$ (minpoly equal $\Rightarrow$ two $J_2$'s), $(x{-}1)^2(x{+}1)^2$ with $M_1^2=I$ (semisimple), $(x{+}1)^2(x^2{-}x{+}1)$ with the $(-1)$ in one $J_2$ | |
| 3 | **corrected** | No ordering gives $M_0M_1M_\infty=I$; the relation that holds exactly is $M_\infty=M_0M_1$, i.e. the printed $M_{\infty}$ is the *inverse* of the monodromy at $\infty$ in the $\gamma_0\gamma_1\gamma_\infty=1$ convention. Harmless (the eigenvalue multiset is inversion-stable and claim 10 is unchanged in either convention), but the note never states its convention. |
| 4 | **verified** | $N_0^2=0$, rank 2; $S_0=J_{\mathbf Z}N_0$ symmetric with charpoly $x^2(x^2-138x+1122)$, so **positive** semidefinite, $\mathrm{SNF}=\mathrm{diag}(1,3,0,0)$; induced form $\begin{psmallmatrix}14&19\\19&26\end{psmallmatrix}$, $\det3$, content $2$, LLL-equivalent to $\begin{psmallmatrix}2&1\\1&2\end{psmallmatrix}\cong A_2$, 6 minimal vectors of norm 2. ($N_0^{\mathsf T}J_{\mathbf Z}=-S_0$: the note's sign convention is the psd one.) |
| 5 | **verified** exactly, matrix and $\mathrm{SNF}=\mathrm{diag}(4,0,0,0)$, rank 1, square-zero | |
| 6 | **verified** exactly; Riemann–Hurwitz $2g-2=-12+(4+3+5)=0$; $c_4=0$, $c_6=864$, $\Delta=-432$, $j=0$ | |
| 7 | **not verified** (Fougeron's formula was not re-derived). This is the single unchecked numerical input; $\delta_1=0$ and $\rho=10$ both rest on $\lambda=1$. And see **D1** below: the twist that repairs $V\simeq\Lambda^2_0W$ must be tracked through this degree computation. |
| 8 | **verified** as an arithmetic consequence; but see **D7** — §7 presents it as an independent check and it is not one ($\delta_1$ is *derived* from the relation). |
| 9 | **verified**, all seven identities: $\chi=\lambda+(g_B{-}1)(g_F{-}1)=1$; $e=4(g_B{-}1)(g_F{-}1)+\delta=10$; $K^2=12\chi-c_2=2$; $p_g=\chi+q-1=1$; $b_1=2$; $b_2=e-2+2b_1=12$; $h^{1,1}=b_2-2p_g=10$; $10\lambda=10=\delta_0+2\delta_1$. |
| 10 | **verified** from the actual integral matrices over $\mathbf Q(\zeta_6)$. The 18 local kernel dimensions are, for $k=0..5$: at $0$ $(0,1,1,0,1,1)$, at $1$ $(2,2,2,2,2,2)$, at $\infty$ $(0,1,0,1,0,1)$, giving $\dim H^1_{\rm par}=4-\sum=(2,0,1,1,1,0)$, total $5$. The same in the inverse convention. Bonus: the algebra generated by the three matrices has $\mathbf Q$-dimension 16, so $W$ is **absolutely irreducible** (§8's assumption). |
| 11 | **argument is sound, conclusion accepted** modulo $\lambda=1$. The step "exactly one rank-2 rational piece contains $H^{2,0}$" is correct (conjugation swaps $H_k\leftrightarrow H_{-k}$; $H_3$ is self-conjugate of dimension 1 so cannot hold both $H^{2,0}$ and $H^{0,2}$), and $\rho=10$ is insensitive to which. |
| 12 | **verified** exactly: $\phi^{-1}(0)=(6)$; $\phi-1=-(x^3+1)^2/(2x^3+1)$ gives $(2,2,2)$; $\phi^{-1}(\infty)$ gives $(3,1,1,1)$ (three roots of $2x^3+1$ plus $x=\infty$ with order 3); $g=0$. **Stronger than the note claims:** reducing the actual matrices mod 2 and letting $\mathrm{Sp}_4(\mathbf F_2)\cong S_6$ act on the six odd theta characteristics gives $T_0\mapsto(6)$, $T_1\mapsto(2,2,2)$, $T_\infty\mapsto(3,1,1,1)$ *in exactly that order*. |
| 13 | **verified** exactly: $v^2-(u^3-1)=\bigl((x^3{+}1)^2-x^6+y^6\bigr)/y^6$ and $(x^3{+}1)^2-x^6+y^6=y^6+2x^3+1$, the relation polynomial itself. |
| 14 | **verified** exactly: both sides $=9\tau^8+228\tau^6+30\tau^4-12\tau^2+1$; $c_4(\widehat E)=16(A^2+12B)$ exactly; $\Delta(\widehat E)=256B(A^2-4B)^2$ factors as $\tau^6(\tau{-}1)^6(\tau{+}1)^6(3\tau{-}1)^2(3\tau{+}1)^2$, whose zero set is exactly Lyons's excluded set $\{0,\pm1,\pm\tfrac13\}$ — an independent corroboration. At $\tau^2=-1/3$: $c_4=0$, $\Delta\ne0$, $j(\widehat E)=0$. |
| 15 | **not checkable** from the note (quoted from Lyons with an unspecified common scalar removed). |

## 3.4 Defects and gaps in N3

- **D1 (serious, and exact).** *"$V\simeq\Lambda^2_0W$" is false as stated.* The
  $\alpha,\beta$ lists of §6 do match the stated ${}_5F_4$ (exponents at $0$ are
  $\{1-b_j\}=(0,0,0,\tfrac13,\tfrac23)=\alpha$, at $\infty$ are $\{a_i\}=\beta$; disjoint,
  $\sum\beta-\sum\alpha=\tfrac32$ so $T_1$ is a reflection; interlacing gives Hodge numbers
  $(1,3,1)$ as claimed).  But $\Lambda^2_0$ of *any* rank-4 symplectic system has local
  eigenvalues $\{1,ts,t/s,s/t,1/(ts)\}$ and therefore carries the eigenvalue $1$ at every
  puncture, while $\beta$ contains no $0$.  Concretely with §2's $W$: at $0$,
  $\Lambda^2_0W$ has exponents $(0,0,0,\tfrac13,\tfrac23)=\alpha$ ✓, but at $\infty$ it has
  $(0,\tfrac13,\tfrac13,\tfrac23,\tfrac23)\ne\beta$ ✗ and at $1$ it is not a reflection ✗.
  **Correction:** $V\simeq\Lambda^2_0W\otimes\chi$ with $\chi$ the rank-1 system with local
  monodromy $(1,-1,-1)$ at $(0,1,\infty)$, i.e. the $(1-z)^{1/2}$ character; with that twist
  all three local exponent lists match.  Hodge numbers are unchanged and $\chi$ dies on $E$
  (ramification 2 and 6 at $1,\infty$), so §§2–5 survive — but §6's Fougeron degree is
  computed *for $V$*, and $\chi$ pulls back to a possibly nontrivial 2-torsion character of
  $\pi_1(E)$, which can shift $\deg F^{2,0}_{V|E}$.  **$\lambda=1$ is the one input the note
  does not verify, and it is now the one input the correction can move.**
- **D2.** The note does **not** determine which rank-2 rational piece carries $H^{2,0}$:
  either $T(X)=H_{\rm triv}$ ($k=0$, the untwisted parabolic cohomology, already living on
  $\mathbf P^1$) or $T(X)=H_{\Phi_3}$ ($k=2,4$, visible only after the sixfold cover).
  Theorem B is insensitive; §13/§19 are not — in the second case the rank-2 transcendental
  core is *not* the Apéry object over $\mathbf P^1$ at all, and the chain
  $L(2,\chi_{-3})\leftarrow$ Apéry extension $\leftarrow V_5$ would land on a Tate piece.
- **D3.** §13/§19 conflate a rank-5 *local system* with a 5-dimensional *cohomology group*;
  they share dimension 5 and Hodge numbers $(1,3,1)$ but the note gives no map between
  them.  The headline principle ("effective arithmetic rank $<$ ODE order") rests on that
  unproved identification.
- **D5.** At $\tau^2=-1/3$, $j(\widehat E_\tau)=0$ but $j(E_\tau)=54000$ (CM by disc $-12$).
  So §15 needs Lyons's Albanese curve to be the 2-isogenous quotient $\widehat E_\tau$, not
  $E_\tau$; that is asserted, not shown.  Also $\tau=\pm i/\sqrt3$ is irrational, so (16.1)
  lives over $\mathbf Q(\sqrt{-3})$, unremarked.
- **D7/D8.** §7's "consistency check" is a derivation, not a check; §3's monodromy
  convention is $M_\infty=M_0M_1$ and is never stated.

## 3.5 What to absorb from N3

- **The dimension vector $(2,0,1,1,1,0)$ and the splitting $H^1(E,j_*W)\simeq T(X)\oplus\mathbf Q(-1)^3$.**
  Independently reproduced here from the integral matrices.  Read as a principle:
  *on a cyclic cover, the parabolic cohomology of an Apéry local system splits by character,
  and the Tate part can be large.*  The obvious follow-up for us: run the identical
  $C_w$-isotypic count on the **level-8 Catalan host** ($I_8I_2I_1I_1$, `HERFURTNER_CLASSIFICATION.md`)
  and on the level-16 cover, and see whether the Catalan $\chi_{-4}$ class survives in a
  rank-2 non-Tate quotient.  If it does, `CATALAN_OBSTRUCTION.md` (2.4) ("$c\le\#\text{cusps}-1$")
  gets a Hodge-theoretic refinement: the number of usable conditional generators is bounded
  by the rank of the non-Tate part, not by the cusp count.
- **The mod-2 $S_6$ cycle-type computation** as a cheap invariant separating hosts.
- **Nothing here bears on the Catalan entry number.** N3 does not evade any of
  (2.1)–(2.6) or supply any of (E1)–(E4); it is a structure theorem for a *different*
  conductor.  Its value to us is methodological.

---

# 4. Scoring against the obstruction framework

**Convention** (as in `CATALAN_OBSTRUCTION.md` §2.2, `CDT_FINDER.md`, and the width law
of `paper/catalan_cdt/main.tex` §7):
$$\text{entry}\ =\ \log|\varphi'(0)|-\tau(\mathbf b;\mathbf e)\ \le\ \tfrac1w\log16+\log|t_2|-k,$$
computed in the coordinate in which the *primary* row is integral; $w$ is the width of the
fold cusp, $t_2$ the surviving singularity after the hypothesis removes the fold, and $k$
the total lcm weight $\sigma$ of the conditional function (so $k=2$ for type $[1..n]^2$,
$k=4$ for $[1..2n]^2$).  Entry $>0$ is *necessary* for CDT to say anything at all
(`CDT_UNPACKED.md` §1); a conditional function is *necessary* for the hypothesis to enter
at all (`CDT_UNPACKED.md` §2).  Benchmarks: CDT's own $L(2,\chi_{-3})$ host $+0.77$;
Zagier $\mathbf E$ / level 8 $-0.61$ ($-0.077$ symmetrised, the best Catalan host known);
$\mu_4$ $-0.563$; level 16 $-1.46$; Hadamard two-host $-11$.

| construction (note) | host: singular points | local exponents | fold? | $\tau$ from measured denominators | ceiling $\log(16^{1/w}|t_2|)$ | **entry** | verdict |
|---|---|---|---|---|---|---|---|
| **cellular row** $C_n$ (N1 §11, N2 §5) | $\{0,\ \varphi^{-5}/16,\ -\varphi^{5}/16,\ \infty\}$, $t_1t_2=-\tfrac1{256}$; no apparent singularities | operator $\theta$-order 6, leading coeff $80x^6(1-176x-256x^2)$; two power-series solutions at $0$ | **yes**, $t_1$ | $[1..2n]^2$, $\sigma=4$ (measured $3.55$ at $n=34$, ratio $0.907\nearrow$) | $5\log\varphi=2.406$ | $\mathbf{-1.594}$ ($-1.14$ at measured $\sigma$; $\mathbf{-4}$ normalised over $\mathbf Q(\sqrt5)$) | negates **(E1)** (irrational $t_2$, $|t_2|=0.693>\tfrac14$) but **falls to (2.2)**: the conductor-4 twist that buys $t_2$ also doubles the lcm width, $+1.02$ against $-2$ |
| **moving EMN** $H=z\Phi$ (N1 §5, N2 §8) | $x=z/2$: $\{0,\tfrac12,1,\infty\}$ | $(0,\tfrac12)$ at $0$ and at $x=1$; $\{0,1\}$ + log at $x=\tfrac12$ (inhomogeneous) | **no** — the $\log$ coefficient is $\tfrac12$, $G$-free | $n[1..2n]$, $\sigma\le2$; measured $1.45$–$1.51$; $\tau^\flat=\tfrac32$ | $\log16=2.773$ *if the fold existed* | $\mathbf{+1.27}$ **counterfactual** | evades **(2.1)–(2.6) entirely** on budget — and then fails at the *prior* step: no conditional function (**F3**).  Not a CDT host, for a new reason |
| **signature-4 moments** $I_n$ (N2 §2) | $t\in\{0,1,\infty\}$, rank 3 | $\mathrm{Sym}^2$ of a signature-4 system | **no** — $G/\pi^2$ is a boundary value; char. roots $(1,1)$, $I_n\asymp(\log n)^2/n$ | n/a (no geometric decay) | — | **n/a** | not an Apéry system; a period relation.  Fails (F3) and has no recessive solution |
| **Gauss-parameter** $D(64x)$ (N1 §17, N2 §3) | $x\in\{0,\tfrac1{64},\infty\}$ | $(0,0)$ + log at $\tfrac1{64}$ | **no** — $16G$ is a regular connection entry, not a residue | $[1..n]^1$, $\sigma=1$ **(verified, $n\le40$)** | $\log(16/64)=-1.386$ | $\mathbf{-2.14}$ | the $64^n$ scaling costs $6\log2$; this is the paper's "conductor-4 hypergeometric realisations carry $16^n$ denominators, costing exactly $4\log2$", here worse.  Falls to (2.2) |
| **EMN radial Padé** $L_n$ (N2 §12) | not a holonomic host; a lattice construction | — | **no** | $\operatorname{den}(256^na_n)$ measured $15.886n$ (bound $L_{12n}^2=24n$) | archimedean gain $\log\tfrac{27}{16}=0.523$ | **n/a**; net $+15.36$/$n$, worthiness $0.033$ | "free integration" is real and worthless; falls to the same lcm cost as everything else |
| **EMN bare central** $I(4n,2n)$ (N2 §11) | as above | — | **no** | $\operatorname{den}(a_n)$ measured $10.483n$, of which $6\log2$ is the $2$-part | $\log64=4.159$ | **n/a**; net $+6.32$/$n$, worthiness $0.397$ | the $2$-adic slope $6$ **is** $\log_264$: the same gain twice — (2.3)/(F5) |
| **conductor-3 surface** (N3) | not a Catalan object | — | — | — | — | — | evades nothing; methodological only |

**Summary against (2.1)–(2.6) / (E1)–(E4).**

- **(2.1) Hodge depth $k=2$.**  Untouched by all three notes.  The one object with $k=1$
  (the Gauss-parameter derivative $D$) carries $G$ as a connection coefficient, not a
  period of the extension at a fold — precisely the escape that (2.1) says is impossible
  for a $k=1$ *realisation of the Apéry pair*.  N2 §3 says this correctly.
- **(2.2) Geometry, $16^{1/w}|t_2|\le4$.**  The cellular row is the **first host in the
  ledger to break the numerical bound $|t_2|\le\tfrac14$** — it has $|t_2|=0.693$ — and
  the width law still defeats it because $k$ rose from $2$ to $4$.  Net: $16^{1/w}|t_2|=11.09$
  against $e^k=e^4=54.6$.  The invariant that is actually pinned is not $|t_2|$ alone but
  $\log(16^{1/w}|t_2|)-k$, and no note moves it above $-1.1$.
- **(2.3) $p$-adic rescue blocked by Calegari.**  Untouched.  N2 §22's $\sigma_2\approx6$,
  $\sigma_3=0$ measurement is about the *rational approximants*' $p$-adic convergence, not
  about the conditional function's slope, so it does not bear on (2.3); the graded twist of
  `ADELIC_HOLONOMY.md` §2.6 needs $(\ast_p)$ for **every** $f_i$ including the conditional
  one, and that is what Calegari's irrationality of $\zeta_2(2)$ forbids.
- **(2.4) $c\le\#\text{cusps}-1$.**  N3 suggests a refinement (rank of the non-Tate
  quotient); untested for Catalan.
- **(2.5) Galois-orbit obstruction to a pure module.**  $\{t_1,t_2\}$ is a single
  $\mathbf Q(\sqrt5)$-orbit.  The *conditional* function is unaffected — its regularity at
  $t_1$ is the vanishing of a transcendental period and is not Galois-equivariant, which is
  Beukers' mechanism.  But the *pure* module is: a $\mathbf Q$-rational holonomic function
  whose singular set is defined over $\mathbf Q$ cannot be singular at $t_2$ and regular at
  $t_1$, so over $\mathbf Q$ there is no pure module on $\mathbf P^1\setminus\{0,t_2,\infty\}$
  and $\varphi$ must avoid $t_1$ as well — the ceiling collapses.  Over $\mathbf Q(\sqrt5)$
  one needs the number-field holonomy bound that `CDT_FINDER.md` §7 calls "the single
  load-bearing unproved input", and the norm accounting gives
  $\log16+\tfrac12\log|N(t_2)|=\log16+\tfrac12\log\tfrac1{256}=0$, hence
  $$\text{entry}_{\mathbf Q(\sqrt5)}=0-4=-4.$$
  So $-1.594$ is the *optimistic* single-place figure and $-4$ the honest normalised one.
  **The cellular row inherits (2.5) as well as (2.2)** — and this is exactly why
  `CDT_FINDER.md` §6 records $|N(t_2)|=1$ as the Apéry-perfect condition: the cellular row
  breaks it by the factor $256$.
- **(2.6) Cuspidal disguise.**  Untouched.
- **(E1)** partially realised and shown insufficient (above).  **(E2)** untouched.
  **(E3)** — a non-Eisenstein motive with period $G$ — is exactly what EMN supplies, and
  N1/N2 are its first serious exploitation; the verdict is that the *budget* becomes
  favourable and the *mechanism* disappears.  **(E4)** untouched.

---

# 5. What to absorb, and what is wrong

## 5.1 Genuinely new and worth absorbing

1. **The cyclotomic cellular Catalan Apéry system** (N1 §11, N2 §5).  Fully verified here:
   $C_n=A_nG+B_n$ exactly, $\pi$-free, $U_n\in\mathbf Z$, one order-2 recurrence for both
   solutions, $\operatorname{den}(V_n)\mid d_{2n}^2$, $\lambda_\pm=88\pm40\sqrt5$.  It is a
   *new row for the ledger* (it appears nowhere in `LEDGER_DUMP.md`, `ROW_LEDGER.md`,
   `SPORADIC_SCAN2.md`, `CATALAN_*`).  **The structural discovery that came out of checking
   it is the important part:** its singular points are $\varphi^{\mp5}/16$, i.e. **Zagier's
   sporadic row $\mathbf D$ / the $X_1(5)$ golden host divided by $16$**, with
   $t_1t_2=-\tfrac1{256}$ instead of the Apéry-perfect $-1$.  So the first genuinely
   non-modular Catalan row we have found lands on a host we had already priced for
   $L(3,\chi_5)$, and the whole cost of carrying $\chi_{-4}$ on it is the factor $16$ —
   paid twice, once in the coordinate and once in the doubled lcm width.  That is a clean
   quantitative statement of "the conductor-4 twist costs $4\log2$" and it should go into
   `CATALAN_OBSTRUCTION.md` §2.2 as the sharpened form of (E1).
2. **The exact reduction $\iint_{[0,1]^2}\frac{x^Ay^B}{(1+x^2y^2)^r}$ to $J^{(r)},K^{(r)}$**
   (N2 §4), with the two closed recurrences.  This is a cheap, exact engine for producing
   Catalan cellular rows in the basis $1,\pi,\log2,G,\pi^2$, and the congruence conditions
   on $(A,B)$ that force the unwanted coordinates to vanish are explicit.  Worth keeping as
   `lattice/sol_notes/01_cellular.gp`: it generates the whole $\mu_4$-cellular family, of
   which $C_n$ is one member, and a scan over numerator shapes is now a finite computation.
3. **$\int_0^1{}_2F_1(\tfrac14,\tfrac34;1;t)^2dt=16G/\pi^2$** with the exact
   $\mathrm{Sym}^2$ operator and the inhomogeneous moment recurrence with inhomogeneity
   $8/\pi^2$ (N2 §2).  Verified to 40 digits.  This puts $G/\pi^2$ *and* $1/\pi^2$ in one
   rank-3 signature-4 period system — a genuinely different Hodge placement for $G$ from
   the weight-3 Eisenstein one, and the only place in the archive where $G/\pi^2$ and
   $1/\pi^2$ appear as two moments of the same system.  It is not an Apéry system (double
   characteristic root $1$), but it is the natural home for a *three-period* Catalan
   statement of CDT type ($1,\pi^2,G$), which `CATALAN_THREE_PERIOD.md` wanted and could
   not find on the modular side.
4. **$d_n[x^n]D(64x)\in\mathbf Z$** where $D=\partial_a\,{}_2F_1(a,1-a;1;t)|_{a=1/4}$
   (N1 §17, N2 §3).  Verified for every $n\le40$, and sharply: the minimal $k$ with
   $\operatorname{den}\mid\mathrm{lcm}(1..k)$ is $\le n$ throughout.  This is the **only
   Catalan-carrying object in the archive with a single lcm layer**.  The mechanism is
   visible: $\binom{4n}{2n}\binom{2n}n$ contains each prime of $(n,4n)$ exactly once, which
   cancels the odd denominators $4k\pm1$ above $n$, and supplies the excess valuation at the
   small primes.  Worth turning into a proved lemma.  Record it as a
   theorem-shaped statement — it is the sharpest available counterexample to "Catalan
   always costs two lcm layers", *at the cost of $G$ appearing as a connection coefficient
   rather than a fold period*.
5. **The reframing of the Catalan obstruction supplied by the moving EMN period.**  For the
   first time we have a Catalan host that passes the budget ($+1.27$) and fails the
   mechanism: $H(z)=G+\tfrac{1-z}2\log(1-z)+O(1-z)$ with a $G$-free branch coefficient.
   The ledger sentence should become: *on every Catalan host constructed, either the
   conformal size is too small for the Hodge depth, or $G$ is a regular value and there is
   no fold.*  Every construction in these three notes, plus $\mu_4$, is on one side or the
   other of that dichotomy.  Finding a host on neither side is the target.
6. **(N3) The $C_6$-isotypic splitting of parabolic cohomology,
   $H^1(E,j_*W)\simeq T(X)\oplus\mathbf Q(-1)^3$ with dimension vector $(2,0,1,1,1,0)$**,
   reproduced here from the integral monodromy.  Method to import: on a cyclic cover, count
   $\dim H^1_{\rm par}(\mathbf P^1,W\otimes L_{\chi_k})=(s-2)\,\mathrm{rk}\,W-\sum_i\dim\ker(T_i\chi_k-1)$
   and split the Tate part off.  Applied to the level-8/level-16 Catalan hosts this would
   give a Hodge-theoretic refinement of (2.4) ($c\le\#\text{cusps}-1$): the number of usable
   conditional generators should be bounded by the rank of the non-Tate quotient.
7. **The closed form of the EMN period family** (N2 §§10–13, and the exact engine built
   here): every $I(m,t)=\iint_\Delta\frac{x^my^m}{(1-x^2-y^2)^{t+1}}$ lies in $\mathbf Q+\mathbf QG$
   because *all* of its $G$-content is the single integral $\int_0^{\pi/4}\log\tan u\,du=-G$;
   the higher harmonics $A_k=-\frac1{2k}\sum_{i<k}\frac{(-1)^i}{2i+1}$ are rational and the
   $\pi$-parts cancel identically.  That is a clean *structural* proof of the period purity
   that the note asserts, and it makes the whole EMN Padé problem exact finite-dimensional
   linear algebra over $\mathbf Q$ (`05_emn_core.py`).
8. **The identity $\sigma_2=6=\log_264$ for the EMN row.**  The measured $2$-adic slope of the
   central family is *numerically equal* to its archimedean saddle exponent.  This is the
   sharpest instance we have of (F5)/(2.3): a large-looking $p$-adic slope that is the
   archimedean gain re-expressed, hence worth nothing in the net.  Add it to
   `CATALAN_OBSTRUCTION.md` (2.3) as a worked example.
9. **(N3) Cheap host invariants:** the mod-2 $\mathrm{Sp}_4(\mathbf F_2)\cong S_6$ cycle
   types of the local monodromies, and the thickness reading of $\mathrm{SNF}(M^e-I)$.

## 5.2 Wrong, unsupported, or misleading

1. **N1 §2: "Viola and Marcovecchio obtain effective approximants with exponent
   approximately $0.6293$."**  No such exponent is in any source we hold; $0.62922$ is
   CDT's contour-loss factor $|\psi'(0)|$ (`CDT_FINDER.md` line 94).  Delete.
2. **N2 §9: the quantitative no-go for the moving EMN period is invalid.**  "Conformal
   radius at most $2$" is the univalent Koebe bound; CDT's method is multivalent and the
   applicable ceiling is $16|t_2|$.  With the fold removed the entry would be $+1.27$, not
   $-0.81$.  The *conclusion* (not a viable CDT host) survives, but for the structural
   reason in the same section, and that reason should be stated as the theorem.
3. **N2 §2–3 omits that the signature-4 moment recurrence has a double characteristic root
   at $1$** ($8\lambda^2-16\lambda+8$), so $I_n\asymp(\log n)^2/n$ and there is no recessive
   solution.  The section reads as though a Diophantine row had been produced; it has not.
4. **N1's "Catalan Realization Hypothesis" and Conjectures A–E are untested here**, and
   Conjecture A (the rank-2 quotient extends to a nonconstant rank-2 variation) is
   contradicted in spirit by the verified computation: the natural deformation $\Phi(z)$ has
   a rank-*three* Gauss–Manin system ($1$, $\arcsin$, $H$), with $G$ at a point where the
   monodromy is $G$-independent.  N1's own stopping condition in Stage 1 fires.
5. **N1 §14's hope that $m=2$** ("the enormous potential advantage … is that one may have
   $m=2$") is a misreading of CDT: small $m$ makes the *bound* easy to satisfy, not hard;
   the contradiction requires *exhibiting more functions than the bound allows*, so a small
   period module is a liability, not an asset.  This should be corrected before any of the
   work packages is executed — it inverts the objective of Work Package H.
6. **N2 §10: the Catalan coefficient is missing a sign.**  The correct closed form is
   $b_{m,t}=(-1)^t4^{-m}\binom m{m/2}\binom mt$; the note's version is off by $-1$ at every
   odd $t$ (all 10 odd-$t$ cases with $m\le8$).  Harmless for the radial family ($t$ even
   there) but not in general.
7. **N2 §22's $v_3=O(1)$ is $O(\log n)$** ($\approx-1.95\log n$ envelope); the substantive
   claim (no $3$-adic slope) survives.  And N2 §10's pole-lowering recurrence **does not
   close on the central family** — applied twice from $(4n,2n)$ it lands on $(4n-4,2n-4)$ —
   so there is no $O(1)$ recurrence for $I_n$ by that route, contrary to the impression
   §11 gives.
8. **N2 §13's "the archimedean side has been solved" is true but misleading.**  It is
   solved by $0.523$ nats/$n$ against a measured companion cost of $15.886$ nats/$n$.  The
   sentence should read: the archimedean side is no longer the binding constraint, and the
   binding constraint is $30\times$ larger than the slack.
9. **N3 §1: "$V\simeq\Lambda^2_0W$" is false as stated**; the correct statement needs the
   quadratic twist $\chi=(1-z)^{1/2}$ (§3.4 D1).  Everything downstream survives except
   possibly $\deg_{\rm par}F^{2,0}=1/6$, hence $\lambda=1$, hence $\delta_1=0$ and
   Theorem B.
10. **N3 §7's "consistency check" is a derivation**, and **N3 §3's monodromy convention is
   $M_\infty=M_0M_1$** (never stated).  **N3 §15** needs Lyons's Albanese to be the
   $2$-isogenous quotient ($j(E_\tau)=54000$, $j(\widehat E_\tau)=0$), and $\tau=\pm i/\sqrt3$
   puts (16.1) over $\mathbf Q(\sqrt{-3})$.

## 5.3 Recommended next computations (in priority order)

1. **Scan the cellular family** with the exact engine of §5.1(2): numerators
   $x^{a}(1-x^4)^{b}y^{a}(1-y^4)^{b}/(1+x^2y^2)^{r}$ and their Jacobi-type relatives, for the
   *joint* objective $\log|\lambda_-|+\sigma$ (Apéry) and $\log(16|t_2|)-\sigma$ (entry).
   The one member checked has $\sigma=4$; the question is whether any member has $\sigma=2$
   while keeping $|t_2|>1/4$.  That is precisely the missing $+2$ nats, and it is now a
   finite search.
2. **Redo the $C_w$-isotypic parabolic count on the level-8 and level-16 Catalan hosts**
   (N3's method) and see whether the $\chi_{-4}$ class is in a rank-2 non-Tate quotient; if
   so, refine (2.4).
3. **Check whether $\int_0^1t^n{}_2F_1(\tfrac14,\tfrac34;1;t)^2dt$ has a *twisted* relative
   with geometric decay** — the double root at $1$ is the only thing wrong with the
   signature-4 realisation, and a $t\mapsto$ rational-pullback that separates the two roots
   would give a three-period $(1,\pi^2,G)$ Apéry system with a genuinely non-Eisenstein
   Hodge placement.
4. Prove the cellular recurrence by creative telescoping (N2's own request); it is verified
   on 35 exact terms here, 61 in the note, but not proved.

5. **Price the EMN lattice properly before any more numerator search.**  The measured
   companion cost is $15.886n$ (radial) and $10.483n$ (bare), against gains $0.523$ and
   $4.159$.  The whole search space of N2 §24 moves the archimedean term by at most
   $\approx0.3$ nats.  The only quantity worth attacking is the **odd** part of
   $\operatorname{den}(a_n)$ — measured $\approx6.32$ nats/$n$ in the bare family after the
   $2$-part is removed — and the sharp statement to prove or refute is: *is
   $\operatorname{den}(256^na_n)$ forced to contain every prime $p\le12n$ to the second
   power?*  E10 shows the minimal $k$ is exactly the largest prime $\le12n$ for all
   $n\le30$, so the answer looks like yes, and if it is a theorem the EMN Padé route is
   closed by $6$ nats and should be retired.

---

## Appendix. Scripts

All under `lattice/sol_notes/`.

| script | what it computes |
|---|---|
| `01_cellular.gp` | exact $C_n=A_nG+B_n$ via the $J^{(r)},K^{(r)}$ engine; $U_n,V_n$ to $n=34$ |
| `02_recur.gp` | recurrence residuals, $d_{2n}^2$ divisibility, characteristic roots |
| `03_moving.py`, `03b_branch.py` | moving-EMN moments, recurrence, PF equation, $H'$ closed form, $\Phi(1)=G$, branch coefficient at $z=1$ |
| `04_sig4.gp`, `04b_sig4num.py` | $\mathrm{Sym}^2$ ODE, $[t^n]F$, $D(64x)$ denominators; $I_n$ moments and the $8/\pi^2$ recurrence |
| `05_emn_core.py`, `05_emn_fast.py`, `05_emn_E*.py`, `05_emn_extra*.py` | exact $I(m,t)$ engine, E0–E12 |
| `06_c3_matrices.gp`, `06_c3_lattice_c7.gp`, `06_c3_curves.gp`, `06_c3_ellj_c9_c11.gp`, `06_c3_mod2_hgeom.py` | conductor-3 monodromy, $A_2$ form, $C_6$ parabolic count, Belyi map, Lyons factorisation |
| `07_types.gp` | denominator type of the moving-EMN period in $x=z/2$ |
| `08_score.gp` | singular points of the cellular row, $\varphi^{\pm5}/16$ identification, entry numbers |
