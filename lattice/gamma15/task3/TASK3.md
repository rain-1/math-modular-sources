# Task 3 — CDT's fourteen functions transported to Zagier's row D on $\Gamma_1(5)$

*Working directory: `lattice/gamma15/task3/`. All scripts and logs are in this directory.
Reference implementation reused: `lattice/cdt_finder/indep_check2.py`, `cdt_bound.py`;
peeling routine from `lattice/hostscan/20_foldreg.gp`. CDT = Calegari–Dimitrov–Tang,
arXiv:2408.15403v2, LaTeX at `papers/cdt/cdt2/L2chi.tex`.*

**Headline.** The measured denominator array on the $\Gamma_1(5)$ host is **identical** to
CDT's $(\ref{typesA})$: $m=14$, $r=2$, $b_1=b_2=2$, $u_1=1$, $u_2=3$,
$\mathbf e=(0,0,1;0,0,0,0,0,0;1,1,1,1,1)$, so
$$\tau^\flat=\tfrac{191}{49},\qquad \tau^\sharp=\tfrac{27}{80},\qquad
\boxed{\tau(\mathbf b;\mathbf e)=\tfrac{16603}{3920}=4.235459183673\ldots}$$
**stands unchanged.** The fourteen functions are $K(y)$-linearly independent in every degree
$\le 5$, certified at a split and at an inert prime of $K=\mathbf Q(\sqrt5)$.

---

## 0. Setting and conventions

Host: $(n+1)^2A_{n+1}=(11n^2+11n+3)A_n+n^2A_{n-1}$, $A=1,3,19,147,1251,11253,\dots$
(equivalently $A_n=\sum_k\binom nk^2\binom{n+k}k$), $F(x)=\sum A_nx^n$.
Recurrence re-verified for all $n\le 209$ against the coefficients produced by the
$q$-expansion peel **[verified, exact]**.

Singularities $t_1=\phi^{-5}$, $t_2=-\phi^5$, $t_1t_2=-1$. **Outer singularity**
$s:=t_2=-\phi^5=\frac{-11-5\sqrt5}{2}$, with
$s^{-1}=-\phi^{-5}=\frac{11-5\sqrt5}{2}$ and $s\cdot s^{-1}=1$: $s\in\mathbf Z[\phi]^\times$
**[proved, exact integer identity]**.

Descent: $w(x)=\dfrac{sx}{x-s}$, $y=x+w(x)=\dfrac{x^2}{x-s}$, branch point $y=4s$,
normalised coordinate $Y:=y/s$ (branch point $Y=4$). Verified as exact $K[[x]]$ identities to
order $42$: $x+w=y$ and $x\,w=s\,y$ **[verified]**.

Because $s$ is a **unit** of $\mathbf Z[\phi]$ and $y=sY$:

* $K(y)=K(Y)$ and the degree filtrations by $y^j$ and $Y^j$ coincide, so the $\deg\le D$
  independence statements of §3 are the same whichever variable is used;
* $\frac{d}{dy}=s^{-1}\frac{d}{dY}$ and $\int\!dy=s\int\!dY$, so derivatives and integrals
  differ from their $Y$-versions by units and neither the span nor the denominator type
  changes;
* if $D_n c_n\in\mathbf Z$ for $B_i(Y)=\sum c_nY^n$, then the *actual function of $y$* is
  $B_i(y/s)=\sum c_ns^{-n}y^n$ and $D_n\,c_n\,s^{-n}\in\mathbf Z[\phi]$.

The last point is **[proved]** ($s^{-n}\in\mathbf Z[\phi]^\times$) and was additionally
**[verified]** coefficient-by-coefficient for all seven $B_i$ and all $n\le 80$
(`t1b_zphi.py` → `t1b_zphi.log`). *All denominator types below are therefore stated in $Y$;
they are literally the same in $y$ over $\mathbf Z[\phi]$.*

---

## 1. (T1) The pure module: $B_1,\dots,B_7$ with measured types

Scripts `t1_pure.py`, `t1b_zphi.py`, `t1c_sharp.py`, `t1d_cbin.py`.
Series built exactly over $\mathbf Q$ to order $96$ (measurements at $n\le 80$):
$$B_1=1,\quad B_2=\sum_{n\ge2}\frac{2(n-2)!\,n!}{(2n)!}Y^n,\quad
B_3=\sum_{n\ge1}\frac{((n-1)!)^2}{(2n)!}Y^n,\quad
B_5=\sum_{n\ge1}\frac{((n-1)!)^2}{(2n-1)!\,(2n-1)}Y^n,$$
$$B_6=\int\frac{B_3}{Y}dY=\sum_{n\ge1}\frac{((n-1)!)^2}{n\,(2n)!}Y^n,\qquad
B_7=\int\frac{B_4}{Y}dY,$$
$$B_4=\operatorname{Sym}^-\operatorname{Li}_2 = 4\sum_{n\ge0}\frac{Y^{n+1}}{16^n}
\sum_{k=0}^n\binom{2k}k\binom{2n-2k}{n-k}\frac1{(2k-1)(2n-2k+1)^2}.$$
$B_4$'s head reproduces CDT's $(\ref{defB})$: $-4Y+\frac49Y^2+\frac{31}{900}Y^3+\frac{389}{88200}Y^4$
**[verified]**, and $B_4$ was independently recomputed to order $20$ by literally
symmetrising $\operatorname{Li}_2$ in $x$ and re-expanding in $Y$ — agreement **[verified]**.

### 1.1 CDT's Lemma (`bdenominators`): status

| $i$ | CDT claim | valid for $n\le80$? | necessary? |
|---|---|---|---|
| $B_1$ | trivial | **VALID** | — |
| $B_2$ | $[1..2n]$ | **VALID** | trivial fails |
| $B_3$ | $[1..2n]\,n$ | **VALID** | $[1..2n]$ alone fails |
| $B_4$ | $[1..2n]^2$ | **VALID** | $[1..n][1..2n]$ and $[1..2n]n^2$ both fail |
| $B_5$ | $[1..2n](2n-1)$ | **VALID** | $[1..2n]$ and $[1..2n]n^2$ both fail |
| $B_6$ | $[1..2n]\,n^2$ | **VALID** | $[1..2n]n$ and $[1..n][1..2n]$ both fail |
| $B_7$ | $[1..2n]^2\,n$ | **VALID** | $[1..2n]^2$ alone fails |

**Every clause of CDT's Lemma is confirmed; none needs correction.** **[verified, $n\le80$]**

**Measured minimal types** in the family $n^e\prod_j[1..b_jn]$ ($b_j\in\{1,2\}$, $e\le2$),
i.e. the sharp entries of the array:

| | minimal admissible types (incomparable alternatives listed) |
|---|---|
| $B_1$ | $1$ |
| $B_2$ | $[1..2n]$ |
| $B_3$ | $[1..2n]\,n$ **or** $[1..n][1..2n]$ |
| $B_4$ | $[1..2n]^2$ |
| $B_5$ | $[1..2n]^2$ |
| $B_6$ | $[1..2n]\,n^2$ **or** $[1..n][1..2n]\,n$ **or** $[1..2n]^2 n$ |
| $B_7$ | $[1..2n]^2\,n$ |

Slack: $\gcd_n\bigl(T_{\rm claimed}(n)/\operatorname{den}(c_n)\bigr)=1$ for $B_2,B_3,B_6$ and
$=2$ for $B_4,B_5,B_7$ — i.e. CDT's stated types for $B_4,B_5,B_7$ carry one spare factor $2$.
A constant factor does not change the array or $\tau$. **[verified, $n\le80$]**

### 1.2 One correction, to CDT's Remark (`central binomial`)

CDT remark that "the $[1,\dots,2n]$ can be relaxed to $n(n-1)\binom{2n}{n}$". Measured
(`t1d_cbin.py`, $2\le n\le 60$), writing $c(n)=n(n-1)\binom{2n}n$:

* $B_2:c(n)$ ✓, $B_3:c(n)\,n$ ✓, $B_5:c(n)(2n-1)$ ✓, $B_6:c(n)\,n^2$ ✓ — the remark holds
  verbatim for the four functions with an explicit factorial coefficient;
* $B_4$: $c(n)^2$ **FAILS**, and so do $c(n)^2n^2$ and $c(n)^2n^4$. What does hold is
  $c(n)\,[1..2n]$ ✓ and $c(n)^2[1..n]$ ✓.
* $B_7$: $c(n)^2n$ **FAILS**; $c(n)[1..2n]\,n$ ✓ holds.

So for the two functions of two-layer type only **one** of the two $[1..2n]$ layers relaxes to
$c(n)$; the literal reading of the Remark is too strong. CDT themselves note the Remark makes
no difference to their setup, and it makes none here either. **[verified, $2\le n\le60$]**

---

## 2. (T2) The conditional side over $K$

### 2.1 Building $H$

`t2_companions.gp` (PARI 2.15) builds, to $q$-precision $222$ and $x$-order $210$:
$x=q\prod_{n\ge1}(1-q^n)^{5(\frac n5)}=q-5q^2+15q^3-30q^4+40q^5-26q^6-\dots$;
$F=\sum A_nx^n$; the four weight-3 Eisenstein series $P_1,\dots,P_4$ for the odd quartic
character $\psi_4$ mod 5; $R_1=P_1+P_2$, $R_2=i(P_1-P_2)$, $R_3=P_3+P_4$, $R_4=i(P_3-P_4)$;
$\Phi_D=R_1/2+R_2$ and $\Phi_{\rm new}=R_3+\phi^5R_4$; and the companions
$B=\bigl(\text{$x$-expansion of }F\cdot D^{-2}\Phi\bigr)$ obtained with the `peel2` routine.

*Because peeling is $\mathbf Q$-linear in the input series, $B_{\rm new}=B_{R_3}+\phi^5B_{R_4}$
is assembled from two rational peels; no arithmetic in $K$ is needed in PARI.*

Checks:
* peel of $F$ itself returns $1,3,19,147,1251,11253$ **[verified]**;
* $B_{D,n}/A_n\to\zeta(2)/5$: at $n=210$ the two agree to all **60 printed digits**
  ($0.328986813369645287294483033329205037843789980241359687547112$) **[verified, 60 digits]**;
* $B_{{\rm new},n}/A_n\to\xi$: at $n=210$,
  $0.655634188840656766330981413872399402411138\ldots$ against the target
  $\xi=0.6556341888406567663309814138723994024111$ — agreement to all **40 digits given**
  **[verified, 40 digits]**;
* $B_D$, $B_{R_3}$, $B_{R_4}$ all have denominator type $[1..n]^2$ over $\mathbf Q$ for
  $n\le100$, and $[1..n]^1$ fails (first at $n=2,4,3$ respectively) **[verified]**.

Heads: $B_D=x+\frac{25}4x^2+\frac{1741}{36}x^3+\dots$, $B_{R_3}=2x+18x^2+166x^3+\dots$,
$B_{R_4}=-\frac12x^2-\frac{113}{18}x^3-\dots$

$H=aA+bB_D+cB_{\rm new}$ with the generic $(a,b,c)=(1,\,-3+\phi,\,5-2\phi)\in\mathbf Z[\phi]^3$.

### 2.2 The symmetrisation recurrence

With $x+w=y$ and $xw=sy$, Newton's identity gives immediately
$$P_0=2,\quad P_1=y,\quad \boxed{P_n(y)=y\,P_{n-1}(y)-s\,y\,P_{n-2}(y)}\qquad
P_n=x^n+w^n. \tag{proved}$$
Setting $Q_n(Y):=s^{-n}P_n(sY)$ the factor $s$ cancels *exactly*:
$Q_0=2$, $Q_1=Y$, $Q_n=YQ_{n-1}-YQ_{n-2}$ — **CDT's own integer polynomials verbatim**.
Hence
$$\operatorname{Sym}^+H=\sum_n h_nP_n(y)=\sum_n\bigl(h_ns^n\bigr)Q_n(Y),\qquad
[Y^N]G=\sum_{n=N}^{2N}\bigl(h_ns^n\bigr)\,[Y^N]Q_n .$$
Verified: $\deg Q_n=n$, $\operatorname{val}Q_n=\lceil n/2\rceil$, integer coefficients, for all
$n\le210$; and $x^n+w^n=s^nQ_n(Y)$ as exact $K[[x]]$ identities for $0\le n\le42$
**[verified]**. As an end-to-end check, $G(Y(x))=H(x)+H(w(x))$ was confirmed as an exact
$K[[x]]$ identity to order $42$ **[verified]**.

Independent calibration of the same code path on CDT's host ($s=1$):
$G_A(y)=2-27y+1014y^2-49536y^3+\dots$, matching CDT's Definition (`defG`) exactly **[verified]**.

On our host, $G_A(Y)=\operatorname{Sym}^+F$ begins
$$2\;-\;\tfrac{4707+2105\sqrt5}{2}\,Y\;+\;\tfrac{38451615+17196085\sqrt5}{2}\,Y^2
\;-\;\tfrac{400066987845+178915396075\sqrt5}{2}\,Y^3+\dots\;\in\;\mathbf Z[\phi][[Y]].$$

### 2.3 Denominator type of $G$ — measured

For $n\le60$, $d_{2n}^2\cdot[Y^n]G\in\mathbf Z[\phi]$, i.e. **$G$ has denominator type
$[1..2n]^2$** — CDT's Lemma (`Gdenominatortype`) transports unchanged. This was measured for
five different $(a,b,c)$:

| $(a,b,c)$ | measured minimal type | $[1..n][1..2n]$ enough? |
|---|---|---|
| $(1,-3+\phi,5-2\phi)$ (generic) | $[1..2n]^2$ | **no** |
| $(1,0,0)$ (=$G_A$) | $1$ (trivial) | yes |
| $(0,1,0)$ | $[1..2n]^2$ | **no** |
| $(0,0,1)$ | $[1..2n]^2$ | **no** |
| $(\phi,\phi^2,\phi^3)$ | $[1..2n]^2$ | **no** |

So **$[1..2n]^2$ is sharp** for $G$ whenever $b$ or $c\ne0$ **[verified, $n\le60$]**.
(The proof is CDT's: $Q_n$ contributes to $Y^N$ only for $N\le n\le 2N$, $h_ns^n$ has
$d_n^2h_ns^n\in\mathbf Z[\phi]$, so $d_{2N}^2[Y^N]G\in\mathbf Z[\phi]$. **[proved]**)

---

## 3. (T3) $K(y)$-linear independence of the fourteen

Script `t3_indep.py` (log `t3_indep.log`); $p$-integrality guard `t3b_pintegral.py`.

The fourteen, in CDT's order:
$$B_1,B_2,B_3;\ B_4,B_5,G,G',G'',G''';\ B_6,B_7,\textstyle\int G\,dy,
\int\frac{G-G(0)}{y}dy,\int\frac{G-G(0)-G'(0)y}{y^2}dy .$$

Two primes near $2^{61}$:

* **SPLIT** $p=2305843009213694009$ ($p\equiv4\bmod5$), $\sqrt5\equiv56964414095014148$,
  $\mathbf Z[\phi]/\mathfrak p=\mathbf F_p$;
* **INERT** $p=2305843009213693967$ ($p\equiv2\bmod5$), $\mathbf Z[\phi]/p=\mathbf F_{p^2}$
  implemented as $\mathbf F_p[t]/(t^2-5)$ (5 verified to be a non-residue).

Guard: no coefficient of any of the fourteen (either host, either convention) has a
denominator divisible by either prime — the reductions are legitimate **[verified]**.

### 3.1 Rank tables

Rank of $\{Y^jf_i:1\le i\le14,\ 0\le j\le D\}$ against $14(D+1)$, series order $L$:

**(b) $\Gamma_1(5)$ host, $K=\mathbf Q(\sqrt5)$, $s=-\phi^5$**

| $D$ | $L$ | rank / needed, SPLIT $\mathbf F_p$ | rank / needed, INERT $\mathbf F_{p^2}$ |
|---|---|---|---|
| 0 | 26 | **14 / 14** | **14 / 14** |
| 1 | 40 | **28 / 28** | **28 / 28** |
| 2 | 54 | **42 / 42** | **42 / 42** |
| 3 | 68 | **56 / 56** | **56 / 56** |
| 4 | 82 | **70 / 70** | **70 / 70** |
| 5 | 96 | **84 / 84** | **84 / 84** |

**(a) Cross-check — CDT's own host ($s=1$), same code path**

| $D$ | $L$ | SPLIT $\mathbf F_p$ | INERT $\mathbf F_{p^2}$ |
|---|---|---|---|
| 0 | 26 | 14 / 14 | 14 / 14 |
| 1 | 40 | 28 / 28 | 28 / 28 |
| 2 | 54 | 42 / 42 | 42 / 42 |
| 3 | 68 | 56 / 56 | 56 / 56 |
| 4 | 82 | 70 / 70 | 70 / 70 |
| 5 | 96 | 84 / 84 | 84 / 84 |

which is exactly the table of `consolidation/CDT_FINDER.md` §2 (14/14, 28/28, 42/42, 56/56,
70/70, 84/84) **[reproduced]**. Both tables are unchanged under either of two readings of
items 13–14 (see §3.3).

### 3.2 What is certified

A full-rank reduction mod a prime of $K$ certifies full rank over $K$. Therefore:

> **[verified]** There is **no relation $\sum_{i=1}^{14}P_i(y)f_i(y)=0$ with
> $P_i\in K[y]$, $\deg P_i\le5$, not all zero** — at *two independent* reductions, one
> where $\mathbf Z[\phi]/\mathfrak p=\mathbf F_p$ (split) and one where
> $\mathbf Z[\phi]/p=\mathbf F_{p^2}$ (inert).

This is a degree-bounded statement. Full $K(y)$-linear independence (all degrees) is **not**
established by this computation; over $\mathbf Q$ CDT prove it by monodromy (their
Lemma `14functions`), and that argument is expected to transport, but that transport is
**[assumed]** here, exactly as in `CDT_FINDER.md` §2.

### 3.3 A coding note on `indep_check2.py`

In `lattice/cdt_finder/indep_check2.py` items 13 and 14 are coded as
`i_([0]+[G[n]/n])` and `i_([0,0]+[G[n]/(n-1)])`. The inner list of item 13 already *is*
$\int\frac{G-G(0)}y dy=\sum_{n\ge1}\frac{G_n}{n}y^n$, so the outer `i_` applies one extra
integration; item 14 is likewise index-shifted relative to
$\int\frac{G-G(0)-G'(0)y}{y^2}dy=\sum_{n\ge1}\frac{G_{n+1}}{n}y^n$. Both CDT's definitions and
`indep_check2.py`'s coding were run here; **both give full rank at every $D\le5$ on both
hosts**, so nothing downstream is affected. Recorded for the record.

---

## 4. (T4) The array and $\tau$

Scripts `t4_array.py` (measurement of all fourteen), `t4_exact.py` (exact $\tau$; calls
`lattice/cdt_finder/cdt_bound.py`).

### 4.1 Measured array on the $\Gamma_1(5)$ host, $n\le50$

Types are stated in CDT's overflow-relaxed family $n^{e}\prod_j[1..b_jn+c_j]$
(their Remark `rem:overflow`), which is needed because differentiating and dividing by powers
of $y$ shifts the index.

| # | function | measured minimal type | layers | $e$ |
|---|---|---|---|---|
| 1 | $B_1$ | $1$ | **0** | 0 |
| 2 | $B_2$ | $[1..2n]$ | **1** | 0 |
| 3 | $B_3$ | $[1..2n]\,n$ | **1** | **1** |
| 4 | $B_4$ | $[1..2n]^2$ | 2 | 0 |
| 5 | $B_5$ | $[1..2n]^2$ (sharper: $[1..2n](2n-1)$) | 2 | 0 |
| 6 | $G$ | $[1..2n]^2$ | 2 | 0 |
| 7 | $G'$ | $[1..2n+1]^2$ | 2 | 0 |
| 8 | $G''$ | $[1..2n+3]^2$ | 2 | 0 |
| 9 | $G'''$ | $[1..2n+5]^2$ | 2 | 0 |
| 10 | $B_6$ | $[1..2n]\,n^2$; also $[1..2n]^2n$ | 2 | **1** |
| 11 | $B_7$ | $[1..2n]^2\,n$ | 2 | **1** |
| 12 | $\int G\,dy$ | $[1..2n]^2\,n$ | 2 | **1** |
| 13 | $\int\frac{G-G_0}{y}dy$ | $[1..2n]^2\,n$ | 2 | **1** |
| 14 | $\int\frac{G-G_0-G_1y}{y^2}dy$ | $[1..2n+2]^2\,n$ | 2 | **1** |

Necessity, measured: for #4–#9 no single-layer type works even after relaxing $c$ up to $6$;
for #10–#14 the factor $n$ is necessary ($[1..2n+6]^2$ alone fails). **[verified, $n\le50$]**

Hence

$$u_1=1\ \ (\text{only }B_1\text{ carries no layer}),\qquad
u_2=3\ \ (B_1,B_2,B_3\text{ carry at most one layer}),\qquad b_1=b_2=2,$$
$$\mathbf e=(0,0,1;\,0,0,0,0,0,0;\,1,1,1,1,1),\qquad \textstyle\sum e_i=6,\quad\max e_i=1 .$$

**The array is IDENTICAL to CDT's $(\ref{typesA})$.**

*The one genuine alternative* is to move $B_6$ into the one-layer column with its own type
$[1..2n]n^2$, giving $u_2=4$, $\sum e=7$, $\max e=2$. Then $\tau^\flat=375/98=3.8265\ldots$
improves but $\tau^\sharp=0.5401785714\ldots$ worsens, and the total is
$\tau=4.366709\ldots$, i.e. **worse by $0.13125$**. CDT's placement is optimal.
**[verified]**

### 4.2 Exact $\tau$

$$\sigma_m=b_1+b_2=4,$$
$$\tau^\flat(\mathbf b)=\sigma_m-\frac{u_1^2b_1+u_2^2b_2}{m^2}
=4-\frac{1\cdot2+9\cdot2}{196}=4-\frac{20}{196}
=\boxed{\tfrac{191}{49}}=3.897959183673\ldots$$
$$I_2^{14}(2)=\tfrac{843}{40}=21.075\ \text{(exact)},\qquad
\tau^\sharp(\mathbf e)=\frac{2}{196}\Bigl(2\cdot6+\tfrac{843}{40}\Bigr)
=\frac{1323}{3920}=\boxed{\tfrac{27}{80}}=0.3375,$$
the minimum being attained on the whole segment $\xi\in[2,\tfrac{13}6]$ (values at
$\xi=2,\ \tfrac{25}{12},\ \tfrac{13}{6}$ all equal $0.3375$).
$$\boxed{\ \tau(\mathbf b;\mathbf e)=\frac{191}{49}+\frac{27}{80}
=\frac{16603}{3920}=4.235459183673469\ldots\ }$$
identical to CDT's value. **[proved, given the measured array; array verified for $n\le50$–$80$]**

---

## 5. Files

| file | role |
|---|---|
| `kfield.py` | $K=\mathbf Q(\sqrt5)$ arithmetic; $\mathbf Z[\phi]$ integrality and least denominator |
| `t1_pure.py` / `.log` | $B_1..B_7$ exact to order 96; CDT Lemma (`bdenominators`) validity; measured minimal types |
| `t1b_zphi.py` / `.log` | $\mathbf Z[\phi]$-integrality of $D_nc_ns^{-n}$, $n\le80$ |
| `t1c_sharp.py` / `.log` | slack $\gcd$s and necessity of each ingredient |
| `t1d_cbin.py` / `.log` | CDT Remark (`central binomial`) tested layer by layer |
| `t2_companions.gp` / `.log` | PARI: $x(q)$, $F$, $\Phi_D$, $R_3$, $R_4$, peels to $x$-order 210; Apéry limits |
| `A_coeffs.txt`, `BD_coeffs.txt`, `BR3_coeffs.txt`, `BR4_coeffs.txt` | exact $x$-expansions, order 210 |
| `t2_G.py` / `.log` | $P_n$ recurrence, $\operatorname{Sym}^+$, $G(Y(x))=H+H\circ w$, denominator type of $G$ |
| `t3_indep.py` / `.log` | rank tables at the split and inert primes, both hosts, both item-13/14 conventions |
| `t3b_pintegral.py` | $p$-integrality guard for the reductions |
| `t4_array.py` / `.log` | measured types of all fourteen; CDT array vs the $u_2=4$ alternative |
| `t4_exact.py` / `.log` | exact $\tau^\flat$, $\tau^\sharp$, $\tau$ via `cdt_bound.py` |

## 6. What is *not* settled

* Full $K(y)$-linear independence in unbounded degree (§3.2) — degrees $\le5$ only. **[assumed]** beyond.
* Denominator types are measured, not proved, for $n$ in the stated ranges; CDT's own proofs
  (explicit factorials for $B_2,B_3,B_5,B_6$; Lemma `etalecover` for $B_4$, $G$; direct
  integration for $B_7$) transport verbatim since $s$ is a unit, but no independent proof is
  written out here.
* Nothing about the entry condition, $\varphi$, the Bost–Charles numerator, or the margin —
  out of scope for this task.
