# The $\chi_{-4}$ Eisenstein class on Atkin–Lehner quotient hosts: the pocket is empty

*Fable (Opus 5), 2026-08-23.  Scripts and logs: `lattice/catalan_al_hosts/`
(`01_hosts.gp`, `04_geom.gp`, `05_rational.gp`, `06_rows.gp`, `07_sing.gp`,
`08_full.gp`, `09_foldreg.gp`, `10_score.py`, `11_po4.gp`; logs `06run.log`,
`07run.log`, `08run.log`, `09run.log`, `09run2.log`).
Conventions from `CDT_FINDER.md` §§1–4, `CDT_NONCONGRUENCE.md` §§3, 5, 10,
`CATALAN_THREE_PERIOD.md` §§1–2, `ADELIC_HOLONOMY.md` §2.6.
Tags: **[proved]** derived here; **[verified]** exact or high-precision
computation in this task; **[estimate]** transported input, flagged at each use.
**No irrationality claim is made anywhere.**

---

## 0. Verdict

| claim | verdict |
|---|---|
| The genus-0 groups $\Gamma_0(N)+W_S$ with $4\mid N$, $N\le120$, are exactly $23$ (list in §1). | **[verified]** §1 |
| For $f\in M_k(\Gamma_0(N),\chi_{-4})$ with $k$ **odd**, $W_Q$ acts with $W_Q^2=\chi_{-4}(-1)\,\chi_{N/Q}(Q)$, so its eigenvalues are $\pm1$ ($Q\equiv1\bmod4$) or $\pm i$ ($4\mid Q$ or $Q\equiv3\bmod4$). Both are legitimate *multiplier systems*; what actually obstructs is **rationality**: the eigenvectors are defined over $\mathbf Q$ **iff every $Q\in S$ is a perfect square**. | **[proved]+[verified]** §2 |
| Consequently the only hosts in the pocket that can carry a row with $A,B\in\mathbf Q[[t]]$ are $\Gamma_0(4)$, $\Gamma_0(4)+4$, $\Gamma_0(8)$, $\Gamma_0(12)$, $\Gamma_0(12)+4$, $\Gamma_0(16)$, $\Gamma_0(16)+16$, $\Gamma_0(20)+4$, $\Gamma_0(36)+4$, $\Gamma_0(36)+4,9$ — **ten** of the twenty-three. | **[verified]** §2 |
| Of these, $\Gamma_0(4)$, $\Gamma_0(4)+4$ and $\Gamma_0(12)+4$ carry **no** row at all: there is no $F\in M_1(N,\chi_{-4})$ with $F\,\mathcal D t$ holomorphic, for any placement of the hauptmodul's pole at a cusp. | **[verified]** §3 |
| Since $4\mid N$ forces $\nu_2(\Gamma_0(N))=\nu_3(\Gamma_0(N))=0$, **every** elliptic point of an AL quotient host is an order-$2$ point coming from a fixed point of an Atkin–Lehner involution. Of the six surviving *proper* quotients, $\Gamma_0(4)+4$ and $\Gamma_0(12)+4$ carry no row and $\Gamma_0(36)+4$ is degenerate; on the remaining three ($\Gamma_0(16)+16$, $\Gamma_0(20)+4$, $\Gamma_0(36)+4,9$) the **fold** (nearest singular point to $t=0$) is such an order-$2$ orbifold point, not a cusp. There is then no Eichler period there, no fold-regular class, and the Apéry limit $b_n/a_n$ is not in $\mathbf Q+\mathbf Q\zeta(2)+\mathbf QG$. **The AL pocket is empty.** | **[verified]** §5 |
| The genuinely new rows are on $\Gamma_0(12)$ itself (trivial $W$): five hauptmodul placements, all carrying $G$ and $\zeta(2)$, with $(\lambda_1,\lambda_2)\in\{(6,4),(4,2),(3,1)\}$. The $(4,2)$ placements have $\operatorname{score}=\log|t_2|-2=-2.6931$, i.e. **$+0.693$ nats better than Zagier E** and $+0.35$ better than the level-16 host. | **[verified]** §4 |
| That gain does not survive: the level-12 hosts have **five** remaining singular points instead of three, no normaliser involution that *removes* the fold, and a hard uniformisation ceiling of at most $\log 4=1.3863$ — Zagier E's exact value — so their entry is no better. | **[verified]** §6 |
| $\lambda_2<1$ occurs **nowhere** in the pocket. The single $\lambda_2=1$ (Apéry-perfect) geometry, $\Gamma_0(12)$ placement 4 with $\Sigma=\{0,\pm\tfrac13,\pm1,\infty\}$, has its fold split into the **conjugate pair** $\pm\tfrac13$ of equal modulus; LLL over $\langle A,B_{\Phi_1},\dots,B_{\Phi_4}\rangle$ finds no combination decaying below $\lambda_1=3$. | **[verified]** §4.3, §7 |
| **Correction to a transported number.** `CATALAN_THREE_PERIOD.md` §5 and `CDT_UNPACKED.md` §4 quote the level-16 host's symmetrised ceiling as $\log(256\cdot\tfrac14)=4.1589$. Computed from the actual $y$-line orbifold ($y=4x^2/(4x+1)$; **both** fixed points of the descending involution are cusps, so the quotient carries **no cone point**) the ceiling is $16$, $\log=2.7726$, and the level-16 entry at the ceiling is $-1.463$, not $-0.077$. Zagier E's $-0.0766$ is unaffected: there the second fixed point $x=\tfrac12$ is a *free* point and does become an order-$2$ cone point. | **[verified]** §6.2 |
| Calibration: Zagier E reproduces $\lambda=(8,4)$, ceiling $\log 64$, entry $-0.0766$, $(S)$ deficit $-0.887$/$-0.632$ with the adelic $\gamma_2$; the level-16 host reproduces $\lambda=(4,2\sqrt2)$, the $4$-dimensional fold-regular space, and all four periods $-\tfrac12P_{\rm in}(2)G$, $\tfrac32(P_{\rm out}(1)-P_{\rm out}(2))\zeta(2)$ of `CATALAN_THREE_PERIOD.md` §2.4 — from an independent construction. | **[verified]** §8 |

One sentence: *widening the Catalan search from three-term rows to every genus-zero
Atkin–Lehner quotient $\Gamma_0(4M)+W$ produces exactly ten candidate hosts, of which
the four proper quotients all put their fold on an order-$2$ orbifold point (so the
hypothesis has nothing to remove) and the rest are level $8$, $12$, $16$; the best new
geometry, on $\Gamma_0(12)$, improves the elementary score by $0.693$ nats but pays two
extra punctures and a lost normaliser descent for it, so **Zagier E's $-0.077$ entry
deficit stands as the best in the whole pocket**.*

---

## 1. The pocket, enumerated

`01_hosts.gp`, `04_geom.gp`.  $\operatorname{genus}(X_0(N)/W_S)=\dim S_2(\Gamma_0(N))^{W_S}$
(joint $+1$-eigenspace, computed with `mfatkininit`); the routine reproduces the classical
Fricke list $\{N: g(X_0(N)+N)=0\}=\{1,\dots,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,29,31,32,35,36,39,41,47,49,50,59,71\}$
exactly, which is the validation. **[verified]**

Since $4\mid N$, $\nu_2(\Gamma_0(N))=0$ (as $4\mid N$) and $\nu_3(\Gamma_0(N))=0$
(as $2\mid N$ and $\bigl(\tfrac{-3}{2}\bigr)=-1$).  So $\Gamma_0(N)$ has **no elliptic
points**, and every elliptic point of a quotient is an order-$2$ point arising from an
interior fixed point of an involution in $W_S$.  Their number follows from the area
formula $e_2=4+\mu'/3-2\nu_\infty'$ with $\mu'=[\mathrm{PSL}_2(\mathbf Z):\Gamma_0(N)]/|W_S|$,
$\nu_\infty'$ = number of $W_S$-orbits on the cusps of $\Gamma_0(N)$ (computed from the
explicit action $\tfrac ac\mapsto W_Q\tfrac ac$ and Diamond–Shurman 3.8.3).

| $N$ | $W_S$ | $\mu'$ | $\#\text{cusps}(\Gamma_0(N))$ | $\#\text{cusps}(\Gamma)$ | $e_2$ | $\#\Sigma$ |
|---|---|---|---|---|---|---|
| 4 | $\{\}$ | 6 | 3 | 3 | 0 | **3** |
| 4 | $+4$ | 3 | 3 | 2 | 1 | **3** |
| 8 | $\{\}$ | 12 | 4 | 4 | 0 | **4** |
| 8 | $+8$ | 6 | 4 | 2 | 2 | **4** |
| 12 | $\{\}$ | 24 | 6 | 6 | 0 | **6** |
| 12 | $+3$ | 12 | 6 | 3 | 2 | 5 |
| 12 | $+4$ | 12 | 6 | 4 | 0 | **4** |
| 12 | $+3,4$ | 6 | 6 | 2 | 2 | **4** |
| 16 | $\{\}$ | 24 | 6 | 6 | 0 | 6 |
| 16 | $+16$ | 12 | 6 | 3 | 2 | 5 |
| 20 | $+4$ | 18 | 6 | 4 | 2 | 6 |
| 20 | $+4,5$ | 9 | 6 | 2 | 3 | 5 |
| 24 | $+8$ | 24 | 8 | 4 | 4 | 8 |
| 24 | $+3,8$ | 12 | 8 | 2 | 4 | 6 |
| 28 | $+7$ | 24 | 6 | 3 | 6 | 9 |
| 28 | $+4,7$ | 12 | 6 | 2 | 4 | 6 |
| 32 | $+32$ | 24 | 8 | 4 | 4 | 8 |
| 36 | $+4$ | 36 | 12 | 8 | 0 | 8 |
| 36 | $+4,9$ | 18 | 12 | 4 | 2 | 6 |
| 44 | $+4,11$ | 18 | 6 | 2 | 6 | 8 |
| 56 | $+7,8$ | 24 | 8 | 2 | 8 | 10 |
| 60 | $+3,4,5$ | 18 | 12 | 2 | 6 | 8 |
| 92 | $+4,23$ | 36 | 6 | 2 | 12 | 14 |

$\#\Sigma=\#\text{cusps}(\Gamma)+e_2$ is the number of singular points of the row's ODE
(the branch locus of $\mathbf H^*\to X(\Gamma)$; the row's $F$ is chosen with divisor on
that locus, see §3).  There is **no genus-0 host with $4\mid N$ beyond $N=92$**, and none
between $60$ and $92$ or between $92$ and $120$. **[verified]**

---

## 2. The rationality obstruction: $Q$ must be a perfect square

> **Proposition 1 [proved].**  Let $\chi=\chi_{-4}$, $Q\,\|\,N$, $4\mid N$, and let
> $W_Q$ act on $M_k(\Gamma_0(N),\chi)$ by $f\mapsto f|_kW_Q$.  Writing
> $\chi=\chi_Q\chi_{N/Q}$, Atkin–Li gives $f|W_Q|W_Q=\chi_Q(-1)\chi_{N/Q}(Q)f$, so
> $$W_Q^2=\begin{cases}-1,& 4\mid Q\ \ (\chi_Q=\chi_{-4}),\\ \chi_{-4}(Q),& Q\ \text{odd},\end{cases}$$
> and the eigenvalues are $\pm1$ if $Q\equiv1\ (4)$ and $\pm i$ otherwise.

An eigenvalue $\pm i$ is **not** an obstruction: it is a multiplier system for the
extended group $\Gamma_0(N)+W_Q$, and a weight-$1$ form with a multiplier system still
gives a rank-$2$ Fuchsian ODE in the quotient hauptmodul.  What obstructs is that
$|_kW_Q$ for **odd** $k$ carries a factor $Q^{k/2}$: on the $V_d$-basis of the
$\chi_{-4}$ Eisenstein spaces the matrix of $W_Q$ is $\zeta\cdot\sqrt Q\cdot(\text{rational})$,
and its eigenvectors are rational **iff $\sqrt Q\in\mathbf Q$**.

**[verified]** (`05_rational.gp`, strict test: the eigenprojector's entries must be
recognised by `bestappr(\cdot,10^6)` to $10^{-25}$ at $60$ digits).  Explicitly, at
$N=8$, $\theta^2|W_8=-i\sqrt2\,V_2\theta^2$ and the eigenvectors are
$\theta^2\mp\sqrt2\,V_2\theta^2$ — the host $\Gamma_0(8)+8$ carries the class only over
$\mathbf Q(\sqrt2)$, exactly like the "$\sqrt{\phantom{x}}$-rows" of
`CDT_NONCONGRUENCE.md` §4.  At $N=16$, $W_{16}$ is $-i$ times a rational matrix with
eigenvalues $\pm1$ and rational eigenvectors.

**Surviving hosts** ($W_S$ generated by $W_Q$'s with $Q\in\{4,9,16\}$):
$$\Gamma_0(4),\ \Gamma_0(4)+4,\ \Gamma_0(8),\ \Gamma_0(12),\ \Gamma_0(12)+4,\
\Gamma_0(16),\ \Gamma_0(16)+16,\ \Gamma_0(20)+4,\ \Gamma_0(36)+4,\ \Gamma_0(36)+4,9 .$$
All of $8+8$, $12+3$, $12+3,4$, $20+4,5$, $24+8$, $24+3,8$, $28+7$, $28+4,7$, $32+32$,
$44+4,11$, $56+7,8$, $60+3,4,5$, $92+4,23$ are killed. **[verified]**

*(Aside, also verified: for $Q\equiv1\bmod4$ and $4\mid N/Q$, $W_Q$ acts **freely** on
$X_0(N)$ — the discriminant $-4Q\equiv12\bmod16$ is not a square mod $16$ — so
$\operatorname{tr}(W_Q\mid S_2)=1$ and $g(X_0(N)/W_Q)=\tfrac{g+1}2\ge1$: no genus-0
quotient by a $Q\equiv1\bmod4$ involution exists at all when $4\mid N$.)*

---

## 3. Building the rows

For a host $\Gamma=\Gamma_0(N)+W_S$ of genus $0$ the data is $(t,F)$ with $t$ the
hauptmodul normalised $t=q+O(q^2)$ and $F\in M_1(\Gamma_0(N),\chi_{-4})$ a joint
$W_S$-eigenvector with $F(0)=1$; the row is $A(t)=F$, $B_\Phi=F\,\mathcal D^{-2}\Phi_+$
with $\Phi$ in the matching eigenspace of $M_3(\Gamma_0(N),\chi_{-4})$.  The two
constraints are
$$\Phi_0:=F\cdot\mathcal D t\in M_3(\Gamma_0(N),\chi_{-4})\ \ (\text{holomorphy}),
\qquad A\in\mathbf Z[[t]] .$$

**Hauptmoduls.**  $t$ is a modular unit with $\operatorname{div}(t)=(\infty_\Gamma)-(P)$
for a cusp $P$ of $\Gamma$; pulling back to $X_0(N)$ and solving Ligozat's linear system
$$\operatorname{ord}_c\Bigl(\prod_{d\mid N}\eta(d\tau)^{r_d}\Bigr)
=\frac{N}{24\gcd(c,N/c)\,c}\sum_{d\mid N}\frac{\gcd(c,d)^2}{d}\,r_d$$
for $r\in\mathbf Q^{\sigma_0(N)}$ gives $t$ from
$\mathcal D\log t=\tfrac1{24}\sum_d r_d\,d\,E_2(d\tau)$, $t=q\exp\bigl(\mathcal D^{-1}(\mathcal D\log t-1)\bigr)$
(`lib2.gp`).  Rational $r_d$ are allowed ($t$ is then a root of an eta quotient; its
$q$-expansion is still in $\mathbf Q[[q]]$).  Every cusp orbit $P\ne\infty_\Gamma$ is tried
— this is the cusp-move family of `CUSP_MOVE_PROGRAM.md`, here enumerated exactly.

**$F$** is then the solution of a linear system: the columns
$\{$coefficients of $F_i\,\mathcal Dt\}$ against the $M_3$ basis, kernel taken to
$q^{240}$ (`06_rows.gp`).

**Result.**  $\Gamma_0(4)$, $\Gamma_0(4)+4$ and $\Gamma_0(12)+4$ admit **no** $F$ for any
pole placement.  (For $\Gamma_0(12)+4$ this is the interesting one: the host has the same
signature as $\Gamma_0(6)$ — $\mu'=12$, $4$ cusps, no elliptic points — i.e. *CDT's own
geometry*, but $M_1(\Gamma_0(12),\chi_{-4})^{W_4=\pm i}$ is one-dimensional and its unique
member does not clear the pole of $\mathcal Dt$ at any of the three candidate cusps.)
**[verified]** `06run.log`.

---

## 4. The rows on $\Gamma_0(8)$, $\Gamma_0(12)$, $\Gamma_0(16)$

`07_sing.gp` fits the rank-$2$ Fuchsian ODE $q_2(t)A''+q_1(t)A'+q_0(t)A=0$ from $200$
coefficients and reads $\Sigma$ off $q_2$; `08_full.gp`, `09_foldreg.gp` measure
$|a_n|^{1/n}$, the fold-regular subspace, $\xi=\lim b_n/a_n$, its identification by
`lindep` against $(1,\zeta(2),G)$, the decay $|b_n-\xi a_n|^{1/n}$, and $k$.

### 4.1 Calibration rows

| host, placement | $F$ | $\Sigma$ | $\lambda_1$ | $\lambda_2$ | $\dim$FR | periods |
|---|---|---|---|---|---|---|
| $\Gamma_0(8)$, pole at $c=2$ (**Zagier E**) | $\theta^2$ | $\{0,\tfrac18,\tfrac14,\infty\}$ | $8$ | $4$ | $2/4$ | $-\tfrac12G$ (inner), $\tfrac38\zeta(2)$ (outer) |
| $\Gamma_0(16)$, pole at $c=2$ | $\theta^2$ | $\{0,\tfrac14,\tfrac12,\tfrac{1\pm i}4,\infty\}$ | $4$ | $2\sqrt2$ | $4/6$ | $-\tfrac1{16}G,\ -\tfrac3{128}G,\ -\tfrac38\zeta(2),\ -\tfrac{15}{32}\zeta(2)$ |
| $\Gamma_0(16)$, pole at $c=1$ | $\theta^2$ | $\{0,-\tfrac14,-\tfrac12,\tfrac{-1\pm i}4,\infty\}$ | $4$ | $2\sqrt2$ | $4/6$ | $-\tfrac12G,\ -\tfrac1{64}G,\ \tfrac38\zeta(2),\ \tfrac9{32}\zeta(2)$ |

The last line is `CATALAN_THREE_PERIOD.md` §1.3–1.4 exactly:
$\Phi=E\Rightarrow\xi=-\tfrac12P_{\rm in}(2)G=-\tfrac12G$;
$\Phi=-\tfrac18V_2E+V_4E\Rightarrow P_{\rm in}(2)=\tfrac1{32}$, $\xi=-\tfrac1{64}G$;
$\Phi=-T+V_2T\Rightarrow\tfrac32(P(1)-P(2))=\tfrac38$;
$\Phi=-T+V_4T\Rightarrow\tfrac9{32}$.  All to $\ge38$ digits. **[verified]**

$\Gamma_0(16)$ also carries a **second** weight-1 form on the same $(\Gamma,t)$,
$F=1-4q^2+4q^4+\cdots=V_2\theta^2$-twisted, with $A=1,0,-4,16,-44,96,\dots$, the same
$\Sigma$, the same $\lambda$'s and the same fold-regular periods.

### 4.2 The new rows: $\Gamma_0(12)$

$\dim M_1(12,\chi_{-4})=2$, $\dim M_3^{\rm Eis}(12,\chi_{-4})=4$
($E,V_3E$ inner; $T,V_3T$ outer).  Five of the six cusp placements carry a row.

| placement (pole at $c=$) | $F$ | $\Sigma$ | $\lambda_1$ | $\lambda_2$ | $\dim$FR | inner period | outer period |
|---|---|---|---|---|---|---|---|
| $1$ | $1-2q-2q^2+6q^3-\cdots$ | $\{0,-\tfrac16,-\tfrac14,-\tfrac13,-\tfrac12,\infty\}$ | $6$ | $4$ | $2/4$ | $\tfrac19G$ | $\tfrac29\zeta(2)$ |
| $2$ | $\theta^2$ | $\{0,\tfrac16,\tfrac14,\tfrac13,\tfrac12,\infty\}$ | $6$ | $4$ | $2/4$ | $-\tfrac1{27}G$ | $-\tfrac49\zeta(2)$ |
| $3$ | $1+2q+2q^2+2q^3+\cdots$ | $\{0,\tfrac14,\pm\tfrac12,1,\infty\}$ | $\mathbf 4$ | $\mathbf 2$ | $2/4$ | $-\tfrac1{27}G$ | $-\tfrac49\zeta(2)$ |
| $6$ | $1+4q^3+4q^6+\cdots$ | $\{0,-\tfrac14,\pm\tfrac12,-1,\infty\}$ | $\mathbf 4$ | $\mathbf 2$ | $2/4$ | $\tfrac19G$ | $\tfrac29\zeta(2)$ |
| $4$ | $1+q+q^2+3q^3+\cdots$ | $\{0,\pm\tfrac13,\pm1,\infty\}$ | $3$ | $(1)$ | $0$ | — | — |

All $A\in\mathbf Z[[t]]$ to $n=234$; measured $|a_n|^{1/n}$ at $n=234$: $5.83\to6$,
$5.87\to6$, $3.907\to4$, $3.876\to4$, $2.925\to3$.  Decay of the fold-regular linear
forms $|b_n-\xi a_n|^{1/n}$ at $n=195$: $3.89\to4$, $3.86\to4$, $1.92\to2$, $1.94\to2$.
Periods identified to $\ge28$ digits by `lindep` with bound $10^{30}$.  $k=2$
(sharp; $k=1$ fails) on every fold-regular class with an integral primitive
representative that was tested. **[verified]**

**The $(4,2)$ placements are the best elementary geometry for Catalan found anywhere**:
$$\operatorname{score}=\log|t_2|-2=\log\tfrac12-2=-2.6931,$$
against $-3.3863$ for Zagier E and $-3.0397$ for the level-16 host.

### 4.3 The Apéry-perfect placement, and why it fails

Placement $4$ has $\Sigma=\{0,\pm\tfrac13,\pm1,\infty\}$ — $\lambda_1\lambda_2=3$,
$|t_2|=1$, i.e. exactly the Fricke-palindromic / Apéry-perfect condition
$\lambda_2^{\rm norm}=1$ that `CDT_NONCONGRUENCE.md` §3.3 identifies as the **only**
geometry where architecture $(S)$ works.  It fails for a different reason: the fold is a
**pair** $\pm\tfrac13$ of equal modulus, so $b_n/a_n$ oscillates and there is no Apéry
limit.  A class regular at *both* folds would decay at rate $1$ and would be a
sensational object; `11_po4.gp` runs LLL on
$\langle A,B_{\Phi_1},\dots,B_{\Phi_4}\rangle$ over a $12$-term window at $n\simeq190$
and every reduced combination still grows at $|c_n|^{1/n}\to3=\lambda_1$.  **No such
class exists in the $\chi_{-4}$ Eisenstein space.** **[verified]**

---

## 5. The proper Atkin–Lehner quotients: the fold is an orbifold point

| host | placement | $\Sigma$ | $\lambda_1$ | $\lambda_2$ | fold is | $\dim$FR | $\lim b_n/a_n$ |
|---|---|---|---|---|---|---|---|
| $\Gamma_0(16)+16$ | $c=2$ | $\{0,\tfrac{3-2\sqrt2}2,\tfrac12,\tfrac{3+2\sqrt2}2,\infty\}$ | $6+4\sqrt2=11.657$ | $2$ | **elliptic** | $1/3$ | not in $\mathbf Q+\mathbf Q\zeta(2)+\mathbf QG$ |
| $\Gamma_0(16)+16$ | $c=4$ | $\{0,\tfrac{\sqrt2-1}4,-\tfrac12,-\tfrac{\sqrt2+1}4,\infty\}$ | $4(\sqrt2+1)=9.657$ | $2$ | **elliptic** | $1/3$ | not identified |
| $\Gamma_0(20)+4$ | $c=1$ | $\{0,-\tfrac15,\tfrac{-1\pm2i}5,-1,\infty\}$ | $5$ | $\sqrt5$ | elliptic | $0/2$ | — |
| $\Gamma_0(20)+4$ | $c=2$ | $\{0,\tfrac15,\tfrac14,\tfrac{2\pm i}{10},\infty\}$ | $5$ | $2\sqrt5$ | elliptic | $0/2$ | — |
| $\Gamma_0(20)+4$ | $c=10$ | $\{0,-\tfrac14,\pm\tfrac i2,1,\infty\}$ | $4$ | $2$ | elliptic | $0/2$ | — |
| $\Gamma_0(36)+4,9$ | $c=2$ | $\{0,\tfrac{2-\sqrt3}2,\tfrac14,1,\tfrac{2+\sqrt3}2,\infty\}$ | $2(2+\sqrt3)=7.464$ | $4$ | elliptic | $0/2$ | — |
| $\Gamma_0(36)+4,9$ | $c=3$ | $\{0,\tfrac1{3+2\sqrt3},\tfrac13,-1,-\tfrac{3+2\sqrt3}3,\infty\}$ | $3+2\sqrt3=6.464$ | $3$ | elliptic | $0/2$ | — |
| $\Gamma_0(36)+4,9$ | $c=6$ | $\{0,-\tfrac13,\pm\tfrac1{2\sqrt3},-\tfrac14,\infty\}$ | $4$ | $2\sqrt3$ | elliptic | $0/2$ | — |
| $\Gamma_0(36)+4$ | all | (degenerate; $A$ algebraic, see below) | — | — | — | — | — |

The mechanism, verified in detail on $\Gamma_0(16)+16$ (`07run.log` + the explicit
descent): $W_{16}$ acts on the level-16 $x$-line by $\sigma(x)=\tfrac1{8x}$ (it swaps the
cusps $0\leftrightarrow\infty$ and $x=-\tfrac14\leftrightarrow x=-\tfrac12$, and swaps the
Galois pair $\tfrac{-1\pm i}4$).  Its fixed points $x=\pm\tfrac1{2\sqrt2}$ are **not**
cusps: they are the two order-$2$ elliptic points of the quotient.  In the quotient
coordinate $t=x/(1+bx+8x^2)$ (a hauptmodul for every $b$; $b=6$ is the placement with the
pole at $\{-\tfrac14,-\tfrac12\}$) they sit at $t=\tfrac1{4\sqrt2+b}$ and
$t=-\tfrac1{4\sqrt2-b}$, and for $b=6$ the first of these, $t=\tfrac1{6+4\sqrt2}=0.0858$,
is **nearer to $0$ than any cusp** (the surviving cusp is at $t=\tfrac12$).

Consequences, all **[verified]**:
* $|a_n|^{1/n}\to6+4\sqrt2$ ($10.18,\,10.80,\,11.17$ at $n=59,118,236$): the elliptic
  point is a genuine $\sqrt{\ }$-branch point of $A$, not an apparent singularity.
* $b_n/a_n$ converges to the ratio of the two leading branch coefficients, which is not a
  period: $30$-digit `lindep` against $(1,\zeta(2),G)$ with bound $10^{7}$ returns nothing
  on any of the three Eisenstein classes, and on the one-dimensional "slow-kernel"
  direction the decay $|b_n-\xi a_n|^{1/n}$ **rises** to $\lambda_1$ rather than falling
  to $\lambda_2$.
* Hence there is no Eichler period polynomial at the fold, no fold-regularity condition,
  no conditional function, and no Apéry mechanism — irrespective of the CDT bound.

$\Gamma_0(36)+4$ is degenerate in a different way: its two $F$'s give $A=1-t$ and
$A=(1+t)^{-2}$-type algebraic series (the fitted ODE has no finite-degree $q_2$), i.e.
$t$ there is not a hauptmodul in the required sense.

> **Structural reading [proved].** $4\mid N\Rightarrow\nu_2=\nu_3=0$ on $\Gamma_0(N)$, so
> *all* the orbifold points of an AL quotient come from AL fixed points, i.e. from CM
> points of discriminant $-4Q$ or $-Q$.  Those are the points of $X_0(N)$ **closest to the
> cusp $\infty$ in the $t$-coordinate** on these hosts, because the quotient map contracts
> the two long cusp-neighbourhoods that used to be $\infty$ and $W_Q\infty$.  The fold of
> an AL quotient host is therefore generically an orbifold point, and the Eisenstein
> extension class has nothing to be regular at.  *Atkin–Lehner descent buys cusps at the
> price of orbifold points, and the Catalan mechanism lives only at cusps.*

---

## 6. Scoring

`10_score.py`.  $\tau(\mathbf b;\mathbf e)=16603/3920=4.235459$ for CDT's $m=14$, $k=2$
inventory (recomputed from `cdt_finder/cdt_bound.py`); contour loss $\log0.62922=-0.4633$;
$\mathrm{BC}=11.845+\log s$.

### 6.1 Ceilings, computed rather than transported

The ceiling is the conformal radius at $t=0$ of the universal (orbifold) cover of
$\Omega=\mathbf P^1\setminus(\Sigma\setminus\{\text{fold}\})$, with $t=0$ counted as a
puncture (legitimate because the row's power series have trivial local monodromy there —
that is exactly what lets $\varphi$ be the covering map extended across $z=0$).  Exact
values used:
$$r\bigl(\mathbf P^1\setminus\{0,s,\infty\}\bigr)=16|s|,\qquad
r\bigl(\mathbf P^1\setminus\{0,a,b\}\bigr)=16\Bigl|\frac{ab}{b-a}\Bigr|,\qquad
r\bigl(\{0,\infty\}\text{ punctures}+\text{cone-}2\text{ at }c\bigr)=64|c| ,$$
the last being the $\Gamma_0(2)$ hauptmodul $-256\Delta(2\tau)/\Delta(\tau)$ and the source
of CDT's $256s=64\cdot(4s)$.  For $\#\Sigma\setminus\{\text{fold}\}>3$ we use monotonicity,
$r(\Omega)\le\min$ over $3$-point subsets — a rigorous upper bound.

### 6.2 The level-16 correction

CDT's normaliser descent is $w(x)=\tfrac{sx}{x-s}$, $y=x+w(x)=\tfrac{x^2}{x-s}$.  On
Zagier E ($s=\tfrac14$) its fixed points are $x=0$ (a cusp) and $x=2s=\tfrac12$, which is
**not** in $\Sigma$; the latter becomes an order-$2$ cone point at $y=4s=1$ and the
quotient orbifold is $\{0,\infty\text{ punctures}\}+\text{cone-}2$, ceiling
$64\cdot1=64=256s$. **[verified]**

On the level-16 host the descending involution is $\sigma(x)=-\tfrac x{4x+1}$,
$y=\tfrac{4x^2}{4x+1}$ (the map used in `CATALAN_THREE_PERIOD.md` §5).  Its fixed points
are $x=0$ and $x=-\tfrac12$, and **both are cusps of $\Gamma_0(16)$**.  The quotient
therefore carries **no cone point**: $\Sigma_y=\{0,-\tfrac12,-1,\infty\}$ with the fold
at $y=\infty$, so
$$\Omega_y=\mathbf P^1\setminus\{0,-\tfrac12,-1\},\qquad
|\varphi'(0)|\le 16\Bigl|\frac{(-\tfrac12)(-1)}{-1+\tfrac12}\Bigr|=16,\qquad \log=2.7726,$$
and the level-16 entry at the ceiling is $2.7726-4.2355=\mathbf{-1.463}$, not the
$-0.077$ obtained by transporting $256|s|$ with $|s|=\tfrac14$.  **This makes the
level-16 host substantially worse than recorded, and leaves Zagier E as the unique host
in the pocket whose symmetrised ceiling reaches $64$.** **[verified]**; it corrects
`CATALAN_THREE_PERIOD.md` §5 and `CDT_UNPACKED.md` §4.2, and should be double-checked
against a direct numerical uniformisation before being used in a paper.

### 6.3 The table

$\lambda_2$ = second distinct singular modulus's reciprocal.  "$\#$rem" $=\#(\Sigma\setminus\{\text{fold}\})$
counting $0$ and $\infty$.  ceil$_D$ = hard unsymmetrised ceiling (rigorous upper bound);
entry$_D$ uses $\tau=16/9$ of `CDT_NONCONGRUENCE.md` §3.1 ($m=3$).  ceil$_S$/entry$_S$ =
symmetrised architecture, computed from the actual quotient orbifold; margin$_S$ uses
CDT's transported contour loss and $\mathrm{BC}$ **[estimate]**.

| host / class | period | $\Sigma$ type | $\lambda_1$ | $\lambda_2$ | $k$ | $\#$rem | score | ceil$_D$ | entry$_D$ | ceil$_S$ | entry$_S$ | margin$_S$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **$\Gamma_0(8)$ (Zagier E)** inner | $-\tfrac12G$ | 4 cusps | $8$ | $4$ | 2 | 3 | $-3.3863$ | $4$ | $-0.391$ | $\mathbf{64}$ | $\mathbf{-0.0766}$ | $-18.02$ |
| $\Gamma_0(8)$ outer | $\tfrac38\zeta(2)$ | 4 cusps | $8$ | $4$ | 2 | 3 | $-3.3863$ | $4$ | $-0.391$ | $64$ | $-0.0766$ | $-18.02$ |
| **$\Gamma_0(16)$** inner | $-\tfrac12G,-\tfrac1{64}G$ | 6 cusps | $4$ | $2\sqrt2$ | 2 | 5 | $-3.0397$ | $\le4$ | $\le-0.391$ | $16$ | $-1.463$ | $-37.8$ |
| $\Gamma_0(16)$ outer | $\tfrac38\zeta(2),\tfrac9{32}\zeta(2)$ | 6 cusps | $4$ | $2\sqrt2$ | 2 | 5 | $-3.0397$ | $\le4$ | $\le-0.391$ | $16$ | $-1.463$ | $-37.8$ |
| $\Gamma_0(12)$ pl. 1 | $\tfrac19G$, $\tfrac29\zeta(2)$ | 6 cusps | $6$ | $4$ | 2 | 5 | $-3.3863$ | $\le4$ | $\le-0.391$ | — | — | — |
| $\Gamma_0(12)$ pl. 2 | $-\tfrac1{27}G$, $-\tfrac49\zeta(2)$ | 6 cusps | $6$ | $4$ | 2 | 5 | $-3.3863$ | $\le4$ | $\le-0.391$ | — | — | — |
| **$\Gamma_0(12)$ pl. 3** | $-\tfrac1{27}G$, $-\tfrac49\zeta(2)$ | 6 cusps | $4$ | $\mathbf 2$ | 2 | 5 | $\mathbf{-2.6931}$ | $\le4$ | $\le-0.391$ | $\le7.11$ | $\le-2.274$ | $\le-49$ |
| **$\Gamma_0(12)$ pl. 6** | $\tfrac19G$, $\tfrac29\zeta(2)$ | 6 cusps | $4$ | $\mathbf 2$ | 2 | 5 | $\mathbf{-2.6931}$ | $\le4$ | $\le-0.391$ | $\le7.11$ | $\le-2.274$ | $\le-49$ |
| $\Gamma_0(12)$ pl. 4 | none | 6 cusps | $3$ | $(1)$ | — | 5 | $(-2.0)$ | — | — | — | — | — |
| $\Gamma_0(16)+16$ (2 rows) | none | 3 cusps + 2 ell | $11.66,\,9.66$ | $2$ | — | 4 | — | — | — | — | — | — |
| $\Gamma_0(20)+4$ (3 rows) | none | 4 cusps + 2 ell | $5,5,4$ | $\sqrt5,2\sqrt5,2$ | — | 5 | — | — | — | — | — | — |
| $\Gamma_0(36)+4,9$ (3 rows) | none | 4 cusps + 2 ell | $2(2{+}\sqrt3),\,3{+}2\sqrt3,\,4$ | $4,3,2\sqrt3$ | — | 5 | — | — | — | — | — | — |
| $\Gamma_0(36)+4$ | degenerate ($A$ algebraic) | 8 cusps | — | — | — | — | — | — | — | — | — | — |

The $\Gamma_0(12)$ ceil$_S$ entry deserves a word.  There **is** an involution of the
$t$-line preserving $\Sigma$, namely $\sigma(t)=\tfrac1{4t}$ (it swaps $0\leftrightarrow\infty$
and the fold $\tfrac14\leftrightarrow1$, and fixes $\pm\tfrac12$) — but it swaps the fold
with another **cusp**, so the symmetrised function $\mathrm{Sym}^+H$ is still singular at
the common image $v=\tfrac45$, and its two fixed points $\pm\tfrac12$ are cusps, so again
no cone point appears.  $\Omega_v=\mathbf P^1\setminus\{0,\tfrac45,1,-1\}$, ceiling
$\le7.11$.  **The descent buys nothing on $\Gamma_0(12)$.** **[verified]**

### 6.4 Reading

* The score ordering is real: $\Gamma_0(12)$ pl. 3/6 $(-2.69)$ $>$ $\Gamma_0(16)$ $(-3.04)$
  $>$ Zagier E $=\Gamma_0(12)$ pl. 1/2 $(-3.39)$.  In the *elementary* (Beukers) and in the
  $(K)/(D)$ architectures of `CDT_NONCONGRUENCE.md` §3 the level-12 $(4,2)$ rows are the
  best Catalan rows known: they beat Zagier E by $0.693$ nats of score.
* But in every architecture that actually reaches the entry test, the hard ceiling is
  controlled by the **whole** remaining singular set, not by $\lambda_2$ alone.  On the
  level-12 rows the pair $\pm\tfrac12$ already forces
  $r\le16\bigl|\tfrac{(\frac12)(-\frac12)}{-1}\bigr|=4$ — numerically identical to
  Zagier E's exact ceiling of $4$ — and the two extra punctures $1$, $\infty$ push it
  strictly below.  The score gain is exactly cancelled.
* $\lambda_2<1$ **nowhere**.  The minimum over the pocket is $\lambda_2=1$ at
  $\Gamma_0(12)$ pl. 4, which has no Apéry limit (§4.3).

---

## 7. What was looked for and not found

* **A row on $\Gamma_0(4)$.**  $\dim M_1(4,\chi_{-4})=1$ ($F=\theta^2$) and
  $\dim M_3(4,\chi_{-4})=2$, whose $a_0=0$ part is $\langle E\rangle$; so
  $\Phi_0=E$ is forced and $t=\mathcal D^{-1}(E/\theta^2)=q+\tfrac43q^3+\cdots$ is not
  integral and not a modular unit.  There is no $\chi_{-4}$ row on $\Gamma_0(4)$ with the
  normalisation $t=q+O(q^2)$.  (The hypergeometric row $\sum\binom{2n}n^2t^{2n}$ does
  appear — as the placement of $\Gamma_0(8)$ whose hauptmodul $t=q-4q^3+14q^5-\cdots$ is
  odd, i.e. as a level-8 object in the $q^2$-coordinate.) **[verified]**
* **A Catalan row with a decaying linear form beyond $\lambda_2=1$.**  None: every
  $\lambda_2$ in the pocket is $\ge1$, and $\lambda_2=1$ occurs once, on a row whose fold
  is a $\pm$ pair.
* **A mixed conditional class $E\pm rT$.**  On every live host the fold-regularity
  functionals *separate* into an inner functional (on the $E$-coordinates) and an outer
  one (on the $T$-coordinates), exactly as `CATALAN_THREE_PERIOD.md` §2.4 found at level
  16: the fold-regular space is $\{\text{inner, period}\in\mathbf QG\}\oplus\{\text{outer,
  period}\in\mathbf Q\zeta(2)\}$, of dimension $2$ at levels $8$ and $12$ and $4$ at
  level $16$.  The three-period hypothesis $a+b\zeta(2)+cG=0$ therefore yields, on each
  host, **one** conditional generator — the mixed class $\mu+rG+s\zeta(2)=0$ — exactly as
  in `CATALAN_THREE_PERIOD.md` Proposition 2, and no new function.  **[verified]** on
  levels 8, 12, 16.
* **Doubly-small ($\xi=0$) classes on the new hosts.**  On $\Gamma_0(12)$ the fold-regular
  space is $2$-dimensional with the two periods $\tfrac19G$ (or $-\tfrac1{27}G$) and
  $\tfrac29\zeta(2)$ (or $-\tfrac49\zeta(2)$), both nonzero on the whole line: the map
  $\Phi\mapsto\xi(\Phi)$ is injective on each of the two $1$-dimensional pieces, so
  **there is no doubly-small class on $\Gamma_0(12)$** — unlike level 16, where the
  fold-regular space is $4$-dimensional and contains a $2$-dimensional doubly-small plane.
  This is a dimension effect: $\dim M_3^{\rm Eis}(12,\chi_{-4})=4$ against $6$ at level 16.
  **[verified]**

---

## 8. Reproduction

```
gp -q lattice/catalan_al_hosts/01_hosts.gp        # genus-0 AL quotients, 4|N
gp -q lattice/catalan_al_hosts/04_geom.gp         # cusps / elliptic points / #Sigma
gp -q lattice/catalan_al_hosts/05_rational.gp     # the perfect-square criterion
gp -q lattice/catalan_al_hosts/06_rows.gp         # hauptmoduls, F, A  -> 06run.log
gp -q lattice/catalan_al_hosts/07_sing.gp         # ODE fit, singular sets -> 07run.log
gp -q lattice/catalan_al_hosts/08_full.gp         # growth, companions -> 08run.log
gp -q lattice/catalan_al_hosts/09_foldreg.gp      # fold-regular spaces, periods -> 09run*.log
gp -q lattice/catalan_al_hosts/11_po4.gp          # LLL test on the (3,1) placement
python3 lattice/catalan_al_hosts/10_score.py      # ceilings, entry, margins
```

Runtimes: all under 20 min on the dev box; nothing needed `snake`.

## 9. Honest ledger

**Computed exactly / verified in this task.**  The genus-0 census and its validation
against the Fricke list; the cusp/elliptic counts; the $W_Q^2$ values and the
perfect-square rationality criterion (strict test at $60$ digits); the hauptmoduls as
Ligozat units and their $q$-expansions; the existence/non-existence of $F$ on each host
and each cusp placement; $A\in\mathbf Z[[t]]$ to $n=234$; the rank-$2$ ODEs and the exact
singular sets; $|a_n|^{1/n}$ and $|b_n-\xi a_n|^{1/n}$ at $n\le236$; the fold-regular
subspaces; the periods to $\ge28$ digits with `lindep`; $k=2$ sharp; the LLL test on
$\Gamma_0(12)$ pl. 4; the conformal ceilings from the exact three-point and
cone-point formulas; $\tau=16603/3920$.

**Cited, not recomputed.**  CDT's Bost–Charles integral $11.845$ and contour loss
$0.62922$; the architecture thresholds of `CDT_NONCONGRUENCE.md` §3; the adelic
$\gamma_2=0.2546$ of `ADELIC_HOLONOMY.md` §4.

**Estimated / flagged.**  Every margin$_S$ column transports CDT's contour shape and
$\mathrm{BC}$ to a host with a *different* orbifold, which is not justified; the numbers
are indicative only.  The ceilings for $\#$rem $>3$ are rigorous **upper** bounds, not
attained values; a genuine numerical uniformisation of the $4$- and $5$-punctured spheres
was not carried out (it would only make the entries worse).  No pure module was
constructed on any new host; no $\mathbf Q(t)$-independence was checked; no
$p$-adic slopes were measured on the new hosts.

**Not done.**  The non-Atkin–Lehner part of the normaliser (the groups $n|h+e,f,\dots$
with $h^2\mid N$, $h\mid24$) was not scanned; those are the only genus-0 congruence hosts
in the $\chi_{-4}$ world that this note does not cover.
