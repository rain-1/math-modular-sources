# The level-60 ζ(7) "near-miss": what it is, and why it is not a positive-score system

*Claude (Opus 5), 2026-08-22. Scripts: `lattice/zeta7_level60/`. Exact PARI arithmetic
throughout; every number below is measured here unless marked `[archive]`.*

Archive provenance: `chatgpt-research-archive/conversations/apery-systems/zeta7-legti-6a8353d9.md`
(source, filtration, genus-two quotient) and
`.../apery-systems-4/continue-modular-proof-6a86a8fb.md` (cocycle kernels, 19/6, log R,
depth-1 obstruction); catalogued in `ARCHIVE_CATALOGUE.md` rows 3, 86, 87.

---

## 0. Verdict

**(b): the claim was conditional on an unproved lift, and the lift fails.** It is not an
arithmetic error — every number in the archive that I could recompute reproduces exactly —
but it is also not a positive-score Apéry system, because **there is no Apéry system**.
Three independent facts, each established here:

1. **No scalar level-60 system exists.** The only genus-zero Atkin–Lehner quotient of
   $X_0(60)$ is the full one, and **no purified weight-8 source is Atkin–Lehner invariant**
   in any sign pattern. So there is no host, no $t$, no $F=\Phi/Dt$, no $A_n$, no $B_n$, no
   recurrence, no $\lambda_2$ and no $k$ to measure. The archive never built one either.
2. **$\log R = 3.4639$ is not a $\lambda_2$.** It is $\log$ of the *largest* of the eight
   branch-point moduli of the genus-two isogeny curve. Using it as an analytic radius
   requires the linear form to continue past **all seven smaller** singular orbits — stated
   as unproved in the archive itself.
3. **Every monodromy-purified source at level 60 costs the full $k=7$.** Computed here:
   the three Atkin–Lehner kernel lines $K_3,K_4,K_5$ all have support depth 1, hence
   denominator entropy exactly 7 — not 19/6, not 7/3. $7 > 3.4639$ by a factor $e^{3.54}\approx34$.

The archive's own endpoint (fallback entropy $11/3=3.6667 > 3.4639$, score $-0.2027$) is
therefore the optimistic bound, not the honest one.

---

## 1. The source (Task 1 — verified exactly)

`lattice/zeta7_level60/source.gp`. Oldform vector on $d\mid 60$,
$\Phi=\sum_d c_d E_8(d\tau)$, $L(\Phi,s)=P(s)\zeta(s)\zeta(s-7)$, $P(s)=\sum_d c_d d^{-s}$;
purity $\iff P(0)=P(2)=P(4)=P(6)=P(8)=0$; $L(\Phi,7)=(-P(7)/2)\zeta(7)$.

* $\dim(\text{purified})=7$ inside the 12-dimensional oldform space;
  $\dim(\text{purified}\cap\{W_{60}=+1\})=\mathbf3$. **[archive: 3 ✓]**
* Basis (archive's, re-verified purified and $W_{60}=+1$):
  $B_0,B_1,B_2$ with $L(B_i,7)/\zeta(7)=\tfrac{2623}{216000},\tfrac{949}{54000},\tfrac{80347}{216000}$.
  *(The archive writes "$L(B_0/480,7)=\tfrac{2623}{216000}\zeta(7)$"; the $/480$ is a
  transcription slip — the value belongs to $B_0$ itself.)*
* **Support-depth filtration** $F_m=\{c_d=0\ \forall d<m\}$:
  $\dim(F_1,F_2,F_3,F_4)=(3,2,1,0)$. **[archive: 3,2,1 ✓]** $B_0$ is the unique depth-3 line.
* **The $p=5$ completion.**
  $\Phi_{60}^{(5)}=(1-5^4V_5)\Phi_{12}$ with $\Phi_{12}$ the level-12 anti-Fricke parent
  $c=(1,-572,11583,-36608,46332,-20736)$ on $d\mid12$. Explicitly on
  $d=(1,2,3,4,5,6,10,12,15,20,30,60)$:
  $$c=(1,-572,11583,-36608,-625,46332,357500,-20736,-7239375,22880000,-28957500,12960000).$$
  Kills $s=0,2,4,6,8$ ✓; $\ L(\Phi_{60}^{(5)},7)=\tfrac{6479}{54000}\zeta(7)$ ✓
  $\bigl(=\tfrac{209}{1728}\cdot\tfrac{124}{125}\zeta(7)$, the $\tfrac{124}{125}$ being the
  Euler factor $1-5^{4-7}\bigr)$.
* **Atkin–Lehner** (via $E_k(d\tau)|W_Q=(Q/a^2)^{k/2}E_k(dQ/a^2\,\tau)$, $a=Q$-part of $d$):
  $W_5=-1$, $W_{12}=-1$, $W_{60}=+1$; **not** an eigenvector for $W_3,W_4,W_{15},W_{20}$. ✓
* Coordinates in the archive's basis: $\Phi_{60}^{(5)}=(\tfrac{13}{3},-\tfrac{52}{3},1)$,
  i.e. the line $\mathbb Q(13,-52,3)=K_5$. **[archive ✓ exactly]**

---

## 2. The host: a no-go (Task 2)

`lattice/zeta7_level60/host_genus.gp`, `allplus.gp`. Genera computed as
$\dim S_2(\Gamma_0(60))^{H}$ via `mfatkininit`.

| quotient | 1 | $W_3$ | $W_4$ | $W_5$ | $W_{12}$ | $W_{15}$ | $W_{20}$ | $W_{60}$ | $\langle W_5,W_{60}\rangle$ | $\langle W_3,W_5\rangle$ | full |
|---|---|---|---|---|---|---|---|---|---|---|---|
| genus | 7 | 4 | 3 | 4 | 4 | 1 | 2 | 3 | **2** | 1 | **0** |

* $X_0(60)/\langle W_5,W_{60}\rangle$ has genus **2** — an independent confirmation of the
  archive's genus-two curve, which is exactly the quotient on which $\Phi_{60}^{(5)}$
  ($W_5=-1$, $W_{60}=+1$) descends.
* The **unique** genus-zero quotient is $X_0(60)^*=X_0(60)/\langle W_3,W_4,W_5\rangle$.
  A scalar Apéry system needs $t$ a Hauptmodul there and $F=\Phi/Dt$ a function there,
  i.e. $\Phi|W_Q=+\Phi$ for **all** $Q\mid60$.
* $\dim(\text{oldform},\ \text{all }W_Q=+1)=2$, but
  $$\boxed{\ \dim(\text{purified},\ W_3=\varepsilon_3,W_4=\varepsilon_4,W_5=\varepsilon_5)=0\quad\text{for all }8\text{ sign patterns.}\ }$$

So no purified weight-8 level-60 source lives on any genus-zero curve. This is the precise
content of the archive's "ζ(7) pushed off genus zero": the object is irreducibly
vector-valued / genus-two, and the paper's score formula
(`paper/sections/03_archimedean.tex`) simply does not apply to it.

---

## 3. What does exist: the level-12 parent, built and measured

`lattice/zeta7_level60/level12_parent.gp`, `level12_rates.gp` (snake, $n\le448$).
Built by the `ZETA5_TWO_ROW.md` recipe (Fricke-normalised host).

* **Host.** $x=\dfrac{\eta(4\tau)^2\eta(12\tau)^2}{\eta(\tau)^2\eta(3\tau)^2}=q+2q^2+5q^3+\cdots$,
  with $x|W_{12}=1/(16x)$; Fricke-invariant $w=x/(1+16x^2)$; normalised
  $$t=\frac{w}{1+2w}=\frac{x}{16x^2+2x+1}=q-15q^3-32q^4+138q^5+768q^6-\cdots$$
* **Rows.** $F=\Phi_{12}/(Dt)$ (weight 6), $A(t)=F$, $B(t)=F\cdot D^{-7}\Phi_{12}$:
  $$A_n:\ 1,\,-443,\,13816,\,-120335,\,1042756,\,-6945002,\,52122064,\dots$$
  $$B_n:\ 0,\,1,\,-\tfrac{57147}{128},\,\tfrac{4302754067}{279936},\,-\tfrac{6365998422665}{35831808},\dots$$

| quantity | measured ($n\le448$) | limit |
|---|---|---|
| $A_n\in\mathbb Z$ | all $n\le448$ | — |
| $B_n/A_n\to$ | `lindep` $=[-1728,209]$ **exact** | $\tfrac{209}{1728}\zeta(7)$ |
| $d_n^7B_n\in\mathbb Z$ | yes; $d_n^6$ fails at $n=2$ | $k=7$ |
| $v_\ell(B_\ell)$, $\ell=11\ldots41$ | $-7$ every time | $k=7$ **sharp** |
| $\log\lvert A_n\rvert/n$ | $2.28986$ | $\log10=2.30259$, $\lambda_1=10$ |
| $\log\lvert B_n-\xi A_n\rvert/n$ | $1.81815$ | $\log6=1.79176$, $\lambda_2=6$ |
| $\lvert B_n/A_n-\xi\rvert^{1/n}$ | $0.62393$ | $6/10=0.6$ |

The singular $t$-values are $1/58$ (image of the reciprocal pair $7\pm4\sqrt3$ under
$t=u/(4u^2+2u+4)$, $u=4x$; apparent — $A$ is regular there), and the two Fricke fixed
points $t=1/10,\,-1/6$. No recurrence of order $\le20$ with polynomial degree $\le14$
exists (mod-$p$ kernel scan, $p=2^{61}-1$, $\ge25$ excess equations): the minimal operator
is large, as expected for a degree-2 host coordinate.

**Honest score.** $\ \text{score}=\log(1/\lambda_2)-k=-\log6-7=\mathbf{-8.79}$. The linear
form *grows* like $6^n$. Even the archive's optimistic level-12 figure
$\lambda_2^{-1}=7+4\sqrt3$ gives $2.6339-7=-4.37$.

*Discrepancy noted:* the archive's level-12 singularities $\{7-4\sqrt3,\,7+4\sqrt3,\,-1\}$
do not match any integral normalisation I could construct; no Möbius map fixing $t=0$
carries $\{1/58,1/10,-1/6\}$ to that set. Its claim
$\limsup|B_n/A_n-C_7|^{1/n}\le7-4\sqrt3$ is not reproduced (measured $\to0.6$). This does
not change any verdict: $\log(7+4\sqrt3)=2.634<7$ either way.

**Normalisation rigidity.** Rescaling $t\mapsto ct$ sends $A_n\mapsto A_n/c^n$ and
$\lambda_i\mapsto\lambda_i/c$, so score $\mapsto$ score $+\log c$; integrality of the row
forces $c^n\mid A_n$, and $A_1=-443$ is prime. **$c=1$: there is no normalisation that
improves the score.**

---

## 4. The 19/6 arithmetic (Task 2, arithmetic half — verified, and correctly derived)

`lattice/zeta7_level60/entropy_test.gp`. Companion rows $B^{(S)}(t)=A(t)\cdot(D^{-7}S)(q(t))$
for each level-60 source $S$, on the integral host $t$ of §3 (a $\Gamma_0(12)\supset\Gamma_0(60)$
coordinate, so this is a legitimate arithmetic probe even though it is not an Apéry pair).

| source | depth | $\max_n\operatorname{den}\bigl(60^7 d_{\lfloor n/3\rfloor}^7 B_n\bigr)$ | sharp plain $k$ | measured $\log\operatorname{den}(B_n)/n$ at $n=230$ | entropy |
|---|---|---|---|---|---|
| $B_0$ | 3 | **1 for all $n\le230$** | 7 | 2.2785 | $7/3=2.3333$ ✓ |
| $B_1$ | 2 | huge | 7 | 3.5724 | $7/2=3.5$ ✓ |
| $B_2$ | 1 | huge | 7 | 6.8750 | $7$ ✓ |
| $\Phi_{60}^{(5)}$ | 1 | huge | 7 | 6.9000 | $7$ ✓ |
| $\Phi_{12}$ (control) | 1 | huge | 7 | 6.9000 | $7$ ✓ |

So the archive's **7/3 claim is true and now verified exactly to $n=230$**:
$60^7d_{\lfloor n/3\rfloor}^7B_n^{(B_0)}\in\mathbb Z$ for every $n$, with no exception.
Reason (elementary, and correct): $D^{-7}E_8(d\tau)=d^{-7}\sum_m\sigma_7(m)m^{-7}q^{dm}$, so
with $d\ge3$ the $n$-dependent denominator only ever sees $m\le n/3$.

The **19/6** is likewise arithmetically correct as stated: for the hypothesised filtered
triple $\bigl(D^{-7}B_0,\ D^{-2}T_1,\ D^{-1}T_2\bigr)$,
$M_N\mid\operatorname{lcm}(d_{\lfloor N/3\rfloor}^7,d_{\lfloor N/2\rfloor}^2,d_N)$ and
$\log M_N=5\psi(N/3)+\psi(N/2)+\psi(N)+o(N)\Rightarrow\tfrac53+\tfrac12+1=\tfrac{19}6$. ✓
With the depth-1 channel replaced by a depth-2 one this becomes $\tfrac73+\tfrac13+1=\tfrac{11}3$. ✓

**But the triple was never realised.** What *is* realisable is measured above, and the
decisive computation is the next one.

---

## 5. The decisive obstruction, computed (`kernels.gp`)

The three Atkin–Lehner cocycle kernel lines of the 3-dimensional purified space
(archive values, re-verified to be the $W_Q$-annihilated lines):

| line | coords in $(B_0,B_1,B_2)$ | $c_1$ | support depth | entropy | $L(\cdot,7)/\zeta(7)$ |
|---|---|---|---|---|---|
| $K_3$ | $(-1,-380,33)$ | 33 | **1** | **7** | $100529/18000$ |
| $K_4$ | $(14723,-112676,8349)$ | 8349 | **1** | **7** | $23476453/18000$ |
| $K_5$ | $(13,-52,3)$ | 3 | **1** | **7** | $6479/18000$ |
| $B_0$ | $(1,0,0)$ | 0 | 3 | $7/3$ | $2623/216000$ |

$K_3,K_4,K_5$ span the whole plane, and **none has $c_1=0$**, so none meets the depth-$\ge2$
plane $F_2$. Therefore:

> **Every monodromy-purified source at level 60 has denominator entropy exactly 7,
> and the unique entropy-$7/3$ source is monodromy-impure at all three cusps $Q=3,4,5$.**

This is the whole story in one line. The archive's programme is precisely the attempt to
buy $B_0$'s arithmetic while paying only a little of $K_Q$'s monodromy; the "Λ⁶ depth-1
lift" was the cheapest such purchase (19/6), and it is obstructed.

---

## 6. The analytic side (verified) and the honest score

`polroots` on the archive's $B_8(S)$ (critical-value polynomial of the $S=x+y$ projection
of the genus-two quotient, node factors removed) reproduces the archive exactly:

$$3125S^8+25000S^7-2065000S^6+11812000S^5+47565456S^4+55270400S^3+25351168S^2+3440640S-400384$$

moduli $\ 0.0719136,\ 0.5286989\,(\times2),\ 1.0483616\,(\times2),\ 13.4743492\,(\times2),\ 31.9419288$,
$\ \log31.9419288=3.463919528286802$ ✓ (archive: $3.463919528\ldots$).

$R=31.94$ is the **largest**: it is an analytic radius only after **seven** singular orbits
cancel. Scoreboard, $\text{score}=\log R-\delta$:

| $\delta$ | provenance | score at $R=31.94$ | score at the next honest radius $13.474$ | score with no cancellation ($0.5287$) |
|---|---|---|---|---|
| $19/6=3.1667$ | hypothesised filtered triple, **depth-1 lift obstructed** | $+0.2973$ | $-0.5659$ | $-3.8040$ |
| $11/3=3.6667$ | archive's own fallback (depth-2) | $-0.2027$ | $-1.0659$ | $-4.3040$ |
| $7/2=3.5$ | best **realisable** purified source ($F_2$, but not monodromy-purified) | $-0.0361$ | $-0.8992$ | $-4.1373$ |
| $\mathbf 7$ | **every monodromy-purified line $K_3,K_4,K_5$** | $\mathbf{-3.5361}$ | $-4.3992$ | $-7.6377$ |
| $7/3=2.3333$ | $B_0$ alone — *not monodromy-purified, no companion exists* | $+1.1306$ | $+0.2677$ | $-2.9707$ |

The only positive entries are the two rows that do not correspond to any constructible
system. Note how near the miss is on the third row: $7/2-\log31.94=0.0361$.

**No irrationality statement for $\zeta(7)$ is claimed or implied anywhere here.**

---

## 7. Files

* `lattice/zeta7_level60/source.gp` — Task 1: purity, dimensions, AL signs, filtration, $K_5$.
* `lattice/zeta7_level60/host_genus.gp` — genera of all AL quotients of $X_0(60)$.
* `lattice/zeta7_level60/allplus.gp` — the no-go: no purified AL eigenvector.
* `lattice/zeta7_level60/level12_parent.gp` — the level-12 parent row (host, $A_n$, $B_n$, $k=7$ sharp).
* `lattice/zeta7_level60/level12_rates.gp` — rates to $n=448$ + mod-$p$ recurrence scan (snake).
* `lattice/zeta7_level60/entropy_test.gp` — the $7/3$, $7/2$, $7$, $19/6$ denominator tests.
* `lattice/zeta7_level60/kernels.gp` — depth of $K_3,K_4,K_5$.
* `lattice/zeta7_level60/level12_rows.txt`, `level12_rates.log` — the rows $A_n,B_n$ to $n=448$ and the run log.
