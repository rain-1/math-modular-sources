# Cooper's three rows: twisted CM traces, the double-pole obstruction, and the master verification

*Round 2, scripts `50_`–`57_`, 2026-09-02.  All computations in PARI/GP 2.15.4 in
`lattice/cooper_congruence/round2/`.  Nothing outside that directory was touched.*

Claims are tagged **[proved]**, **[verified, range]**, **[num, digits]**, **[conjectural]**,
**[refuted, range]**.

---

## 0. Verdict table

| # | claim | status |
|---|---|---|
| 1 | The $s_7$ twisted CM trace $\tau_\chi(m)$ of §2 reproduces independently; $\mathrm{Tr}(m)=i\sqrt3\,\tau_\chi(m)$ with $\tau_\chi(m)\in\mathbf Z$ | **[num, 170 digits]**, $m\le40$ |
| 2 | For $s_7$: every SL$_2(\mathbf Z)$-class of disc $-3m^2$ (imprimitive included) has exactly one $\Gamma_0(7)$-Heegner rep with $B\equiv5m\ (14)$; counts $=\sum_{g^2\mid d}h(d/g^2)$, primitive counts $=h(d)$ | **[verified, $m\le40$]** |
| 3 | $\dfrac{m\,\beta(m)}{(m-\kappa)\,\tau_\chi(m)}\to 3$ for $s_7$ with error $O(R^{-m/2})$, $\kappa=\frac1{2\pi\operatorname{Im}\tau_0}=\frac7{\pi\sqrt3}$ | **[num]**, $m\le40$, $\bigl|(r_2-3)R^{m/2}\bigr|\le0.48$ for $15\le m\le40$ |
| 4 | $\dfrac{\beta(m)}{\tau_\chi(m)}\to 3$ **without** the $\kappa$-correction | **[refuted, $m\le40$]**: $(r_1-3)\cdot m\to-3\kappa=-3.859302267952544\ldots$ to 11 digits |
| 5 | $s_{10}$: the correct Heegner class is $\beta\equiv 6m\ (\mathrm{mod}\ 20)$, **not** $3m$ | **[proved]** ($6^2\equiv-4$, $3^2\not\equiv-4\pmod{40}$) |
| 6 | $s_{10}$: $\mathrm{Tr}_{\chi_{-4}}(m)=2i\,\tau_\chi(m)$, $\tau_\chi(m)\in\mathbf Z$, for all $m\le40$ with $5\nmid m$; untwisted trace is **not** integral | **[num, 150 digits]** |
| 7 | $s_{10}$: same $(1-\kappa/m)$ law with $\kappa=5/\pi$ and $L=2$ ($m$ odd), $L=4$ ($m$ even) | **[num]**, $m\le40$ |
| 8 | $s_{18}$: forms of disc $-36m^2$ whose **content is divisible by 3** have *no* $\Gamma_0(18)$-Heegner representative | **[proved]** (see §4.1) |
| 9 | $s_{18}$: $\mathrm{Tr}_{\chi_{-3}}(m)=\sqrt3\,\tau_\chi(m)$, $\tau_\chi(m)\in\mathbf Z$, for all **odd** $m\le29$; $\chi_{-4}$ and trivial twists are not | **[num, 140 digits]** |
| 10 | $s_{18}$: same $(1-\kappa/m)$ law with $\kappa=3/\pi$ and $L=\frac14$ ($\gcd(m,6)=1$), $L=\frac12$ ($m$ odd, $3\mid m$) | **[num]**, $m\le29$ |
| 11 | **The obstruction is uniform across all three rows**: $\beta$ is not a finite linear combination of twisted CM traces of modular *functions* | **[num, all three rows]** |
| 12 | $n^2\mid\beta(n)$ | **[verified, $n\le12000$]**, all three rows |
| 13 | $\gamma(n)=0\iff7\mid n$ ($s_7$), $\iff5\mid n$ ($s_{10}$), never ($s_{18}$) | **[verified, $n\le12000$]** |
| 14 | $c'(p^k)\equiv\psi(p)c'(p^{k-1})\pmod{p^{2k}}$ | **[verified, $p^k\le12000$, 1493 pairs/row]** |
| 15 | $c'(pq)-\psi(p)c'(q)-\psi(q)c'(p)+\psi(pq)\equiv0\pmod{p^2q^2}$ | **[verified, $pq\le12000$, 3101 pairs/row]** |
| 16 | Conjecture 4.1 in the form $c'(n)=\sum_{d\mid n}\psi(n/d)d^2\gamma(d)$ | **[verified, $n\le12000$]** (tautological given 12, recorded as a consistency check) |
| 17 | $\prod_n(1-q^n)^{\tau_\chi(n)}$ is a meromorphic modular form on $\Gamma_0(N)$ with divisor at cusps + CM points (Borcherds shape) | **[refuted]**, all three rows (§6) |

---

## 1. Set-up recap

| row | $N$ | $\tau_0$ | $D_0$ | $\nu$ | $R=e^{2\pi\operatorname{Im}\tau_0}$ | $\kappa=\frac1{2\pi\operatorname{Im}\tau_0}$ | $\psi$ |
|---|---|---|---|---|---|---|---|
| $s_7$ | 7 | $(5+\sqrt{-3})/14$ | $-3$ | 3 | $2.1756828834\ldots$ | $7/(\pi\sqrt3)=1.2864340893\ldots$ | $\mathbf 1$ |
| $s_{10}$ | 10 | $(3+i)/10$ | $-4$ | 4 | $1.8744560875\ldots$ | $5/\pi=1.5915494309\ldots$ | $\mathbf 1$ |
| $s_{18}$ | 18 | $(3+i)/6$ | $-36$ | 4 | $2.8496539082\ldots$ | $3/\pi=0.9549296585\ldots$ | $\chi_{-3}$ |

$\beta=c'\star(\mu\psi)$ as in `COOPER_CONGRUENCE.md`; $\beta(n)=n^2\gamma(n)$.

**The test function is forced in each row.**  $u|W_N=1/(Cu)$ gives
$H:=1/u-Cu$ with $H|W_N=-H$, matching $\Phi|_4W_N=-\Phi$.
* $s_7$: $X_0(7)$ has 2 cusps, exchanged by $W_7$; $u$ has divisor $(\mathfrak c_7)-(\mathfrak c_1)$, so
  $L(\mathfrak c_1+\mathfrak c_7)$ is 3-dimensional and its $W_7=-1$ eigenspace is 1-dimensional.
  $H=1/u-49u$ spans it. **[proved]**
* $s_{10}$: $u|W_2=u$, $u|W_5=1/(25u)$, so $H=1/u-25u$ is $W_2$-invariant and $W_5$-anti-invariant —
  exactly the eigencharacter of $\Phi$ ($\Phi|_4W_2=+\Phi$, $\Phi|_4W_5=-\Phi$).  The
  Atkin–Lehner group $V=\{1,W_2,W_5,W_{10}\}$ acts *simply transitively* on the four cusps of
  $X_0(10)$, so $L(\sum\mathfrak c)$ ($\dim 5$) $=\mathbf 1\oplus\mathbf Q[V]$ and the
  $(W_2{=}{+}1,W_5{=}{-}1)$ eigenspace is again **1-dimensional**; $H$ spans it.  $u$ has
  degree 2, but its divisor is $(\mathfrak c_5)+(\mathfrak c_{10})-(\mathfrak c_1)-(\mathfrak c_2)$,
  so $1/u-25u$ still has only *simple* poles at the four cusps and no pole in $\mathbf H$.
  **A Hauptmodul is therefore not needed.** **[proved]**
* $s_{18}$: $u|W_9=u$, $u|W_{18}=1/u$, so $H=1/u-u$ is $W_9$-invariant, $W_{18}$- and
  $W_2$-anti-invariant, matching $\Phi|_4W_9=+\Phi$, $\Phi|_4W_{18}=-\Phi$.  $u$ has degree 4
  with divisor $\mathfrak c_2+\mathfrak c_3+\mathfrak c_{3'}+\mathfrak c_{18}-\mathfrak c_1-\mathfrak c_6-\mathfrak c_{6'}-\mathfrak c_9$,
  so $H$ has simple poles at all eight cusps. **[proved]**

---

## 2. Task 1 — $s_7$ reproduced and pushed to $m=40$  (`52_s7.gp`, `52_s7.log`)

**Heegner audit (audit A/B).** For every $m\le40$, with $d=-3m^2$ and $\beta=5m\bmod14$:

* `#redforms(d)` $=\sum_{g^2\mid d,\ 4\mid d/g^2\ \mathrm{or}\ d/g^2\equiv1(4)}h(d/g^2)$ — exact match, $m\le40$;
* the number of *primitive* reduced forms equals `qfbclassno(d)` — exact match, $m\le40$;
* $\sum_Q1/\omega_Q=H(3m^2)$ (`qfbhclassno`) — exact match, $m\le40$;
* **every** reduced form (imprimitive included) has a $\Gamma_0(7)$-Heegner representative with
  $7\mid A$, $B\equiv5m\ (14)$, found with search bound $|p|,|r|\le80$: `#found = #forms` in all 40 cases.

**The trace.** $\mathrm{Tr}(m)=\sum_Q\chi_{-3}(Q)H(\alpha_Q)/\omega_Q=i\sqrt3\,\tau_\chi(m)$ with
$\tau_\chi(m)\in\mathbf Z$ to **170 digits** (residual $\le10^{-159}$ at $m=40$).  Values (`52_tau_s7.txt`):

```
1, -6, 5, 18, -57, 90, 0, -458, 1227, -1138, -2609, 12570, -22061, 0, 104595,
-283366, 274837, 596598, -2922599, 5100094, 0, -24134186, 65485023, -63419650,
-137954013, 674437686, -1176694173, 0, 5569992617, -15111679490, 14632086955,
31835062406, -155631144085, 271539055858, 0, -1285358125434, 3487216342321,
-3376563966126, -7346338747345, 35914151811946
```

This **agrees exactly** with the values you reported for $m\le20$, and fills the two gaps:
at $7\mid m$ the two $\beta$-classes $\pm5m$ merge, $\mathrm{Tr}(m)$ becomes **real** (not
$i\sqrt3\cdot$ real), hence $\tau_\chi(m)=0$ in the $i\sqrt3$ normalisation.  E.g. $m=7$:
$\mathrm{Tr}=-\sqrt3\cdot93.4221\ldots$, $m=14$: $-\sqrt3\cdot111.7274\ldots$.

**The ratio test.**  With $r_1(m)=\beta(m)/\tau_\chi(m)$ and $r_2(m)=\dfrac{m\,\beta(m)}{(m-\kappa)\tau_\chi(m)}$:

| $m$ | $r_1$ | $(r_1-3)\cdot m$ (target $-3\kappa=-3.8593022679525443$) | $r_2$ | $(r_2-3)R^{m/2}$ |
|---|---|---|---|---|
| 27 | 2.857062847884095 | $-3.859302107129434$ | $2.99999996736442962$ | $-0.00118$ |
| 33 | 2.883051446437614 | $-3.859302267558730$ | $3.00000000001241784$ | $+4.62\cdot10^{-6}$ |
| 37 | 2.895694532867466 | $-3.859302283903774$ | $2.99999999955335657$ | $-0.00079$ |
| 39 | 2.901043531594738 | $-3.859302267805231$ | $3.00000000000390610$ | $+1.50\cdot10^{-5}$ |
| 40 | 2.903517456907184 | $-3.859301723712622$ | $3.00000001405811914$ | $+0.079$ |

**Answer to the question asked:** yes, $r_2\to3$ *exactly*, with error $O(R^{-m/2})$:
$\sup_{15\le m\le40}|(r_2-3)R^{m/2}|=0.473$ (attained at $m=18$), and the quantity does not
drift.  By contrast $(r_1-3)m\to-3\kappa$ to **11 significant digits**, i.e. the $1/m$ term is
present and is *exactly* $-3\kappa/m$ with $\kappa=1/(2\pi\operatorname{Im}\tau_0)$.
$|(r_1-3)R^{m/2}|$ over the last ten clean $m$ reaches $5.45\cdot10^5$ — the naive claim is
refuted by 6 orders of magnitude.

---

## 3. Task 2a — $s_{10}$  (`53_s10.gp`, `53_s10.log`)

**Correction to the task's $\beta$.**  For $Q=[A,B,C]$ with $N\mid A$, $B^2\equiv d\pmod{4N}$;
here $4N=40$ and $d=-4m^2$.  $\;(6m)^2+4m^2=40m^2\equiv0$, while $(3m)^2+4m^2=13m^2\not\equiv0$.
So $\beta\equiv6m\pmod{20}$ is the right class; $3m$ is **not** a square root of $-4m^2$ mod 40.
(Directly: $\tau_0=(3+i)/10$ is $\alpha_Q$ for $Q=[10,-6,1]$, $\beta\equiv-6\equiv14$;
$(7+i)/10$ is $[10,-14,5]$, $\beta\equiv6$; the two elliptic points sit in the two classes $\pm6$.)
The valid $\beta$ mod 20 are listed in the log for $m\le6$; for $5\mid m$ there is only **one**
class ($\beta\equiv10$), which is the anomaly.

**Heegner audit.**  `#redforms(-4m^2)` $=\sum_g h$ and `#prim` $=h(-4m^2)$ for all $m\le40$;
**every** reduced form has a $\Gamma_0(10)$-Heegner rep with $10\mid A$, $B\equiv6m\ (20)$
(search bound 80): `#found = #forms` in all 40 cases.

**The trace.**  $\mathrm{Tr}_{\chi_{-4}}(m)=\sum_Q\chi_{-4}(Q)H(\alpha_Q)/\omega_Q = 2i\,\tau_\chi(m)$
with $\tau_\chi\in\mathbf Z$ to 150 digits, for every $m\le40$ with $5\nmid m$.  The **untwisted**
trace is not integral (e.g. $m=6$: $11.2251\ldots-130.8796\ldots i$), so the genus character
$\chi_{-4}$ is the right twist — the exact analogue of $\chi_{-3}$ for $s_7$, and in both rows
the algebraic factor is $\sqrt{D_0}$.  For $5\mid m$ the classes merge and $\mathrm{Tr}$ is real
(e.g. $m=5$: $\mathrm{Tr}=36$).  Values (`53_tau_s10.txt`, 0 = anomalous):

```
2, -3, -10, 14, 0, -33, 94, 76, -542, 0, 1910, -1078, -4134, 6243, 0, -22184,
51170, 48285, -290970, 0, 1022314, -593169, -2220194, 3367780, 0, -11827335,
27408854, 25681474, -155818566, 0, 547481758, -317100624, -1188867022,
1802849133, 0, -6334601378, 14676878282, 13755727527, -83439532098, 0
```

**The ratio test** (`55_ratios.gp`).  The limit is **parity-dependent**:
$L=2$ for $m$ odd, $L=4$ for $m$ even.  (Cause: for even $m$ a positive proportion of the forms
of disc $-4m^2$ have $A$, $C$ and $A+B+C$ all even, so $\chi_{-4}(Q)=0$ and the trace is
halved.)  With that $L$:

| $m$ | $L$ | $(r_1-L)\cdot m$ vs $-L\kappa$ | $r_2$ | $(r_2-L)R^{m/2}$ |
|---|---|---|---|---|
| 29 | 2 | $-3.1830971734$ vs $-3.1830988618$ | $2.00000006160231673$ | $+0.00056$ |
| 34 | 4 | $-6.3661085589$ vs $-6.3661977237$ | $4.00000275128160841$ | $+0.120$ |
| 37 | 2 | $-3.1830985860$ vs $-3.1830988618$ | $2.00000000778887912$ | $+0.00087$ |
| 38 | 4 | $-6.3662649069$ vs $-6.3661977237$ | $3.99999815473545847$ | $-0.282$ |
| 39 | 2 | $-3.1830985656$ vs $-3.1830988618$ | $2.00000000791852303$ | $+0.0017$ |

$\sup$ over the last 10 clean $m$: $|(r_2-L)R^{m/2}|\le0.393$ **with** the $\kappa$-correction,
$|(r_1-L)R^{m/2}|$ up to $2.56\cdot10^4$ **without** it.

---

## 4. Task 2b — $s_{18}$  (`54_s18.gp`, `54_s18.log`)

### 4.1 A genuine obstruction in the Heegner correspondence  **[proved]**

$\tau_0=(3+i)/6$ is $\alpha_Q$ for $Q=[18,-18,5]$ (disc $-36$), so $\beta\equiv18\pmod{36}$;
the phase requirement $e^{2\pi i\beta/36}=(-1)^m$ (matching $q_0^{-m}=(-1)^mR^m$) forces
$\beta\equiv18m\pmod{36}$.

**Not every SL$_2(\mathbf Z)$-class of disc $-36m^2$ has a $\Gamma_0(18)$-Heegner representative.**
If $3\mid\mathrm{cont}(Q)$, write $Q=3Q'$; then $18\mid Q(p,r)=3Q'(p,r)$ needs $6\mid Q'(p,r)$,
and $Q'$ has discriminant $-4(m/1)^2/\ldots$ of the shape $\square+\square$; $3\mid p^2+r^2$
forces $3\mid p$ and $3\mid r$, contradicting $\gcd(p,r)=1$.  Concretely $[3,0,3]$ (disc $-36$)
has no rep, and the first run of `54_s18.gp` flagged exactly the content-divisible-by-3 forms
(20 `NO REP` events for $m\le30$) — this is **not** a search-bound artefact.

The admissible set is therefore $\{Q:3\nmid\mathrm{cont}(Q)\}$, and with that restriction
`#found = #admissible` for **all** $m\le30$ (search bound 80).  All `#redforms` counts match
$\sum_gh(d/g^2)$.

### 4.2 The right twist

Of the three twists tried — trivial, $\chi_{-4}$ (writing $-36m^2=(-4)\cdot 9m^2$) and
$\chi_{-3}$ (writing $-36m^2=(-3)\cdot 12m^2$) — only $\chi_{-3}$ gives an algebraic answer:
$$\mathrm{Tr}_{\chi_{-3}}(m)=\sqrt3\,\tau_\chi(m),\qquad \tau_\chi(m)\in\mathbf Z,$$
to 140 digits, for every **odd** $m\le29$ (plus, accidentally, $m=2$).  This is a pleasing
independent confirmation that $\psi=\chi_{-3}$ is the character of the row.  Values
(`54_tau_s18.txt`, indices 1..30, 0 = not integral):

```
16, -12, 48, 0, 880, 0, 7056, 0, 28656, 0, 464816, 0, 3774160, 0, 15324240, 0,
248882992, 0, 2021062224, 0, 8206044336, 0, 133274809712, 0, 1082261746000, 0,
4394268082224, 0, 71367548709616, 0
```

All entries are divisible by 16.  **Even $m$ are anomalous**: $\beta=18m\equiv0\pmod{36}$, so the
$\pm\beta$ classes merge (exactly as $7\mid m$ for $s_7$ and $5\mid m$ for $s_{10}$, except that
here it happens for *half* of all $m$), and the trace lands outside $\sqrt3\,\mathbf Q$.

### 4.3 The ratio test

$L=\tfrac14$ when $\gcd(m,6)=1$, $L=\tfrac12$ when $m$ is odd and $3\mid m$.

| $m$ | $L$ | $(r_1-L)\cdot m$ vs $-L\kappa$ | $r_2$ | $(r_2-L)R^{m/2}$ |
|---|---|---|---|---|
| 23 | 1/4 | $-0.2387324104964392$ vs $-0.2387324146378430$ | $0.250000000187860768$ | $+3.19\cdot10^{-5}$ |
| 25 | 1/4 | $-0.2387324141594486$ vs $-0.2387324146378430$ | $0.250000000019895736$ | $+9.63\cdot10^{-6}$ |
| 27 | 1/2 | $-0.4774648291617471$ vs $-0.4774648292756860$ | $0.500000000004374684$ | $+6.03\cdot10^{-6}$ |
| 29 | 1/4 | $-0.2387324146803623$ vs $-0.2387324146378430$ | $0.249999999998483895$ | $-5.96\cdot10^{-6}$ |

$(r_1-L)m\to-L\kappa$ to **13 significant digits** — the cleanest of the three rows.
$\sup$ over the last 10 clean $m$: $0.0028$ with the correction, $3.24\cdot10^4$ without.

---

## 5. The conclusion asked for: **the obstruction is uniform** across all three rows

| row | $\kappa$ | $L$ (per class of $m$) | $\sup_{\text{last }10}\bigl|(r_2-L)R^{m/2}\bigr|$ | $\sup_{\text{last }10}\bigl|(r_1-L)R^{m/2}\bigr|$ |
|---|---|---|---|---|
| $s_7$ | $7/(\pi\sqrt3)$ | $3$ | $0.290$ | $5.45\cdot10^5$ |
| $s_{10}$ | $5/\pi$ | $2$ (odd $m$), $4$ (even $m$) | $0.393$ | $2.56\cdot10^4$ |
| $s_{18}$ | $3/\pi$ | $1/4$ ($\gcd(m,6)=1$), $1/2$ ($3\mid m$ odd) | $0.0028$ | $3.24\cdot10^4$ |

**Statement.**  In every row, numerically,
$$\frac{\beta(m)}{\tau_\chi(m)}=L\Bigl(1-\frac{\kappa}{m}\Bigr)+O(R^{-m/2}),\qquad
\kappa=\frac1{2\pi\operatorname{Im}\tau_0},$$
with the $1/m$ coefficient equal to $-L\kappa$ to 11–13 significant digits.  **[num]**

**Why this kills the CM-trace identification.**  From `cooper_sources/REPORT.md` §2.4,
$c(m)=-\nu^2\bigl(q_0^{-m}/g'(u_0)+\overline{\phantom{x}}\bigr)(m-\kappa)+O(R_2^m)$, so
$\beta(m)\sim c(m)/m$ carries the factor $(1-\kappa/m)$: this is precisely the signature of a
**double** pole of $\Phi$ in $\tau$ (a simple pole would give no $1/m$, a triple pole would give a
$m^2$-type factor).  A twisted CM trace of a weight-0 *modular function* $f$ over the Heegner
points of discriminant $D_0m^2$ is
$\mathrm{Tr}_f(m)=\text{const}\cdot R^m\cdot(\text{oscillation})+O(R^{m/2})$ with **no** $1/m$
term — its dominant Heegner point contributes $f$'s $q$-expansion principal part, which has no
$m$-dependence beyond $q^{-m}$.  Any finite linear combination of such traces (over finitely
many discriminants and finitely many functions) inherits this.  Hence:

> **$\beta$ is not a finite $\mathbf Q$-linear combination of twisted CM traces of modular
> functions — in any of the three rows.**  The deficit is exactly $\kappa/m$ with
> $\kappa=1/(2\pi\operatorname{Im}\tau_0)$, i.e. it is the *residue* term $A_1=iA_2/\operatorname{Im}\tau_0$
> of `cooper_sources/REPORT.md` §2.3, not an accident of normalisation. **[num, all three rows]**

What *would* fit is a trace of a **weight-0 real-analytic** object with a $1/\operatorname{Im}\tau$
term (a Maass–Poincaré / harmonic-weak-Maass shadow), or equivalently a trace over the Heegner
points of the non-holomorphic completion; that is the natural next target and is consistent with
`COOPER_CONGRUENCE.md` §4's Shimura–Borcherds reading, which needs a weight-$5/2$ (not $3/2$)
input.

---

## 6. Task 4 — the Borcherds-product test  (`56_borcherds.gp`, `56_borcherds.log`)

**Method.**  $P:=\prod_{n\ge1}(1-q^n)^{\tau_\chi(n)}$; its logarithmic derivative is exact,
$$D\log P=-\sum_{M\ge1}\Bigl(\sum_{d\mid M}d\,\tau_\chi(d)\Bigr)q^M ,$$
so only $\tau_\chi(d)$ for $d\le M$ is needed.  If $P$ is (a rational power of) a meromorphic
modular form on $\Gamma_0(N)$ whose divisor is supported on the cusps and on the CM points
$g(u)=0$ (the only divisor a Borcherds lift of these traces could have), then
$P=c\,q^h\cdot(\text{eta quotient})\cdot u^a g(u)^b$, whence
$$D\log P\in\operatorname{span}_{\mathbf Q}\Bigl\{1\Bigr\}+\operatorname{span}_{\mathbf Q}\bigl\{E_2(d\tau):d\mid N\bigr\}
+\operatorname{span}_{\mathbf Q}\bigl\{F\,u^i/g(u)^j\bigr\},$$
using $D\log q^h=h$, $D\log\eta(d\tau)=\tfrac d{24}E_2(d\tau)$, $D\log R(u)=F\,uR'(u)/R(u)$.
The test is exact linear algebra over $\mathbf Q$ on the $q$-coefficients $q^0,\dots,q^M$
($i=-2..4$, $j=0,1,2$; basis size 23/25/27, equations 41/41/31).

**Positive control.**  With $\tau_\chi(n)\equiv24$ (so $P=\eta^{24}/q$) the routine recovers
$D\log P=E_2(\tau)-1$ exactly, in all three rows.

**Result. [refuted]**  For the actual $\tau_\chi$ the kernel contains **no** vector with a
non-zero $D\log P$ component, in all three rows.  So
$$\prod_n(1-q^n)^{\tau_\chi(n)}\ \text{is not a Borcherds-type product on }\Gamma_0(N)$$
(not even up to a rational power or a multiplier system, both of which the log-derivative test is
blind to).  For $s_{10}$ and $s_{18}$ the exponents at the anomalous $m$ ($5\mid m$, resp. $m$
even) were set to $0$; that is forced for $s_7$ (where $\tau_\chi(7m)=0$ genuinely) and is a
*choice* for the other two — flagged as a caveat, but the $s_7$ result alone is unconditional.

The first log-derivative coefficients $\sum_{d\mid M}d\tau_\chi(d)$ are
$s_7$: $1,-11,16,61,-284,544,1,-3603,\ldots$;
$s_{10}$: $2,-4,-28,52,2,-232,660,660,\ldots$;
$s_{18}$: $16,-8,160,-8,4416,136,49408,-8,\ldots$

---

## 7. Task 3 — the master verification  (`51_master.gp`/`.log` at $n\le3000$; `57_master_big.gp`/`.log` at $n\le12000$)

Total runtime for all three rows at $n\le12000$: **7.4 s** (the eta quotients are computed with
PARI's `eta(q^d+O(q^M))`, which is the pentagonal-number-theorem sparse product; `50_lib.gp`
is a drop-in replacement for `lib.gp` and was checked to agree with it and with the stored data
files coefficient-by-coefficient to $n=300$, `50_sanity.log`).

### 7.1 $n^2\mid\beta(n)$  **[verified, $n\le12000$, exact, all three rows]**
(Previously $n\le1500$.)

### 7.2 Zeros of $\gamma$  **[verified, $n\le12000$]**
$$\gamma_{s_7}(n)=0\iff 7\mid n\quad(1714\text{ values}),\qquad
\gamma_{s_{10}}(n)=0\iff 5\mid n\quad(2400\text{ values}),\qquad
\gamma_{s_{18}}(n)\ne0\ \forall n .$$
Verified to be *exactly* the sets of multiples, not merely contained in them.  This is the
$\gamma$-side shadow of Theorem 2.2 (`COOPER_CONGRUENCE.md`): $\Xi|U_7=\Xi$ and $\Xi|U_5=\Xi$
identically force $\gamma(pd)=0$ through Conjecture 4.1, and conversely.

### 7.3 Sharpness $v_p(\gamma(n))=v_p(\beta(n))-2v_p(n)$, $p\le60$, $n\le12000$

Note $v_p(\beta(n))-2v_p(n)=v_p(\gamma(n))$ identically, so the "strengthened statement"
$v_p(\beta(n))\ge2v_p(n)$ and $n^2\mid\beta(n)$ are the same assertion; the sharpness data is:

* $A_p:=\min_{n,\ \gamma(n)\ne0}v_p(\gamma(n))=\mathbf 0$ for **every** $p\le60$ and every row
  (attained at $n=1$, $\gamma(1)=1$).  So the exponent $2$ in $n^2\mid\beta(n)$ is exactly
  sharp at every prime: nothing more is true.
* $B_p:=\min_{p\mid n,\ \gamma(n)\ne0}v_p(\gamma(n))$ — the interesting refinement, i.e. how much
  more than $p^2$ divides $\beta(n)$ at multiples of $p$:

| row | $B_2$ | $B_3$ | $B_5$ | $B_7$ | $B_p$, $11\le p\le59$ |
|---|---|---|---|---|---|
| $s_7$ | **1** (at $n=2$) | 0 | **1** (at $n=5$) | $+\infty$ ($\gamma$ vanishes) | $0$ for every $p$ |
| $s_{10}$ | 0 | 0 | $+\infty$ | 0 | $0$ for every $p$ |
| $s_{18}$ | 0 | 0 | 0 | 0 | $0$ for every $p$ |

The two $B_p=1$ entries for $s_7$ at $p=2,5$ are exactly the "$p^3$ instead of $p^2$" cells
noted in `cooper_sources/REPORT.md` §0.3; they persist to $n\le12000$.  Everything else is
sharp at $p^2$.

### 7.4 The refined supercongruences of `COOPER_CONGRUENCE.md` §4

* $c'(p^k)\equiv\psi(p)c'(p^{k-1})\pmod{p^{2k}}$: **no failures**, $p^k\le12000$,
  **1493 $(p,k)$ pairs per row** (previously $p\le43$, $k\le3$).
* $c'(pq)-\psi(p)c'(q)-\psi(q)c'(p)+\psi(pq)\equiv0\pmod{p^2q^2}$: **no failures**,
  **3101 pairs $p<q$ with $pq\le12000$ per row**.
* The exact modulus actually reached ($v_p$ of the difference vs the required $2k$), a
  by-product worth recording — it is *not* always sharp:

| row | $p$ | $v_p(c'(p^k)-\psi(p)c'(p^{k-1}))$ for $k=1,2,3,\ldots$ | required |
|---|---|---|---|
| $s_7$ | 2 | $3,5,7,9,11,13,\ldots=2k+1$ | $2k$ |
| $s_7$ | 3 | $2,5,8,11,14,\ldots=3k-1$ | $2k$ |
| $s_7$ | 5 | $3,5,7,9,11=2k+1$ | $2k$ |
| $s_7$ | 7 | $+\infty$ (exact eigenvalue) | $2k$ |
| $s_{10}$ | 2 | $2,5,8,11,\ldots=3k-1$ | $2k$ |
| $s_{10}$ | 5 | $+\infty$ | $2k$ |
| $s_{18}$ | 2 | $2,6,7,9,11,\ldots$ | $2k$ |
| $s_{18}$ | 23 | $2,5$ | $2k$ |
| all | other $p\le53$ | exactly $2k$ | $2k$ |

So the modulus $p^{2k}$ is sharp for every $p$ except $p=2,3,5$ (and $p=7\mid s_7$, $5\mid s_{10}$),
where $v_p$ grows like $3k-1$ or $2k+1$ — a genuinely stronger congruence at the small primes,
uniform in $k$ up to $k=13$.

### 7.5 Conjecture 4.1 consistency
$c'(n)=\sum_{d\mid n}\psi(n/d)d^2\gamma(d)$ holds for all $n\le12000$ in all three rows
(this is equivalent to 7.1 by Möbius inversion; recorded as a check on the code path).

---

## 8. Honest caveats

1. **$s_{18}$, even $m$.**  Half the family is anomalous ($\beta\equiv0\bmod36$); I have *no*
   integral trace there, and hence no ratio data.  The odd-$m$ evidence is very strong
   (13 digits), but the even-$m$ half is a gap.  Trying $\beta\in\{12,24\}\bmod36$ for even $m$
   is possible but the phase $e^{2\pi i\cdot12/36}\ne\pm1$ does not match $q_0^{-m}=(-1)^mR^m$,
   so the ratio would not be real; I did not pursue it.
2. **$s_{18}$, content-3 forms.**  These are genuinely absent from the $\Gamma_0(18)$-Heegner
   set (proved), so my $\tau_\chi$ is a trace over a *proper subset* of $Q_{-36m^2}/\mathrm{SL}_2(\mathbf Z)$.
   That is the correct Heegner set for level 18, but it means the Gross–Kohnen–Zagier bijection
   does **not** hold verbatim at this level/discriminant pair; anyone quoting a "twisted trace"
   theorem here must check the hypotheses.
3. **The constants $L$** are empirical (2, 3, 4, 1/4, 1/2 as tabulated) and I did not derive them
   from $\nu^2/g'(u_0)$ and the count of dominant Heegner points.  Nothing in the conclusion
   depends on their values — only on their being *constant in $m$ within each class*, which is
   what the $(r_2-L)R^{m/2}$ columns establish.
4. **The Borcherds test** allows poles of order $\le2$ at the CM points and $u$-exponents in
   $[-2,4]$.  As argued in §6 this is complete for divisors supported on cusps $\cup$ CM points
   (the divisor is forced to be $a\,\mathrm{div}(u)+b\,\mathrm{div}(g(u))+$ cusps), but if a
   Borcherds lift of these traces had divisor at *other* Heegner points the test would not see it.
5. **Precision.**  All CM evaluations at `realprecision` 140–160; the integrality residuals are
   $\le10^{-140}$ throughout, so the integer identifications are not in doubt.  All divisibility
   claims are exact integer arithmetic.

---

## 9. Files

| file | content |
|---|---|
| `50_lib.gp` | fast drop-in for `lib.gp` (PARI `eta` on a series) |
| `50_sanity.gp/.log` | agreement with `lib.gp` and with the stored data, $n\le300$ |
| `51_master.gp/.log` | Task 3 at $n\le3000$ (writes `50_beta_*.txt`, `50_gamma_*.txt`, `50_cp_*.txt`) |
| `52_s7.gp/.log`, `52_tau_s7.txt` | Task 1: $s_7$ audits + traces, $m\le40$ |
| `53_s10.gp/.log`, `53_tau_s10.txt` | Task 2a: $s_{10}$, $m\le40$ |
| `54_s18.gp/.log`, `54_tau_s18.txt` | Task 2b: $s_{18}$, $m\le30$ |
| `55_ratios.gp/.log` | consolidated ratio test, all three rows |
| `56_borcherds.gp/.log` | Task 4, with positive control |
| `57_master_big.gp/.log` | Task 3 extended to $n\le12000$ |

---

# ADDENDUM (same session): the $\widehat f$ identification for $s_{10}$ and $s_{18}$

*After the above was written, the parallel session identified $\gamma_{s_7}$ as a twisted CM
trace of $\widehat f=Df+f/(2\pi y)$, $f=1/(xF)$ the weight-$(-2)$ weakly holomorphic form with
$f=q^{-1}+O(1)$.  That confirms §5 above from the other side: $\widehat f$'s leading term at the
top Heegner point is $q^{-1}(1-1/(2\pi y))=q^{-1}(1-\kappa/m)$ — exactly the $1/m$ tail that
§5 measured and that no modular **function** can supply.  Scripts `58_`, `59_`.*

## A. Row $s_{10}$: the factor 2 is the Atkin–Lehner involution $W_2$

**A.1 The identity. [num, 130 digits, $m\le30$]**  With $d=-4m^2$, $\beta\equiv6m\ (20)$,
$T(m)=\sum_Q\chi_{-4}(Q)\,\widehat f(\alpha_Q)/\omega_Q$ over the $\Gamma_0(10)$-Heegner classes:
$$\beta_{s_{10}}(m)=i\,T(m)\ (m\ \text{odd}),\qquad \beta_{s_{10}}(m)=2i\,T(m)\ (m\ \text{even}),$$
verified to **130 digits** for every $m\le30$ with $5\nmid m$ (`59_final2.log`), with
`#classes found = #reduced forms` in every case.

**A.2 The mechanism.**  $\Phi|_4W_2=+\Phi\Rightarrow f|W_2=+f\Rightarrow\widehat f$ is
$W_2$-invariant.  I computed the $W_2$-permutation of the Heegner set explicitly
($W_2=\binom{6\ 1}{10\ 2}$, det 2; `59_w2s10.gp`, `.log`):

* **$m$ odd:** $W_2$ preserves the $\beta$-class, acts as an involution on the set of classes,
  and **preserves both $\chi_{-4}$ and $\widehat f$** (checked term by term), with **exactly one
  fixed class** — verified $m=1,3,7,9,11,13,17,19$.  So every free $W_2$-orbit is counted twice
  in $T$.  Consistently, $\beta/(i\,T_{W_2})\to2$ where $T_{W_2}$ is the sum over a
  $W_2$-transversal: $2.0137,\,2.00225,\,2.00061,\,2.00028,\,1.999978,\,1.9999962$ at
  $m=7,9,11,13,17,19$ (the defect is the fixed class, whose term is $O(1)$).
* **$m$ even:** $\chi_{-4}$ **vanishes** exactly on the even-content forms (which exist only for
  even $m$: content 2 needs $-m^2\equiv0,1\bmod4$), and those are exactly the classes whose
  $\Gamma_0(10)$-Heegner fibre is **3, not 1** — the Gross–Kohnen–Zagier bijection genuinely
  fails there (`59_fibre10.log`: three distinct $\widehat f$ values over one $\mathrm{SL}_2(\mathbf Z)$-class,
  for every even-content class at $m=12,14,16$; every content-odd class has fibre 1).  The
  $W_2$-doubling therefore does not occur and $T$ is already "correct".

**A.3 Consequence for $\lambda$.**  $|\lambda|=\nu^2/|g'(u_0)|=16/8=2$ is the value realised for
even $m$; for odd $m$ the sum over $\mathrm{SL}_2(\mathbf Z)$-classes double-counts $W_2$-orbits,
which is why $|\lambda|$ reads as 1 there.  So the parallel session's revised guess (that the
**odd**-$m$ traces are the ones that are twice too large) is **correct**, and the cause is $W_2$.
*Caveat: the factor is exactly 2 empirically; $T=2S_{\text{free}}+t_{\text{fix}}$ is not exactly
$2T_{W_2}$, so "the trace on $X_0(10)/W_2$" is the mechanism but not yet a clean restatement.*

## B. Row $s_{18}$: two fixes, then an exact identity

**B.1 The admissible set. [proved]**  Forms of discriminant $-36m^2$ with $3\mid\mathrm{cont}(Q)$
have **no** $\Gamma_0(18)$-Heegner representative: $Q=3Q'$ needs $6\mid Q'(p,r)$, and
$3\mid p^2+r^2$ forces $3\mid p,3\mid r$, contradicting $\gcd(p,r)=1$.  (First instance:
$[3,0,3]$ at $m=1$.)  With the admissible set $\{3\nmid\mathrm{cont}\}$, **all** classes are
found for every $m\le30$.  This is the whole of the class-count mismatch reported earlier.

**B.2 A bug in `heeg.gp`. [proved]**  `heegrep` tests only **one** completion $(q,s)$ of a first
column $(p,r)$.  The general matrix is $\binom{p\ \ q+tp}{r\ \ s+tr}$, which shifts $B$ by $2At$;
so valid representatives (in particular the ones with **smallest $A$**, i.e. largest
$\operatorname{Im}\alpha$) are silently discarded.  Fixed in `heegmin2` (`59_final3.gp`) by solving
$B_0+2At\equiv\beta\pmod{2N}$ for $t$.  Choosing the **minimal-$A$** representative is what makes
the evaluation of $\widehat f$ numerically safe; the earlier blow-ups (values $\sim10^{36}$) were
`heegrep` returning e.g. $A=5634$ instead of $A=18$ for $[1,0,9]$, so that
$\operatorname{Im}\alpha=3m/A$ was tiny and both routes for $\widehat f$ failed.

**B.3 The genus character must be corrected on imprimitive classes. [num, 135 digits]**
With $\chi_{-3}$ alone the identity failed at $m=7,11,19,21,23$ — and *only* there among odd
$m\le29$; those are exactly the $m$ having a prime factor $p\equiv3\ (4)$ other than 3, i.e. a
prime **inert in $\mathbf Q(i)$** other than the conductor.  The defect is exactly
$\pm8\sqrt3$ (and $-72\sqrt3$ at $m=21$), and `59_s18diag.log` localises it: it is exactly
**twice the content-$p$ part** of the trace.  Since $-36$ is **not fundamental**,
$-36=(-3)\cdot12=(-4)\cdot9$, and the character that works is
$$\boxed{\ \chi^*(Q)=\chi_{-3}(Q)\cdot\Bigl(\tfrac{-4}{\mathrm{cont}(Q)}\Bigr)\ }$$
i.e. the $\chi_{-3}$ genus character of the form times the $\chi_{-4}$ symbol of its content
(trivial on primitive forms, hence invisible at small $m$).  With it, and with
$T^*(m)=\sum_Q\chi^*(Q)\widehat f(\alpha_Q)/\omega_Q$ over the admissible Heegner classes,
$\beta\equiv18m\ (36)$:
$$\beta_{s_{18}}(m)=-\frac{1}{4\sqrt3}\,T^*(m)\quad(3\nmid m),\qquad
\beta_{s_{18}}(m)=-\frac{1}{2\sqrt3}\,T^*(m)\quad(3\mid m),$$
**verified to 135 digits for every odd $m\le29$** (`59_s18fix.log`).  Even $m$ remain anomalous:
$\beta\equiv18m\equiv0\ (36)$, so the $\pm\beta$ classes merge (the $s_{18}$ analogue of $7\mid m$
for $s_7$ and $5\mid m$ for $s_{10}$) — half the family, and still open.

**B.4 $\lambda$ versus the polar coefficient. [refuted for $s_{18}$]**
$|\lambda_{s_{18}}|=\frac1{4\sqrt3}$ (resp. $\frac1{2\sqrt3}$), whereas
$\nu^2/|g'(u_0)|=2/\sqrt3$.  The polar rule is therefore off by a factor **8** (resp. **4**) for
$s_{18}$, although it is exact for $s_7$ ($\sqrt3$) and for $s_{10}$ on even $m$ (2).  The
discrepancy is a power of 2 and depends on $3\mid m$ exactly as $s_{10}$'s depends on $2\mid m$,
so it is presumably the same phenomenon (extra dominant Heegner classes produced by the
Atkin–Lehner group and by the non-fundamental conductor), but I have **not** derived it.

## C. Summary of the constants

| row | $\nu^2/|g'(u_0)|$ | observed $\lambda$ | ratio |
|---|---|---|---|
| $s_7$ | $\sqrt3$ | $i\sqrt3$ | 1 |
| $s_{10}$, $m$ even | 2 | $2i$ | 1 |
| $s_{10}$, $m$ odd | 2 | $i$ | 2 ($W_2$) |
| $s_{18}$, $3\nmid m$ odd | $2/\sqrt3$ | $-1/(4\sqrt3)$ | 8 |
| $s_{18}$, $3\mid m$ odd | $2/\sqrt3$ | $-1/(2\sqrt3)$ | 4 |

## D. Addendum file list

`58_diag10.gp/.log` (per-form s10 diagnostic), `59_fibre10.gp/.log` (Heegner fibre test, s10),
`59_w2s10.gp/.log` ($W_2$ action), `59_s18ident.gp/.log`, `59_s18b.gp/.log` (min-$A$ fix),
`59_fibre18.gp/.log`, `59_final.gp/.log`, `59_final2.gp/.log` (s10 to 130 digits),
`59_final3.gp/.log` (corrected rep search), `59_s18diag.gp/.log` (defect localisation),
`59_s18fix.gp/.log` (**the corrected s18 identity**).
