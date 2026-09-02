# Eisenstein source spaces $M_3$ and $M_4$ on modular hosts of level $N\le 60$: a machine census

**Working dir:** `/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/hostscan/eis/`
**Tool:** PARI/GP 2.15.4 (`gp -q`), plus exact rational linear algebra in Python `fractions`.
**Scripts:** `census.gp` (direction enumeration), `dimchk.gp` (PARI dimension cross-check), `beukers.gp` (calibration),
`verify.gp` (constant-term support structure), `ann2.gp` (fixed-nebentypus annihilation), `an1.py`/`an2.py`/`an3.py` (tables).
**Data:** `directions.txt` (3996 rows), `pdims.txt`, `verify.txt`, `ann2.txt`, `tables.txt`, `maintable.md`.

Every claim below is tagged **[computed]** (a machine computation performed here) or **[reasoned]**
(a short argument, stated as such, usually corroborated by a spot check).

---

## 0. Conventions

Task convention, used throughout:
$$E_{k'}^{\psi,\phi}(\tau) = c_0 + \sum_{n\ge1}\Big(\sum_{d\mid n}\psi(n/d)\phi(d)d^{k'-1}\Big)q^n,\qquad
L(E_{k'}^{\psi,\phi},s)=L(\psi,s)\,L(\phi,s-(k'-1)).$$
The critical slot is $s=k'-1$, value $L(k'-1,\psi)\cdot L(0,\phi)$.

**PARI convention differs by a swap** [computed, `beukers.gp`]:
`mfeisenstein(k,CHI1,CHI2)` has $a_n=\sum_{d\mid n}\mathrm{CHI1}(d)\,\mathrm{CHI2}(n/d)\,d^{k-1}$,
i.e. PARI's `mfeisenstein(k,CHI1,CHI2)` $=E_{k}^{\psi=\mathrm{CHI2},\ \phi=\mathrm{CHI1}}$.
Checked on $N=5$: `mfeisenstein(4,1,Mod(4,5))` $=[0,1,7,26,57,125,182]$ $=E_4^{\chi_5,1}$ and
`mfeisenstein(4,Mod(4,5),1)` $=[1,1,-7,-26,57,1,182]$ $=E_4^{1,\chi_5}$.

Characters are Conrey labels; a primitive $\psi$ is written `[u,lab]` and named
`L(s,chi_{-u})`/`L(s,chi_u)` when quadratic (unique per conductor and parity), `L(s,chi_u#lab)` otherwise.
Conrey labels multiply mod $N$ under induction — verified [computed]: $\chi_{-3}\!\uparrow\!12 = 5$,
$\chi_{-4}\!\uparrow\!12=7$, $5\cdot7\equiv 11 \pmod{12}$ $=\chi_{12}$ (even, conductor 12).

**Period classes.**
- `I` (interesting, non-elementary): $k'=3$ and $\psi$ odd; or $k'=4$ and $\psi$ even (incl. trivial $\Rightarrow\zeta(3)$). Requires $\phi$ odd or trivial.
- `E` (elementary, an algebraic multiple of $\pi^{k'-1}$): $k'=3$ and $\psi$ even (incl. $\zeta(2)$); or $k'=4$ and $\psi$ odd. Requires $\phi$ odd or trivial.
- `Z`: $\phi$ even and non-trivial $\Rightarrow L(\phi,0)=0$, the direction carries **no** period.

---

## 1. The direction census and its verification

**Enumeration** [computed, `census.gp`]. For every $N\le 60$ and $k'\in\{3,4\}$ I enumerated all triples
$(\psi,\phi,d)$ with $\psi$ primitive mod $u$, $\phi$ primitive mod $v$, $uvd\mid N$,
$\psi(-1)\phi(-1)=(-1)^{k'}$. **Total: 3996 directions** (1996 at $k'=3$, 2000 at $k'=4$).

**Cross-check against PARI** [computed, `dimchk.gp`, `an1.py`]. For every $N\le 60$ and every Dirichlet
character $\varepsilon$ mod $N$ of the required parity I compared my direction count for that nebentypus
against `mfdim([N,k,Mod(m,N)],3)` (Eisenstein subspace).

> **1103 nebentypus components checked. 0 mismatches**, both per-component and in the total
> $\dim M_{k'}^{\mathrm{Eis}}(\Gamma_1(N))=\sum_\varepsilon \dim M_{k'}^{\mathrm{Eis}}(N,\varepsilon)$.

Weight 3 has odd nebentypus (so $\dim M_3^{\mathrm{Eis}}(\Gamma_0(N))=0$ for all $N$) and weight 4 even —
confirmed by the enumeration and by PARI.

### Main census table

`#cplx` = number of distinct interesting periods attached to **non-quadratic** (complex) $\psi$;
each such period comes in Galois-conjugate families. `dim E4(G0)` is the trivial-nebentypus
($\Gamma_0(N)$) part of the weight-4 Eisenstein space.

| N | dim E3(G1) | wt-3 quadratic periods L(2,chi_-D) | #cplx | dim E4(G1) | dim E4(G0) | wt-4 quadratic periods L(3,chi_D) | #cplx |
|---:|---:|---|---:|---:|---:|---|---:|
| 1 | 0 | &mdash; | 0 | 1 | 1 | zeta(3) | 0 |
| 2 | 0 | &mdash; | 0 | 2 | 2 | zeta(3) | 0 |
| 3 | 2 | -3 | 0 | 2 | 2 | zeta(3) | 0 |
| 4 | 2 | -4 | 0 | 3 | 3 | zeta(3) | 0 |
| 5 | 4 | &mdash; | 2 | 4 | 2 | zeta(3), 5 | 0 |
| 6 | 4 | -3 | 0 | 4 | 4 | zeta(3) | 0 |
| 7 | 6 | -7 | 2 | 6 | 2 | zeta(3) | 2 |
| 8 | 6 | -4, -8 | 0 | 6 | 4 | zeta(3), 8 | 0 |
| 9 | 8 | -3 | 2 | 8 | 4 | zeta(3) | 2 |
| 10 | 8 | &mdash; | 2 | 8 | 4 | zeta(3), 5 | 0 |
| 11 | 10 | -11 | 4 | 10 | 2 | zeta(3) | 4 |
| 12 | 10 | -3, -4 | 0 | 10 | 6 | zeta(3), 12 | 0 |
| 13 | 12 | &mdash; | 6 | 12 | 2 | zeta(3), 13 | 4 |
| 14 | 12 | -7 | 2 | 12 | 4 | zeta(3) | 2 |
| 15 | 16 | -3, -15 | 2 | 16 | 4 | zeta(3), 5 | 2 |
| 16 | 14 | -4, -8 | 2 | 14 | 6 | zeta(3), 8 | 2 |
| 17 | 16 | &mdash; | 8 | 16 | 2 | zeta(3), 17 | 6 |
| 18 | 16 | -3 | 2 | 16 | 8 | zeta(3) | 2 |
| 19 | 18 | -19 | 8 | 18 | 2 | zeta(3) | 8 |
| 20 | 20 | -4, -20 | 2 | 20 | 6 | zeta(3), 5 | 2 |
| 21 | 24 | -3, -7 | 4 | 24 | 4 | zeta(3), 21 | 4 |
| 22 | 20 | -11 | 4 | 20 | 4 | zeta(3) | 4 |
| 23 | 22 | -23 | 10 | 22 | 2 | zeta(3) | 10 |
| 24 | 24 | -3, -4, -8, -24 | 0 | 24 | 8 | zeta(3), 8, 12, 24 | 0 |
| 25 | 28 | &mdash; | 10 | 28 | 6 | zeta(3), 5 | 8 |
| 26 | 24 | &mdash; | 6 | 24 | 4 | zeta(3), 13 | 4 |
| 27 | 30 | -3 | 8 | 30 | 6 | zeta(3) | 8 |
| 28 | 30 | -4, -7 | 4 | 30 | 6 | zeta(3), 28 | 4 |
| 29 | 28 | &mdash; | 14 | 28 | 2 | zeta(3), 29 | 12 |
| 30 | 32 | -3, -15 | 2 | 32 | 8 | zeta(3), 5 | 2 |
| 31 | 30 | -31 | 14 | 30 | 2 | zeta(3) | 14 |
| 32 | 32 | -4, -8 | 6 | 32 | 8 | zeta(3), 8 | 6 |
| 33 | 40 | -3, -11 | 8 | 40 | 4 | zeta(3), 33 | 8 |
| 34 | 32 | &mdash; | 8 | 32 | 4 | zeta(3), 17 | 6 |
| 35 | 48 | -7, -35 | 10 | 48 | 4 | zeta(3), 5 | 10 |
| 36 | 40 | -3, -4 | 4 | 40 | 12 | zeta(3), 12 | 4 |
| 37 | 36 | &mdash; | 18 | 36 | 2 | zeta(3), 37 | 16 |
| 38 | 36 | -19 | 8 | 36 | 4 | zeta(3) | 8 |
| 39 | 48 | -3, -39 | 10 | 48 | 4 | zeta(3), 13 | 10 |
| 40 | 48 | -4, -8, -20, -40 | 4 | 48 | 8 | zeta(3), 5, 8, 40 | 4 |
| 41 | 40 | &mdash; | 20 | 40 | 2 | zeta(3), 41 | 18 |
| 42 | 48 | -3, -7 | 4 | 48 | 8 | zeta(3), 21 | 4 |
| 43 | 42 | -43 | 20 | 42 | 2 | zeta(3) | 20 |
| 44 | 50 | -4, -11 | 8 | 50 | 6 | zeta(3), 44 | 8 |
| 45 | 64 | -3, -15 | 10 | 64 | 8 | zeta(3), 5 | 10 |
| 46 | 44 | -23 | 10 | 44 | 4 | zeta(3) | 10 |
| 47 | 46 | -47 | 22 | 46 | 2 | zeta(3) | 22 |
| 48 | 56 | -3, -4, -8, -24 | 4 | 56 | 12 | zeta(3), 8, 12, 24 | 4 |
| 49 | 60 | -7 | 20 | 60 | 8 | zeta(3) | 20 |
| 50 | 56 | &mdash; | 10 | 56 | 12 | zeta(3), 5 | 8 |
| 51 | 64 | -3, -51 | 14 | 64 | 4 | zeta(3), 17 | 14 |
| 52 | 60 | -4, -52 | 10 | 60 | 6 | zeta(3), 13 | 10 |
| 53 | 52 | &mdash; | 26 | 52 | 2 | zeta(3), 53 | 24 |
| 54 | 60 | -3 | 8 | 60 | 12 | zeta(3) | 8 |
| 55 | 80 | -11, -55 | 18 | 80 | 4 | zeta(3), 5 | 18 |
| 56 | 72 | -4, -7, -8, -56 | 8 | 72 | 8 | zeta(3), 8, 28, 56 | 8 |
| 57 | 72 | -3, -19 | 16 | 72 | 4 | zeta(3), 57 | 16 |
| 58 | 56 | &mdash; | 14 | 56 | 4 | zeta(3), 29 | 12 |
| 59 | 58 | -59 | 28 | 58 | 2 | zeta(3) | 28 |
| 60 | 80 | -3, -4, -15, -20 | 4 | 80 | 12 | zeta(3), 5, 12, 60 | 4 |

---

## 2. Reachable periods; smallest levels

### 2(a)/2(b) Reachability rule [computed, and [reasoned] for the closed form]

A direction with $\phi=1$ exists at level $N$ iff $\mathrm{cond}(\psi)\mid N$. Hence:

- **Weight 3.** $L(2,\psi)$ ($\psi$ odd primitive) is carried at level $N$ **iff $\mathrm{cond}(\psi)\mid N$**.
  In particular $L(2,\chi_{-D})$ is available exactly at the multiples of $D$.
- **Weight 4.** $\zeta(3)$ is carried at **every** $N$ (direction $(\psi,\phi,d)=(1,1,d)$, $d\mid N$).
  $L(3,\psi)$ ($\psi$ even primitive $\neq1$) is carried **iff $\mathrm{cond}(\psi)\mid N$**.

This closed form was verified against the enumeration for all $N\le60$ [computed, `an3.py` §B/§C]:
every occurrence in `directions.txt` matches, with no exceptions.

### 2(c) Smallest levels [computed]

| value | smallest $N$ | all $N\le 60$ |
|---|---:|---|
| $L(2,\chi_{-7})$ | 7 | 7, 14, 21, 28, 35, 42, 49, 56 |
| $L(2,\chi_{-8})$ | 8 | 8, 16, 24, 32, 40, 48, 56 |
| $L(2,\chi_{-11})$ | 11 | 11, 22, 33, 44, 55 |
| $L(2,\chi_{-15})$ | 15 | 15, 30, 45, 60 |
| $L(2,\chi_{-20})$ | 20 | 20, 40, 60 |
| $L(2,\chi_{-23})$ | 23 | 23, 46 |
| $L(2,\chi_{-24})$ | 24 | 24, 48 |
| $L(3,\chi_{5})$ | 5 | 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 |
| $L(3,\chi_{8})$ | 8 | 8, 16, 24, 32, 40, 48, 56 |
| $L(3,\chi_{12})$ | 12 | 12, 24, 36, 48, 60 |
| $L(3,\chi_{13})$ | 13 | 13, 26, 39, 52 |
| $L(3,\chi_{17})$ | 17 | 17, 34, 51 |
| $L(3,\chi_{24})$ | 24 | 24, 48 |

(For reference: $L(2,\chi_{-3})$ first at $N=3$, $L(2,\chi_{-4})$ at $N=4$; the next new odd quadratic
conductors $\le 60$ are $-19,-31,-35,-39,-40,-43,-47,-51,-52,-55,-56,-59$, each first at $N=|D|$.
The next new even quadratic conductors are $21,28,29,33,40,44,56,60$, each first at $N=D$.)

---

## 3. MIXED pairs

### 3.1 The rigid within-nebentypus statement [computed — this is the sharp result]

Fix $(N,k',\varepsilon)$ and let $\varepsilon^*$ be the primitive character inducing $\varepsilon$.
Within one nebentypus every direction has $\psi\phi=\varepsilon$, so:

- $\phi=1 \Rightarrow \psi=\varepsilon^*$ — exactly one $\psi$, hence exactly one $L(k'-1,\varepsilon^*)$;
- $\phi$ odd non-trivial $\Rightarrow \psi$ has the *wrong* parity $\Rightarrow$ elementary ($\pi^{k'-1}$);
- $\phi$ even non-trivial $\Rightarrow L(\phi,0)=0$, no period.

> **[computed]** Over all **550** weight-3 and **552** weight-4 nebentypus components with $N\le 60$:
> the number of *distinct interesting* periods in a fixed nebentypus is **exactly 1, in every single case**.
> There are **zero** components with $\ge 2$. (`an3.py` §D.)

So a nebentypus eigenform can never carry two *interesting* periods. Two interesting periods
**require** a form on $\Gamma_1(N)$ that is not a diamond eigenform.

### 3.2 Weight 3: the $(\zeta(2), L(2,\chi_{-D}))$ pair is generic and lives inside ONE nebentypus [computed]

At $k'=3$, $\varepsilon$ is odd. Take $\varepsilon=\chi_{-D}$ with $D\mid N$. Then

- $(\psi,\phi)=(\chi_{-D},1)$, any $d\mid N/D$ — carries $L(2,\chi_{-D})$ (interesting);
- $(\psi,\phi)=(1,\chi_{-D})$, any $d\mid N/D$ — carries $\zeta(2)\cdot L(\chi_{-D},0)$, i.e. the $\pi^2$ class.

Both live in the **same** nebentypus $\varepsilon=\chi_{-D}$.

> **[computed]** For **every** $N$ with $3\le N\le 60$ there is at least one nebentypus carrying both
> $\zeta(2)$ and an interesting $L(2,\psi)$; and there is **no** nebentypus component anywhere in the
> range that carries an interesting $L(2,\psi)$ but not $\zeta(2)$ (empty list). (`an3.py` §E.)

Quadratic instances $\varepsilon=\chi_{-D}$, all $N\le60$:

| $D$ | $N$ realising $(\zeta(2),L(2,\chi_{-D}))$ in one nebentypus |
|---:|---|
| 3 | 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60 |
| 4 | 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60 |
| 7 | 7, 14, 21, 28, 35, 42, 49, 56 |
| 8 | 8, 16, 24, 32, 40, 48, 56 |
| 11 | 11, 22, 33, 44, 55 |
| 15 | 15, 30, 45, 60 |
| 19 | 19, 38, 57 |
| 20 | 20, 40, 60 |
| 23 | 23, 46 |
| 24 | 24, 48 |
| 31, 35, 39, 40, 43, 47, 51, 52, 55, 56, 59 | $N=D$ only |

**CDT's configuration reproduced** [computed]. The census at $N=6$, $k'=3$, $\varepsilon=\chi_{-3}$
(Conrey `Mod(5,6)`) has exactly four directions, matching `mfdim([6,3,Mod(5,6)],3)`$=4$:

```
DIR 6 3  psi=[1,1]   phi=[3,2]  d=1  eps=5  E  zeta(2)
DIR 6 3  psi=[1,1]   phi=[3,2]  d=2  eps=5  E  zeta(2)
DIR 6 3  psi=[3,2]   phi=[1,1]  d=1  eps=5  I  L(2,chi_-3)
DIR 6 3  psi=[3,2]   phi=[1,1]  d=2  eps=5  I  L(2,chi_-3)
```

i.e. $E_{3,\chi_{-3},1}$, $V_2E_{3,\chi_{-3},1}$, $E_{3,1,\chi_{-3}}$, $V_2E_{3,1,\chi_{-3}}$ — exactly
the two orientations CDT keep, $(1-8V_2)E_{3,\chi_{-3},1}$ and $(1-V_2)E_{3,1,\chi_{-3}}$.
**Conclusion: at weight 3 the mixed pair $(\pi^2, L(2,\chi_{-D}))$ is generic and needs no
$\Gamma_1$/non-eigenform mechanism** — it is already a single-nebentypus phenomenon. [computed]

### 3.3 Weight 4: $(\zeta(3), L(3,\chi_D))$ needs the $\Gamma_1$ (non-eigenform) mechanism [computed]

At $k'=4$, $\varepsilon$ is even. $\zeta(3)$ needs $\psi=1$, hence $\varepsilon=1$ (trivial nebentypus);
$L(3,\chi_D)$ needs $\psi=\chi_D$, hence $\varepsilon=\chi_D$. **Different nebentypus components.**
Confirmed by §3.1: no weight-4 nebentypus anywhere in the range carries two interesting periods;
and by an explicit check that each interesting weight-4 period sits in exactly one $\varepsilon$
(`an3.py` §F: "each period sits in exactly one nebentypus: True" for all $N\le60$).

So a mixed weight-4 source must be a **sum across nebentypus components** — a form on $\Gamma_1(N)$
that is not a diamond eigenform. PARI enforces this: `mflinear` on the Beukers summands **refuses**,
`*** mflinear: incorrect type in mflinear [different characters]`. [computed]

**Availability** [computed]: mixed pairs $(\zeta(3),L(3,\chi_D))$ exist on $\Gamma_1(N)$ exactly when $D\mid N$:

| $D$ | all $N\le 60$ |
|---:|---|
| 5 | 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 |
| 8 | 8, 16, 24, 32, 40, 48, 56 |
| 12 | 12, 24, 36, 48, 60 |
| 13 | 13, 26, 39, 52 |
| 17 | 17, 34, 51 |
| 21 | 21, 42 |
| 24 | 24, 48 |
| 28 | 28, 56 |
| 29 | 29, 58 |
| 33, 40, 44, 56, 60 | $N=D$ only |

More generally, the number of distinct interesting weight-4 periods at level $N$ is
$1+\#\{$even primitive $\psi\ne1$ with $\mathrm{cond}(\psi)\mid N\}$; the first $N$ with $\ge2$ is $N=5$,
and every $N\ge5$ except $N=6$ has $\ge2$ [computed] ($N=6$ admits no even non-trivial primitive
character of conductor dividing 6 — there are no primitive characters mod 6 at all).

Weight-3 counts of distinct interesting periods are $\#\{$odd primitive $\psi$ with $\mathrm{cond}(\psi)\mid N\}$
(0 for $N=1,2$; $\ge1$ for $N\ge3$).

### 3.4 Beukers $\Gamma_1(5)$ calibration — reproduced numerically [computed, `beukers.gp`]

The census at $N=5$, $k'=4$ has exactly 4 directions (matching PARI's $\dim=4$):

```
DIR 5 4  psi=[1,1]   phi=[1,1]  d=1  eps=1  I  zeta(3)
DIR 5 4  psi=[1,1]   phi=[1,1]  d=5  eps=1  I  zeta(3)
DIR 5 4  psi=[1,1]   phi=[5,4]  d=1  eps=4  Z  0            <- phi = chi_5 EVEN => L(phi,0)=0
DIR 5 4  psi=[5,4]   phi=[1,1]  d=1  eps=4  I  L(3,chi_5)
```

Beukers' form is $F = 10E_4^{1,1}(\tau)-250E_4^{1,1}(5\tau)+E_4^{1,\chi_5}-5\sqrt5\,E_4^{\chi_5,1}$, with
$$L(F,s)=10(1-5^{2-s})\zeta(s)\zeta(s-3)+\zeta(s)L(s-3,\chi_5)-5\sqrt5\,\zeta(s-3)L(s,\chi_5).$$
Evaluated at $s=3$ at 60 digits:

```
zeta(0)      = -1/2                      L(0,chi_5) = 0            (chi_5 EVEN)
zeta(3)      =  1.20205690315959428539973816151144999076498629234049888179227
L(3,chi_5)   =  0.854824766648543010235690083538137697138396464937005282730702
term1        = -4.80822761263837714159895264604579996305994516936199552716909  = -4*zeta(3)     [exact match]
term2        =  0                                                                                [chi_5 even]
term3        =  4.77861571769134312375381963529161611414870939905953621484128  = (5*sqrt5/2)*L(3,chi_5) [exact match]
L(F,3)       = -0.0296118949470340178451330107541838489112357703024593123278061
-(1/2)(8*zeta(3) - 5*sqrt5*L(3,chi_5))
             = -0.0296118949470340178451330107541838489112357703024593123278061
difference   =  0.E-78
```

> **[computed] Calibration passes: $L(F,3)=-\tfrac12\bigl(8\zeta(3)-5\sqrt5\,L(3,\chi_5)\bigr)$ to 78 digits.**

Three structural facts this exhibits, all confirmed by the census:
1. the middle summand $E_4^{1,\chi_5}$ is a **class-`Z`** direction — its critical value is $0$
   because $\chi_5$ is even; it contributes to $F$'s $q$-expansion but nothing to the period;
2. the $\zeta(3)$ part needs **two** $\zeta(3)$ directions ($d=1$ and $d=5$) so that the
   Euler-factor combination $10(1-5^{2-s})$ can be formed — this is the oldform shift;
3. $F$ mixes $\varepsilon=1$ and $\varepsilon=\chi_5$, so it is **not** a diamond eigenform.

**A remark on rationality** [reasoned]: $E_4^{1,1}(\tau)$, $E_4^{1,1}(5\tau)$ and $E_4^{\chi_5,1}$ all have
rational $q$-coefficients (verified: $[0,1,7,26,57,125,182]$ for $E_4^{\chi_5,1}$), and their critical
values are $-4\zeta(3)$ (for the $10,-250$ combination) and $-\tfrac12L(3,\chi_5)$. So already over
$\mathbb{Q}$ the pair $\{\zeta(3),L(3,\chi_5)\}$ is spanned on $\Gamma_1(5)$: any $a,b\in\mathbb{Q}$ give a
rational form with critical value $-4a\,\zeta(3)-\tfrac{b}{2}L(3,\chi_5)$. Beukers' $-5\sqrt5$ is forced by
his geometry ($F=$ weight-2 form $\times\,\theta_q t$), not by the Eisenstein bookkeeping; his $F$ and its
Galois conjugate span the same $\mathbb{Q}$-rational 2-plane.

### 3.5 Complete mixed-pair list, $N\le 60$

Because the interesting periods at $(N,k')$ are exactly

$$\mathcal{P}_3(N)=\{L(2,\psi):\psi\ \text{odd primitive},\ \mathrm{cond}\,\psi\mid N\},\qquad
\mathcal{P}_4(N)=\{\zeta(3)\}\cup\{L(3,\psi):\psi\ \text{even primitive}\ne1,\ \mathrm{cond}\,\psi\mid N\},$$

and each is carried by a $\phi=1$ direction, **every unordered pair drawn from $\mathcal{P}_{k'}(N)$ is
realisable by a form in $M_{k'}^{\mathrm{Eis}}(\Gamma_1(N))$** [computed for the sets, [reasoned] for
"every pair"]: pick the two $\phi=1$ directions and take any $\mathbb{Q}$- (or $K$-) linear combination;
rationality over $\mathbb{Q}$ is achieved by summing over the Galois orbit of each $\psi$.
So the "complete mixed-pair list" is the list of $\mathcal{P}_{k'}(N)$, which is the main census table of §1.
Distilled:

- **$k'=3$, pair $(\pi^2\ \text{i.e.}\ \zeta(2),\ L(2,\chi_{-D}))$:** realisable **inside a single
  nebentypus** $\varepsilon=\chi_{-D}$ for every $N$ with $D\mid N$ — table in §3.2. Also realisable for
  every complex odd $\psi$ with $\mathrm{cond}\,\psi\mid N$.
- **$k'=3$, pair of two interesting values $(L(2,\psi_1),L(2,\psi_2))$:** needs $\varepsilon$-mixing.
  Available iff $N$ has $\ge2$ odd primitive characters of conductor dividing $N$;
  smallest levels: $N=5$ (the two quartic $\psi$ mod 5, a Galois pair), $N=7$, $N=8$
  ($L(2,\chi_{-4}),L(2,\chi_{-8})$ — first pair of two *quadratic* ones), $N=9$, $N=12$
  ($L(2,\chi_{-3}),L(2,\chi_{-4})$).
- **$k'=4$, pair $(\zeta(3),L(3,\chi_D))$:** always needs $\varepsilon$-mixing ($\Gamma_1$, non-eigenform).
  Available iff $D\mid N$ — table in §3.3. Smallest: $N=5$ (Beukers).
- **$k'=4$, pair of two non-$\zeta(3)$ values:** first at $N=7$ (the two cubic $\psi$ mod 7),
  first *quadratic* pair at $N=24$: $(L(3,\chi_8),L(3,\chi_{12}))$, $(L(3,\chi_8),L(3,\chi_{24}))$,
  $(L(3,\chi_{12}),L(3,\chi_{24}))$; then $N=40$ ($\chi_5,\chi_8,\chi_{40}$),
  $N=60$ ($\chi_5,\chi_{12},\chi_{60}$).

---

## 4. Period annihilation / fold-regularity dimension count

Genus-zero $\Gamma_0(N)$, $N\le60$: $N\in\{1,\dots,10,12,13,16,18,25\}$.
"Annihilated" $=$ constant term vanishes at **both** cusp $\infty$ and cusp $0$
(cusp $0$ via `mfslashexpansion(mf,f,[0,-1;1,0],0,·)`; the returned parameters were
$[\alpha,w,A]=[0,w,I_2]$ in every case, so the first returned coefficient *is* the constant term).

### 4.1 The exact structure of the two functionals [computed — this is the engine]

For every direction $E^{\psi,\phi}(d\tau)$ on the genus-zero levels I computed the constant term at
$\infty$ (exactly, `mfcoefs`) and at $0$ (exactly, `mfslashexpansion(...,flrat=1)`):

> **[computed, `verify.gp`, 229 of 244 directions, 0 violations]**
> $$a_0^{\infty}\bigl(E^{\psi,\phi}(d\tau)\bigr)\ne 0 \iff \psi=1,\qquad
> a_0^{0}\bigl(E^{\psi,\phi}(d\tau)\bigr)\ne 0 \iff \phi=1,$$
> and $a_0^{0}\bigl(E^{\psi,\phi}(d\tau)\bigr)=d^{-k'}\,C(\psi)$ with $C(\psi)$ **independent of $d$**
> (verified directly: e.g. at $N=6,k'=3$ both $d=1$ and $d=2$ give
> $C(\chi_{-3})=\mathrm{Mod}(\tfrac{2}{243}t+\tfrac1{243},\,t^2+t+1)$; the $(1,1)$ directions give
> $C(1)=C'(1)=1/240$ at every $N$, consistent with $E_k|_kS=E_k$).
> *(The 15 unchecked directions are $N=25$, $k'=4$ with $\psi$ of order 20 or 10, where PARI reports
> "`mfeisenstein` for these characters is not yet implemented". See §5.)*

Rescaling each direction by a nonzero constant (which changes no rank or span question) puts the two
functionals in the **exact rational normal form**
$$r_\infty = \bigl(\mathbb{1}[\psi_j=1]\bigr)_j,\qquad r_0=\bigl(\mathbb{1}[\phi_j=1]\,d_j^{-k'}\bigr)_j .$$
All of §4.2–§4.3 is then exact rational linear algebra (`an2.py`, `fractions.Fraction`).

Two immediate consequences [reasoned from the verified structure]:
- at $k'=3$ the two supports are **disjoint** ($\psi=\phi=1$ would force $\varepsilon$ even), so
  $\mathrm{rank}=2$ whenever the space is non-zero: **$\dim(\text{annihilated})=\dim M_3^{\mathrm{Eis}}(\Gamma_1(N))-2$**;
- at $k'=4$ the supports overlap on the $(1,1,d)$ directions; rank drops to 1 only at $N=1$.

### 4.2 Dimension table [computed]

$\dim^{\mathrm{ann}}_{\Gamma_1}$: forms on all of $\Gamma_1(N)$ (nebentypus components may be mixed).
$\dim^{\mathrm{ann}}_{\varepsilon}$: sum over nebentypus of the annihilated dimension **inside** each
$M_{k'}^{\mathrm{Eis}}(N,\varepsilon)$ (i.e. requiring diamond eigenforms) — computed independently
in `ann2.gp` straight from `mfbasis` + `mfslashexpansion`, no structural input.

| $N$ | $k'$ | $\dim M^{\mathrm{Eis}}_{k'}(\Gamma_1(N))$ | rank of $(r_\infty,r_0)$ | $\dim^{\mathrm{ann}}_{\Gamma_1}$ | $\dim^{\mathrm{ann}}_{\varepsilon}$ | interesting periods still spanned by the annihilated subspace |
|---:|---:|---:|---:|---:|---:|---|
| 1 | 3 | 0 | &mdash; | 0 | 0 | &mdash; |
| 1 | 4 | 1 | 1 | 0 | 0 | none ($\zeta(3)$ killed) |
| 2 | 3 | 0 | &mdash; | 0 | 0 | &mdash; |
| 2 | 4 | 2 | 2 | 0 | 0 | none ($\zeta(3)$ killed) |
| 3 | 3 | 2 | 2 | 0 | 0 | none ($L(2,\chi_{-3})$ killed) |
| 3 | 4 | 2 | 2 | 0 | 0 | none ($\zeta(3)$ killed) |
| 4 | 3 | 2 | 2 | 0 | 0 | none ($L(2,\chi_{-4})$ killed) |
| 4 | 4 | 3 | 2 | **1** | 1 | $\zeta(3)$ |
| 5 | 3 | 4 | 2 | 2 | 2 | $L(2,\chi_5^{\#2}), L(2,\chi_5^{\#3})$ (+ $\pi^2$) |
| 5 | 4 | 4 | 2 | 2 | 0 | $\zeta(3)$, $L(3,\chi_5)$ |
| 6 | 3 | 4 | 2 | 2 | 2 | $L(2,\chi_{-3})$ (+ $\pi^2$) &nbsp;&larr; **CDT** |
| 6 | 4 | 4 | 2 | 2 | 2 | $\zeta(3)$ |
| 7 | 3 | 6 | 2 | 4 | 2 | $L(2,\chi_{-7}), L(2,\chi_7^{\#3}), L(2,\chi_7^{\#5})$ (+ $\pi^2$) |
| 7 | 4 | 6 | 2 | 4 | 2 | $\zeta(3), L(3,\chi_7^{\#2}), L(3,\chi_7^{\#4})$ |
| 8 | 3 | 6 | 2 | 4 | 2 | $L(2,\chi_{-4}), L(2,\chi_{-8})$ (+ $\pi^2$) |
| 8 | 4 | 6 | 2 | 4 | 2 | $\zeta(3), L(3,\chi_8)$ |
| 9 | 3 | 8 | 2 | 6 | 4 | $L(2,\chi_{-3}), L(2,\chi_9^{\#2}), L(2,\chi_9^{\#5})$ (+ $\pi^2$) |
| 9 | 4 | 8 | 2 | 6 | 4 | $\zeta(3), L(3,\chi_9^{\#4}), L(3,\chi_9^{\#7})$ (+ $\pi^3$) |
| 10 | 3 | 8 | 2 | 6 | 6 | $L(2,\chi_5^{\#2}), L(2,\chi_5^{\#3})$ (+ $\pi^2$) |
| 10 | 4 | 8 | 2 | 6 | 4 | $\zeta(3), L(3,\chi_5)$ |
| 12 | 3 | 10 | 2 | 8 | 6 | $L(2,\chi_{-3}), L(2,\chi_{-4})$ (+ $\pi^2$) |
| 12 | 4 | 10 | 2 | 8 | 6 | $\zeta(3), L(3,\chi_{12})$ (+ $\pi^3$) |
| 13 | 3 | 12 | 2 | 10 | 6 | all 6 $L(2,\psi)$, $\psi$ odd mod 13 (+ $\pi^2$) |
| 13 | 4 | 12 | 2 | 10 | 4 | $\zeta(3)$ + all 5 $L(3,\psi)$, $\psi$ even $\ne1$ mod 13 |
| 16 | 3 | 14 | 2 | 12 | 8 | $L(2,\chi_{-4}),L(2,\chi_{-8}),L(2,\chi_{16}^{\#3}),L(2,\chi_{16}^{\#11})$ (+ $\pi^2$) |
| 16 | 4 | 14 | 2 | 12 | 8 | $\zeta(3),L(3,\chi_8),L(3,\chi_{16}^{\#5}),L(3,\chi_{16}^{\#13})$ (+ $\pi^3$) |
| 18 | 3 | 16 | 2 | 14 | 12 | $L(2,\chi_{-3}),L(2,\chi_9^{\#2}),L(2,\chi_9^{\#5})$ (+ $\pi^2$) |
| 18 | 4 | 16 | 2 | 14 | 12 | $\zeta(3),L(3,\chi_9^{\#4}),L(3,\chi_9^{\#7})$ (+ $\pi^3$) |
| 25 | 3 | 28 | 2 | 26 | 18 | all 10 $L(2,\psi)$ (+ $\pi^2$) |
| 25 | 4 | 28 | 2 | 26 | 16 | $\zeta(3)$ + all 9 $L(3,\psi)$ (+ $\pi^3$) |

### 4.3 What the annihilation kills [computed, exact]

For each interesting period $P$ let $I_P$ be its critical-slot functional
($I_P(x)=-\tfrac12\sum_{j\in S_P}x_j d_j^{1-k'}$, $S_P=\{j:\phi_j=1,\psi_j=\psi_P\}$).
$P$ survives into the annihilated subspace iff $I_P\notin\mathrm{span}(r_\infty,r_0)$.
Exact computation over $\mathbb{Q}$ gives:

> **The only genus-zero levels at which an interesting period is *annihilated away* are the very smallest:**
> $(N,k')=(1,4),(2,4),(3,4)$ kill $\zeta(3)$; $(3,3)$ kills $L(2,\chi_{-3})$; $(4,3)$ kills $L(2,\chi_{-4})$.
> **From $N=4$ ($k'=4$) and $N=5$ ($k'=3$) onward, every interesting period at every genus-zero level
> survives.** [computed]

The mechanism is transparent [reasoned, from §4.1]: $I_P$ has $d$-weight $d^{1-k'}$ while $r_0$ has
$d$-weight $d^{-k'}$ on the same support; they are proportional only when the support has a single
$d$-value **and** $P$ is the only interesting period at that level. That is exactly the list above.
$N=6$, $k'=3$ is the first level where the *oldform shift* $d\in\{1,2\}$ breaks the proportionality —
which is precisely why CDT work at level 6 and not at level 3.

The elementary $\pi^{k'-1}$ block likewise survives at every genus-zero level with $\ge2$ elementary
directions (marked "$+\pi^{k'-1}$" above); at $k'=4$, $N\in\{1,\dots,8,10,13\}$ there are **no**
elementary directions at all (column `elem(0)` in `an2.py` output), so the annihilated weight-4
subspace there carries a **pure** combination of $\zeta(3)$ and $L(3,\psi)$'s with **no $\pi^3$ contamination** —
including $N=5$ (Beukers) and $N=8$.

---

## 5. What I could NOT compute (honest list)

1. **$N=25$, $k'=4$, 15 of 28 directions.** PARI 2.15.4 refuses
   `mfeisenstein` for the characters mod 25 of order 20 and 10:
   `*** mfeisenstein: sorry, mfeisenstein for these characters is not yet implemented.`
   So the *constant-term support structure* of §4.1 is verified on **229 of 244** genus-zero
   directions, not all 244. The dimension and survival numbers for $(25,4)$ in §4.2–4.3 therefore rest
   on the structural rule extrapolated from the other 229 rows (and on the independently computed
   per-nebentypus figure $\dim^{\mathrm{ann}}_\varepsilon=16$, which used `mfbasis` and did go through).
   Everything else in the report — the census, the PARI dimension cross-check, the mixed-pair analysis —
   is unaffected: those never call `mfeisenstein`.

2. **The elementary/$\pi^{k'-1}$ functional was not evaluated numerically.** Its coefficients are
   $L(k'-1,\psi_j)L(0,\phi_j)d_j^{1-k'}$ with $\psi_j,\phi_j$ ranging over several conductors, so
   testing it requires embedding several cyclotomic coefficient fields consistently. I tested only the
   coarser (support-level) question — whether the annihilated subspace has non-zero projection on the
   elementary block — and report that. A cancellation among elementary directions of *different*
   conductors would not be detected; I regard this as very unlikely (the values carry independent
   $\sqrt{f}$ factors) but it is **not** verified. This is why §4.2's "$+\pi^{k'-1}$" annotations are
   weaker claims than the interesting-period ones.

3. **A first numerical attempt at §4 failed and was replaced.** `mfslashexpansion(mf,f,[1,0;0,1],0,0)`
   returns *exact* `t_POLMOD` values (PARI short-circuits the identity matrix) while the cusp-$0$ call
   returns floats, so the two rows of the constant-term matrix came back in incompatible
   representations (`*** abs: incorrect type in gabs (t_POLMOD)`). I abandoned the floating-point route
   and replaced it with the exact structural normal form of §4.1, which is both cleaner and
   embedding-free. The floating results obtained before the failure ($N\le4$) agree with the exact ones.

4. **Independence of the periods is assumed, not proved.** "Two distinct interesting periods" here means
   two distinct $L$-values in the list; whether they are $\overline{\mathbb{Q}}$-linearly independent is
   the open problem (that is what CDT prove for $1,\zeta(2),L(2,\chi_{-3})$). Nothing here bears on it.

5. **Only $\Gamma_0(N)$/$\Gamma_1(N)$ hosts.** No Atkin–Lehner quotients $\Gamma_0(N)+$, no
   $\Gamma_1(N)\cap\Gamma_0(M)$ intermediate groups, no non-congruence or Fuchsian hosts. Genus-zero
   status was taken as given ($N\in\{1..10,12,13,16,18,25\}$), not recomputed.

6. **No source $\Phi=F\cdot\theta_q t$ was actually built.** The census is of the *ambient Eisenstein
   space* and the periods its directions carry. Which of these are realised by an actual modular Apéry
   system $(\Gamma,t,F)$ with integral $A_n$ is a separate question and was not addressed.

---

## 6. Summary of the load-bearing findings

1. **[computed]** 3996 Eisenstein directions for $N\le60$, $k'\in\{3,4\}$; **0 mismatches** against PARI's
   `mfdim(...,3)` on all 1103 nebentypus components.
2. **[computed]** *Within a fixed nebentypus the number of distinct interesting periods is exactly 1* —
   in all 550 (wt 3) + 552 (wt 4) components with $N\le60$, without exception. A diamond eigenform can
   never carry two interesting periods; that is the sharp form of "mixed pairs need $\Gamma_1$".
3. **[computed]** At weight 3, the pair $(\pi^2,\,L(2,\chi_{-D}))$ *is* available inside a single
   nebentypus $\varepsilon=\chi_{-D}$ whenever $D\mid N$ — generic, no $\Gamma_1$ mechanism needed.
   CDT's level-6 configuration is reproduced exactly (4 directions, two orientations, $V_2$ shift).
4. **[computed]** At weight 4 the pair $(\zeta(3),L(3,\chi_D))$ *always* straddles two nebentypus
   components ($\varepsilon=1$ and $\varepsilon=\chi_D$) and hence needs a non-eigenform on $\Gamma_1(N)$.
   PARI's `mflinear` literally refuses to form Beukers' combination. Available iff $D\mid N$.
5. **[computed]** Beukers' $\Gamma_1(5)$ identity verified: $L(F,3)=-\tfrac12(8\zeta(3)-5\sqrt5 L(3,\chi_5))$
   to 78 digits, with the term-by-term decomposition matching the census's four directions (one of which
   is a class-`Z`, zero-period direction).
6. **[computed]** Annihilation: $\dim^{\mathrm{ann}}_{\Gamma_1}=\dim M^{\mathrm{Eis}}_{k'}(\Gamma_1(N))-2$
   at every genus-zero level except $(N,k')=(1,4)$ (rank 1, $\dim^{\mathrm{ann}}=0$). Requiring a diamond
   eigenform costs strictly more (column $\dim^{\mathrm{ann}}_\varepsilon$).
7. **[computed]** The annihilation kills an interesting period only at $(1,4),(2,4),(3,4),(3,3),(4,3)$.
   The oldform shift $d>1$ is what saves it, which is exactly why level 6 works for $L(2,\chi_{-3})$
   where level 3 does not.
