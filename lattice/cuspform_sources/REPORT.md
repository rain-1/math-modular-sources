# Cusp forms as Apéry sources: the Fricke fold is regular exactly on the W_N = −1 eigenspace

*Working directory `lattice/cuspform_sources/`. All scripts and raw outputs listed in §11.
Conventions: `consolidation/COMPANION_ARITHMETIC.md` §1 (companion formula, B = F·D^{−r}Φ),
`consolidation/ASYMPTOTIC_CONSTANTS.md` Thm 3.4 and §5 (weight-two structure, K±),
`lattice/hostscan/REPORT.md` §§4, 8, 9 (the twelve CDT-shape Fricke hosts, the fold/cusp
dichotomy, the sources).  Tags: **[exact]** = closed form / exact rational arithmetic;
**[verified, d digits]** = high-precision computation in this task; **[estimated]**.
No irrationality claim is made — §9 says why none is available.*

---

## 0. Verdict

| claim | verdict |
|---|---|
| **The criterion.** On a Fricke host (x = u/(1+Bu+Cu²), F = D log u, r = 3), a weight-four source Φ ∈ M₄(Γ₀(N)) with a₀(Φ) = 0 is **fold-regular at the near fold ⟺ Φ\|₄W_N = −Φ**. The period-polynomial conditions it must satisfy, L(Φ,2) = 0 and Λ(Φ,1) = −N·Λ(Φ,3), are then **automatic** from the functional equation Λ(s) = ε·Λ(4−s). | **[exact]**, proved §1 |
| The near fold is τ_c = i/√N on **every one of the twelve hosts**: u(i/√N) = 1/√C hence x = 1/λ₁, verified to 60 digits. | **[verified, 60 digits]** §1.4 |
| **No cusp form is fold-regular on either Apéry-perfect host.** dim S₄(Γ₀(6)) = dim S₄(Γ₀(7)) = 1 and both newforms (6.4.a.a, 7.4.a.a) have Fricke eigenvalue **+1**; likewise 5.4.a.a, 8.4.a.a, 9.4.a.a. So dim S₄(Γ₀(N))^{W=−1} = 0 for N = 5,6,7,8,9. | **[verified]** §2 |
| **Five fold-regular (host, cusp form) pairs do exist**, all on the non-perfect hosts N = 10, 12, 18, and every one of them is an **oldform combination** f(τ) − d²f(dτ), never a newform. Their Apéry limits are exactly L(Φ,3) — rational multiples of L(5.4.a.a,3), L(6.4.a.a,3), L(9.4.a.a,3). | **[verified, 12–144 digits]** §3 |
| The rate is **geometric**: A_nξ − B_n ~ K₋·λ₂ⁿ·n^{−3/2}, verified. For W_N = +1 the limit still exists but the convergence is only **O(1/n)** — a clean numerical dichotomy on all twelve hosts. | **[verified]** §§3, 5 |
| **The denominator power is k = 3 in every case: no free integration**, on hosts whose own canonical (meromorphic) source has one. | **[verified, n ≤ 240]** §3 |
| **No Eisenstein admixture can repair a W_N = +1 cusp form.** For ε = +1 the fold-regularity equation reads 2Θ + P = quadratic, so the obstruction is the *whole Eichler integral*, not a scalar; and M₄ = Eis ⊕ S₄ is W_N-stable, so the two +1 components lie in complementary subspaces and cannot cancel. | **[exact]** §5 |
| **Host-independence lemma.** F′(τ_c) = N τ_c F(τ_c) for *any* W_N-antiinvariant weight-two F, hence the Apéry limit of *any* weight-four source on *any* Fricke host of level N is ξ = Θ(q_c) + (2π/√N)·DΘ(q_c), q_c = e^{−2π/√N} — **independent of the host**. For ε = −1 this is a rapidly convergent series for L(Φ,3). | **[exact]** + **[verified, 96 digits]** §6 |
| **A closed form for three of the five sources**: Φ = x·F²·√(1−λ₁x) = Φ₀/√(1−λ₂x), the canonical source with the *far*-fold square root divided out. Equivalently the host's own three-term recurrence applied to the companion leaves the **exact** residual binom(2n+2,n+1)·(λ₂/4)^{n+1}. | **[exact]**, verified n ≤ 177 §7 |
| B_n satisfies **no** three-term recurrence; its minimal recurrence has order 3 (or 4) with characteristic roots {λ₁, λ₂, λ₂} (or {λ₁,λ₁,λ₂,λ₂}). The Casoratian identity K₊K₋(λ₁−λ₂) = κ of `ASYMPTOTIC_CONSTANTS.md` Thm 5.1 therefore **fails**. | **[verified]** §§7, 8 |
| **K₋ identified exactly in 2 of 5 cases**: K₋ = L(6.4.a.a,2)/√π on the Domb host and √3·L(6.4.a.a,2)/√π on Cooper's s₁₈ host — the **central** value of the newform, which is nonzero precisely because the newform's own Fricke sign is +1. Three cases not identified. | **[verified, 50 digits]** §8 |
| **No irrationality.** Only Apéry's own host has log\|λ₂\| + 3 < 0 (= −0.525), and it carries no cusp-form source; every host that does carry one has log\|λ₂\| + 3 ≥ **4.386**. The two conditions are mutually exclusive across the whole family. | **[exact]** §9 |
| **But a genuine mixed pair appears.** The fold-regular source space on the level-18 hosts is **5-dimensional** and carries four independent periods ζ(3), L(3,χ₋₃), L(6.4.a.a,3), L(9.4.a.a,3) on one host; on the level-12 hosts it is 3-dimensional with {ζ(3), L(6.4.a.a,3)}. On Apéry's own host it is **1-dimensional** — the Eisenstein line, ζ(3)/6, and nothing else. | **[verified, 144 digits]** §4 |

**One sentence.** *Cusp forms do give fold-regular Apéry sources, and the criterion is as clean
as it could be — the Fricke eigenvalue must be −1, whereupon the two period-polynomial
conditions come free from the functional equation; but the −1 eigenspace of S₄(Γ₀(N)) is empty
at N = 5,6,7,8,9, so it misses both Apéry-perfect hosts entirely and survives only at
N = 10, 12, 18, where λ₂ is 4, 4, 32, 12 and the d_n³ denominators lose by 4.4 nats — so what
one gets is five genuine Apéry-like linear forms in 1 and L(f,3) for weight-four newforms,
with exact three-term data, geometric decay and (twice) a closed-form constant L(f,2)/√π,
and no irrationality anywhere near.*

---

## 1. The criterion  **[exact]**

### 1.1 Setting

Host: u an eta quotient with W_N-antisymmetric exponent vector (u\|W_N = 1/(Cu)),
x = u/(1+Bu+Cu²), F = D log u ∈ M₂(Γ₀(N)), A_n = [xⁿ]F, λ_{1,2} = B ± 2√C.
Source: Φ ∈ M₄(Γ₀(N)) with a₀(Φ) = 0, Θ := D^{−3}Φ = Σ c(m)m^{−3}qᵐ,
B_n = [xⁿ](F·Θ), companion normalisation B₀ = 0, B₁ = c(1).

Two elementary transformation facts (both immediate from u\|W_N = 1/(Cu)):

* **x is W_N-invariant** — so x descends to X₀(N)+ and has a critical point at each W_N-fixed point (this is the *fold*);
* **F\|₂W_N = −F** and **Dx\|₂W_N = +Dx**, hence the canonical source Φ₀ = F·Dx = x F²√P(x) (Thm 3.4 of `ASYMPTOTIC_CONSTANTS.md`, P(x) = (1−λ₁x)(1−λ₂x)) satisfies **Φ₀\|₄W_N = −Φ₀** — the Eisenstein sources of the ledger are already in the −1 eigenspace.

### 1.2 The fold computation

Write τ = τ_c + s at the near fold τ_c = i/√N; W_N acts as s ↦ −s + O(s²). A = F is
[analytic in x] + (x−t₁)^{1/2}·[analytic], the singular part being the **odd** part of F in s.
So B − ξA = F·(Θ−ξ) is analytic at t₁ **iff F·(Θ−ξ) is even in s, i.e. W_N-invariant as a
plain function**. Substituting F(στ) = −Nτ²F(τ) and Θ(στ) = (Θ\|_{−2}W_N)(τ)/(Nτ²) this reads

> **(★)  Θ\|_{−2}W_N + Θ = ξ·(Nτ² + 1).**

### 1.3 Solving (★)

For the Eichler integral Θ = (2πi)³/2! · ∫_{i∞}^{τ}(τ−v)²Φ(v)dv one has, for Φ\|₄W_N = εΦ,

  Θ\|_{−2}W_N = ε·(Θ + 4π³i·R),  R(τ) = ∫₀^{i∞}(τ−v)²Φ(v)dv = iΛ(1)τ² + 2Λ(2)τ − iΛ(3),

with Λ(j) := ∫₀^∞ y^{j−1}Φ(iy)dy = (j−1)!·L(Φ,j)/(2π)^j. Hence:

* **ε = +1** makes (★) read 2Θ + P = quadratic — impossible unless Φ = 0. *No W_N-invariant source is ever fold-regular.*
* **ε = −1** makes (★) read P(τ) = ξ(Nτ²+1) with P = 4π³i·R, i.e. two conditions
  **Λ(2) = 0** and **Λ(1) = −N·Λ(3)**, and then ξ = 4π³Λ(3) = **L(Φ,3)**.
* Both conditions are **automatic**: the functional equation derived from Φ(στ) = εN²τ⁴Φ(τ) is
  Λ(s) = ε·N^{2−s}·Λ(4−s), which at s = 2 gives Λ(2) = −Λ(2) = 0 and at s = 1 gives
  Λ(1) = −N·Λ(3).

> **Theorem A.** Φ ∈ M₄(Γ₀(N)), a₀(Φ) = 0. The source Φ is fold-regular on every Fricke host
> of level N **iff Φ\|₄W_N = −Φ**, and then the Apéry limit is exactly ξ = L(Φ,3).

(For Eisenstein Φ the integrals need the usual regularisation, but ε = −1 together with
a₀(Φ,∞) = 0 forces a₀(Φ,0) = 0 as well, so both ends converge. Calibration: the Apéry source
Φ_γ = E₄ − 28E₄(2τ) + 63E₄(3τ) − 36E₄(6τ) has Φ_γ\|W₆ = −Φ_γ and
L(Φ_γ,3) = (1 − 28/8 + 63/27 − 36/216)·ζ(3)ζ(0) = ζ(3)/6 ✓, reproduced by the pipeline in
`02_calib.gp` to 60 digits.)

### 1.4 The near fold really is i/√N  **[verified, 60 digits]** (`19_ufold.gp`)

u(i/√N) = 1/√C to 60 digits on all twelve hosts, hence x(i/√N) = 1/(B+2√C) = 1/λ₁, the
**near** singularity. (Evaluating the truncated *x*-series directly at q_c diverges for the
two hosts N = 8, C = 16 and N = 12, C = 1 — q_c lies outside its disc of convergence — which
is why `10_plus.gp` prints nonsense there; going through u is correct and the fold is where
Theorem A puts it. `18_xfold.gp` records the artefact.)

---

## 2. The census  **[verified]** (`01_forms.gp`, `05_space.gp`)

W_N computed on `mfbasis(mfinit([N,4],1))` through `mfatkininit`/`mfatkin`, eigen-split exactly.

| N | dim S₄(Γ₀(N)) | W_N spectrum | **dim S₄^{W=−1}** | a −1 eigenvector |
|---|---|---|---|---|
| 5 | 1 | +1 | **0** | — |
| 6 | 1 | +1 | **0** | — |
| 7 | 1 | +1 | **0** | — |
| 8 | 1 | +1 | **0** | — |
| 9 | 1 | +1 | **0** | — |
| 10 | 3 | +1, +1, −1 | **1** | f₅(τ) − 4f₅(2τ) |
| 12 | 3 | +1, +1, −1 | **1** | f₆(τ) − 4f₆(2τ) |
| 18 | 5 | +1×3, −1×2 | **2** | f₆(τ) − 9f₆(3τ), f₉(τ) − 4f₉(2τ) |

f₅ = 5.4.a.a, f₆ = 6.4.a.a = (η₁η₂η₃η₆)², f₉ = 9.4.a.a = η(3τ)⁸. The five small levels are
exactly the four-point Γ₀(N) of `hostscan` §2 — **the entire Apéry-perfect pocket** — and all
five newforms have Fricke eigenvalue +1, so their central values are not forced to
vanish, and indeed L(f,2) ≠ 0 for each (§8). The
−1 directions that do exist are all **d-oldform** combinations f(τ) − d²f(dτ), whose
W-eigenvalue is −ε(f) by the standard 2×2 computation on span{f(τ), f(dτ)}.

**Consequence for the two Apéry-perfect hosts.** N = 6 (Apéry's ζ(3) row, λ₂ = 17−12√2) and
N = 7 (Cooper's s₇, λ₂ = −1) have one-dimensional S₄ with ε = +1: no cusp-form source there,
and none can be manufactured (§5).

---

## 3. The five fold-regular cusp rows  **[verified]** (`03_cusprows.gp`, `04_foldreg.gp`)

Exact rational arithmetic, A_n and B_n peeled from the q-expansions to n ≤ 240 (n ≤ 280 in
`06_periods.gp`). Host recurrence Σⱼ Pⱼ(n)A_{n+j} = 0 fitted and verified; P₂(n) = (n+2)³ in
every row.

### 3.1 Row 1 — Cooper's s₁₀ host, N = 10, C = 25, B = 6, u = (η₅η₁₀/η₁η₂)²

* λ₁ = 16, λ₂ = −4; A = 1, 2, 18, 164, 1810, 21252, 263844, 3395016, …
* host recurrence: P₀ = −64x³−192x²−188x−60, P₁ = −12x³−54x²−82x−42, P₂ = (x+2)³, char. y²−12y−64
* Φ = f₅(τ) − 4f₅(2τ) = q − 8q² + 2q³ + 24q⁴ − 5q⁵ − 16q⁶ + …
* B = 0, 1, 5, 1460/27, 123425/216, 6106339/900, …  **k = 3** (min k with d_n^k B_n ∈ ℤ, n ≤ 240)
* companion's minimal recurrence: **order 3, degree 4**, char. y³−8y²−112y−256 = (y−16)(y+4)²
* ξ = **L(f₅,3)/2**, agreement of B_n/A_n with L(Φ,3) at n = 240: **144.09 digits**
* rate |λ₂/λ₁| = 1/4; K₋ = 0.20376662415349320849210895209670200943480…

  ξ = **0.3185322739492592814489078985464660740358072155633223741853…**

### 3.2 Row 2 — the Domb host, N = 12, C = 9, B = 10, u = (η₂η₃η₁₂/η₁η₄η₆)⁴

* λ₁ = 16, λ₂ = 4; A = 1, 4, 28, 256, 2716, 31504, 387136, 4951552, …
* host recurrence: P₀ = 64(x+1)³, P₁ = −20x³−90x²−138x−72, P₂ = (x+2)³, char. y²−20y+64
* Φ = f₆(τ) − 4f₆(2τ) = q − 6q² − 3q³ + 12q⁴ + 6q⁵ + 18q⁶ − 16q⁷ − …
* B = 0, 1, 37/4, 818/9, 141587/144, 51588841/4500, …  **k = 3**
* companion's minimal recurrence: **order 3, degree 4**, char. (y−16)(y−4)²
  (Q₀ = −256x⁴−1152x³−1920x²−1408x−384, Q₁ = 144x⁴+992x³+2628x²+3164x+1456,
   Q₂ = −24x⁴−220x³−762x²−1180x−688, Q₃ = x⁴+11x³+45x²+81x+54)
* ξ = **L(f₆,3)/2**, agreement at n = 240: **144.15 digits**
* rate 1/4; **K₋ = L(6.4.a.a, 2)/√π = 0.28757331148922160500069983443231806149557…** [verified, 50 digits]

  ξ = **0.3645567287833946297572856079724189522146819522618335207885…**

### 3.3 Row 3 — N = 12, C = 1, B = 34, deg u = 4

* λ₁ = 36, λ₂ = 32 (λ₂ᶰᵒʳᵐ = 32); A = 1, 12, 252, 6240, 167580, 4726512, …
* host recurrence: P₀ = 1152x³+3456x²+3168x+864, P₁ = −68x³−306x²−466x−240, char. y²−68y+1152
* same Φ = f₆(τ) − 4f₆(2τ); B = 0, 1, 133/4, 9332/9, 1547953/48, …  **k = 3**
* companion: **order 4, degree 6**, char. 4(y−36)²(y−32)²
* ξ = **L(f₆,3)/2** again (the limit is host-independent, §6); agreement at n = 240: **12.2 digits** (= 240·log₁₀(36/32))
* rate 8/9; K₋ = 0.37960819633459956111576090173690623070839… (not identified)

### 3.4 Rows 4, 5 — Cooper's s₁₈ host, N = 18, C = 1, B = 14, deg u = 4

* λ₁ = 16, λ₂ = 12; A = 1, 6, 54, 564, 6390, 76356, 948276, 12132504, …
* host recurrence: P₀ = 192x³+576x²+564x+180, P₁ = −28x³−126x²−194x−102, char. y²−28y+192

| | Φ = f₆(τ) − 9f₆(3τ) | Φ = f₉(τ) − 4f₉(2τ) |
|---|---|---|
| B₁…B₅ | 0, 1, 55/4, 3169/18, 108347/48, 29311423/1000 | 0, 1, 27/2, 171, 17415/8, 112509/4 |
| k | 3 | 3 |
| companion recurrence | order 4, deg 5, char. (y−16)²(y−12)² | order 3, deg 4, char. (y−16)(y−12)² |
| ξ | **2·L(f₆,3)/3** | **L(f₉,3)/2** |
| digits at n = 240 | 29.86 | 29.90 |
| K₋ | **√3·L(6.4.a.a,2)/√π** = 0.49809158640016257987401106018289485845958… | 0.43491056992916815669238152561817236709908… (not identified) |
| ξ (40+ digits) | **0.4860756383778595063430474772965586029529…** | **0.4594751313925473129216738852378661910708…** |

### 3.5 What is checked

For each row `04_foldreg.gp` prints, to 250-digit working precision:
L(Φ,2) = 0 (to 10^{−270}); L(Φ,1) + N·L(Φ,3)/(2π²) = 0 (to 10^{−251}); w_n/w_{n−1} → λ₂;
w_n·n^{3/2}/λ₂ⁿ → K₋; and the digit count of B_n/A_n vs L(Φ,3). **No free integration
anywhere**: k = 3 is minimal in all five rows, although three of the four hosts have a
meromorphic canonical source with k = 2 (`hostscan` §9).

---

## 4. The complete fold-regular source space, and the mixed periods  **[verified, ≤ 167 digits]** (`05_space.gp`, `06_periods.gp`)

V(N) := {Φ ∈ M₄(Γ₀(N)) : Φ\|W_N = −Φ, a₀(Φ) = 0}. Since M₄ = Eis ⊕ S₄ is W_N-stable and cusp
forms have a₀ = 0, V = V_Eis ⊕ S₄^{W=−1}.

| N | dim M₄ | dim V | V_Eis | V_cusp | periods carried (identified by lindep from the row, ≥ 84 digits) |
|---|---|---|---|---|---|
| 5 | 3 | **0** | 0 | 0 | — (the canonical source of this host is meromorphic) |
| 6 | 5 | **1** | 1 | 0 | ζ(3)/6 — **the Eisenstein line, and nothing else** |
| 7 | 3 | **0** | 0 | 0 | — |
| 8 | 5 | **1** | 1 | 0 | 7ζ(3)/32 |
| 9 | 5 | **1** | 1 | 0 | L(3,χ₋₃)/3 |
| 10 | 7 | **2** | 1 | 1 | 9ζ(3)/5, **L(5.4.a.a,3)/2** |
| 12 | 9 | **3** | 2 | 1 | −ζ(3)/24, 11ζ(3)/4, **L(6.4.a.a,3)/2** |
| 18 | 13 | **5** | 3 | 2 | −7ζ(3)/54, 85ζ(3)/54, L(3,χ₋₃)/2, **2L(6.4.a.a,3)/3**, **L(9.4.a.a,3)/2** |

Reading. (i) On **Apéry's own host** the fold-regular source space is one-dimensional: the
Eisenstein source, period ζ(3)/6. There is no second direction of any kind, cuspidal or not —
so no mixed pair, no new period, nothing. Same at N = 8 and N = 9. (ii) The Eisenstein part
carries only **one** interesting period per level (all ζ(3)-multiples at N = 12; ζ(3) and
L(3,χ₋₃) at N = 18) — consistent with the rigidity statement of `hostscan` §10.6.
(iii) **The cusp forms are what creates a genuine mixed configuration**: on the level-12 hosts
the pair {ζ(3), L(6.4.a.a,3)} and on the level-18 host the quadruple
{ζ(3), L(3,χ₋₃), L(6.4.a.a,3), L(9.4.a.a,3)} are all fold-regular on **one** host with **one**
row A_n — exactly the {ζ(2), L(2,χ₋₃)} configuration CDT exploit at Γ₀(6), one weight up.
The obstacle is λ₂ (§9), not the period supply.

---

## 5. The W_N = +1 obstruction, and why no admixture repairs it  **[exact]** + **[verified]** (`10_plus.gp`, `17_pluslimits.gp`)

**Not a scalar.** By §1.3, ε = +1 turns (★) into 2Θ = ξ(Nτ²+1) − P. The left side is a
non-polynomial q-series; there is no finite set of scalar conditions to impose and no
finite-dimensional correction. Equivalently: Φ/Φ₀ is a weight-zero function on X₀(N), and it
is a rational function of x iff ε = −1; for ε = +1 it is *anti*-invariant, so the
inhomogeneity of the row's operator acquires a branch point at the near fold.

**No Eisenstein admixture.** M₄(Γ₀(N)) = Eis ⊕ S₄ is W_N-stable, so for Φ = Φ_E + Φ_S the
+1 component is (Φ_E)⁺ + (Φ_S)⁺ with the two summands in complementary subspaces: they cannot
cancel. A cusp form with ε = +1 is unusable, full stop. (This is the honest negative answer to
the "single scalar condition ⇒ one Eisenstein direction cancels it" hope in the brief.)

**What is measured instead.** For ε = +1 the Apéry limit still exists — B − ξA merely fails to
lose the (x−t₁)^{1/2} branch — and the convergence degrades from geometric to a power law:

| host / Φ | ratio of successive differences of B_n/A_n at n = 200 | expected |
|---|---|---|
| every W_N = +1 form, all 12 hosts | 0.98999 … 0.99001 | ((n−1)/n)² = 0.990025, i.e. **B_n/A_n = ξ + c/n + O(1/n²)** |
| N = 18, W₁₈ = −1 | 0.750051 | λ₂/λ₁ = 3/4 |
| N = 12, C = 1, W₁₂ = −1 | 0.889035 | λ₂/λ₁ = 8/9 |
| N = 10, N = 12 C = 9, W = −1 | exactly 0 at 60 digits | (1/4)²⁰⁰ |

and the limit is **not** L(Φ,3). E.g. 6.4.a.a on all three level-6 hosts gives

  ξ = 0.264718537218080276566224788415170957070513623723586358426586… ,  L(Φ,3) = 0.729113457566789…

(60 digits, `16_limits.gp`). Two of the twelve ε = +1 limits *are* rational multiples of L(Φ,3)
by accident — Φ = f₆+4f₆(2τ) at N = 12 gives ξ = 4L(Φ,3)/9, and one of the level-18 plus forms
gives ξ = L(Φ,3)/2 — which is a useful caution: **an L-value limit does not imply
fold-regularity; fold-regularity is a statement about the rate.** The other ten limits resist
`lindep` against {L(Φ,1), L(Φ,2), L(Φ,3), π³} at 40 digits. **[verified, 60 digits]**

---

## 6. The Apéry limit does not see the host at all  **[exact]** + **[verified, 96 digits]** (`12_checks.gp`)

**Lemma A.** If F ∈ M₂(Γ₀(N)) with F\|₂W_N = −F, then F′(τ_c) = N τ_c F(τ_c) at τ_c = i/√N.
*Proof.* Differentiate F(στ) = −Nτ²F(τ) at the fixed point, using σ(τ_c+s) = τ_c − s + O(s²). ∎

Since the (x−t₁)^{1/2}-coefficient ratio at a fold is d(FΘ)/dF = Θ + F·Θ′/F′, Lemma A gives

> **ξ(Φ) = Θ(τ_c) + Θ′(τ_c)/(N τ_c) = Θ(q_c) + (2π/√N)·(DΘ)(q_c),  q_c = e^{−2π/√N},**

which contains **no reference to x** — the host chooses the *rate* and the *denominators*, never
the limit. Verified: DF(q_c)/F(q_c) = √N/(2π) to 50 digits on all twelve hosts; and the three
level-6 hosts, the two level-8 hosts and the two level-12 hosts each return **identical** limits
for the same Φ, to 60 digits.

Specialising to ε = −1 gives a rapidly convergent series for a critical L-value:

> **L(Φ,3) = Σ_{m≥1} (a(m)/m³)·e^{−2πm/√N} + (2π/√N)·Σ_{m≥1} (a(m)/m²)·e^{−2πm/√N}**
> for every Φ ∈ M₄(Γ₀(N)) with Φ\|W_N = −Φ and a(0) = 0.

Checked against `lfun` for all four −1 eigenforms: agreement to **96 digits** (the truncation
at m = 300). This is the classical Fricke-symmetric expansion, but it is what the fold formula
returns, which is a good consistency check on the whole geometric picture.

---

## 7. The source in closed form, and the recurrence  **[exact]** (`08_recur4.gp`, `09_exact.gp`)

Apply the host's own three-term recurrence to the cusp companion and look at the residual
res(n) := Σⱼ Pⱼ(n)B_{n+j}:

> **res(n) = binom(2n+2, n+1)·(λ₂/4)^{n+1} = [x^{n+1}](1−λ₂x)^{−1/2}** — verified exactly for
> all 0 ≤ n ≤ 177 for the three rows §3.1, §3.2, §3.4-right (Φ = f − 4f(2τ) in each case).

Unwinding through L(B) = p₃(x)·F·Φ/(Dx)³ with p₃ = x³P(x) and Φ₀ = xF²√P this is equivalent to

> **Φ = x·F²·√(1 − λ₁x) = Φ₀ / √(1 − λ₂x)** — verified as an identity of q-series to O(q¹⁶⁰).

So on these three hosts the weight-four **cusp form is the canonical (meromorphic) source with
the far-fold square root divided out** — which is exactly the statement that its W_N-anti-invariant
"odd" factor keeps only the near-fold branch. (The two remaining rows, Φ = f₆−9f₆(3τ) at N = 18
and Φ = f₆−4f₆(2τ) on the deg u = 4 host N = 12, C = 1, have a residual with the same
(1−λ₂x)^{−1/2} leading behaviour but no such clean closed form.)

**Minimal recurrences for B_n** (fitted exactly, verified to reproduce B_n for n ≤ 127):
order 3 degree 4 with characteristic polynomial (y−λ₁)(y−λ₂)² for the three clean rows, order 4
with (y−λ₁)²(y−λ₂)² (N = 18, f₆−9f₆(3τ)) and 4(y−λ₁)²(y−λ₂)² (N = 12, C = 1). **In no case does
B_n satisfy the row's three-term recurrence** — the double root λ₂ is exactly the
(1−λ₂x)^{∓1/2} pair of exponents at the far fold.

**Far singularity.** B − ξA is analytic at the near fold t₁ = 1/λ₁ and has a **square-root
branch** (not a log: both finite singularities of a Fricke host are folds, `hostscan` §8) at the
far fold t₂ = 1/λ₂, with ξA − B = c·(1−λ₂x)^{1/2} + [analytic], c = −2√π·K₋:
c = −0.72233…, −1.01942…, −1.34568…, −1.76569…, −1.54172… for the five rows.

---

## 8. The constants K₊, K₋  **[verified, 50 digits]** (`11_kminus.gp`, `11b_kminus_hi.gp`, `14_kclosed.gp`)

Recurrences iterated exactly to n = 1600, Richardson-extrapolated in 1/n with 22 nodes at
2400-digit working precision. K₊ (measured from A_n) reproduces the weight-two closed form
√N/(2π^{3/2})·√(λ₁/(λ₁−λ₂)) of `ASYMPTOTIC_CONSTANTS.md` Thm 3.4 to **50 digits** on every host —
which calibrates the extrapolator, so K₋ below is good to about the same.

| host | Φ | K₋ | closed form | K₊K₋(λ₁−λ₂) |
|---|---|---|---|---|
| N=10 C=25 B=6 | f₅−4f₅(2τ) | 0.2037666241534932084921089520967020094348 | — | 1.03503070796409682512860016 |
| N=12 C=9 B=10 | f₆−4f₆(2τ) | 0.2875733114892216050006998344323180614956 | **L(6.4.a.a,2)/√π** | 1.23946712184848171267869766 |
| N=12 C=1 B=34 | f₆−4f₆(2τ) | 0.3796081963345995611157609017369062307084 | — | 1.41694386550910601440379377 |
| N=18 C=1 B=14 | f₆−9f₆(3τ) | 0.4980915864001625798740110601828948584596 | **√3·L(6.4.a.a,2)/√π** | 1.51803100074242176837420827 |
| N=18 C=1 B=14 | f₉−4f₉(2τ) | 0.4349105699291681566923815256181723670991 | — | 1.32547456276972086793699040 |

Two remarks.

* The **Casoratian identity fails**: `ASYMPTOTIC_CONSTANTS.md` Thm 5.1 asserts
  K₊K₋(λ₁−λ₂) = κ = ∏(1+d/(cj²)) for two solutions of one three-term recurrence, and B_n is not
  one (§7). The five products above are not κ and are not each other.
* Where K₋ *is* identified it is the **central** value L(f,2) of the underlying newform divided
  by √π (ratio 1 and √3 to 68–70 digits, `14_kclosed.gp`). This is structurally right: K₋ is a
  period-polynomial evaluation at the **far** W_N-fixed point, where the vanishing
  L(Φ,2) = 0 of the source no longer applies and the newform's own non-vanishing central value
  survives. The other three resist `lindep` against {L(f,2), L(f,3), L(f,3)/π, π³} at 45 digits
  (`15_kclosed2.gp`) — presumably genuine modular symbols at the third W_N-fixed point class
  rather than critical values. **[estimated]**

---

## 9. The arithmetic budget: why there is no irrationality statement  **[exact]** (`12_checks.gp`)

Apéry's criterion for the linear form w_n = A_n·ξ − B_n with A_n ∈ ℤ, d_n³B_n ∈ ℤ,
|w_n| ≍ |λ₂|ⁿ n^{−3/2}, d_n³ ≍ e^{3n} needs **log|λ₂| + 3 < 0**.

| host | λ₂ | log\|λ₂\| + 3 | fold-regular cusp source? |
|---|---|---|---|
| **N=6 C=72 B=17 (Apéry)** | 17−12√2 = 0.029437… | **−0.5255** | **no** (S₄^{W=−1} = 0) |
| N=5 C=125 B=22 | 22−10√5 = −0.36068 | +1.980 | no |
| N=8 C=32 B=12 | 12−8√2 = 0.68629 | +2.624 | no |
| N=7 C=49 B=13 (Cooper s₇) | −1 | +3.000 | no |
| N=9 C=27 B=9 | 9−6√3 = −1.39230 | +3.331 | no |
| N=6 C=81, N=6 C=64, **N=10 C=25 (s₁₀)**, **N=12 C=9 (Domb)** | ∓4 | **+4.386** | last two: **yes** |
| **N=18 C=1 (s₁₈)** | 12 | +5.485 | **yes** |
| N=8 C=16 (Catalan) | 16 | +5.773 | no |
| **N=12 C=1** | 32 | +6.466 | **yes** |

The only host inside the budget is Apéry's own, and it is precisely the one with an empty
−1 eigenspace; the best fold-regular cusp row misses by **4.386 nats**. Nothing in the
inventory helps: k = 3 is minimal (§3), so there is no free integration to buy back 1 nat, and
the limit is host-independent (§6), so no re-coordinatisation moves L(f,3) onto a better λ₂ at
the same level. To reach the budget one would need a level with λ₂ᶰᵒʳᵐ = 1 **and**
S₄(Γ₀(N))^{W_N=−1} ≠ 0 — by §2 and `hostscan` §4 no such Fricke host exists at all.

---

## 10. Coverage — what this probe does and does not cover  **(honest ledger)**

Done, with the stated verification:

1. All twelve CDT-shape Fricke hosts of `hostscan` §4.4, every W_N-eigenvector of S₄(Γ₀(N)) on each (16 (host, form) pairs computed in exact arithmetic to n = 240–280).
2. Oldform combinations and Fricke-sign combinations: automatically, since the W_N-eigenspaces of the *full* S₄(Γ₀(N)) are used, not just newforms. Every fold-regular direction found is an oldform combination.
3. Eisenstein + cusp mixtures: settled negatively as a *proof* (§5), not merely a search.
4. The full holomorphic fold-regular source space V(N) including its Eisenstein part, with all periods identified (§4).

Not done:

* **The secondary weight-one / Γ₁(7) question is NOT done.** The weight-one hosts have both
  finite singularities at **cusps** (unipotent monodromy, `hostscan` §8), so §1's fold argument
  does not apply at all: the governing statement there is Theorem 5.2 of
  `ASYMPTOTIC_CONSTANTS.md` (constant term at the cusp × constant term of F), and a cusp form's
  constant terms vanish at every cusp, so the analysis has to be redone from scratch with the
  subleading terms. In addition, `hostscan` §5 records that the only level-7 weight-one
  four-point host has a **complex-conjugate pair** of dominant singularities
  (λ = −13/2 ± (3√3/2)i, λ₂ᶰᵒʳᵐ = 7), so its row has no Apéry limit in the usual sense and the
  CM form η(τ)³η(7τ)³ has nowhere obvious to act. Flagged as open.
* **K₋ in three of five cases** is only [verified, 50 digits], not identified (§8).
* Meromorphic sources (Cooper-type, with poles at CM points) are outside Theorem A's hypotheses:
  the Eichler-integral argument needs Φ holomorphic on ℍ. Whether a meromorphic W_N-anti-invariant
  source can be fold-regular is not addressed here (the canonical ones of the free-integration
  hosts are, by construction, but that is not a consequence of Theorem A).
* Levels beyond the twelve-host census (N ≤ 120 with λ₂ᶰᵒʳᵐ ≤ 1.05 was already swept in
  `hostscan`; a source-side sweep at larger λ₂ᶰᵒʳᵐ would only make §9 worse).
* No attempt at a *p*-adic or Lucas-congruence analysis of the cusp-form companions (the
  companion formula B_n = Σ c(m)m^{−3}e_{n,m} of `COMPANION_ARITHMETIC.md` applies verbatim with
  c(m) = a(m) the Hecke eigenvalues, and hypothesis (a) of its Theorem 5.1 is *false* for a cusp
  form, so the companion Lucas law is expected to fail — untested).

---

## 11. Scripts and data

`gp -q <file>`; every script terminates in seconds (the slowest, `11b_kminus_hi.gp`, is 5 s).

| script | what it does | output |
|---|---|---|
| `lib.gp` | eta/F series, x-peeling, D^{−r}, denominator exponent, recurrence fitting, W_N eigenspaces | — |
| `hosts.gp` | the twelve CDT-shape Fricke hosts (N, C, B, divisors, exponents) | — |
| `01_forms.gp` | W_N matrix and ± eigenspaces of S₄(Γ₀(N)), N = 5,…,18, with q-expansions | `01_forms.out` |
| `02_calib.gp` | calibration: canonical Eisenstein source on Apéry's host ⇒ ξ = ζ(3)/6, k = 3 | `02_calib.out` |
| `03_cusprows.gp` | exact A_n, B_n (n ≤ 180) for every host × every eigenform; k(den); ratios | `03_cusprows.out` |
| `04_foldreg.gp` | the five fold-regular rows: L(Φ,1),L(Φ,2),L(Φ,3), the two FE identities, rate, K₋, digit counts | `04_foldreg.out` |
| `05_space.gp` | dim V(N) = dim{Φ ∈ M₄ : Φ\|W_N = −Φ, a₀ = 0}, split Eis/cusp, with bases | `05_space.out` |
| `06_periods.gp` | all periods of all fold-regular sources, read off the row to ≤ 167 digits and identified by `lindep` | `06_periods.out` |
| `07_recur.gp` | residual of the host's three-term recurrence applied to B | `07_recur.out` |
| `08_recur4.gp` | minimal recurrence of B_n (order/degree/char. poly), residual vs binom(2n,n)(λ₂/4)ⁿ | `08_recur4.out` |
| `09_exact.gp` | exact inhomogeneity for n ≤ 147, and the identity Φ = x F²√(1−λ₁x) to O(q¹⁶⁰) | `09_exact.out` |
| `10_plus.gp` | W_N = +1: Richardson limit vs the fold formula, and the O(1/n) rate diagnostic | `10_plus.out` |
| `11_kminus.gp`, `11b_kminus_hi.gp` | K₊ (vs closed form) and K₋ by iterating to n = 800 / 1600 | `11_kminus.out`, `11b_kminus_hi.out` |
| `12_checks.gp` | Lemma A, the rapidly convergent L(Φ,3) series, the irrationality budget | `12_checks.out` |
| `13_kident.gp`, `14_kclosed.gp`, `15_kclosed2.gp` | identification of K₋ | `*.out` |
| `16_limits.gp` | the Apéry limit of every weight-four cusp form by level, 60 digits | `16_limits.out` |
| `17_pluslimits.gp` | `lindep` on the W_N = +1 limits | `17_pluslimits.out` |
| `18_xfold.gp`, `19_ufold.gp` | the near fold is i/√N: u(i/√N) = 1/√C to 60 digits (and the x-series artefact) | `*.out` |
| `data_*.txt` | exact A_n, B_n vectors written by `03_cusprows.gp` | — |
