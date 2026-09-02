# The machine host sweep: every genus-zero modular Apéry host of level ≤ 120 with
# λ₂ᶰᵒʳᵐ ≤ 1.05, its singular geometry, its Eisenstein source and its CDT score

*Working directory `hostscan/`. All scripts and raw outputs listed in §12.
Conventions: `CDT_FINDER.md` §§1–6 (τ, ceil, entryC, entryR, margin),
`ASYMPTOTIC_CONSTANTS.md` Thm 3.4 (weight-two structure),
`INVENTORY_BOUND.md` (u₁ = 1 forced, so the "best inventory" scenarios are closed),
`NUMBER_FIELD_HOLONOMY.md` (the number-field budget is the **average** over places).
Tags: **[exact]** = closed-form / exact rational arithmetic; **[verified]** = high-precision
computation in this task (digit counts stated); **[cited]**. No irrationality claim is made.*

---

## 0. Verdict

| claim | verdict |
|---|---|
| **λ₂ᶰᵒʳᵐ ≤ 1.05 ⟺ λ₂ is an algebraic unit.** λ₂ is an algebraic integer (integrality of the pure module), so \|N(λ₂)\|^{1/deg} ≥ 1 unless λ₂ = 0. The threshold 1.05 therefore selects exactly the Apéry-perfect hosts. | **[exact]** §1 |
| A **four-singular-point** ("CDT-shape") host requires a genus-zero group with exactly **four special points** (cusps + elliptic points). Over Γ₀(N) these are exactly **N = 5, 6, 7, 8, 9**; adding Γ₁(5) and the Atkin–Lehner quotients gives the complete list of §2. | **[verified]** §2 |
| The **Möbius freedom x ↦ x/(1+βx) is not free**: for β ≠ 0 the pole of x leaves the special locus and a *fifth* singular point appears, so the three-term recurrence and CDT's normaliser descent are both destroyed. The admissible normalisations are the finitely many that place the pole of x on a special point. Enumerated exactly. | **[verified]** §3 |
| The **complete** census of CDT-shape hosts carrying a weight-two form (third-order / Apéry-like rows, "Fricke family") over **all** eta-quotient Fricke parameters with N ≤ 120 and deg u ≤ 4 is **12 hosts**, listed in §4. Exactly **two** have λ₂ᶰᵒʳᵐ = 1: Apéry's own ζ(3) row (k = 3) and Cooper's s₇ (k = 2, period ζ(2)/7). Both are already in the finder's table. | **[verified]** §4 |
| **Five** of those twelve are **not in the project ledger** (they do not appear in `lattice/sporadic_scan2/table.json`, whose coordinate box enumerates eta-quotient t, while the Fricke coordinate x = u/(1+Bu+Cu²) is a rational function of one); all five have λ₂ᶰᵒʳᵐ ≥ 4. Their periods (identified to ≥ 200 digits) are ζ(2)/10 (N = 5), ζ(2)/8 and L(2,χ₋₃)/4 (N = 6), **G/4 — a second Catalan row** (N = 8, λ₂ᶰᵒʳᵐ = 16), 5L(2,χ₋₃)/16 (N = 12). | **[verified]** §4 |
| The **weight-one** (CDT/Zagier) CDT-shape census is nine hosts on Γ₀(5..9) and Γ₁(5). Exactly **three** have λ₂ᶰᵒʳᵐ = 1: Zagier C (= CDT's own, L(2,χ₋₃)/2), Zagier A (ζ(2)/4), Zagier D (ζ(2)/5). No fourth exists at N ≤ 60. | **[verified]** §5 |
| **No new HOST with λ₂ᶰᵒʳᵐ ≤ 1.05 exists**: the four Apéry-perfect k = 2 hosts of the finder's table are the complete list. | **[verified]** §§4–6 |
| **But one of them carries a new PERIOD the finder's census missed.** On Zagier D's host Γ₁(5) (λ₂ᶰᵒʳᵐ = 1, k = 2, margin **+0.0053**) the fold-regular weight-three Eisenstein space is **2-dimensional over ℚ(√5)** and carries, besides ζ(2)/5, the value **ξ = φ⁵·Im L(2,ψ₄) − Re L(2,ψ₄) = 0.655634188840656766…**, where ψ₄ is the odd **quartic** character mod 5 — verified to **128 digits**, with **k = 2**. This is both a new single L-value and a mixed pair {ζ(2), ξ} on one Apéry-perfect host. | **[verified]** §10 |
| Three *additional* configurations with \|λ₂\| = 1 exist off the CDT shape — N = 10 (λ₂ = 1, double), N = 12 (λ₂ = −1, double), N = 18 (λ₂ a primitive 6th root of unity, double) — but they have **five, five and seven** singular points; the far singularity is an interior order-2 orbifold point of the covering group, not one of the two Fricke folds, and the descent has nothing to remove there. The N = 10 one is scan 2's row #109 (k = 3). | **[verified]** §7 |
| **Cusp vs. fold dichotomy.** Weight-one hosts (Γ₀(N), Γ₁(N)) put *both* finite singularities at **cusps** (log monodromy — CDT's own case). Fricke weight-two hosts put both at **order-2 orbifold points** (folds, ℤ/2 monodromy). Cooper's s₇ is the only Apéry-perfect k = 2 host of fold type. | **[verified]** §8 |
| The **M₃/M₄ Eisenstein census** the finder's §8 lists as not done was run (3996 directions, N ≤ 60, k′ ∈ {3,4}; **0 mismatches** against PARI on 1103 nebentypus components). Its independent verdicts: *within one nebentypus there is exactly one interesting period*, in all 1102 components — so mixed pairs always need a non-eigenform on Γ₁(N); at weight 3 the pair (π², L(2,ψ)) lives inside one nebentypus and is generic; at weight 4 (ζ(3), L(3,χ_D)) always straddles two. Its annihilation table gives **dim(fold-regular) = 2 at N = 5, k′ = 3 with periods the two quartic L(2,ψ₄)** — §10 confirmed by a wholly different method. | **[computed]** §10.6, `eis/EIS_REPORT.md` |
| The finder's own numbers are reproduced exactly by this sweep's scorer: Zagier C margin **+0.0053**, Zagier E entryC **−0.0766**, Apéry entryC **−0.435**, Cooper s₁₈ entryR **−1.638**. | **[verified]** §6 |

**One sentence.** *Running the sweep the finder did not run — a machine enumeration of every
genus-zero host of level ≤ 120 with the exact images of all its cusps and elliptic points —
turns up twelve four-point Fricke hosts and nine four-point weight-one hosts, five of them new
to the ledger but all on geometries at least four times worse than CDT's; the Apéry-perfect
pocket is exactly the four hosts the finder already had; and the one thing that is genuinely
new is a **period**, not a host: on Zagier D's Γ₁(5) — margin +0.005, k = 2 — the fold-regular
Eisenstein source space is two-dimensional over ℚ(√5) and its second direction carries
φ⁵·Im L(2,ψ₄) − Re L(2,ψ₄) for the odd quartic character mod 5, a value that is neither ζ(2)
nor a quadratic L(2,χ₋D) and that the finder's character-theoretic census could not see.*

---

## 1. The threshold λ₂ᶰᵒʳᵐ ≤ 1.05 is the unit condition  **[exact]**

`CDT_FINDER.md` §3: integrality of the pure functions Li_j(x/s) = Σ(λ₂x)ⁿ/nʲ forces λ₂ to be
an **algebraic integer**, and the number-field normalisation is λ₂ᶰᵒʳᵐ = |N(λ₂)|^{1/[K:ℚ]}.
A non-zero algebraic integer has |N(λ₂)| ≥ 1, so

  λ₂ᶰᵒʳᵐ ≤ 1.05  ⟺  |N(λ₂)| = 1  ⟺  **λ₂ is a unit**  (the Apéry-perfect condition),

and there is nothing strictly between 1 and √2 = 1.414 (the next value, |N| = 2, deg 2). So
"1.05" and "1" select the same set; the sweep below reports λ₂ᶰᵒʳᵐ for every host and the
next value above 1 that actually occurs is **2** (level-16 weight-one host) and then 4.

Concretely, for a three-term row with characteristic polynomial t² − at + c:
* c a perfect square (λᵢ ∈ ℤ): unit ⟺ λ₂ = ±1;
* c not a square (λᵢ conjugate over a real quadratic field): unit ⟺ **c = λ₁λ₂ = ±1**.

---

## 2. Genus-zero groups and their special points  **[verified]** (`01_groups.gp`)

For every N ≤ 120 and every subgroup S of the Atkin–Lehner group generated by the prime-power
exact divisors, `01_groups.gp` computes genus(X₀(N)/W_S) as the joint (+1)-eigenspace
dimension of {W_Q} on S₂(Γ₀(N)) (`mfatkininit`), the number ν′_∞ of W_S-orbits on the cusps,
and the orbifold defect A := μ′/6 + 2 − ν′_∞ = Σᵢ(1 − 1/mᵢ). It reproduces the classical
Fricke genus-zero list exactly. The number of **special points** (the singular locus of any
weight-w local system on the quotient) is ν′_∞ + #{cone points}.

For the plain groups Γ₀(N) the count is ν_∞ + ν₂ + ν₃:

| N | 1 | 2 | 3 | 4 | **5** | **6** | **7** | **8** | **9** | 10 | 12 | 13 | 16 | 18 | 25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| #Σ | 3 | 3 | 3 | 3 | **4** | **4** | **4** | **4** | **4** | 6 | 6 | 6 | 6 | 8 | 8 |

so **the four-point (CDT-shape) Γ₀(N) hosts are exactly N = 5, 6, 7, 8, 9**, and N ≤ 4 are the
degenerate three-point (hypergeometric) cases. Adding Γ₁(5) (4 cusps, genus 0, no elliptic
points — Zagier D's home) and the four-point Atkin–Lehner quotients
Γ₀(5)+5, Γ₀(6)+2, Γ₀(6)+3, Γ₀(6)+6, Γ₀(6)+{2,3}, Γ₀(8)+8, Γ₀(9)+9, Γ₀(10)+{2,5},
Γ₀(12)+4, Γ₀(12)+{3,4} completes the four-point list at N ≤ 60. Everything else has ≥ 5
singular points, which is enough on its own to break the descent architecture (§3).

`17_allsub.gp` repeats this over **every** subgroup of the Atkin–Lehner group (not only
those generated by prime-power involutions, which is what `01_groups.gp` does and which
misses ⟨W₆⟩ ⊂ ⟨W₂,W₃⟩ and its analogues — the Fricke groups Γ₀(N)+N themselves). Filtering
on ν′_∞ + #cone = 4, using ⌈6A/5⌉ ≤ #cone ≤ 2A, the complete candidate list at N ≤ 60 is

```
N= 5 S=[]            N= 5 S=[5]         N= 6 S=[]          N= 6 S=[2]      N= 6 S=[3]
N= 6 S=[6]           N= 6 S=[2,3,6]     N= 7 S=[]          N= 7 S=[7]      N= 8 S=[]
N= 8 S=[8]           N= 9 S=[]          N= 9 S=[9]         N=10 S=[2]      N=10 S=[5]
N=10 S=[10]          N=10 S=[2,5,10]    N=11 S=[11]        N=12 S=[4]      N=12 S=[3,4,12]
N=13 S=[13]          N=14 S=[2,7,14]    N=15 S=[3,5,15]    N=17 S=[17]     N=18 S=[2,9,18]
N=20 S=[4,5,20]      N=21 S=[3,7,21]    N=22 S=[2,11,22]   N=30 S=[2,3,5,…,30]
```

together with Γ₁(5). Every host of §§4–5 lies in this list, and the list is a *superset*:
e.g. Γ₀(7)+7 has ν′_∞ = 1 and A = 5/3 = 2/3 + 1/2 + 1/2, i.e. one order-3 point (the two
ℤ/3 elliptic points of Γ₀(7) merged by W₇, which is where x = ∞ sits) and the two order-2
W₇-folds — Cooper's s₇ host, Σ = {0, 1/27, −1, ∞}, exactly four points **[verified]**.

---

## 3. The Möbius freedom is **not** free  **[verified]** (`05_rows.gp`, `08_scan3term.gp`, `09_fit.gp`)

`CDT_FINDER.md` and the task statement both invite a search over the integral normalisations
x ↦ x/(1 + βx), β ∈ ℤ, "to maximise the distance of the far singularity". The sweep shows
this search has (essentially) **one point in it** per pole placement:

Let the four special points be P₀ (x = 0, the cusp at i∞), P_pole and P₁, P₂. Writing
μ = 1/x, the shift x ↦ x/(1+βx) translates every μ by β, so λ₁ − λ₂ is invariant and the
pole of x_β sits at x = −1/β. **For β ≠ 0 that point is an ordinary point of the curve, not a
special point**, hence P_pole acquires a *finite* x-value and Σ gains a fifth element: the
row's minimal recurrence lengthens (verified: at Γ₀(10) the coordinate with pole at the cusp
pair {1/2, 1/5} has minimal recurrence of order **4** with quartic coefficients, char.
polynomial y⁴ − 20y³ + 38y² − 20y + 1 = (y² − 18y + 1)(y − 1)², five singular points), and
CDT's normaliser descent — which needs Σ = {0, t₁, t₂, ∞} exactly — no longer applies.

So the admissible normalisations are exactly the finitely many placements of the pole of x on
a special point. With Σ = {P₀, P₁, P₂, P₃} and μᵢ = 1/x_ref(Pᵢ) for a reference Hauptmodul,
the placement "pole at Pⱼ" gives

  **λᵢ = μᵢ − μⱼ  (i ≠ j)**,   an exact and complete description of the Möbius orbit.

Worked example (level 6, reference u = η₂η₆⁵/(η₁⁵η₃), C = 72): μ(cusp 0) = 0,
μ(cusp 1/2) = −9, μ(cusp 1/3) = −8. The three placements give
λ = (−9, −8) → **Zagier F** (λ₂ᶰᵒʳᵐ = 8), λ = (9, 1) → **Zagier C = CDT's host**
(λ₂ᶰᵒʳᵐ = 1), λ = (8, −1) → **Zagier A** (λ₂ᶰᵒʳᵐ = 1). CDT's placement is already the best
of the three; there is no fourth.

**Consequence.** A "β-search" cannot improve any host. The optimisation the finder left open
is empty.

---

## 4. The Fricke family (weight-two form, third-order rows): complete census

### 4.1 The construction and its exhaustive enumeration  **[exact]** (`02_fricke.gp`)

`ASYMPTOTIC_CONSTANTS.md` Thm 3.4: with u = Π_d η(dτ)^{r_d} an eta quotient whose exponent
vector is **W_N-antisymmetric** (r_{N/d} = −r_d, which forces Σr_d = 0 and u|W_N = 1/(Cu)),
ord_∞ u = 1 and C⁻² = Π_d (N/d)^{r_d} ∈ ℚ, one sets

  x = u/(1 + Bu + Cu²),  F = D log u ∈ M₂(Γ₀(N)),  aₙ = [xⁿ]F,  λ₁,₂ = B ± 2√C, c = B² − 4C.

The enumeration is by the (antisymmetric) **cusp-order vector** through the Ligozat matrix
A_{c,d} = N/(24 gcd(c,N/c) c) · gcd(c,d)²/d, with deg u = Σ_c φ(gcd(c,N/c)) max(ord_c,0) ≤ 4.
This is exhaustive and takes seconds. Result: **116 families (N, r, C)** for N ≤ 120,
one of degree 1 for each genus-zero level N ∈ {2,…,10,12,13,16,18} with
C = 4096, 729, 256, 125, 72, 49, 32, 27, 20, 12, 13, 8, 6 respectively.
Calibration: all seven weight-two rows of `ASYMPTOTIC_CONSTANTS.md` §4.1 are reproduced with
their exact (u, B, C) — Apéry (6, 17, 72), ε (8, 12, 32), ζ (9, 9, 27), Cooper s₇ (7, 13, 49),
s₁₀ (10, 6, 25), s₁₈ (18, 14, 1), Domb α (12, 10, 9).

### 4.2 The exact geometry: images of the cusps and folds  **[verified, 60 digits]** (`06_special.gp`, `14_geom.gp`)

u is evaluated with the full Dedekind η at γ(90i) for a cusp matrix γ (error < 10⁻³⁰) and at
the elliptic points τ = (−d+i)/N (d² ≡ −1) and ((1−2d)+i√3)/(2N) (d² − d + 1 ≡ 0); every value
is recognised exactly by `bestappr`. **The two finite singular points of the row are the two
W_N-fixed points u = ±1/√C** (folds), and x = ∞ sits on the W_N-orbit {v, 1/(Cv)} of a
special point, which is what determines B:

  **B = −(C v + 1/v)  must be a rational integer.**

| N | C | special values u ≠ 0, ∞ | admissible B | c = λ₁λ₂ | λ₂ᶰᵒʳᵐ |
|---|---|---|---|---|---|
| 2 | 4096 | ell₂ at −1/64 (W₂-**fixed**) | 128 | 0 | degenerate |
| 3 | 729 | ell₃ at −1/27 (W₃-fixed) | 54 | 0 | degenerate |
| 4 | 256 | cusp 1/2 at −1/16 (W₄-fixed) | 32 | 0 | degenerate |
| 5 | 125 | ell₂ at (−11 ± 2i)/125 | 22 | −16 | 4 |
| **6** | 72 | cusps 1/2, 1/3 at −1/9, −1/8 | **17** | **1** | **1** (Apéry) |
| **7** | 49 | ell₃ at −13/98 ± (3√3/98)i | **13** | **−27** | **1** (Cooper s₇) |
| 8 | 32 | cusps 1/2, 1/4 at −1/8, −1/4 | 12 | 16 | 4 |
| 9 | 27 | cusps 1/3, 2/3 at −1/6 ± (√3/18)i | 9 | −27 | 5.196 |
| 10 | 20 | cusps 1/2, 1/5 at −1/5, −1/4; ell₂ at −1/5 ± i/10 | 9; 8 | 1; −16 | 1; 4 |
| 12 | 12 | cusps at −1/6, −1/4, −1/3, −1/2 | 8; 7 | 16; 1 | 4; 1 |
| 13 | 13 | ell₂ at (−3±2i)/13; ell₃ at −5/26 ± … | 6; 5 | −16; −27 | 4; 5.196 |
| 16 | 8 | cusps at −1/4, (−1±i)/4, −1/2 | 6; 4 | 4; −16 | 2; 4 |
| 18 | 6 | cusps at −1/3, −1/2, −1/4 ± …, −1/2 ± … | 5 | 1 | 1 |

(the three "degenerate" rows have the two folds collapsing onto one point: Γ₀(N)+N then has
only **three** special points and λ₂ = 0.)

### 4.3 Which of these are genuinely CDT-shape  **[verified]** (`08_scan3term.gp`)

A brute-force machine scan over **all 116 families × B ∈ [1, 90]** testing for a three-term
recurrence with cubic coefficients returns exactly **12 hits**, and they are exactly the
(family, B) pairs of §4.2 whose *pole pair and both folds are four distinct points*. The
unit configurations at N = 10 (B = 9), N = 12 (B = 7) and N = 18 (B = 5) are **not** among
them: their coordinate has a fifth (resp. seventh) singular point (§7).

### 4.4 The twelve CDT-shape Fricke rows: arithmetic and periods
**[verified: recurrence exact; k tested to n = 150; ξ to 200–250 digits]** (`10_full.gp`, `13_lindep.gp`)

For each row the minimal three-term recurrence Σⱼ Pⱼ(n)a_{n+j} = 0 (deg Pⱼ = 3) is fitted
exactly and *verified to reproduce the peeled coefficients*; the companion B₀ = 0, B₁ = 1 is
iterated; k = min{k : d_nᵏBₙ ∈ ℤ} is measured; ξ = lim Bₙ/Aₙ is computed at n = 330 and
identified. Calibration: N = 6, C = 72, B = 17 returns Apéry's own row 1, 5, 73, 1445, 33001,
Apéry's recurrence (n+1)³A_{n+1} = (2n+1)(17n²+17n+5)Aₙ − n³A_{n−1}, k = 3, and
ξ = ζ(3)/6 to 250 digits.

| # | N | C | B | deg u | λ₁ | λ₂ | λ₂ᶰᵒʳᵐ | free int. | k | period ξ | in ledger? |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 5 | 125 | 22 | 1 | 22+10√5 | 22−10√5 | 4 | yes | 2 | **ζ(2)/10** | new |
| 2 | 6 | 81 | 14 | 2 | 32 | −4 | 4 | yes | 2 | **ζ(2)/8** | new |
| 3 | 6 | 72 | 17 | 1 | 17+12√2 | 17−12√2 | **1** | no | 3 | ζ(3)/6 | Apéry |
| 4 | 6 | 64 | 20 | 2 | 36 | 4 | 4 | yes | 2 | **L(2,χ₋₃)/4** | new |
| 5 | 7 | 49 | 13 | 1 | 27 | −1 | **1** | yes | 2 | ζ(2)/7 | Cooper s₇ |
| 6 | 8 | 32 | 12 | 1 | 12+8√2 | 12−8√2 | 4 | no | 3 | 7ζ(3)/32 | AZ ε |
| 7 | 8 | 16 | 24 | 2 | 32 | 16 | 16 | yes | 2 | **G/4 (Catalan)** | new |
| 8 | 9 | 27 | 9 | 1 | 9+6√3 | 9−6√3 | 5.196 | no | 3 | L(3,χ₋₃)/3 | AZ ζ |
| 9 | 10 | 25 | 6 | 2 | 16 | −4 | 4 | yes | 2 | ζ(2)/5 | Cooper s₁₀ |
| 10 | 12 | 9 | 10 | 2 | 16 | 4 | 4 | no | 3 | 7ζ(3)/24 | Domb α |
| 11 | 12 | 1 | 34 | 4 | 36 | 32 | 32 | yes | 2 | **5L(2,χ₋₃)/16** | new |
| 12 | 18 | 1 | 14 | 4 | 16 | 12 | 12 | yes | 2 | L(2,χ₋₃)/2 | Cooper s₁₈ |

Every period is an exact rational multiple of a single L-value; no mixed pair occurs anywhere
in this family, as §5.3 of `CDT_FINDER.md` predicts for a trivial-nebentypus host.

**Reading.** (i) The **free integration** (k = 3 → k = 2, and the accompanying weight drop
ζ(3) ↦ ζ(2), L(3,·) ↦ L(2,·)) is far more common than the ledger suggested: it occurs in
**8 of the 12** rows, not only Cooper's three. (ii) Exactly two rows are Apéry-perfect, and
they are the two already known. (iii) Row 7 is a **second Catalan row** (ξ = G/4) with k = 2
— but λ₂ᶰᵒʳᵐ = 16 against Zagier E's 4, so it is 1.386 nats *worse*, not better.
(iv) Row 4 gives L(2,χ₋₃)/4 with k = 2 at λ₂ᶰᵒʳᵐ = 4: a fold-type analogue of CDT's target,
four times too far out.

---

## 5. The weight-one family (CDT's own architecture): complete four-point census
**[verified]** (`14_geom.gp` for the μ-configurations; §3 for the placement rule)

Here x is a Hauptmodul of a four-point genus-zero group, F ∈ M₁(Γ, χ) with χ odd, and both
finite singularities are **cusps**. Using λᵢ = μᵢ − μⱼ (§3) with the μ's read off §4.2:

| group | μ-configuration (μ = 1/x_ref) | pole at | λ₁, λ₂ | λ₂ᶰᵒʳᵐ | period ξ | name |
|---|---|---|---|---|---|---|
| Γ₀(5) | 0, −11∓2i | cusp 0 | −11±2i (complex pair) | 11.18 | complex fold | (AZ η region) |
| Γ₀(6) | 0, −9, −8 | cusp 1/2 | **9, 1** | **1** | L(2,χ₋₃)/2 | **Zagier C = CDT** |
| Γ₀(6) | " | cusp 1/3 | **8, −1** | **1** | ζ(2)/4 | **Zagier A** |
| Γ₀(6) | " | cusp 0 | 9, 8 | 8 | 5L(2,χ₋₃)/8 | Zagier F |
| Γ₀(7) | 0, −13/2 ∓ (3√3/2)i | cusp 0 | −13/2 ± (3√3/2)i | 7 | complex fold | — |
| Γ₀(8) | 0, −8, −4 | cusp 0 or 1/2 | 8, 4 | 4 | **G = Catalan** | **Zagier E** |
| Γ₀(8) | " | cusp 1/4 | 4, −4 | 4 | two dominant sing. | degenerate |
| Γ₀(9) | 0, −9/2 ∓ (3√3/2)i | cusp 0 | −9/2 ± (3√3/2)i | 5.196 | complex fold | Zagier B |
| Γ₁(5) | (4 cusps) | — | (11±5√5)/2 | **1** | ζ(2)/5 | **Zagier D** |

There are exactly **three** Apéry-perfect weight-one hosts, and they are CDT's own host, its
sibling orientation Zagier A, and Zagier D. This is the finder's table §4A verbatim; **the
sweep adds nothing to it**, and the reason is structural: four-point genus-zero groups are
scarce (§2), and at every one of them the three pole placements are exhausted above.

The Atkin–Lehner four-point quotients of §2 carry no rational weight-one row: a weight-one
form needs an odd nebentypus, and `CATALAN_AL_HOSTS.md` §2 shows the W_Q-eigenvectors are
defined over ℚ **iff every Q ∈ S is a perfect square**, which among the four-point quotients
leaves only Γ₀(12)+4 and Γ₀(9)+9, and Γ₀(12)+4 was shown there to carry no row at all
(no F ∈ M₁(12, χ₋₄) with F·𝒟t holomorphic, for any placement) **[cited]**.

---

## 6. Scoring: the ranked table  **[verified]** (`11_score.py`, using `lattice/cdt_finder/cdt_bound.py`)

CDT's inventory, m = 14, p₀ = 7, u = (1,3) for k = 2 and (1,3,5) for k = 3 (the finder's
proportional rule, now *forced* at u₁ = 1 by `INVENTORY_BOUND.md` Cor. 2.2), Σeᵢ = 6, bⱼ = 2;
ceil = log(256/λ₂ᶰᵒʳᵐ), entryC = ceil − τ, entryR = ceil + log 0.62922 − τ,
margin = m·entryR − (11.845 + log s).

```
host                                           w k  lam2^norm          field     tau    ceil  entryC  entryR   margin  period
-----------------------------------------------------------------------------------------------------------------------------
Fricke N=7 C=49 B=13  (Cooper s_7)             2 2     1.0000              Q   4.235   5.545  +1.310  +0.846   +0.005  zeta(2)/7
Gamma_0(6) pole 1/2  = Zagier C                1 2     1.0000              Q   4.235   5.545  +1.310  +0.846   +0.005  L(2,chi_-3)/2  [CDT]
Gamma_0(6) pole 1/3  = Zagier A                1 2     1.0000              Q   4.235   5.545  +1.310  +0.846   +0.005  zeta(2)/4
Gamma_1(5)           = Zagier D                1 2     1.0000  Q(sqrt5) unit   4.235   5.545  +1.310  +0.846   +0.005  zeta(2)/5
Fricke N=5 C=125 B=22                          2 2     4.0000     Q(sqrt125)   4.235   4.159  -0.077  -0.540  -18.017  zeta(2)/10
Fricke N=6 C=81  B=14                          2 2     4.0000              Q   4.235   4.159  -0.077  -0.540  -18.017  zeta(2)/8
Fricke N=6 C=64  B=20                          2 2     4.0000              Q   4.235   4.159  -0.077  -0.540  -18.017  L(2,chi_-3)/4
Fricke N=10 C=25 B=6  (Cooper s_10)            2 2     4.0000              Q   4.235   4.159  -0.077  -0.540  -18.017  zeta(2)/5
Gamma_0(8) pole 0/1/2 = Zagier E               1 2     4.0000              Q   4.235   4.159  -0.077  -0.540  -18.017  G = Catalan
Gamma_0(8) pole 1/4                            1 2     4.0000              Q   4.235   4.159  -0.077  -0.540  -18.017  (degenerate)
Gamma_0(9) pole 0    = Zagier B                1 2     5.1962      Q(sqrt-3)   4.235   3.897  -0.338  -0.801  -21.418  (complex fold)
Fricke N=6 C=72 B=17  (Apery)                  2 3     1.0000       Q(sqrt2)   5.980   5.545  -0.435  -0.898  -24.423  zeta(3)/6
Gamma_0(7) pole 0                              1 2     7.0000      Q(sqrt-3)   4.235   3.599  -0.636  -1.099  -25.292  (complex fold)
Gamma_0(6) pole 0    = Zagier F                1 2     8.0000              Q   4.235   3.466  -0.770  -1.233  -27.027  5 L(2,chi_-3)/8
Gamma_0(5) pole 0                              1 2    11.1803    Q(i sqrt 5)   4.235   3.131  -1.104  -1.568  -31.379  (complex fold)
Fricke N=18 C=1 B=14  (Cooper s_18)            2 2    12.0000              Q   4.235   3.060  -1.175  -1.638  -32.298  L(2,chi_-3)/2
Fricke N=8 C=16 B=24                           2 2    16.0000              Q   4.235   2.773  -1.463  -1.926  -36.038  G/4  (Catalan)
Fricke N=8 C=32 B=12  (AZ eps)                 2 3     4.0000       Q(sqrt2)   5.980   4.159  -1.821  -2.285  -42.445  7 zeta(3)/32
Fricke N=12 C=9 B=10  (Domb)                   2 3     4.0000              Q   5.980   4.159  -1.821  -2.285  -42.445  7 zeta(3)/24
Fricke N=9 C=27 B=9   (AZ zeta)                2 3     5.1962       Q(sqrt3)   5.980   3.897  -2.083  -2.546  -45.846  L(3,chi_-3)/3
Fricke N=12 C=1 B=34                           2 2    32.0000              Q   4.235   2.079  -2.156  -2.619  -45.049  5 L(2,chi_-3)/16
```

Cross-check against `CDT_FINDER.md` §4A: Zagier C margin **+0.0053** ✓, Zagier E entryC
**−0.0766** ✓, Apéry entryC **−0.435** ✓, Cooper s₁₈ entryR **−1.638** ✓, AZ ζ entryR
**−2.546** ✓. m needed for a positive margin: **m ≥ 13.99** on all four unit hosts (CDT's
14 = 13.9938 rounded up), i.e. exactly CDT's own count and no slack.

**No host produced by this sweep improves on CDT's +0.005, and the only new geometries with
entryR > −0.6 are the λ₂ᶰᵒʳᵐ = 4 tier that Zagier E already occupies.**

---

## 7. The three off-shape unit configurations  **[verified]** (`12_fivept.gp`)

The B-scan of §4.2 produced three admissible normalisations whose λ₁λ₂ = 1 but which the
three-term test rejects. Their exact minimal recurrences (fitted and verified) are:

| N | C | B | minimal recurrence | characteristic roots | #Σ |
|---|---|---|---|---|---|
| 10 | 20 | 9 | 5-term, deg 4 | y⁴−20y³+38y²−20y+1 = (y²−18y+1)(y−1)² → 9±4√5, **1, 1** | 5 |
| 12 | 12 | 7 | 5-term, deg 4 | y⁴−12y³−26y²−12y+1 → 7±4√3, **−1, −1** | 5 |
| 18 | 6 | 5 | 7-term, deg 5 | y⁶−12y⁵+24y⁴−34y³+24y²−12y+1 → 5±2√6, **ζ₆, ζ̄₆ (double)** | 7 |

So each has a *genuine* unit far singularity — at x = 1, x = −1, and at a primitive sixth root
of unity respectively — but it is an **interior order-2 orbifold point of the covering group**,
not one of the two Fricke folds, and it comes with two or four extra punctures. The
descending involution of CDT's architecture has nothing to remove there, and the conditional
function would have to be regular at five or seven points, not four.

**Independent confirmation and a correction to my own numbers.** The N = 10 row
aₙ = 1, 3, 25, 267, 3249, 42795, … is `lattice/sporadic_scan2/table.json` **row #109**
(N = 10, w = 2, t the degree-2 eta quotient (η₂η₅/η₁η₁₀)⁶ — the same coordinate written
differently), where it is recorded with **k = 3** and Apéry limit 3.4696821099…·10⁻⁴, and its
five roots 17.944, 1, 1, 1, 0.0557 agree with the exact factorisation above. The N = 12 and
N = 18 configurations are **not** in that table. *For these ≥5-point rows the pair
(B₀,B₁) = (0,1) does not define the canonical companion (L(y_B) = t needs the full
inhomogeneous solution once the recurrence has order > 2), so the k and ξ my `12_fivept.gp`
prints for them are for a non-canonical solution and are not reported here;* scan 2's k = 3
and limit for the N = 10 row are the authoritative ones. The point that survives is purely
geometric and is exact: **a unit far singularity does occur off the CDT shape, and it costs
two extra punctures and (at least) one extra Eichler integration to reach — 1.74 nats of τ
against 0 nats of ceiling.**

---

## 8. Cusp or fold? the dichotomy, and what it does to the scoring  **[verified]**

| family | t₁ and t₂ are… | local monodromy | example |
|---|---|---|---|
| weight one on Γ₀(N)/Γ₁(N) (§5) | **cusps** | unipotent, log-type | CDT's own host: Σ = {0, 1/9, 1, ∞}, all four cusps of Γ₀(6) |
| Fricke weight two (§4) | **folds** = the two W_N fixed points | ℤ/2, √ -type | Cooper s₇: Σ = {0, 1/27, −1, ∞}, with ∞ at the ℤ/3 elliptic pair of Γ₀(7) |

So **Cooper's s₇ is the only Apéry-perfect k = 2 host of fold type** — the three others (Zagier
A, C, D) are CDT's own log-type geometry. This matters for the transported contour, because
CDT's ψ must avoid the preimages of the extra point y(t₁): at a cusp the conditional function
has a log singularity there (CDT's actual situation, loss 0.62922), at a ℤ/2 fold it has a
square-root singularity and *one* preimage suffices in a different way. The sweep does **not**
design a contour for the fold case; that is flagged as an open input, exactly as
`CDT_FINDER.md` §8 estimate (1) is.

A small positive lemma that does fall out: for any **four-point** host the fixed point x = 2s
of the descending involution w(x) = sx/(x−s) is never itself singular (2/λ₂ = 1/λ₁ would need
λ₁ = λ₂/2, impossible for |λ₁| > |λ₂|), so the quotient orbifold does carry its order-2 cone
point and **the ceiling log(256/λ₂ᶰᵒʳᵐ) is valid on every host of §§4–5**. The degradation to
log(16/λ₂ᶰᵒʳᵐ) reported in `CATALAN_AL_HOSTS.md` §6.2 for the level-16 host is a ≥5-point
phenomenon; it does not touch this census.

---

## 9. The sources, and a clean criterion for the free integration  **[verified]** (`15_sources.gp`, `16_polesite.gp`)

For each of the twelve rows the companion source Φ = F·θ_q x (Theorem A of
`paper/sections/02_sources.tex`) was computed to O(q⁶⁰) and tested against
M₄(Γ₀(N)) and M₄^Eis(Γ₀(N)) with `mfinit`/`mftobasis`. Result:

| row | pole pair of x sits at | Φ holomorphic? | Eisenstein coordinates | free int. | k | period slot |
|---|---|---|---|---|---|---|
| N=6 C=72 B=17 (Apéry) | **cusps** 1/2, 1/3 (u = −1/9, −1/8) | **yes** | (1, −28, 63, −36) | no | 3 | ζ(3) |
| N=8 C=32 B=12 (AZ ε) | **cusps** 1/2, 1/4 (u = −1/8, −1/4) | **yes** | (1, −21, 84, −64) | no | 3 | ζ(3) |
| N=9 C=27 B=9 (AZ ζ) | **cusps** 1/3, 2/3 | **yes** | (0, 0, 0, 1) | no | 3 | L(3,χ₋₃) |
| N=12 C=9 B=10 (Domb) | **cusps** 1/2, 1/6 (u = −1/9, −1) | **yes** | (1, −17, −9, 16, 153, −144) | no | 3 | ζ(3) |
| N=5 C=125 B=22 | ℤ/2 elliptic pair | no (triple poles) | — | yes | 2 | ζ(2) |
| N=6 C=81 B=14 | interior pair, u = (−7±i√…)/81 | no | — | yes | 2 | ζ(2) |
| N=6 C=64 B=20 | interior pair, u = −1/4, −1/16 | no | — | yes | 2 | L(2,χ₋₃) |
| N=7 C=49 B=13 (Cooper s₇) | ℤ/3 elliptic pair | no | — | yes | 2 | ζ(2) |
| N=8 C=16 B=24 | interior pair | no | — | yes | 2 | G |
| N=10 C=25 B=6 (Cooper s₁₀) | interior CM pair u = (−3±4i)/25 | no | — | yes | 2 | ζ(2) |
| N=12 C=1 B=34 | interior pair | no | — | yes | 2 | L(2,χ₋₃) |
| N=18 C=1 B=14 (Cooper s₁₈) | interior pair | no | — | yes | 2 | L(2,χ₋₃) |

**Calibration.** The four Eisenstein coordinate vectors reproduce
Table `tab:sources` of `paper/sections/02_sources.tex` exactly and independently:
γ ↦ 1 − 28·2⁻ˢ + 63·3⁻ˢ − 36·6⁻ˢ, ε ↦ (1−2⁻ˢ)(1−4·2⁻ˢ)(1−16·2⁻ˢ) = 1 − 21·2⁻ˢ + 84·4⁻ˢ − 64·8⁻ˢ,
ζ ↦ the primitive (χ₋₃,χ₋₃) newform of level 9 with coefficient exactly 1,
α ↦ 1 − 17·2⁻ˢ − 9·3⁻ˢ + 16·4⁻ˢ + 153·6⁻ˢ − 144·12⁻ˢ. **[verified]**

**New uniform statement (perfect correlation across all twelve rows).**

> **Criterion.** In the Fricke family, Φ = F·θ_q x is a *holomorphic Eisenstein series*
> **iff** the pole pair of x consists of **cusps**. If instead x has its pole at an
> interior (elliptic / CM) point of the covering curve, Φ is *meromorphic* with a triple
> pole there, the row acquires a **free integration** ((n+1) | Aₙ), the denominator exponent
> drops **k = 3 → k = 2**, and the Apéry limit drops one weight:
> ζ(3) ↦ ζ(2), L(3,χ) ↦ L(2,χ).

This generalises Remark `rem:cooper` (which flagged s₇, s₁₀, s₁₈ as having meromorphic
sources) from Cooper's three rows to *all eight* rows of the family with a free integration,
and it gives the geometric reason for the "weight drop" that `WEIGHT_DROP.md` studies. It is
also exactly why **no Fricke host can produce a positive margin with a weight-three period**:
holomorphic source ⇒ k = 3 ⇒ τ = 5.980 > ceiling 5.545 even at λ₂ᶰᵒʳᵐ = 1; and the free
integration that buys k = 2 buys it by destroying the Eisenstein source and dropping the
period to the ζ(2)/L(2,χ) slot.

---

## 10. **The one genuine hit: a new period on an Apéry-perfect host**
**[verified: 128 digits]** (`18_gamma1_5.gp`, `19_g15_periods.gp`, `20_foldreg.gp`, `21_foldreg2.gp`, `22_newperiod.gp`, `23b_knew.gp`)

The finder's Eisenstein census (`CDT_FINDER.md` §5) is carried out through the two
*orientations* with one character trivial, and concludes that weight three offers
{ζ(2)} ∪ {L(2,χ₋_D)} with χ₋_D **quadratic**. On Γ₀(N) that is right. On **Γ₁(N)** it is not:
the nebentypus may be a higher-order character, and then the inner orientation carries the
L-value of that character. The smallest case is also an Apéry-perfect host.

### 10.1 The host

**Zagier D = Apéry's ζ(2) row on Γ₁(5)**, the fourth of the finder's four λ₂ᶰᵒʳᵐ = 1 rows:
x = q∏(1−qⁿ)^{5(n|5)}, F = Σ Aₙxⁿ with Aₙ = Σ_k C(n,k)²C(n+k,k) = 1, 3, 19, 147, 1251, …,
λ₁ = φ⁵ = (11+5√5)/2, λ₂ = −φ⁻⁵, λ₁λ₂ = **−1** (unit), k = 2,
entryC **+1.310**, entryR **+0.846**, **margin +0.0053** — CDT's own score.
*Calibration:* Φ_D = F·θ_q x is computed here and equals Table `tab:sources` row **D**
coefficient by coefficient, c(m) = Σ_{d|m}(ψ₁(d) − 2ψ₂(d))d² = q − 7q² + 19q³ − 23q⁴ + q⁵ + 47q⁶ − 97q⁷ + …,
and its Apéry limit is ζ(2)/5 to **140 digits**. **[verified]**

### 10.2 The weight-three Eisenstein directions on Γ₁(5)

Mod 5 there are four characters: **1** and χ₅ (even), and the two **odd quartic** characters
ψ₄, ψ̄₄ with ψ₄(2) = i. Weight three needs an odd nebentypus, so the admissible pairs (ψ,φ)
with cond ψ · cond φ | 5 are exactly four, and dim M₃^Eis(Γ₁(5)) = 4 = #cusps:

| direction | (ψ, φ) | L(Φ,s) | critical value at s = 2 | interesting? |
|---|---|---|---|---|
| Φ₁ | (**1**, ψ₄) | ζ(s)·L(ψ₄, s−2) | ζ(2)·L(ψ₄,0) | ζ(2) — known irrational |
| Φ₂ | (**1**, ψ̄₄) | ζ(s)·L(ψ̄₄, s−2) | ζ(2)·L(ψ̄₄,0) | ζ(2) |
| **Φ₃** | (**ψ₄**, **1**) | L(ψ₄,s)·ζ(s−2) | **−½ L(2, ψ₄)** | **yes, and new** |
| **Φ₄** | (**ψ̄₄**, **1**) | L(ψ̄₄,s)·ζ(s−2) | **−½ L(2, ψ̄₄)** | **yes, and new** |

L(2,ψ₄) is a genuinely non-critical value: for odd χ the functional equation pairs s = 2 with
s = −1, where L(−1,χ) = −B_{2,χ}/2 = 0 and Γ((s+1)/2) has a pole, so no algebraic evaluation
follows — L(2,ψ₄) is in the same class as CDT's L(2,χ₋₃) and Catalan's G.

### 10.3 Fold-regularity picks out one ℚ(√5)-rational direction in each pair

Working in the real basis R₁ = Φ₁+Φ₂, R₂ = i(Φ₁−Φ₂), R₃ = Φ₃+Φ₄, R₄ = i(Φ₃−Φ₄), the linear
form Bₙ = [xⁿ](F·θ_q^{−2}Φ) was built for each and the convergence of Bₙ/Aₙ measured at
n = 84…86. Two slow (non-fold-regular) modes appear, with distinct rates 0.98449 and 0.98834;
they are killed by exactly one direction in each pair:

* outer: d(R₁) = −2·d(R₂), so **½R₁ + R₂ = Φ_D** is fold-regular — a re-derivation of Zagier D
  from scratch, with the correct coefficients (½ + i, ½ − i) that the "orientation" language
  predicts; **[verified, residual 4·10⁻¹⁵⁵]**
* inner: d(R₃) = −φ⁵·d(R₄) **exactly** (the ratio is constant in n to 130 digits and equals
  −(11+5√5)/2 = −λ₁), so the fold-regular inner direction is

$$\boxed{\ \Phi_{\rm new}\;=\;R_3+\varphi^5R_4\;=\;(1+i\varphi^5)\,E_3^{\psi_4,\mathbf 1}+(1-i\varphi^5)\,E_3^{\bar\psi_4,\mathbf 1}\ }$$

with Fourier coefficients in **ℤ[φ]** (φ⁵ = 5φ+3 is a unit), so the same integral structure as
the host. The fold-regular subspace of M₃^Eis(Γ₁(5)) is therefore **2-dimensional over ℚ(√5)**.

### 10.4 The new period

$$\boxed{\ \xi_{\rm new}\;=\;\lim_{n\to\infty}\frac{B_n}{A_n}\;=\;\varphi^5\,\mathrm{Im}\,L(2,\psi_4)\;-\;\mathrm{Re}\,L(2,\psi_4)\;=\;0.6556341888406567663309814138723994024111\ldots\ }$$

* **verified to 128 digits** (successive n agree to 128 places; the closed form is reproduced
  by the character theory: the critical value of L(Φ_new, s) at s = 2 is
  −½[(1+iφ⁵)L(2,ψ₄) + (1−iφ⁵)L(2,ψ̄₄)] = −(Re L − φ⁵ Im L), with rational prefactor exactly 1);
* `lindep` on [ξ, 1, ζ(2), Re L, Im L, √5, √5·ReL, √5·ImL] returns
  **[2, 0, 0, 2, −11, 0, 0, −5]**, i.e. 2ξ + 2ReL − (11+5√5)·ImL = 0 and **nothing else** —
  no rational part, no ζ(2) admixture, no elementary π² term;
* **denominator exponent k = 2** for both ℚ-components (2R₃+11R₄ and 5R₄) of Φ_new, measured to
  n = 76 — *the same k = 2 as Zagier D and as CDT's own host*. **[verified]**

### 10.5 Completeness: why Γ₁(5) is the only place this can happen

The inner orientation (ψ, **1**) at weight three exists at level N whenever ψ is an **odd**
primitive character with cond ψ | N; its critical value is −½L(2,ψ), interesting for every
odd ψ. Running this over the complete list of four-point genus-zero hosts (§2, §5):

| four-point host | odd characters available (nebentypus) | interesting inner periods | best λ₂ᶰᵒʳᵐ | entryR |
|---|---|---|---|---|
| Γ₀(5) / **Γ₁(5)** | the two **quartic** ψ₄, ψ̄₄ | **L(2,ψ₄) — new** | **1** (Zagier D) | **+0.846** |
| **Γ₀(6)** | χ₋₃ | **L(2,χ₋₃)** | **1** (Zagier C) | **+0.846** |
| Γ₀(7) | χ₋₇ and the odd sextic characters mod 7 | L(2,χ₋₇), L(2,χ₆ mod 7) | 7 | −1.099 |
| Γ₀(8) | χ₋₄, χ₋₈ | **G**, L(2,χ₋₈) | 4 | −0.540 |
| Γ₀(9) | χ₋₃ and the odd sextic characters mod 9 | L(2,χ₋₃), L(2,χ₉ odd) | 5.196 | −0.801 |

so the levels carrying an interesting inner weight-three period are 5, 6, 7, 8, 9, and of
these **only 5 and 6 are Apéry-perfect**. Level 6 is CDT's theorem. Level 5 is §10.4. Every
other interesting weight-two L-value reachable in this architecture — L(2,χ₋₇), L(2,χ₋₈),
Catalan's G, L(2, sextic mod 7 or 9) — sits on a host with λ₂ᶰᵒʳᵐ ∈ {4, 5.196, 7} and hence
entryR ≤ −0.540, i.e. the Catalan tier or worse. **The pocket is exactly two elements wide,
and one of them is already a theorem.**

Note also that Γ₁(5) is the **unique** four-point genus-zero group whose character group is
not 2-torsion: (ℤ/5)^× ≅ ℤ/4, while (ℤ/6)^× ≅ ℤ/2, (ℤ/8)^× ≅ (ℤ/2)², and at levels 7, 9 the
higher-order characters are available but the geometry is not. That is why the finder's
character-theoretic census, which enumerates the quadratic (χ₋_D) directions, sees everything
except this one place.

### 10.5 What this is and is not

**It is** a (host, weight) pair with λ₂ᶰᵒʳᵐ = 1, k = 2, entryR +0.846 and margin **+0.0053** —
numerically identical to CDT's own — whose period is **not** ζ(2), **not** L(2,χ₋₃) and **not**
Catalan's G, but a ℚ(√5)-combination of the values at s = 2 of the two odd **quartic**
characters mod 5. It is simultaneously target (a) (a single new L-value) and target (b) (the
fold-regular source space on this one host is 2-dimensional and carries the **mixed pair**
{ζ(2)/5, ξ_new}, exactly the configuration CDT exploit at Γ₀(6) with {ζ(2), L(2,χ₋₃)}).
It is the level-5 weight-three analogue of Beukers' weight-four ℚ(√5) object
8ζ(3) − 5√5 L(3,χ₅), but one Eichler integration cheaper: **k = 2, not 3** — which is exactly
the 1.74 nats that `CDT_FINDER.md` §6 identifies as the thing the budget cannot absorb.

**It is not** a theorem, or even a construction. Everything that is missing for
X₁(5) Sym² in `CDT_FINDER.md` §7 (C2) is still missing here, minus item 1 (settled by
`NUMBER_FIELD_HOLONOMY.md`) and minus the weight cost:
1. the conditional function H built from a hypothesised relation
   a + b·ζ(2) + c·ξ_new = 0 over ℚ(√5), with its conjugate, and d_n² arithmetic in ℤ[φ];
2. a pure module of type [1..2n]² on **P**¹∖{0, s, ∞} with s = 1/λ₂ = −φ⁵ — and by
   `INVENTORY_BOUND.md` Cor. 2.2 the number of lcm-free functions is **u₁ = 1**, so the entry
   is the CDT-proportional one, +0.846, and m must be ≥ 13.99 exactly as for CDT;
3. ℚ(√5)(y)-linear independence of the m ≈ 14 functions;
4. a contour at **both** real places (the second place is a tax, `NUMBER_FIELD_HOLONOMY.md`
   §3.1) whose loss is no worse than CDT's 0.62922 — untested, and the margin is +0.005, so
   there is no slack whatsoever;
5. the pure module must be holomorphic at both places simultaneously.

The honest summary is that the sweep found **one** (host, weight, period) triple that is
scored exactly like CDT's theorem and whose period is new, and that everything separating it
from a proof is the same list CDT had to complete for χ₋₃ — over a real quadratic field, with
a fourth-order character, and with zero margin to spare.

### 10.6 Independent confirmation from the M₃/M₄ machine census  **[computed]** (`eis/EIS_REPORT.md`)

The Eisenstein census that `CDT_FINDER.md` §8 lists as not done was run as a separate machine
sweep (full report in `eis/EIS_REPORT.md`; 3996 directions for N ≤ 60 — 1996 at k′ = 3, 2000 at k′ = 4 — with
**0 mismatches** against PARI's `mfdim([N,k],3)` on all 1103 nebentypus components). Its
independent findings, and how they meet §10:

* **Annihilation dimensions.** With "annihilated" = constant term zero at both cusps 0 and ∞
  (computed exactly with `mfslashexpansion`), the census gives, at **N = 5, k′ = 3**:
  dim M₃^Eis(Γ₁(5)) = 4, rank of the two constant-term functionals = 2,
  **dim(annihilated) = 2, periods spanned: L(2,χ₅^{#2}), L(2,χ₅^{#3}) together with the π² class**
  — the two quartic characters. That is exactly §10.3's fold-regular subspace, reached by a
  completely different route (exact constant-term linear algebra rather than the
  convergence-rate analysis of the actual Apéry row). Likewise N = 6, k′ = 3 gives
  dim(annihilated) = 2 with L(2,χ₋₃) and π² — CDT's own configuration.
* **A rigidity theorem for periods.** Over all 550 (weight 3) + 552 (weight 4) nebentypus
  components with N ≤ 60, *the number of distinct interesting periods inside one nebentypus is
  exactly 1, without exception*. So a diamond eigenform can never carry two interesting
  periods; a genuine mixed pair always needs a non-eigenform on Γ₁(N).
* **Weight 3 vs weight 4.** At weight 3 the pair (π², L(2,ψ)) lives *inside one nebentypus*
  ε = ψ whenever cond ψ | N — generic. At weight 4 the pair (ζ(3), L(3,χ_D)) *always* straddles
  ε = 1 and ε = χ_D, so it needs the Γ₁ mechanism (PARI's `mflinear` literally refuses to form
  Beukers' combination, "different characters"). Beukers' identity
  L(F,3) = −½(8ζ(3) − 5√5 L(3,χ₅)) was reproduced to 78 digits as calibration.
* **What the annihilation kills.** An interesting period is annihilated away only at
  (N,k′) = (1,4), (2,4), (3,4), (3,3), (4,3); from N = 4 (k′=4) and N = 5 (k′=3) on, every
  interesting period survives. The oldform shift d > 1 is what saves it — "which is exactly why
  CDT work at level 6 and not at level 3".
* **Smallest levels for the new L-values.** (ζ(2), L(2,χ₋_D)) pairs are available at every
  N divisible by D: D = 7 at N = 7, 14, 21, …; D = 8 at N = 8, 16, …; D = 11 at N = 11, 22, …;
  D = 15 at N = 15, 30, …. (ζ(3), L(3,χ_D)) pairs at every N divisible by D: D = 5 at
  N = 5, 10, 15, …; D = 8 at N = 8, 16, …; D = 12 at N = 12, 24, …; D = 13 at N = 13, 26, ….
  **Every one of these levels other than 5, 6 fails the geometry test of §§2, 5 and 10.5.**

So the two sweeps meet exactly: the ambient Eisenstein space offers a great many interesting
periods (every L(2,ψ) for ψ odd with cond ψ | N, every L(3,ψ) for ψ even), *the geometry is
what is scarce*, and the intersection of "interesting new period" with "λ₂ᶰᵒʳᵐ = 1 and k = 2"
is the single entry of §10.4.

Two further findings of that census that bear directly on §10:

* **No π-contamination at weight four on the small levels.** At k′ = 4 and
  N ∈ {1,…,8, 10, 13} there are **no elementary directions at all**, so the annihilated
  weight-four subspace there carries ζ(3) and L(3,ψ) with no algebraic-multiple-of-π³ term —
  including N = 5 (Beukers') and N = 8.
* **The √5 in Beukers' 8ζ(3) − 5√5 L(3,χ₅) is his geometry, not the Eisenstein bookkeeping.**
  E₄^{1,1}(τ), E₄^{1,1}(5τ) and E₄^{χ₅,1} are all ℚ-rational with critical values −4ζ(3) and
  −½L(3,χ₅), so mixed weight-four pairs are already available over **ℚ** on Γ₁(5). By contrast
  the √5 in §10.4 is *forced*: the fold-regularity ratio d(R₃)/d(R₄) = −φ⁵ is an exact
  irrational, so the new weight-three direction is intrinsically over ℚ(√5) — the same
  phenomenon `CDT_FINDER.md` §5 records for the Galois trace on X₁(5).

**Caveats carried over from the census** (its own §5): it enumerates the *ambient* Eisenstein
space and does not build any source Φ = F·θ_q t, does not treat Atkin–Lehner quotients, assumes
nothing about ℚ̄-linear independence of the periods, and could not run `mfeisenstein` for the
order-10 and order-20 characters at N = 25, k′ = 4, so **15 of the 244 genus-zero directions**
have their cusp-0 constant term extrapolated from the structural rule verified on the other
229 rather than computed (the N = 25 rows of its annihilation table therefore rest on that
extrapolation; the direction census, the PARI dimension cross-check and the mixed-pair
analysis never call `mfeisenstein` and are unaffected). It also did not evaluate the
elementary π^{k′−1} functional numerically, so a cancellation among elementary directions of
different conductors would go undetected — the "+π²/+π³" annotations are weaker claims than
the interesting-period ones.

---

## 11. Coverage: what this sweep does and does not cover  **(honest ledger)**

**Covered exhaustively and exactly.**
1. All W_N-antisymmetric eta-quotient Fricke parameters (N, r, C) with N ≤ 120 and deg u ≤ 4
   — 116 families (`02_fricke.gp`). Enumeration is by cusp-order vector through the Ligozat
   matrix, so it is complete in that box, not a heuristic exponent box.
2. For each, all B ∈ ℤ with λ₂ a unit (closed form, `03_bscan.py`), **and** all B ∈ [1, 90]
   tested for a three-term cubic recurrence (`08_scan3term.gp`) — the CDT-shape test.
3. The exact images of every cusp and every elliptic point of Γ₀(N) under the degree-1
   Fricke Hauptmodul, for all thirteen genus-zero levels (`06_special.gp`, `07_admB.gp`,
   `14_geom.gp`), recognised in closed form; hence the complete admissible-B list.
4. Genus and orbifold data for every Γ₀(N)+W_S with N ≤ 120 generated by prime-power AL
   involutions (`01_groups.gp`).
5. Arithmetic (minimal recurrence, free integration, k, Apéry limit to 200–250 digits) of
   the twelve CDT-shape Fricke rows, and their sources against M₄ (`10_full.gp`,
   `15_sources.gp`).

**Covered but by citation, not recomputed here.**
* The weight-one row search over the whole *integral form space* on genus-zero parameters
  (`SPORADIC_SCAN2.md`: 5.5·10⁵ exact bilinear solves, 3979 parameters, 51 levels N ≤ 104,
  weights 1–3) — it found no new Sym²/Sym³ row, and my §5 geometry is consistent with it.
* `CATALAN_AL_HOSTS.md` for the emptiness of the 4 | N Atkin–Lehner pocket.
* CDT's Bost–Charles integral 11.845 and the contour loss 0.62922, transported unchanged
  (this remains the finder's estimate (1), untested here).

**Not covered — explicit gaps.**
* **deg u ≥ 5 Fricke parameters.** The enumeration box is deg u ≤ 4 (scan 2's box too). A
  Fricke parameter of degree 5+ would be a Hauptmodul of a group between Γ₀(N) and its AL
  quotient of index ≥ 5, which does not exist for AL quotients (index is a power of 2), so
  deg u ∈ {1,2,4,8,…}; deg 8 is untested. I expect nothing there (every extra sheet adds
  singular points), but it is not proved.
* **Levels N > 120**, and Fricke parameters with C ∉ ℤ (C rational non-integral) — excluded
  because x would not have an integral q-expansion, but not proved impossible.
* **Non-eta parametrising forms in the Fricke family.** Theorem 3.4 fixes F = D log u. For
  genus-zero levels with dim M₂(Γ₀(N))^{W_N = −1} > 1 there are other candidate F. Scan 2
  covered the whole integral form space at weights 1–3 for N ≤ 104 and found nothing new;
  I did not redo that.
* **Γ₁(N) and non-trivial nebentypus in the Fricke family.** My Fricke sweep is
  trivial-nebentypus throughout (F = D log u has trivial character), which is exactly why no
  even-character L(3,χ) and no mixed pair can appear in §4. The mixed-pair search lives on
  Γ₁(N) and is the subject of §10.
* **No contour was designed** for any host; the Bost–Charles numerator for the fold-type
  geometry (Cooper s₇) is not computed. Flagged in §8.
* The **three off-shape unit rows** of §7: k and ξ not determined canonically here (§7).

## 12. Scripts and data

| file | what it does |
|---|---|
| `01_groups.gp` | genus-0 Γ₀(N)+W_S for N ≤ 120; μ′, cusp orbits, orbifold defect A → `01_groups.out` |
| `02_fricke.gp` | exhaustive antisymmetric eta-quotient Fricke parameters, deg ≤ 4, N ≤ 120 → `02_fricke.out` (116 families) |
| `02b_check.gp` | Ligozat order/degree/C check on Domb and Apéry |
| `03_bscan.py` | closed-form B-scan for unit λ₂ → `03_units.json` (152 (family,B) pairs) |
| `04a_test.gp` | row-builder calibration on Apéry's ζ(3) row |
| `gen05.py`, `05_rows.gp` | rows for all deg ≤ 2 unit hosts → `05_rows.out` |
| `06_special.gp`, `14_geom.gp` | images of every cusp/elliptic point under the 13 degree-1 Fricke Hauptmoduln |
| `07_admB.gp` | admissible B = −(Cv + 1/v) from the special values → `07_admB.out` |
| `08_scan3term.gp` | machine test for a 3-term cubic recurrence, all 116 families × B ∈ [1,90] → `08_scan.out` (12 hits) |
| `09_fit.gp` | higher-order recurrence fits for the off-shape unit hosts |
| `10_full.gp` | recurrence, free integration, k, ξ (250 digits) for the 12 CDT-shape rows → `10_full.out` |
| `11_score.py` | the finder's τ / ceil / entryC / entryR / margin on every host → `11_scored.json` |
| `12_fivept.gp` | the three off-shape unit configurations |
| `13_lindep.gp` | period identification by ratio and `lindep` |
| `15_sources.gp` | Φ = F θ_q x against M₄ and M₄^Eis (Table tab:sources reproduced) |
| `16_polesite.gp` | where the pole pair of x sits (cusp vs interior) |
| `17_allsub.gp` | genus-0 Γ₀(N)+W_S over **all** AL subgroups, N ≤ 60 → `17_allsub.out` |
| `18_gamma1_5.gp` | builds Zagier D's host on Γ₁(5); peel check against Apéry's ζ(2) numbers |
| `19_g15_periods.gp` | Φ_D reproduced against Table tab:sources row D; the four Eisenstein directions |
| `20_foldreg.gp`, `21_foldreg2.gp` | the fold-regular subspace: d(R₁) = −2d(R₂), d(R₃) = −φ⁵d(R₄) |
| `22_newperiod.gp` | **the new period, 128 digits, and its `lindep` identification** |
| `23b_knew.gp` | denominator exponent k = 2 of the new linear form |
| `eis/` | the M₃/M₄ Eisenstein direction census (§10), scripts and `EIS_REPORT.md` |
