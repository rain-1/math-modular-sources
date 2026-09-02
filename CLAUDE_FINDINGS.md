# Findings log (Claude's working notes)

Recovered results, cross-checks, and other findings from mining the archive,
kept separate from Fable's ledger work on `packages/Odd_Zeta_Research_Master_Ledger.xlsx`.
Newest entries on top.

---

## 2026-09-02 — Magnetic sources on Apéry's host: exist, keep k = 2, and are killed by a unit

`consolidation/MAGNETIC_APERY.md` + `lattice/magnetic_apery/`. Complete census of weight-4 meromorphic forms Φ = F²ρ(u) on Γ₀(6): three magnetic families. Φ₋₈ (double pole at the disc −8 CM point, Fricke sign −1, the canonical source of the (C,B) = (81,14) host) transplanted to Apéry's host keeps the free integration (k = 2) and gives B_n/A_n → ζ(2)/8 (185 digits). Theorem U: 1/x = W + B with W a Hauptmodul with integral expansion, so the pole coordinate satisfies |N(W₀+B)| ≥ 1 and the companion has a singularity at |x| ≥ 1 (here x(τ₀) = 1 exactly): A_nξ − B_n ≍ n⁻², no irrationality. A magnetic source only helps its native host (Cooper's poles sit at x = ∞). By-products: free integration ⟺ meromorphy on all 12 Fricke hosts; L(D⁻¹(Δ/E₄²), 2) = 0.0610392510075379919… apparently a new constant.

## 2026-09-02 — Imaginary quadratic fields: the holonomy bound is free, but no host exists

`consolidation/IMAG_QUADRATIC.md` + `lattice/imag_quadratic/`. Over K imaginary quadratic the number-field bound is letter-for-letter the rational one (one complex place, d_v = 2 = n_K; no second place, no hypothesis tax). But a non-real period needs a nebentypus of order ≥ 3, and over all subgroups H ≤ (Z/N)^× with N ≤ 60 the only four-point genus-zero host carrying one is X₁(5) with ψ₄, whose fold-regular ψ₄-component is zero (constant-term functionals have rank 2 there because cond ψ₄ = N). Supply starts at level 7 with ≥ 6 special points; K-rational hosts with conjugate singularities are dead (|λ₁| = |λ₂|) or not three-term. Best margin −18. The would-be prize (L(2,ψ₄) ∉ Q(i)) would have scored exactly CDT's numbers.

## 2026-09-02 — Cooper's magnetic congruence: two cells proved, the rest reduced, and a master conjecture

`consolidation/COOPER_CONGRUENCE.md` + `lattice/cooper_congruence/`.

- **Reduction (proved).** eq:magnetic (Φ|U_{p^n} ≡ (ψ(p)p)^n Φ mod p^{n+2} for all n) ⟺ Ξ|U_p ≡ ψ(p)Ξ mod p² for Ξ = D⁻¹Φ ∈ Z[[q]] ⟺ rad(n)² | β(n), β = c' ⋆ μψ.
- **Theorem.** Φ_{s₇}|U₇ = 7Φ_{s₇} and Φ_{s₁₀}|U₅ = 5Φ_{s₁₀} identically (Fricke sign −1, Atkin–Lehner trace to level N/p, valence formula, two exact coefficients). These are the cells with ε_p = −ψ(p).
- **Mechanism located.** c'(m) = Res_x(η q(x)^{−m}) with η = l(x)dx/(x√P F) = Ξ dq/q; (S) mod p ⟺ Cartier eigenvector C(η) ≡ ψ(p)η ⟺ a_{pj+p−1} ≡ ψ(p)a_j mod p (verified p ≤ 101); a_{p−1} ≡ ψ(p) manufactures χ₋₃ for s₁₈ on the x-side. For ψ = 1 it is a mirror-map-type integrality exp∫η ∈ Z_p[[x]]. Paşol–Zudilin's level-one magnetic forms FAIL (S): the congruence is not generic magnetism.
- **Master conjecture (verified n ≤ 1500, sharp):** n² | β(n); equivalently Φ = D³[Σ_e ψ(e)e⁻²K(q^e)] with K = Σγ(d)q^d integral, i.e. Ξ is an inner weight-3 Eisenstein series whose Lambert kernel q/(1−q) is replaced by an integral non-multiplicative K, and L(Φ,s) = L(s−1,ψ)·Σγ(d)d^{3−s}. Implies (S) and supercongruences such as c'(p²) ≡ ψ(p)c'(p) mod p⁴. In Shimura–Borcherds terms: a(|D|m²) = m³γ(m).

## 2026-09-02 — Cusp periods (one recurrence, many periods) and cusp-form sources

`consolidation/CUSP_PERIODS.md` + `lattice/cusp_periods/`; `consolidation/CUSPFORM_SOURCES.md` + `lattice/cuspform_sources/`.

- **Cusp-period formula.** The Eichler integral of an Eisenstein source splits by orientation, D^{-r}E^{ψ,φ} = Σ_f ψ(f)f^{-r}Λ_φ(q^f) with Λ_φ = N_φ/(1−z^Q) rational; at a cusp a/c the expansion is ρ/t + Π + o(1) with absolutely convergent finite formulas, every cusp period a Q(ζ_M)-combination of L(r,χ), χ mod M = lcm(c, conductors). Validated to 60 digits; the monodromy period at fold-regular cusps equals it (100–135 digits, six rows).
- **Polar part = constant term, exactly:** ρ = (2πi)^{r+1}a₀(Φ|γ)/(r! c^{r+1}) (77 digits). Inner orientation: ρ ∝ ψ(c)L(r+1,ψ), a critical-value condition; Track D's identity Re L(3,ψ₄) = φ⁵ Im L(3,ψ₄) is exactly this condition on Γ₁(5), now proved from generalised Bernoulli numbers.
- **Third regularity condition:** a₀(Φ,∞) = 0 is also necessary (otherwise the inhomogeneity is not rational and is singular at the near cusp); fold-regular space = {a₀(∞)=0} ∩ {a₀(near)=0}. Resolves the hostscan §10.3 slow modes.
- **Orientation dichotomy (one line):** Re[z/(1−z)] = −1/2 on the whole circle, so inner rows have factor +1 at every cusp and are Galois-equivariant across cusps (Φ_new realises ξ and its conjugate σ(ξ), at every cusp of Γ₁(5), 77 digits); outer rows pick up [2^{-r}Λ(1)+(1−2^{-r})Λ(−1)]/Λ(1): −11 (Zagier D), −3 (Zagier A). Cooper's meromorphic rows obey the same dichotomy (−5, −2, +1 matching ψ = 1, 1, χ₋₃).
- **Cusp-form sources on Fricke hosts (r=3):** fold-regular iff Φ|W_N = −Φ, then ξ = L(Φ,3) (period-polynomial conditions automatic). No weight-4 newform of level 5–9 has sign −1, so nothing on the Apéry-perfect hosts; five oldform combinations f − d²f(dτ) give new Apéry-like rows in 1 and L(f,3) for f = 5.4.a.a, 6.4.a.a, 9.4.a.a (e.g. Domb numbers with ξ = L(6.4.a.a,3)/2, rate 1/4, K₋ = L(6.4.a.a,2)/√π, the central value). No irrationality (−4.4 nats).

## 2026-09-02 — Bucket 4 closed: inventory, host sweep, two-variable systems, and the Γ₁(5) theorem attempt

Four tracks (`consolidation/GAMMA15_CLOSURE.md` §6 for the table; `INVENTORY_BOUND.md` §4, `lattice/hostscan/REPORT.md`, `TWO_VARIABLE_HOLONOMY.md` §6, `lattice/gamma15/REPORT.md`).

- **Inventory.** Pólya–Carlson + θ-trick: lcm-free admissible functions are polynomials, so u₁ = 1 on every Apéry-perfect host. CDT's fourteen functions are the complete supply at max eᵢ ≤ 1; the single-layer module is the arcsine tower Σ yⁿ/(nʲC(2n,n)). Entry passes on Catalan and X₁(5), but the best margin over all inventories is −8.04 / −6.99.
- **Hosts.** No new Apéry-perfect host for N ≤ 120. One new period on Zagier D's host Γ₁(5): ξ = φ⁵ Im L(2,ψ₄) − Re L(2,ψ₄) = 0.6556341888…, ψ₄ the odd quartic character mod 5, k = 2, source in Z[φ][[q]].
- **Two variables.** Split-denominator companions do not exist (tail obstruction; Landau ceiling ℓ₁+ℓ₂ ≤ log 256|u₀|).
- **Γ₁(5) theorem (1, ζ(2), ξ independent over Q(√5)).** Does not follow. Same host, k, fourteen functions and τ = 16603/3920 as CDT, entry +0.67, but the averaged bound gives m ≤ 16.15 against 14 (margin −1.44), and that is on a *doubled* hypothesis: the far-cusp period of Zagier D's companion is −11ζ(2)/5, so the second real place needs a relation the first does not imply (an orientation effect: Li₂ parity at q = −1; on the ξ-line equivariance is exact). On the honest hypothesis the margin is −4.40. Second real place responsible for both.
- **By-products.** Two corrections to CDT's Appendix A (published slit parameters miss a preimage pair at |z| = 0.7539; BC of the published φ is 12.12 not 11.845; their constants and theorem are nevertheless right); a Jensen reduction of the Bost–Charles integral (excess over log|φ′(0)| is a pure multivalency term); the identity Re L(3,ψ₄) = φ⁵ Im L(3,ψ₄); a monodromy criterion showing the decisive bad preimages can never be dropped (fixed lines at distinct cusps differ).

## 2026-09-02 — Two censuses: asymptotic constants, and the arithmetic of the companion sequences

**Asymptotic constants** (`consolidation/ASYMPTOTIC_CONSTANTS.md`, scripts `lattice/asymptotic_constants/`).
All seven weight-two rows with a Fricke fold (Domb α, Apéry γ, ε, ζ, Cooper s₇, s₁₀, s₁₈) have
`K = √N/(2π^{3/2}) · √(λ₁/(λ₁−λ₂))`, λ₁,₂ the characteristic roots — proved from the structure
`x = u/(1+Bu+Cu²)`, `F = D log u`, `u|W_N = 1/(Cu)` (verified exactly to O(q^300)). Apéry's
`(1+√2)²/(2^{9/4}π^{3/2})` confirmed to 73 digits; Cooper's rows give `3√3/(4π^{3/2})`, `√2/π^{3/2}`,
`3√2/π^{3/2}`. Above weight two a Chowla–Selberg period Ω^{k−2} (or Ω^k if the companion vanishes at
the fold) survives: level-16 ζ(5) has `K ∈ Q̄·Γ(1/4)⁴π^{−9/2}`, level-12 ζ(7) has
`K ∈ Q̄·Γ(1/3)^{12}π^{−19/2}`, both in closed form. Correction: the level-12 ζ(7) host
`η₄²η₁₂²/(η₁²η₃²)` has degree 2 on X₀(12); the Hauptmodul is `v = 1/h₁₂` (shared with the level-12
ζ(5) system), `x = v(1+3v)/(1+4v)` — which explains the missing order-≤20 recurrence.

**Companion arithmetic** (`consolidation/COMPANION_ARITHMETIC.md`, scripts `lattice/companion_arithmetic/`).
Companion formula `b_n = Σ c(m) m^{−r} e_{n,m}` verified exactly for all twelve Eisenstein rows.
Refined denominator `R_n` is exact at every bad prime but only an upper bound at good primes
(mechanism: `v_p(den b_n) = max(0, r − v_p(c(p) e_{n,p}))` for n/2 < p ≤ n). `v_p(b_n) − v_p(a_n) = v_p(ξ_p)`
constant wherever ξ_p ≠ 0. New exact digit laws `v₂(a_n) = s₂(n)` for s₁₀ and s₁₈ (identical
valuation vectors), with companion laws `v₂(b_n) = s₂(n) − 2⌊log₂n⌋ − 1`. A Lucas law for companions at
every good prime, `β_{np+m} ≡ ψ(p) β_n a_m (mod p)` with `β_n = p^{r⌊log_p n⌋} b_n`, holding in all 199
(row, prime) cells, which determines ψ for Cooper's meromorphic sources: ψ(s₇) = ψ(s₁₀) = 1,
ψ(s₁₈) = χ₋₃ (previously conjectural).

---

## 2026-08-21 — Recovered: full leading asymptotics of the ζ(5) linear form (level-12 CM-isogeny construction)

**Status:** proved exact, up to one constant left numerical (not yet closed form).
**Source:** `packages/phase 1 zeta math package/modular_apery_cm_isogeny_research_notes.txt`, lines 527–594 (sections K "Convergence theorem" and L "Full leading asymptotics"). `C_5 = (11/144) ζ(5)` is defined at line 157.
**Ledger status:** absent from `Odd_Zeta_Research_Master_Ledger.xlsx` — checked all 7 sheets for `cm_isogeny`, `C_5`, `level-12`, `saddle`, `x_+`; no match. This is the "asymptotic formula for the linear form size" that `chats/tasks.txt` flagged as lost.

### The result

Let `W = V - C_5 U`. Fricke invariance makes `W` regular at the nearest branch point
`x_+ = 7 - 4√3`. The only finite singularities of the differential system are
`-1, 7-4√3, 7+4√3`, so `W` is holomorphic in `|x| < 1`, while `U` has a square-root
singularity at `x_+`:

```
c_n ~ K x_+^{-n} n^{-3/2}       (K a nonzero constant)
```

hence `d_n/c_n -> C_5`, with root bound `limsup |d_n/c_n - C_5|^{1/n} <= 7 - 4√3`.

At `x = -1`, the regularized solution has

```
W(x) = (13/480) log^5(1+x) + O(log^4(1+x))
```

so

```
d_n - C_5 c_n  ~  (13/96) (-1)^{n+1} (log n)^4 / n
```

Combined with `c_n ~ -K (7+4√3)^n n^{-3/2}`, this gives the full asymptotic of the
linear form `d_n - C_5 c_n` (linear form in `1` and `ζ(5)`):

```
d_n/c_n - C_5  ~  (13/(96K)) (-1)^n (7-4√3)^n n^{1/2} (log n)^4
```

Numerically, `K ≈ 0.6853718484905344059...`, so `13/(96K) ≈ 0.19758130854792016...`.

**Open piece — CLOSED 2026-09-02** (see `consolidation/ZETA5_K_CLOSED_FORM.md`): `K = Γ(1/3)^12 π^{-17/2} · 3^{8/3}(15√3−23)√((12+7√3)/2) / (2^{37/3}(2205−1273√3)(45+26√3)^{1/3})`, verified to 51 digits, via `K² = 49 x₊ (DE)(τ*)² / (2π Q(x₊)² (−D²x)(τ*))` with `E(τ*) = 0` at the Fricke point.

**Possible next steps (not yet started):**
- ~~Derive the closed form for `K` from the discriminant `-48` CM period.~~ Done 2026-09-02.
- Add this as a formal row to the xlsx ledger (Fable's file) once wording is agreed.
