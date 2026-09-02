# Cusp periods of the Apéry-row Eichler integrals: the closed formula, the
# orientation dichotomy, and the monodromy dictionary

*Working directory `lattice/cusp_periods/`. All scripts and raw outputs listed in §11.
Conventions: `D = q d/dq`; `E_k^{ψ,φ} = c₀ + Σ_m (Σ_{e|m} ψ(m/e) φ(e) e^{k-1}) qᵐ`,
`k = r+1`; "outer" = `(1,φ)`, "inner" = `(ψ,1)`. Sources `Φ = Σ_d c_d E_{r+1}^{ψ,φ}(dτ)`
as in `paper/sections/02_sources.tex` Table `tab:sources`.
Tags: **[exact]** = closed form / exact arithmetic; **[verified, N digits]** = high-precision
computation done here; **[estimated]** = fitted or extrapolated. No irrationality claim is made.*

---

## 0. Verdict

| claim | verdict |
|---|---|
| **(P2)/(P3) hold, in a sharper form than stated.** Grouping `D^{-r}E^{ψ,φ} = Σ_f ψ(f) f^{-r} Λ_φ(q^f)` with `Λ_φ(z)=N_φ(z)/(1-z^Q)` gives, at the cusp `a/c` with `ζ=e(a/c)`, an **absolutely convergent** formula `Π_{a/c} = Σ_d c_d d^{-r} Σ_f ψ(f) f^{-r} w_φ(ζ^{df})`, no Abel regularisation needed, `w_φ(z₀)=N_φ(z₀)/2 − z₀N'_φ(z₀)/Q` at `z₀^Q=1` and `N_φ(z₀)/(1−z₀^Q)` otherwise. Each `f`-sum is periodic, so `Π` is a finite `ℚ(ζ_{P'})`-combination of Hurwitz values `P'^{-r}ζ_H(r, j/P')`, i.e. of `L(r,χ)` for `χ mod P'`. | **[exact]** §1 |
| The task's stated mean-zero condition for the outer Abel formula is **the wrong condition**: the pole condition involves `Li_{r+1}`, the Abel simplification `Li_r`. The Hurwitz form `Σ_j φ(j)Li_r(e(αj))(1/2 − j/P)` is valid unconditionally and is what is implemented. | **[exact]** §1.3 |
| **(P1) is an exact identity, not a proportionality:** `ρ = (2πi)^{r+1} a₀(Φ\|_{r+1}γ) / (r!·c^{r+1})` at every cusp `a/c`, `γ=[a,b;c,d]`. Verified against `mfslashexpansion` for rows **A, C, E, γ, ζ** at every cusp of `Γ₀(N_Φ)`. | **[verified, 77 digits]** §2 |
| For the inner orientation the polar coefficient is `ψ(c)c^{-r-1}L(r+1,ψ)`; for `Φ_new` on `Γ₁(5)` at the cusp `1/2` this is exactly the critical-value condition **Re L(3,ψ₄) = φ⁵ Im L(3,ψ₄)**, which is here **proved exactly** from `L(3,χ) = −i(2π)³τ(χ)B_{3,χ̄}/(12N³)`: the ratio is `(2s₁+s₂)/(2s₂−s₁)` with `s₁/s₂ = 2cos(π/5) = φ`, i.e. `φ³·φ² = φ⁵`. | **[exact]**, checked **[verified, 96 digits]** §3 |
| **The orientation dichotomy is a theorem, with a one-line proof.** For the inner orientation `w(z)=z/(1−z)` has `Re w = −1/2` for **every** `z` on the unit circle (and `w(1):=ζ(0)=−1/2`), so `Re Π_{a/c} = ζ(0)·P(r)·L(r,ψ) = L(Φ,r)` at **every** cusp — orientation factor `+1`, Galois equivariance across cusps automatic — and `Im Π_{a/c} = ½ Σ_d c_d d^{-r} Σ_f ψ(f) f^{-r} cot(π a d f / c)`, an elementary `π^r`-multiple. For the outer orientation `Re Λ_φ` is *not* constant on the circle, and the far-fold factor is `[2^{-r}Λ_φ(1)+(1−2^{-r})Λ_φ(−1)]/Λ_φ(1)` = **−11** for Zagier D and **−3** for Zagier A. | **[exact]**, cot-sum **[verified, 6–17 digits]** §4 |
| **Census (12 table sources + 4 `Γ₁(5)` directions, every cusp).** Every near-fold period is reproduced from scratch: `ζ(2)/4, L(2,χ₋₃)/2, L(2,χ₋₃)/2, ζ(2)/5, G/2, 5L(2,χ₋₃)/8, 7ζ(3)/24, ζ(3)/6, 13ζ(3)/54, 7ζ(3)/32, L(3,χ₋₃)/3, L(3,χ₅)/2`. Every cusp value is `(rational)·L(Φ,r) + (rational)·π^r/√D·i`. Only **A (−3)**, **D (−11)** and **ζ (−1/2, −7/2)** have a factor ≠ 1 — exactly the three sources with `φ ≠ 1`. | **[verified, 60 digits]** §5 |
| **(P4) holds, imaginary parts included.** Adaptive Taylor continuation of the inhomogeneous row ODE reproduces the q-side cusp period at **both** finite singular points for A, C, E, F, D. Dictionary: A `0↦1/8, 1/2↦−1`; C `0↦1/9, 1/3↦1`; E `0↦1/8, 1/4↦1/4`; F `0↦1/9, 1/6↦1/8`; D `0↦φ⁻⁵, 1/2↦−φ⁵`. Domb `α`: `1/12↦1/16, 1/4↦1/4`. | **[verified, 100–135 digits]** §6 |
| For `γ` (Apéry's `ζ(3)`), `ε` and `ζ` the far singularity is a **fold**, not a cusp: the monodromy value there (`ζ(3)/6 + 2πi³/15`, `7ζ(3)/32 + iπ³/24`, `−(7/2)·L(3,χ₋₃)/3`) matches **no** cusp of the source level. `ζ`'s far-fold factor **−7/2** is a genuinely new number that the cusp census cannot see. | **[verified, 100 digits]** §6.2 |
| **The coordinator's discrepancy is settled: both statements were right.** `R₁`, `R₂` on `Γ₁(5)` have zero polar part at cusps `0` and `1/2` and genuine cusp periods `6ζ(2)/5`, `−2ζ(2)/5`; but their constant term at **∞** is nonzero (`−4/5`, `2/5`), so `L B = t(Ψ−a₀)/Φ_can` is **not** a rational function of `t` and acquires a singularity at the **near** fold `t₁`. Measured inhomogeneity grows like `φ^{5n}`; `B_n/A_n` converges only logarithmically. The two slow modes are in the exact ratio `−2` (= ratio of the constant terms at `∞`) and cancel in `Φ_D = ½R₁+R₂`, whose inhomogeneity is exactly `[1,0,0,…]`. **So there is a third regularity condition beyond (P1): `a₀(∞) = 0`.** | **[verified, exact rationals to n=300]** §7 |
| **Galois (task 4).** `Φ_new` is the only `K`-rational direction in the census. Being **inner**, its far cusp realises `σ(ξ)` **exactly** (factor `+1`), at every cusp, while the `ζ(2)`-line of the same host carries `−11`. At weight **four** on `Γ₁(5)` all four Eisenstein directions and the whole fold-regular subspace are `ℚ`-rational: **no `K`-rational direction exists**, confirming that Beukers' `√5` in `8ζ(3) − 5√5 L(3,χ₅)` is geometry, not Eisenstein bookkeeping. | **[verified, 60 digits]** §8 |
| **Bonus (Cooper).** By continuation only (the sources are meromorphic): near-fold periods `ζ(2)/7, ζ(2)/5, L(2,χ₋₃)/2` reproduced to 100 digits, and the far-fold factors are **−5 (s₇), −2 (s₁₀), +1 (s₁₈)**. So the two rows with `ψ = 1` behave like *outer* rows and the one with `ψ = χ₋₃` behaves like an *inner* row — the orientation dichotomy survives meromorphy. `s₁₈`'s far-fold imaginary part `π²/(18√3)` **equals Zagier F's** exactly. | **[verified, 100 digits]** §9 |
| **Caveat found:** `Π_{a/c}` depends on the rational `a/c`, **not** on its `Γ`-class. `Π(1/N) = L(Φ,r) ≠ 0 = Π(∞)` although `1/N ∼ ∞`. The Eichler integral is not modular; the dictionary (P4) picks out one representative, and it is the one the analytic continuation actually reaches. | **[verified, 60 digits]** §5.3 |

**One sentence.** *The cusp period of an Apéry row's Eichler integral has an exact finite closed
form as a `ℚ(ζ_P)`-combination of Hurwitz values; its polar part is literally
`(2πi)^{r+1}a₀/(r!c^{r+1})`; whether it is constant across cusps is decided entirely by whether
the outer character is trivial (`Re[z/(1−z)] ≡ −1/2` on the circle), which is why every inner
row has orientation factor `+1` and Galois equivariance for free while Zagier D gets `−11` and
Zagier A `−3`; and every one of these numbers is the monodromy log-coefficient ratio of the row
ODE, verified to 100–135 digits.*

---

## 1. The formula  **[exact]**  (`lib.gp`)

### 1.1 Derivation

With `Φ = Σ_d c_d E_{r+1}^{ψ,φ}(dτ)` and `Θ = D^{-r}Φ = Σ_{m≥1} c(m) m^{-r} qᵐ`, group the
double divisor sum by the `ψ`-index:

$$D^{-r}E^{\psi,\varphi}=\sum_{f\ge1}\frac{\psi(f)}{f^{r}}\,\Lambda_\varphi(q^{f}),\qquad
\Lambda_\varphi(z)=\sum_{e\ge1}\varphi(e)z^{e}=\frac{N_\varphi(z)}{1-z^{Q}},\quad
N_\varphi(z)=\sum_{a=1}^{Q}\varphi(a)z^{a},$$

`Q` = period of `φ` (`φ = 1`: `Q = 1`, `N(z) = z`, `Λ = z/(1−z)`). Put `q = ζe^{-t}`,
`ζ = e(a/c)`, `t → 0⁺`. For `z₀ = ζ^{df}` with `z₀^Q = 1`,

$$\Lambda_\varphi(z_0e^{-u})=\frac{N_\varphi(z_0)}{Q\,u}+\Bigl[\frac{N_\varphi(z_0)}{2}-\frac{z_0N'_\varphi(z_0)}{Q}\Bigr]+O(u),\qquad u=dft,$$

and for `z₀^Q ≠ 1` the value is just `Λ_φ(z₀)`. Since `Σ_f |ψ(f)|f^{-r} < ∞` for `r ≥ 2`,
dominated convergence gives

$$\boxed{\ \Theta=\frac{\rho}{t}+\Pi_{a/c}+O(t),\quad
\rho=\sum_d c_d d^{-r-1}\!\sum_f \frac{\psi(f)}{f^{r+1}}\,n_\varphi(\zeta^{df}),\quad
\Pi_{a/c}=\sum_d c_d d^{-r}\!\sum_f \frac{\psi(f)}{f^{r}}\,w_\varphi(\zeta^{df})\ }$$

with `w_φ(z₀) = N_φ(z₀)/2 − z₀N'_φ(z₀)/Q` and `n_φ(z₀) = N_φ(z₀)/Q` when `z₀^Q = 1`,
`w_φ(z₀) = N_φ(z₀)/(1−z₀^Q)` and `n_φ(z₀) = 0` otherwise. For `φ = 1` this is `w(1) = −1/2`
— which is exactly `ζ(0)`, the source of the "uniform weight `−1/2`".

### 1.2 The explicit finite formula (P3)  **[exact]**

`w_φ(ζ^{df})` depends on `f` only mod `c_d := c/gcd(c,d)`, and `ψ(f)` only mod `Q_ψ`. With
`M = lcm(c_d, Q_ψ)`,

$$\sum_{f\ge1}\frac{\psi(f)w_\varphi(\zeta^{df})}{f^{s}}
=M^{-s}\sum_{j=1}^{M}\psi(j)\,w_\varphi\!\bigl(e(adj/c)\bigr)\,\zeta_H(s,j/M),$$

exact and finite. Since `ζ_H(r, j/M)` is a `ℚ`-combination of `L(r,χ)`, `χ mod M`, times
`M^r/φ(M)` and Gauss sums, **every cusp period is a `ℚ(ζ_{M})`-linear combination of the
`L(r,χ)` with `χ mod M`, `M = lcm(c, cond φ, cond ψ)`** — which is (P3). Implemented in
`lib.gp:blk/cperiod`.

### 1.3 Where the task statement needed correcting  **[exact]**

The Mellin form `Θ = (1/2πi)∫Γ(s)A(s+r)t^{-s}ds`, `A(s) = Σ c(m)e(am/c)m^{-s}`, shows
`Π_{a/c} = A(r)` and `ρ = Res_{s=r+1}A(s)`. In the outer orientation
`A(s) = P^{r-s}Σ_j φ(j)Li_s(e(αj))ζ_H(s−r, j/P)`, so

* the pole condition is `Σ_j φ(j)Li_{r+1}(e(αj)) = 0` — with `Li_{r+1}`, **not** `Li_r`;
* the value is `Π = Σ_{j=1}^{P} φ(j)Li_r(e(αj))(1/2 − j/P)` **unconditionally**, which reduces
  to the task's Abel formula `−(1/P)Σ j φ(j)Li_r(e(αj))` only when `Σ_j φ(j)Li_r(e(αj)) = 0`.

For Zagier D at both cusps `0` and `1/2` both conditions happen to hold (`Σ_{even}w = Σ_{odd}w = 0`),
which is why the two versions agree there.

### 1.4 Validation  **[verified]**  (`01_check.gp`)

| check | result |
|---|---|
| `Π(Φ_D, 0) = ζ(2)/5`, `Π(Φ_D, 1/2) = −11ζ(2)/5` | agreement to **60 digits** |
| `Π(Φ_new, ·) = ξ`, `Π(Φ'_new, ·) = ξ'` at all four cusps | agreement to **60 digits** |
| Abel summation `Σ_{m≤60000} c(m)m^{-2}ζ^m e^{-mt}`, `t = 0.0025, 0.00125`, one Richardson step | `0.32898681336964528729448303332920347…` vs `ζ(2)/5`, and `−3.61885494706609816023931336662123733…` vs `−11ζ(2)/5`: **32 digits** |

---

## 2. (P1): the polar part is the constant term, exactly  **[verified, 77 digits]**  (`04_polar.gp`)

Approaching `a/c` with `γ = [a,b;c,d]`, `Φ(γτ') = (cτ'+d)^{r+1}(Φ|γ)(τ')` and `τ = a/c+iy`
give `cτ'+d = i/(cy)`; `r` integrations of `(2πi)^{-r}d/dτ` then produce

$$\boxed{\ \rho \;=\; \frac{(2\pi i)^{\,r+1}\,a_0(\Phi|_{r+1}\gamma)}{r!\;c^{\,r+1}}\ }$$

Checked against `mfslashexpansion` (PARI convention: our `E_k^{ψ,φ} = mfeisenstein(k, φ, ψ)`)
at **every** cusp of `Γ₀(N_Φ)` for **A, C, E, γ, ζ**: differences `≤ 5.6·10⁻⁷⁷`. Examples:

| row | cusp | `a₀` | `ρ` computed | `(2πi)^{r+1}a₀/(r!c^{r+1})` |
|---|---|---|---|---|
| A | 1/3 | `−1/8` | `0.5741903088944411143606725…i` | same, diff `1.6·10⁻⁷⁷` |
| C | 1/2 | `−0.06415002990995841828…i` | `−0.9945267882188398388359…` | same, diff `1.8·10⁻⁷⁷` |
| E | 1/2 | `−i/16` | `−0.9689461462593693804836…` | same, diff `2.7·10⁻⁷⁷` |
| γ | 1/2, 1/3 | `−1/9`, `1/4` | `−1.80387205618…`, `0.80172091386…` | same, diff `≤5.4·10⁻⁷⁷` |
| ζ | 1/3, 2/3 | `±0.192450089729875…i` | `±0.617165047242777…i` | same, diff `≤5.6·10⁻⁷⁷` |

In the **inner** orientation the general formula collapses to
`ρ = Σ_d c_d d^{-r-1} ψ(c_d) c_d^{-r-1} L(r+1,ψ)` — automatically zero when
`gcd(c, cond ψ) > 1`, and otherwise a **critical `L(r+1)`-value condition**, exactly as (P1)
predicts. (`Γ₁(5)`: PARI's `mflinear` refuses `Φ_D` because its two components have different
nebentypus; the `Γ₁(5)` constant-term table is instead cross-checked against
`lattice/gamma15/task2/TASK2.md` §4(d), which it reproduces cusp-for-cusp.)

---

## 3. The critical-`L`-value condition on `Γ₁(5)`, proved exactly  **[exact]**  (`03_L3exact.gp`)

`Φ_new = (1+iφ⁵)E₃^{ψ₄,1} + (1−iφ⁵)E₃^{ψ̄₄,1}` is inner, so its polar coefficient at the cusp
`1/c` is `c^{-3}[(1+iφ⁵)ψ₄(c)L(3,ψ₄) + c.c.]`. At `c = 1` fold-regularity is exactly

$$(1+i\varphi^{5})L(3,\psi_4)+(1-i\varphi^{5})L(3,\bar\psi_4)=0
\iff \operatorname{Re}L(3,\psi_4)=\varphi^{5}\operatorname{Im}L(3,\psi_4).$$

**Proof.** From `Σ_{m≥1} sin(2πmx)/m³ = (2π)³B₃(x)/12` and finite Fourier inversion, for
primitive odd `χ mod N`

$$L(3,\chi)=-\,\frac{i\,(2\pi)^{3}\,\tau(\chi)\,B_{3,\bar\chi}}{12\,N^{3}},\qquad
B_{3,\bar\chi}=N^{2}\sum_{b}\bar\chi(b)B_3(b/N).$$

For `ψ₄ mod 5`: `B_{3,ψ̄₄} = (12−6i)/5`, `τ(ψ₄) = −2s₂+2is₁` with `s₁ = sin(2π/5)`,
`s₂ = sin(4π/5)`. Multiplying out,

$$\operatorname{Re}L(3,\psi_4)=\tfrac{16\pi^{3}}{7500}(12s_1+6s_2),\qquad
\operatorname{Im}L(3,\psi_4)=\tfrac{16\pi^{3}}{7500}(12s_2-6s_1),$$

so the ratio is `(2s₁+s₂)/(2s₂−s₁)`. Since `s₁/s₂ = 2cos(π/5) = φ`, this is
`(2φ+1)/(2−φ) = φ³/φ⁻² = φ⁵`. **∎**  (`φ³ = 2φ+1`, `φ⁻² = 2−φ`.)

Numerically: predicted `L(3,ψ₄)` agrees with `lfun` to `5·10⁻⁹⁷`, and
`Re/Im − φ⁵ = −7.5·10⁻⁹⁶`. **[verified, 96 digits]**

**A trap recorded.** The naive generalized-Bernoulli ratio `B_{3,Reψ₄}/B_{3,Imψ₄} = (12/5)/(6/5) = 2`
is **not** `φ⁵`: the Gauss sum does not commute with `Re`/`Im`, so the `L`-ratio is not the
`B`-ratio. Only the full formula above gives `φ⁵`.

At the cusp `1/2` the same expression is `(1/8)[i(1+iφ⁵)L(3,ψ₄) − i(1−iφ⁵)L(3,ψ̄₄)] = −2.76207971863902708052…`,
which is exactly the `ρ` produced by the general machinery — an independent confirmation, and
the `log²` of `gamma15` §4(b).

---

## 4. The orientation dichotomy — a one-line proof  **[exact]**  (`10_cot.gp`)

**Inner (`φ = 1`).** `w(z) = z/(1−z) = 1/(z^{-1}−1)`. For `z = e^{iθ}`,

$$\operatorname{Re}w=\frac{\cos\theta-1}{2-2\cos\theta}=-\tfrac12\ \ \text{for every }\theta\ne0,
\qquad \operatorname{Im}w=\tfrac12\cot(\theta/2),$$

and the regularised value at `z = 1` is `ζ(0) = −1/2` too. Hence at **every** cusp

$$\operatorname{Re}\Pi_{a/c}=-\tfrac12\,L(r,\psi)\!\sum_d c_d d^{-r}=\zeta(0)P(r)L(r,\psi)=L(\Phi,r),$$
$$\operatorname{Im}\Pi_{a/c}=\tfrac12\sum_d c_d d^{-r}\sum_{f}\frac{\psi(f)}{f^{r}}\cot\frac{\pi a d f}{c}.$$

**So for the inner orientation the orientation factor is `+1` at every cusp, and Galois
equivariance across cusps is a triviality.** The cotangent sum was checked directly
(truncated at `f = 4·10⁵`): C at 1/3 `1.4·10⁻⁶`, E at 1/4 `6.2·10⁻⁷`, F at 1/6 `1.2·10⁻⁷`
(second-order rows converge like `1/f²`), and ε at 1/4 `3.9·10⁻¹⁸`, α at 1/3 `1.0·10⁻¹⁷`
(third-order rows, `1/f³`). **[verified, 6–17 digits]**

**Outer (`ψ = 1`).** `Re Λ_φ` is not constant on the unit circle, and at `c = 2`

$$\frac{\Pi_{1/2}}{\Pi_{0}}=\frac{2^{-r}\Lambda_\varphi(1)+(1-2^{-r})\Lambda_\varphi(-1)}{\Lambda_\varphi(1)} .$$

* **Zagier D** (`φ = (1,−2,2,−1,0)`, `Q = 5`, `r = 2`): `Λ(1) = 1/5`, `Λ(−1) = −3`, factor
  `(1/20 − 9/4)/(1/5) = **−11**`.
* **Zagier A** (`φ = χ₋₃`, `Q = 3`, `r = 2`, oldform `1 − V₂`): `Λ(1) = 1/3`, `Λ(−1) = −1`,
  factor `**−3**`.

Both reproduced exactly by the machinery. This is the general form of `gamma15` §3.5.

---

## 5. The census  **[verified, 60 digits]**  (`02_census.gp`, `07_ident.gp`, `11_g15.gp`)

### 5.1 Near-fold periods, re-derived from the cusp formula alone

| row | `r` | `N_Φ` | orientation | `L(Φ,r) = Π(cusp 0)` |
|---|---|---|---|---|
| **A** | 2 | 6 | outer `χ₋₃`, `(1−V₂)` | `ζ(2)/4` |
| **B** | 2 | 36 | inner `χ₋₃` | `L(2,χ₋₃)/2` |
| **C** | 2 | 6 | inner `χ₋₃` | `L(2,χ₋₃)/2` |
| **D** | 2 | 5 | outer `re ψ₄ − 2 im ψ₄` | `ζ(2)/5` |
| **E** | 2 | 8 | inner `χ₋₄` | `G/2` |
| **F** | 2 | 12 | inner `χ₋₃` | `5L(2,χ₋₃)/8` |
| **α** | 3 | 12 | trivial/trivial | `7ζ(3)/24` |
| **γ** | 3 | 6 | trivial/trivial | `ζ(3)/6` |
| **δ** | 3 | 12 | trivial/trivial | `13ζ(3)/54` |
| **ε** | 3 | 8 | trivial/trivial | `7ζ(3)/32` |
| **ζ** | 3 | 9 | `χ₋₃` **both** | `L(3,χ₋₃)/3` |
| **η** | 3 | 20 | inner `χ₅` | `L(3,χ₅)/2` |

Every entry agrees with the ledger's Apéry-limit column. This is an independent
re-derivation of Theorem B's value column from the cusp side alone.

### 5.2 All cusps: fold-regularity, factor, imaginary part

Notation: **fr** = `ρ = 0`; factor = `Re Π / L(Φ,r)`; `Im` in units of `π^r`.

| row | cusp | fr | factor | `Im Π` | | row | cusp | fr | factor | `Im Π` |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 0 | ✓ | 1 | 0 | | α | 0 | ✓ | 1 | 0 |
| A | 1/2 | ✓ | **−3** | 0 | | α | 1/2 | ✗ | 1 | 0 |
| A | 1/3 | ✗ | −1/3 | 0 | | α | 1/3 | ✓ | 1 | `π³/36` |
| A | 1/6 | ✓ | 1 | 0 | | α | 1/4 | ✓ | 1 | `π³/48` |
| C | 0 | ✓ | 1 | 0 | | α | 1/6 | ✗ | 1 | `π³/108` |
| C | 1/2 | ✗ | 1 | 0 | | α | 1/12 | ✓ | 1 | 0 |
| C | 1/3 | ✓ | 1 | `2π²/(9√3)` | | γ | 0 | ✓ | 1 | 0 |
| C | 1/6 | ✓ | 1 | 0 | | γ | 1/2 | ✗ | 1 | 0 |
| E | 0 | ✓ | 1 | 0 | | γ | 1/3 | ✗ | 1 | `π³/27` |
| E | 1/2 | ✗ | 1 | 0 | | γ | 1/6 | ✓ | 1 | 0 |
| E | 1/4 | ✓ | 1 | `π²/16` | | ε | 0 | ✓ | 1 | 0 |
| E | 1/8 | ✓ | 1 | 0 | | ε | 1/2 | ✗ | 1 | 0 |
| F | 0 | ✓ | 1 | 0 | | ε | 1/4 | ✗ | 1 | `π³/64` |
| F | 1/2 | ✗ | 1 | 0 | | ε | 1/8 | ✓ | 1 | 0 |
| F | 1/3 | ✓ | 1 | `π²/(6√3)` | | **ζ** | 0 | ✓ | 1 | 0 |
| F | 1/4 | ✓ | 1 | `π²/(12√3)` | | **ζ** | 1/3 | ✗ | **−1/2** | 0 |
| F | 1/6 | ✓ | 1 | `π²/(18√3)` | | **ζ** | 2/3 | ✗ | **−1/2** | 0 |
| F | 1/12 | ✓ | 1 | 0 | | **ζ** | 1/9 | ✓ | 1 | 0 |

`B` (level 36, 11 cusps), `δ` (level 12) and `η` (level 20) behave identically on the real
part: factor `1` at **every** cusp. Their imaginary parts are
`B`: `{4/27, −4/27, 1/12, 2/27, −2/27, 2/27, 1/27, 2/27, 1/9}·π²/√3`;
`δ`: `{2/81, 7/432, 1/81, 5/324}·π³`;  `η`: `{3/80, 4/125, 2/125, 9/500}·π³/√5`
— full lists in `07_ident.out`.

**On `Γ₁(5)`** (cusps `∞, 0, 1/2, 2/5`, `x = 0, φ⁻⁵, −φ⁵, ∞`):

| direction | ∞ | 0 | 1/2 | 2/5 |
|---|---|---|---|---|
| `Φ_D` (outer) | `0` (fr ✓) | `ζ(2)/5` (fr ✓) | `−11ζ(2)/5` (fr ✓) | `−7ζ(2)/25` (fr ✗) |
| `Φ_new` (inner) | `0` | `ξ` (fr ✓) | `ξ` (fr ✗) | `ξ + 7.0335782…i` (fr ✓) |
| `Φ'_new` (inner) | `0` | `ξ'` (fr ✗) | `ξ'` (fr ✓) | `ξ' + 0.3171086…i` (fr ✓) |
| `R₁ = 2E₃^{1,re ψ₄}` | `0` | `6ζ(2)/5` (fr ✓) | `−6ζ(2)/5` (fr ✓) | `−0.8685251…` (fr ✗) |
| `R₂ = −2E₃^{1,im ψ₄}` | `0` | `−2ζ(2)/5` (fr ✓) | `−8ζ(2)/5` (fr ✓) | `−0.0263189…` (fr ✗) |

`½R₁+R₂ = Φ_D` checks: `½(6/5) − 2/5 = 1/5`, `½(−6/5) − 8/5 = −11/5`. ✓

### 5.3 The representative caveat  **[verified, 60 digits]**

`Π` is a function of the rational number `a/c`, not of its `Γ`-class: `1/N ∼ ∞` under
`Γ₀(N)` (and `1/5 ∼ ∞` under `Γ₁(5)`), yet `Π(1/N) = L(Φ,r)` while `Π(∞) = 0`. The Eichler
integral is not modular and its constant term at a cusp shifts by the period polynomial. §6
shows which representative the geometry actually selects.

---

## 6. (P4): the monodromy dictionary  **[verified, 100–135 digits]**  (`mono.gp`, `05_mono.gp`, `06_mono3.gp`)

### 6.1 Method

The row operator in `x` (second order, `P(x) = 1 − ax + cx²`):
`x P y'' + (P + xP')y' + (cx − b)y = R`, `R = 0` for `A`, `R = 1` for the companion `B`.
Third order (`R3` rows, `P₃ = 1 − 2ax + cx²`, `d`-term included):
`x²P₃ y''' + 3x(1−3ax+2cx²)y'' + (1−(6a+2b)x+(7c+d)x²)y' + (−b+(c+d)x)y = R`.
`mono.gp` shifts every coefficient polynomial to the current point and generates 260 Taylor
coefficients per step, stepping `0.35 ×` the distance to the nearest singularity, at 100–120
digit working precision, starting from 520–620 exact rational coefficients at `x = 0`.
The period is `ℓ(B)/ℓ(A) = ΔB/ΔA` with `Δ = (M−1)`, computed independently from `y`, `y'`
(and `y''`).

### 6.2 Results

| row | singular point | monodromy period | q-side cusp | `|diff|` |
|---|---|---|---|---|
| A | `1/8` (near) | `ζ(2)/4` | `0` | `3.5·10⁻¹³⁵` |
| A | `−1` (far) | `−3ζ(2)/4` | `1/2` | `6.0·10⁻¹³⁴` |
| C | `1/9` (near) | `L(2,χ₋₃)/2` | `0` | `1.6·10⁻¹³⁵` |
| C | `1` (far) | `L(2,χ₋₃)/2 + (2π²/(9√3))i` | `1/3` | `1.1·10⁻¹¹⁹` |
| E | `1/8` (near) | `G/2` | `0` | `3.3·10⁻¹³⁵` |
| E | `1/4` (far) | `G/2 + (π²/16)i` | `1/4` | `6.7·10⁻¹²⁰` |
| F | `1/9` (near) | `5L(2,χ₋₃)/8` | `0` | `4.1·10⁻¹²⁰` |
| F | `1/8` (far) | `5L(2,χ₋₃)/8 + (π²/(18√3))i` | **`1/6`** | `4.0·10⁻¹²⁰` |
| D | `φ⁻⁵` (near) | `ζ(2)/5` | `0` | `1.6·10⁻¹³⁵` |
| D | `−φ⁵` (far) | `−11ζ(2)/5` | `1/2` | `2.9·10⁻¹²¹` |
| α | `1/16` (near) | `7ζ(3)/24` | `0` (and `1/12`) | `<10⁻¹⁰⁰` |
| α | `1/4` (far) | `7ζ(3)/24 + (π³/48)i` | `1/4` | `<10⁻¹⁰⁰` |

So (P4) is confirmed, **including the imaginary parts**, which are the elementary
`π^r/√D`-multiples of §4 and not artefacts.

**Where it fails, and why — folds.** For `γ` (Apéry's `ζ(3)`), `ε` and `ζ` the far singularity
is a `W_N`-fixed point (order-2 orbifold point), not a cusp, and the value matches **no** cusp:

| row | far singularity | monodromy value | nearest cusp value | verdict |
|---|---|---|---|---|
| γ | `17+12√2` | `ζ(3)/6 + (2π³/15)i` | cusp `1/3` gives `ζ(3)/6 + (π³/27)i` | **fold** (ratio `18/5`) |
| ε | `(3+2√2)/4` | `7ζ(3)/32 + (π³/24)i` | cusp `1/4` gives `+ (π³/64)i` | **fold** (ratio `8/3`) |
| ζ | `−0.7182335127930838…` | `−(7/2)·L(3,χ₋₃)/3` | cusps give factors `1` and `−1/2` | **fold**, new factor `−7/2` |

The near-fold values still equal `L(Φ,r) = Π(cusp 0)` in all three cases (Theorem B), so the
*near* fold and the cusp `0` agree numerically even when the geometry differs. That coincidence
is the content of Theorem B's Fricke half, not of (P4).

---

## 7. The `R₁`/`R₂` discrepancy: both statements were right  **[verified, exact rationals, n ≤ 300]**  (`00_disc.gp`, `00b_rhs.gp`)

Building `x(q)`, `F(q)` from the recurrence `(11,3,−1)` alone (Frobenius → nome → reversion, no
modular input) and peeling `B_n = [xⁿ](F·D^{-2}Ψ)` exactly:

| direction | `B₁…B₅` | `B_n/A_n` at `n = 298` | q-side `Π(0)` | defect |
|---|---|---|---|---|
| `Φ_D` | `1, 25/4, 1741/36, …` | `0.328986813369645287294483033329205037843789980241…` | `ζ(2)/5` | `0` to **96 digits** |
| `R₁` | `2, 33/2, 2605/18, …` | `1.66888438460542632…` | `6ζ(2)/5` | `−0.30503649561244540099804782494548498855912…` |
| `R₂` | `0, −2, −24, …` | `−0.50545537893306787…` | `−2ζ(2)/5` | `+0.15251824780622270049902391247274249427956…` |

The two defects are in the **exact** ratio `−2` at every `n` (80 digits), which is the ratio of
the constant terms at `∞`: `c₀(R₁) = −4/5`, `c₀(R₂) = 2/5` (and `½(−4/5) + 2/5 = 0`).

The inhomogeneities `R_n := (n+1)²B_{n+1} − (11n²+11n+3)B_n − n²B_{n−1}` are

| direction | `R₀ … R₅` | ratio `R_{n+1}/R_n` at `n = 10` |
|---|---|---|
| `Φ_D` | `1, 0, 0, 0, 0, 0` | — (exact `[1,0,0,…]`) |
| `R₁` | `2, 16, 156, 1564, 15992, 165780` | `10.70…` → `φ⁵` |
| `R₂` | `0, −8, −78, −782, −7996, −82890` | same (`R₂ = −R₁/2 + δ_{n0}`) |
| `R₃` (inner) | `2, 22, 244, 2706, 30010, 332816` | `11.0901699437494742410…= φ⁵` |
| `R₄` (inner) | `0, −2, −22, −244, −2706, −30010` | `R₄,ₙ = −R₃,ₙ₋₁` |

so **the inhomogeneity of `R₁`, `R₂` is singular at the near fold `t₁ = φ⁻⁵`**. It is also *not*
a rational function of `x`: multiplying by `(1−11x−x²)` or its square does not truncate.

**Why.** `L B_Ψ = C F(Ψ − a₀) = t(Ψ − a₀)/Φ_can` with `Φ_can = F·θ_q t` of weight `r+1`.
`tΨ/Φ_can` is a modular function, hence rational in `t`; but `t·a₀/Φ_can` has weight `−(r+1)`
and is **not**. And `Φ_can` vanishes at the cusp over the near fold (`Φ_D` vanishes at `0`),
which is where the singularity appears. Hence:

> **A third regularity condition.** (P1) fold-regularity (`ρ = 0` at the fold cusp) is
> **necessary but not sufficient** for the companion to have an Apéry limit there. One also
> needs `a₀(Φ, ∞) = 0`, i.e. `Φ` must vanish at the cusp over `x = 0`. On `Γ₁(5)` the outer
> weight-three space has this only on the line `⟨Φ_D⟩` — which is exactly `hostscan` §10.3's
> observation, arrived at from the other side.

At non-fold-regular cusps the q-side constant term is still well defined (it is `A(r)`), but it
is **not** the monodromy period; the coordinator's warning is correct and was respected
throughout §6.

---

## 8. Galois (task 4)  **[verified, 60 digits]**  (`08_galois.gp`)

**Weight three, `Γ₁(5)`, the inner `K`-rational line.** `Φ_new = (1+iφ⁵)E₃^{ψ₄,1} + c.c.`
has coefficients in `ℤ[φ]` and `σ(Φ_new) = Φ'_new`. Being inner, §4 applies verbatim:

$$\Pi_{\mathfrak c}(\Phi_{\rm new})=\xi,\qquad \Pi_{\mathfrak c}(\Phi'_{\rm new})=\xi'=\sigma(\xi)
\quad\text{at \emph{every} cusp } \mathfrak c\in\{0,\,1/2,\,2/5,\,1/5\},$$

to `10⁻⁷⁷`. **The far cusp realises the Galois conjugate exactly, orientation factor `+1`.**
By contrast the `ζ(2)`-line of the same host (Zagier D, outer) has factor `−11`. So the
`gamma15` §3.5 conclusion — equivariance is a theorem on the `ξ`-line and fails only on the
`ζ(2)`-line — is not special to `Γ₁(5)`: it is the inner/outer dichotomy of §4.

**Weight four, `Γ₁(5)` (Beukers' setting).** The four Eisenstein directions
`E₄^{1,1}(τ), E₄^{1,1}(5τ), E₄^{χ₅,1}, E₄^{1,χ₅}` have cusp periods, at all four cusps,

| direction | `Re Π` | `c₀(∞)` |
|---|---|---|
| `E₄^{1,1}(τ)` | `−ζ(3)/2` | `1/240` |
| `E₄^{1,1}(5τ)` | `−ζ(3)/250` | `1/240` |
| `E₄^{χ₅,1}` | `−L(3,χ₅)/2` | `0` |
| `E₄^{1,χ₅}` | `0` (since `L(0,χ₅) = 0` for even `χ`) | `−B_{4,χ₅}/8 = 1` |

All four are `ℚ`-rational, all polar coefficients are `ℚ`-rational, and so is the fold-regular
subspace: **no `K`-rational direction is forced at weight four on `Γ₁(5)`.** This independently
confirms `hostscan` §10.6 — the `√5` in Beukers' `8ζ(3) − 5√5 L(3,χ₅)` comes from his geometry,
not from the Eisenstein bookkeeping — and it isolates weight three as the place where `φ⁵` is
forced (by the fold-regularity ratio `d(R₃)/d(R₄) = −φ⁵`, equivalently §3's `L(3)`-identity).

**No other `K`-rational direction exists in the census**: rows A–F and α–η all have sources
defined over `ℚ` (rows B, C, E, F, ζ, η over `ℚ` with quadratic nebentypus; D over `ℚ` after
combining `ψ₄`, `ψ̄₄`).

---

## 9. Bonus: Cooper's meromorphic rows  **[verified, 100 digits]**  (`09_cooper.gp`)

No q-side formula exists (the sources have double poles at CM points,
`COMPANION_ARITHMETIC.md` §4.5), so this is continuation only.

| row | `(a,b,c,d)` | near fold | period | far fold | period | factor |
|---|---|---|---|---|---|---|
| `s₇` | `(13,4,−27,3)` | `1/27` | `ζ(2)/7` | `−1` | `−5ζ(2)/7` | **−5** |
| `s₁₀` | `(6,2,−64,4)` | `1/16` | `ζ(2)/5` | `−1/4` | `−2ζ(2)/5` | **−2** |
| `s₁₈` | `(14,6,192,−12)` | `1/16` | `L(2,χ₋₃)/2` | `1/12` | `L(2,χ₋₃)/2 + (π²/(18√3))i` | **+1** |

The near-fold values reproduce the ledger exactly (100 digits). The pattern matches the
characters `ψ_{s₇} = ψ_{s₁₀} = 1`, `ψ_{s₁₈} = χ₋₃` determined in `COMPANION_ARITHMETIC.md` §4.5:
**the two rows whose `ψ` is trivial behave like *outer* rows (factor ≠ 1) and the one with
`ψ = χ₋₃` behaves like an *inner* row (factor `+1`)** — the orientation dichotomy survives
meromorphy. `s₁₈`'s far-fold imaginary part `π²/(18√3)` is **numerically identical** to Zagier
F's, to all 100 digits, extending the Conjecture-D family `ξ^B = ξ^C = ξ^{s₁₈}` to the far fold.

---

## 10. Honest ledger — what is not done, and what is uncertain

1. **`ρ` vs constant term on `Γ₁(5)`** was not run through `mfslashexpansion` because PARI's
   `mflinear` refuses to combine forms of different nebentypus. It is cross-checked against
   `gamma15/task2` §4(d) instead (which used `mfslashexpansion` componentwise) and agrees
   cusp-for-cusp, but the *identity* `ρ = (2πi)^{r+1}a₀/(r!c^{r+1})` was verified numerically
   only for the five `Γ₀(N)` sources of §2.
2. **The imaginary parts of `Π(Φ_new)`, `Π(Φ'_new)` at the width-5 cusps `2/5`, `1/5` of
   `Γ₁(5)`** (`7.0335782469643621718…`, `0.3171086774431510891…`, `1.5855433872157554459…`)
   were **not identified**: `lindep` against `π²·{1, 5^{1/4}, √5, 5^{3/4}}` at 40-digit
   tolerance finds nothing. They presumably live in `ℚ(ζ₅,φ)·π²`. They are not needed — those
   cusps map to `x = ∞` and `x = 0`. **[estimated / open]**
3. **Rows `B`, `δ`, `η` have complex-conjugate finite singularities**, so no real far-fold
   monodromy run was done for them; their cusp periods are computed and identified, but the
   cusp ↔ `x` dictionary for them is **not** verified by continuation. **[not done]**
4. **The cotangent-sum identity of §4** was checked by truncated summation only
   (`f ≤ 4·10⁵`), giving 6–7 digits at `r = 2` and 17 digits at `r = 3`. The proof is exact;
   the check is weak at `r = 2` purely because `Σ f^{-2}cot` converges slowly.
5. **The "`Π` depends on the representative"** phenomenon (§5.3) is recorded but not turned
   into a period-polynomial statement. The natural conjecture — that
   `Π(γ(a/c)) − Π(a/c)` is the period polynomial of `Φ` evaluated at the cusp — is untested.
6. **The `−7/2` at row `ζ`'s far fold** and the fold values `2π³/15` (γ), `π³/24` (ε) have no
   q-side derivation here: fold values are not cusp values, and the fold analogue of §1 (a
   `W_N`-twisted Mellin split, as in `thmB_exact/03_fricke_fold.gp`) was **not** developed.
7. No irrationality or independence claim is made anywhere. All `lindep`/`bestappr`
   identifications are at 40–42 digit tolerance against explicitly named bases, and are
   reported as identifications, not proofs.
8. The task's suggestion that the outer polar condition is "mean zero of `φ(j)Li_r(ζ^j)`" is
   **wrong** (§1.3); the implementation uses the unconditional Hurwitz form, so nothing
   downstream depends on the erroneous version.

---

## 11. Scripts and outputs

| script | what it does | output |
|---|---|---|
| `lib.gp` | the twelve table sources + four `Γ₁(5)` directions; `cperiod` (`Π`, `ρ`) by the §1 closed form; `cusplist`; direct Abel summation | — |
| `mono.gp` | generic adaptive Taylor continuation of `Σ_j P_j(x)y^{(j)} = R` | — |
| `00_disc.gp` | `x, F` from the recurrence `(11,3,−1)`; exact `B_n` for `R₁, R₂, Φ_D, Φ_new` to `n = 300`; convergence | `00_disc.out` |
| `00b_rhs.gp` | the inhomogeneities `R_n` of the four `Γ₁(5)` directions | `00b_rhs.out` |
| `01_check.gp` | validation on `Γ₁(5)`; Abel + Richardson cross-check | `01_check.out` |
| `02_census.gp` | `Π`, `ρ`, fold-regularity, `Π/Π(0)` at every cusp of `Γ₀(N_Φ)` for all sources | `02_census.out` |
| `03_L3exact.gp` | the exact `Re L(3,ψ₄) = φ⁵ Im L(3,ψ₄)` proof and its numerical check | `03_L3exact.out` |
| `04_polar.gp` | (P1) `ρ = (2πi)^{r+1}a₀/(r!c^{r+1})` vs `mfslashexpansion` | `04_polar.out` |
| `05_mono.gp` | (P4) for the second-order rows A, C, E, F, D at both folds | `05_mono.out` |
| `06_mono3.gp` | (P4) for γ, α, ε | `06_mono3.out` |
| `07_ident.gp` | exact identification of every cusp period and of `L(Φ,r)` | `07_ident.out` |
| `08_galois.gp` | task 4: the `K`-rational line at weight three; the weight-four census on `Γ₁(5)` | `08_galois.out` |
| `09_cooper.gp` | Cooper `s₇, s₁₀, s₁₈` and row `ζ`'s far fold, by continuation | `09_cooper.out` |
| `10_cot.gp` | the two structural theorems of §4 (cotangent sum; `Λ_φ` on the circle) | `10_cot.out` |
| `11_g15.gp` | the `Γ₁(5)` census rows with identification | `11_g15.out` |
| `12_imid.gp` | attempted identification of the width-5-cusp imaginary parts (negative) | `12_imid.out` |

The `.out` files are present in this directory but are **not** tracked by git
(`*.out` is in the repository `.gitignore`, as for `lattice/hostscan/`); regenerate any of
them with `gp -q <script> > <script-basename>.out`.

Run with `gp -q <file>` from this directory (`05`, `06`, `09` want a large stack; they set
`parisizemax` themselves).
