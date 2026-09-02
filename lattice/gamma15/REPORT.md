# 1, ζ(2), ξ over ℚ(√5): the Calegari–Dimitrov–Tang architecture on Γ₁(5),
# run at **both** real places — with a correction to CDT's own Appendix A

*Working directory `lattice/gamma15/`. All scripts and raw outputs listed in §8.
Conventions: `CDT_FINDER.md` §§1–3 (τ, ceiling, entry, margin), `NUMBER_FIELD_HOLONOMY.md`
(the number-field bound is the **average** over archimedean places), `INVENTORY_BOUND.md`
(u₁ = 1 forced; CDT's fourteen functions are the complete admissible supply at max eᵢ ≤ 1),
`lattice/hostscan/REPORT.md` §10 (the host, the source Φ_new and the period ξ).
Tags: **[exact]** = closed form / exact rational arithmetic; **[verified, N digits]** =
high-precision computation in this task; **[cited]**; **[open]**.
No margin is ever rounded in the favourable direction.*

---

## 0. Verdict

| claim | verdict |
|---|---|
| **The target.** ξ = φ⁵·Im L(2,ψ₄) − Re L(2,ψ₄) = 0.65563418884065676633098141387239940241…, the Apéry limit of the fold-regular inner Eisenstein direction Φ_new on Zagier's row **D** / Γ₁(5); re-derived here from scratch, 142 digits of internal agreement, and \|ξ − (φ⁵ Im L − Re L)\| < 10⁻¹⁴². `lindep` returns **2ξ + 2 Re L − (11+5√5) Im L = 0** and nothing else. | **[verified, 142 digits]** §2 |
| **The number-field bound is the average.** CDT's own Remark `BCboundK` (quoted verbatim in §1.1) gives, for K = ℚ(√5) with two real places, m ≤ (BC₁+BC₂)/(L₁+L₂−2τ) = **B̄C/(L̄−τ)**. This is exactly `NUMBER_FIELD_HOLONOMY.md` Thm 2.2, and it is CDT's, not ours. | **[cited, exact]** §1.1 |
| **τ is unchanged.** The host is Apéry-perfect (λ₂ = −φ⁻⁵ a unit) with k = 2, the same as CDT's; the descent array is CDT's (m = 14, r = 2, b₁ = b₂ = 2, u₁ = 1, u₂ = 3, Σeᵢ = 6, max eᵢ = 1), so **τ = τ♭+τ♯ = 191/49 + 27/80 = 16603/3920 = 4.2354592**. | **[exact]** §1.2, §3 |
| **Correction to CDT Appendix A (i).** Their published slit/lune parameters do **not** satisfy their Lemma A.4.4. The point −1/72 has **six** (not four) non-principal h-preimages inside \|z\| < 77/100 — at \|z\| = 0.401921, 0.592651, **0.753925** (each a conjugate pair); the argument principle gives 7 preimages in that disc. Four slits cannot exclude six points, and with their θᵢ five of the six lie **inside** their Ω (winding number +1). | **[verified]** §1.3 |
| **Correction to CDT Appendix A (ii).** The Bost–Charles integral of *their published* φ is **12.1207**, not 11.845 (two independent evaluators, N up to 16384). | **[verified, 5 digits]** §1.3 |
| **CDT's numbers are nevertheless exactly right.** A systematic *successive-slitting* designer, run with CDT's own R = 77/100, c = 75/10 and shrink 995/1000 and six slits, gives \|ψ′(0)\| = **0.6293401** (CDT 0.6292232680) and BC = **11.8449** (CDT 11.845), bound **13.9906** (CDT 13.9938). Only the transcription of the slit *placement* parameters is wrong; the contour and both of its constants are right. Optimising R in the same family beats them slightly: \|ψ′(0)\| = 0.646136, BC = 12.2007, bound **13.9761**, margin **+0.0209**. | **[verified]** §1.4 |
| **A new exact reduction of the Bost–Charles integral.** For φ = h∘Φ with Φ univalent, BC(φ) = log\|φ′(0)\| + ∫_𝕋 S(z) dμ, S(z) = Σ_{w∈𝔻, φ(w)=φ(z)} log(1/\|w\|). **BC − log\|φ′(0)\| is a pure multivalency term**; it vanishes iff φ is univalent (recovering CDT's univalent case BC = log ρ). The interior solutions are the Γ₀(2)-orbit points of Φ(z) lying in Ω. | **[proved]** §1.3 |
| **The two places are geometrically very different.** t₁t₂ = −1 with t₁ = φ⁻⁵, t₂ = −φ⁵; s := t₂ is forced. At v₁ the extra point is Y₁ = y(t₁)/s = **−φ⁻¹⁵/(5√5) = −6.5574·10⁻⁵**, hyperbolically deep, and the nearest bad h-preimage sits at \|z\| = 0.5360 — **better than CDT's 0.4019**. At v₂ it is Y₂ = **−φ¹⁵/(5√5) = −122.000**, i.e. \|Y₂\|/4\|s\| = 30.5, and the nearest bad preimage sits at \|z\| = **0.2137** — far worse. | **[exact / verified]** §4 |
| **The contours.** Certified (all bad preimages verifiably excluded, principal preimage inside): v₁ \|ψ′(0)\| = **0.7037841**, L₁ = 7.5999530, BC₁ = 15.223765; v₂ \|ψ′(0)\| = **0.3949390**, L₂ = 2.2100944, BC₂ = 6.408897. | **[verified]** §4 |
| **Entry passes.** L̄ = **4.9050237** > σ_m = 4 (Remark BCboundK's holonomicity condition) and L̄ − τ = **+0.6695645** > 0. | **[verified]** §5 |
| **The margin is NEGATIVE.** B̄C = 10.816331, bound **m ≤ 16.1543** against the exhibited m = 14: **margin = 14(L̄−τ) − B̄C = −1.4424**. Transporting CDT's best convexity improvement (a factor 0.97336 on the bound) still leaves **m ≤ 15.72**, margin −1.15. | **[verified]** §5 |
| **Place v₂ is entirely responsible.** With the per-place score s_v = 13 log\|ψ_v′(0)\| − shape_v (margin = 13 log 256 + mean(s_v) − 14τ): s₁ = **−12.1905**, s₂ = **−16.2761**, CDT's own = −12.7700. So v₁ *gains* 0.579 nats over CDT and v₂ *loses* 3.506; halved, that is **−1.463** nats of margin against CDT's +0.021. | **[verified]** §5 |
| **No inventory rescues it.** Adding functions to CDT's fourteen forces max eᵢ ≥ 2, and τ jumps 4.2355 → 4.4482; the bound goes from 16.13 to 23.5 and the margin from −1.44 to −3.96. Every m ≠ 14 is worse. (`INVENTORY_BOUND.md` R1/R4 already proves 14 is the complete supply at max eᵢ ≤ 1.) | **[verified]** §6 |
| **The second place also taxes the HYPOTHESIS, not just the geometry.** The far-cusp periods are computed exactly: the homogeneous-solution obstruction at t₂ gives π(B_D) = **−(11/5)ζ(2)** and π(B′_new) = **ξ′ = −Re L − φ⁻⁵ Im L** (both to 12 digits; the same construction at the near cusp returns ζ(2)/5 and ξ, calibrating it). So regularity of σH at the v₂-fold t₂ is **(ii)** σ(a) − (11σ(b)/5)ζ(2) + σ(c)ξ′ = 0, which is σ applied to (i) *except* in the ζ(2) slot (+1/5 ↦ −11/5): **(i) ⇏ (ii)**. One K-relation gives a conditional function at one place only. The *sources* are nevertheless perfectly Galois-equivariant (L B_new = 2/(1−x/t₂), L B′_new = 2/(1−x/t₁); Φ_new kills the cusp over t₁, Φ′_new the cusp over t₂). | **[verified, 12 digits]** §3.4 |
| **The −11 is explained, and equivariance holds on the ξ-line.** It is the *orientation*, not Galois: the outer (ζ(2)) direction is a sum of dilogarithms, and Li₂ evaluates to ζ(2) at q → +1 but to −ζ(2)/2 on odd d at q → −1, turning the character sum Σ j·w(j) = −1 (factor 1/5) into Σ j·(wg)(j) = 22 (factor −11/5). The inner (ξ) direction is a Lambert series whose cusp weight is −1/2 *uniformly* at both cusps, so there π^{(t₂)}(σΦ) = σ(π^{(t₁)}(Φ)) exactly. | **[proved]** §3.5 |
| **By-product: a new weight-three identity.** The convergence condition behind that computation is **Re L(3,ψ₄) = φ⁵·Im L(3,ψ₄)** — verified to 212 digits, `lindep` = [−2, 11, 5] — the s = 3 shadow of the weight-two fold-regularity ratio −φ⁵. | **[verified, 212 digits]** §3.5 |
| **Both branches are negative.** Assuming the *doubled* hypothesis (i)∧(ii): margin **−1.4424**. Keeping the single, correct hypothesis and letting Ω₂ avoid *every* preimage of Y₂ (principal included, at \|z\| = 0.074191, cap 0.2572): margin **−4.4038**. | **[verified]** §5.2 |
| **The verdict does not rest on any unproved exclusion.** Enforcing *only* the two exclusions proved necessary in §6(b) and letting every other bad preimage be kept or cut, whichever is better, the optimum still cuts them all (they sit near the boundary, cost almost nothing, and lower BC): relaxed margin **−1.4577**, against −1.4416 fully constrained. | **[verified]** §6(b′) |
| **The one loophole, and its closure.** We exclude *every* non-principal h-preimage of Y_v, as CDT do; only the ones at which the continued conditional function is really singular need go. If the single conjugate pair at \|z\| = 0.213693 at v₂ were harmless the margin would be **+0.51**. It is not harmless: those two preimages are the images of z₀ under the parabolic generator (2,±1) of Γ₀(2) at the cusp y = ∞, whose loop on the x-line encircles the outer cusp t₂, so the continuation changes H by a nonzero multiple of u_{t₂}, the homogeneous solution holomorphic at t₂; and **u_{t₁} and u_{t₂} are distinct lines** (numerically \|u_{t₁} ∧ u_{t₂}\| = 0.318 on Γ₁(5), 0.737 on CDT's Γ₀(6); both monodromies unipotent to 10⁻¹³). So the continued H is singular at the fold: the pair must be excluded. | **[verified]** §6(b) |

**One sentence.** *Everything arithmetic about the target is as good as CDT's own — the same
Apéry-perfect host, the same k = 2, the same fourteen functions, the same τ = 16603/3920, a
genuinely new period ξ, Galois-equivariant sources, and (by CDT's own Remark BCboundK) a
legitimate two-place bound whose entry test passes with 0.67 nats to spare — but the second real
place of ℚ(√5) charges twice: it puts the removed fold at 30 orbifold radii instead of CDT's
1/72 and its deepest bad preimage at \|z\| = 0.214 instead of 0.402, which costs 3.5 nats there
and 1.46 on the average and turns CDT's +0.02 into **−1.44**; and its fold-regularity condition
is not the Galois conjugate of the hypothesis (the far-cusp period of B_D is −11ζ(2)/5, not
ζ(2)/5), so a single K-relation buys a conditional function at one place only, and the variant
that survives on the honest hypothesis is worse still, at **−4.40**.*

---

## 1. Calibration [Task 1]

### 1.1 CDT's number-field theorem, quoted exactly  **[cited]**

`papers/cdt/cdt2/L2chi.tex`, Remark `BCboundK` (line 6849 ff.), verbatim:

> The proof immediately gives the following formal generalization to a number field $K$.
> For each $\sigma: K\hookrightarrow \C$, we consider a holomorphic mapping
> $\varphi_\sigma: (\Db, 0) \to (\C,0)$ with $\varphi'_\sigma(0)\neq 0$.
> Assume there exists an $m$-tuple $f_1, \ldots, f_m \in K \llbracket x \rrbracket$ of
> $K(x)$-linearly independent formal functions with denominator types of the form
> $$f_i(x) = a_{i,0} + \sum_{n=1}^{\infty} a_{i,n} \frac{x^n}{n^{e_i} [ 1, \ldots, b_{i,1} \cdot n] \cdots [1,\ldots, b_{i,r} \cdot n]}, \qquad a_{i,n} \in \OL_K,$$
> where $e_i, b_{i,j}$ are the same as in Theorem~\ref{main:elementary form},
> and such that for all $i \in \{1,\ldots,m\}$ and $\sigma:K\hookrightarrow \C$, we have
> $f_i(\varphi_\sigma(z)) \in \C \llbracket z \rrbracket$ convergent on $|z| < 1$. If
> $$\frac{1}{[K:\Q]}\sum_{\sigma:K\hookrightarrow \C} \log |\varphi'_\sigma(0)| > \sigma_m ,$$
> then all $f_i$ are holonomic functions, and
> $$m  \leq  \frac{\sum_{\sigma:K\hookrightarrow \C}  \iint_{\T^2} \log|\varphi_{\sigma}(z)-\varphi_{\sigma}(w)| \, \mv(z) \mv(w) }{ (\sum_{\sigma:K\hookrightarrow \C}  \log{|\varphi_{\sigma}'(0)|}) - [K:\Q](\ovtau(\bb)+ \tau^\sharp(\be)) }.$$
> The convexity improvements also extend in the obvious way.

For K = ℚ(√5) (two real places, d_v = 1, [K:ℚ] = 2) this reads

$$\boxed{\;m\;\le\;\frac{\mathrm{BC}_1+\mathrm{BC}_2}{L_1+L_2-2\tau}\;=\;\frac{\overline{\mathrm{BC}}}{\bar L-\tau},\qquad L_v=\log|\varphi_v'(0)|,\;\;\tau=\tau^\flat(\mathbf b)+\tau^\sharp(\mathbf e),}$$

i.e. exactly the averaged form of `NUMBER_FIELD_HOLONOMY.md` Thm 2.2. The load-bearing input
that `CDT_FINDER.md` §7 (C2) listed as unproved is therefore **not** unproved: it is CDT's own
remark. Note also that CDT's *holonomicity* hypothesis is stated with σ_m = 4, not with τ.

### 1.2 τ reproduced  **[exact]** (`cdt_bound.py`)

| quantity | value | CDT |
|---|---|---|
| σ_m | 4 | 4 |
| τ♭(**b**), u₁ = 1, u₂ = 3, b₁ = b₂ = 2 | **191/49** = 3.8979592 | 191/49 |
| I₂¹⁴(2) | **21.075** exactly | 21.075 |
| τ♯(**e**), **e** = (0,0,1,0⁶,1,1,1,1,1) | **27/80** = 0.3375 at ξ ≈ 2.163 | 27/80, minimiser [2, 13/6] |
| τ(**b**;**e**) | **16603/3920** = 4.2354592 | 16603/3920 |
| CDT's bound 11.845/(log(256·0.6292232680) − τ) | **13.9938** | 13.9938 |

### 1.3 The contour machinery, and two corrections to CDT's Appendix A

`contour.py` implements from scratch: the Γ₀(2) Hauptmodul h(z) = −256 z ∏(1+zⁿ)²⁴ =
−256 Δ(2τ)/Δ(τ) = λ + λ/(λ−1) (CDT eq. `defofh`) and h′; the enumeration of h-preimages by the
Γ₀(2)-orbit (bottom rows (2c,d), gcd = 1, mod ±); CDT's slit map `Slit(z,r)` (their eq. `slit`)
**and its exact inverse**, derived from their own chain
w ↦ v = (1+w)²/w ↦ u = i√((v−A)/(v−4)) ↦ (u−i)/(u+i), A = −(1−r)²/r; and their lune map
`lune(z,c)` (eq. `firstgobble`) with the inverse they print.

`region2.py` builds Φ = (z ↦ Rz)∘T₁∘…∘T_k with T_j(z) = e^{iβ_j}·(Slit or lune)(z). A slit
removes, **from the current region**, the Φ_{j−1}-image of the radial segment at angle β_j+π
from radius r_j out to 1. Given a target point p, one takes ζ = Φ_{j−1}^{−1}(p), β_j = arg ζ − π,
r_j = |ζ|·(over ≤ 1); then p lies **in** the removed segment, and since Ω_k ⊂ Ω_j for k > j,

> **exclusion of every bad point holds by construction, not by a numerical test.**

Conformal radius is exact: R·∏(4r_j/(1+r_j)²)·∏((c²−1)/(c²+1)). A shrink z ↦ ρz (CDT's own
device, ρ = 995/1000) turns the zero-width slits into open removed regions so that the exclusion
can *also* be confirmed by winding numbers; both checks are reported.

**Two Bost–Charles evaluators.**
1. `bcdirect.py`: the model-free double sum (1/N²)Σ log|(φ(z_j)−φ(z_k))/(z_j−z_k)| with the
   diagonal log|φ′|, valid because ∬log|z−w| dμdμ = 0.
2. `bcfast.py`: an exact **Jensen reduction**. φ = h∘Φ has 0 as its only zero in 𝔻, so
   ∫_𝕋 log|φ| dμ = log|φ′(0)|; and Jensen at a = φ(z) gives

   $$\boxed{\;\mathrm{BC}(\varphi)=\log|\varphi'(0)|+\int_{\mathbb T}S(z)\,d\mu(z),\qquad S(z)=\sum_{w\in\mathbb D,\;\varphi(w)=\varphi(z)}\log\frac1{|w|}\;}$$

   with the interior solutions w = Φ^{−1}(η), η running over the Γ₀(2)-orbit points of Φ(z)
   that lie in Ω. **[proved]** So BC − log|φ′(0)| is a pure multivalency ("shape") term,
   zero exactly for univalent φ — which recovers CDT's univalent special case BC = log ρ.
   *Numerical caveat (found and fixed):* Φ^{−1} loses all precision for tiny arguments (the
   composed square roots collapse to u = i exactly), so for |η| < 10⁻⁴ the stable branch
   Φ^{−1}(η) = η/Φ′(0) + O(η²) with one Newton correction is used. Without this the Jensen
   value silently under-reports by up to 0.9 on contours whose boundary passes near a cusp
   where h → 0.

*Agreement of the two evaluators* (`verdict.py`): on concentric discs |Q| < r (host C of
`lattice/catalan_mu4`) to 6 decimals; on the deep-slit place-v₂ designs, Jensen 6.3876 / 6.6229 /
6.8175 against the direct sums' N-extrapolations 6.394 / 6.638 / 6.830 at R = 0.75 / 0.83 / 0.90.

**Correction (i): CDT's Lemma A.4.4 fails for their published parameters.**
`calib_cdt.py`. The h-preimages of −1/72 in the disc are, by increasing modulus,
5.418·10⁻⁵ (principal, the only real one), then the conjugate pairs

| bottom row (2c,d) | \|z\| | arg/2π |
|---|---|---|
| (2, ±1) | 0.401921 | ∓0.453604 |
| (2, ±3) | 0.592651 | ∓0.420115 |
| **(2, ±5)** | **0.753925** | ∓0.428113 |
| (4, ±1) | 0.782767 | ∓0.243767 |

CDT quote 0.782767 correctly but write "there are only 4 such preimages that we need to
exclude"; the pair at 0.753925 is inside R = 77/100 and was missed. The argument principle
confirms **7** solutions of h = −1/72 in |z| < 0.77. Four slits cannot exclude six points, and
with the published θᵢ the slit tips land at arg/2π = −0.4188, −0.3055, −0.2809 while the bad
points sit at ±0.4536, ±0.42011, ±0.42811: **five of the six have winding number +1** about
their contour. **[verified]**

**Correction (ii): BC of the published φ.** 12.120676 (Jensen, N = 16384) and 12.1239 (direct,
N = 4096, still decreasing), against the printed 11.845. **[verified]**

### 1.4 CDT's constants are nevertheless exactly right  **[verified]**

Running the successive-slitting designer with CDT's own base radius, lune and shrink and six
slits reaching 0.1 % past the six bad preimages:

| | this task | CDT |
|---|---|---|
| R, c, shrink | 77/100, 75/10, 995/1000 | 77/100, 75/10, 995/1000 |
| \|ψ′(0)\| | **0.6293401** | 0.6292232680 |
| BC | **11.8449** | 11.845 |
| bound m ≤ | **13.9906** | 13.9938 |
| margin 14(L−τ)−BC | **+0.0080** | +0.0053 |

So CDT's pair (|ψ′(0)|, BC) lies exactly on the curve traced by this family at their own
parameters: the theorem is fine, the appendix's θᵢ are not. Optimising R in the same family:

| R | \|ψ′(0)\| | BC | entry | bound | margin |
|---|---|---|---|---|---|
| 0.77 | 0.6328104 | 11.9133 | 0.852134 | 13.9806 | +0.0166 |
| **0.80** | **0.6461355** | **12.2007** | **0.872972** | **13.9761** | **+0.0209** |
| 0.83 | 0.6582652 | 12.4953 | 0.891571 | 14.0149 | −0.0133 |

---

## 2. The host, the sources and the period ξ — re-verified  **[verified, 142 digits]** (`xi_check.gp`)

Zagier **D** on Γ₁(5): x = q∏(1−qⁿ)^{5(n|5)}, F = ΣAₙxⁿ with
Aₙ = Σ_k C(n,k)²C(n+k,k) = 1, 3, 19, 147, 1251, 11253, 104959, 1004307, …,
(n+1)²A_{n+1} = (11n²+11n+3)Aₙ + n²A_{n−1}; λ₁ = φ⁵, λ₂ = −φ⁻⁵, λ₁λ₂ = **−1** (a unit);
singular points **t₁ = φ⁻⁵ = 0.0901699…** (the fold) and **t₂ = −φ⁵ = −11.0901699…**, k = 2.

ψ₄ = the odd quartic character mod 5, ψ₄(2) = i; E₃^{ψ,χ} = Σ_n (Σ_{d|n} ψ(n/d)χ(d)d²) qⁿ;
R₁ = E(1,ψ₄)+E(1,ψ̄₄), R₂ = i(E(1,ψ₄)−E(1,ψ̄₄)), R₃ = E(ψ₄,1)+E(ψ̄₄,1), R₄ = i(E(ψ₄,1)−E(ψ̄₄,1)).

* **Φ_D = ½R₁ + R₂**, q-expansion 0, 1, −7, 19, −23, 1, 47, −97, 105, −62, … — coefficient by
  coefficient Table `tab:sources` row **D**; companion B_D = F·D_q^{−2}Φ_D peeled in x gives
  Bₙ/Aₙ → **ζ(2)/5** with |Bₙ/Aₙ − ζ(2)/5| < 10⁻¹⁷³ at n = 118. **[verified, 173 digits]**
* **Φ_new = R₃ + φ⁵R₄ = (1+iφ⁵)E₃^{ψ₄,1} + (1−iφ⁵)E₃^{ψ̄₄,1}**, coefficients in ℤ[φ]
  (0, 2, −14.18033988…, 40.18033988…, −58.72135954…, 50, … = 0, 2, −(3+5φ)·… in ℤ[φ]);
  companion B_new gives

  $$\xi=\lim_n\frac{B_{{\rm new},n}}{A_n}=0.6556341888406567663309814138723994024111384591367020284631349608915842669423313859875149423226475996172012313287272190533428707991338394140091\ldots$$

  with 142 digits of agreement between n = 117 and n = 118, and
  |ξ − (φ⁵·Im L(2,ψ₄) − Re L(2,ψ₄))| < 7.2·10⁻¹⁴³. `lindep` on
  [ξ, 1, ζ(2), Re L, Im L, √5, √5·Re L, √5·Im L] returns **[2, 0, 0, 2, −11, 0, 0, −5]**,
  i.e. 2ξ + 2 Re L − (11+5√5) Im L = 0 and **nothing else** — no rational part, no ζ(2)
  admixture, no elementary π² term. **[verified, 142 digits]**

Re L(2,ψ₄) = 0.95871612271688315539193642933117852641597153075829606867240444791039545605986570951600…,
Im L(2,ψ₄) = 0.14556587678508959046170451181198645372080514688909901383422661119255914771670644486343….

This reproduces `lattice/hostscan/REPORT.md` §10.4 by an independent run.

---

## 3. The descent, the pure module over K, and the fourteen functions [Task 3]

**The descent is forced.** The conditional function H = aA + bB_D + cB_new is (under the
hypothesis a + b·ζ(2)/5 + c·ξ = 0, a,b,c ∈ K) regular at the fold t₁ and singular at t₂; so the
pure module lives on **P**¹∖{0, s, ∞} with **s = t₂ = −φ⁵**, and the CDT normaliser descent is
w(x) = sx/(x−s), y = x + w(x) = x²/(x−s), branch point y = 4s, extra point
y(t₁) = t₁²/(t₁−t₂). Since t₂ is a *unit* of ℤ[φ], the polylogarithms Li_j(x/s) = Σ(x/s)ⁿ/nʲ
keep coefficients in ℤ[φ] and the denominator types of the descent are CDT's, with
[1..n]^k ↦ [1..2n]^k. Normalised coordinate **Y := y/s**, so the ℤ/2 point is at Y = 4 and CDT's
own series for B₁,…,B₇ apply verbatim.

**The conditional ODE.** With θ = x d/dx the row operator is
θ² − x(11θ²+11θ+3) − x²(θ+1)², i.e.

$$\mathcal L\,y \;=\; x(1-11x-x^2)\,y'' \;+\; (1-22x-3x^2)\,y' \;-\; (3+x)\,y ,$$

singular at 0, t₁, t₂, ∞ (1−11x−x² = −(x−t₁)(x−t₂)). Verified on the exact series (221
coefficients, `farcusp.py`):

  **L A = 0**,  **L B_D = 1**,  **L B_new = 2/(1 − x/t₂)**,  **L B′_new = 2/(1 − x/t₁)**,

so the conditional generator satisfies L H = b + 2c/(1 − x/t₂) — the exact analogue of CDT's
b + c/(1−x) (their Prop. 11.1.4). More generally L(c₃B₃+c₄B₄) = 2(c₃−c₄x)/(−(x−t₁)(x−t₂)):
the pole at t₁ cancels iff c₃ = c₄t₁, the pole at t₂ iff c₃ = c₄t₂. **[verified]**
Sharp denominators on the x-line: d_n²B_{D,n} ∈ ℤ and d_n²B_{new,n} ∈ ℤ[φ] for all n ≤ 200,
with d_n¹ failing at 199 resp. 197 of the 200 indices — type exactly [1..n]², i.e. **k = 2**
(`task2/`). **[verified, n ≤ 200]**

### 3.1 The measured denominator array  **[verified, exact rationals, n ≤ 80]** (`task3/`)

Every clause of CDT's Lemma `bdenominators` was re-measured on the exact rational series in
Y (n ≤ 80) and every one is **valid and sharp**:

| | measured minimal type | CDT's claim |
|---|---|---|
| B₁ | trivial | trivial ✔ |
| B₂ | [1..2n] | ✔ |
| B₃ | [1..2n]·n | ✔ |
| B₄ | [1..2n]² (and [1..n][1..2n] **fails**) | ✔ |
| B₅ | [1..2n]², sharper [1..2n](2n−1) | ✔ |
| B₆ | [1..2n]n² (a fortiori [1..2n]²n) | ✔ |
| B₇ | [1..2n]²n (and [1..2n]² alone **fails**) | ✔ |

s = −φ⁵ is a unit (s⁻¹ = (11−5√5)/2), and D_n c_n s^{−n} ∈ ℤ[φ] was verified coefficientwise
for all seven, n ≤ 80. For the full fourteen (n ≤ 50) the measured array is
**u₁ = 1** (only B₁), **u₂ = 3** (B₁, B₂, B₃), **b₁ = b₂ = 2**,
**e = (0,0,1; 0⁶; 1,1,1,1,1)**, Σe = 6, max e = 1 — **identical to CDT's**, hence

$$\boxed{\;\tau=\tau^\flat+\tau^\sharp=\tfrac{191}{49}+\tfrac{27}{80}=\tfrac{16603}{3920}\;\text{unchanged.}}$$

The one genuine alternative — demoting B₆ to a single layer with e = 2, giving u₂ = 4,
τ♭ = 375/98 — loses, because τ♯ rises to 0.54018: net **+0.13125** on τ.

### 3.2 The conditional side over K  **[verified]** (`task3/`)

The symmetrisation recurrence for a general s is P₀ = 2, P₁ = y, **Pₙ = y P_{n−1} − s y P_{n−2}**
(x + w = y, x·w = s y), proved and verified as exact K[[x]] identities to order 42; the
substitution Qₙ(Y) = s^{−n}Pₙ(sY) cancels s exactly, so CDT's integer polynomials apply
verbatim in the normalised coordinate. End-to-end, G(Y(x)) = H(x) + H(w(x)) exactly to order 42.
**G = Sym⁺H has sharp denominator type [1..2n]² over ℤ[φ]** (n ≤ 60; [1..n][1..2n] fails), for
the generic (a,b,c) = (1, −3+φ, 5−2φ) and for (0,1,0), (0,0,1), (φ,φ²,φ³) separately.
Apéry limits reconfirmed independently: ζ(2)/5 to 60 digits, ξ to all 40 digits given.

### 3.3 K(y)-linear independence of the fourteen  **[verified]** (`task3/`)

Rank of {Yʲf_i : i ≤ 14, j ≤ D} at a **split** prime p = 2305843009213694009 (≡ 4 mod 5,
ℤ[φ]/𝔭 = 𝔽_p) and at an **inert** prime p = 2305843009213693967 (≡ 2 mod 5, 𝔽_p[t]/(t²−5)):

| deg P_i ≤ | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| rank / needed (both primes) | 14/14 | 28/28 | 42/42 | 56/56 | 70/70 | 84/84 |

series orders 26…96. Certified: **no K(y)-relation of degree ≤ 5**. Unbounded-degree
independence (CDT's Lemma 12.1.1, a monodromy argument) is **[assumed]**, not transported.
The same code path reproduces `CDT_FINDER.md` §2 on CDT's own host exactly.

*Two incidental discrepancies, neither affecting τ:* CDT's Remark `central binomial`
(replacing [1..2n] by n(n−1)C(2n,n)) holds for B₂,B₃,B₅,B₆ but **fails for B₄ and B₇** —
only one of the two layers relaxes; and `lattice/cdt_finder/indep_check2.py` items 13–14 carry
an index shift relative to CDT's definitions (both conventions give full rank).

### 3.4 What the second real place demands of the *hypothesis* — the far-cusp periods
**[verified, 12 digits]** (`farcusp_data.gp`, `farcusp.py`)

At v₂ the pure module is B_i(y/σ(s)) = B_i(y/t₁), so **x = t₁ = σ(s) is where singularities are
allowed** (it maps to y = ∞) and **x = t₂ is the point that must be removed** — t₂ is the
v₂-fold. Two checks that this is the right bookkeeping, and that the *sources* are perfectly
Galois-equivariant:

* the right-hand sides. L(c₃B₃+c₄B₄) = 2(c₃−c₄x)/(−(x−t₁)(x−t₂)); the pole at t₁ cancels iff
  c₃ = c₄t₁, at t₂ iff c₃ = c₄t₂. For Φ_new, (c₃,c₄) = (1, φ⁵) = (1, −t₂) and c₄t₁ = −t₁t₂ = 1:
  the t₁-pole cancels, leaving **L B_new = 2/(1−x/t₂)**, a pole at t₂ = s ↦ y = ∞ — harmless at
  v₁. For the conjugate, (σc₃,σc₄) = (1, −t₁) and σ(c₄)t₂ = −t₁t₂ = 1: the **t₂**-pole cancels,
  leaving **L B′_new = 2/(1−x/t₁)**, a pole at t₁ = σ(s) ↦ y = ∞ — harmless at v₂. Verified to
  8 coefficients against 2(1/t₂)ⁿ resp. 2(1/t₁)ⁿ. **[verified]**
* the cusps. With x(∞, 0, 1/2, 2/5) = (0, φ⁻⁵, −φ⁵, ∞), the cusp over t₁ is 0 and over t₂ is 1/2;
  the constant terms are Φ_D = (0,−1,0,0), Φ_new = (0,0,0,∗), Φ′_new = (0,0,∗,0). So Φ_new
  vanishes at the cusp over t₁ (the v₁-fold) and **Φ′_new vanishes at the cusp over t₂ (the
  v₂-fold)**. **[verified, `task2/`]**

**The periods.** For f ∈ {A, B_D, B′_new} the monodromy difference Δ_f = (M_{t₂}−1)f solves the
homogeneous equation and is M_{t₂}-fixed, hence spans the line ⟨u_{t₂}⟩; the far-cusp periods are
π_f := Δ_f/Δ_A. Computed by numerical continuation around t₂ (base point x = 0.05, 221 exact
series coefficients, DOP853 rtol 10⁻¹²; the two components of Δ agree to 5·10⁻¹⁴, confirming
collinearity). Calibration: the same construction at the **near** cusp t₁ returns
+0.328986813370 = ζ(2)/5 for B_D and +0.655634188841 = ξ for B_new, exactly the Apéry limits
(and an *inconsistent* pair of ratios for B′_new, which is the log² at t₁). Then

$$\boxed{\;\pi_{B_D}^{(t_2)}=-\tfrac{11}{5}\zeta(2)\quad(\text{ratio to }\zeta(2)/5\ \text{is}\ -11.0000000000),\qquad \pi_{B'_{\rm new}}^{(t_2)}=\xi'=-\mathrm{Re}\,L-\varphi^{-5}\mathrm{Im}\,L=-0.971841789639\;}$$

(ratio to ξ′ equal to 1.000000000000). Both ratios are stable to **12 digits** across base
points 0.05, 0.02, −0.03, loop radii 4, 5, 7, path heights 1.5, 2, 3 and tolerances 10⁻¹², 10⁻¹³
(`farcusp_robust.txt`), and were independently recomputed in `task2/` to **174 digits** by
adaptive Taylor continuation at 200-digit precision (residuals |π_D + 11ζ(2)/5| = 2.1·10⁻¹⁷⁴,
|π′ − ξ′| = 1.6·10⁻¹⁷⁴; `lindep` certificates [−5,0,−11,0,0,0,0,0] and [−2,0,0,−2,11,0,0,−5]),
with the near-cusp calibration returning ζ(2)/5 and ξ to 2·10⁻²¹² and a third, wholly different
confirmation to 5–6 digits by Abel summation of the Eichler integral at the cusps
(`task2/19_abel.gp`). So the two conditions the architecture needs are

  **(i)** a + (b/5)·ζ(2) + c·ξ = 0  (the target hypothesis — regularity of H at the v₁-fold t₁);
  **(ii)** σ(a) − (11σ(b)/5)·ζ(2) + σ(c)·ξ′ = 0  (regularity of σH at the v₂-fold t₂).

With ξ = −Re L + φ⁵ Im L and ξ′ = −Re L + σ(φ⁵) Im L, condition (ii) is σ applied to the
K-coefficients of (i) **except** in the ζ(2) slot, where +1/5 becomes −11/5:
  (ii) = σ(i) − (12/5)σ(b)·ζ(2).
Hence **(i) does not imply (ii)**. More: in the K-basis (1, ζ(2)/5, Re L, Im L) the two
coefficient vectors are Λ = (a, b, −c, cφ⁵) and Λ′ = (σa, −11σb, −σc, σ(cφ⁵)), and
Λ′ = μΛ for μ ∈ K^× forces (3rd and 4th entries) σ(φ⁵) = φ⁵ unless c = 0, while Λ′ = Λ forces
b = −11σ(b) ⇒ b = 0 and cφ⁵ = −cφ⁻⁵ ⇒ c = 0. **The two relations are never the same relation
whenever ξ is actually involved.** Splitting Λ = Λ₀ + √5 Λ₁ over ℚ, "(i) and (ii)" says

  Λ₀·v = (6/5)σ(b)ζ(2)  **and**  √5 Λ₁·v = −(6/5)σ(b)ζ(2),  v = (1, ζ(2)/5, Re L, Im L),

two independent ℚ-relations among 1, ζ(2), Re L, Im L, √5, √5ζ(2), √5 Re L, √5 Im L, where one
K-relation supplies one. (In the clean sub-case b = 0 — irrationality of ξ over K — it says the
rational part and the √5 part of the single K-relation a + cξ = 0 must vanish *separately*.)
So one K-linear relation among 1, ζ(2)/5, ξ gives a conditional function at **one place only**.
**[verified, 174 digits]**

### 3.5 Where the −11 comes from — the orientation, not the Galois action
**[proved; character sums verified exactly]** (`task2/20_minus11.gp`, `L3_identity.gp`)

The asymmetry is elementary and lives entirely in the *outer* (ζ(2)) direction.

**Outer orientation (1, ψ): dilogarithms.** With w := Re ψ₄ − 2 Im ψ₄ (period 5, values
1, −2, 2, −1, 0, mean zero), Φ_D = Σ_n (Σ_{d|n} w(d)d²)qⁿ and, summing over n = de,
D⁻²Φ_D = Σ_{d≥1} w(d)·Li₂(q^d). For a P-periodic mean-zero sequence,
lim_{x→1⁻} Σ a(d)x^d = −(1/P)Σ_{j=1}^{P} j·a(j). Hence

* **cusp 0** (q → +1): Li₂(q^d) → ζ(2) for *every* d; Σ_{j≤5} j·w(j) = 1−4+6−4 = **−1**, so the
  factor is −(−1)/5 = **1/5** and π^{(t₁)}(Φ_D) = **ζ(2)/5**;
* **cusp 1/2** (q → −1): Li₂((−1)^d) = ζ(2) for d even but **−ζ(2)/2** for d odd — the weight is
  *parity-dependent*. w·g (g = 1 on evens, −1/2 on odds) has period 10, values
  (−1/2, −2, −1, −1, 0, 1, 1, 2, 1/2, 0), still mean zero, and Σ_{j≤10} j·(wg)(j) = **22**, so the
  factor is −22/10 = **−11/5** and π^{(t₂)}(Φ_D) = **−11ζ(2)/5**.

Both sums reproduced exactly here (`L3_identity.txt`). **Nothing in this involves σ** — which is
precisely why π_D is not σ(ζ(2)/5).

**Inner orientation (ψ, 1): Lambert series, and equivariance is a theorem there.** With
v := 2Re ψ₄ − 2φ⁵ Im ψ₄, D⁻²Φ_new = Σ_{k≥1} (v(k)/k²)·q^k/(1−q^k). Using
1/(e^t−1) = 1/t − 1/2 + t/12 − … with t = 2πyk: at cusp 0 every k contributes weight **−1/2**;
at cusp 1/2 the even k give −1/2 and the odd k give −x^k/(1+x^k) → **−1/2** as well. The weight
is *uniform in k at both cusps*, so the inner period is −½L(2,v), **K-linear in the source's
K-coefficients**, and therefore

$$\pi^{(t_2)}(\sigma\Phi)=\sigma\bigl(\pi^{(t_1)}(\Phi)\bigr)\quad\text{on the inner directions, with no twist.}$$

So the equivariance predicted in §3.4 **is a theorem on the ξ-line**; the failure is confined to
the ζ(2)-line and is a feature of Zagier D's orientation, not of the arithmetic of K.
(The diamond ⟨2⟩ observation — c₀(∞) = 0 is not ⟨2⟩-stable — is the same phenomenon seen
qualitatively; this is the proof.) That 11 = φ⁵ − φ⁻⁵ = L₅ = Tr_{K/ℚ}(φ⁵) also equals −(t₁+t₂)
is **not** claimed to be more than the value of a character sum at level 5.

**A by-product: a weight-three companion to the fold-regularity ratio.**
Both convergence conditions above (Σ_k v(k)/k³ = 0 at cusp 0, Σ_{k even} v′(k)/k³ = 0 at cusp
1/2) reduce to the *same* statement via Re ψ₄(2j) = −Im ψ₄(j), Im ψ₄(2j) = Re ψ₄(j) (because
ψ₄(2) = i), namely

$$\boxed{\ \mathrm{Re}\,L(3,\psi_4)\;=\;\varphi^5\,\mathrm{Im}\,L(3,\psi_4)\ }$$

Verified here independently to **212 digits**: Re L(3,ψ₄) = 0.98819168162405719379746955…,
Im L(3,ψ₄) = 0.08910518834573959516382505…, ratio − φ⁵ = 0 to 2.4·10⁻²¹²; `lindep` on
[Re L₃, Im L₃, √5 Im L₃] returns **[−2, 11, 5]**, i.e. 2Re L₃ = (11+5√5)Im L₃. This is the s = 3
(critical, for odd ψ) shadow of the weight-two fold-regularity ratio d(R₃)/d(R₄) = −φ⁵ of
`hostscan` §10.3 — the cusp constant-term vanishing re-read on the L-value side. By contrast
Φ_new at cusp 1/2 would need Re L₃ = −φ⁻⁵ Im L₃, and −2Im L₃ − 2φ⁵Re L₃ = −22.0966… ≠ 0: that
non-vanishing **is** the log² of §3.4. **[verified, 212 digits]**

Two ways out, both priced in §5.2:
* assume the *doubled* hypothesis (i) ∧ (ii) — two ℚ-linear relations among
  1, √5, ζ(2), √5ζ(2), Re L, √5 Re L, Im L, √5 Im L. The conclusion is then weaker than
  "1, ζ(2), ξ independent over K";
* keep the single hypothesis and make Ω₂ avoid **every** preimage of Y₂, the principal one
  included — then σG pulls back holomorphically whatever it does at Y₂, and no second relation
  is needed. The principal preimage sits at |z| = 0.074191, so the cap drops to 0.257187.

---

## 4. The contours at the two real places [Task 4]  **[exact geometry / verified contours]**

### 4.1 The two configurations, exactly

With t₁t₂ = −1, t₁+t₂ = −11, t₁−t₂ = 5√5, φ¹⁵ = 682 + 305√5 (L₁₅ = 1364, F₁₅ = 610):

| | place **v₁** (√5 ↦ +2.2360679…) | place **v₂** (√5 ↦ −2.2360679…) |
|---|---|---|
| s realised as | t₂ = −11.0901699… | t₁ = +0.0901699… |
| \|s_v\| | φ⁵ = 11.0901699… | φ⁻⁵ = 0.0901699… |
| fold (removed singularity) realised as | t₁ = +0.0901699 (**inner**) | t₂ = −11.0901699 (**outer**) |
| branch point 4s | −44.360680 | +0.360680 |
| extra point y_v = t²/(t−s) | 7.2722464948·10⁻⁴ | −11.000727224649 |
| **Y_v := y_v/s_v** | **(682√5−1525)/25 = −6.5573805732·10⁻⁵** | **−(682√5+1525)/25 = −121.999934426194** |
| log(256\|s_v\|) | 7.951237 | 3.139118 |

**[exact]** Y₁ = −φ⁻¹⁵/(5√5), Y₂ = −φ¹⁵/(5√5); Y₁Y₂ = 1/125, Y₂/Y₁ = φ³⁰; σ(Y₁) = Y₂. The
average of the two ceilings is log(256·|N(s)|^{1/2}) = **log 256 = 5.545177**, the Apéry-perfect
value: the *ceiling* is exactly CDT's. Everything below is the loss below that ceiling.

### 4.2 The bad h-preimages

Ω_v ⊂ 𝔻 must contain **exactly one** h-preimage of Y_v (the principal one, where the conditional
function is regular under the hypothesis); all others must be excluded (CDT's Lemma A.4.4 for
their −1/72). By increasing modulus:

| place | principal preimage | non-principal preimages (\|z\|, conjugate pairs) |
|---|---|---|
| v₁ | z₀ = 2.56146·10⁻⁷ (Im τ₀ = 2.4157, deep in the cusp) | **0.536032**, 0.625433, 0.730538, 0.810740, … |
| v₂ | z₀ = 0.0741911 (Im τ₀ = 0.4156, shallow) | **0.213693**, 0.499025, 0.695729, 0.764482, … |
| *CDT* | *5.418·10⁻⁵* | *0.401921, 0.592651, 0.753925, 0.782767, …* |

This is the whole story: at v₁ the fold is hyperbolically **deeper** than CDT's, so its orbit sits
**further out** and the slits are cheaper; at v₂ it is far shallower, so the orbit comes close to
the centre. The classical extremal domain omitting one point p (radial slit) gives the rigorous
cap |ψ_v′(0)| ≤ 4|p|/(1+|p|)²:

| place | deepest bad preimage | cap on \|ψ_v′(0)\| | cap on L_v |
|---|---|---|---|
| v₁ | 0.536032 | **0.908762** | 7.855565 |
| v₂ | 0.213693 | **0.580274** | 2.594864 |

*A consistency check the unit property provides.* One may equally run the bound in the
normalised coordinate Y = y/s, in which the pure functions have **rational** coefficients and
the ceiling is log 256 at *each* place; the two coordinates rescale φ_v by |s_v|, changing both
L_v and BC_v by log|s_v|. Since s = t₂ is a **unit**, |s₁s₂| = 1 and the two changes cancel in
the average: L̄ and B̄C — hence the bound — are the same in the y- and Y-coordinates.
(Numerically: L̄^{(Y)} = ½(log(256·0.703784) + log(256·0.394939)) = 4.90505 = L̄^{(y)}.)

### 4.3 The trade-off and the optimum

Family: base disc of radius R, one lune of parameter c at the cusp z = 1 (where |h| → ∞), one
curved slit per bad preimage inside |z| < R, order of slits optimised over all permutations
(the successive-slit construction is *exactly* extremal in the symmetric benchmarks
𝔻 minus n equally spaced radial slits, where the closed form (4rⁿ/(1+rⁿ)²)^{1/n} is known:
ratio 1.0000 for n = 2,3,4). Lunes at further cusps (e^{±2πi/3}, e^{±2πik/5}) were tried and
always lose. Scanning R ∈ [0.40, 0.94] and c ∈ {5, 6.5, 7.5, 9, 12} (`two_places.py`), the
objective 14L_v − BC_v is unimodal in R and peaks at

| place | R | c | #bad | \|ψ_v′(0)\| | L_v | BC_v | shape_v = BC−L | 14L−BC |
|---|---|---|---|---|---|---|---|---|
| v₁ | 0.78 | 6.5 | 6 | 0.704596 | 7.601106 | 15.23969 | 7.63859 | **91.17577** |
| v₂ | 0.76 | 7.5 | 6 | 0.395505 | 2.211527 | 6.42751 | 4.21594 | **24.53386** |

**Certified contours** (slit tips 0.05 % past each bad point, global shrink ρ = 0.999, so that
the removed slits are open regions; every bad preimage verified excluded by
|Φ^{−1}(p)| = 1.001001 ≥ 1 *and* by winding number 0, principal preimage verified inside):

| place | \|ψ_v′(0)\| (exact rational, `rad_exact`) | L_v | BC_v | shape_v |
|---|---|---|---|---|
| **v₁** | **0.703784137176** | 7.5999530 | 15.223765 | 7.623812 |
| **v₂** | **0.394939005291** | 2.2100944 | 6.408897 | 4.198803 |

Certification costs only 0.0013 in L̄. BC computed at N = 16384 Fourier nodes, Γ₀(2)-orbit
truncated at |2c|,|d| ≤ 30 (converged: cdmax 14 → 120 changes nothing beyond 10⁻⁵).

---

## 5. The bound  **[verified]**

$$\bar L=\tfrac12(L_1+L_2)=\mathbf{4.9050237},\qquad \overline{\mathrm{BC}}=\tfrac12(\mathrm{BC}_1+\mathrm{BC}_2)=\mathbf{10.816331},\qquad \tau=\tfrac{16603}{3920}=4.2354592 .$$

* **Entry (holonomicity, Remark BCboundK):** L̄ = 4.9050 > σ_m = 4 ✔
* **Entry (the bound's denominator):** **L̄ − τ = +0.6695645 > 0** ✔ — the entry test *passes*,
  with 0.67 nats to spare (CDT's own is +0.8730).
* **Bound:** m ≤ 10.816331/0.6695645 = **16.15428**.
* **Exhibited:** m = 14. **Margin = 14·(L̄−τ) − B̄C = 9.373903 − 10.816331 = −1.44243.**

(The "ideal" contour, slits exactly through the points and no shrink, gives L̄ = 4.9063163,
B̄C = 10.833579, bound 16.14886, margin −1.44158 — the certification is not what costs.)

**Attribution.** Write shape_v = BC_v − L_v and s_v = 13 log|ψ_v′(0)| − shape_v, so that
margin = 13 log 256 + ½(s₁+s₂) − 14τ:

| term | value |
|---|---|
| 13 log 256 | +72.08731 |
| 13·mean log\|ψ′\| (log\|ψ₁′\| = −0.351284, log\|ψ₂′\| = −0.929024) | −8.32200 |
| −mean shape (shape₁ = 7.62381, shape₂ = 4.19880) | −5.91131 |
| −14τ | −59.29643 |
| **margin** | **−1.44243** |

and per place, against CDT's own host in the same accounting (their in-family optimum
|ψ′(0)| = 0.646136, shape = 6.99987, s = −12.76997, margin +0.02091):

| | s_v | vs CDT |
|---|---|---|
| **v₁** | **−12.19050** | **better by +0.57947 nats** |
| **v₂** | **−16.27611** | **worse by −3.50615 nats** |
| average | −14.23331 | −1.46334 nats of margin |

**Place v₂ is entirely responsible for the failure**; place v₁ is genuinely *better* than CDT's
own geometry, because its fold sits deeper in the cusp.

**What v₂ would have to deliver.** Holding v₁ at its optimum, margin 0 needs
14L₂ − BC₂ ≥ 27.4173 against the achieved 24.5324 — a deficit of **2.885 nats**. At unchanged
BC₂ that is |ψ₂′(0)| = 0.4853, against 0.3949 achieved. The crude one-point cap (0.5803) does
not by itself forbid this. What does is §5.1.

### 5.1 The frontier saturates well short of what is needed  **[verified]** (`saturation.txt`)

Cutting **all** the bad preimages (R → 1, up to 134 slits) the conformal radius saturates:

| R | v₁: #bad, \|ψ₁′(0)\|, 14L−BC | v₂: #bad, \|ψ₂′(0)\|, 14L−BC |
|---|---|---|
| 0.94 | 30, 0.804958, 90.2478 | 30, 0.405685, 24.3336 |
| 0.97 | 64, 0.817462, 89.9429 | 68, 0.406501, 24.2636 |
| 0.985 | 130, **0.822756**, 89.7523 | 134, **0.406809**, 24.2260 |

So within this family |ψ₂′(0)| **cannot exceed ≈ 0.4068** — far below the 0.4853 required — and
|ψ₁′(0)| cannot exceed ≈ 0.8228. Even *freezing* the shape terms at their optimal values and
pushing both conformal radii to saturation, the margin only rises to
−1.442 + 13·½(log(0.4068/0.39494) + log(0.82276/0.70378)) = **−0.235**: still negative. Since
the shape terms in fact grow with R (the 14L−BC columns above *fall*), the true optimum over the
family is the decoupled one of §4.3, margin **−1.4416**. The verdict is not marginal.

### 5.2 The two branches, and which hypothesis is being contradicted

| branch | hypothesis assumed | L̄ | B̄C | bound | **margin** |
|---|---|---|---|---|---|
| **A** | (i) ∧ (ii) — the *doubled* hypothesis (§3.4) | 4.90502 | 10.81633 | 16.154 | **−1.4424** |
| **B** | (i) alone, Ω₂ avoiding *all* preimages of Y₂ | 4.59778 | 9.47631 | 26.154 | **−4.4038** |

Branch B (`single_hypothesis.py` → `single_hypothesis.txt`): the extra excluded point at
|z| = 0.074191 caps |ψ₂′(0)| at 0.257187 and the best achieved is 0.213384, L₂ = 1.594459,
BC₂ = 3.71292, 14L₂−BC₂ = 18.6095 against v₁'s 91.1758. So the variant that survives on the
single, *correct* hypothesis is three nats worse still. **Both branches are negative.**

**Convexity.** CDT's convexity improvements take their own bound 13.9938 → 13.730 → 13.7206 →
13.678 → 13.621 (four radii), a factor 0.97336. Remark BCboundK says these "extend in the
obvious way"; transported, our bound becomes **15.72**, still far above 14 (equivalent margin
−1.15). A 13.3 % improvement would be needed, against CDT's realised 2.7 %.

---

## 6. Could anything rescue it?  **[verified]**

**(a) A bigger inventory.** `INVENTORY_BOUND.md` (R1) proves CDT's fourteen are the complete
admissible supply at max eᵢ ≤ 1 on this descent orbifold; going beyond forces max eᵢ ≥ 2 and
τ♯ jumps. With the contour re-optimised for each m from the measured curves:

| added functions (e = 2) | m | τ | entry | bound | margin |
|---|---|---|---|---|---|
| 0 | **14** | 4.23546 | 0.67643 | 16.1337 | **−1.4433** |
| 1 | 15 | 4.44815 | 0.46585 | 23.4914 | −3.9558 |
| 2 | 16 | 4.47487 | 0.44830 | 24.7193 | −3.9089 |
| 4 | 18 | 4.49846 | 0.44381 | 25.7078 | −3.4208 |
| 7 | 21 | 4.49926 | 0.45562 | 25.5920 | −2.0923 |

and every m < 14 is worse (m = 13: −1.906; m = 12: −2.275; m = 10: −3.180; m = 8: −4.285).
**m = 14 with CDT's exact array is the optimum, exactly as on CDT's own host.**

**(b) A better contour at v₂.** The successive-slit construction is extremal in the symmetric
benchmarks; a non-injective ψ cannot beat the Riemann map of the extremal domain (Schwarz); the
exclusion set is forced by the conservative reading of CDT's Lemma A.4.4; lunes are strictly
worse than slits for excluding single points (0.409 vs 0.580 at |p| = 0.2137); extra lunes at the
q = 3, 5 cusps lose.

**The one identified loophole, quantified** (`whatif.py` → `whatif.txt`, `whatif_cdt.txt`).
CDT — and this report — exclude *every* non-principal h-preimage of Y_v. In fact only those at
which the analytically continued conditional function actually acquires its logarithm need be
excluded: the lift is singular at the preimage γ·z₀ iff ρ(γ)⁻¹v ∉ ker(M−1), where M is the local
monodromy at the fold and v the conditional vector — a codimension-one condition, generically
but not necessarily violated. Recomputing the optimum with the k deepest preimages *assumed
harmless*:

| dropped at v₁ | dropped at v₂ | \|ψ₁′\| | \|ψ₂′\| | **margin** |
|---|---|---|---|---|
| 0 | 0 | 0.6912 | 0.3955 | **−1.446** |
| 0 | **2** | 0.6912 | **0.6725** | **+0.512** |
| 0 | 4 | 0.6912 | 0.7227 | +0.651 |
| 2 | 0 | 0.7125 | 0.3955 | −1.403 |
| 4 | 4 | 0.7243 | 0.7227 | +0.711 |

So **everything turns on the single conjugate pair of preimages at \|z\| = 0.213693 at place v₂**
(bottom rows (2,±1), i.e. the images of the principal preimage under
$\left(\begin{smallmatrix}1&0\\ \pm2&1\end{smallmatrix}\right)$): if they are harmless the
margin is **+0.51** and the theorem follows; if not, it is **−1.44**. Dropping preimages at v₁
buys almost nothing (+0.04), because v₁'s geometry is already better than CDT's. The next
paragraph settles it: **they are not harmless.**

**Closure of the loophole** (`monodromy.py` → `monodromy.txt`, `monodromy_robust.txt`).
Let ρ be the monodromy of the rank-4 system ⟨A, A₂, B_D, B_new⟩ on **P**¹∖{0,t₁,t₂,∞}. For any
deck element γ, ρ(γ)H − H is a solution of the *homogeneous* equation (H and ρ(γ)H satisfy the
same inhomogeneous equation with single-valued rational right-hand side), and ρ(γ)H is regular
at the fold iff that difference lies in the line ⟨u_{t₁}⟩ of homogeneous solutions holomorphic
at t₁. For γ = the parabolic generator at the cusp y = ∞ — bottom row (2,±1), which is exactly
the deepest bad pair at **all three** configurations — the fibre of y = x²/(x−s) over y = ∞ is
{s, ∞} (unramified), so the loop lifts to a loop around t₂ and

  ρ(γ)H − H = (M_{t₂} − 1)H ∈ ker(M_{t₂} − 1) ∩ ⟨homogeneous⟩ = ⟨u_{t₂}⟩,

nonzero because H is singular at t₂ (its radius of convergence is exactly |t₂|: that is what
overconvergence *at t₁* means). Hence the (2,±1) preimages are harmless **iff ⟨u_{t₁}⟩ = ⟨u_{t₂}⟩**.
Computing the monodromy numerically (DOP853, rtol 10⁻¹², loops of several radii and base points):

| host | tr M_{t₁} | tr M_{t₂} | \|u_{t₁} ∧ u_{t₂}\| |
|---|---|---|---|
| CDT, x(1−x)(1−9x)y″+(1−20x+27x²)y′+3(3x−1)y = 0 | 2 ± 5·10⁻¹⁴ | 2 ± 5·10⁻¹⁴ | **0.7366** |
| Zagier D, x(1−11x−x²)y″+(1−22x−3x²)y′−(3+x)y = 0 | 2 ± 5·10⁻¹³ | 2 ± 7·10⁻¹³ | **0.3184** |

stable to 6 digits over base points 0.5, 0.3+0.2i, 1 and loop radii 0.03–0.06 / 3–7. Both local
monodromies are unipotent, each with a one-dimensional fixed line, and **the two lines are
distinct**. (Structurally this is forced: the monodromy group is the image of a congruence
group and M_{t₁}, M_{t₂} are parabolic generators at *distinct cusps*, whose fixed points on
∂**H** — hence whose fixed lines — differ.) The refinement is therefore **not available for the
pair that matters**, and the margin stays −1.44. **[verified]**

*Why the question was host-independent.* In **all three** configurations the
deepest non-principal preimages are the images of the principal one under the **same** two group
elements, the bottom rows (2, ±1), i.e. $\left(\begin{smallmatrix}1&0\\ \pm2&1\end{smallmatrix}\right)\in\Gamma_0(2)$
— CDT's at |z| = 0.401921, ours at |z| = 0.536032 (v₁) and 0.213693 (v₂). These are parabolic,
fixing the cusp τ = 0 of Γ₀(2), i.e. the loop around y = ∞; so "is the lift regular at the
(2,±1)-preimage?" is the *same* question about the *same* group element for CDT's row **C** on
Γ₀(6) and for row **D** on Γ₁(5) — both weight-one k = 2 rows on four-cusp genus-zero hosts
descended by the same involution. The same what-if on CDT's **own** host would have taken their
margin from +0.025 to **+0.652**; they won by 0.0053 and still demanded a *unique* preimage in
Lemma A.4.4 — consistent with the computation above, which says the refinement does not exist
for this pair on either host. **What remains genuinely open** is only the status of the
*other*, shallower preimages (bottom rows (4,±1), (6,±1), (2,±3), …): each is a separate
cocycle condition Δ_γ ∈ ⟨u_{t₁}⟩ for the corresponding γ ∈ Γ₀(2). Dropping every one of them
while keeping the proved (2,±1) exclusions is the relaxed problem of §6(b′).

### 6(b′) The relaxed problem: only the proved exclusions  **[verified]** (`relaxed.py` → `relaxed.txt`)

Enforcing **only** the two proved exclusions at each place and letting every other bad preimage
be slit or not, whichever is better (search over R, the lune parameter c, and the set of optional
points with |z| < T):

| place | required exclusions | best R, c, T | #slits | \|ψ_v′\| | L_v | BC_v | 14L−BC |
|---|---|---|---|---|---|---|---|
| v₁ | \|z\| = 0.536032 (×2) | 0.80, 7.5, 0.75 | 6 | 0.725673 | 7.630581 | 15.67637 | 91.15177 |
| v₂ | \|z\| = 0.213693 (×2) | 0.80, 7.5, 0.85 | 8 | 0.398697 | 2.219566 | 6.54823 | 24.52569 |

**relaxed margin = −1.4577**, and since the fully constrained design is itself feasible for the
relaxed problem, the true relaxed optimum lies in [−1.4416, −1.44]. In other words the optimiser,
given permission to keep the unproved bad preimages inside Ω, **declines**: the shallow preimages
sit near the boundary, so slitting them costs almost nothing in |ψ′(0)| and *reduces* BC. Only the
deep exclusions bind, and the deepest — the (2,±1) pair — is the one proved necessary in §6(b).

> **The verdict is therefore independent of every unproved exclusion.** −1.44 either way.

**(c) Not descending at v₂.** Using the Y(2) (λ-)uniformisation instead of Y₀(2) doubles both
log|φ′(0)| and BC while leaving τ alone — CDT's own Basic Remark (equivalently: τ halves and the
other two stay). It is the *same* bound. No gain.

**(d) Choosing s = t₁ instead of t₂.** Both are units, so both give an integral pure module; but
H is regular at t₁ and singular at t₂, so the pure module must live on **P**¹∖{0,t₂,∞}. Forced.

---

## 7. Honest ledger

**Exact (closed form / exact rational arithmetic).**
τ♭ = 191/49, I₂¹⁴(2) = 21.075, τ♯ = 27/80, τ = 16603/3920; the exact conformal radii of all
contours (products of 4r/(1+r)² and (c²−1)/(c²+1) over ℚ); the geometry of both places
(t₁t₂ = −1, Y₁ = −φ⁻¹⁵/(5√5), Y₂ = −φ¹⁵/(5√5), Y₁Y₂ = 1/125, average ceiling log 256); the
Jensen reduction of the Bost–Charles integral; the fact that exclusion holds by construction.

**Verified numerically (digit counts stated).**
ξ to 142 digits and its closed form to 10⁻¹⁴²; ζ(2)/5 as B_D's Apéry limit to 10⁻¹⁷³; the
`lindep` certificate; CDT's |ψ′(0)| = 0.6293401 and BC = 11.8449 against their 0.6292232680 and
11.845; the six bad preimages of −1/72 and the failure of their Lemma A.4.4 for the published
θᵢ (argument principle + winding numbers); BC(their published φ) = 12.1207 by two evaluators;
all bad-preimage tables; the certified contours at both places, with exclusion verified twice;
BC at N = 16384 with the orbit truncation converged; the attribution table; the inventory sweep;
the monodromy matrices at t₁ and t₂ on both hosts and the distinctness of their fixed lines
(§6(b)); the far-cusp periods −11ζ(2)/5 and ξ′ (12 digits here, 174 in `task2/`, with the
near-cusp calibration to 212 digits and an independent Abel-summation cross-check); the
non-proportionality of the v₁ and v₂ relations; branch B's margin −4.4038.

**Cited, not recomputed.**
CDT's Remark BCboundK and their convexity bounds 13.730/13.7206/13.678/13.621; Beukers'
integrality of the companions; `INVENTORY_BOUND.md`'s Theorem 2.1 (u₁ = 1) and its supply counts;
`lattice/hostscan/REPORT.md` §10 for the completeness of the (host, weight, period) census.

**Estimated / open.**
(i) The transported convexity factor 0.97336 is CDT's own relative gain, not a recomputation on
our contours. (ii) The claim that the successive-slit construction is extremal is verified on
symmetric benchmarks (𝔻 minus n equally spaced radial slits, closed form, ratio 1.0000 for
n = 2,3,4), not proved in general. The crude one-point cap alone does **not** settle the sign:
frozen shapes at the caps would give +2.7. What settles it is the *saturation* of §5.1 —
the family's own frontier — together with the exact decoupled optimisation of §4.3.
A domain outside the family "disc ∖ (one lune ∪ slits)" is not excluded by anything proved here.

(iii) The monodromy refinement is **closed for the decisive (2,±1) pair** (§6(b)) and *moot*
for the shallower ones: §6(b′) shows that granting them all changes the margin by less than
0.02, because the optimiser slits them anyway.
(iv) Unbounded-degree K(y)-independence of the fourteen — CDT's Lemma 12.1.1, a monodromy
argument — is **[assumed]**, not transported; only degree ≤ 5 is certified (§3.3).
(v) *(closed)* The factor −11 is explained in §3.5 by an exact character-sum computation
(dilogarithm parity at the two cusps); the numerical coincidence with −(t₁+t₂) is not claimed to
be more than that.
(vi) The identification of the analytic conditions (i), (ii) with fold-regularity at the two
places rests on the standard dictionary "Apéry limit = period at the cusp = log coefficient of
the Eichler integral", verified here at the near cusp to 212 digits; it is not proved.

---

## 8. Scripts and data

| file | what it does |
|---|---|
| `contour.py` | h, h′, Γ₀(2)-orbit preimage enumeration; CDT's Slit/lune maps **and their exact inverses**; CDT's published parameters |
| `wind.py` | winding numbers (independent exclusion check) |
| `region2.py` | the region class Φ = Rz ∘ T₁∘…∘T_k, successive slitting to prescribed points, exact conformal radius |
| `bcfast.py` | Bost–Charles by the Jensen reduction (with the small-argument fix) |
| `bcdirect.py` | Bost–Charles by the model-free double sum |
| `opt2.py` | designer + evaluator for one place (bad points, best slit order, L, BC, verification) |
| `places.py` | the two real places of K, exact data |
| `cdt_bound.py` | CDT's τ♭, τ♯, I_u^v(w), the bound (copied from `lattice/cdt_finder/`) |
| `calib_cdt.py` → `calib_cdt.txt` | Task 1: CDT calibration, the two Appendix-A corrections, the reconciliation |
| `xi_check.gp` → `xi_check.txt` | independent re-verification of ζ(2)/5, ξ, the closed form and the `lindep` certificate |
| `two_places.py` → `two_places.txt` | Task 4: geometry, trade-off curves, optimum, certified contours, the averaged bound |
| `verdict.py` → `verdict.txt` | cross-validation of the two BC evaluators, hard caps, attribution, inventory sweep, convexity |
| `curves.json` | the measured (R, #bad, \|ψ′\|, L, BC) curves at both places |
| `whatif.py` → `whatif.txt`, `whatif_cdt.txt` | the monodromy loophole, quantified at both places and on CDT's own host |
| `monodromy.py` → `monodromy.txt`, `monodromy_robust.txt` | the monodromy computation that CLOSES the loophole for the decisive (2,±1) pair |
| `relaxed.py` → `relaxed.txt` | the relaxed problem: only the proved (2,±1) exclusions enforced |
| `farcusp_data.gp`, `farcusp.py` → `farcusp_robust.txt` | the far-cusp periods −11ζ(2)/5 and ξ′, with the near-cusp calibration and a robustness sweep |
| `L3_identity.gp` → `L3_identity.txt` | Re L(3,ψ₄) = φ⁵ Im L(3,ψ₄) to 212 digits, and the two character sums behind −11/5 |
| `single_hypothesis.py` → `single_hypothesis.txt` | branch B: the variant needing only the single hypothesis |
| `saturation.txt` | the R → 1 saturation of the conformal radii at both places |
| `task2/` | sources, companions, conditional ODE, fold-regularity at both places (agent report) |
| `task3/` | pure module over K, measured denominator array, K(y)-independence (agent report) |

---

## 9. The final statement

**What the construction would refute, if the margin were positive.** Not

> *(T) 1, ζ(2), ξ are linearly independent over K = ℚ(√5)*

but the strictly stronger

> *(T′) there is no (a,b,c) ∈ K³∖{0} with **both** a + b·ζ(2)/5 + c·ξ = 0 **and**
> σ(a) − 11σ(b)·ζ(2)/5 + σ(c)·ξ′ = 0,  where ξ′ = −Re L(2,ψ₄) − φ⁻⁵ Im L(2,ψ₄),*

because a single K-relation supplies a conditional function at the first place only: the
far-cusp period of B_D is −11ζ(2)/5, not ζ(2)/5, so the second place's fold-regularity
condition is *not* the Galois conjugate of the hypothesis (§3.4, verified to 174 digits;
§3.5 proves the discrepancy is the dilogarithm parity of the **outer** orientation and that on
the **inner** (ξ) line equivariance does hold exactly).
Equivalently, (T′) asks for two ℚ-linear relations among 1, ζ(2), Re L, Im L, √5, √5ζ(2),
√5 Re L, √5 Im L where (T) gives one; in the sub-case b = 0 it says the rational part and the
√5 part of the single K-relation a + cξ = 0 must vanish separately.

**And the margins.**

| | hypothesis refuted | L̄ | B̄C | entry | bound | **margin** |
|---|---|---|---|---|---|---|
| **A** | (T′) — the doubled hypothesis | 4.90502 | 10.81633 | +0.66956 | 16.154 | **−1.4424** |
| **B** | (T) — the honest single hypothesis, paid for at v₂ by also excluding the principal preimage of Y₂ | 4.59778 | 9.47631 | +0.36232 | 26.154 | **−4.4038** |

Both entry tests pass; both margins are negative; m = 14 is the optimal inventory in both cases;
the verdict survives dropping every unproved exclusion (§6(b′)) and the transported convexity
improvements (§5). **The theorem does not follow.** The shortfall is **1.44 nats** on the
doubled hypothesis and **4.40 nats** on the target hypothesis, and in both cases the responsible
ingredient is the **second real place of ℚ(√5)**: it puts the removed fold at
|Y₂|/4 = 30.5 orbifold radii rather than CDT's 1/288, its deepest bad h-preimage at |z| = 0.2137
rather than 0.4019, and — separately — it asks for a relation the hypothesis does not give.

**What is left.** (a) The hypothesis tax is confined to the **outer** (ζ(2)) direction: on the
inner (ξ) line equivariance is a theorem (§3.5). A host whose *outer* direction also had uniform
cusp weights — i.e. a Lambert rather than a dilogarithm orientation — would remove the tax
entirely and put branch A's hypothesis back to (T). That is a concrete search criterion, and it
is orthogonal to the geometry. (b) The
place-v₂ geometry is a motivic invariant of the pair (host, s): |Y₂| = φ¹⁵/(5√5) is forced by
t₁t₂ = −1, so no choice of contour or inventory can repair it — a different host would be
needed, and by `lattice/hostscan/REPORT.md` §10.5 there is no other Apéry-perfect k = 2 host
carrying a new period. (c) The two corrections to CDT's Appendix A (§1.3) are worth reporting to
the authors independently of all of this.
