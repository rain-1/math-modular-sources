# Independent audit: the "Catalan irrationality" two-row lattice claim

**Verdict: the claim is FALSE. Catalan's constant is not proved irrational here.**
The rows and the CSVs are exactly correct; the arithmetic in `lattice/catalan_bridge/`
is also correct (I reproduced `log|q_n G - p_n|/n ≈ -0.8` from scratch in PARI with my own
exactly-constructed rows). The failure is **logical, not arithmetical**: in the regime the
other agent entered (`F_n < 0`) the construction degenerates into Dirichlet's theorem, which
produces such pairs for *every* real number, rational ones included, and the step that would
make it a proof — an unconditional lower bound on `|q_n G - p_n|`'s non-vanishing together
with the size bound — is unavailable.

## 1. Exact rows (scripts: `lattice/catalan_audit/rows.gp`, `nest2.gp`)

*Zudilin.* Built `Q_m,P_m` as exact rationals from `Q_0=1,Q_1=7/4,P_0=0,P_1=13/8` and the
stated second-order recurrence, m ≤ 40. Verified: `P_m/Q_m → G` (200+ digits, `Catalan` in
PARI); denominators of `Q_m` are pure powers of 2; `2^{e_m}Q_m ∈ Z` and
`2^{e_m}D_{2m-1}^2 P_m ∈ Z` for all m ≤ 40 with `e_m = min(6m, 4m+3+⌊log2(2m-1)⌋)`;
`log|Q_m|/m` and `log|Q_m G-P_m|/m` drift to ±5 log φ = ±2.4061 (2.2698 / −2.5266 at m=40).

*Nesterenko (4,7).* Built `A_{1,n,j}, A_{2,n,j}` by solving the partial-fraction linear
system for `R_n(t)` exactly over Q (no numerics). Independent confirmation: my `A_{2,n,j}`
equal the paper's closed form (2.11) for every j, n ≤ 4. Then `B_n`, `C_n`,
`V_n = 4^{7n+1}D_{6n}^2 B_n`, `U_n = 4^{7n}D_{6n}^2 C_n`. Verified: `V_n, U_n ∈ Z`;
`J_n = 4B_n G − C_n > 0` and `→ 0`; `log J_n/n → 2 log T^- = −7.3109` (−7.557 at n=10);
`log B_n/n → 2 log T = 7.6507`; `log|V_nG−U_n|/n → E_2 = 14.3931`.

## 2. CSV comparison

`zudilin_rows.csv` (n ≤ 13 checked) and `nesterenko_rows.csv` (n ≤ 10 checked) agree with my
exact rows **digit-for-digit, exactly**. The "mpmath numerical evaluation" suspicion is
**refuted**. `combined_rows.csv` is also fine and *not* the source of the claim: its pairs
give `log|q_nG−p_n|/n ≈ +8.17`, `log|q_n|/n ≈ 23.1`, δ ≈ 0.647 — matching the LCM paper's
0.6589. Nothing anomalous there.

## 3. My own lattice (script: `lattice/catalan_audit/lattice_build.log`, `fresh.gp`)

Using only my exact rows, with `S_n` = the literal prime-window product (3.6) and separately
`S_n = D_{6n}^2`, `L_n = {c : c_1Y_n + c_2U_n ≡ 0 mod S_n}`, index computed directly from an
explicit basis (verified `det = S_n/gcd(S_n,Y_n,U_n)`), short vectors by 2-D Lagrange/LLL and
by the prescribed anisotropic box:

| n | modulus | log|q|/n | log|qG−p|/n |
|---|---|---|---|
| 10 | window | 24.79 | **+10.70** |
| 10 | `D_{6n}^2` | 21.67 | **+7.03** |

`q_nG − p_n` **diverges**, exactly as the papers predict (`F ≈ +8.2`, δ ≈ 0.66 < 1).
This is the honest construction, and it gives no irrationality.

## 4. The bridge run, and why it fails

`lattice/catalan_bridge/boxtest.py` differs by taking `T_n = 2^{v_2(h_n)}` — the **exact
per-n valuation** of the mixed minor `h_n = a_1U_n − a_2Y_n` (`a_r = X_r/S_n`) — and dividing
by `M_n = S_n T_n`. I checked everything sceptically and it all holds up:

* `v_2(h_n)` really is ≈ 24.06 n bits (trivial bound is only ~37 bits); the 2-adic bridge is real.
* `[Z^2 : K_n] ≤ M_n` — I computed it directly via `matkerint`/`mathnf`, not from the formula
  (e.g. n=8: covol = M/11). So the "correlated" index claim is genuinely true.
* The outputs really are integers: verified `q_nM_n = c_1X_n+c_2V_n` in exact integer arithmetic.
* `|q_nG−p_n|` really is `e^{−0.8n}`, confirmed in PARI at 3000 digits from my own exact rows
  (n = 6..13) and from the CSVs (n = 20,40,60,98). **Not a Python bug, not a bad CSV.**

So why is it not a proof?

**(a) The control experiment.** Let `G* = bestappr(G, 10^320)` — a *rational* number with a
320-digit denominator, `|G−G*| = e^{−1472}`. The lattice `K_n` never uses `G`, so the same
`c_n` is produced. At n = 98: `|q_nG − p_n| = e^{−78.051}` and
`|q_nG* − p_n| = e^{−78.050}`, non-zero. The identical numerics are produced by a number that
**is rational**. No contradiction arises because `1/den(G*) = e^{−735.7} ≪ e^{−78}`. The
observation "|ℓ_n| is exponentially small and non-zero" therefore proves nothing at any finite n.

**(b) The pairs are worse than trivial.** A continued-fraction convergent of `G` with
`log|q| = 812` has `|qG−p| = e^{−813}`. The bridge's pair with `log|q| = 1388` has only
`e^{−78}`. Every real number, rational or not, has infinitely many pairs with `|qα−p| < 1/q`.
δ > 1 for individual pairs is not an irrationality criterion.

**(c) The precise broken step.** The papers' Theorem 5.1 hypothesis "`H_n` and `F_n` eventually
bounded below by a **positive** constant" is *violated* here: the agent's `F_n = (E_1+E_2−σ_n)/2`
is ≈ −0.77 because `σ_n = log(M_n)/n ≈ 28.7` exceeds `E_1+E_2 ≈ 27.5`. In the papers' proof
`ℓ_n ≠ 0` is obtained from the output-determinant identity, which only ever says the two
minima's forms are *not both* zero; when the first minimum has `ℓ = 0` (precisely the case if
`G` were rational — the "rationality kernel vector" the papers warn about) the proof switches
to the second vector, which is non-zero but carries **no** `e^{F_n n}` upper bound. So the
construction can deliver "`ℓ_n ≠ 0`" or "`|ℓ_n| ≤ e^{F_n n}`", never both provably. Taking the
Lagrange-shortest vector and *observing* `ℓ_n ≠ 0` closes that gap only numerically.
Equivalently: proving `c_1Λ_1 + c_2Λ_2 ≠ 0` requires knowing `Λ_2/Λ_1 ∉ Q`, i.e. assuming what
is to be proved. **Circular.**

**(d) Secondary, independent gap.** `T_n = 2^{v_2(h_n)}` is the exact per-n valuation, with no
proved lower bound. The 5:8 paper's 2-adic roof (its eq. 4.6/4.8) is proved for Zudilin at 5n
against the *level-8 modular row* at 8n, **not** for Zudilin-at-3n against Nesterenko-(4,7)-at-n.
So `σ_n → 28.7` is itself only empirical; `v_2(h_n)/n` is measurably *decreasing*
(30.5 bits at n=2 → 24.32 at n=98, fitting `24.06 + O(1/n)`), and the whole sign of `F_∞`
turns on whether it stays above the threshold 22.35 bits/n. Nothing proves it does.

> **Update (`P2_HOLONOMIC.md`, 2026-08-23).** This gap is now localised. With
> `xi2 = zeta_2(2)` the common 2-adic limit of both rows, the following hold exactly for
> `4 <= n <= 200` (exact integer arithmetic, xi2 at 7011 2-adic bits):
> `v_2(X_n xi2 - Y_n) = v_2(X_n) + 24n - 1 - 4 s_2(3n)` (**proved**, `ZUDILIN_2ADIC.md`),
> `v_2(V_n xi2 - U_n) = v_2(V_n) + 28n - 2 s_2(n) - 2 s_2(3n) - 1` (**verified**, new), hence
> `v_2(X_n U_n - V_n Y_n) = v_2(X_n) + v_2(V_n) + 24n - 1 - 4 s_2(3n)` and
> `v_2(h_n) = 24n - O(log n)`, explaining the observed decrease exactly. The 2-adic roof is
> therefore **24**, above the threshold 22.35, and what remains to be proved is one clean
> statement about the Nesterenko row alone: `v_2(zeta_2(2) - C_n/(4 B_n)) >= 24n + O(log n)`.
> This does not touch gaps (a)-(c); the irrationality inference remains blocked at (c).

**(e) Modulus double-use.** The suspicion that `M_n` was used both as covolume and division
modulus is, unusually, *not* the error — `[Z^2:K_n] ≈ M_n` is genuinely true here. But the
consequence is that the box area equals the covolume, so Minkowski's first theorem in that box
**is** Dirichlet's theorem. Any construction pushed to `σ_n > E_1+E_2` lands in that regime and
loses all discriminating power.

## Bottom line

The exact rows, the CSVs, the 2-adic divisor, the congruence lattice and the integrality are
all correct. The irrationality inference is not: at `F_n < 0` the method reduces to Dirichlet
and yields the same output for a rational number (demonstrated), and the non-vanishing of
`q_nG − p_n` cannot be established without assuming `G ∉ Q`. Catalan's constant remains open;
the papers' own δ = 0.6589 (LCM) / 0.9025 (5:8), both < 1, are the correct conclusions.

Scripts: `lattice/catalan_audit/{rows.gp, nest.gp, nest2.gp, fresh.gp, check20.gp}` and
exact row dumps `zud_exact.txt`, `nest_exact.txt`.
