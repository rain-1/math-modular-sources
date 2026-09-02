# Findings log (Claude's working notes)

Recovered results, cross-checks, and other findings from mining the archive,
kept separate from Fable's ledger work on `packages/Odd_Zeta_Research_Master_Ledger.xlsx`.
Newest entries on top.

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
