# Packet K — final resolution (2026-08-27, rev. 2 after Sol's correction)

**Verdict: Section K cannot be fixed, and the PTPAH architecture as a whole cannot prove
$G \notin \mathbf Q$.** This is not a failure to find the right lemma; it is closed by
proof at both the scalar and the filtered level. This note records the closure, the new
results produced along the way, and the exact requirements any successor programme must
meet.

**Revision 2:** Sol pointed out that the first version's "Level 3" claim (K11.1
overcredited by $4\log2\,D^2$) was an error: the per-step re-derivation omitted the source
arithmetic Hilbert–Samuel degree $\widehat{\deg}\,\overline E_D = \tfrac m2(L_b\cdot L_b)D^2
= b_\infty D^2$ (for $m=2$), the term the project's own master inequality carries
(`consolidation/ADELIC_HOLONOMY.md`, $\widehat{\deg}\,\overline{E^{\rm ad}_D}$ formula and
the slopes inequality over $\mathcal V_D$). Including it restores K11.1 exactly. §3 is
rewritten accordingly; the honest endpoint deficit is $-0.6313$ (pinned folds), not
$-3.40$. The verdict is unchanged — the corrected sub-bundle ledger (§3.1) is strictly
negative on $(0,1]$ and exactly critical at $x=0$.

Companion documents:
`PACKET_K_RESOLUTION_RESEARCH_HANDOFF_2026-08-27.txt` (the prior audit),
`k_polygon_report.md` + `k_polygon.py` + `k_polygon_output.txt` (the exact computation,
21/21 verification checks, polygons to $D=34$).

---

## 1. Level 1 — the scalar ledger (recap; closed by conservation)

For the fixed Hermite evaluation map, $v_2(\det\psi_D)=D$ exactly (Packet E). Enlarging
the source lattice by the E/Z saturation of index $2^{\Omega_D}$ forces
$$v_2(\det(\psi_D S_D)) = D-\Omega_D$$
identically. Top-determinant bookkeeping of a fixed rational isomorphism is a zero-sum
game under the product formula, so **no scalar rebuild of K can yield $C_{2,\rm eff}>0$**;
the corrected value is $C_{2,\rm eff}=0$ (up to $O(D)/D^2$). The withdrawn claim
$C_{2,\rm eff}=\log 2$ was exactly the size of the uncharged evaluation mass
$\Omega_D - D - 3R_D = D^2 + O(D\log D)$.

Numerically: the Section L contradiction margin was $+0.02165$; with $C_{2,\rm eff}=0$ the
inequality holds with slack $-0.6715$ (or $-0.6313$ after the fold-pinning credit of §3.3).
The hole is $\approx 30\times$ the margin.

## 2. Level 2 — the filtered/polygonal programme (closed by the rigidity theorem)

The handoff's §XIII/§XIX programme asked for the full 2-adic Smith polygon of the
E/Z-adapted map, hoping for interior structure ("misaligned flags") that the endpoint
conceals. The computation answers this, and the answer is a theorem
(`k_polygon_report.md` §5.1):

**Theorem (F1 rigidity).** In the paper's own K4/K5 formulation,
$$N_D = M_D\,\mathcal C_D = 2\,V_D W\,\operatorname{diag}(\Delta_j^{-1}),
\qquad V_D W \in GL_D(\mathbf Z_2),$$
hence the Smith exponents are exactly
$$e_{D,j} = 1 - v_2(\Delta_j) = 1-8j+4s_2(j), \qquad j=1,\dots,D.$$

*Proof.* K4.7 gives $\mathcal C_D\operatorname{diag}(\Delta_j)\in GL_D(\mathbf Z_2)$, so
$\mathcal C_D = W\operatorname{diag}(\Delta_j^{-1})$; E2+E3 (=K5.1) give $M_D = 2V_D$ with
$V_D\in GL_D(\mathbf Z_2)$; Smith exponents are invariant under $GL_D(\mathbf Z_2)$ on
either side. $\square$

Consequences:

1. **No misalignment.** The E/Z filtration and the Hermite filtration are simultaneously
   diagonalisable over $\mathbf Z_2$. The Hankel mass $D$ is spread exactly one unit per
   direction. The polygon is the Packet-A depth sequence translated by $+1$ — it contains
   *no information beyond Packets A and E*.
2. **No two-sided volume.** Positive Smith mass is exactly $0$ (F1) and $\Theta(\log D)$
   in the only mixing formulation found (F2, the reversed-Padé functionals). Contrast
   Packet C, the programme's model case, with both one-sided masses $\approx D^2$. A
   flag-misalignment contradiction needs $\Theta(D^2)$ of interior freedom; there is
   $O(D\log D)$ at most.
3. **Sub-bundle ledger closed.** Writing the slopes inequality for every flag member
   $E'_x = (\text{base block}) \oplus (\text{first } xD \text{ fold directions})$, all
   terms are explicit except the 2-adic net density; the rigidity theorem forces that
   density to $0$ uniformly in $x$ (gain $\sum_{j\le r} v_2(\Delta_j)$ and adverse
   evaluation mass $\sum_{j\le r} (v_2(\Delta_j)-1)$ cancel direction-by-direction). The
   honestly derived gap function $\Gamma(x)$ (§3.1, including the partial source
   Hilbert–Samuel degree) satisfies $\Gamma(0)=0$ (the base block is exactly
   semistable-critical) and $\Gamma(x)<0$ for all $0<x\le 1$, with interior minimum
   $-0.8734$ at $x^*\approx0.687$ and endpoint $\Gamma(1) = -0.6313$. No intermediate
   semistability violation exists.

### 2.1 The Archimedean side cannot make up the difference

Two facts, both new to this resolution:

- **Fold pinning (+0.0402, provable).** J5 minimised over all assignments of the $D$ fold
  directions to pivot orders ($I_*\approx 2.3285$). But the polynomial block has unit
  leading coefficients, so base occupies pivots $0..D-1$ and folds are pinned to
  $D..2D-1$: the correct constant is
  $\int_1^2 g(s)\,ds = 3\log 2 + (\Delta-L)^2/2L \approx 2.3687$.
- **Assignment cap ($\le 0.296$).** The maximum conceivable fold-placement recovery, even
  if folds could be forced to the worst-case positions, is
  $\max_{|F|=1}\int_F g - I_* = \bigl(2d + \Delta^2/2L\bigr) - 2I_* + I_*
  \approx 4.9530 - 2.3285 - 2.3285 \approx 0.296$.

Since $0.296 < 0.6313$, **no pairing of the (rigid) finite polygon against the
base/fold flag can close the deficit**, even granting the adversary's worst case. This
seals the specific pairing route the handoff §XIII proposed.

## 3. Level 3 — the assembly, corrected: K11.1 is derivable (Sol's correction)

Revision 1 of this note claimed K11.1 was overcredited by $4\log2\,D^2$. **That claim is
withdrawn.** The per-step re-derivation omitted the source arithmetic Hilbert–Samuel
degree: for the rank-$2D$ source ($m=2$ channels of polynomial slots $0..D{-}1$, all
normalized in the base metric $\overline L_b$ with $(\overline L_b\cdot\overline L_b) =
\log 16$),
$$\widehat{\deg}\,\overline E_D
 = \frac m2(L_b\cdot L_b)D^2 + o(D^2)
 = b_\infty D^2 = 4\log2\,D^2 .$$
This is the same $\frac m2\widehat T(1,\varphi)D^2$ term the project's proven master
inequality carries on the left (`consolidation/ADELIC_HOLONOMY.md`). The evaluation
heights cancel all-base as before, but this degree survives, and the scalar inequality
becomes $b_\infty + C_2 - 3 \le \Delta - I_*$, which is exactly K11.1
($C_2 + 2b_\infty - 3 + I_* \le T(1)$, using $b_\infty = \Delta = 4\log2$,
$T(1) = 2b_\infty$). The honest endpoint deficit is therefore
$$\Gamma(1) = I_*^{\rm pinned} - 3 = 3\log2 - 3 + \frac{(\Delta-L)^2}{2L}
 = -0.63130\ \ (\text{or } -0.6715 \text{ with the unpinned } I_*).$$

### 3.1 The corrected sub-bundle gap $\Gamma(x)$, all terms derived

For $E'_x$ = base block $\oplus$ first $xD$ fold directions (pivots $0..D{-}1$ base,
$D..(1{+}x)D{-}1$ fold), with $r = xD$:

- source degree: $\widehat{\deg} = 2\log2\,(1+x^2)D^2$
  (base slots $\sum_{j<D} j\log16$, fold slots $\sum_{j<r} j\log16$);
- odd denominators: $-(2x+x^2)D^2$;
- $p=2$: $0$ up to $O(D\log D)$, by the rigidity theorem, uniformly in $x$;
- heights: all-base part $2\log2\,(1-x^2)D^2$, fold conversion
  $\bigl[-2\log2(x+\tfrac{x^2}2) + 4\log2\,x - A(x)\bigr]D^2$, where
  $A(x)=\int_1^{1+x}(\Delta-Ls)_+\,ds$.

Collecting terms ($d=2\log2$, kink at $x_c = \Delta/L - 1 \approx 0.5729$):
$$\Gamma(x)/D^2 =
\begin{cases}
(5\log2 - 1 - \tfrac L2)\,x^2 - (2\log2 + 2 - \Delta + L)\,x
 = 1.5843\,x^2 - 2.3764\,x, & 0\le x\le x_c,\\[2pt]
(5\log2 - 1)\,x^2 - (2\log2+2)\,x + \tfrac{(\Delta-L)^2}{2L}
 = 2.4657\,x^2 - 3.3863\,x + 0.2893, & x_c\le x\le 1.
\end{cases}$$

Properties (verified numerically to machine precision):

- $\Gamma(0) = 0$ **exactly**: the base block alone saturates the slopes inequality — it
  is semistable-critical, as it must be for an unconditional construction. This is a
  strong consistency check on the whole corrected ledger.
- $\Gamma$ strictly decreases from $0$, reaches its minimum $-0.87337$ at
  $x^* = \tfrac{2\log2+2}{2(5\log2-1)} \approx 0.6867$, and recovers to
  $\Gamma(1) = -0.63130$, which coincides exactly with the pinned-fold scalar deficit.
- $\Gamma(x) < 0$ for all $0 < x \le 1$: **no sub-bundle of the flag certifies a
  contradiction** with the currently proven height bounds (I.4/I.5, two radii). Since the
  bounds enter with a minus sign, sharper Archimedean analysis could only raise $\Gamma$;
  the assignment cap of §2.1 limits that headroom to $\le 0.296 < 0.6313$.

Deficit summary (honest assembly, $C_2 = 0$):

| quantity | value per $D^2$ |
|---|---|
| endpoint deficit, unpinned $I_*$ (= K11 corrected) | $-0.6715$ |
| endpoint deficit, pinned folds | $-0.6313$ |
| worst sub-bundle ($x^* \approx 0.687$) | $-0.8734$ |
| max conceivable Archimedean recovery (cap, §2.1) | $+0.296$ |

Every configuration leaves the ledger strictly negative.

## 4. What survives (keep these)

1. **F1 rigidity theorem** (§2) — new, proved, and the sharp closure of §XIII/§XIX.
2. **Exact valuation $v_2(b_n) = 2s_2(n)-1$** for all $n$ (verified to $n=60$) —
   strengthens Packet E's parity theorem E2 from a mod-2 statement to the exact valuation.
3. **Fold pinning** (§2.1) — the flag structure determines the fold pivot set; J's
   distribution-free minimum is not attained by the actual geometry ($+0.0402$).
4. **F2 mixing laws** (conjectural, verified to $D=34$): maximum Smith defect of the
   moving reversed-Padé functionals steps by $+2$ exactly at $D = 2^k-1$ and
   $3\cdot 2^{k-1}-1$ (not at $5$), and the unique positive exponent is
   $e_{D,D} = 2\lfloor\log_2(D+1)\rfloor - 5$. Genuine binary structure, at
   $\Theta(\log D)$ scale.
5. **The sub-bundle ledger formalism** $\Gamma(x)$ — reusable for any successor
   architecture.
6. All items in handoff §XVII (Packets A, B, C, D, E, G exact results; H–J analytic
   results; the fixed-fibre bridge; the contiguous operators).

## 5. Requirements for any successor programme

Any future attempt on $G$ via this style of argument must produce **all three** of:

1. **A non-rigid evaluation pair.** By the rigidity theorem, an architecture can carry a
   $D^2$-scale filtered surplus at $p=2$ only if its evaluation map is *not* of the form
   $c\cdot(\text{unit})\cdot\mathrm{diag}$ relative to its adapted lattices. This is a
   sharp, checkable test — run it before investing in any candidate.
2. **A genuine "$\tau$-resource":** rational subspaces whose evaluation images are deeply
   flag-tangent 2-adically *at low height* — divisibility not paid for coefficient-wise
   through the product formula. Nothing of the kind exists in Packets A–L; the F2 mixing
   ($\Theta(\log D)$) is negative evidence at the natural candidates.
3. **A stated global filtered theorem in the manuscript:** the assembly is now known to
   be derivable (§3, via the Hilbert–Samuel source degree of `ADELIC_HOLONOMY.md`), but
   the manuscript itself should state the slopes inequality with defined lattices,
   metrics, and degree normalizations, so that the ledger is a corollary rather than a
   verbal assembly — revision 1 of this note mis-derived it precisely because the
   manuscript leaves the convention implicit.

## 6. Recommended status statement for the manuscript

> The local Padé, Christoffel–Darboux, reflection, Hermite, odd-prime and Archimedean
> calculations are substantial and several are proved exactly. The adelic assembly
> (Sections K-Finite, K-Global, L) is withdrawn: the $p=2$ source saturation of size
> $\Omega_D$ reappears exactly in the determinant of the fixed Hermite evaluation map
> (conservation), the 2-adic Smith polygon of the E/Z-adapted map is rigid —
> $e_{D,j} = 1 - v_2(\Delta_j)$, with the E/Z and Hermite filtrations simultaneously
> diagonalisable over $\mathbf Z_2$ — so no filtered/polygonal refinement of the ledger
> can recover the deficit: the honestly derived sub-bundle gap satisfies $\Gamma(0)=0$,
> $\Gamma(x)<0$ for $0<x\le1$, with endpoint $-0.63$. No irrationality proof for
> Catalan's constant is established, and the present architecture cannot yield one.

## 7. Reproduction

```
cd /home/ubuntu/code/math-modular-sources/catalan
python3 k_polygon.py     # ~30 s: 21 verification checks + all polygon tables
```

Constants used throughout: $d=2\log2$, $\Delta=4\log2$, $L=\operatorname{arcosh}3$,
$T(1)=8\log2$, $b_\infty=\log16$, $I_*=\log2+6(\log2)^2/L\approx2.3285$,
$\int_1^2 g \approx 2.3687$, $\int_0^2 g = 2d+\Delta^2/2L \approx 4.9530$.
