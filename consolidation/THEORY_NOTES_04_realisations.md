# Theory notes 04 — Apéry limits as realisations of one extension class

*Fable, 2026-08-22 (early). State after Theorems A, B (verified), C, D (C,F at 3), F (Euler criterion, census 14/14). Tags as before.*

## 1. The picture

A modular Apéry row is an **extension** $0\to\mathcal V\to\mathcal E\to\mathbb Q(0)\to0$ of the pure elliptic local system $\mathcal V=\operatorname{Sym}^wH$ by the trivial object, given by the Eisenstein (or cuspidal) source $\Phi=P(V)E$. Its Apéry limits are the *periods of that extension class at the places of $\mathbb Q$*:

| place | what the limit is | mechanism | status |
|---|---|---|---|
| $\infty$ | $\xi_\infty=L(\Phi,w+1)=-\tfrac12P(w+1)\,L(\psi,w+1)$ — the critical value, i.e. the Beilinson regulator of the Eisenstein class | fold lemma: value of the Eichler integral at the Atkin–Lehner fixed point; period polynomial | Theorem B, prefactor verified case by case; uniform proof in progress |
| slope prime $p$ ($p\mid$ level, Euler factor divides $P$) | $\xi_p=-\tfrac12Q(w+1)\,L_p(w+1,\psi\omega^{-w})$ — the Kubota–Leopoldt value, i.e. the $p$-adic (syntomic) regulator of the same class | Calegari's overconvergence: $\Theta$ is a finite oldform combination of the $p$-adic Eisenstein series; $F(\Theta-\xi)$ overconvergent iff $\xi$ = constant term | Theorem F; 14 identifications to ~3000 digits |
| cuspidal source, slope prime | $0$ (no constant term) | Coleman's $\theta^{w+1}$ | Prop. (cuspidal clause); verified on the $L(f,2)$ row |
| outer placement ($\varphi\ne\mathbf 1$), slope prime | $0$ | constant term of $E^{\mathbf 1,\chi}$ at $\infty$ vanishes | verified: A at 2, ζ at 3 |
| good prime $p\nmid N$ | no limit along all $n$; tower limits $\lim_s b_{ap^s}/a_{ap^s}$ exist | **conjecture:** these are periods of the *pure* system — Dwork unit-root Frobenius data — not of the extension; hence Γ_p-values and no rigidity | under test (GOOD_PRIME_TOWERS) |

The archimedean and $p$-adic rational factors differ by exactly the Euler factor at $p$: $r_p=r_\infty/(1-\psi(p)p^{-(w+1)})$ — the standard Euler-factor removal of $p$-adic $L$-functions, now visible in elementary data. Rigidity (Conjecture D) is the statement that two rows with the same $E$ have limits in the same ratio at every place where both have limits; it is a corollary of the table.

## 2. Why the slope/product formula is true (now a theorem-shaped explanation)
$\sigma_p>0$ iff the Euler factor divides $P$ iff (empirically) $p\mid c=\lambda_1\lambda_2$. The equality $\sigma_p=v_p(c)$ is the Casoratian; the *existence* of a limit is Theorem F. The exact $p$-adic radius of $R_\xi$ is $p^{v_p(c)}=|t_1t_2|_p^{-1}$ (radius of $A$ being 1): the conditional function continues $p$-adically to the product of the singular values, as it continues archimedeanly to the second singular value. (E.g. C: $A$ has radius $1/9=|t_1|_\infty$ and $R$ radius $1=|t_2|_\infty$ at $\infty$; $3$-adically $A$ has radius $1$ and $R$ radius $9=|t_1|_3$.) A geometric proof of the exact radius (overconvergence up to the supersingular annuli, Buzzard-style as in Calegari's $2^{12}$) would complete Proposition C from the modular side. [open, minor]

## 3. Consequences for the Diophantine programme
- The common-period lemma of the Catalan two-row method — Zudilin's row and the level-8 row converge $2$-adically to $\zeta_2(2)$ and $\tfrac12\zeta_2(2)$ — is now Theorem F for the modular row plus (Beukers' $p$-adic Padé theory, in progress) for the hypergeometric row. The Nesterenko case is in Lean.
- Which pairs align is decidable from the source table: same $E$, Euler factor dividing both $P$'s at $p$. New aligned pairs predicted and confirmed: Zudilin×E (factor 2), Domb×ε at 2 (ratio $4/3$), B,C,F,$s_{18}$ at 3. ζ(5): the level-16 projector row has $(1-X)\mid P$, so it has a $2$-adic limit $-Q(5)\cdot\tfrac12\zeta_2(5)$; whether Brown–Zudilin's row aligns with it is being tested.
- Nonvanishing past the threshold remains exactly Lemma P2 (irrationality-strength); nothing above changes that, but the realisation picture is the right home for the "adelic Padé" problem: a construction whose $p$-adic contact is independent of its archimedean contact must come from a different extension-class realisation, not from the moment family.

## 4. Open theoretical items, ranked
1. Uniform archimedean Theorem B (period polynomial at the Fricke point) — agent running.
2. Hypotheses (b),(d) of Theorem F; row B descent — agent running.
3. Hypergeometric rows: an Eisenstein-type description at the $p$-adic place (Beukers) — agent running.
4. Good-prime dichotomy — agent running.
5. Syntomic interpretation: identify $\xi_p$ with the syntomic regulator of the Beilinson–Kings Eisenstein class (Bannai–Kings for Dirichlet characters); this would make the table a theorem about one motive rather than a coincidence of two computations. [conceptual; literature]
6. Vector/higher-rank systems (β(4) rank-16 direct image; $L(3,\chi_5)$ Hadamard): Theorem F should extend verbatim (overconvergence is about the source, not the rank); the census of slopes for these systems has not been done.
