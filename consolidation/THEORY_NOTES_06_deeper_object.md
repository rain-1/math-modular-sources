# Theory notes 06 — the object behind the modular forms

*Fable, 2026-08-22. Written for River's question: "what is the deeper object, and what are its applications and implications?" Exposition first, then consequences, then what it says we should do. Tags as before; this note is mostly [lit] + synthesis, with our results as the worked instances.*

## 1. The object, in one paragraph

Every Apéry system we have met — modular, hypergeometric, Cooper's, the square roots, Zudilin's Catalan row, Brown–Zudilin's $\zeta(5)$ row — is the coefficient expansion of a **period of a mixed motive over $\mathbb Q$ with a maximally unipotent degeneration**. Concretely: a one-parameter family of varieties over $\mathbb P^1\setminus S$ (elliptic curves, Calabi–Yau threefolds, the moduli spaces $\mathcal M_{0,n}$, …) gives a **pure** local system $\mathcal V$ of rank $r$ (its cohomology, with Hodge filtration, Gauss–Manin connection, and Frobenius at every good prime); an **extension**
$$0\to\mathcal V\to\mathcal E\to\mathbb Q(0)\to0$$
in the category of mixed motives (or, concretely, of variations of mixed Hodge structure / $F$-isocrystals) is specified by a class $[\nu]\in\operatorname{Ext}^1(\mathbb Q(0),\mathcal V)$; and an Apéry system is a choice of **realisation**: an integral coordinate $t$ at the degenerate point and a de Rham section, in which the periods of $\mathcal V$ become the integers $a_n$ and the period of $\mathcal E$ becomes the rationals $b_n$. The number $\xi$ is the **regulator** of $[\nu]$ — its period at the archimedean place. Modular forms are the case where the family is the universal elliptic curve over a modular curve and $[\nu]$ is Beilinson's Eisenstein class; that case is fully computable, which is why the paper is written in its language. But the structure — pure part, extension class, realisation — is what carries the theorems, and it exists for every row.

(The GPT-era "three layers" were a correct intuition about this: Layer A = $\mathcal V$, Layer B = $[\nu]$, Layer C = the realisation. What was fanciful was treating Layer C's adelic bookkeeping as a theorem; what is real is that the *periods of $[\nu]$ at all places* are theorems.)

## 2. Why it is the right object: each of our results is a statement about one of the three pieces

| result | which piece it is about | the motivic statement |
|---|---|---|
| $a_n\in\mathbb Z$; Lucas/Dwork congruences; unit roots; $\rho^a_s\to1$ along towers | $\mathcal V$, its Frobenius | Dwork crystal structure of a family with a MUM point (Beukers–Vlasenko, *Dwork crystals I–III*): integrality and congruences of the period coefficients are the existence of the Frobenius structure, for *any* family, modular or not |
| $\lambda_1,\lambda_2$; the two fold geometries; the scalar barrier | $\mathcal V$, its monodromy | the singular points are where the family degenerates; $\lambda_2$ is the distance to the next degeneration; "Apéry-perfect" $=$ the two conifold points are exchanged by an involution of the base (Fricke) |
| $k=w+1$; free integrations; roots | the **Hodge depth** of $\mathcal E$ | each integration is one step of the Hodge filtration of the extension; a free integration is a partial splitting; taking the $w$-th root replaces $\operatorname{Sym}^wH$ by $H$ itself — the *same elliptic family*, a smaller motive, one integration |
| Theorem B\*: $\xi_\infty=L(\Phi,w+1)$ | $[\nu]$, archimedean | Beilinson regulator of the Eisenstein class, degenerated at the cusp (Huber–Kings), $=$ Dirichlet $L$-value (Beilinson–Deninger) |
| Theorem F: $\xi_p=$ Kubota–Leopoldt | $[\nu]$, $p$-adic | syntomic regulator of the same class (Huber–Kings 2003, Bannai–Kings 2010) |
| rigidity (Conj. D); alignment primes | $[\nu]$ | a regulator is an invariant of the class, functorial under Hecke operators; two realisations of one class have proportional regulators at every place — this is *why* the two-row method works at all |
| Zudilin's row $\to\zeta_2(2)$; Nesterenko's $\to$ the same | $[\nu]$, different realisation | the hypergeometric rows realise the *same* Eisenstein class in the cellular/Tate world; their $p$-adic limits are forced to agree with the modular row's |
| cuspidal $\Rightarrow\xi_p=0$; outer $\Rightarrow\xi_p=0$; weight drop | $[\nu]$ | the class has no constant term at the cusp, or is an Eisenstein class of lower weight |
| tower trichotomy at good primes | $\mathcal E$'s Frobenius at the cusp | the extension splits $p$-adically iff the Frobenius eigenvalues on $\mathbb Q(0)$ and on the Eisenstein line differ; the coefficient ratios see the splitting, not the regulator |
| Lemma P2; absorption; the rational-$G^*$ control | the realisation | nothing about $\mathcal V$ or $[\nu]$: a statement about a lattice built from *one* realisation, which is why no amount of motivic structure touches it |

Three things in this table are, to my mind, the real content of the week: that the *same* class has a regulator at every place and the Apéry data computes it where the geometry allows (§3 of Notes 05); that the root phenomenon is "replace $\operatorname{Sym}^wH$ by $H$"; and that the two-row method is regulator functoriality.

## 3. What the deeper object predicts — before any computation

Given only the family and the class, the following are determined, and we now know how to read them off:

1. **Arithmetic of $a_n$** (integrality, congruences, $p$-adic unit roots) from the Frobenius structure of $\mathcal V$ — Beukers–Vlasenko's theory does this for any family with a Laurent-polynomial model, which includes all the sporadic rows.
2. **Archimedean rates** $\lambda_1,\lambda_2$ from the singular fibres; whether there is an Apéry-perfect involution.
3. **Hodge depth** $k$, hence the denominator exponent; and whether a root/free integration exists (is the pure part a symmetric power of something smaller? is the class pulled back from a lower weight?).
4. **The period $\xi$** as an $L$-value (Eisenstein class) or as a cusp-form $L$-value or as a non-congruence period, from the type of the class.
5. **Which primes carry $p$-adic limits** — the bad primes of the host at which the oldform polynomial cancels the Euler factor — and their values. This is a motivic invariant; no change of coordinate creates an alignment prime.
6. **Whether two rows align**, i.e. whether they realise the same class: decidable from the source table, and the exact ratio is the Euler-factor-corrected archimedean ratio.

The score $\log(1/\lambda_2)-k$ is then a number attached to (family, class, realisation), and the design rule of §5 is a statement about pairs of realisations of one class. This is what I meant by "investigate at high level before doing the computations": items 1–6 are now genuinely predictable.

## 4. Implications

**(a) The classification question is a motivic one.** "Which numbers have Apéry systems?" becomes "which extension classes of which families admit an integral MUM realisation?" In weight 1 the modular answer is complete (Beauville's six surfaces $=$ Zagier's six rows); our scans say nothing new appears at genus zero up to level 104 in weights $\le3$; Brown's cellular integrals classify the $\zeta$-value side ($\mathcal M_{0,n}$); the CY operators of Almkvist–van Straten–Zudilin are the rank-4 census. A *uniform* classification would be a theorem about integral MUM variations — this is open and, I think, the correct long-term target.

**(b) Non-congruence motives are in play.** $\sqrt{F_{\rm Apéry}}$ is the period of the *elliptic family itself* written in a non-congruence way, and Beukers' Theorem 3 is the irrationality of a critical value of a non-congruence weight-three form. Scholl attached Galois representations (motives) to non-congruence forms; so this is an irrationality result for a critical $L$-value of a Scholl motive — to my knowledge the only such result, and 1987 is early for it. The root theorems say non-congruence realisations cost exactly $\log\lambda$ in score, and CDT's unbounded-denominators theorem is the statement that $\lambda>1$ always. A scan of genus-zero *non-congruence* groups with integral-after-scaling Hauptmoduls, scored by $\log(1/\lambda_2)-k-\log\lambda$, is a well-posed search no one has done; the known positive instance is Beukers'.

**(c) Irrationality proofs: what the object says can and cannot work.**
- *Single rows*: the score is a monodromy/Hodge invariant; for symmetric powers it becomes negative quickly (the scalar barrier — a statement about conifold distances of $\operatorname{Sym}^wH$ against the Hodge depth $w+1$). Roots help by exactly $(w-1)-\log\lambda$ and nothing else.
- *Two rows*: work iff they realise the same class with a common slope prime and one of them decays. The decayer must have $p$-power denominators in $a_n$ ($\kappa_p>0$) — i.e. must come from a realisation whose integrality mechanism is factorial (cellular/hypergeometric), not $q$-expansion (modular). So "cross-world pairs" is a motivic statement: *one class, two worlds*.
- *The threshold*: once past it, the obstruction (P2) is not motivic at all; it lives in the single realisation's lattice. This is why I stopped pushing: no structure of $\mathcal V$ or $[\nu]$ can see it.
- *The way forward is CDT's*: their arithmetic holonomy bounds are precisely the Diophantine geometry of the realisation layer — a holonomic module with an integral structure and a common complex continuation — and they apply to the whole extension (several functions), not to one coefficient ratio. Our sources, slopes and alignments are the inputs their inequality needs; the Catalan level-144 attempt failed their entry test by a computable margin (1.07 vs 1.15), which is the right kind of failure: quantitative, about the realisation.

**(d) Over other fields and in other worlds.** Beukers' Theorem 4 ($8\zeta(3)-5\sqrt5L(3,\chi_5)\notin\mathbb Q(\sqrt5)$) is the same picture for a motive over $\mathbb Q(\sqrt5)$ — the $\eta$ family with its complex fold, which is why its Apéry limit is not real. The $q$-analogues and function-field analogues (Mellit–Vlasenko) are the same object over different bases. The tower-Frobenius trichotomy should have a clean proof in Beukers–Vlasenko's crystal language, which would make Theorem F (currently modular) available for hypergeometric rows uniformly — e.g. it would *prove* Zudilin's row's 2-adic limit without Padé gymnastics.

## 5. What this suggests we do (theory, not sweeps)

1. **Prove the trichotomy and Theorem F in crystal language** (Beukers–Vlasenko): the Frobenius on the Tate fibre at the cusp, and the constant term as the syntomic regulator. Payoff: Theorem F for every row with a Laurent-polynomial model, including hypergeometric ones; rigidity for Zudilin×E as a corollary. [open, well-posed]
2. **Non-congruence realisations**: define the scan (genus-zero non-congruence groups, integral Hauptmoduls after scaling, CDT's $\lambda$), and look for a second Beukers-type theorem. [open; new territory]
3. **One class, two worlds**: for each Dirichlet period with a conductor prime, construct the cellular-world realisation (Brown/Rhin–Viola-type integrals twisted by roots of unity) of the *same* Eisenstein class; its $p$-adic limit is forced, its decay is what we need. This is the motivic phrasing of Problem 7.3 and the only route to new two-row constructions.
4. **Write the three-piece dictionary into the paper's introduction** as the organising statement (it is currently implicit), and move the "three layers" of the old book there as its correct, modest form.

## 6. A sentence to keep
*An Apéry system is a mixed motive with a maximally unipotent degeneration, written in integers; its Apéry limit is the regulator of the extension class at a place where the degeneration has a free boundary; irrationality is a property not of the motive but of one realisation of it, and is decided by the arithmetic holonomy of that realisation.*
