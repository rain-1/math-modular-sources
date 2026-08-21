# Theory notes 02 — true object, representations, prediction, and how far two lattices go

*Claude (Fable), 2026-08-21. Responds to River's questions after draft 1 of the theory overview. Tags: [lit] [proved] [verified] [conj] [open].*

## 1. The score is a search objective, not a theory — and that is fine

River's point stands: a single scalar pipeline proves exactly $\zeta(2)$ and $\zeta(3)$, because the level-5 and level-6 geometry is simultaneously *cellular* (Brown's $\mathcal M_{0,5}$, $\mathcal M_{0,6}$ integrals) and *modular*, with Fricke symmetry forcing $\lambda_1\lambda_2=1$. Everything else has a negative scalar score (the conjectured scalar barrier).

But the two-lattice construction changes what the "score" *is*. For $m$ realisations of one period with data

$$(\lambda_i,\; k_i,\; \{r_{i,p}\}_{p}) \quad=\quad (\text{archimedean rate},\ \text{denominator depth},\ p\text{-adic convergence rates}),$$

the best irrationality exponent obtainable by the successive-minima argument is a computable function $Q(\lambda, k, r)$. That is the number agents can push — like a kissing number — and it is **constant-specific**, because which realisations exist and at which primes they align depends on the constant. So it does not threaten to yield thirty proofs at once; it yields, at most, one per constant whose realisation-space is rich enough. [conj: that $Q$ is the right invariant; the Catalan work computes it for $m=2$.]

## 2. Are modular forms the true object?

No — they are the sector where the theory is *complete and computable*. The true object is, I think:

> an **integral variation of mixed Hodge structure** on $\mathbb P^1\setminus S$ with a maximally-unipotent point at $0$: a pure piece $\mathcal V$ (rank $w{+}1$, a G-operator / Bombieri–Dwork type) and an extension $0\to\mathcal V\to\mathcal E\to\mathbb Q(0)\to0$ whose period at the MUM point is the constant $\xi$. The Apéry pair is the Taylor expansion of a period of $\mathcal E$ in the integral coordinate.

[lit] On the $\zeta$ side this is Brown's *mixed Tate motives over $\mathbb Z$ realised on $\mathcal M_{0,n}$* (cellular integrals; Apéry $=$ $n=5,6$; Brown–Zudilin $\zeta(5)$ lives here natively). On the Dirichlet-$L$ side it is mixed *elliptic/modular* motives: $\mathcal V=\operatorname{Sym}^w H$ of an elliptic family over a modular curve and the extension is a Beilinson–Kings Eisenstein class (elliptic polylogarithm), regulator $=L(w{+}1,\chi)$. Our Eisenstein Source Theorem is the explicit de Rham shadow of exactly this.

What modular forms buy: integrality for free ($q$-expansion principle), Hecke operators to generate families, CM points to evaluate at, explicit Eisenstein extensions. What they cost: the scalar barrier, and no access to the non-symmetric-power motives ($\zeta(5)$ BZ-type, CY operators of rank 4 that are not $\operatorname{Sym}^3$).

Consequence for strategy: the classification we should aim at is of **integral MUM G-operators with an Eisenstein-type extension**, with the modular ones as the solved sub-census. Almkvist–van Enckevort–van Straten–Zudilin's CY tables are the rank-4 census waiting to be run through the source dictionary. [open]

## 3. Are isogenies the right way to find *different* representations?

Partly. Isogenies/Hecke move within one motive's orbit; the extension class changes only by the scalar $P_C(r)$ and the Frobenius structure is isogeny-invariant almost everywhere. So isogenous realisations have *similar* arithmetic — good for building families and predicting new scalar systems (the hidden $\zeta(7)$ parent, the $p=5,7,11$ completions), bad for producing *complementary* lattices.

The Catalan pair that worked is a **cross-world pair**: Zudilin's row is hypergeometric/cellular (integrality from $\Gamma$-factors), the level-8 row is modular (integrality from $q$-expansions). Their denominator defects come from different mechanisms, so they are close to independent — which is what alignment needs. The archive's $p=3$ pair for $\tfrac12L(2,\chi_{-3})$ (Zagier B and C, $\det M_n = 3^{5(n-1)}/n^4$) is a same-world pair and gave a weaker gain.

**Working hypothesis [conj]: two-lattice gain lives at the primes dividing the level/conductor of the extension**, and the valuation growth rate there is controlled by the local structure at that prime: for $p\mid N$ the modular row is an *overconvergent $p$-adic modular form* on the ordinary locus (this is literally Calegari's 2005 mechanism for the irrationality of $\zeta_p(3)$, $p=2,3$, and it is why the Catalan modular row converges 2-adically at a rate like $8m-1-4s_2(m)$), while the hypergeometric row converges $p$-adically by Lucas/Kummer-type laws at its own rate. Different rates $\Rightarrow$ large $v_p$ of the cross determinant.

Evidence consistent with this: Catalan ($\chi_{-4}$, $N=4,8$) aligns at $p=2$; $L(2,\chi_{-3})$ aligns at $p=3$; $\zeta(5)$ (conductor 1) — the attempted $2{+}3$ "fusion" died (defect theorem) and all $p=3$ Catalan slopes were artefacts. Prediction: $L(3,\chi_5)$ should align at $p=5$; $\beta(4)$ at $p=2$; $\zeta(2k+1)$ at no prime — for odd zeta values the two-lattice method is *structurally* unavailable and one is thrown back on analytic perfection (Apéry) or holonomy vectors (CDT). That is a sharp, testable claim and it would explain a large fraction of the archive's dead ends in one sentence.

Other generators of representations, beyond isogeny: algebraic pullbacks (Shvets Belyi, $\Phi_{\rm tr}=R'(x)\Phi_x$), Hadamard products and Galois traces (the $\beta(4)$, $L(4,\chi_{-3})$ pairs), contiguity/acceleration (Zudilin's second Catalan paper — same motive, different marked de Rham section), Hermite–Padé vectors, and CDT's finite-monodromy dihedral covers. The cross-world pairings are the ones to engineer.

## 4. Predicting analytic and arithmetic properties before computing

Given only $(\Gamma, t, f, \Phi)$ — or more generally the operator and extension — the following are computable in advance:

| property | predicted by | status |
|---|---|---|
| $A_n\in\mathbb Z$ | integrality of $t$ and $f$ ($q$-expansion principle) | [lit] |
| $\lambda_1 = 1/\lvert t_1\rvert$ | nearest cusp/elliptic-point image under $t$ | [lit] |
| $\lambda_2 = 1/\lvert t_2\rvert$ and the exact asymptotic of the linear form | next singularity of the *conditional* function $B-\xi A$; local exponents there (connection constants, CM periods) | [proved] case-by-case — e.g. the recovered $\zeta(5)$ formula: $W$ has a $\log^5$ singularity at $x=-1$, $\lvert x\rvert = 1$, so the linear form decays only like $(\log n)^4/n$. The system is period-pure and analytically dead, and one sees it from the singularity at $-1$ without computing a coefficient. |
| denominator depth $k$ | $w+1$ minus free integrations; free $\iff$ source (or its primitive) exact in the integral lattice | [proved] transfer lemma; exactness criterion [conj] |
| Lucas/Dwork congruences, unit-root structure | Frobenius on $\mathcal V$ (Dwork), twisted by $\chi$ | [proved] for the fifteen (H03, H17) |
| $p$-adic convergence rates of $B/A$ | $p\mid N$: overconvergence/canonical subgroup; $p\nmid N$: Lucas towers | [verified] Catalan; [open] in general |
| scalar vs vector | orientation representation of Fricke/translation on the Eisenstein doublet | [proved] examples |

This table *is* "investigate at high level before computing": it is the thing I would build into a tool. [open]

## 5. How far can the two-lattice method go?

Three levels of ambition, in order:

1. **Theory of $Q$ for $m=2$.** Write the construction-quality exponent as an explicit function of $(\lambda_1,\lambda_2,k_1,k_2,r_{1,p},r_{2,p})$ and the sampling ratio (the 5:8 optimality is a special case). Then the Catalan numbers $0.6588$, $0.9025$, $1-\varepsilon$ become evaluations of one formula, and "nonvanishing" is identified as the single remaining hypothesis. [open, near]
2. **Engineering pairs.** For each target with a ramified prime $p\mid N$: list cross-world realisations (one hypergeometric/cellular, one modular), compute their $p$-adic rates, evaluate $Q$. Candidates: $L(2,\chi_{-3})$ at $p=3$ (Zagier C $+$ a hypergeometric row), $\beta(4)$ at $p=2$, $L(3,\chi_5)$ at $p=5$. [open]
3. **$m\ge3$ and the holonomy hybrid.** Several rows give an $m$-dimensional congruence lattice; CDT's quantified bounds (2510.04156) give effective exponents from *analytic* continuation. A combined statement — "$m$ rows with common continuation *and* aligned $p$-adic lattices" — is the natural generalisation of both; nobody has written it down. [conj]

## 6. What this changes in the plan

- Add a strand: **alignment-prime conjecture** (section 3). Cheap first test: compute the $p$-adic convergence rates of $B_n/A_n$ for the twelve modular sporadics at $p\mid N$ versus $p\nmid N$, and for Zudilin's Catalan row at $p=2$ versus $p=3$.
- The theory paper's programme part should state the true object as an integral MUM VMHS with an Eisenstein-type extension, with modular systems as the complete sub-census — not "three layers".
- Keep the scalar score as a diagnostic, promote $Q(\lambda,k,r)$ to the agent objective.
