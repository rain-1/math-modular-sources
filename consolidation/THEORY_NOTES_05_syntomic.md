# Theory notes 05 — the syntomic interpretation of Theorems B\* and F

*Fable, 2026-08-22. A conceptual note: what the Apéry limits are, motivically, and why the $p$-adic ones exist only at slope primes. Tags: [lit] [proved] [verified] [conj] [refuted].*

## 1. The object

Fix a modular Apéry row with Eisenstein source $\Phi=P(V)\,E^{\psi,\varphi}_{w+2}$ on $Y=Y_0(N)$ (or the relevant cover), weight-$w$ form $F=f$, Hauptmodul $t$. The motivic object behind the row is Beilinson's **Eisenstein symbol**
$$\mathrm{Eis}^w\in H^{w+1}_{\mathcal M}\big(\mathcal E^w,\mathbb Q(w+1)\big),$$
the universal elliptic curve's Eisenstein class of weight $w$ (Beilinson; Kings' construction via the elliptic polylogarithm; Huber–Kings 1999 for the degeneration at the cusps). Pushed to the base it is an extension of variations of (mixed) motives over $Y$:
$$0\longrightarrow \operatorname{Sym}^wH(w+1)\longrightarrow\mathcal E\longrightarrow\mathbb Q(0)\longrightarrow0 ,$$
and the oldform polynomial $P(V)$ is a Hecke/degeneracy combination of such classes. The row's companion $B=F\cdot D^{-(w+1)}\Phi$ is the de Rham period of this extension expanded in the integral coordinate $t$: the Eichler integral $\Theta=D^{-(w+1)}\Phi$ is the extension's de Rham class written against the section $F$. [lit: Beukers' interpretation of Apéry; Kings 2015 for the polylog/Eisenstein identification]

## 2. The archimedean regulator (Theorem B\*)

Huber–Kings (Invent. Math. 1999) show that the Eisenstein class **degenerates at the cusps to the cyclotomic (Beilinson–Deninger) elements** in $H^1_{\mathcal M}(\mathbb Q(\mu_N),\mathbb Q(w+1))$. Beilinson's and Deninger's computation of the Deligne regulator of those elements is the classical Dirichlet $L$-value $L(\psi,w+1)$ (up to the explicit rational factor). Theorem B\* is the Apéry shadow of this: the fold of a second-order row is the cusp $0$ and $\xi_\infty=\lim_{y\to0}\Theta(iy)=L(\Phi,w+1)$ by the Mellin integral; for third-order rows the fold is the Fricke fixed point and the period polynomial forces the same value. So

> **$\xi_\infty$ = Deligne regulator of the degenerate Eisenstein class at the cusp**, evaluated against the row's de Rham section. [proved at the level of Theorem B\*; the identification with the literature regulator is a matching of normalisations we have not written out]

## 3. The $p$-adic regulator (Theorem F)

The syntomic (rigid, Bloch–Kato) regulator of the same elements is the **Kubota–Leopoldt value** — Huber–Kings, *Duke* 2003 (Bloch–Kato for Dirichlet characters), via Coleman power series — and, on the modular curve itself, Bannai–Kings (*Amer. J. Math.* 2010) compute the syntomic realisation of the elliptic polylogarithm/Eisenstein class on the **ordinary locus**: it is given by the $p$-adic Eisenstein–Kronecker series, i.e. by Serre/Katz's **$p$-adic Eisenstein series** $G^{(p)}$, whose constant term is the $p$-adic $L$-value. Theorem F is the Apéry shadow of this:

> The Eichler integral $\Theta=D^{-(w+1)}\Phi$ equals $Q(V/\!\cdot)\,G^{(p)}-Q(w+1)\kappa_p$ when the Euler factor divides $P$; the overconvergent function $F\cdot(\Theta-\xi)$ exists iff $\xi=-Q(w+1)\kappa_p$; so **$\xi_p$ = syntomic regulator (constant term of the $p$-adic Eisenstein series) = Kubota–Leopoldt value.** [proved in the census cases; the identification "constant term of $G^{(p)}$ = syntomic regulator of the cyclotomic element" is Bannai–Kings/Huber–Kings]

Corollary (the rigidity principle, Conjecture D): two rows realising the same extension class have the *same* regulators at every place, hence $p$-adic limits in the ratio of their archimedean ones corrected by the Euler factor. This is functoriality of the regulator under Hecke/degeneracy operators, nothing more.

## 4. Why good primes are different (and a refuted shortcut)

At a good prime $p\nmid N$ the syntomic regulator still exists and still equals the $p$-adic $L$-value (Bannai–Kings work on the ordinary locus for any $p$). Yet $b_n/a_n$ has no $p$-adic limit, and neither — **tested and refuted, `lattice/depleted_companion.gp`** — does the $p$-depleted companion $b^{[p]}_n=[t^n]\big(F\cdot D^{-(w+1)}\Phi^{[p]}\big)$, which one might hope converges to $\kappa_p$ at every prime since $F\cdot(\Theta^{[p]}+\kappa_p)=F\,G^{(p)}$ is overconvergent (row C at $p=5,7,11$: increment valuations $\approx0$ for both the full and the depleted companion, $n\le400$).

The reason is geometric and it is the content of the trichotomy (Notes 04, §1):

- **Overconvergent** means analytic on a strict neighbourhood of the ordinary locus, i.e. extending *into* the supersingular annuli.
- At a **slope prime** $p\mid N$ (bad reduction), the component of the ordinary locus containing the cusp $\infty$ is a disc $\{|t|_p\le1\}$ bounded by a *single* supersingular annulus (Calegari's picture for $X_0(2)$ at $p=2$: ordinary locus $\|f\|\le1$ and $\|f\|\ge2^{12}$, one annulus between). Overconvergence of $F(\Theta-\xi)$ is therefore exactly "radius of convergence $>1$ in $t$", which is the $p$-adic Apéry limit.
- At a **good prime**, the supersingular residue discs of $X_0(N)$ lie *on the unit circle* $|t|_p=1$ of the cusp coordinate. An overconvergent function's expansion at the cusp converges on the residue disc $|t|<1$ and on the thin annuli beyond the ordinary discs, but not across the supersingular discs sitting at $|t|=1$: its radius of convergence is exactly $1$. Hence no coefficient-ratio limit — for any companion, depleted or not. The regulator is encoded instead in the **values** of $F\,G^{(p)}$ on the ordinary residue discs (e.g. at Teichmüller/CM points), not in the Taylor coefficients at the cusp; and the tower limits $\Lambda_a$ of the coefficient ratios measure the Frobenius splitting of the extension (Notes 04), which is a different invariant.

So the trichotomy "(good: splits, $\Gamma_p$-assemblies) / (slope, inner: unipotent, KL value) / (slope, outer or cuspidal: 0)" is a statement about **where the supersingular locus sits relative to the cusp**, not about the regulator — which is the same $p$-adic $L$-value throughout. [conj as a general statement; verified 144/144 cells; the geometric explanation is standard rigid geometry of modular curves]

## 5. What this buys

1. A one-sentence description of the whole theory: *an Apéry row is a Hecke combination of Beilinson's Eisenstein symbol written in an integral coordinate; its Apéry limits are the regulators of that class at the places of $\mathbb Q$ where the cusp's ordinary disc has a free boundary.*
2. A prediction about **which primes** can carry a $p$-adic Apéry limit for *any* realisation of a given extension class: exactly the primes of bad reduction of the host curve at which the Euler factor is cancelled by the oldform polynomial — hence the "alignment primes" of the two-row method are motivic invariants of the pair, and no cleverness with coordinates changes them.
3. A route to Problem 7.1(iv) (proving the tower-Frobenius trichotomy): it is the statement that the Frobenius on the fibre of $\mathcal E$ at the cusp is the Frobenius of the Tate curve twisted by $\psi(p)p^{-w}$, i.e. Bannai–Kings' description of the syntomic class at the cusp. [open]
4. What it does *not* buy: any progress on Lemma P2 (Catalan). The regulator picture explains alignment; it has nothing to say about the rationality-kernel vector.

## 6. References to pin down (for the paper's §7 or an appendix)
- Beilinson, *Higher regulators of modular curves* (1986); Deninger, *Higher regulators and Hecke L-series of imaginary quadratic fields*; Huber–Kings, *Degeneration of l-adic Eisenstein classes and of the elliptic polylog*, Invent. Math. 135 (1999); Huber–Kings, *Bloch–Kato conjecture and main conjecture of Iwasawa theory for Dirichlet characters*, Duke 119 (2003); Bannai–Kings, *$p$-adic elliptic polylogarithm, $p$-adic Eisenstein series and Katz measure*, Amer. J. Math. 132 (2010); Kings, *Eisenstein classes, elliptic Soulé elements and the $\ell$-adic elliptic polylogarithm* (2015); Coleman, *Classical and overconvergent modular forms* (1996); Calegari (2005).
