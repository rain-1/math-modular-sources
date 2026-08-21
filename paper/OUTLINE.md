# Modular Apéry systems: sources, slopes, and the two-row design rule

*Working outline — Fable, 2026-08-21. Evidence tags per section: [lit] [proved] [verified] [conj]. Paper = `paper/main.tex`, one file per section in `paper/sections/`.*

## Thesis (abstract in one paragraph)
An Apéry system for a period $\xi$ is controlled by three separate data: the *source* (which L-value), the *singular geometry* (archimedean rates $\lambda_1,\lambda_2$), and the *arithmetic depth* (denominator exponent $k$, and $p$-adic slopes $\sigma_p$). For the fifteen sporadic families and their relatives we (i) prove the companion is the Eichler integral of an explicit Eisenstein/cusp source, so the Apéry limit is a critical $L$-value; (ii) show that the $p$-adic behaviour of $B_n/A_n$ is governed by $\sigma_p=v_p(\lambda_1\lambda_2)+2\kappa_p$, with the product formula tying archimedean failure to $p$-adic gain; (iii) show that $p$-adic Apéry limits are rigid invariants of the extension class at slope primes; (iv) derive a closed-form quality for two-row lattice constructions and a design rule that explains why exactly Catalan and ζ(3) come within 10% of the irrationality threshold and why every other known pair fails; (v) prove that past the threshold no valuation-theoretic device can establish nonvanishing, reducing Catalan's irrationality along this route to a single stated lemma.

## 1. Introduction [me]
- Apéry → Beukers modular picture; the inequality $\lambda_2e^k<1$; the three knobs.
- Summary of results; what is new vs. Beukers/Zagier/Beauville/Almkvist–Zudilin/Cooper/Franke/CDT/Calegari.
- Honesty paragraph: no new irrationality result; the value is structural.

## 2. Modular Apéry systems and their sources [Opus draft from certificates/eisenstein/*.tex + book §3,§5,§8]
- 2.1 Definitions: $(\Gamma,t,f)$, integrality, Picard–Fuchs, companion $B=f\,D^{-(w+1)}\Phi$, $\Phi=f\,Dt$. [lit/proved]
- 2.2 **Theorem A (Source Theorem).** For the twelve modular sporadic pairs $\Phi$ is holomorphic Eisenstein of weight $w+2$ with zero cuspidal projection; explicit table of $\Phi$ and $L(\Phi,s)$. [proved; Lean project exists]
- 2.3 **Theorem B (Apéry limit = critical value).** Fold lemma; limit $=$ rational multiple of $L(\Phi,w+1)$; Eisenstein ⇒ Dirichlet $L$-values, cuspidal ⇒ $L(f,w+1)$ (Domb apparatus, $f_6$; and the new weight-3 row with limit $L(\eta(2\tau)^3\eta(6\tau)^3,2)$). [proved for the fifteen; verified for the new row]
- 2.4 Period annihilation as Hecke linear algebra: barycentric theorem, Gaussian-binomial projectors, Fricke count; relation to Franke. Examples: level-16 ζ(5), level-12 ζ(7) parent (209/1728), β(4) level 24, $L(4,\chi_{-3})$. [proved]
- 2.5 Isogeny/CM: prime completions; CM values → Lambert identities (β(4) example). [proved]

## 3. Archimedean geometry: the score [me + Sonnet table]
- 3.1 $\lambda_1=1/|t_1|$, $\lambda_2=1/|t_2|$; conditional function; Fricke symmetry ⇒ $\lambda_1\lambda_2=\pm1$ (Apéry-perfect).
- 3.2 Denominators: transfer lemma ($k=w+1$), free integrations (Cooper/Bogner). [proved]
- 3.3 **Census table**: every known row with period, $c$, $\lambda_{1,2}$, $k$, score. [verified]
- 3.4 Scalar barrier (conjecture, with E14 census evidence). [conj]

## 4. $p$-adic slopes and extension-class rigidity [me; scripts lattice/]
- 4.1 **Proposition C (slope law).** Casoratian $=-c^{n-1}/n^{w+1}$; $\sigma_p=v_p(c)+2\kappa_p$ when $v_p(a_n)=-\kappa_pn+O(\log n)$. Product formula corollary. [proved]
- 4.2 **Theorem/Conjecture D (rigidity).** Rows with a common period and positive slope at $p$ have $p$-adic limits in the same rational ratio; cross-determinant slope $=\min$. Evidence: B,C,F,$s_{18}$ at 3; Domb,T at 2; negatives; B without archimedean limit. Statement for the Zagier rows as a theorem to prove (Calegari-style). [verified; conj]
- 4.3 Limits of rigidity: unramified tower limits show no relations (padic-apery-limits.tex). [verified]
- 4.4 Slope-free parents: ζ(7) level 24 and all Fricke-purified systems. [verified]

## 5. Two-row constructions: the master formula and design rule [Opus, from ZETA3_TWO_LATTICE §5,§10,§14 + Catalan papers]
- 5.1 The correlated congruence lattice (recall of the Catalan papers, D02 two-minima theorem). [lit-project]
- 5.2 **Theorem E (master formula).** $F,H,\delta$ in closed form; reproduces 0.857914, 0.9025, 0.9010, 0.6589, and the Lean $1-\varepsilon$ at $k^*$. [proved modulo the papers' selection theorem]
- 5.3 **Design rule.** $\delta>1 \iff \log\Lambda_{\rm dec}>k\max(r,1)+r\log\rho_2^{\rm eng}$ (corrected for $\kappa$). Multi-prime version and why it is unrealisable in the census. [proved]
- 5.4 Worked cases: Catalan (Zudilin×E, Zudilin×Nesterenko), ζ(3) Domb×T (0.9010), $L(2,\chi_{-3})$ (four aligned rows, no decayer). Table of budgets. [verified]
- 5.5 Why cross-world pairs: modular rows are integral ($\kappa=0$), hypergeometric rows carry $\kappa_p>0$. [proved from 4.1]

## 6. Past the threshold: what cannot work [me, from CATALAN_AUDIT/DIRECTIONAL/POSITIVITY/EXPLICIT + Sol]
- 6.1 The rationality-kernel vector; rational-$G^*$ control experiment; box = covolume ⇒ Dirichlet. [proved]
- 6.2 Absorption theorems: extra homogeneous congruences, factorisation of $h_n$ (ceiling identity $\gamma_h=k^*\log2+H^*$), extra rows, coefficient-directional divisibility — all absorbed. [proved/verified]
- 6.3 Positivity: both source forms positive; moment-family identity $K_N^{(n)}=K_Z^{(3n)}w^n$; explicit positive $P$ cannot cancel; SOS floor. [proved/verified]
- 6.4 **Lemma P2** as the exact remaining content; proof that it is of irrationality strength. [proved]

## 7. Problems
- Rigidity theorem for Zagier rows (4.2) via 2-/3-adic Eichler integrals.
- Scalar barrier theorem (3.4).
- Hypergeometric partners with $\kappa_p>0$ for $L(f,2)$, $L(2,\chi_{-3})$, ζ(3) (the only way to raise $\delta$).
- New sporadic sequences with large $\Lambda$: genus-zero census beyond eta quotients (growth filter).
- **Adelic Padé**: replace lattice selection by simultaneous real/2-adic Padé for the Catalan moment family; what such a construction must satisfy ($v_2(q_n)<(24-k)n$ on the kernel side).
- Cooper free integrations: exactness criterion.

## Appendices
- A. Verification scripts index (`verify/`, `lattice/`), with what each PASS means.
- B. Ledger cross-reference (which ledger IDs each theorem closes).

## Division of labour
- Fable: §1, §3.1–3.2, §4, §6, §7, final edit.
- Opus agent 1: §2 LaTeX draft (from the proved sources), §5 LaTeX draft.
- Sonnet agent: census table (§3.3, §5.4) as `paper/tables/census.tex` + `budgets.tex`, generated by a script from exact data; appendix A.
