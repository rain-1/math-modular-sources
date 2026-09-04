# Z.-W. Sun, "Catalan's constant is irrational" (arXiv:2609.04176, 3 Sep 2026): the proof fails

*Fable, 2026-09-04. Scripts: `sun_test.gp` (the paper's central inequality evaluated with the true value of G), `sun_padic.gp` (its p-adic lemma tested with fake rational values of G).*

## The architecture

Assume G = a/q. Tails T_m = (-1)^m (G - S_{m-1}) satisfy T_m + T_{m+1} = 1/(2m+1)^2 and 0 < T_m < 1/(2m+1)^2; weighted tails u_m = T_m/(2m+1). With B, S = floor(B/20), N = 2B+S+3, D = 2B, Pi_i = prod_{h=1}^B (2(h+i)+1)^2, the paper forms the N x N determinant qhat_B of the matrix with columns i^r/Pi_i (r < 2B), u_{i+j} (j <= S) and three binomial columns, and shows qhat_B = +- F_B det R[A,J] / prod_i Pi_i with F_B = prod_{r<2B} r! and R the finite-difference residual matrix. Let H_B^min be the denominator of q^S qhat_B; then N_B = q^S H_B^min qhat_B is a nonzero integer (Prop. 3.3, correct). Theorem 9.1 claims log H_B^min + log|qhat_B| <= -delta_0 B^2 + o(B^2) with delta_0 > 0.00966, hence |N_B| < 1 for large B, contradiction.

## The refutation [exact, one line]

Pi_i is odd, q^S R has 2-integral entries (the paper's own Lemma 5.2), so v_2(N_B) >= v_2(F_B) = sum_{r<2B} v_2(r!) = 2B^2 - O(B log B). Hence
  log|N_B| >= v_2(F_B) log 2 = (2 log 2) B^2 (1 + o(1)) = 1.386 B^2 (1+o(1)),
for every B, whatever G is. This contradicts Theorem 9.1's upper bound -0.0097 B^2. So the archimedean ledger of Section 9 is wrong: the power of 2 in F_B sits in the numerator of qhat_B and is never cancelled by the (odd) denominators; the paper's Remark 9.2 believes the 2-adic effect is (19/200) log 2 per B^2, it is 2 log 2 per B^2.

## Numerical confirmation [verified, B <= 200]

With the true Catalan constant (3000 digits), X_B := log|qhat_B| + log Hbar_B, where Hbar_B = prod_{p odd} p^{[A_p - sum_nu m_{p^nu}]_+} is exactly the paper's upper bound for H_B^min (Cor. 5.3), computed with its own minimisation (greedy marginal costs, exact):

| B | S | log|det R| | odd p-adic content M | v_2(F_B) log 2 | X_B | claimed -delta_0 B^2 | X_B/B^2 |
|---|---|---|---|---|---|---|---|
| 20 | 1 | 98.6 | -98.8 | 471.3 | 668.8 | -3.9 | 1.672 |
| 60 | 3 | 1281.6 | -401.8 | 4669.0 | 6352.4 | -34.8 | 1.765 |
| 100 | 5 | 4066.3 | -625.5 | 13286.2 | 17978.1 | -96.6 | 1.798 |
| 200 | 10 | 19020.8 | 202.4 | 54159.7 | 72949.8 | -386.5 | 1.824 |

X_B = log|det R| - M + v_2(F_B) log 2 grows like +1.8 B^2. Even discarding the 2-adic term, log|det R| ~ 0.47 B^2 exceeds the odd p-adic content M ~ 0: the size of the residual determinant is larger, not smaller, than its guaranteed divisibility, so the claimed gains (Lambda_mid, 83/2400) do not beat the archimedean size even on the odd primes.

## What is correct in the paper

The rank theorem (Thm 2.1), the scalar identity (Prop. 3.1), the Cauchy--Binet factorisation (Section 4), and the p-adic lower bound (Lemma 5.1) all check: `sun_padic.gp` verifies Lemma 5.1 at every odd prime for fake rationals a/q with q up to 10^30 at B = 40, 60, 80, no violation. The failure is confined to the archimedean estimate (Prop. 9.1 and its "raw coefficient 39/200"), which is off by about 1.8 B^2.

## Context

The acknowledgment states the proof was produced in conversation with an AI and "passed the verification of ChatGPT 5.6 Solar". For comparison, our own accounting of Catalan's constant in the Calegari--Dimitrov--Tang framework (`consolidation/INVENTORY_BOUND.md` section 4, `consolidation/GAMMA15_CLOSURE.md`) puts the best achievable margin at -8 nats.
