/* ---------------------------------------------------------------------------
   2-adic identification for the level-16 zeta(5) row.
   Predicted (Euler-factor criterion / CONJ_D_PROOF Thm R''):
       P(X) = (1-X)(1-4X)(1-16X)(1-64X)   [X <-> V_2],  w+1 = 5,
       Q(X) = P(X)/(1-X) = (1-4X)(1-16X)(1-64X),  Q(2^-5) = Q(1/32) = -7/16,
       xi_2 = -Q(1/32) * (1/2) zeta_2(5) = (7/32) zeta_2(5).
   zeta_2 = Kubota-Leopoldt L_2(s, 1), computed by Washington Thm 5.11 with F=4.
   Log: lattice/zeta5_two_row/padic_id.log
   ------------------------------------------------------------------------- */
default(parisizemax, 8000000000);
LOG = "lattice/zeta5_two_row/padic_id.log";
W(s) = write(LOG, s);
if(type(KTR)!="t_INT", KTR = 260);       /* truncation of the Washington series */

/* omega(a) for p=2 is the character mod 4:  omega(a) = 1 if a=1 mod 4, -1 if a=3 mod 4 */
om(a) = if(a%4==1, 1, -1);
tw(a) = a/om(a);                          /* <a> = a/omega(a),  = 1 mod 4 */

/* Washington, Introduction to Cyclotomic Fields, Thm 5.11, chi = 1, p = 2, F = 4:
   L_2(s,1) = 1/(4(s-1)) * sum_{a in {1,3}} <a>^{1-s} sum_{j>=0} binom(1-s,j) (4/a)^j B_j  */
Lp2(s, K) = (1/(4*(s-1))) * sum(i=1,2, my(a=[1,3][i]); \
      tw(a)^(1-s) * sum(j=0, K, binomial(1-s, j) * (4/a)^j * bernfrac(j)));

/* --- validation: exact interpolation L_2(1-n,1) = -(1-2^(n-1)) B_n/n for n even --- */
W("=== validation of the 2-adic zeta code (exact rational identities) ===");
ok = 1;
for(m=1, 6, n = 2*m; \
   lhs = Lp2(1-n, n);            /* binom(n,j) vanishes for j>n so the series is finite & exact */ \
   rhs = -(1 - 2^(n-1))*bernfrac(n)/n; \
   if(lhs != rhs, ok = 0); \
   W(Str("  n=", n, "  L_2(1-n,1)=", lhs, "   -(1-2^(n-1))B_n/n=", rhs, "   equal: ", lhs==rhs)));
W(Str("all interpolation checks exact: ", ok));

/* --- zeta_2(5) as a 2-adic number --- */
PR = 2*KTR;
z2of5 = (1/(4*(5-1))) * sum(i=1,2, my(a=[1,3][i]); \
      (tw(a)+O(2^PR))^(1-5) * sum(j=0, KTR, binomial(1-5, j) * ((4/a)+O(2^PR))^j * bernfrac(j)));
W("");
W(Str("zeta_2(5) = L_2(5,1),  K=", KTR, ":"));
W(Str("  v_2 = ", valuation(z2of5,2)));
W(Str("  = ", z2of5 + O(2^40)));

/* the two predicted constants */
xistar   = (7/32) * z2of5;                     /* -Q(1/32) * (1/2) zeta_2(5) */
xinaive  = -(-217/512) * (1/2) * z2of5;        /* -P(1/32) * (1/2) zeta_2(5)  (no depletion) */
W(Str("predicted xi_2 = (7/32) zeta_2(5)  : v_2 = ", valuation(xistar,2), "   = ", xistar+O(2^40)));
W(Str("alt (undepleted) (217/1024) zeta_2(5): v_2 = ", valuation(xinaive,2), " = ", xinaive+O(2^40)));

/* --- the measured 2-adic limit --- */
read("lattice/zeta5_two_row/level16_rows.txt");
NN = #An;
W("");
W(Str("rows loaded, N = ", NN));
for(i=1,6, n = [NN-1, NN-2, NN-11, NN-51, NN-101, NN-201][i]; \
  if(n>2 && An[n+1]!=0, r = Bn[n+1]/An[n+1]; \
    W(Str("  n=", n, "  v_2(B_n/A_n - (7/32)zeta_2(5)) = ", valuation(r - xistar, 2), \
          "     v_2(B_n/A_n) = ", valuation(r,2), \
          "     v_2(B_n/A_n - (217/1024)zeta_2(5)) = ", valuation(r - xinaive,2)))));

/* --- rational-multiple scan: which rational r has B_n/A_n = r*zeta_2(5) ? --- */
W("");
W("--- scan: v_2( B_N/A_N - c*zeta_2(5) ) over candidate rationals c ---");
n = NN-1; r = Bn[n+1]/An[n+1];
cands = [7/32, 217/1024, 7/16, 7/64, 1/2, 1, 7/8, -7/32, 31/32*7/32, 217/512];
for(i=1,#cands, W(Str("  c=", cands[i], "   v_2(diff) = ", valuation(r - cands[i]*z2of5, 2))));

/* --- lindep: express the measured xi_2 against zeta_2(5) --- */
W("");
W(Str("ratio (B_N/A_N)/zeta_2(5) = ", (r/z2of5) + O(2^60)));
W(Str("bestappr of that ratio, viewed 2-adically: v_2(ratio - 7/32) = ", valuation(r/z2of5 - 7/32, 2)));
W("");
W("DONE");
quit;
