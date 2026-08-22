\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/02_decayer_padic.gp
\\ Step 1b/2: does the chi_(-3) hypergeometric decayer have a p-adic limit at BOTH
\\ p = 3 and p = 2, and does it align with the modular rows there?
\\ Predictions from the Euler-factor criterion (EULER_CRITERION.md Theorem F):
\\   r_p = r_inf / E_p(w+1),  E_p(s) = 1 - chi_(-3)(p) p^(-s),  w = 1.
\\   decayer r_inf = 1.  p=3: E_3(2) = 1  -> xi_3 = Lam_3 = zeta_3(2) = 2 xi_3(C).
\\   p=2: chi_(-3)(2) = -1, E_2(2) = 5/4 -> xi_2 = (4/5) Lam_2 = (8/5) xi_2(F).
default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");

NM = 700;
rF = row2(17, 6, 72, NM);
rC = row2(10, 3,  9, NM);
rB = row2( 9, 3, 27, NM);
rE = row2(12, 4, 32, NM);
xi2F = rF[2][NM+1] / rF[1][NM+1];
xi3C = rC[2][NM+1] / rC[1][NM+1];
xi3B = rB[2][NM+1] / rB[1][NM+1];
xi3F = rF[2][NM+1] / rF[1][NM+1];

hdr() = {
  print("row F Cauchy v_2 : ", valuation(rF[2][NM+1]/rF[1][NM+1] - rF[2][NM]/rF[1][NM], 2), "   (3N = ", 3*NM, ")");
  print("row F Cauchy v_3 : ", valuation(rF[2][NM+1]/rF[1][NM+1] - rF[2][NM]/rF[1][NM], 3), "   (2N = ", 2*NM, ")");
  print("row C Cauchy v_3 : ", valuation(rC[2][NM+1]/rC[1][NM+1] - rC[2][NM]/rC[1][NM], 3));
  print("check 4 xi3F - 5 xi3C : v_3 = ", valuation(4*xi3F - 5*xi3C, 3));
  print("check   xi3B -   xi3C : v_3 = ", valuation(xi3B - xi3C, 3));
  print("");
};

scan(pp, qq, KM, step) = {
  print("### chi_-3 decayer  alpha = ", pp, "/", qq, "   (a = ", pp, "k, b = ", qq, "k)");
  print("  k    b   v3(den Q) v2(den Q)   v3(xi_k - 2 xi3C)   v2(xi_k - (8/5) xi2F)   v3incr  v2incr");
  my(prev = 0);
  for(k = 1, KM,
    if(k % step != 0 && k != KM, next);
    my(r = chi3row(pp*k, qq*k), xk = r[2]/r[1], d3, d2, i3, i2);
    d3 = valuation(xk - 2*xi3C, 3);
    d2 = valuation(xk - (8/5)*xi2F, 2);
    i3 = if(prev == 0, 0, valuation(xk - prev, 3));
    i2 = if(prev == 0, 0, valuation(xk - prev, 2));
    prev = xk;
    print("  ", k, "   ", qq*k, "   ", valuation(denominator(r[1]),3), "  ", valuation(denominator(r[1]),2),
          "        ", d3, "               ", d2, "         ", i3, "  ", i2));
  print("");
};

hdr();
scan(2, 1, 40, 8);
scan(3, 2, 30, 6);
scan(5, 3, 24, 6);
scan(4, 3, 18, 6);
scan(5, 4, 16, 4);
scan(7, 4, 16, 4);
scan(9, 5, 14, 4);
scan(1, 1, 30, 6);
scan(1, 2, 30, 6);
quit
