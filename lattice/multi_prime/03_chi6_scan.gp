\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/03_chi6_scan.gp
\\ The conductor-6 well-poised row for L(2,chi_(-3)): denominators, slopes,
\\ and alignment with the modular rows C (p=3) and F (p=2 and p=3).
\\ Predictions (Theorem F):  xi_3 = zeta_3(2) = 2 xi_3(C) = (8/5) xi_3(F)
\\                           xi_2 = (4/5) L_2(2,chi_12) = (8/5) xi_2(F)
\\ (the decayer has r_inf = 1 for L(2,chi_(-3)); E_2(2) = 5/4, E_3(2) = 1).
default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");

NM = 700;
rF = row2(17, 6, 72, NM);
rC = row2(10, 3,  9, NM);
xi2F = rF[2][NM+1] / rF[1][NM+1];
xi3C = rC[2][NM+1] / rC[1][NM+1];

scan6(pp, qq, KM, step) = {
  print("");
  print("### chi_6 row  alpha = ", pp, "/", qq, "   a = ", pp, "k,  b = ", qq, "k");
  print("   k    b  v2(denQ) v3(denQ) otherden | v3(xi-2xi3C) v2(xi-(8/5)xi2F) | v3incr v2incr | log|Q|/b  log|QL-P|/b");
  my(prev = 0, Lc = (zetahurwitz(2,1/3) - zetahurwitz(2,2/3))/9);
  for(k = 1, KM,
    if(k % step != 0 && k != KM, next);
    my(r = chi6row(pp*k, qq*k), qq2 = r[1], pp2 = r[2], xk = pp2/qq2,
       d3, d2, i3, i2, bb = qq*k, den = denominator(qq2), oth);
    oth = den / 2^valuation(den,2) / 3^valuation(den,3);
    d3 = valuation(xk - 2*xi3C, 3);
    d2 = valuation(xk - (8/5)*xi2F, 2);
    i3 = if(prev == 0, 0, valuation(xk - prev, 3));
    i2 = if(prev == 0, 0, valuation(xk - prev, 2));
    prev = xk;
    printf("  %3d %3d   %5d   %5d   %8d |  %6d   %6d |  %5d %5d | %8.4f %10.4f\n",
      k, bb, valuation(den,2), valuation(den,3), oth, d3, d2, i3, i2,
      log(abs(1.0*qq2))/bb, log(abs(1.0*(qq2*Lc - pp2)))/bb));
};

scan6(2, 1, 40, 8);
scan6(3, 2, 30, 6);
scan6(4, 3, 24, 6);
scan6(1, 1, 30, 6);
scan6(1, 2, 24, 6);
scan6(5, 3, 24, 6);
quit
