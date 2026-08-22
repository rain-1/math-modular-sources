\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/04_chi6_rates.gp
\\ Full rate profile of the conductor-6 decayer family, plus high-precision
\\ two-prime alignment against the modular rows C, B, F and controls.
default(parisizemax, 8000000000);
default(realprecision, 200);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");

NM  = 900;
rF  = row2(17, 6, 72, NM);
rC  = row2(10, 3,  9, NM);
rB  = row2( 9, 3, 27, NM);
xi2F = rF[2][NM+1]/rF[1][NM+1];
xi3C = rC[2][NM+1]/rC[1][NM+1];
xi3B = rB[2][NM+1]/rB[1][NM+1];
xi3F = xi2F;
Lc  = (zetahurwitz(2,1/3) - zetahurwitz(2,2/3))/9;

banner() = {
  print("modular reference rows, N = ", NM);
  print("  v_2 Cauchy of F : ", valuation(rF[2][NM+1]/rF[1][NM+1] - rF[2][NM]/rF[1][NM], 2), "  (3N = ", 3*NM, ")");
  print("  v_3 Cauchy of C : ", valuation(rC[2][NM+1]/rC[1][NM+1] - rC[2][NM]/rC[1][NM], 3), "  (2N = ", 2*NM, ")");
  print("  v_3(4 xi3F - 5 xi3C) = ", valuation(4*xi3F - 5*xi3C, 3));
  print("  L(2,chi_-3) = ", Lc);
  print("");
};

\\ full profile for one family member
prof(pp, qq, KM) = {
  my(k = KM, bb = qq*KM, r, qv, pv, den, dp, oth, lamd, lam, kap2, kap3, nu, eta,
     s3, s2, i3, i2, rprev, xk, xprev);
  r = chi6row(pp*k, qq*k); qv = r[1]; pv = r[2];
  rprev = chi6row(pp*(k-1), qq*(k-1)); xprev = rprev[2]/rprev[1];
  xk = pv/qv;
  den = denominator(qv); dp = denominator(pv);
  kap2 = valuation(den,2)/bb; kap3 = valuation(den,3)/bb;
  oth  = den/2^valuation(den,2)/3^valuation(den,3);
  nu   = log(1.0*oth)/bb;
  eta  = log(1.0*lcm(den,dp))/bb;
  lamd = log(abs(1.0*qv))/bb;
  lam  = log(abs(1.0*(qv*Lc - pv)))/bb;
  s3 = valuation(xk - 2*xi3C, 3)/bb;
  s2 = valuation(xk - (8/5)*xi2F, 2)/bb;
  i3 = valuation(xk - xprev, 3)/bb;
  i2 = valuation(xk - xprev, 2)/bb;
  printf("  a/b=%d/%d k=%3d b=%3d | logLam=%8.4f loglam=%8.4f | kap2=%6.3f kap3=%6.3f nu=%7.3f eta=%7.3f | sig2=%6.3f sig3=%6.3f | incr2=%6.3f incr3=%6.3f | w2=%7.3f w3=%7.3f\n",
    pp, qq, k, bb, lamd, lam, kap2, kap3, nu, eta, s2, s3, i2, i3, s2-2*kap2, s3-2*kap3);
  [pp/qq, lamd, lam, kap2, kap3, nu, eta, s2, s3];
};

controls() = {
  my(k = 20, qq = 3, pp = 4, bb = qq*k, r = chi6row(pp*k, qq*k), xk = r[2]/r[1]);
  print("");
  print("### controls at a=", pp*k, " b=", bb, "  (only the predicted rational should give a large valuation)");
  print("  p=2, target c*xi2F :");
  for(i = 1, 10,
    my(cs = [8/5, 1, 2, 4/5, 5/8, 3/2, 5/4, 16/5, 8/3, 5/2][i]);
    print("     c = ", cs, "   v_2 = ", valuation(xk - cs*xi2F, 2)));
  print("  p=3, target c*xi3C :");
  for(i = 1, 8,
    my(cs = [2, 1, 4, 5/2, 8/5, 5/4, 1/2, 3][i]);
    print("     c = ", cs, "   v_3 = ", valuation(xk - cs*xi3C, 3)));
  print("  p=3 against B and F :  v_3(xi - 2 xi3B) = ", valuation(xk - 2*xi3B, 3),
        "   v_3(5 xi - 8 xi3F) = ", valuation(5*xk - 8*xi3F, 3));
  print("  p=2 against F      :  v_2(5 xi - 8 xi2F) = ", valuation(5*xk - 8*xi2F, 2));
};

banner();
print("### rate profile (all rates per unit of b)");
prof(2,1,60);
prof(3,2,40);
prof(5,3,32);
prof(4,3,32);
prof(5,4,28);
prof(1,1,40);
prof(2,3,32);
prof(1,2,32);
prof(1,3,24);
prof(1,4,20);
controls();
quit
