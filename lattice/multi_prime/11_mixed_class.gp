\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/11_mixed_class.gp
\\ The genuinely two-slope-prime members of the sqrt(Domb) cusp-move orbit
\\ (CUSP_MOVE_PROGRAM.md section 4.1 / section 6):
\\   D.2 = (-28,-6,192) : P = -28n^2-22n-6 , Q = 192n^2-96n
\\   D.3 = ( 8, 0,-48)  : P = 8n^2+2n      , Q = -48n^2+24n
\\   D.4 = ( 8, 6,-48)  : P = 8n^2+14n+6   , Q = -48n^2-24n
\\   D.6 = (-28,-12,192): P = -28n^2-34n-12, Q = 192n^2+96n
\\ Their Q(n) carries a linear term, so the Casoratian is not c^n and
\\ sigma_p is not v_p(c): both p=2 and p=3 survive.
\\ Question: does any factorial-world decayer align with them at BOTH primes?
default(parisizemax, 8000000000);
default(realprecision, 300);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");

\\ (n+1)^2 u_(n+1) = P(n) u_n - Q(n) u_(n-1) with P,Q given as coefficient vectors
gen(pv, qv, NN) = {
  my(av = vector(NN+1), bv = vector(NN+1), pp, qq);
  av[1] = 1; av[2] = pv[3]; bv[1] = 0; bv[2] = 1;
  for(n = 1, NN-1,
    pp = pv[1]*n^2 + pv[2]*n + pv[3];
    qq = qv[1]*n^2 + qv[2]*n + qv[3];
    av[n+2] = (pp*av[n+1] - qq*av[n]) / (n+1)^2;
    bv[n+2] = (pp*bv[n+1] - qq*bv[n]) / (n+1)^2);
  [av, bv];
};

NN = 500;
R2 = gen([-28,-22,-6], [192,-96,0], NN);
R3 = gen([  8,  2, 0], [-48, 24,0], NN);
R4 = gen([  8, 14, 6], [-48,-24,0], NN);
R6 = gen([-28,-34,-12],[192, 96,0], NN);
Lc = (zetahurwitz(2,1/3) - zetahurwitz(2,2/3))/9;
Z2 = zeta(2);

show(nm, rr) = {
  my(x = rr[2][NN+1]/rr[1][NN+1], y = rr[2][NN]/rr[1][NN], ld);
  printf("  %-4s  b/a = %.40f\n", nm, 1.0*x);
  ld = lindep([1.0*x, 1, Z2, Lc], 60);
  printf("        lindep[xi,1,zeta2,L]  = %s\n", ld);
  printf("        den(a_n) at n=500 : %s     v2(a)=%d v3(a)=%d\n",
    denominator(rr[1][NN+1]), valuation(rr[1][NN+1],2), valuation(rr[1][NN+1],3));
  for(j = 2, 5, my(n = 100*j, u = rr[2][n+1]/rr[1][n+1], v = rr[2][n]/rr[1][n]);
    printf("        n=%3d  v2(incr)=%5d (sig2~%5.2f)   v3(incr)=%5d (sig3~%5.2f)\n",
      n, valuation(u-v,2), 1.0*valuation(u-v,2)/n, valuation(u-v,3), 1.0*valuation(u-v,3)/n));
};

print("### the four two-prime rows of the sqrt(Domb) orbit");
show("D.2", R2); show("D.3", R3); show("D.4", R4); show("D.6", R6);

print("");
print("### internal p-adic relations (CUSP_MOVE_PROGRAM section 6)");
{my(x2 = R2[2][NN+1]/R2[1][NN+1], x3 = R3[2][NN+1]/R3[1][NN+1],
    x4 = R4[2][NN+1]/R4[1][NN+1], x6 = R6[2][NN+1]/R6[1][NN+1]);
  printf("  v_2(xi2 + 2 xi3) = %5d   v_3(xi2 + 2 xi3) = %5d\n", valuation(x2+2*x3,2), valuation(x2+2*x3,3));
  printf("  v_2(xi4 -   xi6) = %5d   v_3(xi4 -   xi6) = %5d\n", valuation(x4-x6,2), valuation(x4-x6,3));
  print("");
  print("### alignment against the conductor-6 decayer (xi_p = the L(2,chi_-3) avatar at BOTH primes)");
  my(rc = chi6row(60,120), xc = rc[2]/rc[1], rc3 = chi3row(120,72), xc3 = rc3[2]/rc3[1]);
  print("  scan of rationals c with xi = c * xi_c6 :");
  for(i = 1, 14,
    my(cs = [1,2,3,4,1/2,1/3,1/4,3/2,5/2,15/16,15/32,5/8,8/5,2/3][i]);
    printf("    c=%-6s  D.2: v2=%5d v3=%5d | D.3: v2=%5d v3=%5d\n", cs,
      valuation(x2 - cs*xc,2), valuation(x2 - cs*xc,3),
      valuation(x3 - cs*xc,2), valuation(x3 - cs*xc,3)));
  print("  lindep of xi_D2 against the two p-adic constants, p=3, prec 3^60 :");
  printf("    v_3(xi2)=%d  v_3(xi3)=%d  v_3(xi_c6)=%d  v_3(xi_c3)=%d\n",
    valuation(x2,3), valuation(x3,3), valuation(xc,3), valuation(xc3,3));
  print("  best rational relation xi_D2 = c * xi_c6 at p=3 (bestappr on the 3-adic ratio):");
  my(rt = x2/xc);
  printf("    ratio mod 3^40 = %s\n", lift(Mod(numerator(rt),3^40)*Mod(denominator(rt),3^40)^-1));
  printf("    bestappr of the 3-adic ratio: %s\n", bestappr(centerlift(Mod(numerator(rt),3^40)*Mod(denominator(rt),3^40)^-1)/1, 10^8));
  my(rt2 = x2/xc);
  printf("    ratio mod 2^40 = %s   bestappr: %s\n",
    lift(Mod(numerator(rt2),2^40)*Mod(denominator(rt2),2^40)^-1),
    bestappr(centerlift(Mod(numerator(rt2),2^40)*Mod(denominator(rt2),2^40)^-1)/1, 10^8));}
quit
