\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/10_sqrtdomb.gp
\\ The mixed class of CUSP_MOVE_PROGRAM.md: the sqrt(Domb) cusp-move rows
\\   D5 = (20,10,64):  (n+1)^2 u = (20n^2+26n+10) u - 64 n^2 u ,  xi = (15L - 6 zeta(2))/16
\\   D7 = (20, 4,64):  (n+1)^2 u = (20n^2+14n+ 4) u - 64 n^2 u ,  xi = (6 zeta(2) + 15L)/32
\\ Both are reported to have p-adic limits at BOTH p=2 and p=3.
\\ Question: does a factorial-world decayer align with them at both primes?
\\ Theorem F prediction: the zeta(2) placement is the OUTER one (psi = 1,
\\ phi = chi_(-3) odd), so its p-adic avatar is 0 at every prime.  Hence
\\   xi_p(D5) = (15/16) * (L-avatar) ,   xi_p(D7) = (15/32) * (L-avatar),
\\ and the conductor-6 row, whose xi_p IS the L-avatar, must satisfy
\\   16 xi_p(D5) = 15 xi_p(c6)   and   32 xi_p(D7) = 15 xi_p(c6)   at p = 2 AND p = 3.
default(parisizemax, 8000000000);
default(realprecision, 200);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");

\\ general 2nd-order row (n+1)^2 u = (A n^2 + B n + C) u - D n^2 u
gen2(aA, bB, cC, dD, NN) = {
  my(av = vector(NN+1), bv = vector(NN+1));
  av[1] = 1; av[2] = cC; bv[1] = 0; bv[2] = 1;
  for(n = 1, NN-1,
    av[n+2] = ((aA*n^2 + bB*n + cC)*av[n+1] - dD*n^2*av[n]) / (n+1)^2;
    bv[n+2] = ((aA*n^2 + bB*n + cC)*bv[n+1] - dD*n^2*bv[n]) / (n+1)^2);
  [av, bv];
};

NN = 500;
D5 = gen2(20, 26, 10, 64, NN);
D7 = gen2(20, 14,  4, 64, NN);
MA = row2(7,2,-8,NN);
Lc  = (zetahurwitz(2,1/3) - zetahurwitz(2,2/3))/9;
Z2  = zeta(2);

part1() = {
  print("### archimedean check");
  printf("  D5: b/a at n=%d = %.30f   target (15L-6z2)/16 = %.30f\n",
    NN, 1.0*D5[2][NN+1]/D5[1][NN+1], (15*Lc - 6*Z2)/16);
  printf("  D7: b/a at n=%d = %.30f   target (6z2+15L)/32 = %.30f\n",
    NN, 1.0*D7[2][NN+1]/D7[1][NN+1], (6*Z2 + 15*Lc)/32);
  print("");
  print("### integrality / kappa");
  for(j = 1, 5, my(n = 100*j);
    printf("  n=%3d  den(a5)=%s  den(a7)=%s  v2(a5)=%3d v3(a5)=%3d  v2(a7)=%3d v3(a7)=%3d\n",
      n, denominator(D5[1][n+1]), denominator(D7[1][n+1]),
      valuation(D5[1][n+1],2), valuation(D5[1][n+1],3),
      valuation(D7[1][n+1],2), valuation(D7[1][n+1],3)));
  print("");
  print("### measured slopes (increment valuations of b_n/a_n)");
  for(j = 1, 5, my(n = 100*j, x5 = D5[2][n+1]/D5[1][n+1], y5 = D5[2][n]/D5[1][n],
                    x7 = D7[2][n+1]/D7[1][n+1], y7 = D7[2][n]/D7[1][n]);
    printf("  n=%3d  D5: v2=%5d v3=%5d  |  D7: v2=%5d v3=%5d  |  D5 vs D7: v2(2 xi7 + xi5)=%5d\n",
      n, valuation(x5-y5,2), valuation(x5-y5,3), valuation(x7-y7,2), valuation(x7-y7,3),
      valuation(2*x7 - x5, 2)));
};

part2() = {
  print("");
  print("### two-prime alignment against the conductor-6 decayer");
  my(x5 = D5[2][NN+1]/D5[1][NN+1], x7 = D7[2][NN+1]/D7[1][NN+1]);
  for(j = 1, 3,
    my(ab = [[2,1,60],[1,2,60],[1,1,60]][j], r, xc);
    r = chi6row(ab[1]*ab[3], ab[2]*ab[3]); xc = r[2]/r[1];
    printf("  c6 alpha=%d/%d b=%3d :  v_2(16 xi5 - 15 xi_c6) = %5d   v_3(...) = %5d\n",
      ab[1], ab[2], ab[2]*ab[3], valuation(16*x5 - 15*xc, 2), valuation(16*x5 - 15*xc, 3));
    printf("                        v_2(32 xi7 - 15 xi_c6) = %5d   v_3(...) = %5d\n",
      valuation(32*x7 - 15*xc, 2), valuation(32*x7 - 15*xc, 3));
    printf("     controls: v_2(xi5 - xi_c6)=%4d  v_2(16 xi5 - 16 xi_c6)=%4d  v_3(16 xi5 - 14 xi_c6)=%4d  v_2(32 xi7 - 16 xi_c6)=%4d\n",
      valuation(x5 - xc,2), valuation(16*x5 - 16*xc,2), valuation(16*x5 - 14*xc,3), valuation(32*x7 - 16*xc,2)));
  print("");
  print("### and against the conductor-3 decayer (expect p=3 only)");
  my(r = chi3row(5*24, 3*24), xc = r[2]/r[1]);
  printf("  c3 alpha=5/3 b=72 : v_3(16 xi5 - 15 xi_c3) = %5d   v_2(16 xi5 - 15 xi_c3) = %5d\n",
    valuation(16*x5 - 15*xc, 3), valuation(16*x5 - 15*xc, 2));
  print("");
  print("### is the zeta(2) avatar really 0 ?  row A (xi = zeta(2)/4) at p=2 and p=3");
  printf("  v_2(xi_A) = %5d   v_3(xi_A) = %5d   (n=%d)\n",
    valuation(MA[2][NN+1]/MA[1][NN+1],2), valuation(MA[2][NN+1]/MA[1][NN+1],3), NN);
};

part1();
part2();
quit
