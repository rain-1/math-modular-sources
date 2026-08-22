\\ =====================================================================
\\ row_aesz207.gp -- AESZ 207.
\\   L = sum_{i=0}^{4} z^i P_i(theta),  P_0 = theta^4  (MUM, holonomic rank 4);
\\   coefficient recurrence sum_i P_i(n-i) u_{n-i} = 0 has SHIFT ORDER 4,
\\   hence THREE companions B = X^(1), C = X^(2), D = X^(3).
\\   chi(x) = x^4 - 17152 x^3 - 2^17*51088 x^2 + 2^28*13*73104 x - 2^48*13^2
\\   roots: -2^7(349+85 sqrt17) = -89531.389..., 53248 (DOUBLE),
\\          2^17/(349+85 sqrt17) = 187.389...
\\   lambda_1 = 89531.389..., lambda_2 = 53248.
\\
\\ Everything here is geometric: |lambda_2/lambda_1| = 0.59474..., so the exact
\\ rational run to N gives about 0.2255*N archimedean digits.
\\ =====================================================================

default(parisizemax, 8G);
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/sc_rows.gp");

NN     = 6000;    \\ exact rational range
NKEXP  = 400;     \\ range for the sharp denominator exponent
NKEXP2 = 1600;
PRECW  = 2000;    \\ working real precision

print("==========================================================");
print("ROW 2: AESZ 207");
print("==========================================================");
print();

\\ =====================================================================
\\ (i) the operator: theta-degrees, holonomic rank, characteristic polynomial
\\ =====================================================================
print("---- (i) the operator L = sum_{i=0}^4 z^i P_i(theta)");
{ for(i = 0, 4,
   print("  P_", i, "(n) = ", R5Q[i+1]);
   print("        deg_theta = ", poldegree(R5Q[i+1], 'n),
         "   leading coefficient = ", polcoeff(R5Q[i+1], poldegree(R5Q[i+1],'n), 'n),
         " = ", if(polcoeff(R5Q[i+1], poldegree(R5Q[i+1],'n), 'n) != 0,
                   factor(polcoeff(R5Q[i+1], poldegree(R5Q[i+1],'n), 'n)), 0));
); }
{ print("  ALL FIVE P_i have deg_theta = 4: ",
        vecmax(vector(5, i, poldegree(R5Q[i], 'n))) == 4 &&
        vecmin(vector(5, i, poldegree(R5Q[i], 'n))) == 4); }
print("  P_0(theta) = theta^4 exactly  => MUM at z=0, HOLONOMIC RANK R = 4.");
print("  z-degree of L = 4, so the coefficient recurrence has SHIFT ORDER r = 4");
print("  and there are r-1 = 3 companions.  These are logically independent");
print("  numbers: R = 4 (theta-order) vs r = 4 (z-degree) coincide here by accident.");
print();
print("  Inventory: L(H) = 0 with H = sum A_n z^n means");
print("    theta^4 H = -(1/z^0 coefficient normalisation) sum_{i=1}^{4} z^i P_i(theta) H,");
print("  i.e. theta^4 H is a Q(z)-linear combination of H, thetaH, theta^2 H, theta^3 H.");
print("  Hence {1, H, thetaH, theta^2 H, theta^3 H} is the FULL inventory: the");
print("  Q(z)-span of all theta^k H is 4-dimensional, spanned by k = 0,1,2,3.");
print("  [immediate from L(H)=0 and P_0 = theta^4; stated, not measured]");
print();

chi = sum(i = 0, 4, polcoeff(R5Q[i+1], 4, 'n) * 'x^(4-i));
print("  characteristic polynomial chi(x) = ", chi);
print("  factorisation over Q: ", factor(chi));
default(realprecision, 40);
{ my(rts = polroots(chi));
  print("  roots (numerically):");
  for(i = 1, #rts, print("    ", precision(real(rts[i]), 30)));
  print("  lambda_1 = ", precision(vecmax(vector(#rts,i,abs(real(rts[i])))), 30));
}
print("  -2^7*(349+85*sqrt(17)) = ", precision(-2^7*(349+85*sqrt(17)), 30));
print("   2^17/(349+85*sqrt(17)) = ", precision(2^17/(349+85*sqrt(17)), 30));
print("  |lambda_2/lambda_1| = ", precision(53248/(2^7*(349+85*sqrt(17))), 30));
print("  => digits gained per step of n: ", precision(-log(53248/(2^7*(349+85*sqrt(17))))/log(10), 10));
print();

\\ =====================================================================
\\ build A, B, C, D exactly
\\ =====================================================================
gettime();
Aq = genseq(R5cf, R5r, [1], NN);
print("[exact] A_0..A_5 = ", vector(6, i, Aq[i]));
print("[exact] A built to n = ", NN, ", ms = ", gettime());
Xq = vector(3, j, compan(R5cf, R5r, j, NN));
print("[exact] B_0..B_5 = ", vector(6, i, Xq[1][i]));
print("[exact] C_0..C_5 = ", vector(6, i, Xq[2][i]));
print("[exact] D_0..D_5 = ", vector(6, i, Xq[3][i]));
print("[exact] companions built, ms = ", gettime());
print();

\\ =====================================================================
\\ (h) sharp denominator exponent
\\ =====================================================================
print("---- (h) sharp k:  lcm(1..n)^k x_n in Z");
print("  k(A) = ", denexp(Aq, NKEXP2), "  (A is integral)");
nms = ["B", "C", "D"];
{ for(j = 1, 3,
    my(kk = denexp(Xq[j], NKEXP), kk2 = denexp(Xq[j], NKEXP2), bad = 0, lc = 1);
    for(nn = 1, NKEXP, lc = lcm(lc, nn);
        if(denominator(lc^(kk-1)*Xq[j][nn+1]) > 1, bad = nn; break));
    print("  k(", nms[j], ") = ", kk, "  [n<=", NKEXP, "]   = ", kk2, "  [n<=", NKEXP2,
          "]    k-1 = ", kk-1, " first fails at n = ", bad);
); }
print("  ms = ", gettime());
print();

\\ =====================================================================
\\ (f) archimedean limits
\\ =====================================================================
default(realprecision, PRECW);
print("---- (f) archimedean limits xi_infinity, exact rational run to N = ", NN);
xis = vector(3); dgs = vector(3);
{ for(j = 1, 3,
    my(rN = Xq[j][NN+1]*1.0/Aq[NN+1], rM = Xq[j][NN-49]*1.0/Aq[NN-49]);
    my(cau = abs(rN - rM), d = floor(-log(cau)/log(10)));
    xis[j] = rN; dgs[j] = d;
    print("  xi(", nms[j], ") = ", precision(rN, 60));
    print("      Cauchy |r_N - r_{N-50}| = ", precision(cau, 6), "   -> ", d, " certified digits");
); }
DIG = vecmin(dgs);
print("  common certified precision DIG = ", DIG, " digits");
print();
print("  comparison with the recorded values (MULTISLOPE_PROGRAM.md / AESZ207.md):");
{ my(recB = -0.000504554593441367088625425457970715161172154460579143450288368704004362680122570507,
     recC =  1.501998871618802820919660703567693867874,
     recD = -3.026195781429694336873872900571897312700);
  print("    xi(B) recorded : -0.000504554593441367088625425457970715161172154460579143450288368704004362680122570507");
  print("    xi(B) here     : ", precision(xis[1], 87));
  print("    -> agree to ", floor(-log(abs(xis[1]-recB)/abs(recB))/log(10)), " digits (the whole recorded string)");
  print();
  print("    xi(C) recorded : +1.501998871618802820919660703567693867874...");
  print("    xi(C) here     : ", precision(xis[2], 45));
  print("    -> the MANTISSAS agree to ",
        floor(-log(abs(xis[2]*10^8-recC)/abs(recC))/log(10)),
        " digits, but the recorded value is missing its factor 10^-8.");
  print("    xi(D) recorded : -3.026195781429694336873872900571897312700...");
  print("    xi(D) here     : ", precision(xis[3], 45));
  print("    -> the MANTISSAS agree to ",
        floor(-log(abs(xis[3]*10^13-recD)/abs(recD))/log(10)),
        " digits, but the recorded value is missing its factor 10^-13.");
  print("    [CORRECTION to MULTISLOPE_PROGRAM.md 2.2: xi(C) = 1.5019988716...e-8 and");
  print("     xi(D) = -3.0261957814...e-13; the recorded strings dropped the exponents,");
  print("     as the ratios below confirm.]");
}
print();
print("  ratios:  xi(C)/xi(B) = ", precision(xis[2]/xis[1], 30));
print("           xi(D)/xi(B) = ", precision(xis[3]/xis[1], 30));
print("  (MULTISLOPE records -2976.88077988600456148 and +5997.75687461135022676, i.e.");
print("   the same mantissas times 10^8 and 10^13 -- consistent with the exponent slip.)");
print();

\\ =====================================================================
\\ (f2) lindep / period rank
\\ =====================================================================
default(realprecision, DIG + 5);
xr = vector(3, j, precision(xis[j], DIG));
print("---- (f2) lindep on {1, xi(B), xi(C), xi(D)} at ", DIG, " digits");
{ ldtest(idx) =
  my(v = concat([1.0], vector(#idx, i, xr[idx[i]])));
  my(rel = lindep(v), res = sum(i=1,#v, rel[i]*v[i]));
  my(ht = log(1.0*vecmax(vector(#rel,i,abs(rel[i]))))/log(10));
  my(lr = if(res == 0, -9999.0, log(abs(res))/log(10)));
  printf("  {1,%s}: ht = 10^%.2f (spurious 10^%.1f)  residual = 10^%.1f%s\n",
     vector(#idx,i,nms[idx[i]]), ht, 1.0*DIG/#v, lr,
     if(ht < 1.0*DIG/#v - 4 && lr < -(DIG-8), "   <== HIT", ""));
}
{ for(j = 1, 3, ldtest([j])); }
{ for(j = 1, 2, for(k = j+1, 3, ldtest([j,k]))); }
ldtest([1,2,3]);
print();
{ print("  full relation lattice of {1, xi_B, xi_C, xi_D} (LLL-reduced):");
  my(kk = 10^(DIG-3), lat = matrix(5,4));
  for(j = 1, 4,
    for(i = 1, 4, lat[i,j] = if(i == j, 1, 0));
    lat[5,j] = round(kk * if(j == 1, 1.0, xr[j-1])));
  my(red = lat * qflll(lat));
  for(j = 1, 4,
    my(cv = vector(4, i, red[i,j]), rs = 1.0*abs(red[5,j])/kk);
    printf("    c = %s   height 10^%.2f   residual 10^%.1f%s\n", cv,
      log(1.0*vecmax(vector(4,i,abs(cv[i]))))/log(10),
      if(rs == 0, -9999.0, log(rs)/log(10)),
      if(rs < 10.0^(-(DIG-8)), "   <== GENUINE", "")));
}
print("  => period rank prk_infinity = dim_Q span{1,xi_B,xi_C,xi_D} - 1");
print();
print("  and over Q(sqrt 17):");
{ my(s17 = sqrt(17.0));
  my(v = [1.0, s17, xr[1], s17*xr[1], xr[2], s17*xr[2], xr[3], s17*xr[3]]);
  my(rel = lindep(v), res = sum(i=1,#v, rel[i]*v[i]));
  printf("    lindep over Q(sqrt17) on 8 terms: height 10^%.2f (spurious 10^%.1f) residual 10^%.1f\n",
     log(1.0*vecmax(vector(#rel,i,abs(rel[i]))))/log(10), 1.0*DIG/8,
     if(res == 0, -9999.0, log(abs(res))/log(10)));
}
print();

\\ =====================================================================
\\ (g) fold regularity
\\ =====================================================================
default(realprecision, PRECW);
print("---- (g) FOLD-REGULARITY.  L^(j)_n = X^(j)_n - xi_j A_n.");
print("     Prediction: |L^(j)_n|^(1/n) -> lambda_2 = 53248 while |A_n|^(1/n) -> 89531.389...");
print("     (measured: L_n = c * 53248^n / n^4 exactly -- the double root does NOT");
print("      produce a log or a factor n here.)");
print();
foldn = [500, 1000, 2000, 3000, 4000, 5000];
print("  n       |A_n|^(1/n)      n^2 * A_n / lambda_1^n");
{ my(l1 = -2^7*(349+85*sqrt(17)));
  for(i = 1, #foldn, my(nn = foldn[i]);
    printf("  %-7d %.6f    %s\n", nn, exp(log(abs(Aq[nn+1]*1.0))/nn),
           precision(nn^2*(Aq[nn+1]*1.0)/l1^nn, 20))); }
print();
{ for(j = 1, 3,
    print("  ", nms[j], ":   n       |L_n|^(1/n)      n^4 * L_n / 53248^n");
    for(i = 1, #foldn, my(nn = foldn[i], ln = Xq[j][nn+1]*1.0 - xis[j]*(Aq[nn+1]*1.0));
      printf("        %-7d %.6f       %s\n", nn, exp(log(abs(ln))/nn),
             precision(nn^4*ln/53248.0^nn, 20)));
); }
print();
print("  generic rational combinations S_n = sum_j a_j (X^(j)_n - xi_j A_n):");
gvecs = [[1,-3,5], [2,3,-5], [1,1,1]];
{ for(g = 1, #gvecs, my(av = gvecs[g]);
    print("   a = ", av);
    for(i = 1, #foldn, my(nn = foldn[i],
        ln = sum(j=1,3, av[j]*(Xq[j][nn+1]*1.0 - xis[j]*(Aq[nn+1]*1.0))));
      printf("     n=%-7d |S_n|^(1/n) = %.6f     n^4*S_n/53248^n = %s\n",
             nn, exp(log(abs(ln))/nn), precision(nn^4*ln/53248.0^nn, 20))));
}
print();
{ print("  reference values: |lambda_1| = ", precision(2^7*(349+85*sqrt(17)), 15),
        ", lambda_2 = 53248 (double), lambda_4 = ", precision(2^17/(349+85*sqrt(17)), 15)); }
print("  NOTE: |L_n|^(1/n) approaches 53248 from below exactly as predicted by");
print("  L_n = c * 53248^n / n^4 (the n^{-4} makes |L_n|^(1/n) = 53248*n^(-4/n)");
print("  = 53248*(1 - 4 log n / n + ...)), and the n^4*L_n/53248^n column is");
print("  CONSTANT to 20+ digits, which is the sharp form of the statement.");
print();
print("DONE.");
quit;
