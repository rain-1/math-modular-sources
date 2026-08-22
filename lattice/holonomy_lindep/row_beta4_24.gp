\\ =====================================================================
\\ row_beta4_24.gp -- beta(4) at level 24 = Sym^3(Zagier E).
\\   A_n = [t^n] F(t)^3,  F = Zagier E row (12,4,32).
\\   Minimal recurrence R4cf: shift order 6, poly degree 4, MUM rank 4.
\\   chi(x) = (x-4)^3 (x-8)^3  =>  lambda_1 = 8 (TRIPLE), lambda_2 = 4 (TRIPLE).
\\   Five companions X^(1..5); note also X^(0) = A.
\\
\\ METHOD NOTE.  lambda_1 = 8 is a triple root, so X^(j)_n/A_n converges only
\\ like 1/log n and the Cauchy criterion |r_N - r_{N-50}| certifies ~2 digits
\\ even at N = 1200.  We therefore use the exact asymptotic shape, forced by
\\ A = F^3 with F = P(z) log(1-8z) + Q(z) near z = 1/8:
\\
\\    u_n = (8^n/n) * sum_{k>=0} n^{-k} ( c_{k,2} L^2 + c_{k,1} L + c_{k,0} ),
\\    L = log n,      xi_j = c^{(j)}_{0,2} / c^{(A)}_{0,2},
\\
\\ and fit it by exact linear algebra at 450-digit precision.  The model is
\\ validated independently: c^{(A)}_{0,2} must equal 24/pi^3 (since P(1/8) = -2/pi),
\\ and it does, to ~70 digits.
\\ =====================================================================

default(parisizemax, 12G);
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/sc_rows.gp");

NEXACT = 1200;
NKEXP  = 400;
NFLOAT = 1000000;
PRECF  = 450;

Ainit = [1,12,108,880,6876,52752];

\\ ---- memory-light floating sampler: returns u_n for n in the sorted list want
fsamp(Cf, r, ini, Nn, want) = {
  my(buf = vector(r+1, i, 0.0), out = vector(#want), wi = 1, m0 = #ini);
  for(i = 1, m0, for(t = 1, r, buf[t] = buf[t+1]); buf[r+1] = ini[i]*1.0);
  for(nn = m0, Nn,
    my(s = 0.0);
    for(i = 1, r, s += subst(Cf[i+1], 'n, nn) * buf[r+2-i]);
    my(v = -s / subst(Cf[1], 'n, nn));
    for(t = 1, r, buf[t] = buf[t+1]); buf[r+1] = v;
    while(wi <= #want && want[wi] == nn, out[wi] = v; wi++);
  );
  out;
};

\\ ---- fit  n*u_n/8^n ~ sum_{k=0..mm} n^{-k}(c2 L^2 + c1 L + c0)
\\ vals = values at the nodes ns (in that order); returns [c0,c1,c2,...]
fitasy(vals, ns, mm) = {
  my(np = 3*(mm+1), mt = matrix(np,np), rhs = vector(np));
  for(i = 1, np,
    my(nn = ns[i], lg = log(nn*1.0), c = 0);
    for(k = 0, mm, for(j = 0, 2, c++; mt[i,c] = lg^j / nn^k));
    rhs[i] = nn * vals[i] / 8.0^nn;
  );
  matsolve(mt, rhs~);
};

gnodes(mm, rat) = vector(3*(mm+1), i, round(NFLOAT/rat^(i-1)));
agreedig(x, y) = if(x == y, precision(x), floor(-log(abs((x-y)/x))/log(10)));

\\ fit configurations: [m, ratio];  the first is the reference
FITS = [[18, 1.10], [17, 1.10], [16, 1.12], [15, 1.14]];
foldn = [1000, 5000, 20000, 100000, 300000, 1000000];

print("==========================================================");
print("ROW 1: beta(4) at level 24 = Sym^3(Zagier E)");
print("==========================================================");
print();

\\ =====================================================================
\\ 0. exact sequences and structure
\\ =====================================================================
default(realprecision, 60);
gettime();
Aex = r4A(NEXACT);
print("[exact] A_0..A_9 = ", vector(10,i,Aex[i]));
{ print("[exact] R4cf reproduces the series cube up to n=", NEXACT, " : ",
        Aex == genseq(R4cf, R4r, Ainit, NEXACT)); }
{ print("[exact] A also satisfies R4cf from n=1 on (i.e. A = X^(0)) : ",
        Aex == genseq(R4cf, R4r, [1], NEXACT)); }
chi = sum(i=0, R4r, polcoeff(R4cf[i+1], poldegree(R4cf[i+1]), 'n) * 'x^(R4r-i));
print("[exact] chi(x) = ", chi);
print("[exact] chi factors as ", factor(chi));
print("[exact] build ms = ", gettime());
print();
Xex = vector(5, j, compan(R4cf, R4r, j, NEXACT));
print("[exact] X^(1)_0..X^(1)_8 = ", vector(9,i,Xex[1][i]));
print("[exact] X^(5)_0..X^(5)_8 = ", vector(9,i,Xex[5][i]));
print("[exact] companions built, ms = ", gettime());
print();

\\ =====================================================================
\\ (d) sharp denominator exponent
\\ =====================================================================
print("---- (d) sharp k:  lcm(1..n)^k * x_n in Z  for all n <= ", NKEXP);
print("  k(A) = ", denexp(Aex, NKEXP), "   (A is integral)");
{ for(j = 1, 5,
    my(kk = denexp(Xex[j], NKEXP), kk2 = denexp(Xex[j], NEXACT), bad = 0, lc = 1);
    for(nn = 1, NKEXP, lc = lcm(lc, nn);
        if(denominator(lc^(kk-1)*Xex[j][nn+1]) > 1, bad = nn; break));
    print("  k(X^(", j, ")) = ", kk, "  [n<=", NKEXP, "]   = ", kk2, "  [n<=", NEXACT,
          "]   k-1 = ", kk-1, " first fails at n = ", bad);
); }
print("  ms = ", gettime());
print();

\\ =====================================================================
\\ (a) naive Cauchy, then the structured fit
\\ =====================================================================
print("---- (a) NAIVE Cauchy criterion on r_n = X^(j)_n/A_n  [exact ratios, N=", NEXACT, "]");
{ for(j = 1, 5,
    my(rN = Xex[j][NEXACT+1]*1.0/Aex[NEXACT+1],
       rM = Xex[j][NEXACT-49]*1.0/Aex[NEXACT-49]);
    printf("  j=%d  r_N = %.18f   |r_N - r_{N-50}| = %.3e   -> %d digits certified\n",
           j, rN, abs(rN-rM), floor(-log(abs(rN-rM))/log(10)));
); }
print("  -> convergence is O(1/log n): the Cauchy criterion is useless here.");
print();

default(realprecision, PRECF);
print("---- (a) structured asymptotic fit,  N = ", NFLOAT, ", precision ", PRECF, " digits");

\\ union of all nodes we need
allnodes = Set(concat(concat(vector(#FITS, q, gnodes(FITS[q][1], FITS[q][2]))), concat(foldn, [NEXACT])));
allnodes = vecsort(vector(#allnodes, i, eval(allnodes[i])));
nodeidx = Map(); for(i = 1, #allnodes, mapput(nodeidx, allnodes[i], i));
getv(vals, nn) = vals[mapget(nodeidx, nn)];

gettime();
Avals = fsamp(R4cf, R4r, [1], NFLOAT, allnodes);
print("  [float] A sampled at ", #allnodes, " nodes, ms = ", gettime());
{ print("  [float] relative error of A_", NEXACT, " vs exact: ",
        precision(abs(getv(Avals,NEXACT)/(Aex[NEXACT+1]*1.0) - 1), 5)); }

{ fitof(vals, q) = my(mm = FITS[q][1], ns = gnodes(mm, FITS[q][2]));
                  fitasy(vector(#ns, i, getv(vals, ns[i])), ns, mm); }

cAall = vector(#FITS, q, fitof(Avals, q));
cA = cAall[1];
print("  [float] c^{(A)}_{0,2} = ", precision(cA[3], 72));
print("          24/pi^3       = ", precision(24/Pi^3, 72));
print("          agreement = ", agreedig(cA[3], 24/Pi^3), " digits  <== independent validation of the model");
print("  [float] c^{(A)}_{0,1} = ", precision(cA[2], 50));
print("  [float] c^{(A)}_{0,0} = ", precision(cA[1], 50));
print();

xis = vector(5); dig = vector(5); ctri = vector(6); Xfold = vector(5);
ctri[1] = [cA[1], cA[2], cA[3]];
{ for(j = 1, 5,
    my(Xvals = fsamp(R4cf, R4r, concat(vector(j,i,0),[1]), NFLOAT, allnodes));
    print("  [float] X^(", j, ") sampled, ms = ", gettime(),
          "; rel. error at n=", NEXACT, ": ",
          precision(abs(getv(Xvals,NEXACT)/(Xex[j][NEXACT+1]*1.0) - 1), 5));
    my(cc = vector(#FITS, q, fitof(Xvals, q)));
    my(xv = vector(#FITS, q, cc[q][3]/cAall[q][3]));
    my(d = vecmin(vector(#FITS-1, q, agreedig(xv[1], xv[q+1]))));
    xis[j] = xv[1]; dig[j] = d;
    ctri[j+1] = [cc[1][1], cc[1][2], cc[1][3]];
    Xfold[j] = vector(#foldn, i, getv(Xvals, foldn[i]));
    print("  xi_", j, " = ", precision(xv[1], 70));
    print("        cross-fit agreement: ", d, " digits");
); }
Afold = vector(#foldn, i, getv(Avals, foldn[i]));
print();
DIG = vecmin(dig);
print("---- common certified precision DIG = ", DIG, " digits");
print();

\\ =====================================================================
\\ (b) lindep on {1, xi_1..xi_5}
\\ =====================================================================
default(realprecision, DIG + 5);
xr = vector(5, j, precision(xis[j], DIG));
print("---- (b) lindep on {1, xi_1..xi_5} and every subset, at ", DIG, " digits");
{ printf("  spurious level: a (k+1)-term lindep at D=%d digits returns junk of\n", DIG); }
{ printf("  height ~10^(D/(k+1)); a genuine relation of height 10^h leaves residual\n"); }
{ printf("  ~10^(-(D - (k+1)h)).  We call a HIT only if the residual is at full precision.\n"); }
print();
{ ldtest(idx) =
  my(v = concat([1.0], vector(#idx, i, xr[idx[i]])));
  my(rel = lindep(v), res = sum(i=1,#v, rel[i]*v[i]));
  my(ht = log(1.0*vecmax(vector(#rel,i,abs(rel[i]))))/log(10));
  my(lr = if(res == 0, -9999.0, log(abs(res))/log(10)));
  my(spur = 1.0*DIG/#v);
  printf("  {1,xi%s}: rel=%s  ht=10^%.2f (spurious 10^%.1f)  res=10^%.1f%s\n",
     idx, rel~, ht, spur, lr, if(ht < spur - 2 && lr < -(DIG-8), "   <== HIT", ""));
  [rel, ht, lr];
}
{ for(j=1,5, ldtest([j])); }
print();
{ for(j=1,4, for(k=j+1,5, ldtest([j,k]))); }
print();
{ for(j=1,3, for(k=j+1,4, for(l=k+1,5, ldtest([j,k,l])))); }
print();
{ for(j=1,2, for(k=j+1,3, for(l=k+1,4, for(m4=l+1,5, ldtest([j,k,l,m4]))))); }
print();
ldtest([1,2,3,4,5]);
print();
print("  deflation: after imposing the relations found, re-run on the quotient.");
{ my(rel = lindep(concat([1.0], xr)));
  print("  full lindep vector: ", rel~);
  \\ eliminate xi_2 via 2 xi_1 - 3 xi_2 = 0 and retest the rest
  my(rr = ldtest([1,3,4,5]));
  print("  (that is the second relation, if the HIT tag is present)");
}
print();

\\ =====================================================================
\\ (c) constant catalogue
\\ =====================================================================
print("---- (c) identification against a constant catalogue, ", DIG, " digits");
beta4 = sumalt(nn=0, (-1)^nn/(2*nn+1)^4);
{ catnm = ["1","Pi","Pi^2","Pi^3","Pi^4","log2","G","beta(4)","zeta(3)",
           "zeta(2)","Pi*G","Pi^2*G","Pi^3*G","zeta(3)*Pi","zeta(3)/Pi","G/Pi",
           "G/Pi^2","Pi*log2","log2^2","log2^3","zeta(3)/Pi^3","beta(4)/Pi^3",
           "G^3/Pi^3","1/Pi","G^2/Pi^2","sqrt(2)","sqrt(3)","log2/Pi","zeta(3)/Pi^2"]; }
{ catvl = [1.0, Pi, Pi^2, Pi^3, Pi^4, log(2), Catalan, beta4, zeta(3),
           zeta(2), Pi*Catalan, Pi^2*Catalan, Pi^3*Catalan, zeta(3)*Pi,
           zeta(3)/Pi, Catalan/Pi, Catalan/Pi^2, Pi*log(2), log(2)^2, log(2)^3,
           zeta(3)/Pi^3, beta4/Pi^3, Catalan^3/Pi^3, 1/Pi, Catalan^2/Pi^2,
           sqrt(2), sqrt(3), log(2)/Pi, zeta(3)/Pi^2]; }
print("  beta(4) = ", precision(beta4, 40));
{ for(j = 1, 5, my(hits = 0);
   for(c = 1, #catvl,
     my(rel = lindep([xr[j], catvl[c]]), res = 0);
     res = rel[1]*xr[j] + rel[2]*catvl[c];
     if((res == 0 || log(abs(res))/log(10) < -(DIG-8)) && vecmax(abs(rel)) < 10^(DIG/4),
        hits++; printf("  xi_%d : 2-term HIT  %d*xi + %d*%s = 0\n", j, rel[1], rel[2], catnm[c]));
   );
   for(c = 2, #catvl,
     my(rel = lindep([xr[j], 1.0, catvl[c]]), res = 0);
     res = rel[1]*xr[j] + rel[2] + rel[3]*catvl[c];
     if((res == 0 || log(abs(res))/log(10) < -(DIG-8)) && vecmax(abs(rel)) < 10^(DIG/5),
        hits++; printf("  xi_%d : 3-term HIT  %d*xi + %d + %d*%s = 0\n", j, rel[1], rel[2], rel[3], catnm[c]));
   );
   if(hits == 0, print("  xi_", j, " : no identification in the catalogue (2- or 3-term)"));
); }
print();
print("  raw log-cube coefficients c^{(j)}_{0,2} = xi_j * 24/pi^3 :");
{ for(j = 1, 5, print("    j=", j, "  ", precision(xis[j]*24/Pi^3, 45))); }
print();

\\ =====================================================================
\\ (e) fold-regularity
\\ =====================================================================
default(realprecision, PRECF);
print("---- (e) FOLD-REGULARITY.  L^(j)_n = X^(j)_n - xi_j A_n.");
print("     Hypothesis to test: |L^(j)_n|^(1/n) -> lambda_2 = 4 while |A_n|^(1/n) -> 8.");
print();
print("  n           |A_n|^(1/n)");
{ for(i=1,#foldn, printf("  %-10d %.9f\n", foldn[i], exp(log(abs(Afold[i]))/foldn[i]))); }
print();
{ for(j = 1, 5,
    print("  j=", j, ":  n           |L_n|^(1/n)    |L_n/A_n|       log(n)*|L_n/A_n|");
    for(i = 1, #foldn, my(nn = foldn[i], ln = Xfold[j][i] - xis[j]*Afold[i]);
      printf("          %-10d %.9f    %.6e    %.9f\n",
             nn, exp(log(abs(ln))/nn), abs(ln/Afold[i]), log(nn*1.0)*abs(ln/Afold[i])));
); }
print();
print("  exact limits of the last column: log(n)*L^(j)_n/A_n -> kappa_j with");
print("  kappa_j = (c^{(j)}_{0,1} - xi_j c^{(A)}_{0,1}) / c^{(A)}_{0,2}  [from the fit]:");
{ default(realprecision, 40);
  for(j = 1, 5,
    print("    kappa_", j, " = ",
          precision((ctri[j+1][2] - xis[j]*ctri[1][2])/ctri[1][3], 30)));
  default(realprecision, PRECF);
}
print("  i.e.  X^(j)_n - xi_j A_n  ~  kappa_j * A_n / log n  ~  kappa_j * 24/pi^3 * 8^n * log n / n.");
print("  The linear form does NOT drop to radius 1/4: it stays at radius 1/8, damped");
print("  only by the single logarithmic factor 1/log n.  FOLD-REGULARITY FAILS.");
print();
print("  generic rational combinations  S_n = sum_j a_j (X^(j)_n - xi_j A_n):");
gvecs = [[1,-3,5,-7,11], [2,3,-5,7,-11], [1,1,1,1,1]];
{ for(g = 1, #gvecs, my(av = gvecs[g]);
    print("   a = ", av);
    for(i = 1, #foldn, my(nn = foldn[i], ln = sum(j=1,5, av[j]*(Xfold[j][i] - xis[j]*Afold[i])));
      printf("     n=%-10d |S_n|^(1/n) = %.9f   log(n)*|S_n/A_n| = %.9f\n",
             nn, exp(log(abs(ln))/nn), log(nn*1.0)*abs(ln/Afold[i]))));
}
print();

\\ =====================================================================
\\ (e2) the CORRECT fold question for a triple root: is there a rational
\\      combination of A, X^(1..5) that is log-free at z = 1/8 (hence O(4^n))?
\\      That needs the three coefficients (c_{0,2}, c_{0,1}, c_{0,0}) to vanish
\\      simultaneously -- codimension 3, not 1.
\\ =====================================================================
default(realprecision, DIG + 5);
print("---- (e2) codimension-3 fold test.");
print("     For a triple root the local solution space at z=1/8 is spanned by");
print("     y1^3, y1^2 y2, y1 y2^2, y2^3 with y2 = y1 log(1-8z) + hol.  A combination");
print("     has radius 1/4 (i.e. |u_n|^(1/n) -> 4) iff it is log-FREE at z=1/8, i.e.");
print("     iff c_{0,2} = c_{0,1} = c_{0,0} = 0 -- THREE conditions on the six");
print("     sequences A = X^(0), X^(1),...,X^(5).  Real kernel is 3-dimensional;");
print("     the question is whether it contains a rational point.");
print();
\\ integer-relation search for SIMULTANEOUS conditions: find short integer
\\ c in Z^6 with  sum_j c_j * mt[i,j] ~ 0  for every row i.
{ simrel(mt, dg) =
  my(nr = matsize(mt)[1], nc = matsize(mt)[2], kk = 10^dg);
  my(lat = matrix(nc + nr, nc));
  for(j = 1, nc,
    for(i = 1, nc, lat[i,j] = if(i == j, 1, 0));
    for(i = 1, nr, lat[nc+i, j] = round(kk * mt[i,j]));
  );
  my(tt = qflll(lat), red = lat*tt);
  my(best = 1, bs = -1);
  for(j = 1, nc,
    my(tail = vecmax(vector(nr, i, abs(red[nc+i, j]))));
    my(head = vecmax(vector(nc, i, abs(red[i, j]))));
    if(head > 0 && (bs < 0 || tail < bs), bs = tail; best = j);
  );
  my(cv = vector(nc, i, red[i, best]));
  [cv, vector(nr, i, abs(red[nc+i, best]))/kk];
}
{ my(mt = matrix(3, 6, i, j, precision(ctri[j][i], DIG)));
  print("  matrix of (c_{0,0}; c_{0,1}; c_{0,2}) for X^(0)=A, X^(1..5):");
  for(i = 1, 3, print("    row ", i-1, ": ", vector(6, j, precision(mt[i,j], 25))));
  print();
  my(rr = simrel(mt, DIG - 3));
  my(cv = rr[1], resid = rr[2]);
  my(ht = log(1.0*vecmax(vector(6, i, abs(cv[i]))))/log(10));
  print("  3-condition LLL search (kill c_{0,2}, c_{0,1}, c_{0,0} simultaneously):");
  print("    c = ", cv);
  printf("    height = 10^%.2f   spurious level 10^%.1f  (3 conditions, 6 unknowns)\n", ht, 3.0*DIG/6);
  print("    residuals = ", vector(3, i, precision(1.0*resid[i], 6)));
  print("    VERDICT: ", if(ht < 3.0*DIG/6 - 4, "possible rational fold -- investigate",
                            "NO rational fold; the height sits at the spurious level"));
  print();
  my(mt2 = matrix(2, 6, i, j, precision(ctri[j][4-i], DIG)));
  my(r2 = simrel(mt2, DIG - 3));
  my(h2 = log(1.0*vecmax(vector(6, i, abs(r2[1][i]))))/log(10));
  print("  2-condition LLL search (kill c_{0,2} and c_{0,1} only):");
  print("    c = ", r2[1]);
  printf("    height = 10^%.2f   spurious level 10^%.1f\n", h2, 2.0*DIG/6);
  print("    residuals = ", vector(2, i, precision(1.0*r2[2][i], 6)));
  print("    VERDICT: ", if(h2 < 2.0*DIG/6 - 4, "possible rational partial fold",
                            "NO rational partial fold"));
  print();
  my(mt1 = matrix(1, 6, i, j, precision(ctri[j][3], DIG)));
  my(r1 = simrel(mt1, DIG - 3));
  my(h1 = log(1.0*vecmax(vector(6, i, abs(r1[1][i]))))/log(10));
  print("  1-condition LLL search (kill c_{0,2} only) -- this is (b) again:");
  print("    c = ", r1[1], "   height = 10^", precision(h1,4), "  residual = ", precision(1.0*r1[2][1], 6));

  print();
  print("  FULL RELATION LATTICE of {1, xi_1..xi_5} (all LLL-reduced basis vectors,");
  print("  ordered by residual; genuine relations = tiny residual AND small height):");
  my(kk = 10^(DIG-3));
  my(lat = matrix(7, 6));
  for(j = 1, 6,
    for(i = 1, 6, lat[i,j] = if(i == j, 1, 0));
    lat[7,j] = round(kk * precision(ctri[j][3]/ctri[1][3], DIG)));
  my(red = lat * qflll(lat));
  for(j = 1, 6,
    my(cv = vector(6, i, red[i,j]), rs = 1.0*abs(red[7,j])/kk);
    printf("    c = %s   height 10^%.2f   residual 10^%.1f%s\n", cv,
      log(1.0*vecmax(vector(6,i,abs(cv[i]))))/log(10),
      if(rs == 0, -9999.0, log(rs)/log(10)),
      if(rs < 10.0^(-(DIG-12)), "   <== GENUINE", "")));
}
print();
print("---- explicit verification of the two relations at the fit's native precision:");
{ default(realprecision, PRECF);
  my(v1 = 2*xis[1] - 3*xis[2]);
  print("  R1:  2*xi_1 - 3*xi_2                                     = ", precision(v1, 6));
  my(v2 = -16875 - 281250*xis[1] + 2120000*xis[3] - 4218750*xis[4] + 3538944*xis[5]);
  print("  R2: -16875 - 281250 xi_1 + 2120000 xi_3 - 4218750 xi_4 + 3538944 xi_5 = ", precision(v2, 6));
  print("  factorisations of the R2 coefficients:");
  foreach([16875, 281250, 2120000, 4218750, 3538944], cf, print("     ", cf, " = ", factor(cf)));
  print("  (the xi_j carry ~56 correct digits, so a residual of that size means");
  print("   the relation holds to the full accuracy of the data.)");
}
print();
print("DONE.");
quit;
