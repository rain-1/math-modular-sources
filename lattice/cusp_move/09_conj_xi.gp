/* 09_conj_xi.gp -- the Apery limit of the Galois-conjugate placement of the two
   positive-score rows (Apery's zeta(2) row = Zagier D, and Beukers = sqrt(Apery)).
   Does the second realisation approximate the SAME number?                     */
read("lib.gp");
default(realprecision, 220);

{ conjxi(nm, aa, dd, p1, p0, r1, r2, nmax) =
  my(sq = sqrt(aa^2-4*dd), lam, mu, al, pp, qq, av, bv, x1, x2);
  lam = (aa+sq)/2; mu = (aa-sq)/2;            \\ move by mu (the small root)
  al = 1-r1;
  \\ placement at v=mu : lam_move = mu, mu_move = lam
  pp = (lam-2*mu)*nv^2 + (p1-3*mu+2*mu*r1)*nv + (p0-mu+mu*r1);
  qq = (mu^2-mu*lam)*nv^2
     + (-2*mu^2*r1+mu^2+mu*lam*r1-mu*lam*r2+mu*lam-mu*p1)*nv
     + (mu^2*r1^2-mu^2*r1+mu*lam*r1*r2-mu*lam*r1+mu*p1*r1);
  av = seqA([pp,qq], nmax); bv = seqB([pp,qq], nmax);
  x1 = bv[nmax+1]/av[nmax+1];
  x2 = bv[nmax]/av[nmax];
  print("  ", nm, " conjugate placement (v=mu):");
  print("     P# = ", pp);
  print("     Q# = ", qq);
  print("     xi# (n=", nmax, ")  = ", x1);
  print("     Cauchy check |xi#(n)-xi#(n-1)| = ", abs(x1-x2));
  x1; }

print("=== Apery's zeta(2) row, Zagier D (11,3,-1) ===");
xD  = 0; { my(av = seqA([11*nv^2+11*nv+3, -nv^2], 400), bv = seqB([11*nv^2+11*nv+3, -nv^2], 400));
  xD = bv[401]/av[401]*1.0; }
print("  original xi = ", xD, "   zeta(2)/5 = ", zeta(2)/5);
xDc = conjxi("ZagierD", 11, -1, 11, 3, 0, 0, 260);
print("  lindep [xi#, xi, 1] = ", lindep([xDc, xD, 1], 60));
print("  lindep [xi#, zeta(2), 1] = ", lindep([xDc, zeta(2), 1], 60));
print("  xi#/xi = ", xDc/xD, "  bestappr = ", bestappr(xDc/xD, 10^12));

print("");
print("=== Beukers (136,10,4) = sqrt(Apery) ===");
xB = 0; { my(av = seqA([136*nv^2+68*nv+10, 4*(2*nv-1)^2], 400),
             bv = seqB([136*nv^2+68*nv+10, 4*(2*nv-1)^2], 400));
  xB = bv[401]/av[401]*1.0; }
print("  original xi = ", xB);
xBc = conjxi("Beukers", 136, 16, 68, 10, 1/2, 1/2, 260);
print("  lindep [xi#, xi, 1] = ", lindep([xBc, xB, 1], 60));
print("  xi#/xi = ", xBc/xB, "  bestappr = ", bestappr(xBc/xB, 10^12));
print("  lindep [xi#, xi, zeta(3), zeta(2), 1] = ", lindep([xBc, xB, zeta(3), zeta(2), 1], 50));
quit;
