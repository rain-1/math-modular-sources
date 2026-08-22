\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp
\\ Shared row library for the multi-prime two-row lattice scan.
\\ Exact rational arithmetic throughout.
\\ Names avoid PARI builtins (psi, M, Phi, S, cmp, L, W, O, I).

\\ ---- order-2 Zagier rows:  (n+1)^2 u_(n+1) = (a n^2 + a n + b) u_n - c n^2 u_(n-1)
row2(aa, bb, cc, NN) = {
  my(av = vector(NN+1), bv = vector(NN+1));
  av[1] = 1; av[2] = bb; bv[1] = 0; bv[2] = 1;
  for(n = 1, NN-1,
    av[n+2] = ((aa*n^2 + aa*n + bb)*av[n+1] - cc*n^2*av[n]) / (n+1)^2;
    bv[n+2] = ((aa*n^2 + aa*n + bb)*bv[n+1] - cc*n^2*bv[n]) / (n+1)^2);
  [av, bv];
};

\\ ---- order-3 AZ rows: (n+1)^3 u_(n+1) = (2n+1)(a n^2 + a n + b) u_n - c n^3 u_(n-1)
row3(aa, bb, cc, NN) = {
  my(av = vector(NN+1), bv = vector(NN+1));
  av[1] = 1; av[2] = bb; bv[1] = 0; bv[2] = 1;
  for(n = 1, NN-1,
    av[n+2] = ((2*n+1)*(aa*n^2 + aa*n + bb)*av[n+1] - cc*n^3*av[n]) / (n+1)^3;
    bv[n+2] = ((2*n+1)*(aa*n^2 + aa*n + bb)*bv[n+1] - cc*n^3*bv[n]) / (n+1)^3);
  [av, bv];
};

\\ ---- Cooper s_18: (n+1)^3 A = 2(2n+1)(7n^2+7n+3) A - 12 n (16 n^2 - 1) A,  A_1 = 6
cooper18(NN) = {
  my(av = vector(NN+1), bv = vector(NN+1));
  av[1] = 1; av[2] = 6; bv[1] = 0; bv[2] = 1;
  for(n = 1, NN-1,
    av[n+2] = (2*(2*n+1)*(7*n^2+7*n+3)*av[n+1] - 12*n*(16*n^2-1)*av[n]) / (n+1)^3;
    bv[n+2] = (2*(2*n+1)*(7*n^2+7*n+3)*bv[n+1] - 12*n*(16*n^2-1)*bv[n]) / (n+1)^3);
  [av, bv];
};

\\ ---- Zudilin Catalan row, repo normalisation Q_0=1, Q_1=7/4, P_0=0, P_1=13/8
zudrow(MM) = {
  my(qv = vector(MM+1), pv = vector(MM+1), aa, bb, cc);
  qv[1] = 1; qv[2] = 7/4; pv[1] = 0; pv[2] = 13/8;
  for(m = 1, MM-1,
    aa = (2*m+1)^2*(2*m+2)^2*(20*m^2 - 8*m + 1);
    bb = 3520*m^6 + 5632*m^5 + 2064*m^4 - 384*m^3 - 156*m^2 + 16*m + 7;
    cc = (2*m-1)^2*(2*m)^2*(20*m^2 + 32*m + 13);
    qv[m+2] = (bb*qv[m+1] + cc*qv[m]) / aa;
    pv[m+2] = (bb*pv[m+1] + cc*pv[m]) / aa);
  [qv, pv];
};

\\ ---- chi_(-3) well-poised hypergeometric row, ONE_CLASS_TWO_WORLDS Theorem 1.
\\ Numerator  Nab(t) = (2t+b+1) * prod_(j=1..a) (t-j+1)(t+b+j)
\\ Denominator Db(t) = prod_(j=0..b) (t+j+1/3)^2 (t+j+2/3)^2
\\ Returns [Qc, Pc] with  sum_(t>=0) Nab/Db = Qc * L(2,chi_(-3)) - Pc.
chi3row(aa, bb) = {
  my(nn = 2*bb+2, po = vector(nn), av = vector(nn), bv = vector(nn),
     sa1, sa2, sb1, qc, pc, t1, t2, off, base);
  for(j = 0, bb, po[j+1] = -(3*j+1)/3; po[bb+2+j] = -(3*j+2)/3);
  for(i = 1, nn,
    my(p = po[i], val = 2*p + bb + 1, dlg = 2/(2*p + bb + 1), pr = 1, sm = 0);
    for(j = 1, aa, my(f1 = p-j+1, f2 = p+bb+j); val *= f1*f2; dlg += 1/f1 + 1/f2);
    for(i2 = 1, nn, if(i2 != i, my(d = p - po[i2]); pr *= d^2; sm += 1/d));
    av[i] = val/pr;
    bv[i] = av[i]*(dlg - 2*sm));
  if(sum(i = 1, nn, bv[i]) != 0, error("sum B nonzero"));
  sa1 = sum(i = 1, bb+1, av[i]);
  sa2 = sum(i = bb+2, nn, av[i]);
  sb1 = sum(i = 1, bb+1, bv[i]);
  if(sa1 + sa2 != 0, error("zeta2 coefficient nonzero"));
  if(sb1 != 0, error("pi over sqrt3 coefficient nonzero"));
  qc = 9*sa1;
  pc = 0;
  for(cl = 0, 1,
    off = cl*(bb+1); base = if(cl == 0, 1/3, 2/3);
    t1 = 0; t2 = 0;
    for(j = 0, bb,
      my(i = off + j + 1);
      pc += av[i]*t2 + bv[i]*t1;
      t2 += 1/(j+base)^2; t1 += 1/(j+base)));
  [qc, pc];
};

\\ ---- helpers
denfac(x) = factor(denominator(x));
numfac(x) = factor(numerator(x));

\\ ---- NEW: conductor-6 well-poised row for L(2,chi_(-3)).
\\ Poles of order 2 at t = -(j+1/6) and t = -(j+5/6), j = 0..b; same numerator
\\ Nab(t) = (2t+b+1) prod_(j=1..a) (t-j+1)(t+b+j), antisymmetric under t -> -t-b-1,
\\ which swaps the two pole classes because 1/6 + 5/6 = 1.
\\ zeta(2,1/6) - zeta(2,5/6) = 36 * L(2, chi mod 6) = 36*(5/4)*L(2,chi_(-3)) = 45 L.
\\ zeta(2,1/6) + zeta(2,5/6) = 24 zeta(2) ; psi(1/6) - psi(5/6) = -pi*sqrt(3).
\\ Antisymmetry kills the zeta(2) and the pi*sqrt3 parts, leaving Qc*L - Pc.
\\ Returns [Qc, Pc] with sum_(t>=0) Nab/Db = Qc * L(2,chi_(-3)) - Pc.
chi6row(aa, bb) = {
  my(nn = 2*bb+2, po = vector(nn), av = vector(nn), bv = vector(nn),
     sa1, sa2, sb1, qc, pc, t1, t2, off, base);
  for(j = 0, bb, po[j+1] = -(6*j+1)/6; po[bb+2+j] = -(6*j+5)/6);
  for(i = 1, nn,
    my(p = po[i], val = 2*p + bb + 1, dlg = 2/(2*p + bb + 1), pr = 1, sm = 0);
    for(j = 1, aa, my(f1 = p-j+1, f2 = p+bb+j); val *= f1*f2; dlg += 1/f1 + 1/f2);
    for(i2 = 1, nn, if(i2 != i, my(d = p - po[i2]); pr *= d^2; sm += 1/d));
    av[i] = val/pr;
    bv[i] = av[i]*(dlg - 2*sm));
  if(sum(i = 1, nn, bv[i]) != 0, error("chi6: sum B nonzero"));
  sa1 = sum(i = 1, bb+1, av[i]);
  sa2 = sum(i = bb+2, nn, av[i]);
  sb1 = sum(i = 1, bb+1, bv[i]);
  if(sa1 + sa2 != 0, error("chi6: zeta2 coefficient nonzero"));
  if(sb1 != 0, error("chi6: pi sqrt3 coefficient nonzero"));
  qc = 45*sa1;
  pc = 0;
  for(cl = 0, 1,
    off = cl*(bb+1); base = if(cl == 0, 1/6, 5/6);
    t1 = 0; t2 = 0;
    for(j = 0, bb,
      my(i = off + j + 1);
      pc += av[i]*t2 + bv[i]*t1;
      t2 += 1/(j+base)^2; t1 += 1/(j+base)));
  [qc, pc];
};

\\ ---- generic conductor-F well-poised row.
\\ classes: list of residues r/FF with 0 < r < FF, gcd(r,FF)=1, paired r <-> FF-r.
\\ Returns [Qc, Pc] with sum_(t>=0) Nab/Db = Qc * Lval - Pc, where Lval is the
\\ L-value normalised by  const  supplied by the caller.
\\ conductor FF, odd primitive-or-imprimitive character with chi(r) = +1 for r in first
\\ half list, -1 for its mirror.  cst = FF^2 * prod_(p | FF, p not | cond) (1 - chi(p) p^-2).
cfrow(aa, bb, FF, rs, cst) = {
  my(hh = #rs, nn = hh*(bb+1), po = vector(nn), av = vector(nn), bv = vector(nn),
     sgn = vector(nn), tot = 0, qc, pc, t1, t2);
  for(c = 1, hh, for(j = 0, bb,
    po[(c-1)*(bb+1)+j+1] = -(j + rs[c]/FF);
    sgn[(c-1)*(bb+1)+j+1] = if(c <= hh/2, 1, -1)));
  for(i = 1, nn,
    my(p = po[i], val = 2*p + bb + 1, dlg = 2/(2*p + bb + 1), pr = 1, sm = 0);
    for(j = 1, aa, my(f1 = p-j+1, f2 = p+bb+j); val *= f1*f2; dlg += 1/f1 + 1/f2);
    for(i2 = 1, nn, if(i2 != i, my(d = p - po[i2]); pr *= d^2; sm += 1/d));
    av[i] = val/pr;
    bv[i] = av[i]*(dlg - 2*sm));
  if(sum(i = 1, nn, bv[i]) != 0, error("cf: sum B nonzero"));
  \\ chi-weighted A-sum must be the only survivor
  tot = sum(i = 1, nn, sgn[i]*av[i]);
  qc = cst * tot / 2;
  pc = 0;
  for(c = 1, hh,
    t1 = 0; t2 = 0;
    for(j = 0, bb,
      my(i = (c-1)*(bb+1)+j+1, base = rs[c]/FF);
      pc += av[i]*t2 + bv[i]*t1;
      t2 += 1/(j+base)^2; t1 += 1/(j+base)));
  [qc, pc, tot, sum(i = 1, nn, av[i])];
};
