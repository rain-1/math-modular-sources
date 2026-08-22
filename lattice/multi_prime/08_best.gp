\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/08_best.gp
\\ (i) high-precision two-prime alignment certificate for the conductor-6 rows
\\ (ii) empirical sweep over conductor-6 (decayer, engine) pairs, bridge BOTH vs single
default(parisizemax, 8000000000);
default(realprecision, 4000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");

THC = (zetahurwitz(2,1/3) - zetahurwitz(2,2/3))/9;
NMOD = 700;
MC = row2(10,3,9,NMOD); MF = row2(17,6,72,NMOD); MB = row2(9,3,27,NMOD);
xi2F = MF[2][NMOD+1]/MF[1][NMOD+1];
xi3C = MC[2][NMOD+1]/MC[1][NMOD+1];
xi3B = MB[2][NMOD+1]/MB[1][NMOD+1];

cert() = {
  print("### two-prime alignment certificate for the conductor-6 rows");
  print("    predicted  xi_3 = zeta_3(2) = 2 xi_3(C) = 2 xi_3(B) = (8/5) xi_3(F)");
  print("    predicted  xi_2 = (4/5) L_2(2,chi_12) = (8/5) xi_2(F)");
  print("    reference precision: v_3 Cauchy(C) = ", valuation(MC[2][NMOD+1]/MC[1][NMOD+1] - MC[2][NMOD]/MC[1][NMOD],3),
        " , v_2 Cauchy(F) = ", valuation(MF[2][NMOD+1]/MF[1][NMOD+1] - MF[2][NMOD]/MF[1][NMOD],2));
  for(j = 1, 4,
    my(ab = [[2,1,120],[3,2,60],[1,2,60],[1,1,120]][j], r, xk);
    r = chi6row(ab[1]*ab[3], ab[2]*ab[3]); xk = r[2]/r[1];
    printf("  alpha=%d/%d  b=%3d :  v_3(xi - 2 xi3C) = %5d   v_3(xi - 2 xi3B) = %5d   v_2(5 xi - 8 xi2F) = %5d\n",
      ab[1], ab[2], ab[2]*ab[3], valuation(xk - 2*xi3C,3), valuation(xk - 2*xi3B,3), valuation(5*xk - 8*xi2F,2));
    printf("        controls: v_3(xi - xi3C) = %4d   v_3(xi - 4 xi3C) = %4d   v_2(xi - xi2F) = %4d   v_2(xi - 2 xi2F) = %4d\n",
      valuation(xk - xi3C,3), valuation(xk - 4*xi3C,3), valuation(xk - xi2F,2), valuation(xk - 2*xi2F,2)));
};

rawrow(tag, pp, qq, n) = {
  if(tag == "C", return([  MC[1][n+1], 2*MC[2][n+1] ]));
  if(tag == "F", return([5*MF[1][n+1], 8*MF[2][n+1] ]));
  if(tag == "c3", return(chi3row(pp*(n/qq), qq*(n/qq))));
  if(tag == "c6", return(chi6row(pp*(n/qq), qq*(n/qq))));
  error("bad tag");
};

cfg(lab, rws, brid, nn) = {
  my(kk = #rws, Qr = vector(kk), Pr = vector(kk), sq = 1, sp = 1,
     al = vector(kk), YY = vector(kk), XX = vector(kk), tt = 1, mm,
     amat, ker, bas, cov, ee = vector(kk), aa = vector(kk), sg, xs, bx,
     ff, hh, dfor, best = -100, bq = 0, bl = 0, nin = 0, cand);
  for(i = 1, kk,
    my(r = rawrow(rws[i][1], rws[i][2], rws[i][3], rws[i][4]));
    Qr[i] = r[1]; Pr[i] = r[2]; sq = lcm(sq, denominator(r[1])));
  for(i = 1, kk, al[i] = sq*Qr[i]; sp = lcm(sp, denominator(sq*Pr[i])));
  for(i = 1, kk, YY[i] = sq*sp*Pr[i]; XX[i] = sp*al[i]);
  for(b = 1, #brid,
    my(p = brid[b], vm = -1);
    for(i = 1, kk-1, for(j = i+1, kk,
      my(v = valuation(al[i]*YY[j] - al[j]*YY[i], p));
      if(vm < 0 || v < vm, vm = v)));
    if(vm > 0, tt *= p^vm));
  mm = sp*tt;
  amat = matrix(2, kk, i, j, if(i == 1, (mm/tt)*al[j], YY[j]));
  ker = matkermod(amat, mm);
  bas = mathnf(matconcat([ker, mm*matid(kk)]));
  cov = abs(matdet(bas));
  for(i = 1, kk, aa[i] = log(abs(1.0*XX[i]))/nn; ee[i] = log(abs(XX[i]*THC - YY[i]))/nn);
  sg = log(1.0*cov)/nn;
  xs = vector(kk);
  if(kk == 2, xs[1] = 0.5*(sg + ee[2] - ee[1]); xs[2] = sg - xs[1],
     my(av = sum(i=1,kk,ee[i])/kk); for(i = 1, kk, xs[i] = sg/kk + (av - ee[i])));
  ff = xs[1] + ee[1] - sg;  hh = aa[kk] - xs[kk];
  dfor = if(hh <= 0, -99.0, 1 - ff/hh);
  bx = vector(kk, i, exp(xs[i]*nn));
  my(sc = 2^3000, wts = vector(kk, i, sc \ ceil(bx[i])),
     bs = matrix(kk, kk, i, j, wts[i]*bas[i,j]), uu = qflll(bs));
  cand = bas*uu;
  for(j = 1, kk,
    my(v = cand[,j], ok = 1, qv, pv, lv);
    for(i = 1, kk, if(abs(1.0*v[i]) > bx[i], ok = 0));
    qv = sum(i = 1, kk, v[i]*XX[i])/mm; pv = sum(i = 1, kk, v[i]*YY[i])/mm;
    if(qv == 0 || denominator(qv) != 1 || denominator(pv) != 1, next);
    lv = qv*THC - pv; if(lv == 0, next);
    if(ok, nin++; my(dl = 1 - log(abs(lv))/log(abs(1.0*qv)));
      if(dl > best, best = dl; bq = qv; bl = lv)));
  printf("  %-30s n=%3d logT/n=%7.3f sig=%7.3f F=%8.3f H=%8.3f d_for=%8.4f | d_emp=%7.4f  logform/n=%8.3f\n",
    lab, nn, log(1.0*tt)/nn, sg, ff, hh, dfor, if(bq==0,-99.0,best), if(bq==0,0.0,log(abs(bl))/nn));
  if(bq == 0, -99.0, best);
};

sweep() = {
  print("");
  print("### empirical sweep: conductor-6 (engine, decayer) pairs, n = 48");
  my(al6 = [[2,1],[3,2],[5,3],[1,1],[2,3],[1,2],[1,3],[1,4]], nn = 48, bb = -100, bl = "");
  for(ie = 1, #al6, for(id = 1, #al6,
    if(ie == id, next);
    my(la = Str("c6_",al6[ie][1],"/",al6[ie][2]," x c6_",al6[id][1],"/",al6[id][2]),
       d1, d2);
    d2 = cfg(Str(la," [p=3]"), [["c6",al6[ie][1],al6[ie][2],nn], ["c6",al6[id][1],al6[id][2],nn]], [3], nn);
    d1 = cfg(Str(la," [BOTH]"), [["c6",al6[ie][1],al6[ie][2],nn], ["c6",al6[id][1],al6[id][2],nn]], [2,3], nn);
    if(d1 > bb, bb = d1; bl = la)));
  print("  best two-prime d_emp = ", bb, "  at ", bl);
};

bench() = {
  print("");
  print("### benchmark: conductor-3 decayer x modular C, single prime p=3, several sampling ratios");
  for(j = 1, 6,
    my(rr = [1,2,3,4,5,6][j], nn = 36);
    cfg(Str("c3_1/2(n) x C(",rr,"n) [p=3]"), [["C",0,0,rr*nn], ["c3",1,2,nn]], [3], nn));
  for(j = 1, 6,
    my(rr = [1,2,3,4,5,6][j], nn = 36);
    cfg(Str("c3_5/3(n) x C(",rr,"n) [p=3]"), [["C",0,0,rr*nn], ["c3",5,3,nn]], [3], nn));
};

cert();
bench();
sweep();
quit
