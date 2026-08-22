\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/09_control.gp
\\ The rational-surrogate control of CATALAN_AUDIT.md section 4(a) /
\\ paper/sections/05_two_row.tex Remark "the hypothesis F>0 is not decorative".
\\ The congruence lattice never uses Theta, so replacing Theta by a rational
\\ surrogate Th* of astronomically large denominator must reproduce the SAME
\\ numerics.  If it does, the delta > 1 readings carry no information about Theta.
default(parisizemax, 8000000000);
default(realprecision, 4000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");

THC = (zetahurwitz(2,1/3) - zetahurwitz(2,2/3))/9;
THS = bestappr(THC, 10^1200);          \\ rational surrogate, den ~ 10^1200

rawrow(tag, pp, qq, n) = if(tag == "c6", chi6row(pp*(n/qq), qq*(n/qq)), error("tag"));

cfg(lab, rws, brid, nn) = {
  my(kk = #rws, Qr = vector(kk), Pr = vector(kk), sq = 1, sp = 1,
     al = vector(kk), YY = vector(kk), XX = vector(kk), tt = 1, mm, v2h = 0, v3h = 0,
     amat, ker, bas, cov, ee = vector(kk), aa = vector(kk), sg, xs, bx,
     ff, hh, best = -100, bq = 0, bl = 0, bls = 0, nin = 0, cand);
  for(i = 1, kk,
    my(r = rawrow(rws[i][1], rws[i][2], rws[i][3], rws[i][4]));
    Qr[i] = r[1]; Pr[i] = r[2]; sq = lcm(sq, denominator(r[1])));
  for(i = 1, kk, al[i] = sq*Qr[i]; sp = lcm(sp, denominator(sq*Pr[i])));
  for(i = 1, kk, YY[i] = sq*sp*Pr[i]; XX[i] = sp*al[i]);
  v2h = valuation(al[1]*YY[2] - al[2]*YY[1], 2);
  v3h = valuation(al[1]*YY[2] - al[2]*YY[1], 3);
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
  xs = vector(2); xs[1] = 0.5*(sg + ee[2] - ee[1]); xs[2] = sg - xs[1];
  ff = xs[1] + ee[1] - sg;  hh = aa[kk] - xs[kk];
  bx = vector(kk, i, exp(xs[i]*nn));
  my(sc = 2^4000, wts = vector(kk, i, sc \ ceil(bx[i])),
     bs = matrix(kk, kk, i, j, wts[i]*bas[i,j]), uu = qflll(bs));
  cand = bas*uu;
  for(j = 1, kk,
    my(v = cand[,j], ok = 1, qv, pv, lv);
    for(i = 1, kk, if(abs(1.0*v[i]) > bx[i], ok = 0));
    qv = sum(i = 1, kk, v[i]*XX[i])/mm; pv = sum(i = 1, kk, v[i]*YY[i])/mm;
    if(qv == 0 || denominator(qv) != 1 || denominator(pv) != 1, next);
    lv = qv*THC - pv; if(lv == 0, next);
    if(ok, nin++; my(dl = 1 - log(abs(lv))/log(abs(1.0*qv)));
      if(dl > best, best = dl; bq = qv; bl = lv; bls = qv*THS - pv)));
  if(bq == 0, printf("  %-26s n=%3d  NO vector in box\n", lab, nn); return(0));
  printf("  %-26s n=%3d v2(h)/n=%6.3f v3(h)/n=%6.3f sig=%7.3f F=%8.4f | log|q|/n=%7.3f  log|qT-p|/n=%8.4f  d_emp=%7.4f || SURROGATE log|qT*-p|/n=%8.4f  d*=%7.4f\n",
    lab, nn, 1.0*v2h/nn, 1.0*v3h/nn, sg, ff,
    log(abs(1.0*bq))/nn, log(abs(bl))/nn, best,
    log(abs(bls))/nn, 1 - log(abs(bls))/log(abs(1.0*bq)));
};

print("### F<0 candidates: two conductor-6 rows, bridge at BOTH primes");
print("    surrogate Th* = bestappr(Theta, 10^1200), so 1/den(Th*) ~ 1e-1200");
{for(j = 2, 9, my(nn = 12*j);
  cfg("c6_1/3 x c6_1/4 [BOTH]", [["c6",1,3,nn], ["c6",1,4,nn]], [2,3], nn));}
print("");
{for(j = 2, 9, my(nn = 12*j);
  cfg("c6_1/2 x c6_1/3 [BOTH]", [["c6",1,2,nn], ["c6",1,3,nn]], [2,3], nn));}
print("");
print("### same pair, p=3 only (single prime), for comparison");
{for(j = 2, 9, my(nn = 12*j);
  cfg("c6_1/3 x c6_1/4 [p=3] ", [["c6",1,3,nn], ["c6",1,4,nn]], [3], nn));}
print("");
{for(j = 2, 9, my(nn = 12*j);
  cfg("c6_1/3 x c6_1/4 [none]", [["c6",1,3,nn], ["c6",1,4,nn]], [], nn));}
quit
