\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/07_lattice.gp
\\ Step 3: the ACTUAL correlated congruence lattice, exact integers, in the full
\\ architecture of paper/sections/05_two_row.tex 5.1 :
\\   raw rows (Q_i, P_i) with linear form Q_i*Theta - P_i ;
\\   SQ = lcm_i den(Q_i) ,  al_i = SQ*Q_i in Z ;
\\   SP = lcm_i den(SQ*P_i) , Y_i = SQ*SP*P_i in Z , X_i = SP*al_i ;
\\   T   = prod over the bridge primes p of p^(min_(i<j) v_p(al_i Y_j - al_j Y_i)) ;
\\   Mm  = SP*T ;
\\   Kl  = c in Z^k : sum c_i al_i = 0 mod T and sum c_i Y_i = 0 mod Mm ;
\\   q(c) = (c.X)/Mm , p(c) = (c.Y)/Mm  are integers.
\\ Selection is INSIDE the prescribed anisotropic box (ZETA3_TWO_LATTICE.md 8:
\\ unconstrained LLL measures the lattice, not the construction).
default(parisizemax, 8000000000);
default(realprecision, 4000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");

THC = (zetahurwitz(2,1/3) - zetahurwitz(2,2/3))/9;   \\ L(2,chi_-3)
TH3 = zeta(3);
NMOD = 400;
MC = row2(10,3,9,NMOD); MF = row2(17,6,72,NMOD); MB = row2(9,3,27,NMOD);
RD = row3(10,4,64,NMOD); RT = row3(12,4,16,NMOD);

rawrow(tag, pp, qq, n) = {
  if(tag == "C",    return([  MC[1][n+1], 2*MC[2][n+1] ]));
  if(tag == "F",    return([5*MF[1][n+1], 8*MF[2][n+1] ]));
  if(tag == "B",    return([  MB[1][n+1], 2*MB[2][n+1] ]));
  if(tag == "Domb", return([7*RD[1][n+1],24*RD[2][n+1] ]));
  if(tag == "T",    return([7*RT[1][n+1],32*RT[2][n+1] ]));
  if(tag == "c3",   return(chi3row(pp*(n/qq), qq*(n/qq))));
  if(tag == "c6",   return(chi6row(pp*(n/qq), qq*(n/qq))));
  error("bad tag");
};

cfg(lab, rws, brid, nn, theta) = {
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
  for(i = 1, kk, aa[i] = log(abs(1.0*XX[i]))/nn; ee[i] = log(abs(XX[i]*theta - YY[i]))/nn);
  sg = log(1.0*cov)/nn;
  xs = vector(kk);
  if(kk == 2,
    xs[1] = 0.5*(sg + ee[2] - ee[1]); xs[2] = sg - xs[1],
    my(av = sum(i = 1, kk, ee[i])/kk);
    for(i = 1, kk, xs[i] = sg/kk + (av - ee[i])));
  \\ formula prediction (rank-2 theory; decayer = last row)
  ff = xs[1] + ee[1] - sg;
  hh = aa[kk] - xs[kk] + (sg - xs[kk]) - sg + xs[kk];  \\ = aa[kk] - xs[kk]
  hh = aa[kk] - xs[kk];
  dfor = if(hh <= 0, -99.0, 1 - ff/hh);
  bx = vector(kk, i, exp(xs[i]*nn));
  \\ LLL in the rescaled box
  my(sc = 2^3000, wts = vector(kk, i, sc \ ceil(bx[i])),
     bs = matrix(kk, kk, i, j, wts[i]*bas[i,j]), uu = qflll(bs));
  cand = bas*uu;
  for(j = 1, kk,
    my(v = cand[,j], ok = 1, qv, pv, lv);
    for(i = 1, kk, if(abs(1.0*v[i]) > bx[i], ok = 0));
    qv = sum(i = 1, kk, v[i]*XX[i])/mm;
    pv = sum(i = 1, kk, v[i]*YY[i])/mm;
    if(qv == 0, next);
    if(denominator(qv) != 1 || denominator(pv) != 1, print("  !! non-integral output"); next);
    lv = qv*theta - pv;
    if(lv == 0, next);
    if(ok, nin++;
      my(dl = 1 - log(abs(lv))/log(abs(1.0*qv)));
      if(dl > best, best = dl; bq = qv; bl = lv)));
  printf("  %-30s n=%3d logT/n=%7.3f logSP/n=%7.3f sig=%7.3f F=%8.4f H=%8.4f d_for=%7.4f | inbox=%d log|q|/n=%8.3f logform/n=%9.3f d_emp=%7.4f\n",
    lab, nn, log(1.0*tt)/nn, log(1.0*sp)/nn, sg, ff, hh, dfor, nin,
    if(bq == 0, 0.0, log(abs(1.0*bq))/nn), if(bq == 0, 0.0, log(abs(bl))/nn),
    if(bq == 0, -99.0, best));
};

print("### CONTROL: zeta(3) Domb x T at 2:3, single prime p=2  (target delta = 0.90095, F = +1.163)");
{for(j = 1, 8, my(nn = 10*j);
  cfg("Domb(2n) x T(3n) [p=2]", [["Domb",0,0,2*nn], ["T",0,0,3*nn]], [2], nn, TH3));}
print("");

print("### (A) two conductor-6 rows: the two-prime bridge, and what each prime is worth");
{for(j = 3, 12, my(nn = 4*j);
  cfg("c6_1/2 x c6_2/1 [none]", [["c6",2,1,nn], ["c6",1,2,nn]], [],     nn, THC);
  cfg("c6_1/2 x c6_2/1 [p=3] ", [["c6",2,1,nn], ["c6",1,2,nn]], [3],    nn, THC);
  cfg("c6_1/2 x c6_2/1 [p=2] ", [["c6",2,1,nn], ["c6",1,2,nn]], [2],    nn, THC);
  cfg("c6_1/2 x c6_2/1 [BOTH]", [["c6",2,1,nn], ["c6",1,2,nn]], [2,3],  nn, THC);
  print(""));}

print("### (B) benchmark: conductor-3 decayer x modular engine C, single prime p=3");
{for(j = 3, 12, my(nn = 4*j);
  cfg("c3_1/2 x C [p=3]", [["C",0,0,4*nn], ["c3",1,2,nn]], [3], nn, THC));}
print("");

print("### (C) conductor-6 decayer x modular F: two rows, two primes");
{for(j = 3, 12, my(nn = 4*j);
  cfg("c6_1/2 x F [p=3] ", [["F",0,0,nn], ["c6",1,2,nn]], [3],   nn, THC);
  cfg("c6_1/2 x F [BOTH]", [["F",0,0,nn], ["c6",1,2,nn]], [2,3], nn, THC));}
print("");

print("### (D) THREE rows, TWO primes: C bridges p=3, F bridges p=2, c6 is the decayer");
{for(j = 3, 12, my(nn = 4*j);
  cfg("C + F + c6_1/2 [p=3] ", [["C",0,0,nn], ["F",0,0,nn], ["c6",1,2,nn]], [3],   nn, THC);
  cfg("C + F + c6_1/2 [BOTH]", [["C",0,0,nn], ["F",0,0,nn], ["c6",1,2,nn]], [2,3], nn, THC));}
quit
