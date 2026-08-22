\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/06_design_full.gp
\\ Unified design scan for the class L(2,chi_(-3)) at the two primes 2 and 3.
\\ Every row is normalised to  Q_n * Theta - P_n  with Theta = L(2,chi_(-3)),
\\ so all rows share ONE p-adic constant at each slope prime:
\\   Lam_3 = zeta_3(2) ,  Lam_2 = (4/5) L_2(2,chi_12) .
\\ Master formula (paper/sections/05_two_row.tex eq. master), decayer at index n,
\\ engine at index m ~ r n :
\\   Kr = log lcm( Den_eng(m), Den_dec(n) ) / n         (exact lcm of the actual integers)
\\   G  = sum_p min( sig_p_eng * m/n , sig_p_dec ) log p
\\   F  = 1/2 [ Kr + (m/n) logrho_eng + loglam_dec - G ]
\\   H  = F + logLam_dec - loglam_dec ,   delta = 1 - F/H .
default(parisizemax, 8000000000);
default(realprecision, 2500);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");

Lc = (zetahurwitz(2,1/3) - zetahurwitz(2,2/3))/9;
NMOD = 500;
MB = row2( 9,3,27,NMOD); MC = row2(10,3,9,NMOD);
MF = row2(17,6,72,NMOD); MS = cooper18(NMOD);

\\ ROW SPEC: [label, kind, p, q, sig2, sig3]  ; kind 0 = modular, 3 = chi3, 6 = chi6
{ROWS = [
  ["C",   0, 0, 1, 0.00, 2.00],
  ["B",   0, 1, 1, 0.00, 3.00],
  ["F",   0, 2, 1, 3.00, 2.00],
  ["s18", 0, 3, 1, 0.00, 1.00],
  ["c3_2/1", 3, 2, 1, 0.00, 3.00],
  ["c3_5/3", 3, 5, 3, 0.00, 3.00],
  ["c3_3/2", 3, 3, 2, 0.00, 3.00],
  ["c3_1/1", 3, 1, 1, 0.00, 3.00],
  ["c3_1/2", 3, 1, 2, 0.00, 3.00],
  ["c3_1/3", 3, 1, 3, 0.00, 3.00],
  ["c6_2/1", 6, 2, 1, 7.95, 2.90],
  ["c6_3/2", 6, 3, 2, 7.95, 2.90],
  ["c6_1/1", 6, 1, 1, 7.95, 2.90],
  ["c6_1/2", 6, 1, 2, 7.95, 2.90],
  ["c6_1/3", 6, 1, 3, 7.95, 2.90],
  ["c6_1/4", 6, 1, 4, 7.95, 2.90]
];}

\\ (Q_n, P_n) at index n, in the common normalisation
qp(i, n) = {
  my(s = ROWS[i], kd = s[2]);
  if(kd == 0,
    if(s[3] == 0, return([  MC[1][n+1], 2*MC[2][n+1] ]));
    if(s[3] == 1, return([  MB[1][n+1], 2*MB[2][n+1] ]));
    if(s[3] == 2, return([5*MF[1][n+1], 8*MF[2][n+1] ]));
    return([ MS[1][n+1], 2*MS[2][n+1] ]));
  my(qq = s[4], kk = n/qq);
  if(kd == 3, chi3row(s[3]*kk, qq*kk), chi6row(s[3]*kk, qq*kk));
};
\\ legal indices for row i are multiples of step(i)
step(i) = if(ROWS[i][2] == 0, 1, ROWS[i][4]);

\\ cached rows
QCACHE = Map();
qpc(i, n) = {
  my(key = Str(i, ":", n), v);
  if(mapisdefined(QCACHE, key, &v), return(v));
  v = qp(i, n); mapput(QCACHE, key, v); v;
};

\\ decayer profile at index n
dprofc(i, n) = {
  my(r = qpc(i, n));
  [ log(abs(r[1]))/n, log(abs(r[1]*Lc - r[2]))/n, lcm(denominator(r[1]), denominator(r[2])) ];
};

pair(id, nd, ie, verbose) = {
  my(pr = dprofc(id, nd), lLam = pr[1], llam = pr[2], Dd = pr[3],
     s2d = ROWS[id][5], s3d = ROWS[id][6],
     s2e = ROWS[ie][5], s3e = ROWS[ie][6], st = step(ie),
     best = vector(4), bestr = vector(4));
  for(j = 1, 14,
    my(rr = j*0.25, me = st*max(1, round(rr*nd/st)), rreal, kk, lrho, De, g2, g3, gs, ff, hh, dl, qe);
    if(me > NMOD && ROWS[ie][2] == 0, next);
    rreal = 1.0*me/nd;
    qe = qpc(ie, me);
    lrho = log(abs(qe[1]*Lc - qe[2]))/me;
    De = lcm(denominator(qe[1]), denominator(qe[2]));
    kk = log(1.0*lcm(De, Dd))/nd;
    g3 = min(s3e*rreal, s3d)*log(3);
    g2 = min(s2e*rreal, s2d)*log(2);
    gs = [0.0, g3, g2, g2+g3];
    for(w = 1, 4,
      ff = 0.5*(kk + rreal*lrho + llam - gs[w]);
      hh = ff + lLam - llam;
      dl = if(hh <= 0, -1.0, 1 - ff/hh);
      if(dl > best[w], best[w] = dl; bestr[w] = rreal)));
  if(verbose,
    printf("  %-8s(n=%3d) x %-8s : none %7.4f | p3 %7.4f (r=%4.2f) | p2 %7.4f | BOTH %7.4f (r=%4.2f)\n",
      ROWS[id][1], nd, ROWS[ie][1], best[1], best[2], bestr[2], best[3], best[4], bestr[4]));
  [best, bestr];
};

\\ decayer index choices
{NIDX = [0,0,0,0, 60, 120, 120, 60, 120, 120, 60, 80, 60, 80, 84, 96];}

runall() = {
  my(bd = 0, blab = "", bd1 = 0, blab1 = "");
  print("### master formula, exact lcm denominators");
  for(id = 5, #ROWS,
    print("");
    for(ie = 1, #ROWS,
      if(ie == id, next);
      my(res = pair(id, NIDX[id], ie, 1));
      if(res[1][4] > bd, bd = res[1][4]; blab = Str(ROWS[id][1], " x ", ROWS[ie][1], " r=", res[2][4]));
      if(res[1][2] > bd1, bd1 = res[1][2]; blab1 = Str(ROWS[id][1], " x ", ROWS[ie][1], " r=", res[2][2]))));
  print("");
  printf("BEST single-prime (p=3 only) : %.4f   at %s\n", bd1, blab1);
  printf("BEST two-prime   (p=2 and 3) : %.4f   at %s\n", bd,  blab);
};

default(realprecision, 38);
runall();
quit
