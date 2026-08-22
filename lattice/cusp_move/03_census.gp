/* 03_census.gp -- the cusp-move orbit census.
   For every row of the corpus whose characteristic roots are rational we compute the
   full orbit closure under the two cusp moves (and both admissible gauges alpha),
   and for each orbit member: minimal integral rescaling, (a,b,d)-data, lambda_1,
   lambda_2, sharp k, score, archimedean Apery limit and p-adic diagnostics.        */
read("lib.gp"); read("rows.gp");
default(parisize, 2^31);
default(realprecision, 260);
NINT = 200;     \\ integrality / k depth
NARC = 420;     \\ depth for the archimedean limit

{
qroots2(aa, bb) =
  my(ds = aa^2-4*bb, sq);
  if(issquare(ds, &sq), [(aa+sq)/2, (aa-sq)/2], []);
}
{
rrootsQ(rw) =
  my(q2 = polcoef(rw[2],2,nv), q1 = polcoef(rw[2],1,nv), q0 = polcoef(rw[2],0,nv));
  if(q2 == 0, return([]));
  qroots2(-q1/q2, q0/q2);
}
{
moverow(rw, lam, mu, r1, r2) =
  my(p1 = polcoef(rw[1],1,nv), p0 = polcoef(rw[1],0,nv), pp, qq);
  pp = (mu-2*lam)*nv^2 + (p1-3*lam+2*lam*r1)*nv + (p0-lam+lam*r1);
  qq = (lam^2-lam*mu)*nv^2
     + (-2*lam^2*r1+lam^2+lam*mu*r1-lam*mu*r2+lam*mu-lam*p1)*nv
     + (lam^2*r1^2-lam^2*r1+lam*mu*r1*r2-lam*mu*r1+lam*p1*r1);
  [pp, qq];
}
/* scaling-invariant key: (a^2/d, p1/a, p0/a, q1/d, q0/d) */
{
rowkey(rw) =
  my(aa = polcoef(rw[1],2,nv), dd = polcoef(rw[2],2,nv));
  if(aa == 0,
     return([0, dd/dd, polcoef(rw[1],1,nv)^2/dd, polcoef(rw[1],0,nv)^2/dd,
             polcoef(rw[2],1,nv)/dd, polcoef(rw[2],0,nv)/dd]));
  [aa^2/dd, polcoef(rw[1],1,nv)/aa, polcoef(rw[1],0,nv)/aa,
   polcoef(rw[2],1,nv)/dd, polcoef(rw[2],0,nv)/dd];
}
/* all moves out of a row (rational-root case) */
{
allmoves(rw) =
  my(rt = qroots2(polcoef(rw[1],2,nv), polcoef(rw[2],2,nv)), rq = rrootsQ(rw),
     out = List(), nr = 0, jm = 2);
  if(#rt == 0 || #rq == 0, return(Vec(out)));
  if(rq[1] == rq[2], jm = 1);
  for(i = 1, 2,
    for(j = 1, jm,
      nr = moverow(rw, rt[i], rt[3-i], rq[j], rq[3-j]);
      if(polcoef(nr[2],2,nv) != 0, listput(out, nr))));
  Vec(out);
}
/* orbit closure by BFS, up to the scaling equivalence */
{
orbit(rw, maxsz = 40) =
  my(keys = List(), reps = List(), front = List([rw]), nf, ms, kk, seen);
  listput(keys, rowkey(rw)); listput(reps, rw);
  while(#front > 0 && #reps < maxsz,
    nf = List();
    for(i = 1, #front,
      ms = allmoves(front[i]);
      for(j = 1, #ms,
        kk = rowkey(ms[j]); seen = 0;
        for(l = 1, #keys, if(keys[l] == kk, seen = 1; break));
        if(!seen, listput(keys, kk); listput(reps, ms[j]); listput(nf, ms[j]))));
    front = nf);
  Vec(reps);
}
/* minimal integral rescaling of a row (search c over positive rationals p/q, small) */
{
intrescale(rw, nmax) =
  my(v = seqA(rw, min(nmax, 60)), c = minscale(v, 10^6));
  c;
}
{
padicline(rw, av, bv, nmax) =
  my(dd = polcoef(rw[2],2,nv), f, out = "");
  if(dd == 0, return("--"));
  f = factor(abs(numerator(dd)*denominator(dd)));
  for(i = 1, #f~,
    my(p = f[i,1], v1, v2);
    if(av[nmax+1] == 0 || av[nmax] == 0 || av[nmax\2+1] == 0 || av[nmax\2] == 0,
       out = Str(out, " p=", p, ":(a_n vanishes)"); next);
    v1 = bv[nmax+1]/av[nmax+1] - bv[nmax]/av[nmax];
    v2 = bv[nmax\2+1]/av[nmax\2+1] - bv[nmax\2]/av[nmax\2];
    out = Str(out, " p=", p, ":v(", nmax\2, ")=", if(v2==0,"inf",valuation(v2,p)),
              ",v(", nmax, ")=", if(v1==0,"inf",valuation(v1,p))));
  out;
}
{
describe(rw, tagname) =
  my(c, rw2, av, bv, kk, rt, l1, l2, sc, xi, i);
  c = intrescale(rw, 60);
  if(c == 0, print("    [no integral rescaling found]"); return(0));
  rw2 = [c*rw[1], c^2*rw[2]];
  av = seqA(rw2, NINT); bv = seqB(rw2, NINT);
  my(intok = 1);
  for(n = 0, NINT, if(denominator(av[n+1]) != 1, intok = 0; break));
  kk = denomexp(bv, NINT, 6);
  rt = qroots2(polcoef(rw2[1],2,nv), polcoef(rw2[2],2,nv));
  l1 = rt[1]; l2 = rt[2];
  if(abs(l1) < abs(l2), my(tm=l1); l1=l2; l2=tm);
  sc = log(1/abs(l2*1.0)) - kk;
  my(aa2 = seqA(rw2, NARC), bb2 = seqB(rw2, NARC));
  xi = if(abs(l2) < abs(l1) && aa2[NARC+1] != 0, bb2[NARC+1]/aa2[NARC+1]*1.0, 0);
  print("  * ", tagname, "  scale c=", c);
  print("      P = ", rw2[1]);
  print("      Q = ", rw2[2]);
  print("      a=", polcoef(rw2[1],2,nv), " b=", subst(rw2[1],nv,0),
        " d=", polcoef(rw2[2],2,nv), "   I=a^2/d=",
        polcoef(rw2[1],2,nv)^2/polcoef(rw2[2],2,nv));
  print("      integral to n=", NINT, ": ", if(intok,"YES","NO"),
        "   k=", kk, " (sharp)   lam1=", l1, " lam2=", l2,
        "   score=", strprintf("%.5f", sc));
  print("      xi (n=", NARC, ") = ", strprintf("%.40f", xi));
  print("      p-adic:", padicline(rw2, aa2, bb2, 300));
  [rw2, kk, l1, l2, sc, xi];
}
{
for(i = 1, #corpus,
  my(rw = corpus[i], rt, orb);
  rt = qroots2(polcoef(rw[1],2,nv), polcoef(rw[2],2,nv));
  if(#rt == 0, next);
  print("");
  print("################ ", rw[3], "   ", rw[4]);
  orb = orbit([rw[1], rw[2]]);
  print("  orbit size (mod scaling) = ", #orb);
  for(j = 1, #orb, describe(orb[j], Str(rw[3], ".", j)));
);
}
quit;
