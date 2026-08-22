/* 02_orbits.gp -- cusp-move orbits of every row in the corpus.
   Uses the closed formula derived in 01_general_move.py and cross-checks it
   against the coefficient-level (binomial) transform.                     */
read("lib.gp"); read("rows.gp");
default(parisize, 2^30);
NCHK = 40;      \\ cross-check depth for the binomial transform
NINT = 200;     \\ integrality / k depth

/* roots of a quadratic  x^2 - A x + B ; exact, in Q or in Q(sqrt(disc)) */
{
qroots2(aa, bb) =
  my(ds = aa^2-4*bb, sq);
  if(issquare(ds, &sq), return([(aa+sq)/2, (aa-sq)/2]));
  my(y = Mod('y, 'y^2 - aa*'y + bb));
  [y, aa - y];
}
/* roots of Q(n) = q2 n^2 + q1 n + q0 */
{
rrootsQ(rw) =
  my(q2 = polcoef(rw[2],2,nv), q1 = polcoef(rw[2],1,nv), q0 = polcoef(rw[2],0,nv));
  if(q2 == 0, error("deg Q < 2"));
  qroots2(-q1/q2, q0/q2);
}
/* closed-form cusp move; al = 1-r1 branch */
{
moverow(rw, lam, mu, r1, r2) =
  my(p1 = polcoef(rw[1],1,nv), p0 = polcoef(rw[1],0,nv), pp, qq);
  pp = (mu-2*lam)*nv^2 + (p1-3*lam+2*lam*r1)*nv + (p0-lam+lam*r1);
  qq = (lam^2-lam*mu)*nv^2
     + (-2*lam^2*r1+lam^2+lam*mu*r1-lam*mu*r2+lam*mu-lam*p1)*nv
     + (lam^2*r1^2-lam^2*r1+lam*mu*r1*r2-lam*mu*r1+lam*p1*r1);
  [pp, qq];
}
/* companion inhomogeneity  g(t) = t (1-lam t)^{-al}:  g_k = lam^{k-1} binom(al+k-2,k-1) */
ginh(lam, al) = (k) -> lam^(k-1)*binomial(al+k-2, k-1);

{
report(rw, tag) =
  my(rt, rq, lam, mu, r1, r2, nrw, v, w, ok, c, iv, bv, kk, sc, nm);
  nm = rw[3];
  rt = qroots2(rowa(rw), rowd(rw));
  rq = rrootsQ(rw);
  print("");
  print("=== ", nm, "  ", rw[4]);
  print("  P = ", rw[1], "   Q = ", rw[2]);
  print("  char roots lam,mu = ", rt[1], ", ", rt[2],
        "     Q-roots r1,r2 = ", rq[1], ", ", rq[2]);
  v = seqA(rw, NCHK);
  for(i = 1, 2,
    for(j = 1, 2,
      my(lm, ml, ra, rb, al);
      lm = rt[i]; ml = rt[3-i]; ra = rq[j]; rb = rq[3-j]; al = 1-ra;
      if(j == 2 && rq[1] == rq[2], next);
      nrw = moverow(rw, lm, ml, ra, rb);
      w = cmove(v, lm, al);
      ok = 1;
      for(n = 1, NCHK-1,
        if((n+1)^2*w[n+2] != subst(nrw[1],nv,n)*w[n+1] - subst(nrw[2],nv,n)*w[n], ok=0; break));
      print("  -- move lam=", lm, ", al=", al, " : ",
            if(ok, Str("closed formula VERIFIED against binomial transform to n=", NCHK), "MISMATCH"));
      print("     P# = ", nrw[1]);
      print("     Q# = ", nrw[2]);
      print("     char roots# = ", qroots2(polcoef(nrw[1],2,nv), polcoef(nrw[2],2,nv)));
    );
  );
}
for(i = 1, #corpus, report(corpus[i], i));
quit;
