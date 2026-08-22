default(parisize, 4000000000);


/* exact scan over Q for one level.  Input file: level_<N>.txt with
   line 1: divisors
   line 2: number of t, then each t as exponent vector
   then   : number of F, then each F as exponent vector            */

MQ = 56; NROW = 48;

scanlevel(N, divs, TS, FS, outfile) =
{
  my(E = buildE(divs, MQ), res = List());
  my(nF = #FS, Fq = vector(nF, i, etaq(E, FS[i], 0, MQ)));
  for(ti = 1, #TS,
    my(tq = etaq(E, TS[ti], 1, MQ));
    if(polcoeff(tq, 1) != 1, next);
    my(qt = serreverse(tq));
    for(fi = 1, nF,
      my(Ft, g, c, lam, ee, rec, cp, rts, l1, l2, sc);
      Ft = subst(Fq[fi], 'q, qt);
      if(polcoeff(Ft, 0) != 1, next);
      g  = sqrt(Ft);
      c  = vector(NROW+1, i, polcoeff(g, i-1));
      \\ minimal lambda = 2^m
      my(m = 0);
      for(n = 1, NROW, if(c[n+1], my(v = valuation(c[n+1], 2)); if(v < 0, m = max(m, ceil(-v/n)))));
      lam = 2^m;
      \\ any odd prime in the denominators?  (should not happen)
      my(bad = 0);
      for(n = 1, NROW, if(c[n+1] && denominator(c[n+1]) != 2^valuation(denominator(c[n+1]),2), bad = 1; break));
      if(bad, next);
      my(a = vector(NROW+1, i, lam^(i-1)*c[i]));
      rec = fitmin(a, MAXO, MAXD);
      if(rec == 0, next);
      cp = charpol_from(rec[3], rec[1], rec[2]);
      rts = polroots(cp);
      rts = vecsort(vector(#rts, i, abs(rts[i])), , 4);   \\ decreasing
      l1 = rts[1]; l2 = if(#rts >= 2, rts[2], 0);
      sc = if(l2 > 1e-12, log(1/l2) - 2, -99);
      listput(res, [N, ti, fi, lam, rec[1], rec[2], l1, l2, sc]);
    ));
  write(outfile, "");
  for(i = 1, #res, my(r = res[i]);
    write(outfile, r[1], " ", r[2], " ", r[3], " ", r[4], " ", r[5], " ", r[6], " ",
          r[7], " ", r[8], " ", r[9]));
  #res;
};
