default(parisize, 3000000000);

MQ = 42; NROW = 34; NK = 24;

DL = vector(NK+2, i, if(i==1, 1, lcm(vector(i-1,j,j))));   \\ DL[n+1] = lcm(1..n)

\\ companion from a fitted recurrence rec=[o,d,v]; census convention
\\ b_0=..=b_{o-2}=0, b_{o-1}=1, extended by the recurrence for n>=o
compk(rec, NN) = {
  my(o=rec[1], d=rec[2], v=rec[3], b=vector(NN+1));
  for(i=1,o, b[i]=0); b[o]=1;
  my(cf=vector(o+1, j, sum(e=0,d, v[(j-1)*(d+1)+e+1]*'n^e)));
  for(n=o, NN,
    my(s=0);
    for(j=1,o, s += subst(cf[j+1],'n,n)*b[n+1-j]);
    my(lead=subst(cf[1],'n,n));
    if(lead==0, return(-1));
    b[n+1] = -s/lead);
  my(k=0);
  while(k<=8,
    my(ok=1);
    for(n=0,NN, if(denominator(DL[n+1]^k*b[n+1])!=1, ok=0; break));
    if(ok, break); k++);
  [k, b];
};

scanlevel(N, divs, TS, FS, tlo, thi, outfile) =
{
  my(E = buildE(divs, MQ), res = List());
  my(nF = #FS, Fq = vector(nF, i, etaq(E, FS[i], 0, MQ)));
  for(ti = tlo, min(thi,#TS),
    my(tq = etaq(E, TS[ti], 1, MQ));
    if(polcoeff(tq, 1) != 1, next);
    my(qt = serreverse(tq));
    for(fi = 1, nF,
      my(Ft, g, c, lam, rec, cp, rts, l1, l2, sc, kk);
      Ft = subst(Fq[fi], 'q, qt);
      if(polcoeff(Ft, 0) != 1, next);
      g  = sqrt(Ft);
      c  = vector(NROW+1, i, polcoeff(g, i-1));
      my(m = 0, bad = 0);
      for(n = 1, NROW,
        if(c[n+1],
          my(dn = denominator(c[n+1]), v2 = valuation(dn,2));
          if(dn != 2^v2, bad = 1; break);
          if(v2 > 0, m = max(m, ceil(v2/n)))));
      if(bad, next);
      lam = 2^m;
      my(a = vector(NROW+1, i, lam^(i-1)*c[i]));
      rec = fitminp(a, MAXO, MAXD);
      if(rec == 0, next);
      cp = charpol_from(rec[3], rec[1], rec[2]);
      rts = polroots(cp);
      rts = vecsort(vector(#rts, i, abs(rts[i])), , 4);
      l1 = rts[1]; l2 = if(#rts >= 2, rts[2], 0);
      \\ k is only needed where a positive score is conceivable (l2 < 1)
      my(kv = -1);
      if(l2 < 1.5 && l2 > 1e-12,
        kk = compk(rec, NK);
        if(type(kk)=="t_VEC", kv = kk[1]));
      sc = if(l2 > 1e-12, log(1/l2) - if(kv>=0, kv, 2), -99);
      listput(res, [N, ti, fi, lam, rec[1], rec[2], kv, l1, l2, sc]);
    ));
  for(i = 1, #res, my(r = res[i]);
    write(outfile, r[1]," ",r[2]," ",r[3]," ",r[4]," ",r[5]," ",r[6]," ",r[7]," ",r[8]," ",r[9]," ",r[10]));
  #res;
};
