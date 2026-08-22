/* 08_quadratic.gp -- the placements of the irrational-root rows over the real
   quadratic field K = Q(sqrt(a^2-4d)):  the moved rows, the minimal rescaling
   making the sequence an algebraic integer sequence, and the corrected score.  */
read("lib.gp"); read("rows.gp");
default(parisize, 2^30);
default(realprecision, 60);
NQ = 60;

/* z = Mod(u + v*y, y^2 - A y + B) is an algebraic integer iff Tr,N in Z */
{ isalgint(z, aa, dd) =
  my(p = lift(z), u = polcoef(p,0,'y), v = polcoef(p,1,'y));
  (denominator(2*u + v*aa) == 1) && (denominator(u^2 + u*v*aa + v^2*dd) == 1); }

{ report(rw) =
  my(aa = polcoef(rw[1],2,nv), dd = polcoef(rw[2],2,nv), ds, yy, lam, mu,
     q2 = polcoef(rw[2],2,nv), q1 = polcoef(rw[2],1,nv), q0 = polcoef(rw[2],0,nv),
     rr, r1, r2, al, p1, p0, pp, qq, v, c, ok, l1, l2, sq);
  ds = aa^2-4*dd;
  if(ds <= 0 || issquare(ds), return(0));
  print(""); print("=== ", rw[3], "   ", rw[4]);
  yy  = Mod('y, 'y^2 - aa*'y + dd);
  lam = yy; mu = aa - yy;
  sq  = sqrt(ds*1.0);
  rr  = if(q2==0, [0,0], my(dsq = (q1/q2)^2 - 4*q0/q2, s2);
            if(issquare(dsq, &s2), [(-q1/q2+s2)/2, (-q1/q2-s2)/2], error("Q roots irrational")));
  p1 = polcoef(rw[1],1,nv); p0 = polcoef(rw[1],0,nv);
  /* placement at v=mu : move by lam=mu, mu_new=lam  (the one sharing |lambda_2|) */
  for(sel = 1, 2,
    my(lm, ml, num1, num2, al1);
    lm = if(sel==1, lam, mu); ml = if(sel==1, mu, lam);
    for(jj = 1, if(rr[1]==rr[2], 1, 2),
      r1 = rr[jj]; r2 = rr[3-jj]; al1 = 1-r1;
      pp = (ml-2*lm)*nv^2 + (p1-3*lm+2*lm*r1)*nv + (p0-lm+lm*r1);
      qq = (lm^2-lm*ml)*nv^2
         + (-2*lm^2*r1+lm^2+lm*ml*r1-lm*ml*r2+lm*ml-lm*p1)*nv
         + (lm^2*r1^2-lm^2*r1+lm*ml*r1*r2-lm*ml*r1+lm*p1*r1);
      print("  -- placement at v=", if(sel==1,"lam","mu"), ", gauge al=", al1);
      print("     P# = ", pp);
      print("     Q# = ", qq);
      v = seqA([pp,qq], NQ);
      c = 0;
      for(cc = 1, 64,
        ok = 1;
        for(n = 0, NQ, if(!isalgint(cc^n*v[n+1], aa, dd), ok = 0; break));
        if(ok, c = cc; break));
      print("     minimal integer rescaling c with c^n a_n in O_K (n<=", NQ, "): ",
            if(c==0, ">64", c));
      my(ppn, qqn, cc2, anew, dnew, l1n, l2n, dsn);
      cc2 = if(c==0, 1, c);
      ppn = subst(lift(pp), 'y, (aa+sq)/2);
      qqn = subst(lift(qq), 'y, (aa+sq)/2);
      anew = cc2*polcoef(ppn,2,nv); dnew = cc2^2*polcoef(qqn,2,nv);
      dsn = anew^2 - 4*dnew;
      if(dsn >= 0,
        l1n = (anew+sqrt(dsn))/2; l2n = (anew-sqrt(dsn))/2;
        if(abs(l1n) < abs(l2n), my(tm=l1n); l1n=l2n; l2n=tm);
        print("     after rescaling: lam1=", l1n, "  lam2=", l2n,
              "   score(k=2) = ", log(1/abs(l2n))-2),
        print("     after rescaling: complex roots, |lam|=", sqrt(abs(dnew))));
    ));
}
for(i=1,#corpus, report(corpus[i]));
quit;
