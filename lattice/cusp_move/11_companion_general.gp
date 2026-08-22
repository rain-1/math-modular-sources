/* 11_companion_general.gp -- verification of the general companion transform
      B# = (1-lam t)^al * Bhat ,   L Bhat = t (1-lam t)^{-al},
   i.e.  bhat_n from  (n+1)^2 bhat_{n+1} = P(n) bhat_n - Q(n) bhat_{n-1} + g_{n+1},
         g_k = lam^{k-1} binomial(al+k-2, k-1),
   and   b#_n = sum_m bhat_m (-lam)^{n-m} binomial(n+al-1, n-m).
   Tested on all corpus rows with rational characteristic roots, both roots,
   both gauges.                                                              */
read("lib.gp"); read("rows.gp");
default(parisize, 2^30);
NT = 40;

{ qroots2(aa,bb) = my(ds=aa^2-4*bb, sq); if(issquare(ds,&sq), [(aa+sq)/2,(aa-sq)/2], []); }
{ rrootsQ(rw) = my(q2=polcoef(rw[2],2,nv), q1=polcoef(rw[2],1,nv), q0=polcoef(rw[2],0,nv));
  if(q2==0, return([])); qroots2(-q1/q2, q0/q2); }
{ moverow(rw,lam,mu,r1,r2) =
  my(p1=polcoef(rw[1],1,nv), p0=polcoef(rw[1],0,nv), pp, qq);
  pp = (mu-2*lam)*nv^2 + (p1-3*lam+2*lam*r1)*nv + (p0-lam+lam*r1);
  qq = (lam^2-lam*mu)*nv^2
     + (-2*lam^2*r1+lam^2+lam*mu*r1-lam*mu*r2+lam*mu-lam*p1)*nv
     + (lam^2*r1^2-lam^2*r1+lam*mu*r1*r2-lam*mu*r1+lam*p1*r1);
  [pp,qq]; }

{ for(i=1,#corpus,
  my(rw=corpus[i], rt, rq, jm, lm, ml, ra, rb, al, nrw, bh, bt, btrue, av, aw, okA, okB, bn, okN);
  rt = qroots2(polcoef(rw[1],2,nv), polcoef(rw[2],2,nv));
  rq = rrootsQ(rw);
  if(#rt==0 || #rq==0, next);
  jm = if(rq[1]==rq[2], 1, 2);
  av = seqA([rw[1],rw[2]], NT);
  for(s=1,2, for(j=1,jm,
    lm = rt[s]; ml = rt[3-s]; ra = rq[j]; rb = rq[3-j]; al = 1-ra;
    nrw = moverow([rw[1],rw[2]], lm, ml, ra, rb);
    if(polcoef(nrw[2],2,nv)==0, next);
    aw = cmove(av, lm, al);
    okA = 1; for(n=1,NT-1,
      if((n+1)^2*aw[n+2] != subst(nrw[1],nv,n)*aw[n+1] - subst(nrw[2],nv,n)*aw[n], okA=0; break));
    bh = seqBg([rw[1],rw[2]], NT, (k) -> lm^(k-1)*binomial(al+k-2, k-1));
    bt = cmove(bh, lm, al);
    btrue = seqB(nrw, NT);
    okB = 1; for(n=0,NT, if(bt[n+1] != btrue[n+1], okB=0; break));
    bn = cmove(seqB([rw[1],rw[2]], NT), lm, al);
    okN = 1; for(n=0,NT, if(bn[n+1] != btrue[n+1], okN=0; break));
    print(rw[3], "  lam=", lm, " al=", al, " :  A-transform ", if(okA,"OK","FAIL"),
          " ;  B# = cmove(Bhat) ", if(okB,"OK","FAIL"),
          " ;  naive B# = cmove(B) ", if(okN,"(also OK!)","fails as predicted")))));
}
quit;
