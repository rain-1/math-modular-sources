/* 06_padic.gp -- p-adic Apery limits across a cusp-move orbit.
   For each orbit member and each prime p dividing d we measure the Cauchy slope
   v_p(xi_N - xi_{N-1}) and then test, for every pair (i,j) of members with a
   common prime, whether  q1*xi^i = q2*xi^j  for a small rational q1/q2.      */
read("lib.gp"); read("rows.gp");
default(parisize, 2^31);
NP = 400;

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
{ rowkey(rw) = my(aa=polcoef(rw[1],2,nv), dd=polcoef(rw[2],2,nv));
  if(aa==0, return([0, polcoef(rw[1],1,nv)^2/dd, polcoef(rw[1],0,nv)^2/dd,
                    polcoef(rw[2],1,nv)/dd, polcoef(rw[2],0,nv)/dd]));
  [aa^2/dd, polcoef(rw[1],1,nv)/aa, polcoef(rw[1],0,nv)/aa,
   polcoef(rw[2],1,nv)/dd, polcoef(rw[2],0,nv)/dd]; }
{ allmoves(rw) =
  my(rt=qroots2(polcoef(rw[1],2,nv),polcoef(rw[2],2,nv)), rq=rrootsQ(rw),
     out=List(), nr=0, jm=2);
  if(#rt==0 || #rq==0, return(Vec(out)));
  if(rq[1]==rq[2], jm=1);
  for(i=1,2, for(j=1,jm,
      nr = moverow(rw, rt[i], rt[3-i], rq[j], rq[3-j]);
      if(polcoef(nr[2],2,nv)!=0, listput(out,nr))));
  Vec(out); }
{ orbit(rw, maxsz=40) =
  my(keys=List(), reps=List(), front=List([rw]), nf, ms, kk, seen);
  listput(keys,rowkey(rw)); listput(reps,rw);
  while(#front>0 && #reps<maxsz,
    nf=List();
    for(i=1,#front, ms=allmoves(front[i]);
      for(j=1,#ms, kk=rowkey(ms[j]); seen=0;
        for(l=1,#keys, if(keys[l]==kk, seen=1; break));
        if(!seen, listput(keys,kk); listput(reps,ms[j]); listput(nf,ms[j]))));
    front=nf); Vec(reps); }
{ intscale(rw) = my(v=seqA(rw,60)); minscale(v, 10^6); }
{ pcompare(a, b, p, xa, xb) =
  my(va = if(xa==0, 10^9, valuation(xa,p)), vb = if(xb==0, 10^9, valuation(xb,p)), rr);
  if(va > 100, print("     ", a, " vs ", b, " at p=", p, " : xi_a = 0 (v_p >= ", va, ")"); return(0));
  if(vb > 100, print("     ", a, " vs ", b, " at p=", p, " : xi_b = 0 (v_p >= ", vb, ")"); return(0));
  rr = bestappr(xa/xb + O(p^40), 10^6);
  if(type(rr) != "t_INT" && type(rr) != "t_FRAC",
     print("     ", a, " vs ", b, " at p=", p, " : no rational of height <= 10^6"); return(0));
  print("     ", a, " vs ", b, " at p=", p, " : xi_a/xi_b = ", rr,
        "   [v_p(xi_a - r*xi_b) = ", if(xa-rr*xb==0, "inf", valuation(xa-rr*xb,p)), "]"); }

{ for(i=1,#corpus,
  my(rw=corpus[i], rt, orb, xis=List(), pri=List(), nm=rw[3]);
  rt = qroots2(polcoef(rw[1],2,nv), polcoef(rw[2],2,nv));
  if(#rt==0, next);
  print(""); print("################ ", nm, "  ", rw[4]);
  orb = orbit([rw[1],rw[2]]);
  for(j=1,#orb,
    my(c=intscale(orb[j]), rw2, av, bv, dd, ff, lst=List(), lstp=List());
    if(c==0, print("  ",nm,".",j," no integral scaling"); listput(xis,0); listput(pri,[]); next);
    rw2 = [c*orb[j][1], c^2*orb[j][2]];
    av = seqA(rw2, NP); bv = seqB(rw2, NP);
    dd = polcoef(rw2[2],2,nv);
    ff = factor(abs(dd));
    print("  ", nm, ".", j, "  (a,b,d)=(", polcoef(rw2[1],2,nv), ",",
          subst(rw2[1],nv,0), ",", dd, ")  scale=", c);
    for(m=1,#ff~,
      my(p=ff[m,1], x1, x2, s1, s2);
      if(av[NP+1]==0 || av[NP]==0 || av[NP\2+1]==0 || av[NP\2]==0, next);
      x1 = bv[NP+1]/av[NP+1]; x2 = bv[NP]/av[NP];
      s1 = if(x1-x2==0, oo, valuation(x1-x2,p));
      x1 = bv[NP\2+1]/av[NP\2+1]; x2 = bv[NP\2]/av[NP\2];
      s2 = if(x1-x2==0, oo, valuation(x1-x2,p));
      print("        p=", p, "  Cauchy slope v_p(xi_N - xi_{N-1}): N=", NP\2, " -> ", s2,
            ",  N=", NP, " -> ", s1, "   [v_p(d)=", valuation(dd,p), "]");
      if(s1 > 20,
         print("           v_p(xi) = ", if(bv[NP+1]/av[NP+1]==0, "inf",
               valuation(bv[NP+1]/av[NP+1], p)), "   (0 to precision ", s1, " means xi_p = 0)");
         listput(lst, bv[NP+1]/av[NP+1]); listput(lstp, p)));
    listput(xis, [lst, lstp]));
  /* pairwise p-adic comparison */
  print("  -- pairwise p-adic ratios xi_p^(i)/xi_p^(j) where both converge:");
  for(a=1,#orb, for(b=a+1,#orb,
    if(type(xis[a])=="t_VEC" && type(xis[b])=="t_VEC" && #xis[a]==2 && #xis[b]==2,
      for(u=1,#xis[a][2], for(w=1,#xis[b][2],
        if(xis[a][2][u]==xis[b][2][w], pcompare(a,b,xis[a][2][u],xis[a][1][u],xis[b][1][w])))))));
); }
quit;
