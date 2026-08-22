/* 04_periods.gp -- archimedean Apery limits of every orbit member, to enough
   precision to detect rational relations inside an orbit, plus identification
   against the standard constants.                                            */
read("lib.gp"); read("rows.gp");
default(parisize, 2^32);
default(realprecision, 400);

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
    nf = List();
    for(i=1,#front, ms = allmoves(front[i]);
      for(j=1,#ms, kk=rowkey(ms[j]); seen=0;
        for(l=1,#keys, if(keys[l]==kk, seen=1; break));
        if(!seen, listput(keys,kk); listput(reps,ms[j]); listput(nf,ms[j]))));
    front = nf);
  Vec(reps); }

/* archimedean Apery limit to ~DIG digits */
DIG = 60;
{ xival(rw) =
  my(rt=qroots2(polcoef(rw[1],2,nv),polcoef(rw[2],2,nv)), l1, l2, rr, nn, v, w);
  if(#rt==0, return([0,0]));
  l1 = rt[1]; l2 = rt[2];
  if(abs(l1)<abs(l2), rr=l1; l1=l2; l2=rr);
  if(l2==0 || abs(l2)>=abs(l1), return([0,0]));
  rr = log(abs(l1*1.0/l2))/log(10);
  nn = min(6000, max(400, ceil((DIG+20)/rr)));
  v = seqA(rw, nn); w = seqB(rw, nn);
  if(v[nn+1]==0, return([0,0]));
  [w[nn+1]/v[nn+1]*1.0, nn]; }

zt2   = zeta(2);
lchi3 = lfun(-3, 2);
catG  = lfun(-4, 2);
zt3   = zeta(3);
print("constants: zeta(2)=", strprintf("%.20f",zt2), " L(2,chi-3)=", strprintf("%.20f",lchi3), " G=", strprintf("%.20f",catG));

{ for(i=1,#corpus,
  my(rw=corpus[i], rt, orb, xs=List(), nsv=List(), lv);
  rt = qroots2(polcoef(rw[1],2,nv), polcoef(rw[2],2,nv));
  if(#rt==0, next);
  print(""); print("################ ", rw[3], "   ", rw[4]);
  orb = orbit([rw[1],rw[2]]);
  for(j=1,#orb,
    my(z = xival(orb[j]));
    listput(xs, z[1]); listput(nsv, z[2]);
    print("  ", rw[3], ".", j, "  (a,b,d)=(", polcoef(orb[j][1],2,nv), ",",
          subst(orb[j][1],nv,0), ",", polcoef(orb[j][2],2,nv), ")",
          "  P=", orb[j][1], "  Q=", orb[j][2]);
    print("        n=", z[2], "  xi = ", strprintf("%.45f", z[1])));
  /* rational ratios inside the orbit */
  print("  -- ratios xi_j/xi_1 : best rational of height <= 10^12, and its error");
  for(j=1,#orb,
    if(xs[1]!=0 && xs[j]!=0,
      my(rt2 = xs[j]/xs[1], q = bestappr(rt2, 10^12));
      print("     j=",j,"  xi_j/xi_1 = ", strprintf("%.30f", rt2), "  ~ ", q,
            "   err=", strprintf("%.3e", abs(rt2-q*1.0)))));
  print("  -- lindep [xi_j, zeta(2), L(2,chi-3), G, zeta(3)]  and  [xi_j, xi_1, 1]:");
  for(j=1,#orb,
    if(xs[j]!=0,
      lv = lindep([xs[j], zt2, lchi3, catG, zt3], 35);
      print("     j=",j,"  const: ", lv, "    vs xi_1: ",
            if(j>1, lindep([xs[j], xs[1], 1], 35), "--"))));
); }
quit;
