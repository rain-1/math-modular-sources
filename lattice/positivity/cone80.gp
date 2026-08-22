/* lattice/positivity/cone80.gp   --- Task 2 of the positivity programme.
   Exact congruence lattice K_n of the Zudilin(3n) x Nesterenko(4,7)(n) pair,
   its successive minima in the scaled metric, and the minimum of the linear
   form over K_n intersect the POSITIVE CONE (where non-vanishing is a theorem).
   Prepend lattice/positivity/rows_pos.gp.
   Usage:  gp -q -D parisize=1500000000 run.gp   with  KK, NLO, NHI, NSTEP, NRES set.
   Emits CSV:  n,kappa,halflogcov,logl1,logl2,logcone,ratio,logq,logp,cz,cn  */

\p 4000

/* exact cone minimum of  c1*lz + c2*ln  over  {c in L, c>=0} \ {0},
   L given by integer basis matrix BB (columns = basis vectors).
   Sweep i over [-R,R]; for each i the admissible j form an interval
   (exact integer arithmetic), and the linear objective is minimised at an
   endpoint.  Returns [value, c] or [0,0] if the cone is empty in range. */
{
conemin(BB, lz, ln, R) = my(best=0, bv=0);
 for(i=-R, R,
  my(lo=-oo, hi=oo, ok=1);
  for(r=1,2,
   my(aa=i*BB[r,1], bb=BB[r,2]);
   if(bb==0,
     if(aa<0, ok=0),
     my(tt=-aa/bb);
     if(bb>0, if(tt>lo, lo=tt), if(tt<hi, hi=tt))));
  if(ok,
   my(j0=if(lo==-oo, -R, ceil(lo)), j1=if(hi==oo, R, floor(hi)));
   if(j0<=j1,
    for(w=0,1,
     my(j=if(w==0,j0,j1),
        c1=i*BB[1,1]+j*BB[1,2], c2=i*BB[2,1]+j*BB[2,2]);
     if(c1||c2,
      my(v=c1*lz+c2*ln);
      if(best==0 || v<best, best=v; bv=[c1,c2]))))));
 [best,bv];
}

/* LLL reduction of BB in the metric diag(lz,ln) (scaled to integers) */
{
redu(BB, lz, ln) =
 my(EE = ceil(log(1.*vecmax([abs(BB[1,1]),abs(BB[1,2]),abs(BB[2,1]),abs(BB[2,2])]))/log(10)),
    DD = min(precision(1.*lz), precision(1.*ln)),
    off = max(60, min(600, DD - 40 - EE)),
    sc  = 10^(off + ceil(log(1/min(lz,ln))/log(10))));
 BB*qflll([round(BB[1,1]*lz*sc), round(BB[1,2]*lz*sc);
           round(BB[2,1]*ln*sc), round(BB[2,2]*ln*sc)]);
}

{
/* one n, several k: build the exact rows once, then loop over k. */
sweepn(n, klist, GG) =
 localprec(4000);
 my(DD=lcm(vector(6*n,i,i)), SS=DD^2,
    zr=zudrow(n), nr=nestrow(n), XX=zr[1], YY=zr[2], VV=nr[1], UU=nr[2]);
 for(ii=1,#klist,
  my(kk=klist[ii], TT=2^floor(kk*n), MOD=SS*TT);
  my(B0=klat2(XX,YY,VV,UU,MOD), idx=abs(matdet(B0)));
  my(LZ=(XX*GG-YY)/MOD, LN=(VV*GG-UU)/MOD,
     sZ=sign(LZ), sN=sign(LN), lz=abs(LZ), ln=abs(LN));
  my(BB=[sZ*B0[1,1], sZ*B0[1,2]; sN*B0[2,1], sN*B0[2,2]]);
  my(Br=redu(BB, lz, ln));
  my(nrm(v) = sqrt((v[1]*lz)^2 + (v[2]*ln)^2));
  my(v1=[Br[1,1],Br[2,1]], v2=[Br[1,2],Br[2,2]]);
  my(l1=nrm(v1), l2=nrm(v2), w1=v1);
  if(l1>l2, my(tmp=l1); l1=l2; l2=tmp; w1=v2);
  /* does the (unrestricted) shortest vector lie in the closed cone, up to +-? */
  my(incone = if((w1[1]>=0 && w1[2]>=0) || (w1[1]<=0 && w1[2]<=0), 1, 0));
  my(cov=idx*lz*ln);
  my(R=200, r1=conemin(Br,lz,ln,R), r2=conemin(Br,lz,ln,2*R), guard=0);
  while(guard<5 && (r1[1]==0 || r1[1]!=r2[1]), R*=2; r1=r2; r2=conemin(Br,lz,ln,2*R); guard++);
  my(fc=r2[1], cv=r2[2]);
  my(c1=cv[1], c2=cv[2], cz=sZ*c1, cn=sN*c2,
     q=(cz*XX+cn*VV)/MOD, p=(cz*YY+cn*UU)/MOD);
  my(bad = (type(q)!="t_INT") || (type(p)!="t_INT")
         || abs(abs(q*GG-p)-fc) > fc*10.^(-precision(1.*lz)+50));
  printf("%.4f,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%d\n",
    kk, n, log(idx)/n, log(cov)/(2*n), log(l1)/n, log(l2)/n, log(fc)/n,
    fc/l1, log(abs(q))/n, incone, R, guard, if(bad,1,0));
 );
}
