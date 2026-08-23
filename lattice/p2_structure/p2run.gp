/* lattice/p2_structure/p2run.gp   --- per-n structural analysis.
   Prepend: lattice/positivity/rows_pos.gp, lattice/p2_structure/p2core.gp.
   Driver sets GG (Catalan at >=3000 digits), RW = rdrows(...), KLIST, and
   calls anal(n, KLIST, GG, RW, DUMP).   Emits one CSV line per (k,n):

 k,n,logidx,halflogcov,logl1,logl2,skew,incone,logcone,rho,rho2,logq,
 jhi,mrg,ok,cflen,cfmax,cfidx,cfkind,cfpq,logh11,h22,ratio
*/

/* best-approximation classification of the pair (jj,aa) for h12/h11 */
{
cfclass(h11, h12, jj) =
 my(cf=contfrac(h12/h11), qs=vector(#cf), kind=-1, idx=-1, tt=0);
 for(i=1,#cf, qs[i]=if(i==1, 1, if(i==2, cf[2], cf[i]*qs[i-1]+qs[i-2])));
 for(i=1,#cf, if(qs[i]==jj, kind=0; idx=i-1; break));
 if(kind<0,
   for(i=1,#cf-2,
     for(t=1,cf[i+2]-1,
       if(qs[i]+t*qs[i+1]==jj, kind=1; idx=i-1; tt=t; break));
     if(kind>=0, break)));
 [kind, idx, tt, #cf, vecmax(cf), if(idx>=0 && idx+2<=#cf, cf[idx+2], 0)];
}

{
anal(n, klist, GG, RW, DUMP) =
 my(rw=mapget(RW,n), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*n,i,i)), SS=DD^2);
 for(ii=1,#klist,
  my(kk=klist[ii], TT=2^floor(kk*n), MOD=SS*TT);
  my(B0=kfull(XX,YY,VV,UU,MOD), idx=abs(matdet(B0)));
  my(LZ=(XX*GG-YY)/MOD, LN=(VV*GG-UU)/MOD,
     sZ=sign(LZ), sN=sign(LN), lz=abs(LZ), ln=abs(LN));
  my(BB=[sZ*B0[1,1], sZ*B0[1,2]; sN*B0[2,1], sN*B0[2,2]]);
  my(gr=gred(BB,lz,ln), Br=gr[1], l1=gr[2], l2=gr[3], mrg=gr[4]);
  my(cov=idx*lz*ln, sinth=cov/(l1*l2));
  my(v1=[Br[1,1],Br[2,1]], v2=[Br[1,2],Br[2,2]]);
  my(incone = if((v1[1]>=0 && v1[2]>=0) || (v1[1]<=0 && v1[2]<=0), 1, 0));
  my(L0=eight(Br,lz,ln));
  my(jhi=if(L0==0, 8, floor(L0*l1/cov)+2));
  my(sc=conescan(Br,lz,ln,jhi), sc2=conescan(Br,lz,ln,jhi+3));
  my(fc=sc[1], cv=sc[2], l1q=sc[3], rho=fc/l1, rho2=l1q/l1);
  my(stable = (sc[1]==sc2[1]) && (sc[3]==sc2[3]));
  my(cz=sZ*cv[1], cn=sN*cv[2],
     q=(cz*XX+cn*VV)/MOD, p=(cz*YY+cn*UU)/MOD);
  my(ok = (type(q)=="t_INT") && (type(p)=="t_INT") && stable
        && abs(abs(q*GG-p)-fc) < fc*10.^(-precision(1.*lz)+80));
  /* Hermite form of the oriented lattice, and the CF classification of v1 */
  my(HH=mathnf(BB), h11=HH[1,1], h12=HH[1,2], h22=HH[2,2]);
  my(jj=abs(v1[2])/h22, cc=if(jj==0, [-2,-1,0,0,0,0], cfclass(h11,h12,jj)));
  my(jc=abs(cv[2])/h22, dd=if(jc==0, [-2,-1,0,0,0,0], cfclass(h11,h12,jc)));
  printf("%.4f,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%.6f,%d,%.3e,%d,%d,%d,%d,%d,%d,%.6f,%d,%.6f,%d,%d,%d,%d\n",
    kk, n, log(idx)/n, log(cov)/(2*n), log(l1)/n, log(l2)/n, l2/l1, incone,
    log(fc)/n, rho, rho2, log(abs(q))/n, jhi, mrg, if(ok,1,0),
    cc[4], cc[5], cc[2], cc[1], cc[3], log(h11)/n, h22, log(sinth),
    cc[6], dd[1], dd[2], dd[3]);
  if(DUMP!="",
    write(DUMP, kk, " ", n, " ", v1[1], " ", v1[2], " ", v2[1], " ", v2[2],
          " ", cv[1], " ", cv[2], " ", sZ, " ", sN)));
}

/* same for the single-congruence lattice of 06_threshold.tex: c.(Y,U)=0 mod S */
{
analS(n, GG, RW) =
 my(rw=mapget(RW,n), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*n,i,i)), SS=DD^2);
 my(B0=kcong(YY,UU,SS), idx=abs(matdet(B0)));
 my(LZ=(XX*GG-YY)/SS, LN=(VV*GG-UU)/SS,
    sZ=sign(LZ), sN=sign(LN), lz=abs(LZ), ln=abs(LN));
 my(BB=[sZ*B0[1,1], sZ*B0[1,2]; sN*B0[2,1], sN*B0[2,2]]);
 my(gr=gred(BB,lz,ln), Br=gr[1], l1=gr[2], l2=gr[3]);
 my(cov=idx*lz*ln);
 my(v1=[Br[1,1],Br[2,1]]);
 my(incone = if((v1[1]>=0 && v1[2]>=0) || (v1[1]<=0 && v1[2]<=0), 1, 0));
 my(L0=eight(Br,lz,ln), jhi=if(L0==0, 8, floor(L0*l1/cov)+2));
 my(sc=conescan(Br,lz,ln,jhi), fc=sc[1], l1q=sc[3]);
 my(HH=mathnf(BB), h11=HH[1,1], h12=HH[1,2], h22=HH[2,2]);
 my(jj=abs(v1[2])/h22, cc=if(jj==0, [-2,-1,0,0,0], cfclass(h11,h12,jj)));
 printf("S,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,%d\n",
   n, log(idx)/n, log(cov)/(2*n), log(l1)/n, log(l2)/n, l2/l1, incone,
   log(fc)/n, fc/l1, l1q/l1, gcd(gcd(YY,UU),SS)>1, cc[4], cc[5], cc[2], cc[1], cc[3]);
}

/* Gauss-Kuzmin histogram of the partial quotients of h12/h11 */
{
cfhist(n, klist, GG, RW) =
 my(rw=mapget(RW,n), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*n,i,i)), SS=DD^2);
 for(ii=1,#klist,
  my(kk=klist[ii], TT=2^floor(kk*n), MOD=SS*TT);
  my(B0=kfull(XX,YY,VV,UU,MOD));
  my(LZ=(XX*GG-YY)/MOD, LN=(VV*GG-UU)/MOD, sZ=sign(LZ), sN=sign(LN));
  my(BB=[sZ*B0[1,1], sZ*B0[1,2]; sN*B0[2,1], sN*B0[2,2]], HH=mathnf(BB));
  my(cf=contfrac(HH[1,2]/HH[1,1]), h=vector(9));
  for(i=2,#cf, my(a=cf[i]); if(a<=8, h[a]++, h[9]++));
  printf("%.4f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", kk, n, #cf-1,
    h[1],h[2],h[3],h[4],h[5],h[6],h[7],h[8],h[9]));
}
