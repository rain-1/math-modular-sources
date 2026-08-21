\p 3000
{
klat2(X,Y,V,U,M) = my(K=matkerint([X,V,M,0,0,0; Y,U,0,M,0,0]));
  mathnf(matrix(2,#K[1,],i,j,K[i,j]));
}
{
redu(B,lz,ln) = my(sc=10^(400+ceil(log(1/min(lz,ln))/log(10))));
  B*qflll([round(B[1,1]*lz*sc),round(B[1,2]*lz*sc); round(B[2,1]*ln*sc),round(B[2,2]*ln*sc)]);
}
{
conemin(B,lz,ln,R) = my(best=0,bv=0);
 for(i=-R,R, my(lo=-oo,hi=oo,ok=1);
  for(r=1,2, my(a=i*B[r,1],b=B[r,2]);
    if(b==0, if(a<0,ok=0), my(t=-a/b); if(b>0,if(t>lo,lo=t),if(t<hi,hi=t))));
  if(ok, my(j0=if(lo==-oo,-R,ceil(lo)), j1=if(hi==oo,R,floor(hi)));
   if(j0<=j1, for(w=0,1, my(j=if(w==0,j0,j1), cz=i*B[1,1]+j*B[1,2], cn=i*B[2,1]+j*B[2,2]);
     if(cz||cn, my(v=cz*lz+cn*ln); if(best==0||v<best,best=v;bv=[cz,cn]))))));
 [best,bv];
}
{
run(nlist,k,R,GG,tag) = my();
 printf("\n=== %s   k=%.4f  R=%d ===\n",tag,k,R);
 printf("  n | logidx/n | log(cov^.5)/n | logL1/n logL2/n | logFcone/n | Fcone/L1 | logq/n\n");
 for(t=1,#nlist,
  my(n=nlist[t], D=lcm(vector(6*n,i,i)), S=D^2, T=2^floor(k*n), M=S*T,
     zr=zudrow(n), nr=nestrow(n), X=zr[1],Y=zr[2],V=nr[1],U=nr[2]);
  my(B0=klat2(X,Y,V,U,M), idx=abs(matdet(B0)));
  my(LZ=(X*GG-Y)/M, LN=(V*GG-U)/M, sZ=sign(LZ), sN=sign(LN), lz=abs(LZ), ln=abs(LN));
  my(B=[sZ*B0[1,1],sZ*B0[1,2]; sN*B0[2,1],sN*B0[2,2]], Br=redu(B,lz,ln));
  my(nrm(v)=sqrt((v[1]*lz)^2+(v[2]*ln)^2));
  my(l1=nrm([Br[1,1],Br[2,1]]), l2=nrm([Br[1,2],Br[2,2]]));
  if(l1>l2, my(tmp=l1); l1=l2; l2=tmp);
  my(cov=idx*lz*ln);
  my(r=conemin(Br,lz,ln,R), fc=r[1], cv=r[2]);
  my(cz=sZ*cv[1], cn=sN*cv[2], q=(cz*X+cn*V)/M, p=(cz*Y+cn*U)/M);
  my(bad=(type(q)!="t_INT")||(type(p)!="t_INT")||abs(abs(q*GG-p)-fc)>fc*1e-300);
  printf("%4d | %8.4f | %13.4f | %7.4f %7.4f | %10.4f | %8.3f | %8.4f%s\n",
   n, log(idx)/n, log(cov)/(2*n), log(l1)/n, log(l2)/n, log(fc)/n, fc/l1, log(abs(q))/n,
   if(bad," ***BAD***",""));
 );
}
NL=[6,10,14,18,22,26,30,34,38];
run(NL,22.4,600,Catalan,"true G");
Gs=bestappr(Catalan,10^320);
printf("\ncontrol: G*=a/b, log b=%.2f, |G-G*|=e^%.1f\n",log(denominator(Gs)),log(abs(Catalan-Gs)));
run(NL,22.4,600,Gs,"rational G* control");
\q
