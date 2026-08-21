\p 2500
G = Catalan;
/* honest lattice: q=(cZ X+cN V)/M and p=(cZ Y+cN U)/M both integers */
{
klat2(X,Y,V,U,M) = my(K=matkerint([X,V,M,0,0,0; Y,U,0,M,0,0]));
  mathnf(matrix(2,#K[1,],i,j,K[i,j]));
}
{
redu(B,lz,ln) = my(sc=10^(300+ceil(log(1/min(lz,ln))/log(10))));
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
run(nlist,k,R) = my();
 printf("HONEST lattice (q,p both integral).  k=%.4f R=%d   [k*=22.35129, F(k)=%.5f]\n",k,R,(13.0995887908+14.3931452672-12-k*log(2))/2);
 printf("  n | logM/n  logidx/n | logL1/n | logFcone/n  logqcone/n | logFmin/n (any sign)\n");
 for(t=1,#nlist,
  my(n=nlist[t], D=lcm(vector(6*n,i,i)), S=D^2, T=2^floor(k*n), M=S*T,
     zr=zudrow(n), nr=nestrow(n), X=zr[1],Y=zr[2],V=nr[1],U=nr[2]);
  my(B0=klat2(X,Y,V,U,M), idx=abs(matdet(B0)));
  my(LZ=(X*G-Y)/M, LN=(V*G-U)/M, sZ=sign(LZ), sN=sign(LN), lz=abs(LZ), ln=abs(LN));
  my(B=[sZ*B0[1,1],sZ*B0[1,2]; sN*B0[2,1],sN*B0[2,2]], Br=redu(B,lz,ln));
  my(l1=abs(Br[1,1])*lz+abs(Br[2,1])*ln);
  my(mn=+oo); for(i=-40,40,for(j=-40,40,if(i||j, my(v=abs((i*Br[1,1]+j*Br[1,2])*lz-(i*Br[2,1]+j*Br[2,2])*ln*(-1))); )));
  my(fm=+oo); for(i=-40,40,for(j=-40,40,if(i||j,
      my(cz=i*B0[1,1]+j*B0[1,2], cn=i*B0[2,1]+j*B0[2,2], v=abs(cz*LZ+cn*LN)); if(v<fm,fm=v))));
  my(r=conemin(Br,lz,ln,R), fc=r[1], cv=r[2]);
  my(cz=sZ*cv[1], cn=sN*cv[2], q=(cz*X+cn*V)/M, p=(cz*Y+cn*U)/M);
  my(bad = (type(q)!="t_INT") || (type(p)!="t_INT") || abs(abs(q*G-p)-fc)>fc*1e-200);
  printf("%4d | %7.4f %8.4f | %8.4f | %10.4f %11.4f | %10.4f%s\n",
   n, log(M)/n, log(idx)/n, log(l1)/n, log(fc)/n, log(abs(q))/n, log(fm)/n, if(bad," ***BAD***",""));
 );
}
run([4,6,8,10,12,14,16,18,20,24,28,32],22.4,400);
\q
