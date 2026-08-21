\p 3000
{ klat2(X,Y,V,U,M) = my(K=matkerint([X,V,M,0,0,0; Y,U,0,M,0,0]));
  mathnf(matrix(2,#K[1,],i,j,K[i,j])); }
{ redu(B,lz,ln) = my(sc=10^(400+ceil(log(1/min(lz,ln))/log(10))));
  B*qflll([round(B[1,1]*lz*sc),round(B[1,2]*lz*sc); round(B[2,1]*ln*sc),round(B[2,2]*ln*sc)]); }
{ conemin(B,lz,ln,R) = my(best=0,bv=0);
 for(i=-R,R, my(lo=-oo,hi=oo,ok=1);
  for(r=1,2, my(a=i*B[r,1],b=B[r,2]);
    if(b==0, if(a<0,ok=0), my(t=-a/b); if(b>0,if(t>lo,lo=t),if(t<hi,hi=t))));
  if(ok, my(j0=if(lo==-oo,-R,ceil(lo)), j1=if(hi==oo,R,floor(hi)));
   if(j0<=j1, for(w=0,1, my(j=if(w==0,j0,j1), cz=i*B[1,1]+j*B[1,2], cn=i*B[2,1]+j*B[2,2]);
     if(cz||cn, my(v=cz*lz+cn*ln); if(best==0||v<best,best=v;bv=[cz,cn]))))));
 [best,bv]; }
k=22.4; R=600; GG=Catalan;
print("n | signs(lam1 vec in scaled coords) | cone vec (cz,cn) raw signs | log|cz|/n log|cn|/n | v2(cz) v2(cn) | log(cn/cz) | ratio*lamN/lamZ");
{
for(n=4,44,
  my(D=lcm(vector(6*n,i,i)), S=D^2, T=2^floor(k*n), M=S*T,
     zr=zudrow(n), nr=nestrow(n), X=zr[1],Y=zr[2],V=nr[1],U=nr[2]);
  my(B0=klat2(X,Y,V,U,M), idx=abs(matdet(B0)));
  my(LZ=(X*GG-Y)/M, LN=(V*GG-U)/M, sZ=sign(LZ), sN=sign(LN), lz=abs(LZ), ln=abs(LN));
  my(B=[sZ*B0[1,1],sZ*B0[1,2]; sN*B0[2,1],sN*B0[2,2]], Br=redu(B,lz,ln));
  my(nrm(v)=sqrt((v[1]*lz)^2+(v[2]*ln)^2));
  my(v1=[Br[1,1],Br[2,1]], v2=[Br[1,2],Br[2,2]]);
  if(nrm(v1)>nrm(v2), my(t=v1); v1=v2; v2=t);
  my(r=conemin(Br,lz,ln,R), cv=r[2]);
  /* cv is in scaled (u,v) sign coords: cv[1]=sZ*cZ, cv[2]=sN*cN >=0 */
  my(cZ=sZ*cv[1], cN=sN*cv[2]);
  printf("%3d | lam1 signs (%2d,%2d) | cone (sZcZ,sNcN)=(%d,%d) | %8.4f %8.4f | %3d %3d | %9.5f | %9.5f\n",
    n, sign(v1[1]), sign(v1[2]), sign(cv[1]),sign(cv[2]),
    if(cZ,log(abs(cZ))/n,0), if(cN,log(abs(cN))/n,0),
    if(cZ,valuation(cZ,2),-1), if(cN,valuation(cN,2),-1),
    if(cZ&&cN, log(abs(cN*1.0/cZ)),0),
    if(cZ&&cN, log(abs(cN*ln/(cZ*lz))),0) );
); }
\q
