s2(x)=hammingweight(x);
{
klat2(X,Y,V,U,M) = my(K=matkerint([X,V,M,0,0,0; Y,U,0,M,0,0]));
  mathnf(matrix(2,#K[1,],i,j,K[i,j]));
}
{
redu(B,lz,ln) = my(sc=10^(400+ceil(log(1/min(lz,ln))/log(10))));
  B*qflll([round(B[1,1]*lz*sc),round(B[1,2]*lz*sc); round(B[2,1]*ln*sc),round(B[2,2]*ln*sc)]);
}
/* cancellation excess in the q-numerator */
{
run(nlist,k,GG,G2)=
 printf("  n | v2(M) | eps(l1)  eps/n | 24n-1-4s2(3n) | margin(l1) | 24n-eps | v2(q1) | logL1/n\n");
 for(t=1,#nlist, my(n=nlist[t], D=lcm(vector(6*n,i,i)), S=D^2, T=2^floor(k*n), M=S*T,
   zr=zudrow(n), nr=nestrow(n), X=zr[1],Y=zr[2],V=nr[1],U=nr[2]);
  my(B0=klat2(X,Y,V,U,M), LZ=(X*GG-Y)/M, LN=(V*GG-U)/M,
     sZ=sign(LZ),sN=sign(LN),lz=abs(LZ),ln=abs(LN));
  my(B=[sZ*B0[1,1],sZ*B0[1,2]; sN*B0[2,1],sN*B0[2,2]], Br=redu(B,lz,ln));
  my(cz=sZ*Br[1,1], cn=sN*Br[2,1]);
  my(num=cz*X+cn*V, q=num/M, p=(cz*Y+cn*U)/M);
  my(eps=v2(num)-min(v2(cz*X),v2(cn*V)), qZ=24*n-1-4*s2(3*n));
  printf("%4d | %5d | %6d %7.3f | %13d | %10d | %7d | %6d | %8.4f\n",
    n, valuation(M,2), eps, 1.0*eps/n, qZ, v2(G2-p/q), qZ-eps, v2(q), log(abs(cz*lz+cn*ln))/n);
 );
}
G2=G2rat(140)[1];
run([4,6,8,10,12,14,16,20],22.4,Catalan,G2);
print();
/* rank / trivial relations in the moment family */
for(n=1,4, my(N=3*n, A=vector(N+1),B=vector(N+1));
  for(j=0,N, my(r=nestgen(3*n,j)); A[j+1]=r[1]; B[j+1]=r[2]);
  my(Mx=matrix(2,N+1,i,j,if(i==1,A[j],B[j])));
  print("n=",n,"  family size=",N+1,"  rank over Q=",matrank(Mx),
        "  dim of identically-vanishing relation space=",N+1-matrank(Mx)));
\q
