/* the new row (al,be,ga,de,ep,ze) = (234,-39,-78,-27,63,-30):
   xi to 200 digits, then identification attempts.                         */
default(parisizemax, 6000000000);
{
default(realprecision, 210);
N=900;
a=vector(N+2); b=vector(N+2); a[1]=1;a[2]=-78; b[1]=0;b[2]=1;
for(n=1,N, my(P=234*n^2-39*n-78, Q=-27*n^2+63*n-30);
  a[n+2]=(P*a[n+1]-Q*a[n])/(n+1)^2; b[n+2]=(P*b[n+1]-Q*b[n])/(n+1)^2);
W=vector(N+1,i,a[i]*b[i+1]-a[i+1]*b[i]);
S=sum(m=1,N, W[m]/(a[m]*a[m+1]));
xi=S*1.0;
print("xi = ",xi);
cst=[1,Pi,Pi^2,Pi^3,Pi^4,zeta(3),zeta(5),log(2),log(3),log(13),log(127),
     Pi*log(3),Pi^2*log(3),log(2)^2,log(3)^2,Catalan,Pi/sqrt(3),Pi^3/sqrt(3),
     sqrt(3),sqrt(381),sqrt(127),gamma(1/3),gamma(1/6),gamma(1/3)^6,gamma(1/6)^6];
nm=["1","Pi","Pi^2","Pi^3","Pi^4","z3","z5","log2","log3","log13","log127",
    "Pi log3","Pi^2 log3","log2^2","log3^2","G","Pi/r3","Pi^3/r3","r3","r381","r127",
    "gam(1/3)","gam(1/6)","gam(1/3)^6","gam(1/6)^6"];
for(i=1,#cst, my(v=lindep([xi,cst[i]],150)); if(v!=0 && vecmax(abs(v))<10^12, print("  2-term HIT ",nm[i],": ",v~)));
for(i=2,#cst, my(v=lindep([xi,cst[i],1],110)); if(v!=0 && vecmax(abs(v))<10^8, print("  3-term HIT ",nm[i],": ",v~)));
\\ weight 3 and 4 newforms, small levels
for(w=3,4, for(M=1,80,
  my(mf=mfinit([M,w],0), V=mfeigenbasis(mf));
  for(j=1,#V, my(L=lfunmf(mf,V[j])); if(type(L)=="t_VEC", L=L[1]);
    for(s=1,w, my(val=lfun(L,s)); if(abs(imag(val))>1e-40, val=real(val));
      my(v=lindep([xi,real(val)],100));
      if(v!=0 && vecmax(abs(v))<10^8, print("  HIT w=",w," level ",M," form ",j," s=",s,": ",v~))))));
print("battery done");
}
