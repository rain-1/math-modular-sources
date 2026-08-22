default(parisizemax, 6000000000);
/* xi for the level-5 Fricke row to high precision, then identification attempts. */
{
default(realprecision, 260);
N=1200;
a=vector(N+2); b=vector(N+2); a[1]=1;a[2]=6; b[1]=0;b[2]=1;
for(n=1,N, my(P=88*n^2+44*n+6, Q=-64*n^2+64*n-12);
  a[n+2]=(P*a[n+1]-Q*a[n])/(n+1)^2; b[n+2]=(P*b[n+1]-Q*b[n])/(n+1)^2);
W=vector(N+1,i,a[i]*b[i+1]-a[i+1]*b[i]);
S=sum(m=1,N, W[m]/(a[m]*a[m+1]));
xi = S*1.0;
print("xi = ", xi);
print("digits of agreement check (partial sums): ", (S - sum(m=1,N-40, W[m]/(a[m]*a[m+1])))*1.0);
\\ ---- identification battery
cst = [1, Pi, Pi^2, Pi^3, zeta(3), log(2), log(5), log((1+sqrt(5))/2), sqrt(5), 1/sqrt(5)];
nm  = ["1","Pi","Pi^2","Pi^3","zeta(3)","log2","log5","log phi","sqrt5","1/sqrt5"];
for(i=1,#cst, my(v=lindep([xi,cst[i]],120)); if(v!=0 && vecmax(abs(v))<10^12, print("  HIT 2-term ",nm[i],": ",v~)));
\\ weight-3 newform L-values at levels with 5 | M
for(M=5,120, if(M%5==0,
  my(ok=1, mf);
  mf = mfinit([M,3],0);
  my(V=mfeigenbasis(mf));
  for(j=1,#V,
    my(L=lfunmf(mf,V[j]));
    if(type(L)=="t_VEC", L=L[1]);
    for(s=1,3,
      my(val=lfun(L,s));
      if(abs(imag(val))<1e-40, val=real(val));
      my(v=lindep([xi,real(val)],110));
      if(v!=0 && vecmax(abs(v))<10^9, print("  HIT level ",M," form ",j," s=",s,": ",v~)))));
);
print("battery done");
}
