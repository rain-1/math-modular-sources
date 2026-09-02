default(parisize,"4G");
N=400;
Lt=vector(N+1);Lt[1]=1;for(n=1,N,Lt[n+1]=lcm(Lt[n],n));
L(n)=Lt[n+1];
intg(s)=my(r=O(x^(N+2)));for(k=0,N, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
excT(s,T)=my(mx=1);for(n=1,N,my(dn=denominator(polcoeff(s,n)));mx=max(mx,dn/gcd(dn,T(n))));mx;
lg=log(1-x+O(x^(N+2)));
{for(mi=1,2, my(m=[1,2][mi], HA=1/sqrt(1-4*m*x+O(x^(N+2))));
print("### m=",m," N=",N);
print("-- bounded-excess anomalies (kernel t^i log(1-t)^2 (1-t)^-j), excess over [n][n/2] up to N=400:");
for(i=1,2, for(j=0,3,
  my(k=x^i*(1-x+O(x^(N+2)))^(-j)*lg^2, F=HA*intg(HA*k));
  print("   i=",i," j=",j," excess=",excT(F,n->L(n)*L(n\2)))));
print("-- what type works for k=log(1-t)^2/(1-t)^j, j>=1?  (i=0)");
for(j=1,3, my(k=(1-x+O(x^(N+2)))^(-j)*lg^2, F=HA*intg(HA*k));
  print("   j=",j,":  /[n]^2 = ",excT(F,n->L(n)^2),
        "   /[n]^2[n/2] = ",excT(F,n->L(n)^2*L(n\2)),
        "   /[n][n/2][n/3] = ",excT(F,n->L(n)*L(n\2)*L(n\3)),
        "   /[n]^3 = ",excT(F,n->L(n)^3)));
print("-- kernel own types: log(1-t)^2 and log(1-t)^2/(1-t):");
print("   log^2      : /[n]=",excT(lg^2,n->L(n)),"  /[n][n/2]=",excT(lg^2,n->L(n)*L(n\2)));
print("   log^2/(1-t): /[n][n/2]=",excT(lg^2/(1-x+O(x^(N+2))),n->L(n)*L(n\2)));
print("   log/t      : /[n]=",excT(lg/x,n->L(n)),"  /[n+1]=",excT(lg/x,n->L(n+1)),"  /[n]^2=",excT(lg/x,n->L(n)^2)));}
quit;
