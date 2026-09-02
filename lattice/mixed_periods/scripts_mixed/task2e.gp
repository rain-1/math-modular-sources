default(parisize,"4G");
N=400;
Lt=vector(N+3);Lt[1]=1;for(n=1,N+2,Lt[n+1]=lcm(Lt[n],n));
L(n)=Lt[n+1];
intg(s)=my(r=O(x^(N+2)));for(k=0,N, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
excT(s,T)=my(mx=1);for(n=1,N,my(dn=denominator(polcoeff(s,n)));mx=max(mx,dn/gcd(dn,T(n))));mx;
lg=log(1-x+O(x^(N+2)));
{for(mi=1,2, my(m=[1,2][mi], HA=1/sqrt(1-4*m*x+O(x^(N+2))));
print("### m=",m," N=400: kernels t^-1 (1-t)^-j log(1-t)");
for(j=0,3, my(k=lg/(x*(1-x+O(x^(N+2)))^j), F=HA*intg(HA*k));
  print("  j=",j,": kernel /[n]=",excT(k,n->L(n)),"  kernel /[n+1]=",excT(k,n->L(n+1)),
        "   H[k] /[n]^2=",excT(F,n->L(n)^2),"   H[k] /[n][n/2]=",if(excT(F,n->L(n)*L(n\2))>10^12,"BIG",excT(F,n->L(n)*L(n\2))))));}
quit;
