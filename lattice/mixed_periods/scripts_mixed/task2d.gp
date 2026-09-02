default(parisize,"4G");
N=400;
Lt=vector(N+3);Lt[1]=1;for(n=1,N+2,Lt[n+1]=lcm(Lt[n],n));
L(n)=Lt[n+1];
intg(s)=my(r=O(x^(N+2)));for(k=0,N, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
excT(s,T)=my(mx=1);for(n=1,N,my(dn=denominator(polcoeff(s,n)));mx=max(mx,dn/gcd(dn,T(n))));mx;
lg=log(1-x+O(x^(N+2)));
print("kernel log(1-t)/t : /[n]=",excT(lg/x,n->L(n)),"  /[n+1]=",excT(lg/x,n->L(n+1)),"  /[n][n/2]=",excT(lg/x,n->L(n)*L(n\2)));
{for(mi=1,2, my(m=[1,2][mi], HA=1/sqrt(1-4*m*x+O(x^(N+2))));
print("### m=",m);
my(HC=HA*intg(HA*lg/x));
print("  H_C: /[n]^2=",excT(HC,n->L(n)^2), "  /[n][n/2][n/3]=",excT(HC,n->L(n)*L(n\2)*L(n\3)), "  /[n][n/2][n/3][n/4]=",excT(HC,n->L(n)*L(n\2)*L(n\3)*L(n\4)));
for(j=0,3, my(k=(1-x+O(x^(N+2)))^(-j)*lg^2, F=HA*intg(HA*k));
  print("  H[log^2/(1-t)^",j,"]: /[n][n/2]=",excT(F,n->L(n)*L(n\2)),"  /[n][n/2][n/3]=",excT(F,n->L(n)*L(n\2)*L(n\3)),"  /[n]^2[n/2]=",excT(F,n->L(n)^2*L(n\2)),"  /[n]^2=",excT(F,n->L(n)^2)));
for(i=1,2,for(j=0,3, my(k=x^i*(1-x+O(x^(N+2)))^(-j)*lg^2, F=HA*intg(HA*k));
  print("  H[t^",i,"log^2/(1-t)^",j,"]: /[n][n/2]=",excT(F,n->L(n)*L(n\2)),"  /[n][n/2][n/3]=",excT(F,n->L(n)*L(n\2)*L(n\3)))));
);}
quit;
