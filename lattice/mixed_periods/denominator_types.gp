default(parisize,"1G");
N=250;
Lt=vector(N+1);Lt[1]=1;for(n=1,N,Lt[n+1]=lcm(Lt[n],n)); L(n)=Lt[n+1];
intg(s)=my(r=O(x^(N+1)));for(k=0,N-1, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
exc(s,typ)={my(m=1);for(n=1,N-1,my(dn=denominator(polcoeff(s,n)),T=typ(n));m=max(m,dn/gcd(dn,T)));m;}
lg=log(1-x+O(x^(N+2)));
{for(m=1,8, my(HA=1/sqrt(1-4*m*x+O(x^(N+2))), HD=HA*intg(HA*lg/(1-x)), HB=HA*intg(HA/(1-x)), HC=HA*intg(HA*lg/x));
  print("m=",m,": HD excess over [n][n/2]: ",exc(HD,n->L(n)*L(n\2)), "   HB excess over [n]: ",exc(HB,n->L(n)), "   HC excess over [n][n/2]: ",exc(HC,n->L(n)*L(n\2)), "  HC over [n]^2: ",exc(HC,n->L(n)^2)));}
\\ structure at m=1
HA=1/sqrt(1-4*x+O(x^(N+2))); HD=HA*intg(HA*lg/(1-x)); HC=HA*intg(HA*lg/x);
print("m=1 HD: n, den, [n][n/2]/den:");
for(n=1,24, my(d=denominator(polcoeff(HD,n)));print1("(",n,",",d,",",L(n)*L(n\2)/d,") "));print();
print("m=1 HC: n, den, den/([n][n/2]) :");
for(n=1,24, my(d=denominator(polcoeff(HC,n)));print1("(",n,",",d,",",d/gcd(d,L(n)*L(n\2)),") "));print();
print("m=1 HC numerators*[n]^2 mod small primes? factor of den/([n][n/2]) at n=2^k: ");
for(k=1,7, my(n=2^k,d=denominator(polcoeff(HC,n)));print1(n,":",factor(d/gcd(d,L(n)*L(n\2)))," "));print();
quit;
