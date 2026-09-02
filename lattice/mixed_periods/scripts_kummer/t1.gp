default(parisize,"1G");
N=150;
Lt=vector(N+1);Lt[1]=1;for(n=1,N,Lt[n+1]=lcm(Lt[n],n)); L(n)=Lt[n+1];
intg(s)=my(r=O(x^(N+1)));for(j=0,N-1, r+=polcoeff(s,j)*x^(j+1)/(j+1));r;
exc(s,typ)={my(m=1);for(n=1,N-1,my(dn=denominator(polcoeff(s,n)),T=typ(n));m=max(m,dn/gcd(dn,T)));m;}
isint(s)={for(n=0,N-1,if(denominator(polcoeff(s,n))!=1,return(0)));1}
lg=log(1-x+O(x^(N+2)));
{for(k=2,6, for(m=1,2, my(H=(1-k^k*m*x+O(x^(N+2)))^(-1/k));
  my(ok=vector(k-1,j,isint(H^j)));
  my(HB=H*intg(H^(k-1)/(1-x)), HD=H*intg(H^(k-1)*lg/(1-x)), HL=H*intg(H^(k-1)*lg));
  print("k=",k," m=",m,"  H^j integral for j=1..k-1: ",ok, "   HB excess/[n]: ",exc(HB,n->L(n)), "  HD excess/[n][n/2]: ",exc(HD,n->L(n)*L(n\2)), "  HL: ",exc(HL,n->L(n)*L(n\2)));
  \\ other sector a=2,b=k-2 (k>=3)
  if(k>=3, my(HB2=H^2*intg(H^(k-2)/(1-x)), HD2=H^2*intg(H^(k-2)*lg/(1-x))); print("      sector (2,",k-2,"): HB2/[n]: ",exc(HB2,n->L(n)),"  HD2/[n][n/2]: ",exc(HD2,n->L(n)*L(n\2))));
));}
quit;
