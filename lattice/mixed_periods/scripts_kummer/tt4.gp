Lt=vector(20);Lt[1]=1;for(n=1,19,Lt[n+1]=lcm(Lt[n],n)); L(n)=if(n<=0,1,Lt[n+1]);
exc(s,typ)={my(mx=1,arg=-1);for(n=0,3,my(dn=denominator(polcoeff(s,n)),T=typ(n),e=dn/gcd(dn,T));if(e>mx,mx=e;arg=n));[mx,arg];}
tB(n)=L(n);
s=1/2+x/3+O(x^5);
print(exc(s,tB));
quit;
