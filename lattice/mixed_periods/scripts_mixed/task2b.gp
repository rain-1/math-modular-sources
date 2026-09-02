default(parisize,"4G");
N=200;
Lt=vector(N+1);Lt[1]=1;for(n=1,N,Lt[n+1]=lcm(Lt[n],n));
L(n)=Lt[n+1];
TY(n,k)=if(k==0,1,if(k==1,L(n),if(k==2,L(n)*L(n\2),L(n)^2)));
intg(s)=my(r=O(x^(N+2)));for(k=0,N, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
exc(s,k,lo,hi)=my(mx=1);for(n=lo,hi,my(dn=denominator(polcoeff(s,n)));mx=max(mx,dn/gcd(dn,TY(n,k))));mx;
cls(s)=my(v=vector(4,k,exc(s,k-1,1,N)));if(v[1]==1,"1",if(v[2]==1,"[n]",if(v[3]==1,"[n][n/2]",if(v[4]==1,"[n]^2","WORSE"))));
sz(e)=if(e==1,"1",Str("10^",round(log(e*1.0)/log(10.))));
lg=log(1-x+O(x^(N+2)));
{for(mi=1,2, my(m=[1,2][mi], HA=1/sqrt(1-4*m*x+O(x^(N+2))), HA2=1/(1-4*m*x+O(x^(N+2))));
print("### m=",m,"  N=",N);
print("i j l s | kernel_class | H[k]_class | H'[k]_class | H''[k]_class | exc12@n<=50 | exc12@n<=100 | exc12@n<=200 | exc2@n<=200");
for(i=-1,2, for(j=0,3, for(l=0,2, for(s=0,1,
   if(i==-1 && l==0, next);
   my(k = x^i*(1-x+O(x^(N+2)))^(-j)*lg^l*HA2^s);
   my(F1 = HA*intg(HA*k), F2 = HA^3*intg(HA*k), F3 = (HA/(1-x))*intg(HA*k));
   print(i," ",j," ",l," ",s," | ",cls(k)," | ",cls(F1)," | ",cls(F2)," | ",cls(F3)," | ",
     sz(exc(F1,2,1,50))," | ",sz(exc(F1,2,1,100))," | ",sz(exc(F1,2,1,200))," | ",sz(exc(F1,3,1,200)));
)))));}
quit;
