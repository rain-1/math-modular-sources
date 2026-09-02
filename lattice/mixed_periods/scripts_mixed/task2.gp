default(parisize,"4G");
N=200;
Lt=vector(N+1);Lt[1]=1;for(n=1,N,Lt[n+1]=lcm(Lt[n],n));
L(n)=Lt[n+1];
TY(n,k)=if(k==0,1,if(k==1,L(n),if(k==2,L(n)*L(n\2),L(n)^2)));
intg(s)=my(r=O(x^(N+2)));for(k=0,N, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
exc(s,k)=my(mx=1);for(n=1,N,my(dn=denominator(polcoeff(s,n)));mx=max(mx,dn/gcd(dn,TY(n,k))));mx;
cls(s)=my(v=vector(4,k,exc(s,k-1)));if(v[1]==1,"1",if(v[2]==1,"[n]",if(v[3]==1,"[n][n/2]",if(v[4]==1,"[n]^2","WORSE"))));
firstfail(s,k)=my(f=0);for(n=1,N,my(dn=denominator(polcoeff(s,n)));if(dn/gcd(dn,TY(n,k))>1,f=n;break));f;
lg=log(1-x+O(x^(N+2)));
{for(mi=1,2, my(m=[1,2][mi], HA=1/sqrt(1-4*m*x+O(x^(N+2))), HA2=1/(1-4*m*x+O(x^(N+2))));
print("############ m=",m," (D'=",4*m-1,")  N=",N);
print("outer | i | j | l | s | class | exc/[n][n/2] | exc/[n]^2 | 1st fail [n][n/2]");
for(oi=1,3,
 my(OUT=[HA,HA^3,HA/(1-x)][oi], onm=["H_A","H_A^3","H_A/(1-x)"][oi]);
 for(i=-1,2, for(j=0,3, for(l=0,2, for(s=0,1,
   if(i==-1 && l==0, next);
   my(k = x^i*(1-x+O(x^(N+2)))^(-j)*lg^l*HA2^s);
   my(F = OUT*intg(HA*k));
   my(e12=exc(F,2), e2=exc(F,3));
   print(onm," | ",i," | ",j," | ",l," | ",s," | ",cls(F)," | ",if(e12>10^12,"BIG",e12)," | ",if(e2>10^12,"BIG",e2)," | ",firstfail(F,2));
 )))));
);}
quit;
