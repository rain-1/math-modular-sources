default(parisize,"4G");
N=200;
Lt=vector(N+2);Lt[1]=1;for(n=1,N+1,Lt[n+1]=lcm(Lt[n],n));
L(n)=if(n<=0,1,Lt[n+1]);
tB(n)=L(n);
intg(s)={my(r=O(x^(N+2)));for(j=0,N, r+=polcoeff(s,j)*x^(j+1)/(j+1));r;}
exc(s,tf)={my(mx=1,arg0=-1);for(n=0,N,my(dn=denominator(polcoeff(s,n)),T=tf(n),e=dn/gcd(dn,T));if(e>mx,mx=e;arg0=n));[mx,arg0];}
lg=log(1-x+O(x^(N+3)));
print("kernel log(1-t): excess of H^(a)[log(1-t)] against type [1..n] alone");
{for(k=2,6, for(m=1,3, forstep(sg=1,-1,-2,
  my(H=(1-sg*k^k*m*x+O(x^(N+3)))^(-1/k));
  my(v=vector(k-1,a,my(e=exc(H^a*intg(H^(k-a)*lg),tB));Str(a,":",e[1],"@",e[2])));
  print("  k=",k," m=",m," fam=",if(sg==1,"C","R"),"  ",v))));}
quit;
