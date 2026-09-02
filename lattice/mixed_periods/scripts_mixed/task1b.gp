default(parisize,"2G");
N=400;
Lt=vector(N+1);Lt[1]=1;for(n=1,N,Lt[n+1]=lcm(Lt[n],n));
L(n)=Lt[n+1];
intg(s)=my(r=O(x^(N+2)));for(k=0,N, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
lg=log(1-x+O(x^(N+2)));
HA=1/sqrt(1-4*x+O(x^(N+2)));
HC=HA*intg(HA*lg/x);
print("n, log(e_n)/n  and  e_n / prod_{2n/3<p<=n} p");
{for(i=1,12, my(n=[10,20,40,80,120,160,200,240,280,320,360,400][i],
   dn=denominator(polcoeff(HC,n)), e=dn/gcd(dn,L(n)*L(n\2)),
   P=1); forprime(p=2*n\3+1,n,P*=p);
   print("  n=",n,"  log(e_n)/n=",log(e*1.0)/n,"  e_n/P = ",e/P, "  (P=prod primes in (2n/3,n])"));}
print();
print("Is e_n always divisible by prod_{2n/3<p<=n} p, and what is the cofactor? n=1..60:");
{for(n=1,60, my(dn=denominator(polcoeff(HC,n)), e=dn/gcd(dn,L(n)*L(n\2)), P=1);
   forprime(p=2*n\3+1,n,P*=p);
   print1("(",n,",",if(e%P==0,e/P,"NOTDIV"),") "));}
print();
quit;
