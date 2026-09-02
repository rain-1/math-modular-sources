default(parisize,"2G");
N=400;
Lt=vector(N+1);Lt[1]=1;for(n=1,N,Lt[n+1]=lcm(Lt[n],n));
L(n)=Lt[n+1];
TY(n,k)=if(k==0,1,if(k==1,L(n),if(k==2,L(n)*L(n\2),L(n)^2)));
intg(s)=my(r=O(x^(N+2)));for(k=0,N, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
exc(s,k)=my(mx=1);for(n=1,N,my(dn=denominator(polcoeff(s,n)));mx=max(mx,dn/gcd(dn,TY(n,k))));mx;
excs(s,k,K)=vector(K,n,my(dn=denominator(polcoeff(s,n)));dn/gcd(dn,TY(n,k)));
lg=log(1-x+O(x^(N+2)));
{for(j=1,10,
  my(ms=[1,2,3,4,5,6,7,11,17,41], m=ms[j],
     HA=1/sqrt(1-4*m*x+O(x^(N+2))),
     HB=HA*intg(HA/(1-x)),
     HC=HA*intg(HA*lg/x),
     HD=HA*intg(HA*lg/(1-x)));
  print("m=",m," D'=",4*m-1);
  print("  HA exc/1     = ",exc(HA,0));
  print("  HB exc/1     = ",exc(HB,0),"   exc/[n] = ",exc(HB,1),"   exc/[n][n/2] = ",exc(HB,2));
  print("  HD exc/[n]   = ",exc(HD,1),"   exc/[n][n/2] = ",exc(HD,2),"   exc/[n]^2 = ",exc(HD,3));
  print("  HC exc/[n][n/2] = ",exc(HC,2),"   exc/[n]^2 = ",exc(HC,3));
);}
print();
HA=1/sqrt(1-4*x+O(x^(N+2)));
HC=HA*intg(HA*lg/x);
HD=HA*intg(HA*lg/(1-x));
print("m=1 HC excess over [n][n/2], n=1..40:");
print(excs(HC,2,40));
print("m=1 HC excess over [n][n/2], n=41..80:");
print(vector(40,i,my(n=40+i,dn=denominator(polcoeff(HC,n)));dn/gcd(dn,TY(n,2))));
print("m=1 HC excess factored, n=1..48 (only e>1):");
for(n=1,48,my(dn=denominator(polcoeff(HC,n)),e=dn/gcd(dn,TY(n,2)));if(e>1,print1("n=",n,": ",e," = ",factor(e),"  ")));
print();
print("m=1 HC excess at n=2^k and nearby, k=1..8:");
for(k=1,8,my(n=2^k,dn=denominator(polcoeff(HC,n)),e=dn/gcd(dn,TY(n,2)));print("  n=",n," e=",e," factored ",factor(e)));
print();
print("m=1 HD table: n, den(a_n), L(n)L(n/2), L(n)L(n/2)/den");
for(n=1,30,my(d=denominator(polcoeff(HD,n)));print("  ",n,"  ",d,"  ",TY(n,2),"  ",TY(n,2)/d));
quit;
