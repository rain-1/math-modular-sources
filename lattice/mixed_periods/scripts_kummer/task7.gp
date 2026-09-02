default(parisize,"4G");
N=200;
Lt=vector(N+2);Lt[1]=1;for(n=1,N+1,Lt[n+1]=lcm(Lt[n],n));
L(n)=if(n<=0,1,Lt[n+1]);
tB(n)=L(n);
tD(n)=L(n)*L(n\2);
intg(s)={my(r=O(x^(N+2)));for(j=0,N, r+=polcoeff(s,j)*x^(j+1)/(j+1));r;}
exc(s,tf)={my(mx=1,arg0=-1);for(n=0,N,my(dn=denominator(polcoeff(s,n)),T=tf(n),e=dn/gcd(dn,T));if(e>mx,mx=e;arg0=n));[mx,arg0];}
isint(s)={for(n=0,N,if(denominator(polcoeff(s,n))!=1,return(n)));-1;}
rad(k)=my(r=1);fordiv(k,d,if(isprime(d),r*=d));r;
lg=log(1-x+O(x^(N+3)));
print("=== 7a. integrality of (1 -+ N x)^{-j/k}, N = k*rad(k)*m, j=1..k-1 (-1 = integral for n<=",N,") ===");
{for(k=2,7, my(K=k*rad(k)); for(m=1,4, forstep(sg=1,-1,-2,
  my(H=(1-sg*K*m*x+O(x^(N+3)))^(-1/k));
  print("  k=",k," k*rad(k)=",K," N=",K*m," fam=",if(sg==1,"M","P"),"  ",vector(k-1,j,isint(H^j))))));}
print();
print("=== 7b. SHARPNESS: for each k, smallest N>0 with (1-Nx)^{-j/k} integral for all j=1..k-1 (scan N=1..k*rad(k)) ===");
{for(k=2,7, my(K=k*rad(k), good=List());
  for(NN=1,K, my(H=(1-NN*x+O(x^62))^(-1/k), ok=1);
    for(j=1,k-1, for(n=0,60, if(denominator(polcoeff(H^j,n))!=1, ok=0;break(2))));
    if(ok, listput(good,NN)));
  print("  k=",k,"  k*rad(k)=",K,"   N in 1..",K," that work (minus fam): ",Vec(good)));}
print();
{for(k=2,7, my(K=k*rad(k), good=List());
  for(NN=1,K, my(H=(1+NN*x+O(x^62))^(-1/k), ok=1);
    for(j=1,k-1, for(n=0,60, if(denominator(polcoeff(H^j,n))!=1, ok=0;break(2))));
    if(ok, listput(good,NN)));
  print("  k=",k,"  k*rad(k)=",K,"   N in 1..",K," that work (plus fam): ",Vec(good)));}
print();
print("=== 7c. named counterexamples ===");
{my(v=[[3,3,1],[3,3,-1],[4,4,1],[4,4,-1],[4,2,1],[5,5,1],[6,6,1],[6,12,1],[6,18,1]]);
 for(i=1,#v, my(k=v[i][1],NN=v[i][2],sg=v[i][3], H=(1-sg*NN*x+O(x^62))^(-1/k));
   print("  (1",if(sg==1,"-","+"),NN,"x)^{-j/",k,"}: first non-integral n for j=1..k-1: ",vector(k-1,j,isint(H^j))));}
print();
print("=== 7d. denominator types at general N=k*rad(k)*m (excess vs [1..n] for B, [1..n][1..n/2] for D,Lk) ===");
{for(k=3,6, my(K=k*rad(k)); for(m=1,3, forstep(sg=1,-1,-2,
  my(H=(1-sg*K*m*x+O(x^(N+3)))^(-1/k));
  my(v=vector(k-1,a,my(Hk=H^(k-a),Ha=H^a,B=Ha*intg(Hk/(1-x)),Dd=Ha*intg(Hk*lg/(1-x)),Ll=Ha*intg(Hk*lg));
     Str("a",a,":B=",exc(B,tB)[1],",D=",exc(Dd,tD)[1],",L=",exc(Ll,tD)[1])));
  print("  k=",k," N=",K*m," fam=",if(sg==1,"M","P"),"  ",v))));}
quit;
