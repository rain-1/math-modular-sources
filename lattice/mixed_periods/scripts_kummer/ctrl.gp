default(parisize,"2G");
N=120;
Lt=vector(N+2);Lt[1]=1;for(n=1,N+1,Lt[n+1]=lcm(Lt[n],n));
L(n)=if(n<=0,1,Lt[n+1]);
tB(n)=L(n);
tD(n)=L(n)*L(n\2);
one(n)=1;
intg(s)={my(r=O(x^(N+2)));for(j=0,N, r+=polcoeff(s,j)*x^(j+1)/(j+1));r;}
exc(s,tf)={my(mx=1,arg0=-1);for(n=0,N,my(dn=denominator(polcoeff(s,n)),T=tf(n),e=dn/gcd(dn,T));if(e>mx,mx=e;arg0=n));[mx,arg0];}
isint(s)={for(n=0,N,if(denominator(polcoeff(s,n))!=1,return(n)));-1;}
lg=log(1-x+O(x^(N+3)));
print("=== CONTROL A: wrong scaling c=m*x (not k^k m x): H^j should FAIL integrality ===");
{for(k=3,5, for(m=1,2, my(H=(1-m*x+O(x^(N+3)))^(-1/k));
  print("  k=",k," m=",m,"  first non-integral coeff index of H^j, j=1..k-1: ",vector(k-1,j,isint(H^j)))));}
print("=== CONTROL B: near-miss scaling c=k^k m -1 (=D) ===");
{for(k=3,5, for(m=1,2, my(H=(1-(k^k*m-1)*x+O(x^(N+3)))^(-1/k));
  print("  k=",k," m=",m,"  first non-integral coeff index of H^j: ",vector(k-1,j,isint(H^j)))));}
print("=== SHARPNESS: is D-type series of type [1..n] ALONE?  (excess vs L(n) only; >1 means L(n/2) genuinely needed) ===");
{for(k=3,4, for(m=1,2, my(H=(1-k^k*m*x+O(x^(N+3)))^(-1/k));
  for(a=1,k-1, my(Hk=H^(k-a),Ha=H^a);
    my(Dd=Ha*intg(Hk*lg/(1-x)), Ll=Ha*intg(Hk*lg));
    my(e1=exc(Dd,tB), e2=exc(Ll,tB), e3=exc(Dd,one));
    print("  k=",k," m=",m," a=",a,"  D vs [1..n]: ",e1[1],"@",e1[2],"   Lk vs [1..n]: ",e2[1],"@",e2[2],"   D vs 1 (raw den): ",e3[1],"@",e3[2]);
  )));}
print("=== SHARPNESS: is B-type series integral (type [1] )? ===");
{for(k=3,4, for(m=1,2, my(H=(1-k^k*m*x+O(x^(N+3)))^(-1/k));
  for(a=1,k-1, my(B=H^a*intg(H^(k-a)/(1-x))); my(e=exc(B,one));
    print("  k=",k," m=",m," a=",a,"  B raw denominator excess vs 1: ",e[1],"@n=",e[2]))));}
quit;
