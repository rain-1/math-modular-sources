default(parisize,"6G");
N=300;
Lt=vector(N+2);Lt[1]=1;for(n=1,N+1,Lt[n+1]=lcm(Lt[n],n));
L(n)=if(n<=0,1,Lt[n+1]);
tB(n)=L(n);
tD(n)=L(n)*L(n\2);
intg(s)={my(r=O(x^(N+2)));for(j=0,N, r+=polcoeff(s,j)*x^(j+1)/(j+1));r;}
exc(s,tf)={my(mx=1,arg0=-1);for(n=0,N,my(dn=denominator(polcoeff(s,n)),T=tf(n),e=dn/gcd(dn,T));if(e>mx,mx=e;arg0=n));[mx,arg0];}
isint(s)={for(n=0,N,if(denominator(polcoeff(s,n))!=1,return(n)));-1;}
lg=log(1-x+O(x^(N+3)));
run(k,m,sgn)={
  my(c=k^k*m, H=(1-sgn*c*x+O(x^(N+3)))^(-1/k));
  my(bad=vector(k-1,j,isint(H^j)));
  print("k=",k," m=",m," fam=",if(sgn==1,"C","R"),"  H^j integral j=1..k-1 (-1=integral): ",bad);
  for(a=1,k-1,
    my(Hk=H^(k-a), Ha=H^a);
    my(B=Ha*intg(Hk/(1-x)), Dd=Ha*intg(Hk*lg/(1-x)), Ll=Ha*intg(Hk*lg));
    my(eB=exc(B,tB), eD=exc(Dd,tD), eL=exc(Ll,tD));
    print("    a=",a,"  B exc[1..n]=",eB[1],"@n=",eB[2],"   D exc[1..n][1..n/2]=",eD[1],"@n=",eD[2],"   Lk exc[1..n][1..n/2]=",eL[1],"@n=",eL[2]);
  );
}
{for(k=2,7, for(m=1,3, run(k,m,1); run(k,m,-1)); print(""));}
quit;
