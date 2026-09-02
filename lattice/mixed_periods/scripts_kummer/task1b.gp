default(parisize,"4G");
N=20;
Lt=vector(N+2);Lt[1]=1;for(n=1,N+1,Lt[n+1]=lcm(Lt[n],n)); L(n)=if(n<=0,1,Lt[n+1]);
intg(s)=my(r=O(x^(N+2)));for(j=0,N, r+=polcoeff(s,j)*x^(j+1)/(j+1));r;
exc(s,typ)={my(mx=1,arg=-1);for(n=0,N,my(dn=denominator(polcoeff(s,n)),T=typ(n),e=dn/gcd(dn,T));if(e>mx,mx=e;arg=n));[mx,arg];}
isint(s)={for(n=0,N,if(denominator(polcoeff(s,n))!=1,return(n)));-1}
lg=log(1-x+O(x^(N+3)));
tB(n)=L(n); tD(n)=L(n)*L(n\2);
run(k,m,sgn)={
  my(c=k^k*m, H=(1-sgn*c*x+O(x^(N+3)))^(-1/k));
  my(bad=vector(k-1,j,isint(H^j)));
  my(res=List());
  for(a=1,k-1,
    my(Hk=H^(k-a), Ha=H^a);
    my(B=Ha*intg(Hk/(1-x)), Dd=Ha*intg(Hk*lg/(1-x)), Ll=Ha*intg(Hk*lg));
    listput(res,[a,exc(B,tB),exc(Dd,tD),exc(Ll,tD)]);
  );
  print("k=",k," m=",m," fam=",if(sgn==1,"C","R"),"  H^j integral 1..k-1 (-1=ok, else first bad n): ",bad);
  for(i=1,#res, my(r=res[i]); print("    a=",r[1],"  B[1..n] exc=",r[2][1],"@n=",r[2][2],"   D[1..n][1..n/2] exc=",r[3][1],"@n=",r[3][2],"   L[1..n][1..n/2] exc=",r[4][1],"@n=",r[4][2]));
};
{run(2,1,1);}
quit;
