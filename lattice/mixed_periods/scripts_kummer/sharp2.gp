default(parisize,"4G");
N=200;
Lt=vector(N+2);Lt[1]=1;for(n=1,N+1,Lt[n+1]=lcm(Lt[n],n));
L(n)=if(n<=0,1,Lt[n+1]);
intg(s)={my(r=O(x^(N+2)));for(j=0,N, r+=polcoeff(s,j)*x^(j+1)/(j+1));r;}
bad(s)={my(v=List());for(n=0,N,my(dn=denominator(polcoeff(s,n)),e=dn/gcd(dn,L(n)));if(e>1,listput(v,[n,e])));Vec(v);}
lg=log(1-x+O(x^(N+3)));
print("all (n,excess) with excess>1 vs type [1..n], kernel log(1-t), C family");
{for(k=3,6, for(m=1,3,
  my(H=(1-k^k*m*x+O(x^(N+3)))^(-1/k));
  for(a=1,k-1, my(b=bad(H^a*intg(H^(k-a)*lg)));
    if(#b>0, print("  k=",k," m=",m," a=",a,"  -> ",b)))));}
quit;
