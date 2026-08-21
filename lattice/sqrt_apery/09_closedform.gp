default(parisizemax,6000000000);
N=14; a=vector(N+2); a[1]=1;a[2]=10;
{for(n=1,N, a[n+2]=((136*n^2+68*n+10)*a[n+1]-4*(2*n-1)^2*a[n])/(n+1)^2);}
tgt=vector(9,i,a[i]); print("target a_n: ",tgt);
v=vector(9,i,a[i]/binomial(2*(i-1),i-1)); print("v_n=a_n/binom(2n,n): ",v);
\\ brute search: a_n = sum_k  prod of binomials with small exponents
B=[ (n,k)->binomial(n,k), (n,k)->binomial(n+k,k), (n,k)->binomial(2*k,k), (n,k)->binomial(2*n-2*k,n-k), (n,k)->binomial(2*n,2*k), (n,k)->binomial(n,k)^2, (n,k)->binomial(n+k,k)^2, (n,k)->binomial(2*k,k)^2, (n,k)->binomial(2*n-2*k,n-k)^2, (n,k)->binomial(2*k,k)*binomial(2*n-2*k,n-k), (n,k)->4^k, (n,k)->4^(n-k), (n,k)->16^k, (n,k)->2^k, (n,k)->(-1)^k, (n,k)->binomial(2*n-k,k), (n,k)->binomial(n,k)*binomial(n+k,k) ];
nm=["C(n,k)","C(n+k,k)","C(2k,k)","C(2n-2k,n-k)","C(2n,2k)","C(n,k)^2","C(n+k,k)^2","C(2k,k)^2","C(2n-2k,n-k)^2","C(2k,k)C(2n-2k,n-k)","4^k","4^(n-k)","16^k","2^k","(-1)^k","C(2n-k,k)","C(n,k)C(n+k,k)"];
found=0;
{for(i1=1,#B, for(i2=i1,#B, for(i3=i2,#B,
  my(ok=1);
  for(n=0,6, my(s=sum(k=0,n, B[i1](n,k)*B[i2](n,k)*B[i3](n,k)));
     if(s!=tgt[n+1], ok=0; break));
  if(ok, found++; print("HIT: ",nm[i1]," * ",nm[i2]," * ",nm[i3])))));}
print("triples found: ",found);
\\ pairs only
{for(i1=1,#B, for(i2=i1,#B,
  my(ok=1);
  for(n=0,6, my(s=sum(k=0,n, B[i1](n,k)*B[i2](n,k)));
     if(s!=tgt[n+1], ok=0; break));
  if(ok, print("HIT2: ",nm[i1]," * ",nm[i2]))));}
\q
