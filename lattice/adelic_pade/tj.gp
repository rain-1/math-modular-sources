{
g=G2rat(40); G2=g[1]; PREC=g[2];
print("G2 mod 2^",PREC);
for(n=1,4,
  N=3*n; A=vector(N+1); B=vector(N+1); T=vector(N+1);
  for(j=0,N, r=nestgen(3*n,j); A[j+1]=r[1]; B[j+1]=r[2];
     T[j+1]=v2(A[j+1]*G2-B[j+1]));
  print("n=",n);
  for(j=0,N, print("  j=",j," v2A=",v2(A[j+1])," v2B=",v2(B[j+1])," tau=",T[j+1]," tau-v2A=",T[j+1]-v2(A[j+1])));
  /* bridge identity check on all pairs, in cleared-denominator integers */
  bad=0; mn=1000;
  for(j0=0,N, for(j1=j0+1,N,
    d=lcm(lcm(denominator(A[j0+1]),denominator(B[j0+1])),lcm(denominator(A[j1+1]),denominator(B[j1+1])));
    a0=d*A[j0+1];b0=d*B[j0+1];a1=d*A[j1+1];b1=d*B[j1+1];
    h=a0*b1-a1*b0;
    pred=min(v2(a1)+v2(a0*G2-b0), v2(a0)+v2(a1*G2-b1));
    if(v2(h)!=pred, bad++; if(bad<4,print("  MISMATCH ",j0," ",j1," ",v2(h)," vs ",pred)));
    if(v2(h)-2*v2(d)<mn, mn=v2(h)-2*v2(d)) ));
  print("  pairs mismatching ultrametric-equality: ",bad, "   min v2(h)/n normalised(after removing 2v2(d)): ",mn/n);
);
}
\q
