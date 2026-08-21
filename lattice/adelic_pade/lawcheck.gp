/* Conjectured closed forms, m=3n, M=m+j:
   v2(A_{n,j}) = 3 - 2*(m+M) + s2(m) + s2(M)
   v2(B_{n,j}) = v2(A)-1
   tau_j := v2(A G2 - B) = 5 - v2(A)
   => v2(G2 - B/A) = 5 - 2 v2(A) = 4*(m+M) - 1 - 2 s2(m) - 2 s2(M) = 24n+4j-1-2s2(3n)-2s2(3n+j) */
s2(x)=hammingweight(x);
{
g=G2rat(46); G2=g[1];
print("G2 mod 2^",g[2]);
bad=0;
for(n=1,5, m=3*n;
  for(j=0,3*n, MM=m+j; r=nestgen(m,j); A=r[1]; B=r[2];
    pa = 3-2*(m+MM)+s2(m)+s2(MM);
    if(v2(A)!=pa, bad++; print("A FAIL n=",n," j=",j," ",v2(A)," vs ",pa));
    if(v2(B)!=pa-1, bad++; print("B FAIL n=",n," j=",j));
    if(v2(A*G2-B)!=5-pa, bad++; print("TAU FAIL n=",n," j=",j," ",v2(A*G2-B)," vs ",5-pa));
  ); print("n=",n," done, quality at j=0,n,3n: ",[4*(m+m)-1-4*s2(m), 4*(m+m+n)-1-2*s2(m)-2*s2(m+n), 4*(m+m+3*n)-1-2*s2(m)-2*s2(2*m)]);
);
print("total failures: ",bad);
}
\q
