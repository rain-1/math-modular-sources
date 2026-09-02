\\ 26_mahler.gp -- (3e) K(q) = A(q) K(q^p) + B(q) mod p ;  and the p-kernel of gamma mod p.
read("lib.gp");
G  = [read("20_gamma_s7.txt"), read("20_gamma_s10.txt"), read("20_gamma_s18.txt")];
N  = #G[1];
print("N = ", N);
print();
print("=== K(q) = A(q) K(q^p) + B(q) mod p,  deg A, deg B <= 60 ===");
{
my(D=60);
for(k=1,3,
 for(ip=1,6,
  my(p=prime(ip), L=N+1, v, Kp, A, ker, good, rel, W, ok, R, c0);
  v = vector(L); for(n=1,N, v[n+1]=G[k][n]%p);
  Kp = vector(L); for(n=1,(L-1)\p, Kp[p*n+1] = v[n+1]);
  R = 2*(D+1)+15; if(R>L, R=L);
  A = matrix(R, 2*(D+1)+1);
  for(i=0,D,
    for(n=0,R-1,
      A[n+1,i+1]   = if(n-i>=0, Kp[n-i+1], 0);
      A[n+1,D+2+i] = if(n==i, 1, 0)));
  for(n=0,R-1, A[n+1,2*(D+1)+1] = -v[n+1]);
  ker = matker(A*Mod(1,p));
  good = 0;
  for(j=1,#ker, if(lift(ker[2*(D+1)+1,j])%p!=0, good=j; break));
  if(good==0,
    print("  ",NAM[k]," p=",p,":  NO relation with deg A,B <= ",D)
  ,
    rel = lift(ker[,good]);
    W = vector(L);
    for(i=0,D,
      my(ca=rel[i+1], cb=rel[D+2+i]);
      if(ca!=0, for(n=i,L-1, W[n+1]=(W[n+1]+ca*Kp[n-i+1])%p));
      if(cb!=0, W[i+1]=(W[i+1]+cb)%p));
    c0=rel[2*(D+1)+1];
    for(n=1,L, W[n]=(W[n]-c0*v[n])%p);
    ok=1; for(n=1,L, if(W[n]%p!=0, ok=0; break));
    print("  ",NAM[k]," p=",p,":  relation found (kerdim ",#ker,"), verified to q^",N,"? ",ok);
  );
 );
);
}
print();
print("=== p-kernel: # distinct subsequences n -> gamma(p^s n + r) mod p, s=0,1,2 ===");
{
for(k=1,3,
 for(ip=1,6,
  my(p=prime(ip), S=List(), cnt=vector(3), W);
  for(s=0,2,
    my(ps=p^s);
    if(ps*100+ps > N, cnt[s+1]=-1; next);
    for(r=0,ps-1,
      W = vector(100, i, my(idx=ps*i+r); if(idx>=1 && idx<=N, lift(Mod(G[k][idx],p)), -1));
      listput(S,W));
    cnt[s+1]=#Set(Vec(S));
  );
  print("  ",NAM[k]," p=",p,":  cumulative distinct kernel elements, depth 0,1,2 = ",cnt,
        "   (total subsequences 1+p+p^2 = ",1+p+p^2,")");
 );
);
}
quit;
