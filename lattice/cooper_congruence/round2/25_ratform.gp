\\ 25_ratform.gp -- explicit minimal rational form of K mod p where it exists,
\\ periodic patterns, and periodicity mod p^k.
read("lib.gp");
G3 = [read("22_gamma_s7.txt"), read("22_gamma_s10.txt"), read("22_gamma_s18.txt")];
N  = #G3[1];
print("N = ", N);

\\ minimal (P,Q) with Q*K = P mod (p, q^L), deg <= d ; returns [d,P,Q] or 0
{ minrat(v,p,L,dmax) =
  my(A, ker, rel, P, Q, KK);
  KK = vector(L, i, v[i]%p);          \\ coefficients of q^0..q^(L-1)
  for(d=1,dmax,
    my(R = 2*d+2+10); if(R>L, R=L);
    A = matrix(R, 2*d+2);
    for(i=0,d, for(n=0,R-1,
      A[n+1,i+1]      = if(n-i>=0, KK[n-i+1], 0);     \\ q^i * K
      A[n+1,d+2+i]    = if(n==i, -1, 0)));            \\ -q^i
    ker = matker(A*Mod(1,p));
    if(#ker>0,
      rel = lift(ker[,1]);
      Q = sum(i=0,d, rel[i+1]*'q^i);
      P = sum(i=0,d, rel[d+2+i]*'q^i);
      return([d,P,Q]));
  );
  0;
}

\\ verify Q*K - P = 0 mod (p, q^L)
{ ver(v,p,L,P,Q) =
  my(KK = vector(L,i,v[i]%p), CQ=Vec(Q,-(poldegree(Q)+1)), CP=Vec(P,-(poldegree(P)+1)), W=vector(L));
  for(i=1,#CQ, if(CQ[i]!=0, for(n=1,L-i+1, W[n+i-1]=(W[n+i-1]+CQ[i]*KK[n])%p)));
  for(i=1,#CP, W[i]=(W[i]-CP[i])%p);
  for(n=1,L, if(W[n]%p!=0, return(n-1)));
  -1;
}

print();
print("=== minimal rational form of K mod p ===");
{
for(k=1,3, for(ip=1,3, my(p=prime(ip), v, r);
  v = vector(N+1); for(n=1,N, v[n+1]=G3[k][n]);
  r = minrat(v,p,800,60);
  if(r==0, print("  ",NAM[k]," p=",p,": none with deg<=60 (found from 800 terms)"),
    print("  ",NAM[k]," p=",p,":  deg=",r[1]);
    print("      P = ", lift(r[2]));
    print("      Q = ", lift(r[3]));
    print("      verified to q^",N,"? ", if(ver(v,p,N+1,r[2],r[3])==-1,"YES",concat("NO, first failure n=",ver(v,p,N+1,r[2],r[3])))));
));
}

print();
print("=== periodic pattern of gamma(n) mod p^k (pure period; n=1..) ===");
{
for(k=1,3, for(ip=1,3, my(p=prime(ip));
  for(e=1,6, my(m=p^e, T=0);
    for(t=1,3000, if(t*2>N, break);
      my(ok=1); for(n=1,N-t, if((G3[k][n]-G3[k][n+t])%m!=0, ok=0; break));
      if(ok, T=t; break));
    if(T>0, print("  ",NAM[k]," mod ",p,"^",e,"=",m,": period ",T,
                  if(e==1, concat("   pattern gamma(1..T) = ", vector(T,i,lift(Mod(G3[k][i],m)))), "")),
            print("  ",NAM[k]," mod ",p,"^",e,"=",m,": no period <= ",N\2); break));
));
}
quit;
