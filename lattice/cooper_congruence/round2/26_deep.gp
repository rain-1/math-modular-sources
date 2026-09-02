\\ 26_deep.gp -- deeper mod-p structure: high-degree rationality, larger algebraic
\\ boxes, the Mahler/Cartier relation K(q)=A(q)K(q^p)+B(q), and the p-kernel.
read("lib.gp");
G  = [read("22_gamma_s7.txt"), read("22_gamma_s10.txt"), read("22_gamma_s18.txt")];
CP = [read("22_cp_s7.txt"),    read("22_cp_s10.txt"),    read("22_cp_s18.txt")];
N  = #G[1];
print("N = ", N);
KIND = ["K","Xi","B"];
{ getser(k,kind) = my(v=vector(N+1));
  if(kind==1, for(n=1,N, v[n+1]=G[k][n]));
  if(kind==2, for(n=1,N, v[n+1]=CP[k][n]));
  if(kind==3, for(n=1,N, v[n+1]=n^2*G[k][n]));
  v; }
{ mul(a,b,p) = my(L=#a, c=vector(L));
  for(i=1,L, if(a[i]%p!=0, my(ai=a[i]%p); for(j=1,L-i+1, c[i+j-1]=(c[i+j-1]+ai*b[j])%p)));
  c; }

\\ generic Hermite-Pade: columns are q^i*S^j, i<=D, j<=E; find kernel from R rows,
\\ verify against all N+1 rows.  returns [kerdim, verified?]
{ hp(v,p,D,E) =
  my(L=#v, pw=vector(E+1), R, A, ker, rel, W, ok);
  pw[1]=vector(L); pw[1][1]=1;
  for(j=1,E, pw[j+1]=mul(pw[j], vector(L,i,v[i]%p), p));
  R=(D+1)*(E+1)+15; if(R>L,R=L);
  A=matrix(R,(D+1)*(E+1));
  for(j=0,E, for(i=0,D, my(col=j*(D+1)+i+1);
    for(n=0,R-1, A[n+1,col]=if(n-i>=0, pw[j+1][n-i+1], 0))));
  ker=matker(A*Mod(1,p));
  if(#ker==0, return([0,0]));
  rel=lift(ker[,1]); W=vector(L);
  for(j=0,E, for(i=0,D, my(cc=rel[j*(D+1)+i+1]);
    if(cc!=0, for(n=i,L-1, W[n+1]=(W[n+1]+cc*pw[j+1][n-i+1])%p))));
  ok=1; for(n=1,L, if(W[n]%p!=0, ok=0; break));
  [#ker,ok]; }

print();
print("=== rationality, deg <= 250, p = 3..43, N = ", N, " ===");
{ for(kind=1,3, for(k=1,3, for(ip=2,14, my(p=prime(ip), r=hp(getser(k,kind),p,250,1));
    if(r[1]>0, print("  ",KIND[kind]," ",NAM[k]," p=",p,"  kerdim=",r[1]," verified=",r[2]))))); }
print("  (only s18/p=3 expected)");

print();
print("=== algebraic boxes (D,E) = (20,20), (60,4), (10,30) for K, p=3,5,7,11,13 ===");
{ for(bx=1,3, my(BX=[[20,20],[60,4],[10,30]][bx], D=BX[1], E=BX[2]);
   for(k=1,3, for(ip=2,6, my(p=prime(ip), r=hp(getser(k,1),p,D,E));
     if(r[1]>0, print("  K ",NAM[k]," p=",p," (D,E)=(",D,",",E,")  kerdim=",r[1]," verified=",r[2]))))); }
print("  (only s18/p=3 expected)");

print();
print("=== Mahler/Cartier: K(q) = A(q) K(q^p) + B(q) mod p, deg A,B <= 60 ===");
{
for(k=1,3, for(ip=1,6, my(p=prime(ip), v=getser(k,1), L=N+1, Kp, A, ker, rel, D=60, R);
  Kp = vector(L); for(n=1,(L-1)\p, Kp[p*n+1] = v[n+1]%p);   \\ K(q^p)
  R = 2*(D+1)+15; if(R>L,R=L);
  A = matrix(R, 2*(D+1)+0);
  for(i=0,D, for(n=0,R-1,
    A[n+1,i+1]     = if(n-i>=0, Kp[n-i+1], 0);              \\ q^i K(q^p)
    A[n+1,D+2+i]   = if(n==i, 1, 0)));                      \\ q^i
  \\ target: K(q); solve A*z = K  ->  kernel of [A | -K]
  my(AA = matconcat([A, matrix(R,1,x,y,-v[x])]));
  ker = matker(AA*Mod(1,p));
  \\ need a kernel vector with last coord nonzero
  my(good=0);
  for(j=1,#ker, if(lift(ker[2*(D+1)+1,j])%p!=0, good=j; break));
  if(good==0, print("  ",NAM[k]," p=",p,": NO relation with deg<=",D),
    rel = lift(ker[,good]); my(c0=rel[2*(D+1)+1], W=vector(L), ok=1);
    for(i=0,D, my(ca=rel[i+1], cb=rel[D+2+i]);
      if(ca!=0, for(n=i,L-1, W[n+1]=(W[n+1]+ca*Kp[n-i+1])%p));
      if(cb!=0, W[i+1]=(W[i+1]+cb)%p));
    for(n=1,L, W[n]=(W[n]-c0*v[n])%p);
    for(n=1,L, if(W[n]%p!=0, ok=0; break));
    print("  ",NAM[k]," p=",p,": relation found, kerdim=",#ker,", verified to q^",N,"? ",ok)));
}

print();
print("=== p-kernel: how many distinct subsequences n -> gamma(p^s n + r) mod p ? ===");
{
for(k=1,3, for(ip=1,6, my(p=prime(ip), S=List(), cnt=vector(3));
  for(s=0,2, my(ps=p^s, LL=N\ps);
    if(LL<80, break);
    for(r=0,ps-1, my(w=vector(80,i,lift(Mod(G[k][ps*(i-1)+r+ if(r==0,ps,0)],p))));
      \\ subsequence n=1.. : index ps*n+r ; skip if index 0
      w = vector(80, i, my(idx=ps*i+r); if(idx>=1 && idx<=N, lift(Mod(G[k][idx],p)), -1));
      listput(S,w));
    cnt[s+1]=#Set(Vec(S)));
  print("  ",NAM[k]," p=",p,"  cumulative #distinct kernel elements at depth 0,1,2: ",cnt)));
}
quit;
