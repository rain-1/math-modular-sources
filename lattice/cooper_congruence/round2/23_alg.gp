\\ 23_alg.gp -- is K (resp. Xi, B) mod p algebraic / rational over F_p(q)?
\\ Hermite-Pade style: kernel of the matrix of coefficients of q^i * S^j, checked
\\ against far more coefficients than unknowns.
read("lib.gp");
G  = [read("20_gamma_s7.txt"), read("20_gamma_s10.txt"), read("20_gamma_s18.txt")];
CP = [read("20_cp_s7.txt"),    read("20_cp_s10.txt"),    read("20_cp_s18.txt")];
N  = #G[1];
print("N = ", N);

\\ series (as vector of coefficients of q^0..q^N) for row k, kind 1=K,2=Xi,3=B
{ getser(k,kind) =
  my(v = vector(N+1));
  if(kind==1, for(n=1,N, v[n+1]=G[k][n]));
  if(kind==2, for(n=1,N, v[n+1]=CP[k][n]));
  if(kind==3, for(n=1,N, v[n+1]=n^2*G[k][n]));
  v;
}
KIND = ["K","Xi","B"];

\\ truncated product of two coefficient-vectors mod p
{ mul(a,b,p) = my(L=#a, c=vector(L));
  for(i=1,L, if(a[i]!=0, my(ai=a[i]); for(j=1,L-i+1, c[i+j-1]=(c[i+j-1]+ai*b[j])%p)));
  c;
}

\\ algebraic relation search: deg_q <= D, deg_S <= E
{ algtest(v,p,D,E,verb) =
  my(L=#v, pw=vector(E+1), R, A, ker, rel, ok, W);
  pw[1]=vector(L); pw[1][1]=1;
  for(j=1,E, pw[j+1]=mul(pw[j], vector(L,i,v[i]%p), p));
  R = (D+1)*(E+1) + 20;      \\ rows used to FIND the relation
  if(R>L, R=L);
  A = matrix(R,(D+1)*(E+1));
  for(j=0,E, for(i=0,D,
    my(col=j*(D+1)+i+1);
    for(n=0,R-1, A[n+1,col] = if(n-i>=0, pw[j+1][n-i+1], 0));
  ));
  ker = matker(A*Mod(1,p));
  if(#ker==0, return([0,0]));
  \\ verify the first kernel vector against ALL L coefficients
  rel = lift(ker[,1]);
  W = vector(L);
  for(j=0,E, for(i=0,D,
    my(cc = rel[j*(D+1)+i+1]);
    if(cc!=0, for(n=i,L-1, W[n+1] = (W[n+1] + cc*pw[j+1][n-i+1])%p));
  ));
  ok = 1; for(n=1,L, if(W[n]%p!=0, ok=0; break));
  [#ker, ok];
}

print();
print("=== rationality: relation sum_i c_i q^i * S^j = 0, j<=1 (i.e. S = P/Q), deg <= D ===");
{
for(kind=1,3,
 for(k=1,3,
  for(ip=1,14, my(p=prime(ip), v=getser(k,kind), r);
    r = algtest(v,p,120,1,0);
    if(r[1]>0, print("  ",KIND[kind]," ",NAM[k]," p=",p,"  RATIONAL? kerdim=",r[1]," verified_to_N=",r[2]));
  );
 );
);
}
print("  (no line printed above = no rational relation with deg<=120 for any row/p)");
print();
print("=== algebraicity: deg_q <= 40, deg_S <= 6 ===");
{
for(kind=1,3,
 for(k=1,3,
  for(ip=1,14, my(p=prime(ip), v=getser(k,kind), r);
    r = algtest(v,p,40,6,0);
    if(r[1]>0, print("  ",KIND[kind]," ",NAM[k]," p=",p,"  ALGEBRAIC? kerdim=",r[1]," verified_to_N=",r[2]));
  );
 );
);
}
print("  (no line printed above = no algebraic relation with (D,E)=(40,6) for any row/p)");
quit;
