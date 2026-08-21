\\ Empirical measurement of the two-lattice construction quality by LLL.
\\ No asymptotic theory: build the actual congruence lattice, reduce, measure.
read("/home/ubuntu/code/math-modular-sources/lattice/zeta3_lattice/rows.gp");
default(realprecision,2000);
allocatemem(2*10^9);
XAN=6.7069026; SIG=17.3177662;  \\ theory anisotropy and covolume rate
D=rowsI(10,4,64); T=rowsI(12,4,16);
bb(R,n)=R[2][n+1]/(n!)^3;          \\ b_n as exact rational
z3=zeta(3);
dl=vector(400); {my(d=1);for(i=1,400,d=lcm(d,i);dl[i]=d);}
run(n)={
  my(m=2*n,k=3*n, S=dl[k]^3);
  my(X1=S*7*D[1][m+1], Y1=S*24*bb(D,m), X2=S*7*T[1][k+1], Y2=S*32*bb(T,k));
  if(denominator(Y1)!=1||denominator(Y2)!=1, error("non-integral"));
  my(a1=7*D[1][m+1], a2=7*T[1][k+1]);
  my(Del=3*T[1][k+1]*bb(D,m)-4*D[1][m+1]*bb(T,k));
  my(v=valuation(56*Del,2), Tn=2^v, M=S*Tn);
  \\ lattice K = {c in Z^2 : a.c = 0 mod Tn, Y.c = 0 mod M}
  my(sol=matsolvemod([a1,a2;Y1,Y2],[Tn,M]~,[0,0]~,1));
  my(B=sol[2]);  \\ columns = basis of K
  \\ embed: for c in K, q=(c.X)/M, e=(c.X)*z3-(c.Y) all over M
  \\ Measure the CONSTRUCTION, not the best lattice vector:
  \\ reduce K in the anisotropic box |c1|<=e^{x n}, |c2|<=e^{(sig-x)n}
  my(W=round(exp((SIG-2*XAN)*n)));   \\ relative weight on c1
  my(Bw=[W*B[1,1],W*B[1,2]; B[2,1],B[2,2]]);
  my(U=qflll(Bw), best=[0,0,-oo]);
  for(j=1,2,for(i=1,1,
    my(c=B*U[,j]);
    my(q=(c[1]*X1+c[2]*X2)/M, p=(c[1]*Y1+c[2]*Y2)/M);
    if(q!=0, my(e=abs(q*z3-p));
      if(e!=0 && (best[3]==-oo || log(abs(q*1.))>best[1]),
        best=[log(abs(q*1.)),log(e),1-log(e)/log(abs(q*1.))]))));
  [best[1]/n, best[2]/n, best[3], v-12*n, log(M*1.)/n];
};
print("n |  logq/n (H_emp) | log|err|/n (F_emp) | delta_emp | v2 deficit | log M/n");
print("theory (n->oo):      11.7392            1.1627         0.90095                    17.3178");
{for(j=2,20, my(n=5*j, r=run(n));
  printf("%3d |  %10.4f | %10.4f | %10.6f | %5d | %10.4f\n", n, r[1],r[2],r[3],r[4],r[5]));}
\q
