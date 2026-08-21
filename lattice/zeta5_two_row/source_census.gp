/* Census of trivial-character weight-6 zeta(5) Eisenstein-oldform sources that satisfy
   the p-adic Euler-factor criterion.  See header of the earlier version.        */
LOG = "lattice/zeta5_two_row/source_census.log";
W(s) = write(LOG, s);
NMAX = 200;
g0(N) = my(f=factor(N), idx=N*prod(i=1,#f[,1],1+1/f[i,1]), \
   nu2=if(N%4==0,0,prod(i=1,#f[,1],1+kronecker(-1,f[i,1]))), \
   nu3=if(N%9==0,0,prod(i=1,#f[,1],1+kronecker(-3,f[i,1]))), \
   ninf=sumdiv(N,d,eulerphi(gcd(d,N/d)))); \
   if(N==1, 0, 1 + idx/12 - nu2/4 - nu3/3 - ninf/2);
prim(v) = my(n=#v, den=1, g=0, w); for(j=1,n, den=lcm(den,denominator(v[j]))); \
   w = vector(n, j, v[j]*den); for(j=1,n, g=gcd(g,w[j])); if(g>0, w=vector(n,j,w[j]/g)); w;

W(Str("=== zeta(5) weight-6 trivial-character source census, N <= ", NMAX, " ==="));
W("purification: P(0)=P(2)=P(4)=P(6)=0 and P(5) != 0.");
W("slope at p (Euler-factor criterion): sum_a c_{p^a m} = 0 for every prime-to-p class m.");
W("");
{
for(N = 1, NMAX,
  D = divisors(N); nd = #D;
  if(nd < 5, next);
  Mp = matrix(4, nd, i, j, (D[j])^(-2*(i-1)));
  K = matker(Mp);
  if(#K == 0, next);
  v5 = vector(nd, j, D[j]^(-5));
  if(vecmax(vector(#K, c, if(v5*K[,c]!=0,1,0))) == 0, next);
  fp = factor(N)[,1];
  for(ip = 1, #fp,
    p = fp[ip];
    ms = Set(vector(nd, j, D[j]/p^valuation(D[j],p)));
    Me = matrix(#ms, nd, i, j, if(D[j]/p^valuation(D[j],p) == ms[i], 1, 0));
    K2 = matker(matconcat([Mp; Me]));
    if(#K2 == 0, next);
    good = [];
    for(c = 1, #K2, if(v5*K2[,c] != 0, good = concat(good,[c])));
    if(#good == 0, next);
    W(Str("N=", N, "  genus X_0(N)=", g0(N), "  p=", p, "  dim(purified)=", #K,
          "  dim(purified & slope)=", #K2, "  divisors=", D));
    for(ic = 1, #K2,
      vv = prim(Vec(K2[,ic]));
      W(Str("    c_d = ", vv, "    P(5) = ", v5*vv~, "    L(Phi,5) = ", -(v5*vv~)/2, " * zeta(5)")))));
}
W("");
W("DONE"); quit;
