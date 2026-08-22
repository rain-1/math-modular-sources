
\\ 20_family.gp -- the family (n+1)^2 u_{n+1} = 4(2n+1) u_n - c n^2 u_{n-1}.
\\ (a) which c make the A-row integral?  (b) deep check of the members c = -4^m:
\\ integrality, k, kappa_2, sigma_2, archimedean rate, nonvanishing, rationality,
\\ algebraicity of xi_2.
\\ Run: gp -q /home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/20_family.gp
default(parisizemax,"12G");
default(realprecision,40);

arow(c,N) = { my(u=vector(N+2)); u[1]=0; u[2]=1; for(n=0,N-1, u[n+3]=((8*n+4)*u[n+2]-c*n^2*u[n+1])/(n+1)^2); vector(N+1,j,u[j+1]); }
brow(c,N) = { my(u=vector(N+2)); u[1]=0; u[2]=0; u[3]=1; for(n=1,N-1, u[n+3]=((8*n+4)*u[n+2]-c*n^2*u[n+1])/(n+1)^2); vector(N+1,j,u[j+1]); }
okint(c,N) = { my(u0=0,u1=1,t,d); for(n=0,N-1, t=(8*n+4)*u1-c*n^2*u0; d=(n+1)^2; if(t%d, return(0)); u0=u1; u1=t/d); 1; }

print("(a) c with |c| <= 30000 for which the A-row is integral to n = 60:");
{
my(L=List());
for(c=-30000,30000, if(c!=0 && okint(c,60), listput(L,c)));
print("    ", Vec(L));
print("    (these are exactly c = -4^m, m = 2..7)");
}

print("");
print("(b) deep check of c = -4^m, exact to n = 400:");
{
my(N=400);
for(m=2,9,
  my(c=-4^m, aa=arow(c,N), bb=brow(c,N));
  my(ok=1); for(n=0,N, if(denominator(aa[n+1])!=1, ok=0; break));
  my(kmax=0, d=1, v2b=0);
  for(n=1,N, d=lcm(d,n);
    my(x=bb[n+1], od=d/2^valuation(d,2), kk=0);
    if(x!=0, if(valuation(x,2)<v2b, v2b=valuation(x,2)));
    while(denominator(x*od^kk)!=1 && kk<20, kk++); if(kk>kmax,kmax=kk));
  my(xi=bb[N+1]/aa[N+1]);
  my(v1=valuation(xi-bb[201]/aa[201],2), v2=valuation(xi-bb[301]/aa[301],2));
  my(sg=(v2-v1)/100.);
  my(kap=-(valuation(aa[N+1],2)-valuation(aa[201],2))/200.);
  my(l1=2^m, arch=log(abs(aa[N+1]*1.))/N);
  my(nonvan = !(bb[201]/aa[201]==bb[301]/aa[301] && bb[301]/aa[301]==bb[401]/aa[401]));
  printf("m=%d c=%9d: A int=%d  k=%d  min v_2(B_n)=%d  kappa_2=%6.3f  sigma_2=%7.3f  arch=%8.5f (log l1=%8.5f)  nonvanishing=%d  S=%+9.6f  theta=%9.6f\n",
     m, c, ok, kmax, v2b, kap, sg, arch, log(l1*1.), nonvan,
     2*m*log(2)-2-log(l1*1.), 2*m*log(2)/(2+log(l1*1.)));
);
}

print("");
print("(c) rationality / algebraicity of xi_2, m = 3..7 (2-adic lindep):");
{
my(N=400);
for(m=3,7,
  my(c=-4^m, aa=arow(c,N), bb=brow(c,N), xi=bb[N+1]/aa[N+1]);
  my(PR=2*m*300-60, x=xi+O(2^PR));
  my(l1=lindep([1+O(2^PR), x]));
  my(l4=lindep(vector(5,j,x^(j-1))));
  printf("m=%d: max|coeff| of lindep[1,xi] = %d digits ; of lindep[1,xi,..,xi^4] = %d digits\n",
     m, #Str(vecmax(abs(Vec(l1)))), #Str(vecmax(abs(Vec(l4)))));
);
print("   (hundreds of digits = NO relation: xi_2 is neither rational nor algebraic of degree <= 4)");
}
quit;
