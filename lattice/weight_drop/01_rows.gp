default(realprecision,150);
NMAX = 400;
{par(ab, N) = my(a=ab[1],b=ab[2],c=ab[3],d=ab[4], A=vector(N+1), B=vector(N+1));
  A[1]=1; A[2]=b; B[1]=0; B[2]=1;
  for(n=1,N-1,
    A[n+2] = ((2*n+1)*(a*n^2+a*n+b)*A[n+1] + n*(c*n^2+d)*A[n])/(n+1)^3;
    B[n+2] = ((2*n+1)*(a*n^2+a*n+b)*B[n+1] + n*(c*n^2+d)*B[n])/(n+1)^3;
  ); [A,B];}
{root(p1,p0,N) = my(A=vector(N+1),B=vector(N+1));
  A[1]=1; A[2]=subst(p1,'n,0); B[1]=0; B[2]=1;
  for(n=1,N-1,
    A[n+2] = (subst(p1,'n,n)*A[n+1] - subst(p0,'n,n)*A[n])/(n+1)^2;
    B[n+2] = (subst(p1,'n,n)*B[n+1] - subst(p0,'n,n)*B[n])/(n+1)^2;
  ); [A,B];}
cauchy(u,v,N) = vector(N+1, k, sum(i=0,k-1, u[i+1]*v[k-i]));
rows = [["s7", [13,4,27,-3], 26*'n^2+13*'n+2, -3*(3*'n-1)*(3*'n-2), 1], ["s10", [6,2,64,-4], 24*'n^2+12*'n+2, -4*(8*'n-3)*(8*'n-5), 2], ["s18", [14,6,-192,12], 56*'n^2+28*'n+6, 12*(8*'n-3)*(8*'n-5), 2]];
{for(r=1,#rows,
  my(nm=rows[r][1], PB=par(rows[r][2],NMAX), RT=root(rows[r][3],rows[r][4],NMAX), lam=rows[r][5]);
  my(A=PB[1],B=PB[2],a=RT[1],b=RT[2]);
  print("=== ", nm, "   lambda=", lam);
  print("  A = ", vector(6,k,A[k]), "   a = ", vector(6,k,a[k]));
  print("  B = ", vector(6,k,B[k]), "   b = ", vector(6,k,b[k]));
  my(aa=cauchy(a,a,NMAX), ab=cauchy(a,b,NMAX));
  print("  max|a*a - lam^n A_n| = ", vecmax(apply(abs,vector(NMAX+1,k, aa[k] - lam^(k-1)*A[k]))));
  my(dl = vector(NMAX+1,k, lam^(k-1)*B[k] - ab[k]));
  print("  delta_n := lam^n B_n - (a*b)_n = ", vector(8,k,dl[k]));
  print("  xi_par  = ", B[NMAX+1]/A[NMAX+1]*1.0);
  print("  xi_root = ", b[NMAX+1]/a[NMAX+1]*1.0);
  print("  ratio   = ", (B[NMAX+1]/A[NMAX+1])/(b[NMAX+1]/a[NMAX+1])*1.0);
  print("  delta_n/(a*a)_n -> ", dl[NMAX+1]/aa[NMAX+1]*1.0);
);}
print("targets: zeta(2)/7=",zeta(2)/7," zeta(2)/5=",zeta(2)/5);
quit;
