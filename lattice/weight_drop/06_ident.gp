default(parisizemax,4000000000);
default(realprecision,120);
NMAX=900;
{par(ab,N)=my(a=ab[1],b=ab[2],c=ab[3],d=ab[4],A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]+n*(c*n^2+d)*A[n])/(n+1)^3;
             B[n+2]=((2*n+1)*(a*n^2+a*n+b)*B[n+1]+n*(c*n^2+d)*B[n])/(n+1)^3;);[A,B];}
{root(p1,p0,N)=my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=subst(p1,'n,0);B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=(subst(p1,'n,n)*A[n+1]-subst(p0,'n,n)*A[n])/(n+1)^2;
             B[n+2]=(subst(p1,'n,n)*B[n+1]-subst(p0,'n,n)*B[n])/(n+1)^2;);[A,B];}
rows=[["s7",[13,4,27,-3],26*'n^2+13*'n+2,-3*(3*'n-1)*(3*'n-2)],["s10",[6,2,64,-4],24*'n^2+12*'n+2,-4*(8*'n-3)*(8*'n-5)],["s18",[14,6,-192,12],56*'n^2+28*'n+6,12*(8*'n-3)*(8*'n-5)]];
L2m3 = lfun(-3,2);
targ=[zeta(2)/7, zeta(2)/5, L2m3/2];
{for(r=1,3, my(PB=par(rows[r][2],NMAX),RT=root(rows[r][3],rows[r][4],NMAX));
  my(xp=PB[2][NMAX+1]/PB[1][NMAX+1]*1.0, xr=RT[2][NMAX+1]/RT[1][NMAX+1]*1.0);
  print("=== ",rows[r][1]);
  print("  xi_par - target = ", xp - targ[r]);
  print("  xi_root = ",xr);
  print("  ratio xi_par/xi_root = ",xp/xr);
  print("  algdep(ratio,4) = ",algdep(xp/xr,4),"  algdep(...,6)=",algdep(xp/xr,6));
  print("  lindep([1,xp,xr]) = ",lindep([1,xp,xr]));
  print("  lindep([xp,xr]) = ",lindep([xp,xr]));
);}
print();
print("--- root periods vs candidate constants ---");
{my(RT=root(rows[3][3],rows[3][4],NMAX), xr=RT[2][NMAX+1]/RT[1][NMAX+1]*1.0);
 print("  s18 root xi = ",xr);
 print("  lindep([xr,1,zeta(2),L2m3,zeta(3),Pi^2,Pi^3]) = ",lindep([xr,1,zeta(2),L2m3,zeta(3),Pi^2,Pi^3]));
 print("  algdep(xr,8) = ",algdep(xr,8));}
quit;
