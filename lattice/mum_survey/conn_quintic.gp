/* Connection coefficients: quintic mirror, MUM Frobenius basis at z=0
   vs the local solutions at the conifold z=1.
   L = th^4 - z*prod_{j=1..4}(th + j/5),  th = z d/dz.                */
read("cyops.gp");
default(realprecision, 120);
NM = 2000;   \\ MUM series length
NC = 400;    \\ conifold series length

/* ---- MUM Frobenius: R_n(rho) = prod_j (j/5+rho)_n / ((1+rho)_n)^4 ---- */
An = vector(NM+1); r1=vector(NM+1); r2=vector(NM+1); r3=vector(NM+1);
{ my(a=1.0, S1=0.0,S2=0.0,S3=0.0);
  An[1]=1.0; r1[1]=0; r2[1]=0; r3[1]=0;
  for(n=1,NM,
    my(f=1.0);
    for(j=1,4, f *= (j/5 + n-1));
    a *= f/n^4;  An[n+1]=a;
    /* s_k(n) = (-1)^(k+1)/k [ sum_j sum_{i<n} (j/5+i)^-k - 4 sum_{i<=n} i^-k ] */
    my(t1=0.0,t2=0.0,t3=0.0);
    for(j=1,4, my(x=1.0*(j/5+n-1)); t1+=1/x; t2+=1/x^2; t3+=1/x^3);
    S1 += t1 - 4/(1.0*n); S2 += t2 - 4/(1.0*n)^2; S3 += t3 - 4/(1.0*n)^3;
    my(s1=S1, s2=-S2/2, s3=S3/3);
    r1[n+1]=a*s1; r2[n+1]=a*(s2+s1^2/2); r3[n+1]=a*(s3+s1*s2+s1^3/6);
  );
}
Phi(k,z) = my(V=[An,r1,r2,r3][k+1], s=0.0, zz=1.0); \
  for(n=0,NM, s += V[n+1]*zz; zz*=z); s;
ytil(j,z) = my(lg=log(z), s=0.0); for(i=0,j, s += lg^i/i! * Phi(j-i,z)); s;

/* ---- conifold local solutions (exact rationals, then numeric) ---- */
polh = prod(j=1,4,'th+j/5);
CC = vector(5,i,if(i==5,1,0));
for(k=0,4, CC[k+1] -= 't*polcoef(polh,k,'th));
P = thetaToD(CC); Ps = shiftD(P,1); RD = recData(Ps);
mkexp(e,u2) = { my(imp=vector(NC+1,i,0), ini=vector(NC+1,i,0));
  for(i=1,e+1, imp[i]=1); ini[e+1]=1;
  if(e==1 && u2!=0, imp[3]=1; ini[3]=u2);
  solveSeries(RD,NC,ini,imp,[0],0)[1]; }
g0 = mkexp(0,0);
V  = mkexp(1,-7/10);
g2 = mkexp(2,0);
ev(U,s) = my(t=0.0, ss=1.0); for(n=0,NC, t += U[n+1]*ss; ss*=s); t;

/* ---- match at 4 points ---- */
pts = [45/100, 60/100, 75/100, 90/100];
M = matrix(4,4); rhsV = vector(4); rhsg0=vector(4); rhsg2=vector(4);
{ for(a=1,4, my(z=1.0*pts[a], s=z-1);
    for(j=0,3, M[a,j+1] = ytil(j,z));
    rhsV[a]=ev(V,s); rhsg0[a]=ev(g0,s); rhsg2[a]=ev(g2,s));
}
alV  = matsolve(M, rhsV~);
alg0 = matsolve(M, rhsg0~);
alg2 = matsolve(M, rhsg2~);
print("alpha(V)  = ", alV~);
print("alpha(g0) = ", alg0~);
print("alpha(g2) = ", alg2~);
print();
LL = log(3125.0); z2=zeta(2); z3=zeta(3);
print("alpha3*4*Pi^2/sqrt(5) = ", alV[4]*4*Pi^2/sqrt(5));
{ my(b=[alV[3]/alV[4], alV[2]/alV[4], alV[1]/alV[4]]);
  print("beta2 = ",b[1]);  print("   lindep[b2,1,LL] = ",lindep([b[1],1,LL]));
  print("beta1 = ",b[2]);  print("   lindep[b1,1,LL,LL^2,z2] = ",lindep([b[2],1,LL,LL^2,z2]));
  print("beta0 = ",b[3]);  print("   lindep[b0,1,LL,LL^2,LL^3,z2,z2*LL,z3] = ",lindep([b[3],1,LL,LL^2,LL^3,z2,z2*LL,z3]));
}
print();
print("--- g0 in units of alpha3(V) ---");
{ my(c=alg0); print("  ratios: ",[c[1]/c[4],c[2]/c[4],c[3]/c[4]]);
  print("  c4*4*Pi^2/sqrt(5) = ", c[4]*4*Pi^2/sqrt(5)); }
print("--- g2 ---");
{ my(c=alg2); print("  ratios: ",[c[1]/c[4],c[2]/c[4],c[3]/c[4]]);
  print("  c4*4*Pi^2/sqrt(5) = ", c[4]*4*Pi^2/sqrt(5)); }
quit
