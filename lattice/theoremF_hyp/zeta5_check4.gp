/* ===========================================================================
   zeta5_check4.gp -- exact Picard-Fuchs operator for the level-16 zeta(5)
   A-row (order 5, degree 38) and the 2-adic location of its singular points;
   plus the archimedean/2-adic singularity data for the level-12 h_12 host,
   whose dominant singularity is t = 7-4*sqrt(3) (a 2-adic UNIT).
   =========================================================================== */
default(parisizemax, 24000000000);
\p 60
if(type(MM)!="t_INT", MM = 420);
M = MM;
q = 'q; ee(d) = eta(q^d + O(q^(M+2)));

mkrow(cv, dv, T) =
{ my(g, PHI, DT, F, QT, A, s);
  g = vector(M);
  for(n=1,M, s=0; for(i=1,#dv, if(n%dv[i]==0, s += cv[i]*sigma(n/dv[i],5))); g[n]=s);
  PHI = sum(n=1,M, g[n]*q^n) + O(q^(M+1));
  DT = q*deriv(T,q); F = PHI/DT; QT = serreverse(T);
  A = subst(F,q,QT);
  vector(M, i, polcoeff(A, i-1));
}
pfrow(V, R, D, n) =
{ my(L=List(), m, co);
  for(k=0,R, for(j=0,D, m = n-j+k;
    if(m<0 || m>#V-1, co=0, co = prod(i=0,k-1, n-j+k-i)*V[m+1]); listput(L,co)));
  Vec(L);
}
pfexact(name, V, R, D) =
{ my(NN=#V, nv, nhi, K, VV, P, den, gg, P5, ff, vmin, rts);
  nv=(R+1)*(D+1); nhi=min(NN-1, nv+80);
  K = matker(matconcat(vector(nhi+1, i, pfrow(V,R,D,i-1))~));
  print("--- ", name, ": exact order-",R," degree-",D," fit, kernel dim = ", #K);
  if(#K==0, return(0));
  VV = K[,1];
  P = vector(R+1, k, sum(j=0,D, VV[(k-1)*(D+1)+j+1]*'t^j));
  den=1; for(k=1,R+1, den=lcm(den, denominator(content(P[k]))));
  P = vector(R+1, k, P[k]*den);
  gg=0; for(k=1,R+1, gg=gcd(gg, content(P[k])));
  if(gg>0, P = vector(R+1, k, P[k]/gg));
  P5 = P[R+1];
  print("    deg P_k, k=0..",R," : ", vector(R+1,k,poldegree(P[k])));
  print("    LEADING coefficient P_",R,"(t) factored over Q:");
  print("      ", factor(P5));
  print("    Wronskian/Casoratian:  W'/W = -P_",R-1,"/P_",R,";  P_",R-1," factored:");
  print("      ", factor(P[R]));
  print("    2-adic: t-multiplicity at 0 = ", valuation(P5,'t));
  ff = factorpadic(P5/'t^valuation(P5,'t), 2, 30);
  vmax = -10^9;
  for(i=1, matsize(ff)[1],
    my(f=ff[i,1], e=ff[i,2], d=poldegree(f),
       v0=valuation(polcoeff(f,0),2), vd=valuation(polcoeff(f,d),2), vv=(v0-vd)/d);
    print("      irreducible /Q_2 factor of degree ", d, " (mult ", e, "):  v_2(its roots) = ", vv,
          "   |t|_2 = 2^", -vv);
    if(vv>vmax, vmax=vv));
  /* rho = min_s |s|_2 = 2^(-max_s v_2(s)) ; sigma_2 = log_2 rho = -max_s v_2(s) */
  print("    ==> max over singular t != 0 of v_2(t) = ", vmax, "  (the 2-adically NEAREST singularity)");
  print("    ==> largest singularity-free disc {|t|_2 <= rho}:  rho = 2^", -vmax);
  print("    ==> PREDICTED SLOPE  sigma_2 = log_2 rho = ", -vmax);
  rts = polroots(P5);
  print("    |complex roots| sorted: ", vecsort(vector(#rts,i,abs(rts[i]))));
  P;
}
XX  = q*ee(2)*ee(16)^2/(ee(1)^2*ee(8)) + O(q^(M+1)); T16 = XX/(8*XX^2+2*XX+1);
c16 = [1,-85,1428,-5440,4096]; d16 = [1,2,4,8,16];
A16 = mkrow(c16,d16,T16);
pfexact("level-16", A16, 5, 38);
print("");
print("--- level-12 h_12: archimedean singularity check");
HH  = (ee(1)^3*ee(4)*ee(6)^2/(ee(2)^2*ee(3)*ee(12)^3))/q + O(q^(M+1));
TX  = HH/((HH+3)*(HH+4));
c12 = [1,-113,567,112,-1863,1296]; d12 = [1,2,3,4,6,12];
A12X = mkrow(c12,d12,TX);
print("   a_n/a_{n+1} at n=", M-2, " : ", A12X[M-1]*1.0/A12X[M], "     7-4*sqrt(3) = ", 7-4*sqrt(3));
print("   log|a_n|/n at n=", M-1, " : ", log(abs(A12X[M]*1.0))/(M-1), "   log(7+4sqrt3) = ", log(7+4*sqrt(3)));
print("   lambda^2-14 lambda+1 over Q_2: ", factorpadic('x^2-14*'x+1, 2, 20));
print("   v_2 of its roots: Newton polygon of x^2-14x+1 : v(1)=0, v(-14)=1, v(1)=0  => both roots are 2-adic UNITS");
print("   so the singular point t = 7-4sqrt3 = 1/(7+4sqrt3) has |t|_2 = 1  (ON the unit circle)");
print("DONE"); quit;
