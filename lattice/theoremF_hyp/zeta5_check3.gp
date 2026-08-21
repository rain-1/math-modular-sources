/* ===========================================================================
   zeta5_check3.gp -- the HOST side.  Fit the order-5 Picard-Fuchs ODE
        sum_{k=0}^{5} P_k(t) A^{(k)}(t) = 0
   satisfied by the A-row generating function A(t)=sum a_n t^n of each
   zeta(5) host, then read off:
     * the singular points  = roots of P_5(t)  (t=0 excluded: it is the cusp),
     * their 2-adic valuations, hence the radius rho = min_{s != 0} |s|_2
       of the largest disc {|t|_2 <= rho'} free of singularities,
     * the predicted slope  sigma_2 = log_2 rho,
     * the Wronskian (Casoratian) exponent data  W'/W = -P_4/P_5.
   Set MM before reading (default 420).
   =========================================================================== */
default(parisizemax, 24000000000);
default(threadsizemax, 8000000000);
\p 60
if(type(MM)!="t_INT", MM = 420);
M = MM;
q = 'q; ee(d) = eta(q^d + O(q^(M+2)));
PP = 2305843009213693951;

mkrow(cv, dv, T) =
{ my(g, PHI, DT, F, QT, A, s);
  g = vector(M);
  for(n=1,M, s=0; for(i=1,#dv, if(n%dv[i]==0, s += cv[i]*sigma(n/dv[i],5))); g[n]=s);
  PHI = sum(n=1,M, g[n]*q^n) + O(q^(M+1));
  DT = q*deriv(T,q); F = PHI/DT; QT = serreverse(T);
  A = subst(F,q,QT);
  vector(M, i, polcoeff(A, i-1));
}

/* row n of the linear system, unknowns c_{k,j}, k=0..R, j=0..D            */
pfrow(V, R, D, n) =
{ my(L=List(), idx, m, co);
  for(k=0,R, for(j=0,D,
    m = n - j + k;
    if(m < 0 || m > #V-1, co = 0,
       co = prod(i=0,k-1, n-j+k-i) * V[m+1]);
    listput(L, co)));
  Vec(L);
}

pfscan(name, V, R) =
{ my(NN=#V, Vm, s);
  Vm = vector(NN, i, Mod(V[i], PP));
  print("=== ", name, ":  order-", R, " PF fit, kernel dim vs degree D  (terms to n=", NN-1, ")");
  for(D=2, 40,
    my(nv=(R+1)*(D+1), nhi=min(NN-1, nv+60));
    if(nhi+1 >= nv+30,
      my(K = #matker(matconcat(vector(nhi+1, i, pfrow(Vm,R,D,i-1))~)));
      if(K>0, print("    D=", D, " : kernel dim ", K, "   [equations ", nhi+1, ", unknowns ", nv, "]"))));
}

pfexact(name, V, R, D) =
{ my(NN=#V, nv, nhi, K, VV, P, den, gg, P5, P4, rts, ff);
  nv = (R+1)*(D+1); nhi = min(NN-1, nv+80);
  K = matker(matconcat(vector(nhi+1, i, pfrow(V,R,D,i-1))~));
  print("--- ", name, ": exact order-", R, " degree-", D, " fit, kernel dim = ", #K);
  if(#K==0, return(0));
  VV = K[,1];
  P = vector(R+1, k, sum(j=0,D, VV[(k-1)*(D+1)+j+1]*'t^j));
  den = 1; for(k=1,R+1, den = lcm(den, denominator(content(P[k]))));
  P = vector(R+1, k, P[k]*den);
  gg = 0; for(k=1,R+1, gg = gcd(gg, content(P[k])));
  if(gg>0, P = vector(R+1, k, P[k]/gg));
  P5 = P[R+1]; P4 = P[R];
  print("    deg P_k, k=0..",R," : ", vector(R+1, k, poldegree(P[k])));
  print("    P_",R,"(t) = ", P5);
  print("    P_",R,"(t) FACTORED over Q: ", factor(P5));
  print("    Wronskian:  W'/W = -P_",R-1,"/P_",R,",   P_",R-1,"(t) = ", P4);
  print("    2-adic factorisation of P_",R,"(t)  [each factor: deg, v_2 of its roots]:");
  ff = factorpadic(P5/'t^valuation(P5,'t), 2, 30);
  print("      t-multiplicity at 0 (the cusp): ", valuation(P5,'t));
  for(i=1, matsize(ff)[1],
    my(f=ff[i,1], e=ff[i,2], d=poldegree(f),
       v0=valuation(polcoeff(f,0),2), vd=valuation(polcoeff(f,d),2));
    print("      deg ", d, " (mult ", e, ")  v_2(root) = ", (v0-vd)*1.0/d, "    ", f));
  my(vmin = 10^9);
  for(i=1, matsize(ff)[1],
    my(f=ff[i,1], d=poldegree(f), v0=valuation(polcoeff(f,0),2), vd=valuation(polcoeff(f,d),2), vv=(v0-vd)/d);
    if(vv < vmin, vmin = vv));
  print("    ==> nearest 2-adic singularity has v_2(t) = ", vmin, "  i.e. |t|_2 = 2^", -vmin);
  print("    ==> PREDICTED SLOPE sigma_2 = ", vmin, "   (rho = 2^", vmin, ")");
  print("    archimedean singularities |t| : ", vecsort(vector(poldegree(P5), i, 0)));
  rts = polroots(P5);
  print("    complex roots of P_",R,": ", vecsort(vector(#rts, i, abs(rts[i]))));
  P;
}

WW  = q*(ee(1)*ee(12)/(ee(3)*ee(4)))^4 + O(q^(M+1));
TD  = WW/(1+WW)^2;
HH  = (ee(1)^3*ee(4)*ee(6)^2/(ee(2)^2*ee(3)*ee(12)^3))/q + O(q^(M+1));
TX  = HH/((HH+3)*(HH+4));
XX  = q*ee(2)*ee(16)^2/(ee(1)^2*ee(8)) + O(q^(M+1));
T16 = XX/(8*XX^2 + 2*XX + 1);
c12 = [1,-113,567,112,-1863,1296]; d12 = [1,2,3,4,6,12];
c16 = [1,-85,1428,-5440,4096];     d16 = [1,2,4,8,16];

A12D = mkrow(c12, d12, TD);
A12X = mkrow(c12, d12, TX);
A16  = mkrow(c16, d16, T16);
print("rows built to n=", M-1);
pfscan("level-12 Domb", A12D, 5);
pfscan("level-12 h_12", A12X, 5);
pfscan("level-16     ", A16 , 5);
print("");
print("SCAN DONE");
quit;
