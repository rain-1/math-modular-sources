/* ===========================================================================
   zeta5_check2.gp -- the HOST side: minimal recurrence, characteristic
   polynomial and 2-adic Newton polygon for the level-12 zeta(5) A-rows
   (Domb and h_12 coordinates) and, for contrast, the level-16 A-row.

   The characteristic polynomial's roots lambda_i are the reciprocals of the
   singular points of the t-line Picard-Fuchs operator; the p-adic radius of
   the overconvergent disc is rho = min_i |1/lambda_i|_p = p^{min_i v_p(lambda_i)},
   and the predicted slope is sigma_p = log_p rho = min_i v_p(lambda_i).
   Set MM before reading (default 300).
   =========================================================================== */
default(parisizemax, 24000000000);
default(threadsizemax, 8000000000);
\p 60
if(type(MM)!="t_INT", MM = 300);
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

mk(V,r,D,nlo,nhi) = matconcat(vector(nhi-nlo+1, i, my(n=nlo+i-1); vector((r+1)*(D+1), c, my(j=(c-1)\(D+1), k=(c-1)%(D+1)); n^k * V[n-j+1]))~);

scan(name, An) =
{ my(NN=#An, Am, s, best=0);
  Am = vector(NN, i, Mod(An[i], PP));
  print("=== ", name, ":  mod-p kernel-dimension scan, rows to n=", NN-1);
  for(rr=4,16,
    s="";
    for(DD=2,12,
      my(nv=(rr+1)*(DD+1), nhi=min(NN-1, rr+1+nv+40));
      if(nhi-rr-1 >= nv+25, s=Str(s,"  D=",DD,":",#matker(mk(Am,rr,DD,rr+1,nhi))), s=Str(s,"  D=",DD,":-")));
    print("   r=",rr,s));
}

recon(name, An, RR, DDD) =
{ my(nv, nhi, K, VV, cf, den, gg, lead, cpol, trail, res, lam, NP);
  nv = (RR+1)*(DDD+1); nhi = min(#An-1, RR+1+nv+60);
  K = matker(mk(An, RR, DDD, RR+1, nhi));
  print("--- ", name, "  exact reconstruction at (r,D)=(",RR,",",DDD,"), kernel dim = ", #K);
  if(#K==0, return(0));
  VV = K[,1];
  cf = vector(RR+1, j, sum(k=0,DDD, VV[(j-1)*(DDD+1)+k+1]*'x^k));
  den = 1; for(j=1,RR+1, den = lcm(den, denominator(content(cf[j]))));
  cf = vector(RR+1, j, cf[j]*den);
  gg = 0; for(j=1,RR+1, gg = gcd(gg, content(cf[j])));
  if(gg>0, cf = vector(RR+1, j, cf[j]/gg));
  res = vecmax(vector(30, i, my(n=#An-31+i); abs(sum(j=0,RR, subst(cf[j+1],'x,n)*An[n-j+1]))));
  print("    residual on the last 30 indices: ", res);
  print("    deg p_k, k=0..r : ", vector(RR+1, j, poldegree(cf[RR+2-j])));
  lead = vector(RR+1, j, polcoeff(cf[j], DDD));
  cpol = sum(j=0, RR, lead[j+1]*'y^(RR-j));
  print("    characteristic polynomial (top-degree coefficients), factored:");
  print("      ", factor(cpol));
  print("    c = prod(lambda_i) = (-1)^r p_r-top / p_0-top   ->  lead[r+1]/lead[1] = ", lead[RR+1]/lead[1]);
  print("    v_2 of that ratio = ", if(lead[RR+1]!=0 && lead[1]!=0, valuation(lead[RR+1]/lead[1],2), "n/a"));
  print("    2-adic factorisation of the characteristic polynomial:");
  my(cp = cpol/'y^valuation(cpol,'y));
  print("      y-multiplicity at 0: ", valuation(cpol,'y));
  my(FF = factorpadic(cp, 2, 25));
  for(i=1, matsize(FF)[1],
     my(f=FF[i,1], e=FF[i,2], d=poldegree(f), v0=valuation(polcoeff(f,0),2), vd=valuation(polcoeff(f,d),2));
     print("      factor deg ", d, " exp ", e, "  v_2(roots) = ", (v0-vd)*1.0/d, "   ", f));
  cf;
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
print("A12D[1..9] = ", vector(9,i,A12D[i]));
print("A12X[1..9] = ", vector(9,i,A12X[i]));
print("A16 [1..9] = ", vector(9,i,A16[i]));
print("");
scan("level-12 Domb", A12D);
print("");
scan("level-12 h_12", A12X);
print("");
print("DONE-SCAN");
quit;
