/* ===========================================================================
   zeta5_check1.gp  --  Which hypothesis of Theorem F fails for the
   level-12 zeta(5) host?   (set MM before reading; default 420)

   Source  Phi_12 = sum c_d E_6(d tau),  c=(1,-113,567,112,-1863,1296) on
   d=(1,2,3,4,6,12).  Mellin polynomial P(s)=sum c_d d^-s.
     (a) E_2(s)=1-2^-s divides P(s)         -> checked symbolically
     (c) t*j in Z_(2)[[t]], unit const term -> checked on the t-expansion
     (d) v_2(a_n)=O(log n), sigma_2>0       -> measured to n=MM
   Predicted xi* = -Q(5)*kappa_2, kappa_2=(1/2)zeta_2(5).
   Rows: level-12 Domb, level-12 h_12, level-16 (contrast).
   =========================================================================== */
default(parisizemax, 24000000000);
\p 120
read("lattice/euler_criterion/lp.gp");
if(type(MM)!="t_INT", MM = 420);
M = MM;
PR = M + 320;
q = 'q;
ee(d) = eta(q^d + O(q^(M+2)));
dd(n) = if(n<1, 1, lcm(vector(n,i,i)));

print("############ PART 1: hypothesis (a), symbolic ############");
{
  my(cv=[1,-113,567,112,-1863,1296], dv=[1,2,3,4,6,12], P, Q, R);
  P = sum(i=1,#dv, cv[i] * 'X^valuation(dv[i],2) * 'Z^valuation(dv[i],3));
  print("P(X,Z) = ", P, "      (X=2^-s, Z=3^-s)");
  R = divrem(P, 1-'X, 'X);
  print("  P/(1-X) = ", R[1], "     remainder = ", R[2]);
  print("  (a) at p=2 HOLDS ? ", R[2]==0);
  print("  (a) at p=3 [(1-Z)|P] ? ", divrem(P,1-'Z,'Z)[2]==0);
  QQ5 = subst(subst(R[1],'X,1/32),'Z,1/243);
  PP5 = subst(subst(P,'X,1/32),'Z,1/243);
  print("  Q(s=5) = ", QQ5, "     P(s=5) = ", PP5, "     L(Phi,5) = ", -PP5/2, "*zeta(5)");
}
print("");
print("############ PART 1b: kappa_2, xi* ############");
Z25 = Lp(2, triv, 5, PR);
print("  v_2(zeta_2(5)) = ", valuation(Z25,2));
{
  for(n=1, 6,
    my(k=2*n, lhs, rhs);
    lhs = Lp(2,triv,1-k,60);
    rhs = -(1-2^(k-1))*bernfrac(k)/k;
    print("    L_2(1-",k,",1) = ",rhs," matches interpolation: ", valuation(lhs-rhs,2)>50));
}
KAP2   = Z25/2;
XISTAR = -QQ5*KAP2;
XI16   = 7/32*Z25;
print("  xi*(level 12) = -Q(5)*kappa_2 = ", -QQ5/2, "*zeta_2(5),   v_2 = ", valuation(XISTAR,2));
print("  xi*(level 16) = 7/32*zeta_2(5),                v_2 = ", valuation(XI16,2));
print("");

print("############ PART 2: build the rows ############");
buildrow(name, cv, dv, T) =
{ my(g, PHI, THETA, DT, F, QT, A, B, An, Bn, badA, badB, mA, s, TT);
  g = vector(M);
  for(n=1,M, s=0; for(i=1,#dv, if(n%dv[i]==0, s += cv[i]*sigma(n/dv[i],5))); g[n]=s);
  PHI   = sum(n=1,M, g[n]*q^n) + O(q^(M+1));
  THETA = sum(n=1,M, (g[n]/n^5)*q^n) + O(q^(M+1));
  DT = q*deriv(T,q); F = PHI/DT; QT = serreverse(T);
  A = subst(F,q,QT); B = subst(F*THETA,q,QT);
  An = vector(M,i,polcoeff(A,i-1)); Bn = vector(M,i,polcoeff(B,i-1));
  badA=0; badB=0; mA=1;
  for(n=0,M-1, if(denominator(An[n+1])!=1, badA++; mA=lcm(mA,denominator(An[n+1])));
               if(denominator(dd(n)^5*Bn[n+1])!=1, badB++));
  TT = truncate(T + O(q^(M+1)));
  print("--- ", name);
  print("      t = ", truncate(T+O(q^7)), " ...   integral&monic: ",
        (denominator(content(TT))==1) && (polcoeff(T,1)==1) && (polcoeff(T,0)==0));
  print("      A_0..A_8 = ", vector(9,i,An[i]));
  print("      #A_n notin Z: ", badA, " (lcm denom ",mA,")   #d_n^5 B_n notin Z: ", badB);
  print("      log|A_n|/n at n=",M-1," : ", log(abs(An[M]*1.0))/(M-1));
  [An, Bn];
}

WW  = q*(ee(1)*ee(12)/(ee(3)*ee(4)))^4 + O(q^(M+1));
TD  = WW/(1+WW)^2;
HH  = (ee(1)^3*ee(4)*ee(6)^2/(ee(2)^2*ee(3)*ee(12)^3))/q + O(q^(M+1));
TX  = HH/((HH+3)*(HH+4));
XX  = q*ee(2)*ee(16)^2/(ee(1)^2*ee(8)) + O(q^(M+1));
T16 = XX/(8*XX^2 + 2*XX + 1);

c12 = [1,-113,567,112,-1863,1296]; d12 = [1,2,3,4,6,12];
c16 = [1,-85,1428,-5440,4096];     d16 = [1,2,4,8,16];

R12D = buildrow("level-12 Domb  t=w/(1+w)^2", c12, d12, TD);
R12X = buildrow("level-12 h_12  t=h/((h+3)(h+4))", c12, d12, TX);
R16  = buildrow("level-16       t=x/(8x^2+2x+1)", c16, d16, T16);
print("");

print("############ PART 3: v_2(a_n)  --  hypothesis (d), H1 ############");
vprof(name, An) =
{ my(v, mx=0, nmx=0, mn=1000000, rmax=0.0, nr=0, rmin=1000.0, nrmin=0, r, IDX, tl);
  v = vector(M-1, n, if(An[n+1]!=0, valuation(An[n+1],2), -1));
  for(n=1,M-1, if(An[n+1]!=0, if(v[n]>mx, mx=v[n]; nmx=n); if(v[n]<mn, mn=v[n])));
  IDX = [50,100,150,200,250,300,350,400,M-1];
  print("--- ", name);
  print("   v_2(a_n), n=1..40 : ", vector(40,n,v[n]));
  print("   n=", IDX, " : ", vector(#IDX, i, if(IDX[i]<=M-1, v[IDX[i]], -999)));
  print("   min v_2 = ", mn, "    max v_2 = ", mx, " at n=", nmx, "    log_2(M) = ", log(M*1.0)/log(2));
  for(n=20,M-1, if(An[n+1]!=0, r=1.0*v[n]/n;
      if(r>rmax, rmax=r; nr=n); if(r<rmin, rmin=r; nrmin=n)));
  print("   20<=n<=",M-1,":   limsup-est max v_2(a_n)/n = ", rmax, " at n=", nr,
        "    liminf-est min v_2(a_n)/n = ", rmin, " at n=", nrmin);
  tl = max(1, M-120);
  print("   tail ", tl, "<=n<=", M-1, ":  max v_2/n = ", vecmax(vector(M-tl, i, 1.0*v[tl+i-1]/(tl+i-1))),
        "   mean v_2 = ", (1.0/(M-tl))*sum(i=1,M-tl, v[tl+i-1]));
  print("   max v_2(a_n)/log_2(n) over 2<=n<=",M-1," = ",
        vecmax(vector(M-2, n, 1.0*v[n+1]*log(2)/log(n+1.0))));
  print("   fine tail v_2(a_n), n=",M-21,"..",M-1," : ", vector(20, i, v[M-21+i]));
}
vprof("level-12 Domb", R12D[1]);
vprof("level-12 h_12", R12X[1]);
vprof("level-16     ", R16[1]);
print("");

print("############ PART 4: v_2(b_n - xi* a_n) and slope increments ############");
resid(name, An, Bn, XI) =
{ my(IDX, n, r, rr, inc, va);
  print("--- ", name, "    v_2(xi*) = ", valuation(XI,2));
  print("     n | v_2(b_n - xi* a_n) | v_2(a_n) | v_2(b_n/a_n - xi*) | v_2(b_n/a_n - b_{n-1}/a_{n-1})");
  IDX = [10,20,40,60,80,100,140,180,220,260,300,340,380,M-1];
  IDX = select(x->x<=M-1, IDX);
  for(i=1,#IDX, n=IDX[i];
    if(n<=M-1 && An[n+1]!=0,
      r  = valuation(Bn[n+1]-XI*An[n+1], 2);
      va = valuation(An[n+1],2);
      rr = valuation(Bn[n+1]/An[n+1]-XI, 2);
      inc = if(n>1 && An[n]!=0, valuation(Bn[n+1]/An[n+1]-Bn[n]/An[n],2), -999);
      print("   ", n, "  |  ", r, "  |  ", va, "  |  ", rr, "  |  ", inc)));
}
resid("level-12 Domb, xi*=(1/6)zeta_2(5)", R12D[1], R12D[2], XISTAR);
resid("level-12 h_12, xi*=(1/6)zeta_2(5)", R12X[1], R12X[2], XISTAR);
resid("level-16, xi*=(7/32)zeta_2(5)",     R16[1],  R16[2],  XI16);
print("");

print("############ PART 5: hypothesis (c): t*j in Z_(2)[[t]] with unit constant term ############");
JJ = ellj(q + O(q^(M+1)));
jcheck(name, T) =
{ my(QT, TJ, cf, bad=0, mn=0, v, K);
  K = min(60, M-2);
  QT = serreverse(T);
  TJ = subst(T*JJ, q, QT);
  cf = vector(K, i, polcoeff(TJ, i-1));
  print("--- ", name);
  print("   (t j)(t), first 8 coefficients: ", vector(8,i,cf[i]));
  print("   constant term = ", cf[1], "   2-adic unit: ", valuation(cf[1],2)==0);
  for(i=1,K, if(denominator(cf[i])!=1, bad++));
  for(i=1,K, if(cf[i]!=0, v=valuation(cf[i],2); if(v<mn, mn=v)));
  print("   # of first ",K," coefficients not in Z: ", bad, "    min v_2 = ", mn);
  print("   => (c) as literally stated: ", if(bad==0 && valuation(cf[1],2)==0, "HOLDS", "FAILS"));
}
jcheck("level-12 Domb", TD);
jcheck("level-12 h_12", TX);
jcheck("level-16     ", T16);
print("");
print("DONE");
quit;
