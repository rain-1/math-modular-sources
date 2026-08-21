/* ============================================================================
   The level-16 zeta(5) modular Apery row.
     source  Phi = -(1/504)[E6 - 85E6(2t) + 1428E6(4t) - 5440E6(8t) + 4096E6(16t)]
             = C^full_{5,2} E_6,  C^full_{5,2} = (1-V2)(1-4V2)(1-16V2)(1-64V2);
             weight 6, level 16, trivial character; kills s=0,2,4,6; L(Phi,5)=(217/1024)zeta(5).
     host    x = eta(2t)eta(16t)^2/(eta(t)^2 eta(8t))  (degree-1 Gamma_0(16) Hauptmodul,
             x|W16 = 1/(8x));  Fricke-invariant w = x/(1+8x^2);  t = w/(1+2w) = x/(8x^2+2x+1).
     rows    F = Phi/(Dt) (weight 4);  A(t)=F=sum A_n t^n;  B(t)=F*D^{-5}Phi=sum B_n t^n.
   Log: lattice/zeta5_two_row/level16.log
   ========================================================================= */
default(parisizemax, 12000000000);
\p 900
LOG = "lattice/zeta5_two_row/level16.log";
W(s) = write(LOG, s);
if(type(MM)!="t_INT", MM = 200);
M = MM;
dd(n) = if(n<1, 1, lcm(vector(n,i,i)));
IDX = [10,25,50,75,100,150,200,300,M-1];

q = 'q;
ee(d) = eta(q^d + O(q^(M+2)));      /* prod_{n>=1}(1-q^{dn}) = eta(d tau)/q^{d/24} */
E1 = ee(1); E2 = ee(2); E8 = ee(8); E16 = ee(16);
X  = q * E2*E16^2/(E1^2*E8) + O(q^(M+1));
T  = X/(8*X^2 + 2*X + 1);
W(Str("=== level-16 zeta(5) row, q-order M=", M, " ==="));
W(Str("x = ", truncate(X + O(q^9))));
W(Str("t = ", truncate(T + O(q^9))));

cv = [1,-85,1428,-5440,4096]; dv = [1,2,4,8,16];
g = vector(M);
for(n=1, M, s=0; for(i=1,5, if(n%dv[i]==0, s += cv[i]*sigma(n/dv[i],5))); g[n]=s);
PHI   = sum(n=1,M, g[n]*q^n) + O(q^(M+1));
THETA = sum(n=1,M, (g[n]/n^5)*q^n) + O(q^(M+1));
W(Str("Phi = ", truncate(PHI + O(q^9))));

DT = q*deriv(T, q);
F  = PHI/DT;
W(Str("F = Phi/Dt = ", truncate(F + O(q^9))));

QT = serreverse(T);
A  = subst(F, q, QT);
B  = subst(F*THETA, q, QT);
An = vector(M, i, polcoeff(A, i-1));
Bn = vector(M, i, polcoeff(B, i-1));
W(Str("A_0..A_10 = ", vector(11,i,An[i])));
W(Str("B_0..B_6  = ", vector(7,i,Bn[i])));

badA=0; badB=0; mA=1; mB=1;
for(n=0,M-1, if(denominator(An[n+1])!=1, badA++; mA=lcm(mA,denominator(An[n+1]))); \
             if(denominator(dd(n)^5*Bn[n+1])!=1, badB++; mB=lcm(mB,denominator(dd(n)^5*Bn[n+1]))));
W(Str("# n<",M," with A_n not in Z: ", badA, "   lcm denom(A_n): ", mA));
W(Str("# n<",M," with d_n^5 B_n not in Z: ", badB, "   lcm of residual denoms: ", mB));

z5 = zeta(5); L = 217/1024*z5;
W("");
W(Str("--- archimedean; predicted L(Phi,5) = (217/1024) zeta(5) = ", L*1.0));
for(i=1,#IDX, n=IDX[i]; if(n<=M-1 && n>0 && An[n+1]!=0, \
  W(Str("  n=",n,"  B_n/A_n=", Bn[n+1]/An[n+1]*1.0, "   diff=", (Bn[n+1]/An[n+1]-L)*1.0, \
        "   |diff|^(1/n)=", abs(Bn[n+1]/An[n+1]-L)^(1.0/n)))));
W("");
W("--- growth rates: log|A_n|/n (predict log(4sqrt2+2)=2.03574), log|B-LA|/n (predict log4=1.38629) ---");
for(i=1,#IDX, n=IDX[i]; if(n<=M-1 && n>0 && An[n+1]!=0, \
  W(Str("  n=",n,"  log|A|/n=", log(abs(An[n+1]*1.0))/n, "   log|B-LA|/n=", log(abs(Bn[n+1]*1.0-L*An[n+1]))/n))));

W("");
W("--- v_p(A_n) (kappa test) ---");
for(j=1,4, p=[2,3,5,7][j]; s=""; \
  for(i=1,#IDX, n=IDX[i]; if(n<=M-1 && n>0 && An[n+1]!=0, s=Str(s," ",valuation(An[n+1],p)))); \
  W(Str("  p=",p,":",s)));
W("");
W("--- SLOPE: v_p( B_n/A_n - B_{n-1}/A_{n-1} ) ---");
for(j=1,4, p=[2,3,5,7][j]; s=""; \
  for(i=1,#IDX, n=IDX[i]; if(n<=M-1 && n>1 && An[n+1]!=0 && An[n]!=0, s=Str(s," ",valuation(Bn[n+1]/An[n+1]-Bn[n]/An[n],p)))); \
  W(Str("  p=",p,":",s)));
W("--- fine v_2 increments, last 40 n ---");
s=""; for(n=max(2,M-41), M-1, if(An[n+1]!=0 && An[n]!=0, s=Str(s," ",valuation(Bn[n+1]/An[n+1]-Bn[n]/An[n],2))));
W(s);

write("lattice/zeta5_two_row/level16_rows.txt", "An = ", An);
write("lattice/zeta5_two_row/level16_rows.txt", "Bn = ", Bn);
W("");
W("DONE");
quit;
