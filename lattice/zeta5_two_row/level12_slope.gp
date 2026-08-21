/* ---------------------------------------------------------------------------
   The NEW level-12 zeta(5) source with a 2-adic slope:
       c_d = (1, -113, 567, 112, -1863, 1296)  over d = (1,2,3,4,6,12)
       P(0)=P(2)=P(4)=P(6)=0,  P(5) = -31/96,  L(Phi,5) = (31/192) zeta(5).
   Euler-factor criterion at p=2 holds: 1-113+112 = 0 and 567-1863+1296 = 0.
   Prediction:  P(X,Z) = (1-X) Q(X,Z),  Q = (1-112X) + Z(567-1296X),
                Q(2^-5, 3^-5) = -1/3,   xi_2 = -Q(s=5)*(1/2)zeta_2(5) = (1/6) zeta_2(5).
   Conjecture D check:  (1/6)/(31/192) = 32/31 = (7/32)/(217/1024)  -- same ratio as level 16.
   Built on both level-12 coordinates:
     (a) Domb    w = (eta1 eta12/(eta3 eta4))^4,  t = w/(1+w)^2
     (b) h_12    h = eta1^3 eta4 eta6^2/(eta2^2 eta3 eta12^3),  x = h/((h+3)(h+4))
   Log: lattice/zeta5_two_row/level12_slope.log
   ------------------------------------------------------------------------- */
default(parisizemax, 6000000000);
\p 900
LOG = "lattice/zeta5_two_row/level12_slope.log";
W(s) = write(LOG, s);
if(type(MM)!="t_INT", MM = 200);
M = MM;
dd(n) = if(n<1, 1, lcm(vector(n,i,i)));
q = 'q; ee(d) = eta(q^d + O(q^(M+2)));
IDX = [10,25,50,75,100,150,200,250,M-1];

cv = [1,-113,567,112,-1863,1296]; dv = [1,2,3,4,6,12];
g = vector(M);
for(n=1,M, s=0; for(i=1,#dv, if(n%dv[i]==0, s += cv[i]*sigma(n/dv[i],5))); g[n]=s);
PHI   = sum(n=1,M, g[n]*q^n) + O(q^(M+1));
THETA = sum(n=1,M, (g[n]/n^5)*q^n) + O(q^(M+1));
W(Str("=== NEW level-12 zeta(5) row (31/192), M=", M, " ==="));
W(Str("Phi = ", truncate(PHI + O(q^10))));

WW = q*(ee(1)*ee(12)/(ee(3)*ee(4)))^4 + O(q^(M+1));
TD = WW/(1+WW)^2;
HH = (ee(1)^3*ee(4)*ee(6)^2/(ee(2)^2*ee(3)*ee(12)^3))/q + O(q^(M+1));
TX = HH/((HH+3)*(HH+4));

doit(name, T, tag) = { \
  my(DT, F, QT, A, B, An, Bn, badA, badB, mA, L, n, s, p, r); \
  DT = q*deriv(T,q); F = PHI/DT; QT = serreverse(T); \
  A = subst(F,q,QT); B = subst(F*THETA,q,QT); \
  An = vector(M,i,polcoeff(A,i-1)); Bn = vector(M,i,polcoeff(B,i-1)); \
  W(""); W(Str("#### coordinate: ", name)); \
  W(Str("  A_0..A_8 = ", vector(9,i,An[i]))); \
  W(Str("  B_0..B_4 = ", vector(5,i,Bn[i]))); \
  badA=0; badB=0; mA=1; \
  for(n=0,M-1, if(denominator(An[n+1])!=1, badA++; mA=lcm(mA,denominator(An[n+1]))); \
               if(denominator(dd(n)^5*Bn[n+1])!=1, badB++)); \
  W(Str("  # A_n not in Z: ", badA, " (lcm denom ", mA, ")   # d_n^5 B_n not in Z: ", badB)); \
  L = 31/192*zeta(5); \
  for(i=1,#IDX, n=IDX[i]; if(n<=M-1 && n>0 && An[n+1]!=0, r=Bn[n+1]/An[n+1]; \
     W(Str("    n=",n,"  B/A-L=",(r-L)*1.0,"  |.|^(1/n)=",abs(r-L)^(1.0/n), \
           "  log|A|/n=",log(abs(An[n+1]*1.0))/n, \
           "  log|B-LA|/n=",log(abs(Bn[n+1]*1.0-L*An[n+1]))/n)))); \
  W("  v_p(A_n):"); \
  for(j=1,4, p=[2,3,5,7][j]; s=""; for(i=1,#IDX, n=IDX[i]; if(n<=M-1 && An[n+1]!=0, s=Str(s," ",valuation(An[n+1],p)))); W(Str("    p=",p,":",s))); \
  W("  SLOPE v_p(B_n/A_n - B_{n-1}/A_{n-1}):"); \
  for(j=1,4, p=[2,3,5,7][j]; s=""; for(i=1,#IDX, n=IDX[i]; if(n<=M-1 && n>1 && An[n+1]!=0 && An[n]!=0, s=Str(s," ",valuation(Bn[n+1]/An[n+1]-Bn[n]/An[n],p)))); W(Str("    p=",p,":",s))); \
  write(Str("lattice/zeta5_two_row/level12_rows_",tag,".txt"), Str("An",tag," = "), An); \
  write(Str("lattice/zeta5_two_row/level12_rows_",tag,".txt"), Str("Bn",tag," = "), Bn); \
  0; };

doit("Domb  t = w/(1+w)^2", TD, "D");
doit("h_12  x = h/((h+3)(h+4))", TX, "X");
W("");
W("DONE"); quit;
