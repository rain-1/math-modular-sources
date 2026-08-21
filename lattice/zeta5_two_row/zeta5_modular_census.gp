/* ---------------------------------------------------------------------------
   Census of trivial-character weight-6 zeta(5) modular Apery rows:
     level 16  P(X)=(1-X)(1-4X)(1-16X)(1-64X)                L(Phi,5)=(217/1024)z5
     level 12- (anti-Fricke)  (1,-176,2079,-4928,4752,-1728)  L(Phi,5)=(11/144)z5
     level 12+ (Domb)         (1,-104, 351,  832,-2808, 1728) L(Phi,5)=(25/144)z5
   For each: build A_n, B_n via  F = Phi/(Dt),  A(t)=F,  B(t)=F*D^{-5}Phi,
   and measure integrality, archimedean rates, v_p(A_n), and the p-adic slope.
   Also evaluates the p-adic Euler-factor criterion:  for each prime-to-p class m,
   sum_{a} c_{p^a m}  must vanish  (<=> (1-V_p)-divisibility of the p-part).
   Log: lattice/zeta5_two_row/census.log
   ------------------------------------------------------------------------- */
default(parisizemax, 12000000000);
\p 600
LOG = "lattice/zeta5_two_row/census.log";
W(s) = write(LOG, s);
if(type(MM)!="t_INT", MM = 150);
M = MM;
dd(n) = if(n<1, 1, lcm(vector(n,i,i)));
q = 'q;
ee(d) = eta(q^d + O(q^(M+2)));

/* --- Euler-factor criterion ------------------------------------------------ */
crit(cv, dv, p) = my(ms = Set(vector(#dv, i, dv[i]/p^valuation(dv[i],p))), out = []); \
  for(i=1,#ms, my(m = ms[i], s = 0); \
     for(j=1,#dv, if(dv[j]/p^valuation(dv[j],p) == m, s += cv[j])); \
     out = concat(out, [[m, s]])); out;

/* --- the row builder ------------------------------------------------------- */
build(name, cv, dv, T, Lfac) = \
{ my(g, PHI, THETA, DT, F, QT, A, B, An, Bn, badA, badB, mA, L, IDX, n, s, p, r, ratefit); \
  W(""); W(Str("################  ", name, "  ################")); \
  W(Str("  source c_d = ", cv, "   over d = ", dv)); \
  for(pp=1,2, p = [2,3][pp]; W(Str("  Euler-factor criterion at p=", p, \
      " (need sum of coeffs = 0 in every prime-to-p class): ", crit(cv,dv,p)))); \
  g = vector(M); \
  for(n=1, M, s = 0; for(i=1,#dv, if(n%dv[i]==0, s += cv[i]*sigma(n/dv[i],5))); g[n] = s); \
  PHI   = sum(n=1,M, g[n]*q^n) + O(q^(M+1)); \
  THETA = sum(n=1,M, (g[n]/n^5)*q^n) + O(q^(M+1)); \
  DT = q*deriv(T, q); \
  F  = PHI/DT; \
  W(Str("  t = ", truncate(T + O(q^8)))); \
  W(Str("  F = Phi/Dt = ", truncate(F + O(q^8)))); \
  QT = serreverse(T); \
  A  = subst(F, q, QT); B = subst(F*THETA, q, QT); \
  An = vector(M, i, polcoeff(A, i-1)); Bn = vector(M, i, polcoeff(B, i-1)); \
  W(Str("  A_0..A_8 = ", vector(9,i,An[i]))); \
  badA = 0; badB = 0; mA = 1; \
  for(n=0,M-1, if(denominator(An[n+1])!=1, badA++; mA = lcm(mA,denominator(An[n+1]))); \
               if(denominator(dd(n)^5*Bn[n+1])!=1, badB++)); \
  W(Str("  # A_n not in Z: ", badA, " (lcm denom ", mA, ")    # d_n^5 B_n not in Z: ", badB)); \
  L = Lfac*zeta(5); \
  IDX = [10,25,50,75,100,M-1]; \
  W(Str("  predicted archimedean limit ", Lfac, "*zeta(5)")); \
  for(i=1,#IDX, n = IDX[i]; if(n<=M-1 && An[n+1]!=0, r = Bn[n+1]/An[n+1]; \
     W(Str("    n=", n, "  B/A - L = ", (r-L)*1.0, "   |.|^(1/n)=", abs(r-L)^(1.0/n), \
           "   log|A|/n=", log(abs(An[n+1]*1.0))/n)))); \
  W("  v_p(A_n) at n=10,25,50,75,100,M-1:"); \
  for(j=1,4, p = [2,3,5,7][j]; s = ""; \
     for(i=1,#IDX, n = IDX[i]; if(n<=M-1 && An[n+1]!=0, s = Str(s," ",valuation(An[n+1],p)))); \
     W(Str("    p=",p,":",s))); \
  W("  SLOPE  v_p(B_n/A_n - B_{n-1}/A_{n-1}):"); \
  for(j=1,4, p = [2,3,5,7][j]; s = ""; \
     for(i=1,#IDX, n = IDX[i]; if(n<=M-1 && n>1 && An[n+1]!=0 && An[n]!=0, \
        s = Str(s," ",valuation(Bn[n+1]/An[n+1]-Bn[n]/An[n],p)))); \
     W(Str("    p=",p,":",s))); \
  [An, Bn]; };

W(Str("=== zeta(5) modular row census, M=", M, " ==="));

/* ---- level 16 ---- */
X16 = q * ee(2)*ee(16)^2/(ee(1)^2*ee(8)) + O(q^(M+1));
T16 = X16/(8*X16^2 + 2*X16 + 1);
R16 = build("LEVEL 16  (C^full_{5,2} E_6)", [1,-85,1428,-5440,4096], [1,2,4,8,16], T16, 217/1024);

/* ---- level 12, Domb coordinate  w = (eta1 eta12/(eta3 eta4))^4,  t = w/(1+w)^2 ---- */
W12 = q * (ee(1)*ee(12)/(ee(3)*ee(4)))^4 + O(q^(M+1));
T12 = W12/(1+W12)^2;
R12p = build("LEVEL 12 (+, Domb, 25/144)", [1,-104,351,832,-2808,1728], [1,2,3,4,6,12], T12, 25/144);
R12m = build("LEVEL 12 (-, anti-Fricke, 11/144) on the SAME Domb coordinate", \
             [1,-176,2079,-4928,4752,-1728], [1,2,3,4,6,12], T12, 11/144);

/* ---- level 12, h_12 coordinate: h = eta1^3 eta4 eta6^2/(eta2^2 eta3 eta12^3), x=h/((h+3)(h+4)) ---- */
H12 = (ee(1)^3*ee(4)*ee(6)^2/(ee(2)^2*ee(3)*ee(12)^3))/q + O(q^(M+1));
X12 = H12/((H12+3)*(H12+4));
R12x = build("LEVEL 12 (-, anti-Fricke, 11/144) on the h_12/x coordinate", \
             [1,-176,2079,-4928,4752,-1728], [1,2,3,4,6,12], X12, 11/144);

/* ---- level 18 source, for the criterion only ---- */
W("");
W("### LEVEL 18 source (1,-136,918,-7344,12393,-5832) over d|18 -- criterion only");
W(Str("  p=2: ", crit([1,-136,918,-7344,12393,-5832],[1,2,3,6,9,18],2)));
W(Str("  p=3: ", crit([1,-136,918,-7344,12393,-5832],[1,2,3,6,9,18],3)));
W("");
W("DONE");
quit;
