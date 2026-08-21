/* The level-16 zeta(5) source C^full_{5,2}E_6 on SEVERAL host coordinates,
   looking for a second row with a 2-adic slope (a possible alignment partner).
   Log: lattice/zeta5_two_row/level16_variants.log */
default(parisizemax, 6000000000);
\p 600
LOG = "lattice/zeta5_two_row/level16_variants.log";
W(s) = write(LOG, s);
if(type(MM)!="t_INT", MM = 150);
M = MM; dd(n)=if(n<1,1,lcm(vector(n,i,i)));
q='q; ee(d) = eta(q^d + O(q^(M+2)));
X = q*ee(2)*ee(16)^2/(ee(1)^2*ee(8)) + O(q^(M+1));
cv=[1,-85,1428,-5440,4096]; dv=[1,2,4,8,16];
g = vector(M); for(n=1,M, s=0; for(i=1,#dv, if(n%dv[i]==0, s += cv[i]*sigma(n/dv[i],5))); g[n]=s);
PHI = sum(n=1,M,g[n]*q^n)+O(q^(M+1)); THETA = sum(n=1,M,(g[n]/n^5)*q^n)+O(q^(M+1));
IDX = [10,25,50,75,100,M-1];
z5 = zeta(5);

doit(name, T) = { my(DT,F,QT,A,B,An,Bn,badA,badB,mA,n,s,p,lim,ld); \
  DT = q*deriv(T,q); F = PHI/DT; QT = serreverse(T); \
  A = subst(F,q,QT); B = subst(F*THETA,q,QT); \
  An = vector(M,i,polcoeff(A,i-1)); Bn = vector(M,i,polcoeff(B,i-1)); \
  W(""); W(Str("#### ", name)); W(Str("  t = ", truncate(T+O(q^7)))); \
  W(Str("  A_0..A_6 = ", vector(7,i,An[i]))); \
  badA=0; badB=0; mA=1; \
  for(n=0,M-1, if(denominator(An[n+1])!=1, badA++; mA=lcm(mA,denominator(An[n+1]))); \
               if(denominator(dd(n)^5*Bn[n+1])!=1, badB++)); \
  W(Str("  # A_n not in Z: ",badA," (lcm ",mA,")   # d_n^5 B_n not in Z: ",badB)); \
  lim = Bn[M]/An[M]; ld = lindep([lim*1.0, z5, 1], 40); \
  W(Str("  B_n/A_n at n=",M-1," = ", lim*1.0)); \
  W(Str("  lindep([B/A, zeta(5), 1]) = ", ld, "   => B/A ~ ", if(ld[1]!=0, -ld[2]/ld[1], "?"), " * zeta(5) + ", if(ld[1]!=0, -ld[3]/ld[1], "?"))); \
  W(Str("  log|A_n|/n = ", log(abs(An[M]*1.0))/(M-1))); \
  W("  v_p(A_n):"); \
  for(j=1,3, p=[2,3,5][j]; s=""; for(i=1,#IDX, n=IDX[i]; if(n<=M-1 && An[n+1]!=0, s=Str(s," ",valuation(An[n+1],p)))); W(Str("    p=",p,":",s))); \
  W("  SLOPE v_p(B_n/A_n - B_{n-1}/A_{n-1}):"); \
  for(j=1,3, p=[2,3,5][j]; s=""; for(i=1,#IDX, n=IDX[i]; if(n<=M-1 && n>1 && An[n+1]!=0 && An[n]!=0, s=Str(s," ",valuation(Bn[n+1]/An[n+1]-Bn[n]/An[n],p)))); W(Str("    p=",p,":",s))); \
  0; };

W(Str("=== level-16 zeta(5) source on several host coordinates, M=", M, " ==="));
doit("plain Gamma_0(16) Hauptmodul  x", X);
doit("Fricke quotient  t = x/(8x^2+2x+1)  [the canonical row]", X/(8*X^2+2*X+1));
doit("Fricke invariant  w = x/(1+8x^2)", X/(1+8*X^2));
doit("Catalan-notes  u = x(1+3x)/(1+4x)", X*(1+3*X)/(1+4*X));
doit("Catalan-notes  s = x(1+2x)/(1+4x)^2", X*(1+2*X)/(1+4*X)^2);
doit("shifted  x/(1+4x)", X/(1+4*X));
W(""); W("DONE"); quit;
