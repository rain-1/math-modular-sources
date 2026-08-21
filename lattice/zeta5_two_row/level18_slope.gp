/* level-18 zeta(5) source with the p=3 Euler-factor criterion:
   c_d = (1,-112,566,-1184,-567,1296) over d=(1,2,3,6,9,18); P(5)=-242/729; L=(121/729)zeta(5).
   3-part classes: m=1 (d=1,3,9): 1+566-567=0 ; m=2 (d=2,6,18): -112-1184+1296=0.  OK.
   Predicted xi_3 = -Q(s=5)*(1/2) zeta_3(5),  Q = P/(1-Z), Z = 3^-s.
   Log: lattice/zeta5_two_row/level18.log */
default(parisizemax, 6000000000);
\p 500
LOG = "lattice/zeta5_two_row/level18.log";
W(s) = write(LOG, s);
if(type(MM)!="t_INT", MM = 140);
M = MM; dd(n) = if(n<1,1,lcm(vector(n,i,i)));
q='q; ee(d) = eta(q^d + O(q^(M+2)));
N = 18; DV = divisors(N);
ordc(r,c) = (N/(24*gcd(c^2,N))) * sum(i=1,#DV, gcd(c,DV[i])^2*r[i]/DV[i]);
valid(r) = if(sum(i=1,#DV,r[i])!=0, 0, \
  if(sum(i=1,#DV,DV[i]*r[i])%24!=0, 0, \
  if(sum(i=1,#DV,(N/DV[i])*r[i])%24!=0, 0, \
  if(!issquare(prod(i=1,#DV,DV[i]^r[i])) && !issquare(1/prod(i=1,#DV,DV[i]^r[i])), 0, 1))));
found = List();
{ forvec(r = vector(#DV,i,[-9,9]), if(valid(r), my(o=vector(#DV,i,ordc(r,DV[i]))); \
   if(o[#DV]==1 && sum(i=1,#DV, if(o[i]<0,-o[i],0))==1, listput(found,[r,o])))); }
W(Str("=== level 18: degree-1 eta-quotient Hauptmoduls found: ", #found));
for(i=1,min(#found,6), W(Str("  r=",found[i][1]," cusp orders=",found[i][2])));

cv = [1,-112,566,-1184,-567,1296]; dv = [1,2,3,6,9,18];
g = vector(M);
for(n=1,M, s=0; for(i=1,#dv, if(n%dv[i]==0, s += cv[i]*sigma(n/dv[i],5))); g[n]=s);
PHI = sum(n=1,M,g[n]*q^n)+O(q^(M+1)); THETA = sum(n=1,M,(g[n]/n^5)*q^n)+O(q^(M+1));
IDX = [10,25,50,75,100,M-1];

doit(name, T) = { my(DT,F,QT,A,B,An,Bn,badA,mA,n,s,p); \
  DT = q*deriv(T,q); F = PHI/DT; QT = serreverse(T); \
  A = subst(F,q,QT); B = subst(F*THETA,q,QT); \
  An = vector(M,i,polcoeff(A,i-1)); Bn = vector(M,i,polcoeff(B,i-1)); \
  W(""); W(Str("#### ", name)); W(Str("  A_0..A_7 = ", vector(8,i,An[i]))); \
  badA=0; mA=1; for(n=0,M-1, if(denominator(An[n+1])!=1, badA++; mA=lcm(mA,denominator(An[n+1])))); \
  W(Str("  # A_n not in Z: ", badA, "  lcm denom ", mA)); \
  W(Str("  log|A_n|/n at n=", M-1, ": ", log(abs(An[M]*1.0))/(M-1))); \
  W("  v_p(A_n):"); \
  for(j=1,4, p=[2,3,5,7][j]; s=""; for(i=1,#IDX, n=IDX[i]; if(n<=M-1 && An[n+1]!=0, s=Str(s," ",valuation(An[n+1],p)))); W(Str("    p=",p,":",s))); \
  W("  SLOPE v_p(B_n/A_n - B_{n-1}/A_{n-1}):"); \
  for(j=1,4, p=[2,3,5,7][j]; s=""; for(i=1,#IDX, n=IDX[i]; if(n<=M-1 && n>1 && An[n+1]!=0 && An[n]!=0, s=Str(s," ",valuation(Bn[n+1]/An[n+1]-Bn[n]/An[n],p)))); W(Str("    p=",p,":",s))); \
  W(Str("  B_n/A_n at n=",M-1,": ", Bn[M]/An[M]*1.0)); \
  W(Str("  (121/729) zeta(5) = ", 121/729*zeta(5)*1.0)); 0; };

{ if(#found > 0, my(r = found[1][1], hh = q^(sum(i=1,#DV,DV[i]*r[i])/24)*prod(i=1,#DV, ee(DV[i])^r[i]) + O(q^(M+1)));
  W(Str("  chosen Hauptmodul q-expansion: ", truncate(hh + O(q^8))));
  doit("plain Gamma_0(18) Hauptmodul u", hh);
  doit("Fricke-type quotient u/(1+a u)^2 with a=1", hh/(1+hh)^2);
  doit("quotient u/(1+3u+ ... ) heuristic: u/(1+u)^2 scaled", hh/(1-hh)^2)); }
W(""); W("DONE"); quit;
