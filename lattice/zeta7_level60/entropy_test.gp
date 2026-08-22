\\ ================================================================
\\ Test of the archive's level-60 denominator-entropy claims.
\\ Sources at level 60 (weight 8, purified, W_60=+1):
\\   B0 = the unique support-depth-3 source (c_1=c_2=0), L(B0,7)=(2623/216000)zeta(7)
\\   Phi5 = (1-5^4 V_5)Phi_12,  c_1=1,       L(Phi5,7)=(6479/54000)zeta(7)
\\ Companion rows built on the SAME integral host t as the level-12 parent
\\ (t is a modular function for Gamma_0(12) hence for Gamma_0(60)):
\\   A(t) = Phi12/(Dt),  B^{(S)}(t) = A(t) * (D^{-7} S)(q(t))
\\ Claim (archive): denominator of B^{(B0)}_n divides const * d_{floor(n/3)}^7,
\\ i.e. entropy 7/3 rather than 7; and a "filtered triple" would give 19/6.
\\ ================================================================
default(parisize,8000000000); default(realprecision,60);
N = 240;
et(k) = eta(q^k + O(q^(N+2)));
dl = divisors(60);
d12=[1,2,3,4,6,12]; c12=[1,-572,11583,-36608,46332,-20736];
B0  =[0,0,2673,-73216,398125,-497664,-3840000,13208832,-14478750,5280000,0,0];
B1  =[0,33,0,-39424,276250,-375921,-2900625,9165312,-7796250,0,1670625,0];
B2  =[1,0,0,-402688,3062500,-4313088,-33280000,101606400,-79633125,0,0,12960000];
Phi5=[1,-572,11583,-36608,-625,46332,357500,-20736,-7239375,22880000,-28957500,12960000];
{coefs(c) = vector(N, m, sum(i=1,#dl, if(m%dl[i]==0, c[i]*sigma(m\dl[i],7), 0)));}
{g12 = vector(N, m, sum(i=1,6, if(m%d12[i]==0, c12[i]*sigma(m\d12[i],7), 0)));}
Phi12 = sum(m=1,N, g12[m]*q^m) + O(q^(N+1));
x = q*et(4)^2*et(12)^2/(et(1)^2*et(3)^2);
t = x/(16*x^2+2*x+1);
Dt = q*deriv(t,q);
Qt = serreverse(t);
A  = subst(Phi12/Dt, q, Qt);
M  = N-10;
An = vector(M+1,i,polcoeff(A,i-1));
{dn = vector(M+1); dn[1]=1; for(n=1,M, dn[n+1]=lcm(dn[n],n));}
{test(nm, c) =
  my(gg=coefs(c), Th, BB, Bn, k, worst, wn, v, Mn);
  Th = sum(m=1,N, (gg[m]/m^7)*q^m) + O(q^(N+1));
  BB = A*subst(Th,q,Qt);
  Bn = vector(M+1,i,polcoeff(BB,i-1));
  print("--- source ", nm, " ---");
  print("  B_1..B_3 = ", vector(3,i,Bn[i+1]));
  k = 0;
  while(k<12, my(ok=1); for(n=1,M, if(denominator(dn[n+1]^k*Bn[n+1])!=1, ok=0;break)); if(ok, break); k++);
  print("  sharp plain exponent: d_n^", k, " B_n in Z  (d_n^",k-1," fails)");
  worst=0; wn=0;
  for(n=1,M, v = denominator(dn[n\3+1]^7*Bn[n+1]); if(v>worst, worst=v; wn=n));
  print("  max denom of d_{floor(n/3)}^7 B_n = ", worst, " at n=", wn, if(worst>1, Str("  = ",factor(worst)), "   ALWAYS INTEGRAL"));
  worst=0; wn=0;
  for(n=1,M, v = denominator(60^7*dn[n\3+1]^7*Bn[n+1]); if(v>worst, worst=v; wn=n));
  print("  max denom of 60^7 d_{floor(n/3)}^7 B_n = ", worst, " at n=", wn, if(worst>1, Str("  = ",factor(worst)), "   ALWAYS INTEGRAL"));
  worst=0; wn=0;
  for(n=1,M, Mn = lcm(lcm(dn[n\3+1]^7, dn[n\2+1]^2), dn[n+1]); v = denominator(60^7*Mn*Bn[n+1]); if(v>worst, worst=v; wn=n));
  print("  max denom of 60^7 lcm(d_{n/3}^7,d_{n/2}^2,d_n) B_n = ", worst, " at n=", wn, if(worst>1, Str("  = ",factor(worst)), "   ALWAYS INTEGRAL"));
  print("  log denom(B_n)/n at n=60,120,180,",M,": ", vector(4,i, my(n=[60,120,180,M][i]); log(1.0*denominator(Bn[n+1]))/n));
  print("     (7 = 7.0,  19/6 = 3.16667,  7/3 = 2.33333)");
  ;}
test("B0 (depth 3)", B0);
test("B1 (depth 2)", B1);
test("B2 (depth 1)", B2);
test("Phi60^(5) = (1-5^4V_5)Phi12", Phi5);
test("Phi12 (level-12 parent, control)", [1,-572,11583,-36608,0,46332,0,-20736,0,0,0,0]);
\q
