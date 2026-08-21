/* ---------------------------------------------------------------------------
   Cross-row 2-adic alignment test for zeta(5):
     row X = level-16 modular row   (A_m, B_m),  B/A -> (217/1024) zeta(5)
     row Y = Brown-Zudilin cellular row (Q_n, P_n), P/Q -> 1 * zeta(5)
   Cross determinant with the archimedean factors built in (THEORY_NOTES_03 §2):
     Delta_{m,n} = r * A_m * P_n - r' * Q_n * B_m ,   r = 217/1024, r' = 1.
   Aligned  <=>  v_2(Delta) grows linearly along some sampling ratio.
   Log: lattice/zeta5_two_row/align.log
   ------------------------------------------------------------------------- */
default(parisizemax, 8000000000);
\p 100
LOG = "lattice/zeta5_two_row/align.log";
W(s) = write(LOG, s);
if(type(NB)!="t_INT", NB = 200);

c0(n) = 2*(2*n+1)*(41218*n^3-48459*n^2+20010*n-2871)*(n+1)^5;
c1(n) = -(97604224*n^9 + 178061760*n^8 + 72005308*n^7 - 48634688*n^6 - 39076836*n^5 + 2622730*n^4 + 7581006*n^3 + 920112*n^2 - 543402*n - 120582);
c2(n) = -2*n*(3874492*n^8 - 2617900*n^7 - 3144314*n^6 + 2947148*n^5 + 647130*n^4 - 1182926*n^3 + 115771*n^2 + 170716*n - 44541);
c3(n) = n*(41218*n^3+75195*n^2+46746*n+9898)*(n-1)^5;
run(v0,v1,v2) = { my(v=vector(NB+1)); v[1]=v0; v[2]=v1; v[3]=v2; \
  for(n=2, NB-1, v[n+2] = -(c1(n)*v[n+1] + c2(n)*v[n] + c3(n)*v[n-1])/c0(n) ); v; };
Q = run(1,21,2989); Ph = run(0,101/4,344923/96); P = run(0,87/4,1190161/384);

read("lattice/zeta5_two_row/level16_rows.txt");
MX = #An - 1;
W(Str("=== zeta(5) cross-row alignment, BZ n<=", NB-1, ", level-16 m<=", MX, " ==="));
W("  archimedean:  B_m/A_m -> (217/1024) zeta(5);   P_n/Q_n -> zeta(5)");
W("  2-adic slopes: level-16 sigma_2 ~ 1 (measured);  Brown-Zudilin sigma_2 = 0 (measured)");

W("");
W("--- v_2( 217*A_m*P_n - 1024*Q_n*B_m ) on a grid of sampling ratios m = round(rho*n) ---");
{ for(k=1,7, my(rho = [1/4,1/2,1,3/2,2,3,4][k], s="");
    for(i=1,6, my(n = [20,40,60,90,120,NB-1][i], m = round(rho*n));
      if(m>=1 && m<=MX && n>=1, s = Str(s, "  n=",n,",m=",m,":", valuation(217*An[m+1]*P[n+1] - 1024*Q[n+1]*Bn[m+1], 2))));
    W(Str("  rho=m/n=", rho, s))); }

W("");
W("--- control: the same test for the two ALIGNED zeta(3) rows would give ~min(sigma)*n.");
W("--- scan of rational multipliers c in  v_2( c*A_m*P_n - Q_n*B_m ), m=n ---");
{ my(cands = [217/1024, 1, 1024/217, 7/32, 32/7, 217/512, 1/2, 2, 7/8, 1/1024, 1024]);
  for(j=1,#cands, my(cc = cands[j], s="");
    for(i=1,5, my(n = [20,40,60,90,120][i]);
      if(n<=MX && n<=NB-1, s = Str(s, " n=",n,":", valuation(cc*An[n+1]*P[n+1] - Q[n+1]*Bn[n+1], 2))));
    W(Str("  c=", cc, s))); }

W("");
W("--- for reference: v_2 of each factor separately (m=n) ---");
{ for(i=1,5, my(n=[20,40,60,90,120][i]);
    W(Str("  n=", n, "  v2(A_n)=", valuation(An[n+1],2), "  v2(B_n)=", valuation(Bn[n+1],2),
          "  v2(Q_n)=", valuation(Q[n+1],2), "  v2(P_n)=", valuation(P[n+1],2)))); }

W("");
W("--- BZ 2-adic Cauchy test again, dense, n=100..140 : v_2(P_n/Q_n - P_{n-1}/Q_{n-1}) ---");
{ my(s=""); for(n=100,min(140,NB-1), s=Str(s," ",valuation(P[n+1]/Q[n+1]-P[n]/Q[n],2))); W(s); }
W("--- BZ: v_2(Q_n) dense n=100..140 ---");
{ my(s=""); for(n=100,min(140,NB-1), s=Str(s," ",valuation(Q[n+1],2))); W(s); }
W("--- BZ: v_2(numerator of P_n) minus v_2(denominator) dense n=100..140 ---");
{ my(s=""); for(n=100,min(140,NB-1), s=Str(s," ",valuation(P[n+1],2))); W(s); }

W("");
W("DONE"); quit;
