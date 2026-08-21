/* Brown-Zudilin zeta(5) cellular row, arXiv:2210.03391v3 Sect.2 (totally symmetric case).
   Exact rational construction of Q_n, Phat_n, P_n from the published third-order
   recursion + initial data; arithmetic / archimedean / p-adic diagnostics.
   Log: lattice/zeta5_two_row/bz_row.log                                            */
default(parisizemax, 8000000000);
\p 6000
LOG = "lattice/zeta5_two_row/bz_row.log";
W(s) = write(LOG, s);
if(type(BZN)!="t_INT", BZN = 600);
N = BZN;

c0(n) = 2*(2*n+1)*(41218*n^3-48459*n^2+20010*n-2871)*(n+1)^5;
c1(n) = -(97604224*n^9 + 178061760*n^8 + 72005308*n^7 - 48634688*n^6 - 39076836*n^5 + 2622730*n^4 + 7581006*n^3 + 920112*n^2 - 543402*n - 120582);
c2(n) = -2*n*(3874492*n^8 - 2617900*n^7 - 3144314*n^6 + 2947148*n^5 + 647130*n^4 - 1182926*n^3 + 115771*n^2 + 170716*n - 44541);
c3(n) = n*(41218*n^3+75195*n^2+46746*n+9898)*(n-1)^5;

run(v0,v1,v2) = { my(v=vector(N+1)); v[1]=v0; v[2]=v1; v[3]=v2;
  for(n=2, N-1, v[n+2] = -(c1(n)*v[n+1] + c2(n)*v[n] + c3(n)*v[n-1])/c0(n) ); v; };

Q  = run(1, 21, 2989);
Ph = run(0, 101/4, 344923/96);
P  = run(0, 87/4, 1190161/384);

W(Str("=== Brown-Zudilin zeta(5) row, N=", N, " ==="));

/* closed form (7), as transcribed and validated */
Qcf(n) = sum(k1=0,n, sum(k2=0,n, binomial(n+k1,n)*binomial(n,k1)^2*binomial(n+k2,n)*binomial(n,k2)^2*binomial(n+k1+k2,n)));
{ my(ok=1); for(n=0,7, if(Q[n+1]!=Qcf(n), ok=0));
  W(Str("closed form (7) [sum binom(n+k1,n)binom(n,k1)^2 * (k2) * binom(n+k1+k2,n)] matches recurrence for n<=7: ", ok)); }
W(Str("Q_0..Q_5 = ", vector(6,i,Q[i])));

d(n) = if(n<1, 1, lcm(vector(n,i,i)));
{ my(badQ=0, badP=0, badPh=0, bad12=0);
  for(n=0,N,
    if(denominator(Q[n+1])!=1, badQ++);
    if(denominator(d(n)^5*P[n+1])!=1, badP++);
    if(denominator(12*d(n)^5*P[n+1])!=1, bad12++);
    if(denominator(d(2*n)^2*Ph[n+1])!=1, badPh++));
  W(Str("# n<=",N," with Q_n not in Z          : ", badQ));
  W(Str("# n<=",N," with d_n^5 P_n not in Z    : ", badP));
  W(Str("# n<=",N," with 12 d_n^5 P_n not in Z : ", bad12));
  W(Str("# n<=",N," with d_2n^2 Phat_n not in Z: ", badPh)); }

W("");
W("--- exact 2-adic and 3-adic content of denom(P_n) vs 5*floor(log2 n) ---");
{ for(n=3,20, W(Str("  n=",n," v2(denom P_n)=",valuation(denominator(P[n+1]),2),
        "  5*floor(log2 n)=",5*floor(log(n)/log(2)+1e-9),
        "  v3(denom)=",valuation(denominator(P[n+1]),3))) ); }

W("");
W("--- KAPPA test: is v_2(Q_n) or v_2(P_n) linear in n?  (kappa_p rate) ---");
{ for(i=1,8, my(n=[10,50,100,200,300,400,500,N-1][i]); if(n<=N,
    W(Str("  n=",n,"  v2(Q_n)=",valuation(Q[n+1],2),
          "  v3(Q_n)=",valuation(Q[n+1],3),
          "  v5(Q_n)=",valuation(Q[n+1],5),
          "  v7(Q_n)=",valuation(Q[n+1],7),
          "  v2(P_n)=",valuation(P[n+1],2),
          "  v2(Ph_n)=",valuation(Ph[n+1],2))) )); }

W("");
W("--- archimedean (BZ predict: log|Q|/n -> 6.38364071, log|I'|/n -> -2.47237372, log|I|/n -> -5.29756135) ---");
z5 = zeta(5); z3 = zeta(3); z2 = zeta(2);
{ for(i=1,8, my(n=[10,50,100,200,300,400,500,N-1][i]); if(n<=N && n>2,
    my(Ip = Q[n+1]*z5-P[n+1], Ipp = Q[n+1]*z3-Ph[n+1], Ifull = 2*Ip+4*Ipp*z2);
    W(Str("  n=",n," log|Q|/n=",log(abs(Q[n+1]*1.0))/n,
          " log|Qz5-P|/n=",log(abs(Ip))/n,
          " log|Qz3-Ph|/n=",log(abs(Ipp))/n,
          " log|I|/n=",log(abs(Ifull))/n)) )); }

W("");
PL = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,557];
W("--- SLOPE TEST: v_p( P_n/Q_n - P_{n-1}/Q_{n-1} ) ---");
{ for(j=1,#PL, my(p=PL[j], s="");
    for(i=1,8, my(n=[10,50,100,200,300,400,500,N-1][i]); if(n<=N && n>0,
      s=Str(s," ",valuation(P[n+1]/Q[n+1]-P[n]/Q[n],p))));
    W(Str("  p=",p,":",s))); }
W("--- same for Phat/Q (the zeta(3) direction) ---");
{ for(j=1,4, my(p=[2,3,5,7][j], s="");
    for(i=1,8, my(n=[10,50,100,200,300,400,500,N-1][i]); if(n<=N && n>0,
      s=Str(s," ",valuation(Ph[n+1]/Q[n+1]-Ph[n]/Q[n],p))));
    W(Str("  p=",p,":",s))); }

W("");
W("--- v_p( 2x2 minor  Q_n P_{n-1} - Q_{n-1} P_n ) ---");
{ for(j=1,4, my(p=[2,3,5,7][j], s="");
    for(i=1,8, my(n=[10,50,100,200,300,400,500,N-1][i]); if(n<=N && n>0,
      s=Str(s," ",valuation(Q[n+1]*P[n]-Q[n]*P[n+1],p))));
    W(Str("  p=",p,":",s))); }

W("");
W("--- v_p( 3x3 Casoratian det[Q,Ph,P] at (n,n-1,n-2) ); predicted product-of-roots -1/4 => v_2 ~ 2n ---");
Cas(n) = matdet([Q[n+1],Ph[n+1],P[n+1]; Q[n],Ph[n],P[n]; Q[n-1],Ph[n-1],P[n-1]]);
{ for(j=1,4, my(p=[2,3,5,7][j], s="");
    for(i=1,8, my(n=[10,50,100,200,300,400,500,N-1][i]); if(n<=N && n>1,
      s=Str(s," ",valuation(Cas(n),p))));
    W(Str("  p=",p,":",s))); }

W("");
W("--- fine v_2 increments of P_n/Q_n, n=N-30..N-1 ---");
{ my(s=""); for(n=max(3,N-30),N-1, s=Str(s," ",valuation(P[n+1]/Q[n+1]-P[n]/Q[n],2))); W(s); }
W("--- fine v_2(Q_n), n=N-30..N-1 ---");
{ my(s=""); for(n=max(3,N-30),N-1, s=Str(s," ",valuation(Q[n+1],2))); W(s); }

W("");
W("--- 2-adic tower n=a*2^e : v_2(P_n/Q_n - P_m/Q_m) between consecutive tower members ---");
{ for(j=1,4, my(a=[1,3,5,7][j], s="", pn=0);
    for(e=0,12, my(n=a*2^e); if(n<=N-1 && n>0,
      if(pn>0, s=Str(s,"  (",pn,"->",n,"):",valuation(P[n+1]/Q[n+1]-P[pn+1]/Q[pn+1],2)));
      pn=n));
    W(Str("  a=",a,":",s))); }

W("");
W("--- 2-adic Newton polygon of char poly 4L^3-2368L^2-188L+1 ---");
W(Str("  v2 of coeffs (const..lead): ",[valuation(1,2),valuation(188,2),valuation(2368,2),valuation(4,2)]));
W(Str("  factorisation over Q_2 (2-adic, prec 20): ", factorpadic(4*x^3-2368*x^2-188*x+1, 2, 20)));
W(Str("  factorisation over Q_3: ", factorpadic(4*x^3-2368*x^2-188*x+1, 3, 10)));

W("");
W("DONE");
quit;
