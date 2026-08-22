default(parisizemax, 6000000000);
default(threadsizemax, 2000000000);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row3_slopes.log";
W(s) = write(LOG, s);
read("/home/ubuntu/code/math-modular-sources/lattice/zeta5_two_row/level16_rows.txt");
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row3_rec18.txt");
QQ = QROW3B; RR = 18;
QC = vector(RR+1, i, Vecrev(QQ[i]));
evq(i,n) = { my(v=QC[i+1], s=0, t=1); for(j=1,#v, s+=v[j]*t; t*=n); s };
NB = 400;
Q0 = vector(NB+1, i, evq(0, i-1));    /* Q0[n+1] = Q_0(n) */
PI = vector(NB+1); PI[1] = 1;
{ for(n=1, NB, PI[n+1] = PI[n]*Q0[n+1]); }
/* w_n = PI_n u_n ;  w_n = -sum_{i=1}^{RR} Q_i(n) * (PI_{n-1}/PI_{n-i}) * w_{n-i} */
run(seed) = { my(cur = vector(NB+1), sm, pr);
  cur[seed+1] = PI[seed+1];
  for(n = seed+1, NB, sm = 0; pr = 1;
    for(i = 1, min(RR,n),
      if(i>1, pr *= Q0[n-i+2]);
      if(cur[n-i+1] != 0, sm += evq(i,n)*pr*cur[n-i+1]));
    cur[n+1] = -sm);
  cur};
gettime();
WA = run(0);
W(Str("=== level-16 zeta(5): companions of the exact order-18 operator, N=", NB, " ==="));
W(Str("  A built ", gettime(), " ms"));
W(Str("  check A-companion(seed 0) reproduces An: ", vector(8, i, WA[i]/PI[i]) == vector(8,i,An[i])));
CN = ["X1","X2","X3"];
WS = vector(3, s, run(s));
{ for(s=1,3, W(Str("  ", CN[s], "_0..7 = ", vector(8, i, WS[s][i]/PI[i])))); }
W(Str("  B_0..7 (modular)  = ", vector(8,i,Bn[i])));
/* is X1 proportional to Bn ? */
{ my(rat = vector(10, i, if(WS[1][i+1]!=0 && Bn[i+1]!=0, (WS[1][i+1]/PI[i+1])/Bn[i+1], 0)));
  W(Str("  X1_n / B_n for n=0..9: ", rat)); }
dn = vector(NB+1); dn[1]=1; { for(n=1,NB, dn[n+1]=lcm(dn[n],n)); }
W("");
W("--- sharp lcm-denominator exponent k (n <= 150) ---");
{ my(kf); for(k=0,12, my(ok=1); for(n=1,150, if(denominator(dn[n+1]^k*Bn[n+1])!=1, ok=0;break)); if(ok, kf=k; break));
  W(Str("  modular B_n: k = ", kf));
  for(s=1,3, kf=-1; for(k=0,14, my(ok=1); for(n=1,150, if(denominator(dn[n+1]^k*WS[s][n+1]/PI[n+1])!=1, ok=0;break)); if(ok, kf=k; break));
    W(Str("  ", CN[s], ": k = ", kf))); }
W("");
W("--- slopes v_p(r_n - r_{n-1}), p = 2,3,5,7 ---");
SM = [50,100,200,300,400];
{ for(j=1,4, my(p=[2,3,5,7][j]);
    my(line = Str("  p=", p, "   modular B: "));
    for(i=1,#SM, my(n=SM[i]); if(n<=399, my(d=Bn[n+1]/An[n+1]-Bn[n]/An[n]); line=Str(line, if(d==0,"Z",valuation(d,p)), " ")));
    W(line);
    for(s=1,3, my(l2 = Str("        ", CN[s], ": "));
      for(i=1,#SM, my(n=SM[i], d = WS[s][n+1]/WA[n+1] - WS[s][n]/WA[n]); l2 = Str(l2, if(d==0,"Z",valuation(d,p)), " "));
      W(l2))); }
W("");
W("--- v_2(xi) and Cauchy precision at n = 400 ---");
{ for(s=1,3, my(r1 = WS[s][NB+1]/WA[NB+1], r0 = WS[s][NB]/WA[NB]);
    W(Str("  ", CN[s], "  v_2(r_N) = ", valuation(r1,2), "   v_2(r_N - r_{N-1}) = ", valuation(r1-r0,2)))); }
write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row3_data.txt", Str("R3A = ", WA[NB+1]));
{ for(s=1,3, write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row3_data.txt", Str("R3X",s," = ", WS[s][NB+1]))); }
W("DONE"); quit;
