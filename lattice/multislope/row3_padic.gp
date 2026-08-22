default(parisizemax, 6000000000);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row3_padic.log";
W(s) = write(LOG, s);
read("/home/ubuntu/code/math-modular-sources/lattice/zeta5_two_row/level16_rows.txt");
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row3_rec18.txt");
QQ = QROW3B; RR = 18; QC = vector(RR+1, i, Vecrev(QQ[i]));
evq(i,n) = { my(v=QC[i+1], s=0, t=1); for(j=1,#v, s+=v[j]*t; t*=n); s };
NB = 2600; PREC = 30000;
Z2 = 1 + O(2^PREC);
runp(init) = { my(cur = vector(NB+1), sm);
  for(i=1, #init, cur[i] = init[i]*Z2);
  for(n = #init, NB, sm = 0;
    for(i=1, min(RR,n), sm += evq(i,n)*cur[n-i+1]);
    cur[n+1] = -sm/evq(0,n));
  cur};
seedvec(s) = { my(v = vector(s+1, i, 0)); v[s+1] = 1; v };
gettime();
PA = runp(vector(RR+1, i, An[i]));
PB = runp(vector(RR+1, i, Bn[i]));
W(Str("=== level-16 zeta(5): 2-adic run to N=", NB, ", working precision 2^", PREC, "  (", gettime(), " ms) ==="));
NM = ["B(modular)","X1","X2","X3"];
PS = [PB, runp(seedvec(1)), runp(seedvec(2)), runp(seedvec(3))];
W("--- v_2(r_n - r_{n-1}) at sampled n (slope test) ---");
{ for(s=1,4, my(line = Str("  ", NM[s], ": "));
    for(i=1,6, my(n=[400,800,1200,1800,2200,NB][i], d = PS[s][n+1]/PA[n+1] - PS[s][n]/PA[n]);
      line = Str(line, " n=", n, ":", valuation(d,2)));
    W(line)); }
XI = vector(4, s, PS[s][NB+1]/PA[NB+1]);
W("");
{ for(s=1,4, W(Str("  ", NM[s], "  v_2(xi_2) = ", valuation(XI[s],2), "   Cauchy precision = ", valuation(XI[s]-PS[s][NB]/PA[NB], 2)))); }
/* Kubota-Leopoldt zeta_2 */
KTR = 1100; PR2 = 2200;
om(a) = if(a%4==1, 1, -1);
tw(a) = a/om(a);
z2f(s) = { (1/(4*(s-1))) * sum(i=1,2, my(a=[1,3][i]); (tw(a)+O(2^PR2))^(1-s) * sum(j=0, KTR, binomial(1-s,j)*((4/a)+O(2^PR2))^j*bernfrac(j))) };
ZT = vector(6, i, z2f(i+1));
W("");
W(Str("  zeta_2(m) valuations m=2..7: ", vector(6,i,valuation(ZT[i],2))));
W("");
W("--- VALIDATION: is xi_2(B) = (7/32) zeta_2(5)? ---");
{ my(d = XI[1] - (7/32)*ZT[4]); W(Str("  v_2( xi_2(B) - (7/32) zeta_2(5) ) = ", valuation(d,2), "   (large = confirmed)")); }
CUT = 1900;
tst(nm, v) = { my(L = lindep(v), h = vecmax(vector(#L,i,#Str(abs(L[i])))));
  W(Str("  ", nm, " -> max coeff ", h, " digits", if(h<30, Str("   ", L), "  [spurious]"))); L };
W("");
W(Str("--- lindep against zeta_2(m), at ", CUT, " 2-adic digits ---"));
{ for(s=1,4, for(i=1,6, tst(Str(NM[s], " vs zeta_2(", i+1, ")"), [XI[s]+O(2^CUT), ZT[i]+O(2^CUT)]))); }
W("");
W("--- are the four 2-adic limits Q-linearly dependent? ---");
tst("[1,xiB,xiX1,xiX2,xiX3]", [1+O(2^CUT), XI[1]+O(2^CUT), XI[2]+O(2^CUT), XI[3]+O(2^CUT), XI[4]+O(2^CUT)]);
tst("[xiB,xiX1]            ", [XI[1]+O(2^CUT), XI[2]+O(2^CUT)]);
tst("[xiB,xiX2]            ", [XI[1]+O(2^CUT), XI[3]+O(2^CUT)]);
tst("[xiB,xiX3]            ", [XI[1]+O(2^CUT), XI[4]+O(2^CUT)]);
tst("[xiX1,xiX2]           ", [XI[2]+O(2^CUT), XI[3]+O(2^CUT)]);
tst("[1,xiB,xiX1]          ", [1+O(2^CUT), XI[1]+O(2^CUT), XI[2]+O(2^CUT)]);
W("");
W("--- xi_2(X_s) in Q*zeta_2(5) + Q*zeta_2(3) ? and Q + Q*zeta_2(5) ? ---");
{ for(s=2,4,
    tst(Str(NM[s], " [xi, z2(3), z2(5)]"), [XI[s]+O(2^CUT), ZT[2]+O(2^CUT), ZT[4]+O(2^CUT)]);
    tst(Str(NM[s], " [xi, 1, z2(5)]   "), [XI[s]+O(2^CUT), 1+O(2^CUT), ZT[4]+O(2^CUT)])); }
W("DONE"); quit;
