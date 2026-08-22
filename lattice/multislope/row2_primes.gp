default(parisizemax, 6000000000);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row2_primes.log";
W(s) = write(LOG, s);
th='th;
P0 =  th^4;
P1 = -2^4*(1072*th^4 - 17824*th^3 - 10888*th^2 - 1976*th - 145);
P2 = -2^17*(51088*th^4 + 116368*th^3 - 45264*th^2 - 14228*th - 1397);
P3 =  2^28*13*(73104*th^4 + 1536*th^3 - 488*th^2 + 384*th + 97);
P4 = -2^44*13^2*(2*th+1)^4;
PL = [P0,P1,P2,P3,P4];
Qp = vector(5, i, subst(PL[i], th, 'n-(i-1)));
QC = vector(5, i, Vecrev(Qp[i]));
evq(i,n) = { my(v=QC[i+1], s=0, t=1); for(j=1,#v, s+=v[j]*t; t*=n); s };
NB = 600;
run(seed) = { my(cur=vector(NB+1), sm); cur[seed+1]=(seed!)^4;
  for(n=seed+1, NB, sm=0;
    if(cur[n]  !=0, sm += evq(1,n)*cur[n]);
    if(n>=2, if(cur[n-1]!=0, sm += evq(2,n)*(n-1)^4*cur[n-1]));
    if(n>=3, if(cur[n-2]!=0, sm += evq(3,n)*((n-1)*(n-2))^4*cur[n-2]));
    if(n>=4, if(cur[n-3]!=0, sm += evq(4,n)*((n-1)*(n-2)*(n-3))^4*cur[n-3]));
    cur[n+1] = -sm); cur};
WA = run(0); WB = run(1); WC = run(2); WD = run(3);
W(Str("=== AESZ207: slope scan over all primes p < 60, N = ", NB, " ==="));
W("  reported: v_p(r_n - r_{n-1}) at n = 200, 400, 600  (a slope exists only if these grow ~ linearly)");
{ my(nm=["B","C","D"], ws=[WB,WC,WD]);
  forprime(p=2, 60,
    my(line = Str("  p=", p, ": "));
    for(s=1,3,
      my(str="");
      for(k=1,3, my(n=[200,400,600][k], d = ws[s][n+1]/WA[n+1] - ws[s][n]/WA[n]);
        str = Str(str, if(d==0,"Z",valuation(d,p)), if(k<3,",","")));
      line = Str(line, nm[s], "[", str, "]  "));
    W(line)); }
W("DONE"); quit;
