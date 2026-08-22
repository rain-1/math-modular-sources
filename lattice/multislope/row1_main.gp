default(parisizemax, 6000000000);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_main.log";
W(s) = write(LOG, s);
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_rec.txt");  /* QROW1 */
Qp = QROW1; RR = #Qp - 1;   /* RR = 6 */
qc(i, n) = subst(Qp[i+1], 'n, n);
NB = if(type(NBIG)=="t_INT", NBIG, 800);

/* scaled sequences: w_n = (n!)^4 u_n ; w_n = -sum_{i=1}^{6} Q_i(n) ((n-1)!/(n-i)!)^4 w_{n-i} */
fac4 = vector(NB+1, i, ((i-1)!)^4);
buildseq(seed) = {                 /* seed = index s with u_s = 1, u_j = 0 for j<s */
  my(w = vector(NB+1));
  w[seed+1] = fac4[seed+1];
  for(n = seed+1, NB,
    my(s = 0);
    for(i = 1, min(RR, n),
      if(w[n-i+1] != 0,
         s += qc(i, n) * ((((n-1)!)/((n-i)!))^4) * w[n-i+1]));
    w[n+1] = -s);
  w};

WA = buildseq(0);
W("=== Row 1 sequences (scaled by (n!)^4) ===");
Anv = vector(NB+1, i, WA[i]/fac4[i]);
W(Str("A_0..A_8 = ", vector(9,i,Anv[i])));
bad=0; for(i=1,NB+1, if(denominator(Anv[i])!=1, bad++));
W(Str("A_n integral for n<=", NB, "? ", bad==0));

CN = ["B","C","D","E","F"];
WX = vector(5, s, buildseq(s));
for(s=1,5, my(xv = vector(9, i, WX[s][i]/fac4[i])); W(Str(CN[s], "_0..B_8 = ", xv)));

/* denominator exponent k */
dn = vector(NB+1); dn[1]=1; for(n=1,NB, dn[n+1]=lcm(dn[n],n));
W("");
W("--- sharp lcm-denominator exponent k (n<=200) ---");
for(s=1,5, my(kf = -1); \
  for(k=0, 12, my(ok=1); \
    for(n=1, min(200,NB), my(v = dn[n+1]^k * WX[s][n+1]/fac4[n+1]; ); if(denominator(v)!=1, ok=0; break)); \
    if(ok, kf=k; break)); \
  W(Str("  ", CN[s], ": k = ", kf)));

/* p-adic slopes */
W("");
W("--- v_p(A_n)/n  (unit-root slope mu_1) ---");
for(j=1,4, my(p=[2,3,5,7][j], s=""); for(i=1,6, my(n=[100,200,400,600,700,NB][i]); if(n<=NB, s=Str(s,"  n=",n,": ",valuation(Anv[n+1],p),"(",valuation(Anv[n+1],p)/n*1.0,")"))); W(Str(" p=",p,s)));

W("");
W("--- SLOPE sigma_p: v_p(x_{n}/a_{n} - x_{n-1}/a_{n-1}) ---");
for(j=1,4, my(p=[2,3,5,7][j]); W(Str(" == p = ", p, " ==")); \
 for(s=1,5, my(str=""); \
   for(i=1,7, my(n=[100,200,300,400,500,600,NB][i]); \
     if(n<=NB && Anv[n+1]!=0 && Anv[n]!=0, \
        my(d = WX[s][n+1]/WA[n+1] - WX[s][n]/WA[n]); \
        str = Str(str, "  n=", n, ":", if(d==0,"ZERO",valuation(d,p))))); \
   W(Str("   ", CN[s], str))));

/* fine 2-adic slope fit */
W("");
W("--- fine 2-adic increments, n = NB-10 .. NB ---");
for(s=1,5, my(str=""); for(n=NB-10, NB, my(d = WX[s][n+1]/WA[n+1] - WX[s][n]/WA[n]); str=Str(str," ",if(d==0,"Z",valuation(d,2)))); W(Str("  ",CN[s],":",str)));

/* archimedean */
default(realprecision, 300);
W("");
W("--- archimedean ratios x_n/a_n ---");
for(s=1,5, my(r1 = WX[s][NB+1]/WA[NB+1], r2 = WX[s][NB-49]/WA[NB-49]); \
  W(Str("  ", CN[s], "  n=", NB, ": ", r1*1.0)); \
  W(Str("        n=", NB-50, ": ", r2*1.0)); \
  W(Str("        |diff| = ", abs(r1-r2)*1.0)));

/* save 2-adic + arch data */
write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_data.txt", Str("NB1 = ", NB));
write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_data.txt", Str("WA1 = ", WA));
for(s=1,5, write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_data.txt", Str("WX1", s, " = ", WX[s])));
W("DONE");
quit;
