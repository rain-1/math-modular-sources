/* AESZ 207: A, B, C, E (companions seeded at 1,2,3); p-adic analysis */
default(parisizemax, 6000000000);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row2_padic.log";
W(s) = write(LOG, s);
th = 'th;
P0 =  th^4;
P1 = -2^4  * (1072*th^4 - 17824*th^3 - 10888*th^2 - 1976*th - 145);
P2 = -2^17 * (51088*th^4 + 116368*th^3 - 45264*th^2 - 14228*th - 1397);
P3 =  2^28*13 * (73104*th^4 + 1536*th^3 - 488*th^2 + 384*th + 97);
P4 = -2^44*13^2 * (2*th+1)^4;
PL = [P0,P1,P2,P3,P4];
/* recurrence: sum_{i=0}^{4} P_i(n-i) u_{n-i} = 0, i.e. Q_i(n) := P_i(n-i) */
Qp = vector(5, i, subst(PL[i], th, 'n - (i-1)));
RR = 4;
qc(i,n) = subst(Qp[i+1], 'n, n);
W("=== AESZ 207 ===");
for(i=0,4, W(Str("  Q_",i,"(n) = P_",i,"(n-",i,") = ", Qp[i+1])));
cp = sum(i=0,4, polcoeff(Qp[i+1],4,'n)*'x^(4-i));
W(Str("char poly = ", cp));
W(Str("  factored: ", factor(cp)));
default(realprecision, 40);
W(Str("  roots: ", polroots(cp)));
W(Str("  product of roots = ", polcoeff(cp,0,'x)));
W(Str("  v_2 of product = ", valuation(polcoeff(cp,0,'x),2), "   13-part: ", valuation(polcoeff(cp,0,'x),13)));
/* Newton polygon at p=2,13,17 : need v_p of each root */
W("");
W("--- Newton polygon of char poly ---");
for(j=1,3, my(p=[2,13,17][j]); \
  W(Str(" p=",p,"  Newton slopes: ", newtonpoly(cp, p))));
W(Str(" 53248 = 2^12*13 : v_2=",valuation(53248,2)," v_13=",valuation(53248,13)));
q17 = 'x^2 + 89344*'x - 2^24;
W(Str(" quadratic factor x^2+89344x-2^24, disc = ", 89344^2+4*2^24, " = ", factor(89344^2+4*2^24)));
s17 = sqrt(17 + O(2^60));
W(Str(" sqrt17 in Q_2 (17 = 1 mod 8): ", s17));
r1 = -2^7*(349 + 85*s17); r2 = -2^7*(349 - 85*s17);
W(Str(" v_2(-2^7(349+85 sqrt17)) = ", valuation(r1,2), "    v_2(-2^7(349-85 sqrt17)) = ", valuation(r2,2)));
W(Str(" so v_2 of the 4 roots = {", valuation(r1,2), ", ", valuation(r2,2), ", 12, 12}"));

/* --- sequences, scaled w_n = (n!)^4 u_n --- */
NB = if(type(NBIG)=="t_INT", NBIG, 1600);
buildseq(seed) = {my(w = vector(NB+1)); w[seed+1] = ((seed)!)^4;
  for(n=seed+1, NB, my(s=0);
    for(i=1, min(RR,n), if(w[n-i+1]!=0, s += qc(i,n) * ((((n-1)!)/((n-i)!))^4) * w[n-i+1]));
    w[n+1] = -s); w};
fac4 = vector(NB+1, i, ((i-1)!)^4);
WA = buildseq(0);
CN = ["B","C","D"];
WX = vector(3, s, buildseq(s));
W("");
W(Str("A_0..A_5 = ", vector(6,i,WA[i]/fac4[i])));
for(s=1,3, W(Str(CN[s],"_0..: ", vector(7,i,WX[s][i]/fac4[i]))));
bad=0; for(i=1,NB+1, if(denominator(WA[i]/fac4[i])!=1, bad++));
W(Str("A_n integral for n<=",NB,"? ", bad==0));

dn = vector(NB+1); dn[1]=1; for(n=1,NB, dn[n+1]=lcm(dn[n],n));
W("");
W("--- sharp lcm-denominator exponent k (tested n<=200) ---");
for(s=1,3, my(kf=-1); for(k=0,10, my(ok=1); \
    for(n=1,min(200,NB), if(denominator(dn[n+1]^k*WX[s][n+1]/fac4[n+1])!=1, ok=0;break)); \
    if(ok,kf=k;break)); \
  W(Str("  ", CN[s], ": smallest k with d_n^k x_n in Z is k = ", kf)));

W("");
W("--- v_p(A_n) ---");
SAMP = [50,100,200,400,800,1200,NB];
for(j=1,6, my(p=[2,3,5,7,13,17][j], s=""); for(i=1,#SAMP, my(n=SAMP[i]); if(n<=NB, s=Str(s," n=",n,":",valuation(WA[n+1]/fac4[n+1],p)))); W(Str(" p=",p,s)));

W("");
W("--- SLOPE: v_p(x_n/a_n - x_{n-1}/a_{n-1}) ---");
for(j=1,6, my(p=[2,3,5,7,13,17][j]); W(Str(" == p=",p," ==")); \
 for(s=1,3, my(str=""); for(i=1,#SAMP, my(n=SAMP[i]); if(n<=NB && n>1, \
    my(d = WX[s][n+1]/WA[n+1] - WX[s][n]/WA[n]); str=Str(str,"  n=",n,":",if(d==0,"Z",valuation(d,p))))); \
  W(Str("   ",CN[s],str))));

W("");
W("--- 2-adic: v_2(xi) and increments at the top ---");
for(s=1,3, my(rN = WX[s][NB+1]/WA[NB+1], rM = WX[s][NB]/WA[NB]); \
  W(Str("  ",CN[s],"  v_2(x_N/a_N) = ", valuation(rN,2), "   Cauchy precision v_2(r_N - r_{N-1}) = ", valuation(rN-rM,2))));

write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row2_data.txt", Str("NB2 = ", NB));
write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row2_data.txt", Str("R2A = ", WA[NB+1]));
for(s=1,3, write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row2_data.txt", Str("R2X",s," = ", WX[s][NB+1])));
W("DONE");
quit;
