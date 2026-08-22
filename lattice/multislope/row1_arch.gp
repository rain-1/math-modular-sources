default(parisizemax, 6000000000);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_arch.log";
W(s) = write(LOG, s);
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_rec.txt");
Qp = QROW1; RR = 6;
QC = vector(RR+1, i, Vecrev(Qp[i]));
ev(i,n) = { my(v=QC[i+1], s=0, t=1); for(j=1,#v, s+=v[j]*t; t*=n); s };
NB = if(type(NBIG)=="t_INT", NBIG, 3000);
run(seed) = {
  my(cur = vector(NB+1), sm, ff);
  cur[seed+1] = (seed!)^4;
  for(n = seed+1, NB,
    sm = 0; ff = 1;
    for(i = 1, min(RR,n),
      if(i>1, ff *= (n-i+1)^4);
      if(cur[n-i+1] != 0, sm += ev(i,n)*ff*cur[n-i+1]));
    cur[n+1] = -sm);
  cur};
W(Str("=== Row 1 (Sym^3 Zagier E, t-coordinate) archimedean, N = ", NB, " ==="));
gettime(); WA = run(0);
W(Str("A built in ", gettime(), " ms; A_1..A_5 = ", vector(5,i,WA[i+1]/(i!)^4)));
default(realprecision, 400);
/* asymptotic exponent of A_n:  A_n ~ C 8^n n^alpha  */
W("");
W("--- diagnosing A_n ~ C * 8^n * n^alpha ---");
for(i=1,4, my(n=[500,1000,2000,NB][i], u, v); u = (WA[n+1]/(n!)^4)*1.0/8^n; v = (WA[n]/((n-1)!)^4)*1.0/8^(n-1); W(Str("  n=",n,"  A_n/8^n = ",u,"   log ratio/log(n/(n-1)) = ", log(u/v)/log(n/(n-1.0)))));
CN = ["B","C","D","E","F"];
RT = vector(5);
for(s=1,5, my(wx = run(s)); RT[s] = vector(NB+1, i, if(WA[i]!=0, (wx[i]*1.0)/WA[i], 0)); wx=0);
W("");
W("--- ratios r_n = x_n/a_n and increments ---");
for(s=1,5, W(Str("  ", CN[s], ":  r_N = ", RT[s][NB+1])); \
  for(i=1,3, my(n=[NB, NB-500, NB-1000][i]); W(Str("      n=",n," r_n=",RT[s][n+1], "   n^2*(r_n-r_{n-1}) = ", n^2*(RT[s][n+1]-RT[s][n])))));
/* Richardson: fit r_n = L + sum_{j=1}^{J} c_j/n^j using J+1 points */
rich(rv, J, N) = {
  my(mt, rhs, idx);
  idx = vector(J+1, i, N - 20*(i-1));
  mt = matrix(J+1, J+1, a, b, if(b==1, 1.0, 1.0/idx[a]^(b-1)));
  rhs = vector(J+1, a, rv[idx[a]+1])~;
  (matsolve(mt, rhs))[1] };
W("");
W("--- Richardson extrapolation (fit L + sum_{j<=J} c_j n^-j) ---");
for(s=1,5, my(str=Str("  ",CN[s],":")); for(J=2,12,2, str=Str(str,"\n      J=",J,": ", rich(RT[s], J, NB))); W(str));
/* targets */
default(realprecision, 60);
bet4 = 4^(-4)*(zetahurwitz(4,1/4) - zetahurwitz(4,3/4));
W("");
W(Str("beta(4) = ", bet4, "  (expect 0.98894455174110533...)"));
GG = Catalan;
tg = [1.0, GG, bet4, Pi^4, Pi^3, zeta(3), Pi*zeta(3), GG*Pi^2];
tn = ["1","G","beta4","pi^4","pi^3","zeta3","pi*zeta3","G*pi^2"];
W("");
W("--- lindep of each extrapolated limit against {1,G,beta4,pi^4,pi^3,zeta3,pi*zeta3,G*pi^2} (40 digits) ---");
default(realprecision, 40);
for(s=1,5, my(LL = rich(RT[s], 10, NB)); \
  W(Str("  ", CN[s], "  L ~ ", LL)); \
  W(Str("      lindep([L,1,G,beta4]) = ", lindep([LL, 1.0, GG, bet4]))); \
  W(Str("      lindep([L,1,G])       = ", lindep([LL, 1.0, GG]))); \
  W(Str("      lindep([L,1])         = ", lindep([LL, 1.0]))); \
  W(Str("      bestappr(L, 10^12)    = ", bestappr(LL, 10^12))));
W("DONE");
quit;
