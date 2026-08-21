/* Fit a polynomial-coefficient recurrence to the level-16 zeta(5) A-row.
   Log: lattice/zeta5_two_row/fitrec16.log */
default(parisizemax, 8000000000);
LOG = "lattice/zeta5_two_row/fitrec16.log";
W(s) = write(LOG, s);
read("lattice/zeta5_two_row/level16_rows.txt");
NN = #An;
W(Str("=== recurrence fit, level-16 zeta(5) A-row, n up to ", NN-1, " ==="));

mkmat(V, r, D, nlo, nhi) = matconcat(vector(nhi-nlo+1, i, my(n=nlo+i-1); \
   vector((r+1)*(D+1), c, my(j=(c-1)\(D+1), k=(c-1)%(D+1)); n^k * V[n-j+1]))~);

{
FOUND = 0;
for(rr = 2, 26, if(FOUND, next); \
 for(DD = 1, 9, if(FOUND, next); \
  nv = (rr+1)*(DD+1); \
  nhi = min(NN-1, rr + 2 + 2*nv + 20); \
  if(nhi - (rr+2) + 1 < nv + 15, next); \
  KK = matker(mkmat(An, rr, DD, rr+2, nhi)); \
  if(#KK == 1, FOUND = 1; RR = rr; DDD = DD; VV = KK[,1]; \
     W(Str("  minimal fit: order r=", rr, ", degree D=", DD, ", equations=", nhi-rr-1)))));

if(FOUND,
  cf = vector(RR+1, j, sum(k=0,DDD, VV[(j-1)*(DDD+1)+k+1]*x^k));
  den = 1; for(j=1,RR+1, den = lcm(den, denominator(content(cf[j]))));
  cf = vector(RR+1, j, cf[j]*den);
  gg = 0; for(j=1,RR+1, gg = gcd(gg, content(cf[j])));
  cf = vector(RR+1, j, cf[j]/gg);
  W(Str("  p_0..p_r (variable x = n): ", cf));
  lead = vector(RR+1, j, polcoeff(cf[j], DDD));
  W(Str("  leading x^", DDD, " coefficients: ", lead));
  cpol = sum(j=0, RR, lead[j+1]*y^(RR-j));
  W(Str("  characteristic polynomial: ", cpol));
  W(Str("  factored: ", factor(cpol)));
  W(Str("  roots: ", polroots(cpol)));
  W(Str("  product of roots  c = (-1)^r p_r/p_0 = ", (-1)^RR*lead[RR+1]/lead[1]));
  W(Str("  v_2(c) = ", if(lead[1]!=0 && lead[RR+1]!=0, valuation((-1)^RR*lead[RR+1]/lead[1],2), "n/a")));
  W(Str("  verify recurrence on n=", NN-40, "..", NN-1, ": max residual = ",
        vecmax(vector(40, i, my(n=NN-41+i); abs(sum(j=0,RR, subst(cf[j+1],x,n)*An[n-j+1]))))));
  W(Str("  same recurrence on B-row residual (inhomogeneity) n=", NN-8, "..", NN-1, ": ",
        vector(8, i, my(n=NN-9+i); sum(j=0,RR, subst(cf[j+1],x,n)*Bn[n-j+1]))));
, W("  NO fit found in the searched range (order<=26, degree<=9)"));
}
W("");
W("DONE"); quit;
