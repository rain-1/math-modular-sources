/* Fit a polynomial-coefficient recurrence to the level-16 zeta(5) A-row.
   Strategy: locate (order r, degree D) modulo a 61-bit prime (fast), then
   reconstruct the coefficients exactly over Q on the smallest consistent system.
   Log: lattice/zeta5_two_row/fitrec16.log */
default(parisizemax, 6000000000);
default(threadsizemax, 4000000000);
LOG = "lattice/zeta5_two_row/fitrec16.log";
W(s) = write(LOG, s);
read("lattice/zeta5_two_row/level16_rows.txt");
NN = #An;
W(Str("=== recurrence fit, level-16 zeta(5) A-row, rows to n=", NN-1, " ==="));
PP = 2305843009213693951;   /* 2^61-1 */
Am = vector(NN, i, Mod(An[i], PP));

mk(V, r, D, nlo, nhi) = matconcat(vector(nhi-nlo+1, i, my(n=nlo+i-1); \
   vector((r+1)*(D+1), c, my(j=(c-1)\(D+1), k=(c-1)%(D+1)); n^k * V[n-j+1]))~);

FOUND = 0;
{ for(rr = 2, 30, if(FOUND, next); \
   for(DD = 1, 10, if(FOUND, next); \
    nv = (rr+1)*(DD+1); \
    nhi = min(NN-1, rr + 1 + nv + 40); \
    if(nhi - rr - 1 < nv + 25, next); \
    KK = matker(mk(Am, rr, DD, rr+1, nhi)); \
    if(#KK >= 1, FOUND = 1; RR = rr; DDD = DD; \
      W(Str("  mod-p hit: order r=", rr, ", degree D=", DD, ", unknowns=", nv,
            ", equations=", nhi-rr, ", kernel dim=", #KK))))); }

{ if(!FOUND, W("  no fit with order<=30, degree<=10"),
  nv = (RR+1)*(DDD+1);
  nhi = min(NN-1, RR + 1 + nv + 60);
  K = matker(mk(vector(NN,i,An[i]), RR, DDD, RR+1, nhi));
  W(Str("  exact kernel dim over Q: ", #K));
  if(#K >= 1,
    VV = K[,1];
    cf = vector(RR+1, j, sum(k=0,DDD, VV[(j-1)*(DDD+1)+k+1]*x^k));
    den = 1; for(j=1,RR+1, den = lcm(den, denominator(content(cf[j]))));
    cf = vector(RR+1, j, cf[j]*den);
    gg = 0; for(j=1,RR+1, gg = gcd(gg, content(cf[j])));
    if(gg>0, cf = vector(RR+1, j, cf[j]/gg));
    W(Str("  p_0..p_r (x = n):"));
    for(j=1,RR+1, W(Str("    p_",j-1,"(x) = ", cf[j])));
    lead = vector(RR+1, j, polcoeff(cf[j], DDD));
    cpol = sum(j=0, RR, lead[j+1]*y^(RR-j));
    W(Str("  leading coefficients: ", lead));
    W(Str("  characteristic polynomial: ", cpol));
    W(Str("  factored: ", factor(cpol)));
    W(Str("  roots: ", polroots(cpol)));
    if(lead[1]!=0,
      cc = (-1)^RR*lead[RR+1]/lead[1];
      W(Str("  product of roots c = (-1)^r p_r/p_0 = ", cc, "   v_2(c) = ", valuation(cc,2))));
    W(Str("  residual on n=", NN-30, "..", NN-1, ": ",
      vecmax(vector(30, i, my(n=NN-31+i); abs(sum(j=0,RR, subst(cf[j+1],x,n)*An[n-j+1]))))));
    W(Str("  same operator on B-row (inhomogeneity), n=", NN-5, "..", NN-1, ": ",
      vector(5, i, my(n=NN-6+i); sum(j=0,RR, subst(cf[j+1],x,n)*Bn[n-j+1])))))); }
W(""); W("DONE"); quit;
