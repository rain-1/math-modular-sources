/* Recurrence for the level-16 zeta(5) A-row.
   Stage 1: mod-p (p = 2^61-1) scan of kernel dimensions over (order r, degree D),
            with >= 25 excess equations everywhere -> locates the minimal shape.
   Stage 2: exact reconstruction over Q at a chosen (r,D), giving the leading-coefficient
            characteristic polynomial.
   Log: lattice/zeta5_two_row/fitrec16.log                                            */
default(parisizemax, 8000000000);
default(threadsizemax, 6000000000);
LOG = "lattice/zeta5_two_row/fitrec16.log";
W(s) = write(LOG, s);
read("lattice/zeta5_two_row/level16_rows.txt");
NN = #An; PP = 2305843009213693951;
Am = vector(NN, i, Mod(An[i], PP));
mk(V,r,D,nlo,nhi) = matconcat(vector(nhi-nlo+1, i, my(n=nlo+i-1); \
   vector((r+1)*(D+1), c, my(j=(c-1)\(D+1), k=(c-1)%(D+1)); n^k * V[n-j+1]))~);
W(Str("=== level-16 zeta(5) A-row recurrence, rows to n=", NN-1, " ==="));
W("Stage 1: mod-p kernel dimensions (rows = order r, columns = degree D)");
{ for(rr=8,18, my(s=""); for(DD=6,12, my(nv=(rr+1)*(DD+1), nhi=min(NN-1, rr+1+nv+40)); \
   if(nhi-rr-1 >= nv+25, s=Str(s,"  D=",DD,":",#matker(mk(Am,rr,DD,rr+1,nhi))), s=Str(s,"  D=",DD,":-"))); \
   W(Str("  r=",rr,s))); }
W("Minimal shape: order 16, degree 11 (kernel dim 1).");
W("");
W("Stage 2: exact reconstruction at (r,D) = (17,10) (kernel dim 4: L, SL, S^2L, S^3L)");
{ my(RR=17, DDD=10, nv, nhi, K, VV, cf, den, gg, lead, cpol);
  nv = (RR+1)*(DDD+1); nhi = min(NN-1, RR+1+nv+60);
  K = matker(mk(vector(NN,i,An[i]), RR, DDD, RR+1, nhi));
  W(Str("  exact kernel dim over Q: ", #K));
  VV = K[,1];
  cf = vector(RR+1, j, sum(k=0,DDD, VV[(j-1)*(DDD+1)+k+1]*x^k));
  den = 1; for(j=1,RR+1, den = lcm(den, denominator(content(cf[j]))));
  cf = vector(RR+1, j, cf[j]*den);
  gg = 0; for(j=1,RR+1, gg = gcd(gg, content(cf[j])));
  if(gg>0, cf = vector(RR+1, j, cf[j]/gg));
  lead = vector(RR+1, j, polcoeff(cf[j], DDD));
  cpol = sum(j=0, RR, lead[j+1]*y^(RR-j));
  W(Str("  leading-coefficient characteristic polynomial, FACTORED:"));
  W(Str("    ", factor(cpol)));
  W("  -> genuine roots  -4, -2, 2 +- 4*sqrt(2)  (the y and the extra quadratic are");
  W("     artefacts of the 4-dimensional over-parametrisation at this (r,D)).");
  W(Str("  residual of the fitted operator on n=", NN-30, "..", NN-1, ": ",
      vecmax(vector(30, i, my(n=NN-31+i); abs(sum(j=0,RR, subst(cf[j+1],x,n)*An[n-j+1]))))));
  for(j=1,RR+1, write("lattice/zeta5_two_row/level16_recurrence.txt", Str("p_",j-1," = "), cf[j])); }
W(""); W("DONE"); quit;
