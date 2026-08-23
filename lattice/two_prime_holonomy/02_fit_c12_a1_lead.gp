default(parisizemax, 8000000000);
\\ 02_fit_c12_a1_lead.gp -- for c12rowB(b,b,b) the minimal joint operator is
\\ (r,d) = (4,72), nullity 1.  Reconstructing all 365 coefficients over Q is
\\ expensive and unnecessary: only the degrees deg c_i and their leading
\\ coefficients decide the Newton polygon / balance polynomial / singularities.
\\ So we reconstruct just those, multi-modularly.
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/02_fit_recur.gp");
read(concat(outdir, "c12_a1_long.seq"));
SQ = [QV, PV];  NA = #QV - 1;
{
my(r = 4, d = 72, uu = (r+1)*(d+1), nl, p0 = nextprime(2^62), kk, vv,
   dg, lead, acc, modu = 1, cand, prev = 0, pos, nrows = 0);
nl = vector(ceil(uu/2) + 30, k, k);
print("rows = ", 2*#nl, "   unknowns = ", uu);
kk = matker(buildmat(SQ, r, d, nl, p0));
print("nullity at (r,d)=(", r, ",", d, ") = ", #kk);
vv = kk[,1];
\\ degrees, from the mod-p kernel vector
dg = vector(r+1, i, my(e = -1);
   for(j = 0, d, if(lift(vv[(i-1)*(d+1)+j+1]) != 0, e = j)); e);
print("deg c_i = ", dg);
pos = (r+1-1)*(d+1) + dg[r+1] + 1;   \\ leading coefficient of c_r : normalise to 1
acc = vector(r+1);
for(c = 1, 60,
  p0 = nextprime(p0+1);
  kk = matker(buildmat(SQ, r, d, nl, p0));
  if(#kk != 1, next);
  vv = kk[,1] / kk[pos,1];
  for(i = 1, r+1, acc[i] = lift(chinese(Mod(acc[i], modu),
                                        vv[(i-1)*(d+1)+dg[i]+1])));
  modu *= p0;
  cand = vector(r+1, i, bestappr(Mod(acc[i], modu)));
  if(cand == prev, print("stabilised after ", c, " primes"); break);
  prev = cand);
lead = prev;
print("leading coeff of c_i (normalised, c_4 leading = 1) = ", lead);
my(bal = sum(i = 1, r+1, lead[i]*'x^(i-1)));
print("balance polynomial sum_i lead_i rho^i = ", bal);
print("factored = ", factor(bal));
print("roots rho = ", polroots(bal));
print("Newton polygon slopes (deg c_(i+1) - deg c_i) = ",
      vector(r, i, dg[i+1] - dg[i]));
print("=> u_(n+1)/u_n ~ rho * n^", dg[1]-dg[2]);
my(dd = vecmax(dg), pgf = 0);
for(i = 1, r+1, if(dg[i] == dd, pgf += lead[i]*'x^(r-i+1)));
print("Pgf(x) (coefficient of theta^", dd, " in the ODE) = ", pgf);
if(poldegree(pgf,'x) <= 0,
   print("Pgf constant => NO finite nonzero singular point; the generating",
         " functions are ENTIRE."),
   print("finite nonzero singular points = ", polroots(pgf)));
}
quit;
