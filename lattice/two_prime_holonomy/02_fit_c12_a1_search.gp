default(parisizemax, 8000000000);
\\ 02_fit_c12_a1_search.gp -- hunt for the minimal joint (Q,P) operator of the
\\ conductor-12 product form at alpha = 1, c12rowB(b,b,b), using b <= 340.
\\ Coarse scan first (step 10), then refine downwards.
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/02_fit_recur.gp");
read(concat(outdir, "c12_a1_long.seq"));
SQ = [QV, PV];  NA = #QV - 1;
nul(r, d, nl) = #matker(buildmat(SQ, r, d, nl, nextprime(2^61)));
{
foreach([3, 4], r,
  my(nl = vector(NA - 10 - r, k, k), dmax = floor((2*#nl - 20)/(r+1)) - 1, hit = 0);
  print("--- r=", r, "  rows=", 2*#nl, "  dmax allowed = ", dmax);
  forstep(d = 40, dmax, 10,
    my(nu = nul(r, d, nl));
    print("    d=", d, " nullity=", nu);
    if(nu > 0, hit = d; break));
  if(hit,
    for(k = 1, 9,
      my(d = hit - k, nu = nul(r, d, nl));
      print("    refine d=", d, " nullity=", nu);
      if(nu == 0, print("    => MINIMAL for r=", r, " is d=", d+1); break)),
    print("    no hit up to d=", dmax)));
}
quit;
