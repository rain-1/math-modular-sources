\\ 08_odeorder.gp -- the order-degree curve of the Hadamard module.
\\ A recurrence of order r and degree d corresponds to an ODE in theta=x d/dx of
\\ order d and degree r; so min over r of d = the differential (holonomic) order.
default(parisizemax, 12000000000);
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read("/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/lib_fit.gp");
read(concat(outdir, "had.gp"));
N = #HZZ; PP = nextprime(2^61);
S4 = [HZZ, HZB, HNB, HBB];

print("== order-degree curve, joint module (4 sequences) ==");
{ my(nl = vector(N-8, k, k), best = 10^6);
  for(r = 4, 24,
    my(nlr = vector(N-8-r, k, k), dmax = (4*(N-8-r)-6)\(r+1) - 1, res);
    dmax = min(dmax, 200);
    res = hh_mind(S4, r, 0, dmax, nlr, PP);
    print("  r=", r, "  min d = ", res[1], "  (nullity ", res[2], ")");
    if(res[1] > 0 && res[1] < best, best = res[1]));
  print("  --> minimal differential order over r<=24 : D = ", best); }

print("");
print("== order-degree curve, W alone ==");
{ my(best = 10^6);
  for(r = 4, 24,
    my(nlr = vector(N-8-r, k, k), dmax = ((N-8-r)-6)\(r+1) - 1, res);
    dmax = min(dmax, 200);
    res = hh_mind([WW], r, 0, dmax, nlr, PP);
    print("  r=", r, "  min d = ", res[1], "  (nullity ", res[2], ")");
    if(res[1] > 0 && res[1] < best, best = res[1]));
  print("  --> W: minimal differential order over r<=24 : D = ", best); }
\q
