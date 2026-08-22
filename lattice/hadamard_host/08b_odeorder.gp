default(parisizemax, 12000000000);
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read("/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/lib_fit.gp");
read(concat(outdir, "had.gp"));
N = #HZZ; PP = nextprime(2^61);
S4 = [HZZ, HZB, HNB, HBB];
print("== high-order end of the order-degree curve (joint module) ==");
{ for(k = 1, 6, my(r = [24,30,36,44,52,60][k],
    nlr = vector(N-8-r, kk, kk), dmax = (4*(N-8-r)-6)\(r+1) - 1, res);
    dmax = min(dmax, 60);
    res = hh_mind(S4, r, 0, dmax, nlr, PP);
    print("  r=", r, "  min d = ", res[1], "  (nullity ", res[2], ")")); }
\q
