\\ 05_fitprod.gp -- the minimal joint operator of the four Hadamard products,
\\ and the (possibly smaller) operator of W.
default(parisizemax, 12000000000);
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read("/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/lib_fit.gp");
read(concat(outdir, "had.gp"));
N = #HZZ;
PP = nextprime(2^61);
S4 = [HZZ, HZB, HNB, HBB];

print("== order search, joint on the four products, n = 1..", N-8, " ==");
{ for(r = 2, 5,
    my(nl = vector(N-8-r, k, k), res);
    res = hh_mind(S4, r, 0, min(300, (4*(N-8-r)-6)\(r+1) - 1), nl, PP);
    print("  -> r=", r, " minimal d = ", res[1], " nullity ", res[2])); }
\q
