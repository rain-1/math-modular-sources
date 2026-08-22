\\ 06_operator.gp -- exact minimal operators and singular sets on the Hadamard host.
default(parisizemax, 12000000000);
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read("/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/lib_fit.gp");
read(concat(outdir, "had.gp"));
N = #HZZ; PP = nextprime(2^61);
S4 = [HZZ, HZB, HNB, HBB];

\\ ---- W alone: is it annihilated by a smaller operator? ----
print("== W alone: order search ==");
{ for(r = 2, 4,
    my(nl = vector(N-8-r, k, k), res);
    res = hh_mind([WW], r, 0, min(300, ((N-8-r)-6)\(r+1) - 1), nl, PP);
    print("  -> W: r=", r, " minimal d = ", res[1], " nullity ", res[2])); }

\\ ---- exact order-4 joint operator ----
print("");
print("== exact reconstruction of the joint order-4 operator (r=4, d=96) ==");
r = 4; d = 96;
nl = vector(N-8-r, k, k);
ek = hh_ker(S4, r, d, nl, 80);
print("   primes used = ", ek[3]);
cs = hh_clean(hh_tocoefs(ek[1], r, d));
print("   verify (all four products, n=1..", N-r, "): failures = ", hh_verify(S4, cs, r, 1, N-r));
print("   verify (W,            n=1..", N-r, "): failures = ", hh_verify([WW], cs, r, 1, N-r));
hh_sing(cs, r);
write(concat(outdir,"op4.gp"), "OP4 = ", cs, ";");
print("");
print("   c_0 factored = ", factor(cs[1]));
print("   c_4 factored = ", factor(cs[5]));
\q
