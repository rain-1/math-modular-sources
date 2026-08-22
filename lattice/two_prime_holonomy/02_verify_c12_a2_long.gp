default(parisizemax, 8000000000);
\\ 02_verify_c12_a2_long.gp -- independent check: refit the minimal joint operator
\\ for c12rowB(2b,2b,b) from the b<=200 data ONLY, then verify it exactly over Q
\\ on the long sequence b <= 340 (i.e. ~150 completely fresh indices).
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/02_fit_recur.gp");
{
my(r = 3, d = 60, seqsS, seqsL, nl, ek, cs, den, gg, NL);
read(concat(outdir, "c12_a2.seq"));      seqsS = [QV, PV];
nl = vector(189, k, k);
ek = exactker(seqsS, r, d, nl, 40);
cs = tocoefs(ek[1], r, d);
den = lcm(vector(#ek[1], t, denominator(ek[1][t])));
cs = vector(r+1, i, cs[i]*den);
gg = 0; for(i = 1, r+1, gg = gcd(gg, content(cs[i])));
cs = vector(r+1, i, cs[i]/gg);
read(concat(outdir, "c12_a2_long.seq")); seqsL = [QV, PV]; NL = #QV - 1;
print("fit used n <= 189 (b <= 192);  verifying exactly on n = 1 .. ", NL - r);
print("total failures over Q = ", verify(seqsL, cs, r, 1, NL - r));
print("failures on FRESH n = 200..", NL - r, " : ", verify(seqsL, cs, r, 200, NL - r));
}
quit;
