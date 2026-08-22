\\ 10_folds.gp -- which singular points each function is regular at.
\\ Measured by the exponential growth rate of its coefficients (= 1/|nearest singularity|).
default(parisizemax, 12000000000);
\p 6000
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read(concat(outdir, "rows_raw.gp"));
read(concat(outdir, "rows_full.gp"));
read(concat(outdir, "had.gp"));
N = #HZZ; G = Catalan;
phi = (1+sqrt(5))/2; TT = (3303+437*sqrt(57))/144; TM = (3303-437*sqrt(57))/144;
s1 = phi^(-15); s2 = -phi^15; t1 = TT^(-2); t2 = TM^(-2);
{ printf("singular points of A_Z : s1 = %.10e  s2 = %.10e\n", s1, s2);
  printf("singular points of A_N : t1 = %.10e  t2 = %.10e\n", t1, t2);
  printf("products: s1t1=%.10e  s1t2=%.10e  s2t1=%.10e  s2t2=%.10e\n",
    s1*t1, s1*t2, s2*t1, s2*t2); }

\\ true conditional functions (G exact, not a test rational)
CZ = vector(N, n, AZ[n]*G - BZ[n]);        \\ Zudilin fold-regular form
CN = vector(N, n, AN[n]*G - BN[n]);        \\ Nesterenko fold-regular form
CONDX = vector(N, n, AZ[n]*CN[n]);         \\ A_Z (*) (B_N - G A_N), up to sign
COND2 = vector(N, n, AN[n]*CZ[n]);         \\ A_N (*) (B_Z - G A_Z)
DBL   = vector(N, n, CZ[n]*CN[n]);         \\ doubly conditional
print("");
print("  growth log|c_n|/n  (n = 100,200,300,400)   vs   -log|nearest singularity|");
{ my(nms = ["A_Z(*)A_N","W","A_Z(*)(B_N-G A_N)","A_N(*)(B_Z-G A_Z)","(B_Z-G A_Z)(*)(B_N-G A_N)"],
     sq = [HZZ, WW, CONDX, COND2, DBL],
     pr = [-log(abs(s1*t1)), -log(abs(s2*t1)), -log(abs(s1*t2)), -log(abs(s2*t1)), 0]);
  for(s = 1, 5,
    printf("  %-28s", nms[s]);
    for(k = 1, 4, my(n = 100*k); printf(" %9.5f", log(abs(sq[s][n]))/n));
    printf("   pred %9.5f\n", pr[s])); }
{ printf("\n  reference rates: -log|s1t1|=%.5f  -log|s1t2|=%.5f  -log|s2t1|=%.5f  -log|s2t2|=%.5f\n",
    -log(abs(s1*t1)), -log(abs(s1*t2)), -log(abs(s2*t1)), -log(abs(s2*t2))); }
\q
