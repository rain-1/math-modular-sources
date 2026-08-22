\\ 04_hadamard.gp -- the Hadamard host: build the products, fit their minimal
\\ recurrences, read off rank + singular sets, check Hadamard's theorem.
default(parisizemax, 12000000000);
\p 6000
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read("/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/lib_fit.gp");
read(concat(outdir, "rows_raw.gp"));
read(concat(outdir, "rows_full.gp"));
NMAX = #AZ;
print("NMAX = ", NMAX, "   #BN = ", #BN);

\\ sanity: the extended companion still gives the Nesterenko linear form rate
{ for(k = 1, 8, my(n = 50*k);
    printf("  n=%3d  log|a^N G - b^N|/n = %+.6f   (target %+.6f)\n",
      n, log(abs(AN[n]*Catalan - BN[n]))/n, 2*log((3303-437*sqrt(57))/144))); }

HZZ = vector(NMAX, n, AZ[n]*AN[n]);          \\ A_Z (*) A_N
HZB = vector(NMAX, n, AZ[n]*BN[n]);          \\ A_Z (*) B_N
HNB = vector(NMAX, n, AN[n]*BZ[n]);          \\ A_N (*) B_Z
HBB = vector(NMAX, n, BZ[n]*BN[n]);          \\ B_Z (*) B_N
WW  = vector(NMAX, n, HZB[n] - HNB[n]);      \\ W

write(concat(outdir,"had.gp"), "HZZ = ", HZZ, ";");
write(concat(outdir,"had.gp"), "HZB = ", HZB, ";");
write(concat(outdir,"had.gp"), "HNB = ", HNB, ";");
write(concat(outdir,"had.gp"), "HBB = ", HBB, ";");
write(concat(outdir,"had.gp"), "WW  = ", WW,  ";");
print("wrote out/had.gp");

print("");
print("### archimedean growth  log|c_n|/n ###");
{ my(names = ["A_Z(*)A_N","A_Z(*)B_N","A_N(*)B_Z","B_Z(*)B_N","W"],
     seqs = [HZZ,HZB,HNB,HBB,WW]);
  print("   n:      100      200      300      400");
  for(s = 1, 5,
    printf("  %-10s %8.5f %8.5f %8.5f %8.5f\n", names[s],
      log(abs(seqs[s][100]))/100, log(abs(seqs[s][200]))/200,
      log(abs(seqs[s][300]))/300, log(abs(seqs[s][400]))/400)); }
{ my(l1z=15*log((1+sqrt(5))/2), l2z=-15*log((1+sqrt(5))/2),
     l1n=2*log((3303+437*sqrt(57))/144), l2n=2*log((3303-437*sqrt(57))/144));
  printf("  predicted: A_Z(*)A_N %.5f ; W  max(%.5f, %.5f) = %.5f\n",
    l1z+l1n, l1z+l2n, l1n+l2z, max(l1z+l2n, l1n+l2z)); }

\\ The minimal operator has (order,degree) = (4,96); a single sequence gives only
\\ ~390 equations, which is not enough for 5*97 unknowns.  The order search is done
\\ jointly on all four products in 05_fitprod.gp, and the exact operator in 06_operator.gp.
\q
