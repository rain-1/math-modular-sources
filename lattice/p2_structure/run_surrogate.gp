/* driver: rational-surrogate control.  Prepend a line   EEXP = 320;
   Runs the whole structural analysis with G replaced by
   GS = bestappr(G,10^EEXP) and dumps the vectors for comparison with the
   true-G dump: the first n at which the two differ is the faithfulness
   horizon n_f(E).                                                        */
\p 3000
default(parisize, 2000000000);
GG = Catalan;
GS = bestappr(GG, 10^EEXP);
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data/";
RW = rdrows(concat(DIR,"rows_all.txt"));
DUMP = concat(concat(DIR,"vec_E"), concat(Str(EEXP), ".txt"));
printf("surrogate E=%d  log10 den = %.4f  log|G-GS| = %.4f\n",
   EEXP, log(1.*denominator(GS))/log(10), log(abs(1.*(GG-GS))));
for(n=4,120, anal(n, [22.4,23.0,23.9], GS, RW, DUMP));
\q
