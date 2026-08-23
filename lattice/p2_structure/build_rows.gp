/* lattice/p2_structure/build_rows.gp
   Cache the exact integer rows of the Zudilin(3n) x Nesterenko(4,7)(n) pair.
   Prepend lattice/positivity/rows_pos.gp.  Set NLO, NHI, OUTF in the driver.
   Writes one line per n:   n X Y V U    (space separated exact integers)
   X_n = 2^{e_{3n}} D_{6n}^2 Q_{3n},  Y_n = 2^{e_{3n}} D_{6n}^2 P_{3n},
   V_n = 4^{7n+1} D_{6n}^2 B_n,       U_n = 4^{7n} D_{6n}^2 C_n.          */

{
buildrows(NLO, NHI, OUTF) =
 for(n=NLO, NHI,
   my(t0=getabstime(), zr=zudrow(n), nr=nestrow(n));
   if(denominator(zr[1])!=1 || denominator(zr[2])!=1 ||
      denominator(nr[1])!=1 || denominator(nr[2])!=1,
      error("non-integral row at n=", n));
   write(OUTF, n, " ", zr[1], " ", zr[2], " ", nr[1], " ", nr[2]);
   printf("built n=%d  ms=%d\n", n, getabstime()-t0));
}
