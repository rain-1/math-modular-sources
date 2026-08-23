/* lattice/p2_holonomic/adelic.gp
   The 2-adic side: the common 2-adic limit xi2 of the two rows, the 2-adic
   linear forms, and the adelic quality of the output pairs.
   Prepend lattice/positivity/rows_pos.gp, lattice/p2_structure/p2core.gp,
   lattice/p2_holonomic/hcore.gp.                                           */

default(parisizemax, 8000000000);

/* xi2 := lim_m P_m/Q_m in Q_2  ( = zeta_2(2), ZUDILIN_2ADIC.md, proved ).
   Returns [xi as t_PADIC, prec in bits, the exact rational P_m/Q_m].       */
{
xi2at(MTOP) =
 my(ZZ=zud(MTOP), r=ZZ[2][MTOP+1]/ZZ[1][MTOP+1],
    pr=8*MTOP-1-4*hammingweight(MTOP)-8);
 [r + O(2^pr), pr, r, ZZ];
}

/* v_2 of the 2-adic linear form of a row entry pair (A,B) against xi:
   v_2(A*xi - B), computed exactly from the rational representative.        */
vform(AA,BB,rat) = valuation(AA*rat - BB, 2);
