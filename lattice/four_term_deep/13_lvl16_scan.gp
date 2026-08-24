/* General search for the minimal recurrence satisfied by a_n.  (fixed transposes) */
default(parisize, 3000000000);
L  = readvec("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl16_A.txt");
NA = #L - 1;
a(n) = L[n+1];
print("a_0..a_",NA);
/* mkmat(r,D,NF): NF x (r+1)(D+1) matrix, unknowns = coeffs of P_s(n), s=0..r,
   relation sum_{s=0..r} P_s(n) a_{n+1-s} = 0, equations n = r..r+NF-1        */
{mkmat(r,D,NF) = matrix(NF, (r+1)*(D+1), i, j,
   my(n = i+r-1, s = (j-1)\(D+1), e = (j-1)%(D+1));
   n^(D-e) * a(n+1-s));}
{ for(r=2,7,
    for(D=1,5,
      my(NU=(r+1)*(D+1), NF=NU+45);
      if(NF+r+2 > NA, next);
      my(kd = #matker(mkmat(r,D,NF)));
      print("  terms=",r+1," deg=",D," unknowns=",NU," eqs=",NF," ker dim = ",kd)));
}
