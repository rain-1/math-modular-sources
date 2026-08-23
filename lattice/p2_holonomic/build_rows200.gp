/* build_rows200.gp -- extend the exact row cache to n = 200 using the fitted
   order-3 recurrence for B_n, C_n and Zudilin's own order-2 recurrence.
   Prepend lattice/positivity/rows_pos.gp, lattice/p2_structure/p2core.gp,
   lattice/p2_holonomic/rowrec.gp.                                          */
default(parisize, 6000000000);
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data/";
OUT = "/home/ubuntu/code/math-modular-sources/lattice/p2_holonomic/data/rows_n200.txt";
NHI = 200;
RW  = rdrows(concat(DIR,"rows_all.txt"));
BV = vector(120); CV = vector(120);
{
for(n=4,120, my(rw=mapget(RW,n), DD=lcm(vector(6*n,i,i)), SS=DD^2);
  BV[n]=rw[3]/(4^(7*n+1)*SS); CV[n]=rw[4]/(4^(7*n)*SS));
}
AB = vector(117,i,BV[i+3]); AC = vector(117,i,CV[i+3]);
RB = fitrec2(AB,3,14); RC = fitrec2(AC,3,14);
print("recurrences equal for B and C: ", RB==RC);
{ my(bad=0); for(i=3,116, if(sum(j=0,3, subst(RB[j+1],'m,i)*AB[i+1-j])!=0, bad++));
  print("B relation failures on the cache: ", bad); }
{ my(bad=0); for(i=3,116, if(sum(j=0,3, subst(RC[j+1],'m,i)*AC[i+1-j])!=0, bad++));
  print("C relation failures on the cache: ", bad); }
AB = extend(AB, RB, NHI-3);  AC = extend(AC, RC, NHI-3);
ZZ = zud(3*NHI);
{
for(n=4,NHI,
  my(DD=lcm(vector(6*n,i,i)), SS=DD^2, m=3*n,
     XX=2^ee(m)*SS*ZZ[1][m+1], YY=2^ee(m)*SS*ZZ[2][m+1],
     VV=4^(7*n+1)*SS*AB[n-3], UU=4^(7*n)*SS*AC[n-3]);
  if(n<=120, my(rw=mapget(RW,n));
     if([XX,YY,VV,UU]!=rw, print("MISMATCH at n=",n)));
  if(type(XX)!="t_INT" || type(YY)!="t_INT" || type(VV)!="t_INT" || type(UU)!="t_INT",
     print("NONINTEGER at n=",n));
  write(OUT, n, " ", XX, " ", YY, " ", VV, " ", UU));
}
RECF = "/home/ubuntu/code/math-modular-sources/lattice/p2_holonomic/data/nest_recurrence.txt";
write(RECF, "/* order-3 P-recursion annihilating BOTH B_n and C_n, the Nesterenko (4,7) Catalan row. */");
write(RECF, "/* With A[i] the value at n = i+3, the relation is  sum_{j=0..3} R[j+1](m=i) * A[i+1-j] = 0. */");
write(RECF, "/* Fitted from n = 4..80 only; exact on n = 81..120; used here to build n <= 200. */");
write(RECF, RB);
print("done");
\q
