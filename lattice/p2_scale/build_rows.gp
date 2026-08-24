/* lattice/p2_scale/build_rows.gp
   Extend the exact row cache X_n, Y_n, V_n, U_n to n = NHI (default 1000),
   reusing the order-3 degree-14 P-recursion of lattice/p2_holonomic/rowrec.gp
   for the Nesterenko (4,7) entries B_n, C_n and Zudilin's own order-2
   recurrence for Q_m, P_m.  Everything is exact rational/integer arithmetic.
   Checks performed on the way:
     * the fitted operator annihilates B and C on the whole n <= 120 cache;
     * the rebuilt rows agree digit-for-digit with the cached n <= 200 rows;
     * X_n, Y_n, V_n, U_n are integers at every n (a strong non-trivial test
       of the extension: the recurrence works with rationals throughout).
   Prepend: lattice/positivity/rows_pos.gp, lattice/p2_structure/p2core.gp,
            lattice/p2_holonomic/rowrec.gp.
   No builtin names shadowed (no psi, M, Phi, S, cmp).                      */

default(parisize, 6000000000);
default(parisizemax, 13000000000);

/* D_m = lcm(1..m) by prime powers -- O(m/log m) big multiplications */
dlcm(mm) = prod(i=1, #(my(pv=primes([2,mm]))), pv[i]^logint(mm,pv[i]));
{
dlcm(mm) = my(pv=primes([2,mm]));
  prod(i=1,#pv, pv[i]^logint(mm,pv[i]));
}

NHI  = if(type(NHI0)=="t_INT", NHI0, 1000);
DIRA = "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data/";
DIRB = "/home/ubuntu/code/math-modular-sources/lattice/p2_holonomic/data/";
OUT  = "/home/ubuntu/code/math-modular-sources/lattice/p2_scale/data/rows_scale.txt";

RW  = rdrows(concat(DIRA,"rows_all.txt"));       /* n <= 120, from the O(n^4.1) solve */
RW2 = rdrows(concat(DIRB,"rows_n200.txt"));      /* n <= 200 */

BV = vector(120); CV = vector(120);
{
for(n=4,120, my(rw=mapget(RW,n), DD=dlcm(6*n), SS=DD^2);
  BV[n]=rw[3]/(4^(7*n+1)*SS); CV[n]=rw[4]/(4^(7*n)*SS));
}
AB = vector(117,i,BV[i+3]); AC = vector(117,i,CV[i+3]);
RB = fitrec2(AB,3,14); RC = fitrec2(AC,3,14);
print("[rows] one operator for both B and C: ", RB==RC);
{ my(bad=0); for(i=3,116, if(sum(j=0,3, subst(RB[j+1],'m,i)*AB[i+1-j])!=0, bad++);
                          if(sum(j=0,3, subst(RB[j+1],'m,i)*AC[i+1-j])!=0, bad++));
  print("[rows] operator failures on the n<=120 cache (B and C): ", bad); }

AB = extend(AB, RB, NHI-3);  AC = extend(AC, RC, NHI-3);
ZZ = zud(3*NHI);
BADI = 0; BADM = 0;
{
for(n=4,NHI,
  my(DD=dlcm(6*n), SS=DD^2, mm=3*n,
     XX=2^ee(mm)*SS*ZZ[1][mm+1], YY=2^ee(mm)*SS*ZZ[2][mm+1],
     VV=4^(7*n+1)*SS*AB[n-3], UU=4^(7*n)*SS*AC[n-3]);
  if(type(XX)!="t_INT" || type(YY)!="t_INT" || type(VV)!="t_INT" || type(UU)!="t_INT",
     BADI++; print("[rows] NONINTEGER at n=",n));
  if(n<=200, my(rw=mapget(RW2,n)); if([XX,YY,VV,UU]!=rw, BADM++; print("[rows] MISMATCH at n=",n)));
  write(OUT, n, " ", XX, " ", YY, " ", VV, " ", UU);
  if(n%100==0, print("[rows] n=",n," ok, log10 X = ",sizedigit(XX))));
}
print("[rows] noninteger count = ", BADI, "   mismatch vs cached n<=200 = ", BADM);
print("[rows] done, NHI = ", NHI);
\q
