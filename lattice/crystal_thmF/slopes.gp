/* slopes.gp -- p-adic Frobenius-slope census for Apery rows (crystal criterion).
   mu1 = rate of v_p(a_n)/n ; muC = rate of v_p(Casoratian)/n ;
   sigma_p = muC - 2*mu1  (tail rate of b_n/a_n).  Limit exists iff sigma_p > 0. */
default(parisizemax, 6000000000);
read("lattice/euler_criterion/rows.gp");
N = 300;
{ ROWS = [
 ["A     (7,2,-8)  R2", rowR2(7,2,-8,N)],
 ["B     (9,3,27)  R2", rowR2(9,3,27,N)],
 ["C     (10,3,9)  R2", rowR2(10,3,9,N)],
 ["D z2  (11,3,-1) R2", rowR2(11,3,-1,N)],
 ["E     (12,4,32) R2", rowR2(12,4,32,N)],
 ["F     (17,6,72) R2", rowR2(17,6,72,N)],
 ["gamma z3 (17,5,1)R3", rowR3(17,5,1,N)],
 ["alpha Domb      R3", rowR3(10,4,64,N)],
 ["delta (7,3,81)  R3", rowR3(7,3,81,N)],
 ["eps   (12,4,16) R3", rowR3(12,4,16,N)],
 ["zeta  (9,3,-27) R3", rowR3(9,3,-27,N)],
 ["eta   (11,5,125)R3", rowR3(11,5,125,N)],
 ["s18 Cooper      R3", rowS18(N)],
 ["cusp L(f,2) N=12  ", rowCusp(N)],
 ["Zudilin Catalan   ", rowZud(N)]
]; }
PR = [2,3,5,7,11,13];
vp(x,p) = if(x==0, 0, valuation(x,p));
print("row                 |  p |   mu1  |   muC  | sigma_p | verdict");
{ for(r=1,#ROWS,
  my(nm=ROWS[r][1], A=ROWS[r][2][1], B=ROWS[r][2][2], Cas=vector(N));
  for(n=1,N-1, Cas[n]=A[n+2]*B[n+1]-A[n+1]*B[n+2]);
  for(j=1,#PR,
    my(p=PR[j], n1=200, n2=N-2);
    my(m1=(vp(A[n2+1],p)-vp(A[n1+1],p))/(n2-n1));
    my(mC=(vp(Cas[n2],p)-vp(Cas[n1],p))/(n2-n1));
    my(sg=mC-2*m1, vd);
    vd = if(sg>0, "LIMIT", if(sg==0, "flat(no limit)", "DIVERGES"));
    print(nm, " | ", p, " | ", m1, " | ", mC, " | ", sg, " | ", vd));
  print("--------------------+----+--------+--------+---------+-------")); }
