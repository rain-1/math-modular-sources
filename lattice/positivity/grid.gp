/* lattice/positivity/grid.gp  --- Task 3 (a): the single-kernel grid.
   For 0 <= j <= m <= MMAX evaluate the Catalan-world Beukers moment
      mom(m,j) = I(a,b,a,b,a+1),  a = m+j, b = m,
   which is  K_Z^{(m-3j... )}: writing s = m-3j, t = j one has
      mom(m,j) = int int K_Z^{(s)} * v^t,   v = x^4(1-x)^3 y^4 (1-y)^3 (1-xy)^{-4},
   i.e. the (s,t) grid of PRODUCTS of a Zudilin kernel of index s with t copies
   of the Nesterenko shape.  Zudilin(m) = (m,0); Nesterenko(n) = (3n,n).
   Every mom(m,j) > 0 (positive integrand), so den*mom is an admissible
   POSITIVE-CONE linear form  q*G - p  with q,p in Z.
   Output CSV: m,j,s,logden_m,logval_m,obj_m,obj_n,logq_m,delta            */
\p 500
{
gridrow(m) =
 localprec(500);
 for(j=0,m,
  my(r=mom(m,j), AA=r[1], BB=r[2],
     dd=lcm(denominator(AA),denominator(BB)),
     val=AA*Catalan-BB, q=dd*AA, p=dd*BB, obj=dd*val);
  printf("%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
    m, j, m-3*j, log(dd)/m, log(val)/m, log(obj)/m, 3*log(obj)/m,
    log(abs(q))/m, -log(obj)/log(abs(q)));
 );
}
