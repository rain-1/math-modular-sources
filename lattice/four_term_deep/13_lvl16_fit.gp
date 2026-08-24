/* Fit the five-term recurrence to a_n = [x^n] A(x) for the level-16 Catalan host. */
default(parisize, 2000000000);
AF = "/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl16_A.txt";
L  = readvec(AF);
NA = #L - 1;           /* a_0 .. a_NA  ; index a(n) = L[n+1] */
a(n) = L[n+1];
print("have a_0..a_", NA);

/* --- ansatz 1: (n+1)^2 a_{n+1} = P a_n - Q a_{n-1} + R a_{n-2} - T5 a_{n-3}
       unknowns u = [a,b,c, d,e,f, g,h,j, k,l,m]
       row: [n^2 a_n, n a_n, a_n, -n^2 a_{n-1}, -n a_{n-1}, -a_{n-1},
             n^2 a_{n-2}, n a_{n-2}, a_{n-2}, -n^2 a_{n-3}, -n a_{n-3}, -a_{n-3}]
       rhs: (n+1)^2 a_{n+1}                                                     */
rowvec(n) = [n^2*a(n), n*a(n), a(n), -n^2*a(n-1), -n*a(n-1), -a(n-1), n^2*a(n-2), n*a(n-2), a(n-2), -n^2*a(n-3), -n*a(n-3), -a(n-3)]~;
NFIT = 30;
MM = matconcat(vector(NFIT, i, rowvec(i+2)));   /* n = 3..NFIT+2 */
RHS = vectorv(NFIT, i, my(n=i+2); (n+1)^2*a(n+1));
print("system size: ", matsize(MM), "  rank = ", matrank(MM));
print("nullspace dim of the 12-col matrix = ", #matker(MM));
SOL = 0;
{ if(matrank(MM)==12 && matrank(matconcat([MM,RHS]))==12,
     SOL = matsolve(MM[,1..12]~ , 0);  /* placeholder */
  );}
/* solve least-structured: use first 12 independent rows */
{
  my(Mt = MM~, Rt = RHS);   /* Mt is NFIT x 12 */
  my(sel = List(), Cur = matrix(0,12));
  for(i=1,NFIT,
     my(Try = matconcat([Cur; Mt[i,]]));
     if(matrank(Try)>matrank(Cur), Cur = Try; listput(sel,i));
     if(#sel==12, break));
  print("selected rows (n values): ", vector(#sel,i,sel[i]+2));
  if(#sel<12, print("*** LOUD: cannot reach rank 12 -- ansatz underdetermined"),
    my(Ms = matrix(12,12,i,j, Mt[sel[i],j]), Rs = vectorv(12,i, Rt[sel[i]]));
    SOL = matsolve(Ms, Rs);
    print("SOL = ", SOL~);
  );
}
{ if(type(SOL)=="t_COL",
  my(u=SOL, aa=u[1],bb=u[2],cc=u[3],dd=u[4],ee=u[5],ff=u[6],gg=u[7],hh=u[8],jj=u[9],kk=u[10],ll=u[11],mm=u[12]);
  print("\n(a,b,c) = ",[aa,bb,cc]);
  print("(d,e,f) = ",[dd,ee,ff]);
  print("(g,h,j) = ",[gg,hh,jj]);
  print("(k,l,m) = ",[kk,ll,mm]);
  /* full verification */
  my(bad=0);
  for(n=3,NA-1,
    my(lhs=(n+1)^2*a(n+1),
       rhs=(aa*n^2+bb*n+cc)*a(n) - (dd*n^2+ee*n+ff)*a(n-1) + (gg*n^2+hh*n+jj)*a(n-2) - (kk*n^2+ll*n+mm)*a(n-3));
    if(lhs!=rhs, bad=n; break));
  print("\nrecurrence holds for ALL n=3..",NA-1,"? ", bad==0, if(bad, Str("  FIRST FAILURE n=",bad), ""));
  print("(a,d,g,k) = ", [aa,dd,gg,kk], "   predicted (10,40,80,64) or (-10,40,-80,64)");
  print("char poly lam^4 - a lam^3 + d lam^2 - g lam + k factors: ", factor(lam^4-aa*lam^3+dd*lam^2-gg*lam+kk));
  print("roots (numeric): ", polroots(lam^4-aa*lam^3+dd*lam^2-gg*lam+kk));
  /* class invariants */
  Rc = 1 - aa*t + dd*t^2 - gg*t^3 + kk*t^4;
  Sc = 1 - (aa+bb)*t + (3*dd+ee)*t^2 - (5*gg+hh)*t^3 + (7*kk+ll)*t^4;
  Vc = -cc + (dd+ee+ff)*t - (4*gg+2*hh+jj)*t^2 + (9*kk+3*ll+mm)*t^3;
  Tc = 1 - bb*t + (dd+ee)*t^2 - (2*gg+hh)*t^3 + (3*kk+ll)*t^4;
  print("\nRc = ",Rc);
  print("Sc = ",Sc);
  print("Vc = ",Vc);
  print("Tc = Sc - t Rc' = ",Tc, "   check: ", Sc - t*deriv(Rc,t) - Tc);
  print("Rc factors: ", factor(Rc));
  default(realprecision,40);
  my(rts = polroots(Rc), Rp = deriv(Rc,t));
  for(i=1,#rts, my(ti=rts[i]);
     print("  t_",i," = ",ti,"   rho = ", -subst(Tc,t,ti)/(ti*subst(Rp,t,ti))));
  /* T5 roots */
  my(T5 = kk*nn^2+ll*nn+mm);
  print("\nT5(n) = ", T5, "   factors: ", factor(T5));
  my(sr = polroots(T5));
  print("  s roots = ", sr, "   delta_infinity = |s2-s1| = ", abs(sr[2]-sr[1]));
  print("  exponents at infinity 3-s: ", [3-sr[1], 3-sr[2]]);
  print("  T5(n0)=0 for integer n0>=1? ", vector(12,n, [n, kk*n^2+ll*n+mm]));
)}
