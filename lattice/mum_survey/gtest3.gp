read("gammap.gp");
read("../euler_criterion/lp.gp");
Lpom(p, m, s, PR) =
{ my(q=if(p==2,4,p), F=q, J, tot, br, inn, x, tt, c);
  J = ceil((PR+5)/valuation(F,p)) + 3; setbern(J);
  tot = O(p^PR);
  for(a=1,F, if(a%p==0, next);
    c = if(p==2, if(a%4==1,1,-1)^m, teichmuller(a+O(p^(PR+8)))^m);
    br = (a+O(p^(PR+8)))/tch(a,p,PR+8);
    x = (F/a)+O(p^(PR+8)); inn=O(p^(PR+8)); tt=1+O(p^(PR+8));
    for(j=0,J, inn += binomial(1-s,j)*tt*bern(j); tt*=x);
    tot += c*br^(1-s)*inn);
  tot/(F*(s-1));
}
/* how many digits of X agree with the rational r ? */
agree(X,r,p) = my(d=X-r); if(d==0, "exact", valuation(d,p));
{
for(pi=1,7,
  my(p=[2,3,5,7,11,13,17][pi], PR=200, K=40);
  my(c=gpTaylor(p,K,PR));
  my(L2=2*c[2], L3=6*c[3], L4=24*c[4], L5=120*c[5], L7=5040*c[7]);
  my(z3=Lpom(p,-2,3,80), z5=Lpom(p,-4,5,80), z7=Lpom(p,-6,7,80));
  print("p=",p, "  [K=",K," PR=",PR,"]");
  print("   v_p(L2)=",valuation(L2,p),"   v_p(L4)=",valuation(L4,p), "   (0 up to Vandermonde loss)");
  print("   L3 / L_p(3,om^-2) + 2  -> v_p = ", valuation(L3/z3+2,p));
  print("   L5 / L_p(5,om^-4) + 24 -> v_p = ", valuation(L5/z5+24,p));
  print("   L7 / L_p(7,om^-6) + 720-> v_p = ", valuation(L7/z7+720,p));
);
}
quit
