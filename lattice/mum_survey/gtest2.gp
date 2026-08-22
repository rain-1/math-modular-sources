read("gammap.gp");
read("../euler_criterion/lp.gp");
/* L_p(s, omega^m) via Washington Thm 5.11, character values Teichmuller powers */
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
{
for(pi=1,6,
  my(p=[2,3,5,7,11,13][pi], PR=70, K=16);
  my(c=gpTaylor(p,K,PR), L3=6*c[3], L5=120*c[5]);
  my(z3t=Lp(p,triv,3,45), z3w=Lpom(p,-2,3,45));
  print("p=",p);
  print("   v(c2)=",valuation(c[2],p)," v(c4)=",valuation(c[4],p),"  (expect ~precision-loss => 0)");
  print("   L3/L_p(3,triv)   = ", L3/z3t);
  print("   L3/L_p(3,om^-2)  = ", L3/z3w);
);
}
quit
