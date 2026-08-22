/* Kubota-Leopoldt L_p(s, chi*omega^m) via Washington Thm 5.11.
   chi encoded as [f,[chi(1..f)]] (rational values).                */
setbern(J) = if(BCACHE[1] < J, BCACHE = [J, vector(J+1,i,bernfrac(i-1))]);
BCACHE = [0,0];
bern(j) = BCACHE[2][j+1];
chival(chi,a) = chi[2][ ((a-1) % chi[1]) + 1 ];
tch(a,p,PR) = if(p==2, if(a%4==1,1,-1), teichmuller(a+O(p^PR)));
LpG(p, chi, m, s, PR) =
{ my(q=if(p==2,4,p), f=chi[1], F=lcm(q,f), J, tot, br, inn, x, tt, c);
  J = ceil((PR+5)/valuation(F,p)) + 3; setbern(J);
  tot = O(p^PR);
  for(a=1,F, if(a%p==0, next);
    c = chival(chi,a); if(c==0, next);
    c *= tch(a,p,PR+8)^m;
    br = (a+O(p^(PR+8)))/tch(a,p,PR+8);
    x = (F/a)+O(p^(PR+8)); inn = O(p^(PR+8)); tt = 1+O(p^(PR+8));
    for(j=0,J, inn += binomial(1-s,j)*tt*bern(j); tt *= x);
    tot += c*br^(1-s)*inn);
  tot/(F*(s-1));
}
triv  = [1,[1]];
chim3 = [3,[1,-1,0]];
chim4 = [4,[1,0,-1,0]];
