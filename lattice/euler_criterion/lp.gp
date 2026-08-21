/* ================================================================
   lp.gp  --  Kubota-Leopoldt p-adic L-functions L_p(s,chi) at
   integer s, via Washington, "Introduction to Cyclotomic Fields",
   Thm 5.11:

     L_p(s,chi) = 1/(F(s-1)) * sum_{a=1,(a,p)=1}^{F} chi(a)<a>^{1-s}
                    * sum_{j>=0} binom(1-s,j) (F/a)^j B_j ,

   F any common multiple of q (= p, or 4 if p=2) and cond(chi),
   <a> = a/omega(a).  For integer s the binomial coefficients are
   exact integers, so no p-adic precision is lost in them.

   A character is encoded as [f, [chi(1),...,chi(f)]].
   ================================================================ */

triv   = [1,[1]];
chim3  = [3,[1,-1,0]];            /* chi_{-3}, odd, cond 3   */
chim4  = [4,[1,0,-1,0]];          /* chi_{-4}, odd, cond 4   */
chi12  = [12,[1,0,0,0,-1,0,-1,0,0,0,1,0]];  /* chi_{-3}chi_{-4}, even */
chi5   = [5,[1,-1,-1,1,0]];       /* Legendre mod 5, even    */

chival(chi,a) = chi[2][ ((a-1) % chi[1]) + 1 ];

/* Teichmuller in Washington's normalisation (p=2: the character mod 4) */
tch(a,p,PR) = if(p==2, if(a%4==1, 1, -1), teichmuller(a + O(p^PR)));

BCACHE = [0,0];   /* [Jmax, vector of B_j, j=0..Jmax] */
bern(j) = BCACHE[2][j+1];
setbern(J) = if(BCACHE[1] < J, BCACHE = [J, vector(J+1,i,bernfrac(i-1))]);

Lp(p, chi, s, PR) =
{ my(q, f, F, J, tot, br, inn, x, tt, c);
  q = if(p==2, 4, p);  f = chi[1];  F = lcm(q,f);
  J = ceil((PR+5)/valuation(F,p)) + 3;
  setbern(J);
  tot = O(p^PR);
  for(a=1, F,
    if(a % p == 0, next);
    c = chival(chi,a); if(c == 0, next);
    br  = (a + O(p^(PR+8))) / tch(a,p,PR+8);
    x   = (F/a) + O(p^(PR+8));
    inn = O(p^(PR+8)); tt = 1 + O(p^(PR+8));
    for(j=0, J,
      inn += binomial(1-s,j) * tt * bern(j);
      tt  *= x);
    tot += c * br^(1-s) * inn);
  tot/(F*(s-1));
}

/* generalised Bernoulli number B_{n,chi} (exact rational) */
Bchi(chi,n) = my(f=chi[1]); f^(n-1)*sum(a=1,f, chival(chi,a)*subst(bernpol(n),'x,a/f));

/* product character */
chimul(c1,c2) = my(f=lcm(c1[1],c2[1])); [f, vector(f,a, chival(c1,a)*chival(c2,a))];
