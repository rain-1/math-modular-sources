/* ================================================================
   criterion.gp -- the Euler-factor divisibility criterion, checked
   symbolically on the twelve-family source table of
   paper/sections/02_sources.tex (Theorem A).

   Source polynomial P as a Dirichlet polynomial  P(s) = sum_d c_d d^{-s};
   we record it as a polynomial in the formal variables X_d = d^{-s}
   written multiplicatively, i.e. as a function of s.
   Euler factor at p:   E_p(s) = 1 - psi(p) p^{-s},  psi = character
   carried by the *depleted* variable (= character on the co-divisor of
   the weight-(w+2) source; "inner" placement => psi = chi, "outer" => psi = 1).
   Criterion: E_p | P  in Q[p^{-s}]  <=>  P vanishes at p^{-s} = 1/psi(p).
   ================================================================ */

/* P given as a vector of [d, c_d] pairs; evaluate at s */
Pev(P, s) = sum(i=1,#P, P[i][2]*P[i][1]^(-s));
/* P as a polynomial in y = p^{-s}, other divisors evaluated at s=w+1 */
Ppoly(P, p, w) =
{ my(y='y, t=0);
  for(i=1,#P,
    my(d=P[i][1], vp=valuation(d,p), dp=d/p^vp);
    t += P[i][2] * dp^(-(w+1)) * y^vp);
  t;
}

/* families: [name, w, P (list of [d,c_d]), psi(2), psi(3), psi(5)] */
{FAM = [
 ["A", 1, [[1,1],[2,-1]],                                            1, 1, 1],
 ["B", 1, [[1,1],[2,-6],[4,-8]],                                    -1, 0, 1],
 ["C", 1, [[1,1],[2,-8]],                                           -1, 0, 1],
 ["D", 1, [[1,1]],                                                   1, 1, 1],
 ["E", 1, [[1,1],[2,-8]],                                            0,-1, 1],
 ["F", 1, [[1,1],[2,-7],[4,-8]],                                    -1, 0, 1],
 ["alpha", 2, [[1,1],[2,-17],[3,-9],[4,16],[6,153],[12,-144]],       1, 1, 1],
 ["gamma", 2, [[1,1],[2,-28],[3,63],[6,-36]],                        1, 1, 1],
 ["delta", 2, [[1,1],[2,-14],[3,-1],[4,16],[6,14],[12,-16]],         1, 1, 1],
 ["eps",   2, [[1,1],[2,-21],[4,84],[8,-64]],                        1, 1, 1],
 ["zeta",  2, [[1,1]],                                              -1, 0,-1],
 ["eta",   2, [[1,1],[2,-14],[4,-16]],                              -1,-1, 0]
];}
/* psi = character on the depleted variable e:
   inner placement (sigma_chi):  psi = chi ;  outer (tilde sigma): psi = 1.
   A, D outer -> psi = 1 ; alpha,gamma,delta,eps trivial characters -> psi = 1;
   B,C,F inner chi_{-3}; E inner chi_{-4}; zeta inner chi_{-3}; eta inner chi_5. */

{
print("=== Euler-factor criterion vs. the slope set {p : p | c} ===");
CVAL = ["A",-8;"B",27;"C",9;"D",-1;"E",32;"F",72;"alpha",64;"gamma",1;"delta",81;"eps",16;"zeta",-27;"eta",125];
for(i=1,#FAM,
  my(f=FAM[i], nm=f[1], w=f[2], P=f[3], cc=0);
  for(j=1,matsize(CVAL)[1], if(CVAL[j,1]==nm, cc=CVAL[j,2]));
  print1("  ",nm,"  w=",w,"  P(",w+1,")=",Pev(P,w+1),"  c=",cc,"   ");
  for(k=1,3,
    my(p=[2,3,5][k], psi=f[3+k], pol, ok, Q);
    pol = Ppoly(P,p,w);
    if(psi==0,
      ok = 1; Q = pol,   /* Euler factor is 1 */
      ok = (subst(pol,'y, 1/psi)==0);
      if(ok, Q = pol/(1 - psi*'y), Q = 0));
    print1("p=",p,": crit=",if(ok,"YES","no"),
           if(ok, Str(" Q=", subst(Q,'y,p^(-(w+1)))), ""),
           " | p|c: ", if(cc%p==0,"YES","no"), "   "));
  print());
}
quit
