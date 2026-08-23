/* lib2.gp -- Hauptmoduls of Gamma_0(N)+W as (roots of) eta quotients, and the
   chi_{-4} rows they carry.
   Ligozat: ord at the cusp of denominator c (in the local parameter q_h) of
     prod_d eta(d tau)^{r_d}   is   (N/(24 gcd(c,N/c) c)) * sum_d gcd(c,d)^2 r_d/d.
   For a hauptmodul t with div(t) = (cusp N) - (pole), t is a modular unit, so
   the r_d solve a linear system; rational r_d are allowed (t is then a root of
   an eta quotient, still with a q-expansion in Q[[q]]).
   D log t = (1/24) sum_d r_d d E_2(d tau),  t = q exp( D^{-1}(D log t - 1) ).
*/

ligoz(nn, cc, dd) = (nn/(24*gcd(cc,nn/cc)*cc)) * gcd(cc,dd)^2/dd;

/* returns the exponent vector r (indexed by divisors of N) for the unit with
   prescribed orders ordv at the cusps (indexed by divisors), or 0 */
{unitexp(nn, ordv) =
  my(dv = divisors(nn), k = #dv);
  my(mm = matrix(k,k,i,j, ligoz(nn, dv[i], dv[j])));
  my(sol = matsolve(mm, ordv~));
  sol~;
}

E2ser(nn) = 1 - 24*sum(m=1,nn-1, sigma(m)*x^m) + O(x^nn);
{Vop(s, dd, nn) = my(r=O(x^nn)); for(m=0,(nn-1)\dd, r += polcoeff(s,m)*x^(m*dd)); r;}

{hauptmod(nn, rvec, nq) =
  my(dv = divisors(nn), lg = O(x^nq));
  for(j=1,#dv, if(rvec[j]!=0, lg += rvec[j]*dv[j]*Vop(E2ser(nq), dv[j], nq)));
  lg = lg/24;
  if(abs(polcoeff(lg,0)-1)>0, error("normalisation: constant term of Dlog t = ", polcoeff(lg,0)));
  my(ll = sum(m=1,nq-1, polcoeff(lg,m)/m*x^m) + O(x^nq));
  x*exp(ll);
}
