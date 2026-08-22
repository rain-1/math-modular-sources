/* ============================================================
   lattice/acf_one_surface/lib.gp
   Shared utilities: eta products, Zagier rows, Eisenstein series.
   Variable for q-expansions: 'q .   (avoid PARI builtins as names)
   ============================================================ */

NQ = 260;          \\ working q-precision

/* Euler product E(k) = prod_{n>=1} (1 - q^{kn}) , via PARI's eta on a series */
{ etaprod(k, prec=NQ) = eta('q^k + O('q^prec)); }

/* Zagier row: (n+1)^2 u_{n+1} = (aa*(n^2+n)+bb) u_n - dd*n^2 u_{n-1}, u_0=1 */
{ zrow(aa,bb,dd,N) =
  my(u = vector(N+2)); u[1] = 1;
  for(n=0, N,
      my(Pn = aa*(n^2+n)+bb, Qn = dd*n^2, prv = if(n==0, 0, u[n]));
      u[n+2] = (Pn*u[n+1] - Qn*prv)/(n+1)^2);
  vector(N+1, k, u[k]);
}

/* character chi_{-3} */
{ chim3(n) = my(r = n % 3); if(r==0, 0, if(r==1, 1, -1)); }

/* Ess = E_{3,chi_{-3},1} : coefficient  sum_{d|m} chi(m/d) d^2 , no constant term */
{ Ess(prec=NQ) =
  my(v = vector(prec)); \\ v[m] = coeff of q^m
  for(m=1, prec-1, my(s=0); fordiv(m, d, s += chim3(m/d)*d^2); v[m]=s);
  Ser(concat([0], vector(prec-1, m, v[m])), 'q, prec);
}

/* Ett = E_{3,1,chi_{-3}} : coefficient sum_{d|m} chi(d) d^2 , constant term set to 0
   (the true constant term is L(-2,chi_{-3})/2 ; every source we build has zero
    constant term, and all our combinations P(V_2) with P(1)=0 kill it.)      */
{ Ett(prec=NQ) =
  my(v = vector(prec));
  for(m=1, prec-1, my(s=0); fordiv(m, d, s += chim3(d)*d^2); v[m]=s);
  Ser(concat([0], vector(prec-1, m, v[m])), 'q, prec);
}

/* V_d f (tau) = f(d tau) */
{ Vop(d, f) = my(prec = serprec(f,'q)); subst(f, 'q, 'q^d) + O('q^prec); }

/* D = q d/dq */
{ Dop(f) = 'q * deriv(f, 'q); }

/* Dinv2 : formal D^{-2} on series with zero constant term */
{ Dinv2(f) =
  my(prec = serprec(f,'q), v = Vec(f + O('q^prec)), val = valuation(f,'q));
  my(w = vector(prec-1));
  for(m=1, prec-1, w[m] = polcoef(f, m, 'q)/m^2);
  Ser(concat([0], w), 'q, prec);
}

/* twist q -> -q  (tau -> tau + 1/2) */
{ twist(f) = my(prec = serprec(f,'q)); subst(f, 'q, -'q) + O('q^prec); }

/* substitute a power series s (val>=1) into a coefficient vector: sum c[n] s^n */
{ evalrow(cv, s, prec) =
  my(r = 0*'q + O('q^prec), pw = 1 + O('q^prec));
  for(n=1, #cv, r += cv[n]*pw; pw = pw*s + O('q^prec));
  r;
}

/* t-expansion of a q-series G w.r.t. a coordinate t = q + O(q^2):
   returns the series G(t) = sum c_n t^n (variable is still 'q, now meaning t) */
{ texp(G, t) =
  my(prec = min(serprec(G,'q), serprec(t,'q)));
  subst(G + O('q^prec), 'q, serreverse(t + O('q^prec)));
}
