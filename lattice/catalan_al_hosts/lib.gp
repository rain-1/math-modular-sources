/* lib.gp -- shared machinery for the chi_{-4} Eisenstein rows on genus-0
   Atkin-Lehner quotient hosts.   Conventions of CATALAN_THREE_PERIOD.md.
   Avoids the builtin names psi, M, Phi, S, cmp.  */

NQ = 260;                              /* number of q-coefficients kept      */

chim(n) = if(n%2==0, 0, if(n%4==1, 1, -1));

/* exact real matrix of W_Q on a modular-form space, together with the global
   phase (1 or I) that W_Q carries; returns [Rmat, phase] with W_Q = phase*Rmat */
{alreal(mfx, qq) =
  my(al = mfatkininit(mfx, qq), mm = al[2]/al[3]);
  my(re = 0.);
  for(i=1,matsize(mm)[1], for(j=1,matsize(mm)[2], re = max(re, abs(real(mm[i,j])))));
  my(ph, rr);
  if(re > 1e-20, ph = 1; rr = real(mm), ph = I; rr = real(mm/I));
  [matrix(matsize(rr)[1],matsize(rr)[2],i,j, bestappr(rr[i,j],10^20)), ph];
}

/* joint eigenspace of the W_Q (Q in qs) with sign vector sgn (entries +-1),
   returned as a matrix whose columns are a basis in terms of mfbasis coords */
{jointeig(mfx, qs, sgn) =
  my(d = mfdim(mfx));
  if(d==0, return(matrix(0,0)));
  my(pr = matid(d));
  for(j=1,#qs,
    my(ar = alreal(mfx, qs[j])[1]);
    pr = pr * (matid(d) + sgn[j]*ar)/2;
  );
  matimage(pr);
}

/* q-expansion vector (length nn+1, index 1 = constant term) of a linear
   combination of mfbasis(mfx) with coefficient vector cv */
{qexp(mfx, cv, nn) =
  my(bb = mfbasis(mfx), s = vector(nn+1));
  for(i=1,#bb, if(cv[i]!=0, my(c = mfcoefs(bb[i], nn)); s += cv[i]*c));
  s;
}

vec2ser(v) = sum(n=1,#v, v[n]*x^(n-1)) + O(x^#v);

/* D^{-1} and D^{-2} on a q-series with zero constant term (positive part) */
Dinv(s)  = my(n=serprec(s,x)); sum(k=1,n-1, polcoeff(s,k)/k * x^k) + O(x^n);
Dinv2(s) = my(n=serprec(s,x)); sum(k=1,n-1, polcoeff(s,k)/k^2 * x^k) + O(x^n);

/* fit a recurrence  sum_{j=0}^{r} p_j(n) a_{n+j} = 0,  deg p_j <= dg,
   to the coefficient vector av (av[n+1] = a_n).  Returns the coefficient
   matrix rows (p_j) or 0 if none. */
{fitrec(av, r, dg) =
  my(nvars = (r+1)*(dg+1), rows = List(), nmax = #av-1-r);
  my(nstart = 1);
  for(n=nstart, nmax,
    my(row = vector(nvars));
    for(j=0,r, for(e=0,dg, row[j*(dg+1)+e+1] = n^e * av[n+j+1]));
    listput(rows, row);
  );
  my(mm = Mat(Col(Vec(rows)))); mm = matconcat(Vec(rows)~);
  my(ker = matker(mm));
  if(matsize(ker)[2]==0, return(0));
  ker;
}

/* characteristic (singular) polynomial from a fitted recurrence kernel vector:
   leading behaviour a_n ~ lam^n  =>  sum_j (coeff of n^dg in p_j) lam^j = 0 */
{charpolyof(kv, r, dg) =
  my(p = 0);
  for(j=0,r, p += kv[j*(dg+1)+dg+1] * x^j);
  p;
}
