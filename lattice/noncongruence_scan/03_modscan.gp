/* Modular scan: rows a_n = lam^n [t^n] sqrt(F) for eta-quotient pairs (t,F).
 *   t : weight 0, ord_inf = 1        (parameter / Hauptmodul)
 *   F : weight 2, holomorphic, F(0)=1
 * g = sqrt(F) is a weight-1 form on an index-<=2 subgroup; by the classical
 * Sym^2 fact it satisfies a SECOND-order ODE in t, so the row is Apery-like
 * with k=2 (ROOT_ROWS Thm R3).  score = log(1/|lam_2|) - 2, lam_2 = lam*|1/t_2|.
 *
 * Stage 1 (this file): mod-p fit of the minimal recurrence of the UNSCALED row
 * c_n = [t^n]sqrt(F) (in Z[1/2]); lift the leading symbol; char roots; and the
 * 2-adic lambda from a 2-adic run.  Output: level, t, F, lambda, order, degree,
 * |lam1|, |lam2|, score.
 */
default(parisize, 2000000000);

MQ   = 60;          \\ q-precision
NROW = 52;          \\ number of row terms used
P1   = 2305843009213693951;   \\ 2^61-1
P2   = 1152921504606846883;   \\ another 60-bit prime
MAXO = 4;  MAXD = 3;

\\ ---- eta products mod m -------------------------------------------------
etaprod(d, M) = { my(s = 1 + O('q^M)); forstep(n=d, M-1, d, s *= (1 - 'q^n)); s; };

\\ build the vector of E_d for the divisors
buildE(divs, M) = vector(#divs, i, etaprod(divs[i], M));

\\ eta quotient with exponents r and q-order o
etaq(E, r, o, M) = { my(s = 'q^o + O('q^M)); for(i=1,#r, if(r[i], s *= E[i]^r[i])); s; };

\\ ---- fit minimal recurrence over a field --------------------------------
fitmin(A, maxo, maxd) = {
  my(L = #A - 1);   \\ A[1..L+1] = a_0..a_L
  for(o = 1, maxo,
    for(d = 0, maxd,
      my(nc = (o+1)*(d+1), rows = List(), r);
      for(n = o, L,
        r = vector(nc);
        for(j = 0, o, for(e = 0, d, r[j*(d+1)+e+1] = n^e * A[n+1-j]));
        listput(rows, r));
      if(#rows >= nc + 6,
        my(Mx = matconcat(Vec(rows)~), K = matker(Mx));
        if(#K == 1, return([o, d, K[,1]])))));
  0;
};

\\ characteristic roots from the fitted vector v (coefficients of n^d)
charpol_from(v, o, d) = {
  my(c = vector(o+1, j, v[(j-1)*(d+1)+d+1]));
  sum(j=1, o+1, c[j]*'L^(o+1-j));
};


\\ ---- two-stage fit: mod-p gate, then exact over Q ------------------------
PGATE = 2147483647;   \\ 2^31-1

fitminp(A, maxo, maxd) = {
  my(L = #A - 1, Ap = vector(#A, i, Mod(A[i], PGATE)));
  for(o = 1, maxo,
    for(d = 0, maxd,
      my(nc = (o+1)*(d+1), rows = List(), r);
      if(L - o + 1 < nc + 6, next);
      for(n = o, L,
        r = vector(nc);
        for(j = 0, o, for(e = 0, d, r[j*(d+1)+e+1] = Mod(n,PGATE)^e * Ap[n+1-j]));
        listput(rows, r));
      my(Mx = matconcat(Vec(rows)~), K = matker(Mx));
      if(#K == 1,
        my(rows2 = List(), r2);
        for(n = o, L,
          r2 = vector(nc);
          for(j = 0, o, for(e = 0, d, r2[j*(d+1)+e+1] = n^e * A[n+1-j]));
          listput(rows2, r2));
        my(M2 = matconcat(Vec(rows2)~), K2 = matker(M2));
        if(#K2 == 1, return([o, d, K2[,1]])))));
  0;
};
