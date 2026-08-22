/* Elliptic-surface test for a second-order Apery-like row.
 *
 *  row:  (n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1},  P = A n^2 + be n + B,
 *        Q = C (M n - j1)(M n - j2).
 *  L = th^2 - t P(th) + t^2 Q(th+1);  Frobenius basis y0 = sum u_n t^n,
 *  y1 = y0 log t + g,  g = sum g_n t^n.  Nome q = t exp(g/y0).
 *
 *  The projective monodromy of L is a Fuchsian group G commensurable with
 *  PSL_2(Z).  G is CONJUGATE INTO PSL_2(Z) (equivalently: L is the
 *  Picard-Fuchs operator of an elliptic surface over P^1_t) iff, for some
 *  positive integer h (the width of the cusp t=0), the function
 *  J(h*tau) is a RATIONAL function of t.  Its degree is then the index mu.
 *
 *  jtest(M,j1,j2,A,B,C, mumax, prec) returns [h, degnum, degden] on success,
 *  0 on failure.
 */

NS = 80;            \\ series precision in q / t

/* j-function q-expansion, 1/q + 744 + ... */
jser(N) = {
  my(E4, Dl);
  E4 = 1 + 240*sum(n=1,N, sigma(n,3)*x^n) + O(x^(N+1));
  Dl = x*prod(n=1,N, (1-x^n+O(x^(N+1)))^24);
  E4^3/Dl;
}

JS = 0;

nome(M,j1,j2,A,B,C,N) = {
  my(be, de, ep, ze, u, g, y0, gg, q);
  be = A*(2*M-j1-j2)/(2*M);
  de = C*M^2; ep = -C*M*(j1+j2); ze = C*j1*j2;
  Pv(n) = A*n^2 + be*n + B;
  Qv(n) = de*n^2 + ep*n + ze;
  dPv(n) = 2*A*n + be;
  dQv(n) = 2*de*n + ep;
  u = vector(N+1); g = vector(N+1);
  u[1] = 1; g[1] = 0;                       \\ index n+1 <-> u_n
  for(n=1, N,
    my(t1 = Pv(n-1)*u[n], t2 = if(n>=2, Qv(n-1)*u[n-1], 0));
    u[n+1] = (t1 - t2)/n^2;
    my(s = dPv(n-1)*u[n] + Pv(n-1)*g[n]
         - if(n>=2, dQv(n-1)*u[n-1] + Qv(n-1)*g[n-1], 0) - 2*n*u[n+1]);
    g[n+1] = s/n^2;
  );
  y0 = sum(n=0,N, u[n+1]*x^n) + O(x^(N+1));
  gg = sum(n=0,N, g[n+1]*x^n) + O(x^(N+1));
  q = x*exp(gg/y0);
  q;
}

/* try to write the Laurent series F (in x) as U(x)/V(x), deg <= dmax */
ratfit(F, dmax, N) = {
  my(v, sh, Fs, rows, mat, ker, U, V);
  v = valuation(F, x);
  Fs = F;                              \\ Laurent series
  \\ unknowns: U_0..U_dmax, V_0..V_dmax ; equation V*F - U = 0
  \\ use coefficients of x^v .. x^(v+2*dmax+4)
  my(neq = 2*dmax+6);
  if(v + neq > N-2, neq = N-2-v);
  mat = matrix(neq, 2*dmax+2);
  for(i=1, neq,
    my(e = v + i - 1);
    for(jj=0, dmax, mat[i, jj+1] = polcoeff(Fs, e-jj));      \\ V_jj * F
    for(jj=0, dmax, mat[i, dmax+2+jj] = if(e==jj, -1, 0));   \\ -U_jj
  );
  ker = matker(mat);
  if(#ker == 0, return(0));
  \\ pick a kernel vector with V != 0
  for(c=1, #ker,
    V = sum(jj=0,dmax, ker[,c][jj+1]*x^jj);
    U = sum(jj=0,dmax, ker[,c][dmax+2+jj]*x^jj);
    if(V != 0,
      if(polcoeff(V*Fs - U, 0) === 0 || truncate(V*Fs - U + O(x^(N-4))) == 0,
        return([U, V]));
      \\ generic check: all computed coefficients vanish
      my(R = V*Fs - U);
      my(ok = 1);
      for(e = valuation(R,x), N-6, if(polcoeff(R,e) != 0, ok=0; break));
      if(ok, return([U,V]));
    );
  );
  0;
}

jtest(M,j1,j2,A,B,C, mumax=12, N=NS) = {
  my(q, tq, res, Jh, Jt);
  if(JS==0, JS = jser(N));
  q = nome(M,j1,j2,A,B,C,N);
  tq = serreverse(q);                  \\ t as a series in q
  res = List();
  for(h=1, mumax,
    Jh = subst(JS, x, x^h);            \\ J(h tau), series in q
    Jt = subst(Jh, x, q);              \\ hmm: need J as series in t
    );
  0;
}
