/* Morita p-adic Gamma: Gamma_p(n) = (-1)^n prod_{1<=j<n, p|/j} j  for n>=1 */
Gp(n,p) = { my(pr=1); for(j=1,n-1, if(j%p, pr*=j)); (-1)^n*pr; }

/* Taylor coefficients of  log( Gamma_p(1+x)/Gamma_p(1) )  at x=0,
   valid on |x|_p <= 1/q (q=p for p odd, 4 for p=2).
   Evaluate at x = q*m, m=1..K, solve Vandermonde.                */
gpTaylor(p, K, PR) =
{ my(q=if(p==2,4,p), V=matrix(K,K), rhs=vector(K));
  for(m=1,K,
    my(g = Gp(1+q*m,p)/Gp(1,p) + O(p^PR));
    rhs[m] = log(g);
    for(k=1,K, V[m,k] = (q*m)^k));
  matsolve(V, rhs~);
}
