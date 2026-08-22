/* Elliptic-surface test.  See HERFURTNER_CLASSIFICATION.md sec. 4.
 * jtest(M,j1,j2,A,B,C) -> [h, mu] if J(h*tau) is a rational function of t of
 * degree mu <= mumax (so the projective monodromy of L is conjugate into
 * PSL_2(Z) and L is the Picard-Fuchs operator of an elliptic surface over
 * P^1_t, with a cusp of width h at t = 0), and 0 otherwise.
 */
default(realprecision,38);
NS = 56;

jser(N) = { my(E4, Dl);
  E4 = 1 + 240*sum(n=1,N, sigma(n,3)*x^n) + O(x^(N+1));
  Dl = x*prod(n=1,N, (1-x^n+O(x^(N+1)))^24);
  E4^3/Dl; }

JS = jser(NS);
JJ = x*JS;          /* power series, JJ = 1 + 744 x + 196884 x^2 + ... */

nome(M,j1,j2,A,B,C,N) = {
  my(be,de,ep,ze,u,g,y0,gg);
  be = A*(2*M-j1-j2)/(2*M);
  de = C*M^2; ep = -C*M*(j1+j2); ze = C*j1*j2;
  u = vector(N+2); g = vector(N+2); u[1]=1; g[1]=0;
  for(n=1, N+1,
    my(P0 = A*(n-1)^2 + be*(n-1) + B, Q0 = de*(n-1)^2 + ep*(n-1) + ze,
       dP = 2*A*(n-1) + be,           dQ = 2*de*(n-1) + ep);
    u[n+1] = (P0*u[n] - if(n>=2, Q0*u[n-1], 0))/n^2;
    g[n+1] = (dP*u[n] + P0*g[n] - if(n>=2, dQ*u[n-1] + Q0*g[n-1], 0) - 2*n*u[n+1])/n^2;
  );
  y0 = sum(n=0,N, u[n+1]*x^n) + O(x^(N+1));
  gg = sum(n=0,N, g[n+1]*x^n) + O(x^(N+1));
  x*exp(gg/y0);
}

/* W a power series (valuation 0) with W = x^h * Jh(t); look for polys U,V of
   degree <= d with V*W = x^h*U. */
ratfit(W, h, d, N) = {
  my(neq, mat, ker, U, V, R);
  neq = 2*d + 6; if(neq > N-6, neq = N-6);
  mat = matrix(neq, 2*d+2);
  for(i=1, neq, my(e=i-1);
    for(k=0, d, mat[i,k+1]   = if(e-k>=0, polcoeff(W, e-k), 0));
    for(k=0, d, mat[i,d+2+k] = if(e == h+k, -1, 0)));
  ker = matker(mat);
  for(c=1, #ker,
    V = sum(k=0,d, ker[,c][k+1]*x^k);
    U = sum(k=0,d, ker[,c][d+2+k]*x^k);
    if(V != 0 && polcoeff(V,0) != 0 || V != 0,
      R = V*W - x^h*U;
      my(ok=1, emax = serprec(W,x)-3); if(emax > neq+8, emax = neq+8);
      for(e=0, emax, if(polcoeff(R,e) != 0, ok=0; break));
      if(ok, return([U,V]))));
  0;
}

GAMS = [1,2,3,4,5,6,8,9,12,16,18,24,25,27,32,36,48,54,64,100,125,128,216,500];

jtest(M,j1,j2,A,B,C, mumax=12, N=NS) =
{
  my(q, xq, W, fit, gl, qh, gam, U, V, gg);
  q = nome(M,j1,j2,A,B,C,N);
  xq = x/q;
  for(h = 1, mumax,
    qh = q^h;
    gl = List();
    for(i = 1, #GAMS,
      listput(gl,  GAMS[i]^h); listput(gl, -GAMS[i]^h);
      listput(gl,  1/GAMS[i]^h); listput(gl, -1/GAMS[i]^h)
    );
    for(gi = 1, #gl,
      gam = gl[gi];
      W = subst(JJ, x, gam*qh) * xq^h / gam;
      fit = ratfit(W, h, mumax, N);
      if(fit != 0,
        U = fit[1]; V = fit[2]; gg = gcd(U,V); U = U/gg; V = V/gg;
        return([h, max(poldegree(U),poldegree(V)), gam, U, V])
      )
    )
  );
  0;
}
