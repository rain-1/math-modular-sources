/* Zudilin's Catalan row (arXiv:math/0201024, Thm 1) and its identification
   with Beukers' Theta-Pade denominators at the moving point x = -n+1/2. */
pz(n) = 20*n^2-8*n+1;
qz(n) = 3520*n^6+5632*n^5+2064*n^4-384*n^3-156*n^2+16*n+7;
zudrow(N) =
{ my(u=vector(N+2), v=vector(N+2));
  u[1]=1; u[2]=7/4; v[1]=0; v[2]=13/8;
  for(n=1,N,
    u[n+2] = (qz(n)*u[n+1] + (2*n-1)^2*(2*n)^2*pz(n+1)*u[n]) / ((2*n+1)^2*(2*n+2)^2*pz(n));
    v[n+2] = (qz(n)*v[n+1] + (2*n-1)^2*(2*n)^2*pz(n+1)*v[n]) / ((2*n+1)^2*(2*n+2)^2*pz(n)));
  [u,v];
}
/* Beukers Pade denominator q_n(x) = sum_k C(n,k) C(-x,k) C(k-x,k) */
bq(n,xx) = sum(k=0,n, binomial(n,k)*binomial(xx,k)*binomial(xx+k,k));  /* with xx = -x */
