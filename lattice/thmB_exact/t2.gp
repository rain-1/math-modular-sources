default(parisize, 2*10^9);
NT = 400;
frob(r,a,b,c,NT) = {
  my(A=vector(NT+1)); A[1]=1; A[2]=b;
  if(r==2, for(n=1,NT-1, A[n+2] = ((a*n^2+a*n+b)*A[n+1] - c*n^2*A[n])/(n+1)^2),
           for(n=1,NT-1, A[n+2] = ((2*n+1)*(a*n^2+a*n+b)*A[n+1] - c*n^3*A[n])/(n+1)^3));
  A;
}
build(r,a,b,c,NT) = {
  my(A=frob(r,a,b,c,NT), y0=Ser(A,'t,NT+1), P, K, g, qs, tq, Pq, Kq, Fq);
  if(r==2, P = 1 - a*'t + c*'t^2 + O('t^(NT+1)); K = P*y0^2,
           P = 1 - 2*a*'t + c*'t^2 + O('t^(NT+1)); K = y0*sqrt(P));
  g = intformal((1/K - 1)/'t);
  qs = 't*exp(g); tq = serreverse(qs);
  \\ Kq = theta_q t / t  as q-series
  Kq = ('t*deriv(tq,'t))/tq;
  if(r==2, Pq = 1 - a*tq + c*tq^2; Fq = sqrt(Kq/Pq),
           Pq = 1 - 2*a*tq + c*tq^2; Fq = Kq/sqrt(Pq));
  [tq, Fq];
}
gettime();
BB = build(3,17,5,1,NT);
print("gamma NT=",NT," time ",gettime()," ms");
print("t coeff 400: ", polcoeff(BB[1],400));
print("F coeff 1..8: ", vector(8,i,polcoeff(BB[2],i)));
print("F coeff 400: ", polcoeff(BB[2],400));
BB2 = build(2,17,6,72,NT); print("F-family time ",gettime()," ms  F[1..8]=",vector(8,i,polcoeff(BB2[2],i)));
quit;
