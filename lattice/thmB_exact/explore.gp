default(realprecision, 120);
NT = 50;

frob(r,a,b,c,NT) = {
  my(A=vector(NT+1)); A[1]=1; A[2]=b;
  if(r==2,
    for(n=1,NT-1, A[n+2] = ((a*n^2+a*n+b)*A[n+1] - c*n^2*A[n])/(n+1)^2),
    for(n=1,NT-1, A[n+2] = ((2*n+1)*(a*n^2+a*n+b)*A[n+1] - c*n^3*A[n])/(n+1)^3));
  A;
}

\\ returns [tq, Fq] as q-series
build(r,a,b,c,NT) = {
  my(A=frob(r,a,b,c,NT), y0=Ser(A,'t,NT+1), P, K, iv, g, qs, tq, Fq);
  if(r==2, P = 1 - a*'t + c*'t^2 + O('t^(NT+1)); K = P*y0^2,
           P = 1 - 2*a*'t + c*'t^2 + O('t^(NT+1)); K = y0*sqrt(P));
  iv = 1/K - 1;
  g = intformal(iv/'t);        \\ integral of (1/K-1) dt/t
  qs = 't*exp(g);
  tq = serreverse(qs);
  Fq = subst(y0, 't, tq);
  [tq, Fq];
}

test(nm,r,a,b,c) = {
  my(BB=build(r,a,b,c,NT), tq=BB[1], Fq=BB[2]);
  print(nm,"  t(q) = ", tq + O(q^9));
  print(nm,"  F(q) = ", Fq + O(q^9));
}
test("gamma",3,17,5,1);
test("alpha",3,10,4,64);
test("C",2,10,3,9);
test("D",2,11,3,-1);
quit;
