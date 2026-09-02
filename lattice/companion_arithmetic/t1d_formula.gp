default(parisize,"4G");
read("lib.gp");
read("src.gp");
NT = 60;
/* q as a series in t, and F = y_0(t); exact rational series */
{ param(R,NT) =
  my(r=R[2], a=R[3], b=R[4], c=R[5], AB=genrow(R,NT+2), A=AB[1], y0, P, K, g, qs);
  y0 = Ser(vector(NT+2,j,A[j]), 't, NT+2);
  P = if(r==2, 1 - a*'t + c*'t^2 + O('t^(NT+2)), 1 - 2*a*'t + c*'t^2 + O('t^(NT+2)));
  K = if(r==2, P*y0^2, y0*sqrt(P));
  g = intformal((1/K - 1)/'t);
  qs = 't*exp(g);
  [qs, y0];
}
{ for(i=1,12,
  my(R=ROWS[i], r=R[2], cf=CMS[i], pq=param(R,NT), qs=pq[1], F=pq[2]);
  my(AB=genrow(R,NT), B=AB[2], E=matrix(NT,NT), pw=F, allint=1, badq=0, ok=1, firstbad=0, S);
  /* q(t) integrality */
  for(m=1,NT, if(denominator(polcoeff(qs,m,'t))!=1, badq=m; break));
  for(m=1,NT, pw = pw*qs;  /* pw = F * q^m */
    for(n=m,NT, E[n,m] = polcoeff(pw,n,'t); if(denominator(E[n,m])!=1, allint=0)));
  for(n=1,NT, S = sum(m=1,n, cf(m)*E[n,m]/m^r);
    if(S != B[n+1], ok=0; if(firstbad==0, firstbad=n)));
  print(R[1], " r=", r, " : q(t) in Z[[t]] to t^", NT, ": ", if(badq==0,"YES",Str("NO, first non-integral coeff at t^",badq)), " | e_{n,m} in Z for all 1<=m<=n<=", NT, ": ", if(allint,"YES","NO"), " | b_n = sum_m c(m) m^-r e_{n,m} for all n<=", NT, ": ", if(ok,"YES",Str("NO, first n=",firstbad))));
}
quit;
