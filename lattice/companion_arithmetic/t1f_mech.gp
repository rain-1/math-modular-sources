default(parisize,"4G");
read("lib.gp");
read("src.gp");
NT = 120;
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
  my(R=ROWS[i], r=R[2], cf=CMS[i], pq=param(R,NT), qs=pq[1], F=pq[2], AB=genrow(R,NT), B=AB[2]);
  my(pw=F, E=matrix(NT,NT), ok=1, tot=0, bad=List());
  for(m=1,NT, pw=pw*qs; for(n=m,NT, E[n,m]=polcoeff(pw,n,'t)));
  for(n=3,NT, forprime(p=n\2+1, n,
    my(pred = max(0, r - if(cf(p)*E[n,p]==0, 10^6, valuation(cf(p)*E[n,p],p))), act = max(0,-valuation(B[n+1],p)));
    tot++;
    if(pred!=act, ok=0; if(#bad<5, listput(bad,[n,p,pred,act])))));
  print(R[1], " r=", r, " : large-prime mechanism v_p(den b_n) = max(0, r - v_p(c(p) e_{n,p})) for n/2<p<=n, 3<=n<=", NT, " : ", if(ok, Str("HOLDS (", tot, " cases)"), Str("FAILS, first: ", Vec(bad)))));
}
quit;
