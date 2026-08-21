/* Fit a linear ODE  sum_{i<=rr} p_i(t) theta^i y = 0, deg p_i <= dm, to series y (a polynomial in x). */
{
odefit(y, rr, dm) =
  my(NT, terms, lab, M);
  NT = poldegree(y);
  terms = List(); lab = List();
  for(i=0,rr, for(j=0,dm,
     listput(terms, x^j*sum(m=0,NT, m^i*polcoeff(y,m)*x^m));
     listput(lab,[i,j])));
  terms = Vec(terms); lab = Vec(lab);
  M = matrix(NT-dm-1, #terms, a, b, polcoeff(terms[b], a+dm-1));
  [matker(M), lab];
}
{
showrel(K, lab) =
  my(v);
  for(c=1, matsize(K)[2], v = K[,c];
    print("  relation ", c, ":");
    for(b=1,#lab, if(v[b]!=0, print("    ", v[b], " * t^", lab[b][2], " * theta^", lab[b][1]))));
}
