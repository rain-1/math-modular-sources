/* mono.gp -- generic numerical analytic continuation of an inhomogeneous linear
   ODE  sum_j P_j(x) y^{(j)} = R(x)  by adaptive Taylor stepping.
   Rows of type (R2): (n+1)^2 A_{n+1} = (a n^2+a n+b)A_n - c n^2 A_{n-1}
      operator  x(1-ax+cx^2) y'' + (1-2ax+3cx^2) y' + (cx-b) y = R,  R = 1 for B.
   Rows of type (R3): (n+1)^3 A_{n+1} = (2n+1)(a n^2+a n+b)A_n - n(c n^2+d)A_{n-1}
      operator  x^2 P_3 y''' + 3x(1-3ax+2cx^2) y'' + (1-(6a+2b)x+7cx^2) y'
                 + (-b+cx) y = R,  P_3 = 1-2ax+cx^2   (d = 0 rows only).       */

MT = 260;                       /* Taylor terms per step */
fall(m, j) = prod(l=0, j-1, m-l);
/* shift a polynomial (coefficient vector, ascending) to x0 */
shft(v, x0) = my(n=#v); vector(n, i, sum(k=i, n, v[k]*binomial(k-1,i-1)*x0^(k-i)));
/* one Taylor step: PL = list of shifted coefficient vectors P_0..P_ord,
   iv = [y, y', ..., y^{(ord-1)}] at x0, Rv = RHS Taylor coefs at x0, h = step */
stp(PL, iv, Rv, h) = {
  my(ord = #PL-1, NB = MT+ord, cc = vector(NB+1), acc, top, res);
  for(j=1, ord, cc[j] = iv[j]/(j-1)!);
  for(n=0, MT,
    acc = 0;
    for(j=0, ord, my(Pv = PL[j+1]);
      for(i=0, #Pv-1,
        if(j==ord && i==0, next());
        my(idx = n-i+j);
        if(idx < 0 || idx > NB-1, next());
        if(Pv[i+1]==0, next());
        acc += Pv[i+1]*fall(idx, j)*cc[idx+1]));
    top = PL[ord+1][1]*fall(n+ord, ord);
    cc[n+ord+1] = (Rv[n+1] - acc)/top);
  res = vector(ord);
  for(j=0, ord-1, my(s=0); forstep(k=NB, j, -1, s = s*h + fall(k,j)*cc[k+1]); res[j+1] = s);
  res;
}
/* radius to nearest singularity */
radf(SG, x) = vecmin(vector(#SG, i, abs(x-SG[i])));
/* walk along a straight segment, updating the local data */
walk(PC, SG, RC, iv, x0, x1) = {
  my(cur = x0, st = iv, rem, h, dir, PL, Rv, ord = #PC-1);
  while(abs(x1-cur) > 1e-70,
    rem = abs(x1-cur); dir = (x1-cur)/rem;
    h = min(rem, 0.35*radf(SG,cur));
    PL = vector(#PC, i, shft(PC[i], cur));
    Rv = RC(cur);
    st = stp(PL, st, Rv, h*dir);
    cur = cur + h*dir);
  st;
}
/* follow a polyline */
runp(PC, SG, RC, iv, verts) = my(st=iv); for(j=2,#verts, st = walk(PC,SG,RC,st,verts[j-1],verts[j])); st;
/* evaluate a power series and its derivatives at x */
evs(cf, x, ord) = vector(ord, j, my(s=0); forstep(k=#cf-1, j-1, -1, s = s*x + fall(k,j-1)*cf[k+1]); s);
