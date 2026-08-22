default(parisizemax, 8G);
default(realprecision, 40);
read("/home/ubuntu/code/math-modular-sources/lattice/mum_survey/ops.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/sc_rows.gp");
chipoly(op) = {
  my(dz = op[3], Ps = op[4], dmax = 0);
  for(i=1, #Ps, dmax = max(dmax, poldegree(Ps[i], 'X)));
  sum(i=0, dz, polcoeff(Ps[i+1], dmax, 'X) * 'x^(dz-i));
};
\\ recurrence coefficients Cf[i+1](n) = P_i(n-i)
mkcf(op) = { my(dz=op[3], Ps=op[4]); vector(dz+1, i, subst(Ps[i], 'X, 'n-(i-1))); };
NDEN = 200;
{
for(t=1, #OPS,
  my(op = OPS[t]);
  if(op[3] < 3, next);
  my(cp = chipoly(op));
  my(rr = polroots(cp), m = vector(#rr, i, abs(rr[i])));
  my(idx = vecsort(m,,5));
  my(a1 = m[idx[1]], a2 = m[idx[2]]);
  if(abs(a2-a1) < 1e-20*a1, next);   \\ non-simple max modulus
  my(rat = a2/a1);
  if(rat > 0.97, next);
  my(dz = op[3], cf = mkcf(op));
  my(ks = vector(dz-1));
  for(j=1, dz-1,
    my(x = compan(cf, dz, j, NDEN));
    ks[j] = denexp(x, NDEN);
  );
  my(sc = -log(rat));
  print(t,",",op[1],",",op[2],",",dz,",",a1,",",a2,",",rat,",",sc,",",ks);
);
}
quit;
