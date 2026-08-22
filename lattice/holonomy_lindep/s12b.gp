default(parisizemax, 8G);
default(realprecision, 40);
read("/home/ubuntu/code/math-modular-sources/lattice/mum_survey/ops.gp");
chipoly(op) = {
  my(dz = op[3], Ps = op[4], dmax = 0);
  for(i=1, #Ps, dmax = max(dmax, poldegree(Ps[i], 'X)));
  sum(i=0, dz, polcoeff(Ps[i+1], dmax, 'X) * 'x^(dz-i));
};
{
for(t=1, #OPS,
  my(op = OPS[t]);
  if(op[3] < 3, next);
  my(cp = chipoly(op));
  my(rr = polroots(cp), m = vector(#rr, i, abs(rr[i])));
  my(idx = vecsort(m,,5));
  my(l1 = rr[idx[1]], a1 = m[idx[1]], l2 = rr[idx[2]], a2 = m[idx[2]]);
  my(sim = if(abs(a2-a1) < 1e-20*a1, 0, 1));
  my(r1 = if(abs(imag(l1)) < 1e-20*max(1,a1), 1, 0));
  my(r2 = if(abs(imag(l2)) < 1e-20*max(1,a2), 1, 0));
  my(rat = a2/a1);
  my(sc = -log(rat));
  print(t,",",op[1],",",op[2],",",op[3],",",a1,",",a2,",",rat,",",sim,",",r1,",",r2,",",sc);
);
}
quit;
