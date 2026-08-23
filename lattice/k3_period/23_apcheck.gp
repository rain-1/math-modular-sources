default(parisizemax, 8000000000);
v = readvec("/home/ubuntu/code/math-modular-sources/lattice/k3_period/out/ap_ext.log");
mf32 = mfinit([32,3,-8],0); gg = mfeigenbasis(mf32)[1];
cf = mfcoefs(gg, 400);
bad = 0; nn = 0;
{
for(i=1,#v,
  my(l = v[i]);
  if(type(l) != "t_VEC" && type(l) != "t_INT", next));
}
s = Vec(readstr("/home/ubuntu/code/math-modular-sources/lattice/k3_period/out/ap_ext.log"));
{
for(i=1,#s,
  my(w = Vec(Str(s[i])), parts);
  parts = strsplit(s[i], " ");
  if(#parts != 2, next);
  my(p = eval(parts[1]), a = eval(parts[2]));
  if(!isprime(p), next);
  nn++;
  if(cf[p+1] != a, bad++; print("MISMATCH p=", p, " ours=", a, " form=", cf[p+1])));
}
print("primes checked: ", nn, "   mismatches: ", bad);
quit;
