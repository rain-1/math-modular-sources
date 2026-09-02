\\ 09_gamma_probe.gp -- structure probes on gamma: Dirichlet inverse, growth, and
\\ cross-row comparisons.
read("lib.gp");
M = 200;
{ dirinv(f) = my(N=#f, g=vector(N)); g[1] = 1/f[1];
  for(n=2,N, my(s=0); fordiv(n,d, if(d<n, s += f[n/d]*g[d])); g[n] = -s/f[1]);
  g; }
GG = vector(3);
{
for(k=1,3,
  my(cp=CPvec(k,M), b, g);
  b = Bvec(k,cp); g = vector(M,n,b[n]/n^2); GG[k]=g;
  print("row ",NAM[k]);
  print("  gamma(1..20)      = ", vector(20,i,g[i]));
  print("  gamma^{-1}(1..20) = ", vector(20,i,dirinv(g)[i]));
);
}
\\ growth: gamma(n) * n^2 / R^n  should oscillate boundedly; check |gamma(n)|^{1/n} -> R
{
my(Rp=[exp(Pi*sqrt(3)/7), exp(Pi/5), exp(Pi/3)]);
default(realprecision,30);
for(k=1,3,
  print("row ",NAM[k],"  predicted R = ",Rp[k]);
  for(j=1,4, my(n=M-4+j); print("    n=",n,"  |gamma(n)|^(1/n) = ", abs(GG[k][n]*1.)^(1/n)));
);
}
\\ cross-row: is gamma of one row equal to beta/c'/gamma of another (lindep on 30 terms)?
{
my(v);
for(k=1,3, for(l=1,3,
  if(k==l, next);
  v = lindep(concat(vector(12,i,GG[k][i]*1.), vector(12,i,GG[l][i]*1.)));
  print("lindep(gamma_",NAM[k],", gamma_",NAM[l],") -> nontrivial? ", if(#v>0, "see", "-"));
));
}
quit;
