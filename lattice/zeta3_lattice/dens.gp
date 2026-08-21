read("/home/ubuntu/code/math-modular-sources/lattice/zeta3_lattice/rows.gp");
D=rowsI(10,4,64); T=rowsI(12,4,16);
print("true denominator growth  k_eff = log(den b_n)/n   (d_n^3 would give 3)");
{my(N=400);for(j=1,10,my(n=40*j);
  my(fD=n!^3, gD=gcd(D[2][n+1],fD), dD=fD/gD, fT=n!^3, gT=gcd(T[2][n+1],fT), dT=fT/gT);
  print("n=",n,"  Domb k_eff=",log(dD*1.)/n,"   T k_eff=",log(dT*1.)/n,
        "   [log d_n^3/n=",3*log(lcm(vector(n,i,i))*1.)/n,"]"));}
print();
print("does d_n^3 b_n in Z hold to n=400? (both rows)");
{my(d=1,okD=1,okT=1);for(n=1,400,d=lcm(d,n);
  if((d^3*D[2][n+1])%(n!^3)!=0,okD=0);if((d^3*T[2][n+1])%(n!^3)!=0,okT=0));
 print("  Domb ",okD,"  T ",okT);}
print();
print("is 2^j d_n^3 b_n integral for smaller exponent? test k=2 with extra lcm-free factor:");
{my(d=1);for(n=1,200,d=lcm(d,n));
 print("  (informational only)");}
\q
