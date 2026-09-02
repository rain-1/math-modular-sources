\\ 21_gammap.gp -- gamma(p) mod p, all rows, p <= 43 (and p <= 200 for structure hunting).
read("lib.gp");
G = [read("gamma_s7.txt"), read("gamma_s10.txt"), read("gamma_s18.txt")];
M = #G[1];
print("data length M = ", M);
PR = select(x->isprime(x), vector(M,i,i));
print("primes up to ", M, ": ", #PR);
print();
print("== gamma(p) mod p ==");
{
for(k=1,3,
  print("row ", NAM[k]);
  print("  p            : ", vector(14,i,PR[i]));
  print("  gamma(p)     : ", vector(14,i,G[k][PR[i]]));
  print("  gamma(p) mod p: ", vector(14,i,lift(Mod(G[k][PR[i]],PR[i]))));
  print("  psi(p)       : ", vector(14,i,psival(k,PR[i])));
  print("  zero?        : ", select(p->G[k][p]%p==0, vector(#PR,i,PR[i])));
  print();
);
}
print("== full list p <= 600, gamma(p) mod p ==");
{
for(k=1,3,
  print("row ", NAM[k], ":");
  for(i=1,#PR, my(p=PR[i]); print1("[",p,",",lift(Mod(G[k][p],p)),"] "));
  print(); print();
);
}
quit;
