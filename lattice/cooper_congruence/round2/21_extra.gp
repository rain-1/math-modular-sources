\\ 21_extra.gp -- follow-ups: sign of gamma_s18, p | gamma(p), refined periodicity.
read("lib.gp");
G  = [read("22_gamma_s7.txt"), read("22_gamma_s10.txt"), read("22_gamma_s18.txt")];
N  = #G[1];
print("N = ", N);
print();
print("=== sign pattern ===");
{
my(bad);
bad = select(n->(G[3][n]*(-1)^(n-1) <= 0), vector(N,i,i));
print("  s18: sign(gamma(n)) = (-1)^(n-1) for all n<=",N,"? ", if(#bad==0,"YES",concat("NO, first at ",bad[1])));
for(k=1,2,
  bad = select(n->(n%LEV[k]!=0 && G[k][n]==0), vector(N,i,i));
  print("  ",NAM[k],": gamma(n)=0 only at multiples of the special prime? ", if(#bad==0,"YES","NO"));
);
}
print();
print("=== primes p <= N with p | gamma(p) ===");
{
my(np = #select(x->isprime(x), vector(N,i,i)), hs = sum(i=1,N, if(isprime(i),1.0/i,0)));
for(k=1,3,
  my(S = select(p->(isprime(p) && G[k][p]%p==0), vector(N,i,i)));
  print("  ",NAM[k],": ", S, "   (",#S," of ",np," primes; heuristic sum 1/p = ", strprintf("%.2f",hs)," )");
);
}
print();
print("=== is gamma(n) mod 2 multiplicative in n ? ===");
{
for(k=1,3,
  my(f=0,t=0);
  for(a=2,60, for(b=2,N\a, if(gcd(a,b)==1, t++;
     if((G[k][a*b]-G[k][a]*G[k][b])%2!=0, f++))));
  print("  ",NAM[k],": coprime pairs (a<=60) failing gamma(ab)=gamma(a)gamma(b) mod 2 : ",f,"/",t);
);
}
print();
print("=== s18: gamma mod 9, 27 restricted to n coprime to 3 : periodic ? ===");
{
my(ML=[9,27]);
for(mi=1,2,
  my(m=ML[mi], found=0);
  for(T=1,300,
    my(ok=1);
    for(n=1,N-T, if(n%3!=0 && (n+T)%3!=0 && (G[3][n]-G[3][n+T])%m!=0, ok=0; break));
    if(ok, print("  s18 mod ",m," on n coprime to 3: period ",T); found=1; break);
  );
  if(found==0, print("  s18 mod ",m," on n coprime to 3: no period <= 300"));
);
}
print();
print("=== distribution of gamma(p) mod p, 4 bins ===");
{
for(k=1,3,
  my(bins=vector(4), tot=0);
  forprime(p=11,N, my(r=lift(Mod(G[k][p],p)), b=1+(4*r)\p); if(b>4,b=4); bins[b]++; tot++);
  print("  ",NAM[k],": bins ",bins," of ",tot,"  chi^2 = ",
        strprintf("%.2f", sum(i=1,4,(bins[i]-tot/4.0)^2/(tot/4.0))), "   [3 dof, 5% crit 7.81]");
);
}
quit;
