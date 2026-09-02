\\ 21_laws6000.gp -- the mod-2 / mod-3 laws and the zero set, re-verified at n <= 6000.
read("lib.gp");
G  = [read("20b_gamma_s7.txt"), read("20b_gamma_s10.txt"), read("20b_gamma_s18.txt")];
N  = #G[1];
print("N = ", N);
print();
{
my(bad);
bad = select(n->((G[1][n]%2==1) != (gcd(n,14)==1)), vector(N,i,i));
print("  L1 s7 : gamma odd <=> gcd(n,14)=1                : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->((G[2][n]%2==1) != (n%4!=0 && n%5!=0)), vector(N,i,i));
print("  L2 s10: gamma odd <=> 4 nmid n and 5 nmid n      : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->((G[3][n]%2==1) != (n%4!=0 && n%3!=0)), vector(N,i,i));
print("  L3 s18: gamma odd <=> 4 nmid n and 3 nmid n      : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->(lift(Mod((-1)^(n-1)*G[3][n],3)) != if(n%3==0,2,1)), vector(N,i,i));
print("  L4 s18: (-1)^(n-1)gamma(n) = 1 (3 nmid n), -1 (3|n) mod 3 : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->(G[3][n]*(-1)^(n-1)<=0), vector(N,i,i));
print("  L5 s18: sign(gamma(n)) = (-1)^(n-1)              : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->((G[1][n]==0) != (n%7==0)), vector(N,i,i));
print("  L6 s7 : gamma(n)=0 <=> 7|n                       : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->((G[2][n]==0) != (n%5==0)), vector(N,i,i));
print("  L7 s10: gamma(n)=0 <=> 5|n                       : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->(G[3][n]==0), vector(N,i,i));
print("  L8 s18: gamma(n) never 0                         : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->(n%7!=0 && valuation(G[1][n],2)!=if(n%2==0,1,0)), vector(N,i,i));
print("  L9 s7 : v_2(gamma(n)) = [n even] for 7 nmid n    : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->(n%5!=0 && valuation(G[2][n],2) < valuation(n,2)-1), vector(N,i,i));
print("  L10 s10: v_2(gamma(n)) >= v_2(n)-1               : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->(n%7!=0 && valuation(G[1][n],3) < valuation(n,3)-1), vector(N,i,i));
print("  L11 s7 : v_3(gamma(n)) >= v_3(n)-1               : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->(valuation(G[3][n],3)!=0), vector(N,i,i));
print("  L12 s18: 3 nmid gamma(n)                         : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->(valuation(G[3][n],2)!=valuation(G[3][if(n>24,n-24,n+24)],2)), vector(N,i,i));
print("  L13 s18: v_2(gamma(n)) depends only on n mod 24  : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
}
print();
print("=== periodicity of gamma mod p^k at N=",N," (recheck) ===");
{
my(PP=[2,2,2,3], EE=[1,2,3,1]);
for(k=1,3,
  for(i=1,4,
    my(p=PP[i], e=EE[i], m=p^e, T=0);
    for(t=1,1000, my(ok=1);
      for(n=1,N-t, if((G[k][n]-G[k][n+t])%m!=0, ok=0; break));
      if(ok, T=t; break));
    print("  ",NAM[k]," mod ",m,": ", if(T>0, concat("period ",T), "no period <= 1000"));
  );
);
}
print();
print("=== s7 mod 4 restricted to n odd; s10/s18 mod 8 restricted ===");
{
my(bad, T);
T=0; for(t=1,600, my(ok=1); for(n=1,N-t, if(n%2==1 && (n+t)%2==1 && (G[1][n]-G[1][n+t])%4!=0, ok=0;break)); if(ok,T=t;break));
print("  s7  mod 4 on odd n: ", if(T>0,concat("period ",T),"no period <= 600"));
T=0; for(t=1,600, my(ok=1); for(n=1,N-t, if(n%2==1 && (n+t)%2==1 && (G[2][n]-G[2][n+t])%8!=0, ok=0;break)); if(ok,T=t;break));
print("  s10 mod 8 on odd n: ", if(T>0,concat("period ",T),"no period <= 600"));
T=0; for(t=1,600, my(ok=1); for(n=1,N-t, if(n%2==1 && (n+t)%2==1 && (G[3][n]-G[3][n+t])%8!=0, ok=0;break)); if(ok,T=t;break));
print("  s18 mod 8 on odd n: ", if(T>0,concat("period ",T),"no period <= 600"));
}
quit;
