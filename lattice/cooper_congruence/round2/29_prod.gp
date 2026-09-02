\\ 29_prod.gp -- (4) the sequence e_n = beta(n)/n = n*gamma(n) and the q-product g,
\\ plus explicit verification of the closed-form mod-2 / mod-3 laws for gamma.
read("lib.gp");
G  = [read("22_gamma_s7.txt"), read("22_gamma_s10.txt"), read("22_gamma_s18.txt")];
CP = [read("22_cp_s7.txt"),    read("22_cp_s10.txt"),    read("22_cp_s18.txt")];
N  = #G[1];
print("N = ", N);

print();
print("=== closed-form mod-2 / mod-3 laws for gamma, verified n<=",N," ===");
{
my(bad);
bad = select(n->((G[1][n]%2==1) != (gcd(n,14)==1)), vector(N,i,i));
print("  s7 : gamma(n) odd <=> gcd(n,14)=1        : ", if(#bad==0,"HOLDS",concat("FAILS first at ",bad[1])));
bad = select(n->((G[2][n]%2==1) != (n%4!=0 && n%5!=0)), vector(N,i,i));
print("  s10: gamma(n) odd <=> 4 nmid n and 5 nmid n : ", if(#bad==0,"HOLDS",concat("FAILS first at ",bad[1])));
bad = select(n->((G[3][n]%2==1) != (n%4!=0 && n%3!=0)), vector(N,i,i));
print("  s18: gamma(n) odd <=> 4 nmid n and 3 nmid n : ", if(#bad==0,"HOLDS",concat("FAILS first at ",bad[1])));
bad = select(n->(lift(Mod((-1)^(n-1)*G[3][n],3)) != if(n%3==0,2,1)), vector(N,i,i));
print("  s18: (-1)^(n-1) gamma(n) = 1 mod 3 if 3 nmid n, = -1 mod 3 if 3|n : ", if(#bad==0,"HOLDS",concat("FAILS first at ",bad[1])));
bad = select(n->(n%7!=0 && valuation(G[1][n],2)!=(if(n%2==0,1,0))), vector(N,i,i));
print("  s7 : v_2(gamma(n)) = 0 (n odd), = 1 (n even), for 7 nmid n : ", if(#bad==0,"HOLDS",concat("FAILS first at ",bad[1])));
bad = select(n->(valuation(G[2][n]+if(G[2][n]==0,1,0),2) < valuation(n,2)-1 && n%5!=0), vector(N,i,i));
print("  s10: v_2(gamma(n)) >= v_2(n)-1           : ", if(#bad==0,"HOLDS",concat("FAILS first at ",bad[1])));
bad = select(n->(n%7!=0 && valuation(G[1][n],3) < valuation(n,3)-1), vector(N,i,i));
print("  s7 : v_3(gamma(n)) >= v_3(n)-1           : ", if(#bad==0,"HOLDS",concat("FAILS first at ",bad[1])));
}
print();
print("  s18: v_2(gamma(n)) as a function of n mod 24, n=1..24 : ", vector(24,n,valuation(G[3][n],2)));

print();
print("=== (4) e_n = n*gamma(n) : periodicity of e_n mod p^k ===");
{
for(k=1,3,
  for(ip=1,3,
    my(p=prime(ip));
    for(e=1,3,
      my(m=p^e, T=0);
      for(t=1,600,
        my(ok=1);
        for(n=1,N-t, if((n*G[k][n]-(n+t)*G[k][n+t])%m!=0, ok=0; break));
        if(ok, T=t; break);
      );
      if(T>0, print("  ",NAM[k]," e_n mod ",m,": period ",T), print("  ",NAM[k]," e_n mod ",m,": no period <= 600"); break);
    );
  );
);
}

print();
print("=== the q-product g = prod (1-q^n)^{-e_n} = exp(sum c'(m)q^m/m), s7 & s10 ===");
{
my(M=1200);
for(k=1,2,
  my(ser, gg, cf, bad, T);
  ser = sum(m=1,M, CP[k][m]*'q^m/m) + O('q^(M+1));
  gg  = exp(ser);
  cf  = Vec(gg, -(M+1));
  bad = select(i->denominator(cf[i])!=1, vector(M+1,i,i));
  print("  ",NAM[k],": g integral to q^",M,"? ", if(#bad==0,"YES",concat("NO at ",bad[1]-1)));
  print("    g(1..12) = ", vector(13,i,cf[i]));
  for(ip=1,3,
    my(p=prime(ip), TT=0);
    for(t=1,400,
      my(ok=1);
      for(n=1,M+1-t, if((cf[n]-cf[n+t])%p!=0, ok=0; break));
      if(ok, TT=t; break);
    );
    print("    coefficients of g mod ",p,": ", if(TT>0, concat("periodic, period ",TT), "not periodic with period <= 400"));
  );
  print("    g mod 2, coefficients 0..59: ", vector(60,i,lift(Mod(cf[i],2))));
);
}
quit;
