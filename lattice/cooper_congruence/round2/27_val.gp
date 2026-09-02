\\ 27_val.gp -- zero sets, divisibility sets, p-adic valuations of gamma.
read("lib.gp");
G  = [read("22_gamma_s7.txt"), read("22_gamma_s10.txt"), read("22_gamma_s18.txt")];
N  = #G[1];
print("N = ", N);
print();
print("=== exact zeros: is { n : gamma(n)=0 } exactly (special prime)*Z ? ===");
{
for(k=1,3,
  my(Z, g);
  Z = select(n->G[k][n]==0, vector(N,i,i));
  if(#Z==0, print("  ",NAM[k],":  NO zeros for n<=",N),
    g = gcd(Z);
    print("  ",NAM[k],":  #zeros=",#Z,"  gcd=",g,"  equals g*Z? ", if(#Z==N\g && Z==vector(N\g,i,g*i), "YES","NO"));
  );
);
}
print();
print("=== S_p = { n<=N : p | gamma(n) } : size, density, gcd, is it exactly gcd*Z ? ===");
{
for(k=1,3,
  for(ip=1,14,
    my(p=prime(ip), S, g, ex);
    S = select(n->G[k][n]%p==0, vector(N,i,i));
    if(#S==0,
      print("  ",NAM[k]," p=",p,":  EMPTY")
    ,
      g = gcd(S);
      ex = (g>1 && #S==N\g && S==vector(N\g,i,g*i));
      print("  ",NAM[k]," p=",p,":  #S=",#S,"  density=",strprintf("%.4f",1.0*#S/N),"  gcd(S)=",g,"  S = gcd*Z? ",if(ex,"YES","NO"));
    );
  );
);
}
print();
print("=== min_{ v_p(n)=j } v_p(gamma(n)),  j=0..5,  all rows, p<=43 ===");
{
for(k=1,3,
  for(ip=1,14,
    my(p=prime(ip), mn=vector(6,i,-1), ct=vector(6));
    for(n=1,N,
      if(G[k][n]!=0,
        my(w=valuation(n,p)+1);
        if(w<=6, ct[w]++; my(v=valuation(G[k][n],p)); if(mn[w]<0 || v<mn[w], mn[w]=v));
      );
    );
    print("  ",NAM[k]," p=",p,":  ",vector(6,i,if(ct[i]>0,mn[i],-9)),"   counts ",ct);
  );
);
}
print();
print("=== is v_p(gamma(n)) periodic in n (period<=120, ignoring gamma=0 slots)? ===");
{
for(k=1,3,
  for(ip=1,3,
    my(p=prime(ip), ok=0);
    for(T=1,120,
      my(good=1);
      for(n=1,N-T,
        if(G[k][n]!=0 && G[k][n+T]!=0 && valuation(G[k][n],p)!=valuation(G[k][n+T],p), good=0; break);
      );
      if(good, print("  ",NAM[k]," p=",p,":  v_p(gamma) periodic, period ",T); ok=1; break);
    );
    if(ok==0, print("  ",NAM[k]," p=",p,":  v_p(gamma) NOT periodic with period <= 120"));
  );
);
}
quit;
